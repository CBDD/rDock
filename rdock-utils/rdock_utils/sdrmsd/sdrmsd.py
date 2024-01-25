import functools
import logging
import math

from openbabel import pybel

from .superpose3d import MolAlignmentData, Superpose3D
from .types import AutomorphismRMSD, CoordsArray, PoseMatchData, SDRMSDData

logger = logging.getLogger("SDRMSD")


class SDRMSD:
    def __init__(
        self,
        reference_sdf: str,
        input_sdf: str,
        fit: bool,
        threshold: float,
        output_filename: str,
    ) -> None:
        self.reference_sdf = reference_sdf
        self.input_sdf = input_sdf
        self.fit = fit
        self.threshold = threshold
        self.output_filename = output_filename

    def run(self) -> None:
        # Find the RMSD between the crystal pose and each docked pose
        docked_poses = pybel.readfile("sdf", self.input_sdf)
        self.display_fit_message()
        data = SDRMSDData()

        # Read crystal pose
        crystal_pose = self.get_crystal_pose()
        crystal_atoms = len(crystal_pose.atoms)

        for i, docked_pose in enumerate(docked_poses, start=1):
            atoms_number = self.process_docked_pose(docked_pose)

            if atoms_number != crystal_atoms:
                data.skipped.append(i)
                continue

            rmsd_result = self.calculate_rmsd(crystal_pose, docked_pose)
            pose_match_data = PoseMatchData(i, docked_pose, data)
            self.handle_pose_matching(rmsd_result, pose_match_data)

        if self.output_filename:
            self.process_and_save_selected_molecules(data)

        if data.skipped:
            logger.warning(f"SKIPPED input molecules due to the number of atom mismatch: {data.skipped}")

    def get_automorphism_rmsd(self, target: pybel.Molecule, molecule: pybel.Molecule) -> AutomorphismRMSD:
        """
        Use Automorphism to reorder target coordinates to match reference coordinates atom order
        for correct RMSD comparison. Only the lowest RMSD will be returned.

        Returns:
        If fit=False: bestRMSD	(float)					RMSD of the best matching mapping.
        If fit=True:  (bestRMSD, molecCoordinates)	(float, npy.array)	RMSD of best match and its molecule fitted coordinates.
        """
        superposer = Superpose3D(MolAlignmentData(molecule), MolAlignmentData(target))
        result = superposer.automorphism_rmsd(self.fit)
        return result

    def update_coordinates(self, obmol: pybel.Molecule, new_coordinates: CoordsArray) -> None:
        """
        Update OBMol coordinates. new_coordinates is a numpy array.
        """
        for i, atom in enumerate(obmol):
            atom.OBAtom.SetVector(*new_coordinates[i])

    def get_best_matching_pose(self, pose_match_data: PoseMatchData) -> tuple[int | None, float]:
        docked_pose = pose_match_data.docked_pose
        molecules_dict = pose_match_data.sdrmsd_data.molecules_dict
        get_rmsd = functools.partial(self.get_automorphism_rmsd, target=docked_pose)
        poses_rmsd = ((index, get_rmsd(molecule)[0]) for index, molecule in molecules_dict.items())
        filtered_by_threshold = (t for t in poses_rmsd if t[1] < self.threshold)
        return min(filtered_by_threshold, key=lambda t: t[1], default=(None, math.inf))

    def process_and_save_selected_molecules(self, data: SDRMSDData) -> None:
        with pybel.Outputfile("sdf", self.output_filename, overwrite=True) as output_sdf:
            for i in sorted(data.out_dict.keys()):
                molecule, rmsd_result = data.out_dict[i]
                # Get the number of matches in the thresholding operation
                population_value = data.population.get(i, 1)
                self.save_molecule_with_rmsd(output_sdf, molecule, rmsd_result, population_value)

    def save_molecule_with_rmsd(
        self, output_sdf: pybel.Outputfile, molecule: pybel.Molecule, rmsd: float, population_value: int
    ) -> None:
        new_data = pybel.ob.OBPairData()
        new_data.SetAttribute("RMSD")
        new_data.SetValue(f"{rmsd:.3f}")

        if self.threshold:
            pop_data = pybel.ob.OBPairData()
            pop_data.SetAttribute("Population")
            pop_data.SetValue(f"{population_value}")
            molecule.OBMol.CloneData(pop_data)

        molecule.OBMol.CloneData(new_data)
        output_sdf.write(molecule)

    def get_crystal_pose(self) -> pybel.Molecule:
        """
        Read crystal pose file and remove hydrogen atoms. Returns crystal pose molecule.
        """
        crystal_pose = next(pybel.readfile("sdf", self.reference_sdf))
        crystal_pose.removeh()
        return crystal_pose

    def display_fit_message(self) -> None:
        message = "FIT" if self.fit else "NOFIT"
        print(f"POSE\tRMSD_{message}")

    def process_docked_pose(self, docked_pose: pybel.Molecule) -> int:
        """
        Remove hydrogen atoms and return the number of atoms in docked pose molecule.
        """
        docked_pose.removeh()
        return len(docked_pose.atoms)

    def calculate_rmsd(self, crystal: pybel.Molecule, docked_pose: pybel.Molecule) -> float:
        """
        Perform RMSD calculations and update coordinates if required.
        """
        rmsd, fitted_coords = self.get_automorphism_rmsd(crystal, docked_pose)

        if self.fit:
            if fitted_coords is not None:
                self.update_coordinates(docked_pose, fitted_coords)
            else:
                logger.warning("Automorphism failed. skipping alignment")

        return rmsd

    def handle_pose_matching(self, rmsd: float, pose_match_data: PoseMatchData) -> None:
        """
        Function to handle pose matching and filtering based on 'threshold' parser argument.
        """
        if self.threshold:
            match_pose, best_match_value = self.get_best_matching_pose(pose_match_data)
            if match_pose is not None:
                self.print_matching_info(pose_match_data, match_pose, best_match_value)
            else:
                self.save_and_print_info(rmsd, pose_match_data)
        else:
            self.save_and_print_info(rmsd, pose_match_data)

    def print_matching_info(self, pose_match_data: PoseMatchData, match_pose: int, best_match_value: float) -> None:
        index = pose_match_data.pose_index
        population = pose_match_data.sdrmsd_data.population
        logger.info(f"Pose {index} matches pose {match_pose} with {best_match_value:.3f} RMSD")
        population[match_pose] += 1

    def save_and_print_info(self, rmsd: float, pose_match_data: PoseMatchData) -> None:
        """
        Function to save and print information based on 'out' parser argument.
        """
        index = pose_match_data.pose_index
        docked_pose = pose_match_data.docked_pose
        out_dict = pose_match_data.sdrmsd_data.out_dict
        molecules_dict = pose_match_data.sdrmsd_data.molecules_dict
        population = pose_match_data.sdrmsd_data.population

        if self.output_filename:
            out_dict[index] = (docked_pose, rmsd)
        print(f"{index}\t{rmsd:.2f}")
        molecules_dict[index] = docked_pose
        population[index] = 1
