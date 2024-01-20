import functools
import logging
import math

import numpy
from openbabel import pybel

from .superpose3d import Superpose3D
from .types import AutomorphismRMSD, Coordinate, FloatArray, PoseMatchData, SDRMSDData

logger = logging.getLogger("SDRMSD")


class SDRMSD:
    def __init__(
        self,
        reference_sdf: str,
        input_sdf: str,
        fit: bool,
        threshold: float,
        out: bool,
    ) -> None:
        self.reference_sdf = reference_sdf
        self.input_sdf = input_sdf
        self.fit = fit
        self.threshold = threshold
        self.out = out

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

        if self.out:
            output_sdf = pybel.Outputfile("sdf", self.out, overwrite=True)
            self.process_and_save_selected_molecules(output_sdf, data)

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
        mappings = pybel.ob.vvpairUIntUInt()
        lookup = [i for i, _ in enumerate(target)]
        success = pybel.ob.FindAutomorphisms(target.OBMol, mappings)

        molecule_coordinates = [atom.coords for atom in molecule]
        target_coordinates = [atom.coords for atom in target]  # TODO: ADAPTED TO ORIGINAL
        raw_mappose = self.map_to_crystal(target, molecule)
        # index_to_target_coordinates = dict(enumerate(target_coordinates))

        # take to superpose3d v v v
        mappose = numpy.array(raw_mappose)
        sorted_indices = numpy.argsort(mappose[:, 0])
        mappose_result = mappose[sorted_indices][:, 1]
        pose_coordinates = numpy.array(molecule_coordinates)[mappose_result]
        rmsd_result = math.inf

        # Loop through automorphisms
        for mapping in mappings:
            automorph_coords = [None] * len(target_coordinates)  # TODO: ADAPTED TO ORIGINAL
            for x, y in mapping:  # TODO: ADAPTED TO ORIGINAL
                automorph_coords[lookup.index(x)] = target_coordinates[lookup.index(y)]  # TODO: ADAPTED TO ORIGINAL

            mapping_rmsd = self.rmsd(pose_coordinates, automorph_coords)

            # Update result if the current mapping has a lower RMSD
            if mapping_rmsd < rmsd_result:
                rmsd_result = mapping_rmsd

            # Additional fitting if fit=True
            if self.fit:
                superposer = Superpose3D(numpy.array(automorph_coords))
                superpose_result = superposer.superpose_3D(numpy.array(pose_coordinates))
                fitted_pose, fitted_rmsd, _ = superpose_result

                # Update result if the fitted RMSD is lower
                if fitted_rmsd < rmsd_result:
                    rmsd_result = fitted_rmsd

        # ^ ^ ^
        return (rmsd_result, fitted_pose) if self.fit else rmsd_result

    def rmsd(self, all_coordinates_1: list[Coordinate], all_coordinates_2: list[Coordinate]) -> float:
        """
        Find the root mean square deviation between two lists of 3-tuples.
        """
        deviation = sum(
            self.squared_distance(atom_1, atom_2) for atom_1, atom_2 in zip(all_coordinates_1, all_coordinates_2)
        )
        return math.sqrt(deviation / len(all_coordinates_1))

    def squared_distance(self, coordinates_1: Coordinate, coordinates_2: Coordinate) -> float:
        """
        Find the squared distance between two 3-tuples.
        """
        return sum((a - b) ** 2 for a, b in zip(coordinates_1, coordinates_2))

    def get_automorphism_rmsd_2(self, target: pybel.Molecule, molecule: pybel.Molecule) -> AutomorphismRMSD:
        """
        Use Automorphism to reorder target coordinates to match reference coordinates atom order
        for correct RMSD comparison. Only the lowest RMSD will be returned.

        Returns:
        If fit=False: bestRMSD	(float)					RMSD of the best matching mapping.
        If fit=True:  (bestRMSD, molecCoordinates)	(float, npy.array)	RMSD of best match and its molecule fitted coordinates.
        """
        mappings = pybel.ob.vvpairUIntUInt()
        lookup = [i for i, _ in enumerate(target)]
        success = pybel.ob.FindAutomorphisms(target.OBMol, mappings)

        molecule_coordinates = [atom.coords for atom in molecule]
        target_coordinates = [atom.coords for atom in target]  # TODO: ADAPTED TO ORIGINAL
        raw_mappose = self.map_to_crystal(target, molecule)
        # index_to_target_coordinates = dict(enumerate(target_coordinates))

        # take to superpose3d v v v
        superposer = Superpose3D()

    def map_to_crystal(self, xtal: pybel.Molecule, pose: pybel.Molecule) -> tuple[int, int]:
        """
        Some docking programs might alter the order of the atoms in the output (like Autodock Vina does...)
        this will mess up the rmsd calculation with OpenBabel.
        """
        query = pybel.ob.CompileMoleculeQuery(xtal.OBMol)
        mapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
        mapping_pose = pybel.ob.vvpairUIntUInt()
        exit = mapper.MapUnique(pose.OBMol, mapping_pose)
        return mapping_pose[0]

    def update_coordinates(self, obmol: pybel.Molecule, new_coordinates: FloatArray) -> None:
        """
        Update OBMol coordinates. new_coordinates is a numpy array.
        """
        for i, atom in enumerate(obmol):
            atom.OBAtom.SetVector(*new_coordinates[i])

    def get_best_matching_pose(self, pose_match_data: PoseMatchData) -> tuple[int | None, float]:
        docked_pose = pose_match_data.docked_pose
        molecules_dict = pose_match_data.sdrmsd_data.molecules_dict
        get_rmsd = functools.partial(self.get_automorphism_rmsd, target=docked_pose)

        poses_rmsd = ((index, get_rmsd(molecule)) for index, molecule in molecules_dict.items())
        filtered_by_threshold = (t for t in poses_rmsd if t[1] < self.threshold)

        return min(filtered_by_threshold, key=lambda t: t[1], default=(None, math.inf))

    def process_and_save_selected_molecules(self, output_sdf: pybel.Outputfile, data: SDRMSDData):
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

        molecule.OBMol.CloneData(new_data)  # Add new data
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

    def calculate_rmsd(self, crystal: pybel.Molecule, docked_pose: pybel.Molecule) -> AutomorphismRMSD:
        """
        Perform RMSD calculations and update coordinates if required.
        """

        if self.fit:
            rmsd_result, fitted_result = self.get_automorphism_rmsd(crystal, docked_pose)
            self.update_coordinates(docked_pose, fitted_result)
        else:
            rmsd_result = self.get_automorphism_rmsd(crystal, docked_pose)

        return rmsd_result

    def handle_pose_matching(self, rmsd_result: AutomorphismRMSD, pose_match_data: PoseMatchData) -> None:
        """
        Function to handle pose matching and filtering based on 'threshold' parser argument.
        """
        if self.threshold:
            match_pose, best_match_value = self.get_best_matching_pose(pose_match_data)
            if match_pose is not None:
                self.print_matching_info(pose_match_data, match_pose, best_match_value)
            else:
                self.save_and_print_info(rmsd_result, pose_match_data)
        else:
            self.save_and_print_info(rmsd_result, pose_match_data)

    def print_matching_info(self, pose_match_data: PoseMatchData, match_pose: int, best_match_value: float) -> None:
        index = pose_match_data.pose_index
        population = pose_match_data.sdrmsd_data.population
        logger.info(f"Pose {index} matches pose {match_pose} with {best_match_value:.3f} RMSD")
        population[match_pose] += 1

    def save_and_print_info(self, rmsd_result: AutomorphismRMSD, pose_match_data: PoseMatchData) -> None:
        """
        Function to save and print information based on 'out' parser argument.
        """
        index = pose_match_data.pose_index
        docked_pose = pose_match_data.docked_pose
        out_dict = pose_match_data.sdrmsd_data.out_dict
        molecules_dict = pose_match_data.sdrmsd_data.molecules_dict
        population = pose_match_data.sdrmsd_data.population

        if self.out:
            out_dict[index] = (docked_pose, rmsd_result)
        print(f"{index}\t{rmsd_result:.2f}")
        molecules_dict[index] = docked_pose
        population[index] = 1
