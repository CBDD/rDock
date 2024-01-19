import functools
import logging
import math
import sys
from dataclasses import dataclass, field

import numpy
from openbabel import pybel

logger = logging.getLogger("sdrmsd")

Coordinate = tuple[float, float, float]
SingularValueDecomposition = tuple[numpy.ndarray[float], numpy.ndarray[float], numpy.ndarray[float]]
AutomorphismRMSD = float | tuple[float, numpy.ndarray[float]]
Superpose3D = tuple[numpy.ndarray[float], float, numpy.ndarray[float]] | tuple[numpy.ndarray[float], float]


@dataclass
class SDRMSDData:
    skipped: list[int] = field(default_factory=list)
    molecules_dict: dict[int, pybel.Molecule] = field(default_factory=dict)  # Save all poses with their dockid
    population: dict[int, int] = field(default_factory=dict)  # Poses to be written
    out_dict: dict[int, tuple[pybel.Molecule, AutomorphismRMSD]] = field(default_factory=dict)


@dataclass
class PoseMatchData:
    pose_index: int
    docked_pose: pybel.Molecule
    sdrmsd_data: SDRMSDData


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

    def superpose_3D(
        self,
        reference: numpy.ndarray[float],
        target: numpy.ndarray[float],
        weights: numpy.ndarray[float] | None = None,
        reference_mask: numpy.ndarray[bool] | None = None,
        target_mask: numpy.ndarray[bool] | None = None,
        return_rotation_matrix: bool = False,
    ) -> Superpose3D:
        """superpose_3D performs 3d superposition using a weighted Kabsch algorithm : http://dx.doi.org/10.1107%2FS0567739476001873 & doi: 10.1529/biophysj.105.066654
        definition : superpose3D(reference, target, weights,refmask,target_mask)
        @parameter 1 :  reference - xyz coordinates of the reference structure (the ligand for instance)
        @type 1 :       float64 numpy array (nx3)
        ---
        @parameter 2 :  target - theoretical target positions to which we should move (does not need to be physically relevant.
        @type 2 :       float 64 numpy array (nx3)
        ---
        @parameter 3:   weights - numpy array of atom weights (usuallly between 0 and 1)
        @type 3 :       float 64 numpy array (n)
        @parameter 4:   mask - a numpy boolean mask for designating atoms to include
        Note reference and target positions must have the same dimensions -> n*3 numpy arrays where n is the number of points (or atoms)
        returns:
        - Tuple containing new coordinates and RMSD (default behavior).
        OR
        - Tuple containing new coordinates, RMSD, and rotation matrix (if return_rotation_matrix is True).
        """
        weights = weights or 1.0
        reference_mask = reference_mask or numpy.ones(len(reference), "bool")
        target_mask = target_mask or numpy.ones(len(target), "bool")
        # First get the centroid of both states
        reference_centroid = numpy.mean(reference[reference_mask] * weights, axis=0)
        # Print reference_centroid
        reference_centered_coords = reference - reference_centroid
        # Print reference_centered_coords
        target_centroid = numpy.mean(target[target_mask] * weights, axis=0)
        target_centered_coords = target - target_centroid
        # Print target_centered_coords
        # The following steps come from : http://www.pymolwiki.org/index.php/OptAlign#The_Code and http://en.wikipedia.org/wiki/Kabsch_algorithm
        # Initial residual, see Kabsch.
        E0 = numpy.sum(
            numpy.sum(
                reference_centered_coords[reference_mask] * reference_centered_coords[reference_mask] * weights, axis=0
            ),
            axis=0,
        ) + numpy.sum(
            numpy.sum(target_centered_coords[target_mask] * target_centered_coords[target_mask] * weights, axis=0),
            axis=0,
        )
        reference_tmp = numpy.copy(reference_centered_coords[reference_mask])
        target_tmp = numpy.copy(target_centered_coords[target_mask])
        # print reference_centered_coords[reference_mask]
        # single value decomposition of the dotProduct of both position vectors
        try:
            V, S, Wt = self.perform_svd(reference_tmp, target_tmp, weights)
        except numpy.linalg.LinAlgError:
            warning_msg = "Couldn't perform the Single Value Decomposition, skipping alignment"
            logger.warning(warning_msg)
            print(warning_msg, file=sys.stderr)
            return (reference, 0)
        # we already have our solution, in the results from SVD.
        # we just need to check for reflections and then produce
        # the rotation.  V and Wt are orthonormal, so their det's
        # are +/-1.
        reflect = float(numpy.linalg.det(V) * numpy.linalg.det(Wt))
        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]
        rmsd = E0 - (2.0 * sum(S))
        rmsd = numpy.sqrt(abs(rmsd / len(reference[reference_mask])))  # get the rmsd
        # U is simply V*Wt
        U = numpy.dot(V, Wt)  # get the rotation matrix
        # rotate and translate the molecule
        new_coords = numpy.dot((reference_centered_coords), U) + target_centroid  # translate & rotate
        # new_coords=(reference_centered_coords + target_centroid)
        # print U
        return (new_coords, rmsd, U) if return_rotation_matrix else (new_coords, rmsd)

    def get_automorphism_rmsd(self, target: pybel.Molecule, molecule: pybel.Molecule) -> AutomorphismRMSD:
        """
        Use Automorphism to reorder target coordinates to match reference coordinates atom order
        for correct RMSD comparison. Only the lowest RMSD will be returned.

        Returns:
        If fit=False: bestRMSD	(float)					RMSD of the best matching mapping.
        If fit=True:  (bestRMSD, molecCoordinates)	(float, npy.array)	RMSD of best match and its molecule fitted coordinates.
        """
        fit = self.fit
        mappings = pybel.ob.vvpairUIntUInt()
        lookup = [i for i, _ in enumerate(target)]
        success = pybel.ob.FindAutomorphisms(target.OBMol, mappings)

        target_coordinates = [atom.coords for atom in target]  # TODO: ADAPTED TO ORIGINAL
        # index_to_target_coordinates = dict(enumerate(target_coordinates))
        mappose = numpy.array(self.map_to_crystal(target, molecule))
        sorted_indices = numpy.argsort(mappose[:, 0])
        mappose_result = mappose[sorted_indices][:, 1]
        molecule_coordinates = [atom.coords for atom in molecule]
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
            if fit:
                fitted_pose, fitted_rmsd = self.superpose_3D(
                    numpy.array(automorph_coords), numpy.array(pose_coordinates)
                )

                # Update result if the fitted RMSD is lower
                if fitted_rmsd < rmsd_result:
                    rmsd_result = fitted_rmsd

        return (rmsd_result, fitted_pose) if fit else rmsd_result

    def perform_svd(
        self,
        reference_tmp: numpy.ndarray[float],
        target_tmp: numpy.ndarray[float],
        weights: numpy.ndarray[float] | int,
    ) -> SingularValueDecomposition | None:
        try:
            dot_product = numpy.dot(numpy.transpose(reference_tmp), target_tmp * weights)
            svd_result = numpy.linalg.svd(dot_product)
        except numpy.linalg.LinAlgError:
            svd_result = self._handle_svd_linalg_error(reference_tmp, target_tmp)
        return svd_result

    @staticmethod
    def _handle_svd_linalg_error(
        reference_tmp: numpy.ndarray[float], target_tmp: numpy.ndarray[float]
    ) -> SingularValueDecomposition | None:
        try:
            dot_product = numpy.dot(numpy.transpose(reference_tmp), target_tmp)
            svd_result = numpy.linalg.svd(dot_product)
            return svd_result
        except numpy.linalg.LinAlgError:
            raise

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

    def update_coordinates(self, obmol: pybel.Molecule, new_coordinates: numpy.ndarray[float]) -> None:
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
