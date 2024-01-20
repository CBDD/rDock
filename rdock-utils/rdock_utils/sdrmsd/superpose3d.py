import logging
import math

import numpy
from numpy.linalg import LinAlgError
from numpy.linalg.linalg import SVDResult
from openbabel import pybel

from .types import AutomorphismRMSD, Coordinate, FloatArray

Superpose3DResult = tuple[FloatArray, float, FloatArray] | tuple[FloatArray, float]

logger = logging.getLogger("Superpose3D")


class Superpose3D:
    def __init__(self, reference: FloatArray) -> None:
        self.reference = reference

    def superpose_3D(
        self,
        target: FloatArray,
        weights: FloatArray | float = 1.0,
        reference_mask: numpy.ndarray[bool] | None = None,
        target_mask: numpy.ndarray[bool] | None = None,
        return_rotation_matrix: bool = False,
    ) -> Superpose3DResult:
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
        reference_mask = reference_mask or numpy.ones(len(self.reference), "bool")
        target_mask = target_mask or numpy.ones(len(target), "bool")
        # First get the centroid of both states
        reference_centroid = numpy.mean(self.reference[reference_mask] * weights, axis=0)
        # Print reference_centroid
        reference_centered_coords = self.reference - reference_centroid
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

        result = self.perform_svd(reference_tmp, target_tmp, weights)
        if result is None:
            logger.warning("Couldn't perform the Single Value Decomposition, skipping alignment")
            return (self.reference, 0, None)

        V, S, Wt = result

        # we already have our solution, in the results from SVD.
        # we just need to check for reflections and then produce
        # the rotation.  V and Wt are orthonormal, so their det's
        # are +/-1.
        reflect = float(numpy.linalg.det(V) * numpy.linalg.det(Wt))
        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]
        rmsd = E0 - (2.0 * sum(S))
        rmsd = numpy.sqrt(abs(rmsd / len(self.reference[reference_mask])))  # get the rmsd
        # U is simply V*Wt
        rotation_matrix = numpy.dot(V, Wt)  # get the rotation matrix
        # rotate and translate the molecule
        new_coords = numpy.dot((reference_centered_coords), rotation_matrix) + target_centroid  # translate & rotate
        # new_coords=(reference_centered_coords + target_centroid)
        # print U
        if not return_rotation_matrix:
            rotation_matrix = None
        return (new_coords, rmsd, rotation_matrix)

    def perform_svd(
        self,
        reference_tmp: FloatArray,
        target_tmp: FloatArray,
        weights: FloatArray | float,
    ) -> SVDResult | None:
        weights_values = [weights, 1.0]

        for weight in weights_values:
            try:
                dot_product = numpy.dot(numpy.transpose(reference_tmp), target_tmp * weight)
                svd_result = numpy.linalg.svd(dot_product)
                return svd_result
            except LinAlgError as e:
                logger.exception(e)

        return None

    def automorphism_rmsd(
        self,
        fit: bool,
        raw_mappose: tuple[int, int],
        molecule_coordinates: list[tuple],
        mappings: pybel.ob.vvpairUIntUInt,
        target_coordinates: list[tuple],
        lookup,
    ) -> AutomorphismRMSD:
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
            if fit:
                superpose_result = self.superpose_3D(numpy.array(pose_coordinates))
                fitted_pose, fitted_rmsd, _ = superpose_result

                # Update result if the fitted RMSD is lower
                if fitted_rmsd < rmsd_result:
                    rmsd_result = fitted_rmsd

        return (rmsd_result, fitted_pose) if fit else rmsd_result

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
