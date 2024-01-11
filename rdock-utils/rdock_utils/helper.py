import argparse
import logging
import math
import sys

import numpy
from numpy.linalg import LinAlgError
from numpy.typing import ArrayLike
from openbabel import pybel

logger = logging.getLogger("sdrmsd")

Coordinate = tuple[float, float, float]
SingularValueDecomposition = tuple[ArrayLike[float], ArrayLike[float], ArrayLike[float]]
RMSDResult = float | tuple[float, ArrayLike[float]]


class Helper:
    def __init__(
        self,
        reference_sdf: argparse.FileType,
        input_sdf: argparse.FileType,
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
        reference: ArrayLike[float],
        target: ArrayLike[float],
        weights: ArrayLike[float] | None = None,
        reference_mask: ArrayLike[bool] | None = None,
        target_mask: ArrayLike[bool] | None = None,
        return_rotation_matrix: bool = False,
    ) -> tuple[ArrayLike[float], float, ArrayLike[float]] | tuple[ArrayLike[float], float]:
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
        except LinAlgError:
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

    def perform_svd(
        self, reference_tmp: ArrayLike[float], target_tmp: ArrayLike[float], weights: ArrayLike[float] | int
    ) -> SingularValueDecomposition | None:
        try:
            dot_product = numpy.dot(numpy.transpose(reference_tmp), target_tmp * weights)
            svd_result = numpy.linalg.svd(dot_product)
        except LinAlgError:
            svd_result = self._handle_svd_linalg_error(reference_tmp, target_tmp)
        return svd_result

    def _handle_svd_linalg_error(
        self, reference_tmp: ArrayLike[float], target_tmp: ArrayLike[float]
    ) -> SingularValueDecomposition | None:
        try:
            dot_product = numpy.dot(numpy.transpose(reference_tmp), target_tmp)
            svd_result = numpy.linalg.svd(dot_product)
            return svd_result
        except LinAlgError:
            raise

    def rmsd(self, all_coordinates_A: list[Coordinate], all_coordinates_B: list[Coordinate]) -> float:
        """
        Find the root mean square deviation between two lists of 3-tuples.
        """
        deviation = sum(
            self.squared_distance(atom_A, atom_B) for atom_A, atom_B in zip(all_coordinates_A, all_coordinates_B)
        )
        return math.sqrt(deviation / len(all_coordinates_A))

    def squared_distance(self, coordinates_A: Coordinate, coordinates_B: Coordinate) -> float:
        """
        Find the squared distance between two 3-tuples.
        """
        return sum((a - b) ** 2 for a, b in zip(coordinates_A, coordinates_B))

    def map_to_crystal(self, xtal: pybel.Molecule, pose: pybel.Molecule) -> tuple[int, int]:
        """
        Some docking programs might alter the order of the atoms in the output (like Autodock Vina does...)
        this will mess up the rmsd calculation with OpenBabel.
        """
        query = pybel.ob.CompileMoleculeQuery(xtal.OBMol)
        mapper: pybel.ob.OBIsomorphismMapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
        mapping_pose = pybel.ob.vvpairUIntUInt()
        exit = mapper.MapUnique(pose.OBMol, mapping_pose)
        return mapping_pose[0]

    def update_coordinates(self, obmol: pybel.Molecule, new_coordinates: ArrayLike[float]) -> None:
        """
        Update OBMol coordinates. new_coordinates is a numpy array.
        """
        for i, atom in enumerate(obmol):
            atom.OBAtom.SetVector(*new_coordinates[i])
