import logging
import math
from dataclasses import dataclass

import numpy
from numpy.linalg import LinAlgError
from numpy.linalg.linalg import SVDResult
from openbabel import pybel

from .types import AutomorphismRMSD, BoolArray, Coordinate, FloatArray

Superpose3DResult = tuple[FloatArray, float, FloatArray | None]

logger = logging.getLogger("Superpose3D")


@dataclass
class MolAlignmentData:
    molecule: pybel.Molecule
    mask: BoolArray | None = None
    weights: FloatArray | float = 1.0

    # @functools.lru_cache(1)
    def get_mask(self) -> BoolArray:
        return self.mask if self.mask is not None else numpy.ones(len(self.coords()), "bool")

    # @functools.lru_cache(2)
    def coords(self, use_mask: bool = False):
        if not use_mask:
            return numpy.array([atom.coords for atom in self.molecule])
        return self.coords()[self.get_mask()]

    # @functools.lru_cache(2)
    def centroid(self, use_mask: bool = False):
        return numpy.mean(self.coords(use_mask) * self.weights, axis=0)

    # @functools.lru_cache(2)
    def centered_coords(self, use_mask: bool = False, mask_centroid: bool = False):
        return self.coords(use_mask) - self.centroid(mask_centroid)


class Superpose3D:
    def __init__(self, target: MolAlignmentData, source: MolAlignmentData) -> None:
        self.target = target
        self.source = source

    def superpose_3D(self, return_rotation_matrix: bool = False) -> Superpose3DResult:
        # The following steps come from:
        #   - http://www.pymolwiki.org/index.php/OptAlign#The_Code
        #   - http://en.wikipedia.org/wiki/Kabsch_algorithm

        centered_target = self.target.centered_coords(mask_centroid=True)
        centered_source = self.source.centered_coords(mask_centroid=True)

        result = self.perform_svd(centered_source, centered_target)
        if result is None:
            logger.warning("Couldn't perform the Single Value Decomposition, skipping alignment")
            return (self.source.coords(), 0, None)

        # check for reflections and then produce the rotation.
        # V and Wt are orthonormal, so their det's are +/-1.
        V, S, Wt = result
        reflect = numpy.linalg.det(V) * numpy.linalg.det(Wt)
        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        source_residual = numpy.sum(numpy.sum(centered_source * centered_source * self.source.weights, axis=0), axis=0)
        target_residual = numpy.sum(numpy.sum(centered_target * centered_target * self.target.weights, axis=0), axis=0)
        E0 = source_residual + target_residual
        aux_rmsd = E0 - (2.0 * sum(S))
        rmsd = numpy.sqrt(abs(aux_rmsd / len(self.source.coords(use_mask=True))))

        # rotate and translate the molecule
        rotation_matrix = numpy.dot(V, Wt)
        new_coords = numpy.dot((self.source.centered_coords()), rotation_matrix) + self.target.centroid()

        if not return_rotation_matrix:
            rotation_matrix = None

        return (new_coords, rmsd, rotation_matrix)

    def perform_svd(self, centered_source: FloatArray, centered_target: FloatArray) -> SVDResult | None:
        weights_values = [self.target.weights, 1.0]

        for weight in weights_values:
            try:
                dot_product = numpy.dot(numpy.transpose(centered_source), centered_target * weight)
                svd_result = numpy.linalg.svd(dot_product)
                return svd_result

            except LinAlgError as e:
                logger.exception(e)

        return None

    def map_to_crystal(self) -> tuple[tuple[int, int], ...]:
        """
        Some docking programs might alter the order of the atoms in the output (like Autodock Vina does...)
        this will mess up the rmsd calculation with OpenBabel.
        """
        query = pybel.ob.CompileMoleculeQuery(self.target.molecule.OBMol)
        mapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
        mapping_pose = pybel.ob.vvpairUIntUInt()
        exit = mapper.MapUnique(self.source.molecule.OBMol, mapping_pose)
        result: tuple[tuple[int, int], ...] = mapping_pose[0]
        return result

    def automorphism_rmsd(self, fit: bool) -> AutomorphismRMSD:
        mappings = pybel.ob.vvpairUIntUInt()
        raw_mappose = self.map_to_crystal()
        success = pybel.ob.FindAutomorphisms(self.target.molecule.OBMol, mappings)
        mappose = numpy.array(raw_mappose)
        sorted_indices = numpy.argsort(mappose[:, 0])
        mappose_result = mappose[sorted_indices][:, 1]
        pose_coordinates = numpy.array(self.source.coords())[mappose_result]
        lookup = [i for i, _ in enumerate(self.source.coords())]
        rmsd_result = math.inf

        # Loop through automorphisms
        target_coords = self.target.coords()
        for mapping in mappings:
            indices = ((lookup.index(x), lookup.index(y)) for x, y in mapping)
            automorph_coords = [target_coords[j] for _, j in sorted(indices)]
            mapping_rmsd = self.rmsd(pose_coordinates, automorph_coords)

            if mapping_rmsd < rmsd_result:
                rmsd_result = mapping_rmsd

            if fit:
                superpose_result = self.superpose_3D()
                fitted_pose, fitted_rmsd, _ = superpose_result

                if fitted_rmsd < rmsd_result:
                    rmsd_result = fitted_rmsd

        return (rmsd_result, fitted_pose) if fit else (rmsd_result, None)

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