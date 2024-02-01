import logging
import math
from dataclasses import dataclass
from typing import Iterable

import numpy
from numpy.linalg import LinAlgError
from numpy.linalg.linalg import SVDResult
from openbabel import pybel

from .types import (
    AutomorphismRMSD,
    BoolArray,
    CoordsArray,
    FloatArray,
    MatchIds,
    Matrix3x3,
    Superpose3DResult,
    Vector3D,
)

logger = logging.getLogger("Superpose3D")


def update_coordinates(molecule: pybel.Molecule, new_coordinates: CoordsArray) -> pybel.Molecule:
    for i, atom in enumerate(molecule):
        atom.OBAtom.SetVector(*new_coordinates[i])

    return molecule


@dataclass
class MolAlignmentData:
    molecule: pybel.Molecule
    mask: BoolArray | None = None
    weights: FloatArray | float = 1.0

    def get_mask(self) -> BoolArray:
        return self.mask if self.mask is not None else numpy.ones(len(self.coords()), "bool")

    def coords(self, use_mask: bool = False) -> CoordsArray:
        mask = self.get_mask() if use_mask else None
        return self._masked_coords(mask)

    def _masked_coords(self, mask: Iterable[int] | None) -> CoordsArray:
        coords = numpy.array([atom.coords for atom in self.molecule])
        return coords[mask] if mask is not None else coords

    def centroid(self, use_mask: bool = False) -> Vector3D:
        return numpy.mean(self.coords(use_mask) * self.weights, axis=0)  # type: ignore

    def centered_coords(self, use_mask: bool = False, mask_centroid: bool = False) -> CoordsArray:
        return self.coords(use_mask) - self.centroid(mask_centroid)


class Superpose3D:
    def __init__(self, target: MolAlignmentData, source: MolAlignmentData) -> None:
        self.target = target
        self.source = source

    def align(self, ref_mask: BoolArray | None = None, target_mask: BoolArray | None = None) -> Superpose3DResult:
        ref_use_mask = ref_mask is not None
        target_use_mask = target_mask is not None
        ref = self.source if not ref_use_mask else MolAlignmentData(self.source.molecule, ref_mask)
        target = self.target if not target_use_mask else MolAlignmentData(self.target.molecule, target_mask)
        try:
            # Get rotation matrix
            ref_align_selection = ref.centered_coords(use_mask=ref_use_mask, mask_centroid=True)
            target_align_selection = target.centered_coords(use_mask=target_use_mask, mask_centroid=True)
            rmsd, rotation_matrix = self.get_rotation_matrix(ref_align_selection, target_align_selection)

            # Apply transformation
            all_coords = ref.centered_coords(use_mask=False, mask_centroid=True)
            translation_vector = target.centroid(target_use_mask)
            new_coords = self.apply_transformation(all_coords, rotation_matrix, translation_vector)

            if ref_use_mask:
                target_coordinates = target.coords(use_mask=False)
                rmsd = self.rmsd(new_coords, target_coordinates)

            return (new_coords, rmsd, rotation_matrix)

        except ValueError:
            logger.warning("Couldn't perform the Single Value Decomposition, skipping alignment")
            return (ref.coords(), 0, numpy.identity(3))

    def apply_transformation(
        self, coords: CoordsArray, rotation_matrix: Matrix3x3, translation_vector: Vector3D
    ) -> CoordsArray:
        return numpy.dot(coords, rotation_matrix) + translation_vector

    def get_rotation_matrix(self, source: CoordsArray, target: CoordsArray) -> tuple[float, Matrix3x3]:
        # The following steps come from:
        #   - http://www.pymolwiki.org/index.php/OptAlign#The_Code
        #   - http://en.wikipedia.org/wiki/Kabsch_algorithm
        svd = self.perform_svd(source, target)

        if svd is None:
            raise ValueError

        # check for reflections and then produce the rotation.
        # V and Wt are orthonormal, so their det's are +/-1.
        V, S, Wt = svd
        reflect = numpy.linalg.det(V) * numpy.linalg.det(Wt)

        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        source_residual = numpy.sum(numpy.sum(source * source * self.source.weights, axis=0), axis=0)
        target_residual = numpy.sum(numpy.sum(target * target * self.target.weights, axis=0), axis=0)
        E0 = source_residual + target_residual
        aux_rmsd = E0 - (2.0 * sum(S))
        rmsd = numpy.sqrt(abs(aux_rmsd / len(self.source.coords(use_mask=True))))

        # rotate and translate the molecule
        rotation_matrix = numpy.dot(V, Wt)

        return (rmsd, rotation_matrix)

    def perform_svd(self, centered_source: CoordsArray, centered_target: CoordsArray) -> SVDResult | None:
        weights_values = [self.target.weights, 1.0]

        for weight in weights_values:
            try:
                dot_product = numpy.dot(numpy.transpose(centered_source), centered_target * weight)
                svd_result: SVDResult | None = numpy.linalg.svd(dot_product)
                return svd_result

            except LinAlgError as e:
                logger.exception(e)

        return None

    def map_to_crystal(self) -> MatchIds:
        query = pybel.ob.CompileMoleculeQuery(self.target.molecule.OBMol)
        mapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
        mapping_pose = pybel.ob.vvpairUIntUInt()
        mapper.MapUnique(self.source.molecule.OBMol, mapping_pose)
        result: MatchIds = mapping_pose[0]
        return result

    def get_pose_rmsd(self, fit: bool, pose_coordinates: CoordsArray, mapping: MatchIds) -> AutomorphismRMSD:
        target_coords = self.target.coords()
        coords_mask = [j for _, j in sorted(mapping)]
        automorph_coords = target_coords[coords_mask]
        rmsd = self.rmsd(pose_coordinates, automorph_coords)
        fitted_pose = None

        if fit:
            superpose_result = self.align()
            fitted_pose, fitted_rmsd, _ = superpose_result
            rmsd = min(fitted_rmsd, rmsd)

        return (rmsd, fitted_pose)

    def automorphism_rmsd(self, fit: bool) -> AutomorphismRMSD:
        mappings = pybel.ob.vvpairUIntUInt()
        raw_mappose = self.map_to_crystal()
        pybel.ob.FindAutomorphisms(self.target.molecule.OBMol, mappings)
        mappose = numpy.array(raw_mappose)
        sorted_indices = numpy.argsort(mappose[:, 0])
        mappose_result = mappose[sorted_indices][:, 1]
        pose_coordinates = numpy.array(self.source.coords())[mappose_result]
        return min(
            (self.get_pose_rmsd(fit, pose_coordinates, mapping) for mapping in mappings),
            key=lambda t: t[0],
            default=(math.inf, None),
        )

    def rmsd(self, all_coordinates_1: CoordsArray, all_coordinates_2: CoordsArray) -> float:
        differences = all_coordinates_2 - all_coordinates_1
        deviation = sum(delta.dot(delta) for delta in differences)
        return math.sqrt(deviation / len(all_coordinates_1))
