from dataclasses import dataclass, field
from typing import Any

# from typing import Literal
import numpy
from openbabel import pybel

Coordinate = tuple[float, float, float]
FloatArray = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
BoolArray = numpy.ndarray[Any, numpy.dtype[numpy.bool_]]
CoordsArray = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
AutomorphismRMSD = tuple[float, CoordsArray | None]
Vector3D = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
Matrix3x3 = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
SingularValueDecomposition = tuple[Matrix3x3, Vector3D, Matrix3x3]
Superpose3DResult = tuple[CoordsArray, float, Matrix3x3]

## Shape support for type hinting is not yet avaialable in numpy
## let's keep this as a guide for numpy 2.0 release
# Coordinate = tuple[float, float, float]
# FloatArray = numpy.ndarray[Literal["N"], numpy.dtype[float]]
# BoolArray = numpy.ndarray[Literal["N"], numpy.dtype[bool]]
# CoordsArray = numpy.ndarray[Literal["N", 3], numpy.dtype[float]]
# AutomorphismRMSD = tuple[float, CoordsArray | None]
# Vector3D = numpy.ndarray[Literal[3], numpy.dtype[float]]
# Matrix3x3 = numpy.ndarray[Literal[3, 3], numpy.dtype[float]]
# SingularValueDecomposition = tuple[Matrix3x3, Vector3D, Matrix3x3]
# Superpose3DResult = tuple[CoordsArray, float, Matrix3x3]


@dataclass
class SDRMSDData:
    skipped: list[int] = field(default_factory=list)
    molecules_dict: dict[int, pybel.Molecule] = field(default_factory=dict)  # Save all poses with their dockid
    population: dict[int, int] = field(default_factory=dict)  # Poses to be written
    out_dict: dict[int, tuple[pybel.Molecule, float]] = field(default_factory=dict)


@dataclass
class PoseMatchData:
    pose_index: int
    docked_pose: pybel.Molecule
    sdrmsd_data: SDRMSDData
