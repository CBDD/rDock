from typing import Any

import numpy

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
