from typing import Any

import numpy
import numpy.typing

# TODO: Review common types for all rdock_utils scripts
# SDRMSD types
FloatArray = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
CoordsArray = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
AutomorphismRMSD = tuple[float, CoordsArray | None]
Vector3D = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
Matrix3x3 = numpy.ndarray[Any, numpy.dtype[numpy.float64]]
SingularValueDecomposition = tuple[Matrix3x3, Vector3D, Matrix3x3]
Superpose3DResult = tuple[CoordsArray, float, Matrix3x3]
AtomsMapping = tuple[tuple[int, int], ...]

# RBHTFinder types
SDReportArray = numpy.ndarray[list[int | str | float], numpy.dtype[numpy.object_]]
Array1DFloat = numpy.typing.NDArray[numpy.float_]
Array2DFloat = numpy.typing.NDArray[numpy.float_]
Array3DFloat = numpy.typing.NDArray[numpy.float_]
Array1DStr = numpy.typing.NDArray[numpy.str_]
Array1DInt = numpy.typing.NDArray[numpy.int_]
ColumnNamesArray = Array1DStr | list[str]
InputData = tuple[SDReportArray, ColumnNamesArray]
MinScoreIndices = dict[int, Array1DInt]
FilterCombination = tuple[float, float]

## Shape support for type hinting is not yet avaialable in np
## let's keep this as a guide for np 2.0 release
# FloatArray = numpy.ndarray[Literal["N"], numpy.dtype[float]]
# BoolArray = numpy.ndarray[Literal["N"], numpy.dtype[bool]]
# CoordsArray = numpy.ndarray[Literal["N", 3], numpy.dtype[float]]
# AutomorphismRMSD = tuple[float, CoordsArray | None]
# Vector3D = numpy.ndarray[Literal[3], numpy.dtype[float]]
# Matrix3x3 = numpy.ndarray[Literal[3, 3], numpy.dtype[float]]
# SingularValueDecomposition = tuple[Matrix3x3, Vector3D, Matrix3x3]
# Superpose3DResult = tuple[CoordsArray, float, Matrix3x3]
