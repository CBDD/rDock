from typing import Any

import numpy as np
import numpy.typing as npt

# TODO: Review common types for all rdock_utils scripts
# SDRMSD types
FloatArray = np.ndarray[Any, np.dtype[np.float64]]
CoordsArray = np.ndarray[Any, np.dtype[np.float64]]
AutomorphismRMSD = tuple[float, CoordsArray | None]
Vector3D = np.ndarray[Any, np.dtype[np.float64]]
Matrix3x3 = np.ndarray[Any, np.dtype[np.float64]]
SingularValueDecomposition = tuple[Matrix3x3, Vector3D, Matrix3x3]
Superpose3DResult = tuple[CoordsArray, float, Matrix3x3]
AtomsMapping = tuple[tuple[int, int], ...]

# RBHTFinder types
SDReportArray = np.ndarray[list[int | str | float], np.dtype[np.object_]]
Array1DFloat = npt.NDArray[np.float_]
Array2DFloat = npt.NDArray[np.float_]
Array3DFloat = npt.NDArray[np.float_]
Array1DStr = npt.NDArray[np.str_]
Array1DInt = npt.NDArray[np.int_]
ColumnNamesArray = Array1DStr | list[str]
InputData = tuple[SDReportArray, ColumnNamesArray]
MinScoreIndices = dict[int, Array1DInt]
FilterCombination = tuple[float, float]

## Shape support for type hinting is not yet avaialable in np
## let's keep this as a guide for np 2.0 release
# FloatArray = np.ndarray[Literal["N"], np.dtype[float]]
# BoolArray = np.ndarray[Literal["N"], np.dtype[bool]]
# CoordsArray = np.ndarray[Literal["N", 3], np.dtype[float]]
# AutomorphismRMSD = tuple[float, CoordsArray | None]
# Vector3D = np.ndarray[Literal[3], np.dtype[float]]
# Matrix3x3 = np.ndarray[Literal[3, 3], np.dtype[float]]
# SingularValueDecomposition = tuple[Matrix3x3, Vector3D, Matrix3x3]
# Superpose3DResult = tuple[CoordsArray, float, Matrix3x3]
