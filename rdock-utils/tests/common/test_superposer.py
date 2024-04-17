import numpy
import pytest
from openbabel import pybel

from rdock_utils.common import CoordsArray, MolAlignmentData

from ..conftest import FIXTURES_FOLDER

SDTETHER_FIXTURES_FOLDER = FIXTURES_FOLDER / "sdtether"
REF_FILE = str(SDTETHER_FIXTURES_FOLDER / "ref.sdf")
INPUT_FILE = str(SDTETHER_FIXTURES_FOLDER / "ref.sdf")
EXPECTED_OUTPUT_FILE_1 = str(SDTETHER_FIXTURES_FOLDER / "out.sdf")


@pytest.mark.parametrize(
    "match, expected",
    [
        pytest.param(
            (4, 5, 6),
            numpy.array([[0.14633, -1.2568, 0.0], [-0.42316, 0.0222, 0.0], [0.27683, 1.2346, 0.0]]),
            id="match_1",
        ),
        pytest.param(
            (6, 5, 13),
            numpy.array([[0.92310, 0.71123, 0.0], [0.22310, -0.50116, 0.0], [-1.14619, -0.21006, 0.0]]),
            id="match_2",
        ),
    ],
)
def test_get_centroid(match: tuple[int], expected: CoordsArray):
    mol = next(pybel.readfile("sdf", REF_FILE))
    mask = numpy.array(match)
    data = MolAlignmentData(mol, mask=mask)
    result_coords = data.centered_coords(use_mask=True, mask_centroid=True)
    assert numpy.allclose(result_coords, expected, atol=0.0001)
