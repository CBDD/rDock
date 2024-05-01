from ..conftest import FIXTURES_FOLDER
from .helper import get_data, get_data_fast_mode

SDSORT_FIXTURES_FOLDER = FIXTURES_FOLDER / "sdsort"

BASIC_INPUT_FILE = str(SDSORT_FIXTURES_FOLDER / "basic_input.sdf")
FAST_MODE_INPUT_FILE = str(SDSORT_FIXTURES_FOLDER / "fast_mode_input.sdf")
FAST_MODE_INPUT_WITH_ID_FILE = str(SDSORT_FIXTURES_FOLDER / "fast_mode_input_with_id.sdf")


BASIC_INPUT_DATA = get_data(BASIC_INPUT_FILE, "SCORE")
BASIC_INPUT_DATA_DIFFERENT_FIELD = get_data(BASIC_INPUT_FILE, "test_field")
FAST_MODE_INPUT_DATA = get_data_fast_mode(FAST_MODE_INPUT_WITH_ID_FILE, "SCORE")
FAST_MODE_INPUT_DATA_DIFFERENT_FIELD = get_data_fast_mode(FAST_MODE_INPUT_WITH_ID_FILE, "test_field")
