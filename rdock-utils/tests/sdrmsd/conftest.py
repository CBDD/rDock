import pytest

from ..conftest import FIXTURES_FOLDER

SDRMSD_FIXTURES_FOLDER = FIXTURES_FOLDER / "sdrmsd"
REF_FILE = str(SDRMSD_FIXTURES_FOLDER / "ref.sdf")
INPUT_FILE = str(SDRMSD_FIXTURES_FOLDER / "input.sdf")
ALIGNED_FILE = str(SDRMSD_FIXTURES_FOLDER / "aligned-nofit.sdf")
ALIGNED_FIT_FILE = str(SDRMSD_FIXTURES_FOLDER / "aligned-fit.sdf")


@pytest.fixture(scope="session")
def aligned_molecules_raw_data() -> bytes:
    with open(ALIGNED_FILE, "rb") as f:
        return f.read()

@pytest.fixture(scope="session")
def aligned_fit_molecules_raw_data() -> bytes:
    with open(ALIGNED_FIT_FILE, "rb") as f:
        return f.read()
