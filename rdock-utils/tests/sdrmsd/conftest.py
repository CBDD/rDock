import pytest

from ..conftest import FIXTURES_FOLDER

SDRMSD_FIXTURES_FOLDER = FIXTURES_FOLDER / "sdrmsd"
REF_FILE = str(SDRMSD_FIXTURES_FOLDER / "ref.sdf")
INPUT_FILE = str(SDRMSD_FIXTURES_FOLDER / "input.sdf")


@pytest.fixture(scope="session")
def aligned_nofit_filename() -> str:
    return str(SDRMSD_FIXTURES_FOLDER / "aligned-nofit.sdf")


@pytest.fixture(scope="session")
def aligned_fit_filename() -> str:
    return str(SDRMSD_FIXTURES_FOLDER / "aligned-fit.sdf")


@pytest.fixture(scope="session")
def aligned_molecules_raw_data(aligned_nofit_filename: str) -> bytes:
    with open(aligned_nofit_filename, "rb") as f:
        return f.read()


@pytest.fixture(scope="session")
def aligned_fit_molecules_raw_data(aligned_fit_filename: str) -> bytes:
    with open(aligned_fit_filename, "rb") as f:
        return f.read()
