import pytest

from ..conftest import FIXTURES_FOLDER

RBHTFINDER_FIXTURES_FOLDER = FIXTURES_FOLDER / "rbhtfinder"

INPUT_FILE = str(RBHTFINDER_FIXTURES_FOLDER / "input_tabs.txt")
THRESHOLD_FILE = str(RBHTFINDER_FIXTURES_FOLDER / "threshold.txt")
EXPECTED_OUTPUT_FILE = str(RBHTFINDER_FIXTURES_FOLDER / "output.txt")


@pytest.fixture
def file_path(tmp_path):
    output_path = tmp_path / "output.txt"
    return output_path


@pytest.fixture
def argv(file_path):
    return [
        "-i",
        INPUT_FILE,
        "-o",
        str(file_path),
        "-t",
        THRESHOLD_FILE,
        "-f",
        "column=4,steps=3,min=-10.0,max=0.0,interval=5.0",
        "column=6,steps=5,min=1.0,max=5.0,interval=5.0",
        "--max-time",
        "1",
        "--min-perc",
        "1.0",
        "-v",
        "5",
        "--header",
    ]


def get_file_content(file: str) -> list[str]:
    with open(file, "r") as f:
        return f.readlines()
