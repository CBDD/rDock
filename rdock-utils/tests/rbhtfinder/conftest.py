from pathlib import Path

import pytest

from ..conftest import FIXTURES_FOLDER

RBHTFINDER_FIXTURES_FOLDER = FIXTURES_FOLDER / "rbhtfinder"

INPUT_FILE = str(RBHTFINDER_FIXTURES_FOLDER / "input.txt")
EXPECTED_THRESHOLD_FILE = str(RBHTFINDER_FIXTURES_FOLDER / "threshold.txt")
EXPECTED_OUTPUT_FILE = str(RBHTFINDER_FIXTURES_FOLDER / "output.txt")


@pytest.fixture
def output_temp(tmp_path: Path) -> Path:
    output_path = tmp_path / "output.txt"
    return output_path


@pytest.fixture
def threshold_temp(tmp_path: Path) -> Path:
    threshold_path = tmp_path / "threshold.txt"
    return threshold_path


@pytest.fixture
def argv(output_temp: Path, threshold_temp: Path) -> list[str]:
    return [
        "-i",
        INPUT_FILE,
        "-o",
        str(output_temp),
        "-t",
        str(threshold_temp),
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


def get_file_content(file: str | Path) -> str:
    with open(file, "r") as f:
        return f.read()
