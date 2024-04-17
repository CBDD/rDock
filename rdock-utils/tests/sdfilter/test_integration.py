from io import StringIO

import pytest

from rdock_utils.common import read_molecules
from rdock_utils.sdfilter.main import main

from .conftest import INPUT_FILE


def generate_titles(ids: list[int]) -> list[str]:
    return [f"MOL{i}" for i in ids]


def test_do_nothing():
    with pytest.raises(SystemExit):
        main()


@pytest.mark.parametrize(
    "args, expected_titles",
    [
        pytest.param(
            ["-f", "$test_field <= 2", INPUT_FILE],
            generate_titles([1, 2, 3, 6]),
            id="molecule field filter",
        ),
        pytest.param(
            ["-f", "$_COUNT == 1", "-s", "test_field", INPUT_FILE],
            generate_titles([1, 3, 4, 5]),
            id="summary field filter",
        ),
        pytest.param(["-f", "$_REC == 3", INPUT_FILE], generate_titles([3]), id="record filter"),
    ],
)
def test_basic_run(args: list[str], expected_titles: list[str], capsys: pytest.CaptureFixture):
    main(args)
    captured = capsys.readouterr()
    output = StringIO(captured.out)
    molecules = read_molecules(output)
    titles = [m.title for m in molecules]
    assert titles == expected_titles
