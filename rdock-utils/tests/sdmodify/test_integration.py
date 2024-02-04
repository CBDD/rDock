from io import StringIO

import pytest

from rdock_utils.common import read_molecules
from rdock_utils.sdmodify import main

from .conftest import INPUT_FILE


def test_do_nothing():
    with pytest.raises(SystemExit):
        main()


@pytest.mark.parametrize(
    "args, expected_titles",
    [
        pytest.param(
            ["-f", "test_field", INPUT_FILE],
            (["0.0", "0.0", "2.0", "3.0", "4.0", "0.0"]),
            id="molecule field filter",
        ),
        pytest.param(
            ["-f", "_REC", INPUT_FILE],
            list(map(str, range(1, 7))),
            id="molecule field filter with output file",
        ),
    ],
)
def test_basic_run(args: list[str], expected_titles: list[str], capsys: pytest.CaptureFixture):
    main(args)
    captured = capsys.readouterr()
    input = StringIO(captured.out)
    molecules = read_molecules(input)
    titles = [m.title for m in molecules]
    assert titles == expected_titles
