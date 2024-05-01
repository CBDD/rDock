from io import StringIO

import pytest

from rdock_utils.sdsort.parser import get_config
from rdock_utils.sdsort.sdsort import SDSort

from .conftest import (
    BASIC_INPUT_DATA,
    BASIC_INPUT_DATA_DIFFERENT_FIELD,
    BASIC_INPUT_FILE,
    FAST_MODE_INPUT_DATA,
    FAST_MODE_INPUT_DATA_DIFFERENT_FIELD,
    FAST_MODE_INPUT_FILE,
)
from .helper import generate_result


@pytest.mark.parametrize(
    "args, expected_result",
    [
        pytest.param(
            [BASIC_INPUT_FILE],
            generate_result([1, 2, 6, 5, 4, 3], BASIC_INPUT_DATA),
            id="Default sort (Lexicographic by '-f' argument)",
        ),
        pytest.param(
            ["-f", "test_field", BASIC_INPUT_FILE],
            generate_result([6, 1, 3, 4, 5, 2], BASIC_INPUT_DATA_DIFFERENT_FIELD),
            id="Default sort different field",
        ),
        pytest.param(
            ["-r", BASIC_INPUT_FILE],
            generate_result([3, 4, 5, 6, 2, 1], BASIC_INPUT_DATA),
            id="Reverse sort",
        ),
        pytest.param(
            ["-n", BASIC_INPUT_FILE],
            generate_result([2, 1, 6, 4, 3, 5], BASIC_INPUT_DATA),
            id="Numeric sort",
        ),
        pytest.param(
            ["-n", "-r", BASIC_INPUT_FILE],
            generate_result([5, 3, 4, 6, 1, 2], BASIC_INPUT_DATA),
            id="Reverse numeric sort",
        ),
        pytest.param(
            ["-s", FAST_MODE_INPUT_FILE],
            generate_result([2, 1, 4, 3, 5, 6, 8, 9, 7, 10], FAST_MODE_INPUT_DATA, True),
            id="Default fast mode",
        ),
        pytest.param(
            ["-s", FAST_MODE_INPUT_FILE],
            generate_result([2, 1, 4, 3, 5, 6, 8, 9, 7, 10], FAST_MODE_INPUT_DATA, True),
            id="Default fast mode",
        ),
        pytest.param(
            ["-s", "-f", "test_field", FAST_MODE_INPUT_FILE],
            generate_result([2, 1, 3, 4, 5, 6, 7, 8, 9, 10], FAST_MODE_INPUT_DATA_DIFFERENT_FIELD, True),
            id="Default fast mode sort different field",
        ),
        pytest.param(
            ["-s", "--group-key", "id", FAST_MODE_INPUT_FILE],
            generate_result([2, 4, 3, 8, 1, 9, 7, 5, 6, 10], FAST_MODE_INPUT_DATA, True),
            id="Default fast mode sort different group key",
        ),
    ],
)
def test_basic_run(args: list[str], expected_result: list[str]):
    config = get_config(args)
    output = StringIO()
    sdsort = SDSort(config, output)

    sdsort.run()
    output.seek(0)
    output_lines = output.getvalue().strip().splitlines()
    titles = [line for line in output_lines if line.startswith("MOL")]
    field_values = [
        output_lines[i + 1].strip()
        for i, line in enumerate(output_lines)
        if line.startswith(f">  <{config.sorting_field}>")
    ]
    result = list(zip(titles, field_values))

    assert result == expected_result
