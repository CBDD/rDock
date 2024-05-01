import itertools
from io import StringIO
from typing import Iterable

import pytest

from rdock_utils.common import read_molecules
from rdock_utils.sdsort.parser import get_config
from rdock_utils.sdsort.sdsort import SDSort

from .conftest import BASIC_INPUT_FILE, FAST_MODE_INPUT_FILE, FAST_MODE_INPUT_WITH_ID_FILE


def get_data(path: str, key_field: str = "_TITLE1", value_field: str = "SCORE") -> list[tuple[str, str]]:
    with open(path, "r") as f:
        molecules = read_molecules(f)
        return [(m.get_field(key_field), m.get_field(value_field)) for m in molecules]


def sorted_grouped_values(
    values: Iterable[str],
    reverse: bool = False,
    numeric: bool = False,
) -> list[tuple[str, str]]:
    convert = float if numeric else str
    return list(
        itertools.chain.from_iterable(
            sorted(values, key=lambda x: convert(x[1]), reverse=reverse)
            for _, values in itertools.groupby(values, key=lambda x: x[0])
        )
    )


@pytest.mark.parametrize("numeric", [pytest.param(False, id="lexicographic"), pytest.param(True, id="numeric")])
@pytest.mark.parametrize("reverse", [pytest.param(False, id="ascending"), pytest.param(True, id="descending")])
@pytest.mark.parametrize(
    "sorting_field",
    [pytest.param(["-f", "test_field"], id="custom field"), pytest.param([], id="default field")],
)
def test_sdsort_basic(numeric: bool, reverse: bool, sorting_field: list[str]):
    args = sorting_field + (["-n"] if numeric else []) + (["-r"] if reverse else []) + [BASIC_INPUT_FILE]
    config = get_config(args)
    output = StringIO()
    sdsort = SDSort(config, output)
    data = get_data(BASIC_INPUT_FILE, value_field=config.sorting_field)
    convert = float if numeric else str
    expected_result = sorted(data, key=lambda x: convert(x[1]), reverse=reverse)

    sdsort.run()

    output.seek(0)
    molecules = read_molecules(output)
    result = [(m.get_field(config.group_key), m.get_field(config.sorting_field)) for m in molecules]

    assert result == expected_result


@pytest.mark.parametrize("numeric", [pytest.param(False, id="lexicographic"), pytest.param(True, id="numeric")])
@pytest.mark.parametrize("reverse", [pytest.param(False, id="ascending"), pytest.param(True, id="descending")])
@pytest.mark.parametrize(
    "sorting_field",
    [pytest.param(["-f", "test_field"], id="custom field"), pytest.param([], id="default field")],
)
def test_sdsort_fast_mode(numeric: bool, reverse: bool, sorting_field: list[str]):
    args = ["-s"] + sorting_field + (["-n"] if numeric else []) + (["-r"] if reverse else []) + [FAST_MODE_INPUT_FILE]
    config = get_config(args)
    output = StringIO()
    sdsort = SDSort(config, output)
    data = get_data(FAST_MODE_INPUT_FILE, value_field=config.sorting_field)
    expected_result = sorted_grouped_values(data, reverse=reverse, numeric=numeric)

    sdsort.run()

    output.seek(0)
    molecules = read_molecules(output)
    result = [(m.get_field(config.group_key), m.get_field(config.sorting_field)) for m in molecules]

    assert result == expected_result


@pytest.mark.parametrize("numeric", [pytest.param(False, id="lexicographic"), pytest.param(True, id="numeric")])
@pytest.mark.parametrize("reverse", [pytest.param(False, id="ascending"), pytest.param(True, id="descending")])
@pytest.mark.parametrize(
    "group_key",
    [pytest.param(["--group-key", "_TITLE2"], id="custom group field"), pytest.param([], id="default group field")],
)
def test_sdsort_by_group_key(numeric: bool, reverse: bool, group_key: list[str]):
    filename = FAST_MODE_INPUT_WITH_ID_FILE
    args = ["-s"] + group_key + (["-n"] if numeric else []) + (["-r"] if reverse else []) + [filename]
    config = get_config(args)
    output = StringIO()
    sdsort = SDSort(config, output)
    data = get_data(filename, key_field=config.group_key, value_field=config.sorting_field)
    expected_result = sorted_grouped_values(data, reverse=reverse, numeric=numeric)

    sdsort.run()

    output.seek(0)
    molecules = read_molecules(output)
    result = [(m.get_field(config.group_key), m.get_field(config.sorting_field)) for m in molecules]

    assert result == expected_result
