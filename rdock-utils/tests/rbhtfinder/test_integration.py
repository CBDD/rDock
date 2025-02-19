from pathlib import Path
from typing import Callable

import pytest

from rdock_utils.rbhtfinder.main import main as rbhtfinder_main
from rdock_utils.rbhtfinder_original_copy import main as rbhtfinder_old_main

from .conftest import EXPECTED_OUTPUT_FILE, EXPECTED_THRESHOLD_FILE, get_file_content

parametrize_main = pytest.mark.parametrize(
    "main",
    [
        pytest.param(rbhtfinder_old_main, id="Original version Python 3"),
        pytest.param(rbhtfinder_main, id="Improved version Python 3.12"),
    ],
)


@parametrize_main
def test_do_nothing(main: Callable[[list[str]], None]):
    with pytest.raises(SystemExit):
        main()


@parametrize_main
def test_integration(main: Callable[[list[str]], None], output_temp: Path, threshold_temp: Path, argv: list[str]):
    main(argv)
    output = get_file_content(output_temp)
    threshold = get_file_content(threshold_temp)
    expected_output = get_file_content(EXPECTED_OUTPUT_FILE)
    expected_threshold = get_file_content(EXPECTED_THRESHOLD_FILE)
    assert output == expected_output
    assert threshold == expected_threshold
