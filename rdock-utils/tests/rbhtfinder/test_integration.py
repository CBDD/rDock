from pathlib import Path
from typing import Callable

import pytest

from rdock_utils.rbhtfinder.main import main as rbhtfinder_main
from rdock_utils.rbhtfinder_original_copy import main as rbhtfinder_old_main

from .conftest import EXPECTED_OUTPUT_FILE, get_file_content

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
def test_integration(main: Callable[[list[str]], None], file_path: Path, argv: list[str]):
    main(argv)
    result = get_file_content(file_path)
    expected_result = get_file_content(EXPECTED_OUTPUT_FILE)
    assert result == expected_result
