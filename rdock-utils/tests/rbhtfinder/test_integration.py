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
def test_do_nothing(main):
    with pytest.raises(SystemExit):
        main()


@parametrize_main
def test_integration(main, file_path, argv):
    main(argv)
    result = get_file_content(file_path)
    expected_result = get_file_content(EXPECTED_OUTPUT_FILE)
    assert result == expected_result
