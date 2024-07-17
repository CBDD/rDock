import pytest

from rdock_utils.rbhtfinder.main import main

from .conftest import EXPECTED_OUTPUT_FILE, RESULT_OUTPUT_FILE


def test_do_nothing():
    with pytest.raises(SystemExit):
        main()


def test_integration():
    # result = main()
    with open(EXPECTED_OUTPUT_FILE, "r") as expected_file, open(RESULT_OUTPUT_FILE, "r") as result_file:
        assert result_file.readlines() == expected_file.readlines()
