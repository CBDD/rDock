from tempfile import NamedTemporaryFile

import pytest

from rdock_utils.sdrmsd_OOP import main as sdrmsd_main
from rdock_utils.sdrmsd_original import main as sdrmsd3_main
from tests.sdrmsd.conftest import INPUT_FILE, REF_FILE

parametrize_main = pytest.mark.parametrize(
    "main",
    [
        pytest.param(sdrmsd_main, id="Improved version (Python 3.12)"),
        pytest.param(sdrmsd3_main, id="Original version (Python 3)"),
    ],
)


@parametrize_main
def test_do_nothing(main):
    with pytest.raises(SystemExit):
        main()


@parametrize_main
def test_basic_run(main):
    args = [REF_FILE, INPUT_FILE]
    main(args)


@parametrize_main
def test_no_fit_resulsts(main, aligned_molecules_raw_data: bytes):
    with NamedTemporaryFile() as tmp:
        args = ["-o", tmp.name, REF_FILE, INPUT_FILE]
        main(args)
        with open(tmp.name, "rb") as f:
            output = f.read()
        assert output == aligned_molecules_raw_data


@parametrize_main
def test_fit_resulsts(main, aligned_fit_molecules_raw_data: bytes):
    with NamedTemporaryFile() as tmp:
        args = ["-f", "-o", tmp.name, REF_FILE, INPUT_FILE]
        main(args)
        with open(tmp.name, "rb") as f:
            output = f.read()
        assert output == aligned_fit_molecules_raw_data
