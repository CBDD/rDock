from pathlib import Path
from unittest import mock

import pytest

from rdock_utils.common import read_molecules
from rdock_utils.sdsplit import main, sanitize_args

from .conftest import INPUT_FILE


@pytest.mark.parametrize(
    "args, expected_molecule_count",
    [
        pytest.param(["-2", INPUT_FILE], [2, 2, 2], id="Deprecated (2 mols per file)"),
        pytest.param(["-r", "2", INPUT_FILE], [2, 2, 2], id="Non-deprecated (2 mols per file)"),
        pytest.param(["-r", "3", INPUT_FILE], [3, 3], id="Non-deprecated (3 mols per file)"),
        pytest.param(["-r", "4", INPUT_FILE], [4, 2], id="Non-deprecated (4 mols per file)"),
    ],
)
def test_basic_run(args: list[str], expected_molecule_count: list[int], tmp_path: Path):
    output_root = str(tmp_path / "test_")
    output_args = ["-o", output_root]
    full_args = output_args + args
    main(full_args)

    for i, molecule_count in enumerate(expected_molecule_count):
        current_file = Path(f"{output_root}{i}.sd")
        assert current_file.is_file()
        with open(current_file, "r") as f:
            assert sum(1 for _ in read_molecules(f)) == molecule_count


@pytest.mark.parametrize(
    "args, gives_warning",
    [
        pytest.param(["-4"], True, id="Deprecated args"),
        pytest.param(["-r", "4"], False, id="accepted args"),
    ],
)
@mock.patch("rdock_utils.sdsplit.logger")
def test_sanitize_args(m_logger: mock.Mock, args: list[str], gives_warning: bool):
    sanitize_args(args)
    assert (m_logger.warning.call_count > 0) == gives_warning
