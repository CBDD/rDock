from tempfile import NamedTemporaryFile

import pytest
from openbabel import pybel

from rdock_utils.sdtether.sdtether import main as sdtether_main
from rdock_utils.sdtether_original import main as sdtether_old_main
from tests.sdtether.conftest import EXPECTED_OUTPUT_FILE_1, INPUT_FILE, REF_FILE

parametrize_main = pytest.mark.parametrize(
    "main",
    [
        pytest.param(sdtether_main, id="Improved version (Python 3.12)"),
        pytest.param(sdtether_old_main, id="Original version (Python 3)"),
    ],
)


@parametrize_main
def test_do_nothing(main):
    with pytest.raises(SystemExit):
        main()


@parametrize_main
def test_basic_run(main):
    with NamedTemporaryFile() as tmp:
        args = [REF_FILE, INPUT_FILE, tmp.name, "cnc"]
        main(args)
        assert compare_sdf_files(tmp.name, EXPECTED_OUTPUT_FILE_1)


def atoms_are_equal(atom_1: pybel.Atom, atom_2: pybel.Atom) -> bool:
    if atom_1.atomicnum != atom_2.atomicnum:
        return False

    if atom_1.coords != atom_2.coords:
        return False

    if atom_1.formalcharge != atom_2.formalcharge:
        return False

    if atom_1.type != atom_2.type:
        return False

    return True


def molecules_are_equal(molecule_1: pybel.Molecule, molecule_2: pybel.Molecule) -> bool:
    molecule_1_data = {k: v for k, v in molecule_1.data.items()}
    molecule_2_data = {k: v for k, v in molecule_2.data.items()}

    if molecule_1_data != molecule_2_data:
        return False

    return all(atoms_are_equal(atom_1, atom_2) for atom_1, atom_2 in zip(molecule_1.atoms, molecule_2.atoms))


def compare_sdf_files(filename_1: str, filename_2: str) -> bool:
    molecules_1 = pybel.readfile("sdf", filename_1)
    molecules_2 = pybel.readfile("sdf", filename_2)
    return all(molecules_are_equal(molecule_1, molecule_2) for molecule_1, molecule_2 in zip(molecules_1, molecules_2))
