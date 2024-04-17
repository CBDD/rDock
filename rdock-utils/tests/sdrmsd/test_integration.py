from tempfile import NamedTemporaryFile

import pytest
from openbabel import pybel

from rdock_utils.sdrmsd.main import main as sdrmsd_main
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
def test_no_fit_resulsts(main, aligned_nofit_filename: str, capsys: pytest.CaptureFixture[str]):
    with NamedTemporaryFile() as tmp:
        args = ["-o", tmp.name, REF_FILE, INPUT_FILE]
        main(args)
        assert compare_sdf_files(tmp.name, aligned_nofit_filename)
        assert capsys.readouterr().out == "POSE\tRMSD_NOFIT\n1\t0.00\n2	6.69\n3	0.48\n4	0.12\n5	0.14\n"


@parametrize_main
def test_fit_resulsts(main, aligned_fit_filename: str):
    with NamedTemporaryFile() as tmp:
        args = ["-f", "-o", tmp.name, REF_FILE, INPUT_FILE]
        main(args)
        assert compare_sdf_files(tmp.name, aligned_fit_filename)


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
