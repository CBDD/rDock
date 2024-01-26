from openbabel import pybel

from .types import CoordsArray


def update_coordinates(molecule: pybel.Molecule, new_coordinates: CoordsArray) -> pybel.Molecule:
    for i, atom in enumerate(molecule):
        atom.OBAtom.SetVector(*new_coordinates[i])

    return molecule
