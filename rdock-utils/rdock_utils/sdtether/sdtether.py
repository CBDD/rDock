import itertools
import logging
from typing import Iterable

import numpy
from openbabel import pybel

from rdock_utils.common import CoordsArray, MolAlignmentData, Superpose3D, update_coordinates
from rdock_utils.sdtether.parser import SDTetherConfig

logger = logging.getLogger("SDTether")


def atoms_string(atom_ids: Iterable[int], batch_size: int = 35) -> str:
    str_ids = map(str, atom_ids)
    batches = itertools.batched(str_ids, batch_size)
    comma_formatted = (",".join(batch) for batch in batches)
    atoms_string = ",\n".join(comma_formatted) + "\n"
    return atoms_string


class SDTether:
    def __init__(self, config: SDTetherConfig) -> None:
        self.refsdf = config.reference_filename
        self.molsdf = config.input_filename
        self.outsdf = config.output_filename
        self.smarts = pybel.Smarts(config.smarts)
        self.ref_mol = next(pybel.readfile("sdf", self.refsdf))
        self.ref_matches = self.get_matches(self.ref_mol)
        self.ref_align_data = MolAlignmentData(self.ref_mol)

    def run(self):
        self.show_matches_info()

        with pybel.Outputfile("sdf", self.outsdf, overwrite=True) as out:
            molSupp = pybel.readfile("sdf", self.molsdf)
            for i, mol in enumerate(molSupp, start=1):
                logger.info(f"Processing molecule {i}")
                print(f"## Molecule {i}")
                self.process_molecule(mol, out)

        print("DONE")

    def show_matches_info(self):
        num_matches = len(self.ref_matches)

        if num_matches == 0:
            print("No match")
            message = "No match found in the reference structure and the SMARTS string given. Please check it."
            raise RuntimeError(message)

        elif num_matches == 1:
            print("Match")

        elif num_matches > 1:
            print(
                "More than one match in the reference molecule for the SMARTS string given. "
                "Will tether each input molecule all possible ways."
            )

    def process_molecule(self, mol: pybel.Molecule, out: pybel.Outputfile) -> None:
        mol_align_data = MolAlignmentData(mol)
        mol_matches = self.get_matches(mol)
        superposer = Superpose3D(mol_align_data, self.ref_align_data)

        for mol_match_index, mol_match in enumerate(mol_matches):
            mol_mask = numpy.array(mol_match) - 1
            for ref_match_index, ref_match in enumerate(self.ref_matches):
                ref_mask = numpy.array(ref_match) - 1
                new_coords, _, _ = superposer.align(mol_mask, ref_mask)
                selection_rmsd = superposer.rmsd(new_coords[mol_mask], self.ref_align_data.coords()[ref_mask])
                self.write_aligned_molecule(mol, mol_match, new_coords, out)
                print(
                    f"\tBest RMSD reached (match {mol_match_index}, refmatch {ref_match_index}): {selection_rmsd:.4f}"
                )

    def get_matches(self, molecule: pybel.Molecule) -> list[tuple[int, ...]]:
        match_ids = self.smarts.findall(molecule)
        return match_ids

    def write_aligned_molecule(
        self,
        mol: pybel.Molecule,
        match_ids: list[tuple[int, ...]],
        new_coords: CoordsArray,
        out: pybel.Outputfile,
    ) -> None:
        new_mol = mol.clone
        update_coordinates(new_mol, new_coords)
        newData = pybel.ob.OBPairData()
        newData.SetAttribute("TETHERED ATOMS")
        newData.SetValue(atoms_string(match_ids))
        new_mol.OBMol.DeleteData("TETHERED ATOMS")  # Remove Previous DATA
        new_mol.OBMol.CloneData(newData)  # Add new data
        out.write(new_mol)
