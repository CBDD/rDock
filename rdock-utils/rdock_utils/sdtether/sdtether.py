import itertools
import logging

import numpy
from openbabel import pybel

from rdock_utils.common import AutomorphismRMSD, BoolArray, MolAlignmentData, Superpose3D, update_coordinates
from rdock_utils.sdtether.parser import SDTetherConfig, get_config

logger = logging.getLogger("SDTether")


def atoms_string(id_list: list[int], batch_size: int = 35) -> str:
    str_ids = map(str, id_list)
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
        self.ref_matches_align_data = self._get_ref_matches_align_data()

    def _get_ref_matches_align_data(self) -> list[MolAlignmentData]:
        ref = next(pybel.readfile("sdf", self.refsdf))
        refMatchIds = self.get_matches(ref, raise_error=True)
        ref_matches_align_data = self.get_matches_align_data(ref, refMatchIds)
        return ref_matches_align_data

    def match_to_mask(self, match: list[int], num_atoms: int) -> BoolArray:
        return numpy.array([i in match for i in range(1, num_atoms + 1)], dtype=numpy.bool_)

    def get_matches_align_data(self, molecule: pybel.Molecule, matches: list[list[int]]) -> list[MolAlignmentData]:
        masks = (self.match_to_mask(match, len(molecule.atoms)) for match in matches)
        align_data_list = [MolAlignmentData(molecule, mask) for mask in masks]
        return align_data_list

    def get_matches(self, molecule: pybel.Molecule, raise_error: bool = False) -> list[list[int]]:
        match_ids = self.smarts.findall(molecule)
        num_matches = len(match_ids)
        if num_matches == 0:
            logger.info("No match")
            if raise_error:
                raise RuntimeError(
                    "No match found in the reference structure and the SMARTS string given. Please check it."
                )
        elif num_matches == 1:
            logger.info("Match")
        elif num_matches > 1:
            logger.info(
                f"More than one match in the molecule '{molecule.title}' for the SMARTS string given. "
                f"Will tether each input molecule all possible ways ({num_matches})."
            )
        return match_ids

    def run(self):
        # Do the same for molecule in molsdf
        with pybel.Outputfile("sdf", self.outsdf, overwrite=True) as out:
            molSupp = pybel.readfile("sdf", self.molsdf)
            # pybel.ob.OBForceField_FindForceField("MMFF94")
            for i, mol in enumerate(molSupp, start=1):
                logger.info(f"Processing molecule {i}")
                self.process_molecule(mol, out)

        logger.info("DONE")

    def process_molecule(self, mol: pybel.Molecule, out: pybel.Outputfile) -> None:
        # find all matches for the smarts query in the input molecule
        mol_matches = self.get_matches(mol)
        for mol_match in mol_matches:
            mol_alignment_data = MolAlignmentData(mol, self.match_to_mask(mol_match, len(mol.atoms)))
            alignement_data_pairs = itertools.product([mol_alignment_data], self.ref_matches_align_data)
            superposers = (Superpose3D(target, source) for source, target in alignement_data_pairs)
            for superposer in superposers:
                alignment_result = superposer.automorphism_rmsd(fit=True)
                self.write_aligned_molecule(mol, mol_match, alignment_result, out)

    def write_aligned_molecule(
        self,
        mol: pybel.Molecule,
        match_ids: list[tuple[int, ...]],
        aligned_data: AutomorphismRMSD,
        out: pybel.Outputfile,
    ) -> None:
        rmsd, coords = aligned_data
        update_coordinates(mol, coords)
        newData = pybel.ob.OBPairData()
        newData.SetAttribute("TETHERED ATOMS")
        newData.SetValue(atoms_string(match_ids))

        mol.OBMol.DeleteData("TETHERED ATOMS")  # Remove Previous DATA
        mol.OBMol.CloneData(newData)  # Add new data

        out.write(mol)


def main(argv: list[str] | None = None):
    config = get_config(argv)
    sdtether = SDTether(config)
    sdtether.run()
    # run
    # Read reference pose and get atom list matching smarts query
    # if more than 1 match, take the first one


if __name__ == "__main__":
    main()
