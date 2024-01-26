import itertools
import logging
import math
import sys

import numpy
from nptyping import NDArray
from openbabel import pybel

from rdock_utils.common.superpose3d import MolAlignmentData
from rdock_utils.common.types import BoolArray, CoordsArray
from rdock_utils.sdtether.parser import SDTetherConfig, get_config

logger = logging.getLogger("SDTether")


def superpose3D(ref, target, weights=None, refmask=None, targetmask=None, returnRotMat=False):
    """superpose3D performs 3d superposition using a weighted Kabsch algorithm : http://dx.doi.org/10.1107%2FS0567739476001873 & doi: 10.1529/biophysj.105.066654
    definition : superpose3D(ref, target, weights,refmask,targetmask)
    @parameter 1 :  ref - xyz coordinates of the reference structure (the ligand for instance)
    @type 1 :       float64 numpy array (nx3)
    ---
    @parameter 2 :  target - theoretical target positions to which we should move (does not need to be physically relevant.
    @type 2 :       float 64 numpy array (nx3)
    ---
    @parameter 3:   weights - numpy array of atom weights (usuallly between 0 and 1)
    @type 3 :       float 64 numpy array (n)
    @parameter 4:   mask - a numpy boolean mask for designating atoms to include
    Note ref and target positions must have the same dimensions -> n*3 numpy arrays where n is the number of points (or atoms)
    Returns a set of new coordinates, aligned to the target state as well as the rmsd
    """
    if weights is None:
        weights = 1.0
    if refmask is None:
        refmask = numpy.ones(len(ref), "bool")
    if targetmask is None:
        targetmask = numpy.ones(len(target), "bool")
    # first get the centroid of both states
    ref_centroid = numpy.mean(ref[refmask] * weights, axis=0)
    # print ref_centroid
    refCenteredCoords = ref - ref_centroid
    # print refCenteredCoords
    target_centroid = numpy.mean(target[targetmask] * weights, axis=0)
    targetCenteredCoords = target - target_centroid
    # print targetCenteredCoords
    # the following steps come from : http://www.pymolwiki.org/index.php/OptAlign#The_Code and http://en.wikipedia.org/wiki/Kabsch_algorithm
    # Initial residual, see Kabsch.
    E0 = numpy.sum(
        numpy.sum(refCenteredCoords[refmask] * refCenteredCoords[refmask] * weights, axis=0), axis=0
    ) + numpy.sum(
        numpy.sum(targetCenteredCoords[targetmask] * targetCenteredCoords[targetmask] * weights, axis=0), axis=0
    )
    reftmp = numpy.copy(refCenteredCoords[refmask])
    targettmp = numpy.copy(targetCenteredCoords[targetmask])
    # print refCenteredCoords[refmask]
    # single value decomposition of the dotProduct of both position vectors
    try:
        dotProd = numpy.dot(numpy.transpose(reftmp), targettmp * weights)
        V, S, Wt = numpy.linalg.svd(dotProd)
    except Exception:
        try:
            dotProd = numpy.dot(numpy.transpose(reftmp), targettmp)
            V, S, Wt = numpy.linalg.svd(dotProd)
        except Exception:
            print("Couldn't perform the Single Value Decomposition, skipping alignment", sys.stderr)
        return ref, 0
    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    rmsd = E0 - (2.0 * sum(S))
    rmsd = numpy.sqrt(abs(rmsd / len(ref[refmask])))  # get the rmsd
    # U is simply V*Wt
    U = numpy.dot(V, Wt)  # get the rotation matrix
    # rotate and translate the molecule
    new_coords = numpy.dot((refCenteredCoords), U) + target_centroid  # translate & rotate
    # new_coords=(refCenteredCoords + target_centroid)
    # print U
    if returnRotMat:
        return U, ref_centroid, target_centroid, rmsd
    return new_coords, rmsd


def squared_distance(coordsA, coordsB):
    """Find the squared distance between two 3-tuples"""
    sqrdist = sum((a - b) ** 2 for a, b in zip(coordsA, coordsB))
    return sqrdist


def rmsd(allcoordsA, allcoordsB):
    """Find the RMSD between two lists of 3-tuples"""
    deviation = sum(squared_distance(atomA, atomB) for (atomA, atomB) in zip(allcoordsA, allcoordsB))
    return math.sqrt(deviation / float(len(allcoordsA)))


def mapToCrystal(xtal, pose):
    """Some docking programs might alter the order of the atoms in the output (like Autodock Vina does...)
    this will mess up the rmsd calculation with OpenBabel"""
    query = pybel.ob.CompileMoleculeQuery(xtal.OBMol)
    mapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
    mappingpose = pybel.ob.vvpairUIntUInt()
    exit = mapper.MapUnique(pose.OBMol, mappingpose)
    return mappingpose[0]


def takeCoords(obmol):
    """Take coordinates of an OBMol as a numpy array"""
    return numpy.array([atom.coords for atom in obmol])


def updateCoords(obmol, newcoords):
    "Update OBMol coordinates. newcoords is a numpy array"
    for i, atom in enumerate(obmol):
        atom.OBAtom.SetVector(*newcoords[i])


def prepareAtomString(idlist):
    s = ""
    n = len(idlist)
    for i, id in enumerate(idlist):
        s += "%i" % id
        if (i + 1) == n:
            s += "\n"
        elif (i + 1) % 35 == 0:
            s += ",\n"
        else:
            s += ","
    return s


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

        print("DONE")

    def process_molecule(
        self,
        mol: pybel.Molecule,
        out: pybel.Outputfile,
    ) -> None:
        # find all matches for the smarts query in the input molecule
        # do the product of all matches with the reference matches
        # for each match pair, align the input molecule to the reference molecule
        # keep the best alignment
        # add the tethered atoms property to the input molecule
        # write the molecule to the output file
        
        # smells like itertools spirit
        
        mol.OBMol.DeleteNonPolarHydrogens()
        molMatchAllIds = self.get_matches(mol)
        numMatchs = len(molMatchAllIds)
        numRefMatchs = len(self.ref_matches_align_data)

        if numMatchs == 0:
            return

        # If more than one match, write an output of the same molecule for each match
        # Start a default bestcoord and rmsd for later looping for each pose

        # func 3
        mol_align_data = MolAlignmentData(mol)
        molCoords = mol_align_data.coords()
        for imatch, molMatchIds in enumerate(molMatchAllIds):
            mol_match_align_data = MolAlignmentData(mol, self.match_to_mask(molMatchIds))
            # loop in each matching way for the input molecule
    
            molMatchCoords = mol_match_align_data.coords(use_mask=True)  # ref
            # return MolAlignmentData(molMatchCoords)
            # func 3

            # func 4
            # Loop over the reference matches
            # Align: Get rotation matrix between the two sets of coords
            # Apply rotation to the whole target molecule
            # MolAlignmentData
            
            bestCoordPerMatch = list(itertools.repeat([0] * numMatchs, numRefMatchs))
            bestRMSPerMatch = list(itertools.repeat([math.inf] * numMatchs, numRefMatchs))
            for ir, ref_match_align_data in enumerate(self.ref_matches_align_data):
                refMatchCoord = ref_match_align_data.coords(use_mask=True)
                # aligned_data = MolAlignmentData(molMatchCoords)
                # superpose = Superpose3D(aligned_data, refMatchCoord)
                # refMatchCoord --> target
                rotMat, targetCentroid, refCentroid, rmsd = superpose3D(
                    molMatchCoords, refMatchCoord, returnRotMat=True
                )
                if rmsd < bestRMSPerMatch[ir][imatch]:
                    newcoords = numpy.dot((molCoords - targetCentroid), rotMat) + refCentroid
                    bestRMSPerMatch[ir][imatch] = rmsd
                    bestCoordPerMatch[ir][imatch] = newcoords

            # func 4

            # func 5
            # Finally update molecule coordinates with the best matching coordinates found
            # change molecule coordinates, set TETHERED ATOMS property and save
            for i in range(numMatchs):
                for irefmatch in range(numRefMatchs):
                    bestCoord = bestCoordPerMatch[irefmatch][i]
                    bestRMS = bestRMSPerMatch[irefmatch][i]
                    print("\tBest RMSD reached (match %d, refmatch %d): %s" % (i, irefmatch, bestRMS))
                    molMatchID = molMatchAllIds[i]

                    updateCoords(mol, bestCoord)

                    newData = pybel.ob.OBPairData()
                    newData.SetAttribute("TETHERED ATOMS")
                    newData.SetValue(prepareAtomString(molMatchID))

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
