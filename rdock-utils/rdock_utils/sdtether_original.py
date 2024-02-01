import math
import sys

import numpy
from openbabel import pybel


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
    refCenteredCoords = ref - ref_centroid
    target_centroid = numpy.mean(target[targetmask] * weights, axis=0)
    targetCenteredCoords = target - target_centroid
    
    # the following steps come from : http://www.pymolwiki.org/index.php/OptAlign#The_Code and http://en.wikipedia.org/wiki/Kabsch_algorithm
    # Initial residual, see Kabsch.
    E0 = numpy.sum(
        numpy.sum(refCenteredCoords[refmask] * refCenteredCoords[refmask] * weights, axis=0), axis=0
    ) + numpy.sum(
        numpy.sum(targetCenteredCoords[targetmask] * targetCenteredCoords[targetmask] * weights, axis=0), axis=0
    )
    reftmp = numpy.copy(refCenteredCoords[refmask])
    targettmp = numpy.copy(targetCenteredCoords[targetmask])
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


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        argv = ["sdtether"] + argv
    argv = argv or sys.argv
    if len(argv) != 5:
        sys.exit("USAGE: %s reference.sdf input.sdf output.sdf 'SMARTS'" % argv[0])

    refsdf = argv[1]
    molsdf = argv[2]
    outsdf = argv[3]
    smarts = pybel.Smarts(argv[4])

    # Read reference pose and get atom list matching smarts query
    # if more than 1 match, take the first one
    ref = next(pybel.readfile("sdf", refsdf))
    refMatchIds = smarts.findall(ref)
    numRefMatchs = len(refMatchIds)

    if not numRefMatchs:
        sys.exit("No match found in the reference structure and the SMARTS string given. Please check it.")

    if numRefMatchs > 1:
        print(
            "More than one match in the reference molecule for the SMARTS string given. Will tether each input molecule all possible ways."
        )

    refIndxPerMatch = [numpy.array(rmi) - 1 for rmi in refMatchIds]

    # Take coordinates for the reference matched atoms
    refCoords = takeCoords(ref)
    refMatchCoords = [numpy.take(refCoords, refIndx, axis=0) for refIndx in refIndxPerMatch]

    # Do the same for molecule in molsdf
    out = pybel.Outputfile("sdf", outsdf, overwrite=True)
    molSupp = pybel.readfile("sdf", molsdf)
    # ff = pybel.ob.OBForceField.FindForceField("MMFF94")
    for i, mol in enumerate(molSupp):
        print(f"## Molecule {i+1}")
        mol.OBMol.DeleteNonPolarHydrogens()
        molMatchAllIds = smarts.findall(mol)
        numMatchs = len(molMatchAllIds)

        if numMatchs == 0:
            print("No_Match")
            continue
        elif numMatchs == 1:
            print("Match")
        elif numMatchs > 1:
            print(f"Multiple_Match SMART Matches for this molecule ({numMatchs})")

        # If more than one match, write an output of the same molecule for each match
        # Start a default bestcoord and rmsd for later looping for each pose
        bestCoordPerMatch = [[0 for i in range(numMatchs)] for i in range(numRefMatchs)]
        bestRMSPerMatch = [[999 for i in range(numMatchs)] for i in range(numRefMatchs)]

        # Will do a randomrotorsearch to find conformer with the lower rmsd when superposing
        # At least 20 when possible
        # ff.Setup(mol.OBMol)
        # numats = mol.OBMol.NumAtoms()
        # numrot = mol.OBMol.NumRotors()
        # print "Atoms: %i, Rotors: %i"%(numats, numrot)
        # geomopt = 300
        # genconf = 100
        # increase iterations if bigger molecule or bigger number of rotatable bonds
        # for allowing better sampling
        # if numats > 40 and numrot > 5:
        # geomopt = 300
        # genconf = 150
        # if numats > 55 and numrot > 7:
        # genconf = 100
        # geomopt = 500
        # print "\tDoing conformational search with WeightedRotorSearch (%i, %i)..."%(genconf, geomopt),
        # ff.SteepestDescent(500, 1.0e-4)
        # ff.WeightedRotorSearch(genconf,geomopt)
        # ff.ConjugateGradients(500, 1.0e-6)
        # ff.GetConformers(mol.OBMol)
        # numconf = mol.OBMol.NumConformers()
        numconf = 1
        # print "%i conformers generated"%numconf
        if numconf > 1:
            # Doing conf search
            # for i in range(numconf):
            #    mol.OBMol.SetConformer(i)
            #    confCoords = takeCoords(mol)
            # 	print 'coord:',confCoords[0,:]
            #
            #    for imatch, molMatchIds in enumerate(molMatchAllIds):
            #        molMatchIndx = numpy.array(molMatchIds) - 1
            #        confMatchCoords = numpy.take(confCoords, molMatchIndx, axis=0)
            #
            #        # Align: Get rotation matrix between the two sets of coords
            #        # Apply rotation to the whole target molecule
            #        rotMat, targetCentroid, refCentroid, rmsd = superpose3D(confMatchCoords, refMatchCoords, returnRotMat=True)
            #        if rmsd < bestRMSPerMatch[imatch]:
            #            newcoords = numpy.dot((confCoords - targetCentroid), rotMat) + refCentroid
            #            bestRMSPerMatch[imatch] = rmsd
            #            bestCoordPerMatch[imatch] = newcoords
            #   if bestrms < 0.01: break
            pass
        else:
            molCoords = takeCoords(mol)
            for imatch, molMatchIds in enumerate(molMatchAllIds):
                # loop in each matching way for the input molecule
                molMatchIndx = numpy.array(molMatchIds) - 1
                molMatchCoords = numpy.take(molCoords, molMatchIndx, axis=0)

                # Loop over the reference matches
                # Align: Get rotation matrix between the two sets of coords
                # Apply rotation to the whole target molecule
                for ir, refMatchCoord in enumerate(refMatchCoords):
                    rotMat, targetCentroid, refCentroid, rmsd = superpose3D(
                        molMatchCoords, refMatchCoord, returnRotMat=True
                    )
                    if rmsd < bestRMSPerMatch[ir][imatch]:
                        newcoords = numpy.dot((molCoords - targetCentroid), rotMat) + refCentroid
                        bestRMSPerMatch[ir][imatch] = rmsd
                        bestCoordPerMatch[ir][imatch] = newcoords

            # Finally update molecule coordinates with the best matching coordinates found
            # change molecule coordinates, set TETHERED ATOMS property and save
            for imatch in range(numMatchs):
                for irefmatch in range(numRefMatchs):
                    bestCoord = bestCoordPerMatch[irefmatch][imatch]
                    bestRMS = bestRMSPerMatch[irefmatch][imatch]
                    print(f"\tBest RMSD reached (match {imatch}, refmatch {irefmatch}): {bestRMS:.4f}")
                    molMatchID = molMatchAllIds[imatch]
                    updateCoords(mol, bestCoord)
                    newData = pybel.ob.OBPairData()
                    newData.SetAttribute("TETHERED ATOMS")
                    newData.SetValue(prepareAtomString(molMatchID))
    
                    mol.OBMol.DeleteData("TETHERED ATOMS")  # Remove Previous DATA
                    mol.OBMol.CloneData(newData)  # Add new data
                    out.write(mol)

    out.close()

    print("DONE")


if __name__ == "__main__":
    main()
