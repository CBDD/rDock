import argparse
import math
import os
import sys
from typing import TextIO

import numpy
from numpy.typing import ArrayLike
from openbabel import pybel

Coordinate = tuple[float, float, float]


def superpose_3D(
    ref: ArrayLike[float],
    target: ArrayLike[float],
    weights: ArrayLike[float] | None = None,
    refmask: ArrayLike[bool] | None = None,
    targetmask: ArrayLike[bool] | None = None,
    returnRotMat: bool = False,
) -> tuple[ArrayLike[float], float, ArrayLike[float]] | tuple[ArrayLike[float], float]:
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
            print("Couldn't perform the Single Value Decomposition, skipping alignment", file=sys.stderr)
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
        return new_coords, rmsd, U
    return new_coords, rmsd


def squared_distance(coordsA: Coordinate, coordsB: Coordinate) -> float:
    """Find the squared distance between two 3-tuples"""
    sqrdist = sum((a - b) ** 2 for a, b in zip(coordsA, coordsB))
    return sqrdist


def rmsd(allcoordsA: list[Coordinate], allcoordsB: list[Coordinate]) -> float:
    """Find the RMSD between two lists of 3-tuples"""
    deviation = sum(squared_distance(atomA, atomB) for (atomA, atomB) in zip(allcoordsA, allcoordsB))
    return math.sqrt(deviation / float(len(allcoordsA)))


def map_to_crystal(xtal: pybel.Molecule, pose: pybel.Molecule) -> list[tuple[int, int]]:
    """
    Some docking programs might alter the order of the atoms in the output (like Autodock Vina does...)
    this will mess up the rmsd calculation with OpenBabel
    """
    query = pybel.ob.CompileMoleculeQuery(xtal.OBMol)
    mapper: pybel.ob.OBIsomorphismMapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
    mappingpose = pybel.ob.vvpairUIntUInt()
    exit = mapper.MapUnique(pose.OBMol, mappingpose)
    return mappingpose[0]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="SDRMSD",
        usage="%(prog)s [options] reference.sdf input.sdf",
        description="Superpose molecules before RMSD calculation",
        epilog=(
            "Arguments:\n"
            "   reference.sdf   SDF file with the reference molecule.\n"
            "   input.sdf       SDF file with the molecules to be compared to reference.\n"
        ),
    )
    #TODO: add the two required arguments: reference.sdf and input.sdf
    parser.add_argument(
        "-f",
        "--fit",
        dest="fit",
        action="store_true",
        default=False,
        help="Superpose molecules before RMSD calculation",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        dest="threshold",
        action="store",
        type=float,
        help=(
            "Discard poses with RMSD < THRESHOLD with respect previous poses "
            "which were not rejected based on the same principle. A Population "
            "SDField will be added to output SD with the population number."
        ),
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="outfilename",
        metavar="FILE",
        default=False,
        help=(
            "If declared, write an output SDF file with the input molecules with "
            "a new sdfield <RMSD>. If the molecule was fitted, the fitted molecule coordinates will be saved."
        ),
    )

    # Check we have two positional arguments
    if len(sys.argv) < 3:
        parser.error("Incorrect number of positional arguments. Use -h or --help options to print help.")

    return parser.parse_args()


def update_coordinates(obmol: pybel.Molecule, new_coordinates: ArrayLike[float]) -> None:
    "Update OBMol coordinates. new_coordinates is a numpy array"
    for i, atom in enumerate(obmol):
        atom.OBAtom.SetVector(*new_coordinates[i])


def get_automorphism_RMSD(
    target: pybel.Molecule, molec: pybel.Molecule, fit: bool = False
) -> float | tuple[float, ArrayLike[float]]:
    """
    Use Automorphism to reorder target coordinates to match ref coordinates atom order
    for correct RMSD comparison. Only the lowest RMSD will be returned.

    Returns:
      If fit=False: 	bestRMSD	(float)					RMSD of the best matching mapping.
      If fit=True:	(bestRMSD, molecCoordinates)	(float, npy.array)	RMSD of best match and its molecule fitted coordinates.
    """
    mappings = pybel.ob.vvpairUIntUInt()
    bitvec = pybel.ob.OBBitVec()
    lookup = []
    for i, atom in enumerate(target):
        lookup.append(i)
    success = pybel.ob.FindAutomorphisms(target.OBMol, mappings)
    targetcoords = [atom.coords for atom in target]
    mappose = numpy.array(map_to_crystal(target, molec))
    mappose = mappose[numpy.argsort(mappose[:, 0])][:, 1]
    posecoords = numpy.array([atom.coords for atom in molec])[mappose]
    resultrmsd = 999999999999
    for mapping in mappings:
        automorph_coords = [None] * len(targetcoords)
        for x, y in mapping:
            automorph_coords[lookup.index(x)] = targetcoords[lookup.index(y)]
        mapping_rmsd = rmsd(posecoords, automorph_coords)
        if mapping_rmsd < resultrmsd:
            resultrmsd = mapping_rmsd
            fitted_result = False
        if fit:
            fitted_pose, fitted_rmsd = superpose_3D(numpy.array(automorph_coords), numpy.array(posecoords))
            if fitted_rmsd < resultrmsd:
                resultrmsd = fitted_rmsd
                fitted_result = fitted_pose

    if fit:
        return (resultrmsd, fitted_pose)
    else:
        return resultrmsd


def save_molecule_with_RMSD(outsdf: TextIO, molec: pybel.Molecule, rmsd: float, population: bool = False) -> None:
    newData = pybel.ob.OBPairData()
    newData.SetAttribute("RMSD")
    newData.SetValue("%.3f" % rmsd)

    if population:
        popData = pybel.ob.OBPairData()
        popData.SetAttribute("Population")
        popData.SetValue("%i" % population)
        molec.OBMol.CloneData(popData)
    molec.OBMol.CloneData(newData)  # Add new data
    outsdf.write(molec)


if __name__ == "__main__":

    def main() -> None:
        args = parse_arguments()
        xtal = sys.argv[1]
        poses = sys.argv[2]

        if not os.path.exists(xtal) or not os.path.exists(poses):
            sys.exit("Input files not found. Please check the path given is correct.")
        
        

    (opts, args) = parse_arguments()

    xtal = args[0]
    poses = args[1]

    if not os.path.exists(xtal) or not os.path.exists(poses):
        sys.exit("Input files not found. Please check the path given is correct.")

    fit = opts.fit
    outfname = opts.outfilename
    threshold = opts.threshold

    # Read crystal pose
    crystal = next(pybel.readfile("sdf", xtal))
    crystal.removeh()
    crystalnumatoms = len(crystal.atoms)

    # If outfname is defined, prepare an output SDF sink to write molecules
    if outfname:
        outsdf = pybel.Outputfile("sdf", outfname, overwrite=True)

    # Find the RMSD between the crystal pose and each docked pose
    dockedposes = pybel.readfile("sdf", poses)
    if fit:
        print("POSE\tRMSD_FIT")
    else:
        print("POSE\tRMSD_NOFIT")

    skipped = []
    moleclist = {}  # Save all poses with their dockid
    population = {}  # Poses to be written
    outlist = {}
    for docki, dockedpose in enumerate(dockedposes):
        dockedpose.removeh()
        natoms = len(dockedpose.atoms)
        if natoms != crystalnumatoms:
            skipped.append(docki + 1)
            continue
        if fit:
            resultrmsd, fitted_result = get_automorphism_RMSD(crystal, dockedpose, fit=True)
            update_coordinates(dockedpose, fitted_result)
        else:
            resultrmsd = get_automorphism_RMSD(crystal, dockedpose, fit=False)

        if threshold:
            # Calculate RMSD between all previous poses
            # Discard if rmsd < FILTER threshold
            if moleclist:
                match = None
                bestmatchrmsd = 999999
                for did, prevmol in moleclist.items():
                    tmprmsd = get_automorphism_RMSD(prevmol, dockedpose)
                    if tmprmsd < threshold:
                        if tmprmsd < bestmatchrmsd:
                            bestmatchrmsd = tmprmsd
                            match = did

                if match is not None:
                    # Do not write this one
                    # sum one up to the matching previous molecule id
                    print(f"Pose {docki + 1} matches pose {match + 1} with {bestmatchrmsd:.3f} RMSD", file=sys.stderr)
                    population[match] += 1
                else:
                    # There's no match. Print info for this one and write to outsdf if needed
                    # Save this one!
                    if outfname:
                        outlist[docki] = (dockedpose, resultrmsd)
                    print(f"{docki + 1}\t{resultrmsd:.2f}")
                    moleclist[docki] = dockedpose
                    population[docki] = 1
            else:
                # First molecule in list. Append for sure
                moleclist[docki] = dockedpose
                population[docki] = 1
                if outfname:
                    outlist[docki] = (dockedpose, resultrmsd)
        else:
            # Just write the best rmsd found and the molecule to outsdf if demanded
            if outfname:
                save_molecule_with_RMSD(outsdf, dockedpose, resultrmsd)
            print(f"{docki + 1}\t{resultrmsd:.2f}")

    if outlist:
        # Threshold applied and outlist needs to be written
        for docki in sorted(outlist.keys()):
            molrmsd = outlist[docki]
            # Get the number of matches in the thresholding operation
            pop = population.get(docki)
            if not pop:
                pop = 1
            # Save molecule
            save_molecule_with_RMSD(outsdf, molrmsd[0], molrmsd[1], pop)

    if skipped:
        print(f"SKIPPED input molecules due to the number of atom mismatch: {skipped}", file=sys.stderr)
