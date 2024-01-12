import argparse
import logging
import math
import sys

import numpy
from numpy.linalg import LinAlgError
from numpy.typing import ArrayLike
from openbabel import pybel

logger = logging.getLogger("sdrmsd")

Coordinate = tuple[float, float, float]
SingularValueDecomposition = tuple[ArrayLike[float], ArrayLike[float], ArrayLike[float]]
RMSDResult = float | tuple[float, ArrayLike[float]]


def superpose_3D(
    reference: ArrayLike[float],
    target: ArrayLike[float],
    weights: ArrayLike[float] | None = None,
    reference_mask: ArrayLike[bool] | None = None,
    target_mask: ArrayLike[bool] | None = None,
    return_rotation_matrix: bool = False,
) -> tuple[ArrayLike[float], float, ArrayLike[float]] | tuple[ArrayLike[float], float]:
    """superpose_3D performs 3d superposition using a weighted Kabsch algorithm : http://dx.doi.org/10.1107%2FS0567739476001873 & doi: 10.1529/biophysj.105.066654
    definition : superpose3D(reference, target, weights,refmask,target_mask)
    @parameter 1 :  reference - xyz coordinates of the reference structure (the ligand for instance)
    @type 1 :       float64 numpy array (nx3)
    ---
    @parameter 2 :  target - theoretical target positions to which we should move (does not need to be physically relevant.
    @type 2 :       float 64 numpy array (nx3)
    ---
    @parameter 3:   weights - numpy array of atom weights (usuallly between 0 and 1)
    @type 3 :       float 64 numpy array (n)
    @parameter 4:   mask - a numpy boolean mask for designating atoms to include
    Note reference and target positions must have the same dimensions -> n*3 numpy arrays where n is the number of points (or atoms)
    returns:
    - Tuple containing new coordinates and RMSD (default behavior).
      OR
    - Tuple containing new coordinates, RMSD, and rotation matrix (if return_rotation_matrix is True).
    """
    weights = weights or 1.0
    reference_mask = reference_mask or numpy.ones(len(reference), "bool")
    target_mask = target_mask or numpy.ones(len(target), "bool")
    # First get the centroid of both states
    reference_centroid = numpy.mean(reference[reference_mask] * weights, axis=0)
    # Print reference_centroid
    reference_centered_coords = reference - reference_centroid
    # Print reference_centered_coords
    target_centroid = numpy.mean(target[target_mask] * weights, axis=0)
    target_centered_coords = target - target_centroid
    # Print target_centered_coords
    # The following steps come from : http://www.pymolwiki.org/index.php/OptAlign#The_Code and http://en.wikipedia.org/wiki/Kabsch_algorithm
    # Initial residual, see Kabsch.
    E0 = numpy.sum(
        numpy.sum(
            reference_centered_coords[reference_mask] * reference_centered_coords[reference_mask] * weights, axis=0
        ),
        axis=0,
    ) + numpy.sum(
        numpy.sum(target_centered_coords[target_mask] * target_centered_coords[target_mask] * weights, axis=0), axis=0
    )
    reference_tmp = numpy.copy(reference_centered_coords[reference_mask])
    target_tmp = numpy.copy(target_centered_coords[target_mask])
    # print reference_centered_coords[reference_mask]
    # single value decomposition of the dotProduct of both position vectors
    try:
        V, S, Wt = perform_svd(reference_tmp, target_tmp, weights)
    except LinAlgError:
        warning_msg = "Couldn't perform the Single Value Decomposition, skipping alignment"
        logger.warning(warning_msg)
        print(warning_msg, file=sys.stderr)
        return (reference, 0)
    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(numpy.linalg.det(V) * numpy.linalg.det(Wt))
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    rmsd = E0 - (2.0 * sum(S))
    rmsd = numpy.sqrt(abs(rmsd / len(reference[reference_mask])))  # get the rmsd
    # U is simply V*Wt
    U = numpy.dot(V, Wt)  # get the rotation matrix
    # rotate and translate the molecule
    new_coords = numpy.dot((reference_centered_coords), U) + target_centroid  # translate & rotate
    # new_coords=(reference_centered_coords + target_centroid)
    # print U
    return (new_coords, rmsd, U) if return_rotation_matrix else (new_coords, rmsd)


def perform_svd(
    reference_tmp: ArrayLike[float], target_tmp: ArrayLike[float], weights: ArrayLike[float] | int
) -> SingularValueDecomposition | None:
    try:
        dot_product = numpy.dot(numpy.transpose(reference_tmp), target_tmp * weights)
        svd_result = numpy.linalg.svd(dot_product)
    except LinAlgError:
        svd_result = _handle_svd_linalg_error(reference_tmp, target_tmp)
    return svd_result


def _handle_svd_linalg_error(
    reference_tmp: ArrayLike[float], target_tmp: ArrayLike[float]
) -> SingularValueDecomposition | None:
    try:
        dot_product = numpy.dot(numpy.transpose(reference_tmp), target_tmp)
        svd_result = numpy.linalg.svd(dot_product)
        return svd_result
    except LinAlgError:
        raise


def squared_distance(coordinates_A: Coordinate, coordinates_B: Coordinate) -> float:
    """
    Find the squared distance between two 3-tuples.
    """
    return sum((a - b) ** 2 for a, b in zip(coordinates_A, coordinates_B))


def rmsd(all_coordinates_A: list[Coordinate], all_coordinates_B: list[Coordinate]) -> float:
    """
    Find the root mean square deviation between two lists of 3-tuples.
    """
    deviation = sum(squared_distance(atom_A, atom_B) for atom_A, atom_B in zip(all_coordinates_A, all_coordinates_B))
    return math.sqrt(deviation / len(all_coordinates_A))


def map_to_crystal(xtal: pybel.Molecule, pose: pybel.Molecule) -> tuple[int, int]:
    """
    Some docking programs might alter the order of the atoms in the output (like Autodock Vina does...)
    this will mess up the rmsd calculation with OpenBabel.
    """
    query = pybel.ob.CompileMoleculeQuery(xtal.OBMol)
    mapper: pybel.ob.OBIsomorphismMapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
    mapping_pose = pybel.ob.vvpairUIntUInt()
    exit = mapper.MapUnique(pose.OBMol, mapping_pose)
    return mapping_pose[0]


def get_parser() -> argparse.ArgumentParser:
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
    parser.add_argument(
        "reference",
        type=argparse.FileType("r"),
        help="Path to the SDF file with the reference molecule.",
    )
    parser.add_argument(
        "input",
        type=argparse.FileType("r"),
        help="Path to the SDF file with the molecules to be compared to reference.",
    )
    parser.add_argument(
        "-f",
        "--fit",
        action="store_true",
        default=False,
        help="Superpose molecules before RMSD calculation",
    )
    parser.add_argument(
        "-t",
        "--threshold",
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
        default=False,
        metavar="FILE",
        help=(
            "If declared, write an output SDF file with the input molecules with "
            "a new sdfield <RMSD>. If the molecule was fitted, the fitted molecule coordinates will be saved."
        ),
    )
    return parser


def update_coordinates(obmol: pybel.Molecule, new_coordinates: ArrayLike[float]) -> None:
    """
    Update OBMol coordinates. new_coordinates is a numpy array.
    """
    for i, atom in enumerate(obmol):
        atom.OBAtom.SetVector(*new_coordinates[i])


def get_automorphism_rmsd(target: pybel.Molecule, molecule: pybel.Molecule, fit: bool = False) -> RMSDResult:
    """
    Use Automorphism to reorder target coordinates to match reference coordinates atom order
    for correct RMSD comparison. Only the lowest RMSD will be returned.

    Returns:
      If fit=False: bestRMSD	(float)					RMSD of the best matching mapping.
      If fit=True:  (bestRMSD, molecCoordinates)	(float, npy.array)	RMSD of best match and its molecule fitted coordinates.
    """
    mappings = pybel.ob.vvpairUIntUInt()
    bit_vector = pybel.ob.OBBitVec()
    success = pybel.ob.FindAutomorphisms(target.OBMol, mappings)
    target_coordinates = [atom.coords for atom in target]
    index_to_target_coordinates = dict(enumerate(target_coordinates))
    mappose = numpy.array(map_to_crystal(target, molecule))
    sorted_indices = numpy.argsort(mappose[:, 0])
    mappose_result = mappose[sorted_indices][:, 1]
    molecule_coordinates = [atom.coords for atom in molecule]
    pose_coordinates = numpy.array(molecule_coordinates)[mappose_result]
    rmsd_result = math.inf

    # Loop through automorphisms
    for mapping in mappings:
        automorph_coords = [index_to_target_coordinates[i] for i in mapping]
        mapping_rmsd = rmsd(pose_coordinates, automorph_coords)

        # Update result if the current mapping has a lower RMSD
        if mapping_rmsd < rmsd_result:
            rmsd_result = mapping_rmsd

        # Additional fitting if fit=True
        if fit:
            fitted_pose, fitted_rmsd = superpose_3D(numpy.array(automorph_coords), numpy.array(pose_coordinates))

            # Update result if the fitted RMSD is lower
            if fitted_rmsd < rmsd_result:
                rmsd_result = fitted_rmsd

    return (rmsd_result, fitted_pose) if fit else rmsd_result


def save_molecule_with_rmsd(
    output_sdf: pybel.Outputfile, molecule: pybel.Molecule, rmsd: float, population_value: int
) -> None:
    new_data = pybel.ob.OBPairData()
    new_data.SetAttribute("RMSD")
    new_data.SetValue(f"{rmsd:.3f}")

    pop_data = pybel.ob.OBPairData()
    pop_data.SetAttribute("Population")
    pop_data.SetValue(f"{population_value}")
    molecule.OBMol.CloneData(pop_data)

    molecule.OBMol.CloneData(new_data)  # Add new data
    output_sdf.write(molecule)


# ------------------------- MANDANGA AFTER NAME == MAIN STARTS HERE -------------------------


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    reference_sdf = args.reference
    input_sdf = args.input
    fit = args.fit
    threshold = args.threshold
    out = args.out

    # Read crystal pose
    crystal_pose = get_crystal_pose(reference_sdf)
    crystal_atoms = len(crystal_pose.atoms)

    # Find the RMSD between the crystal pose and each docked pose
    docked_poses = pybel.readfile("sdf", input_sdf)

    display_fit_message(fit)

    skipped = []
    molecules_dict = {}  # Save all poses with their dockid
    population = {}  # Poses to be written
    out_dict = {}

    for i, docked_pose in enumerate(docked_poses, start=1):
        atoms_number = process_docked_pose(docked_pose)

        if atoms_number != crystal_atoms:
            skipped.append(i)
            continue

        rmsd_result = calculate_rmsd(crystal_pose, docked_pose, fit)
        handle_pose_matching(out, i, docked_pose, rmsd_result, threshold, molecules_dict, population, out_dict)

    if out:
        output_sdf = pybel.Outputfile("sdf", out, overwrite=True)  # TODO: Check the second argument
        process_and_save_selected_molecules(output_sdf, out_dict, population)

    if skipped:
        print(f"SKIPPED input molecules due to the number of atom mismatch: {skipped}", file=sys.stderr)


def display_fit_message(fit: bool) -> None:
    message = "FIT" if fit else "NOFIT"
    print(f"POSE\tRMSD_{message}")


def get_crystal_pose(reference_sdf: argparse.FileType) -> pybel.Molecule:
    """
    Read crystal pose file and remove hydrogen atoms. Returns crystal pose molecule.
    """
    crystal_pose = next(pybel.readfile("sdf", reference_sdf))
    crystal_pose.removeh()
    return crystal_pose


def process_docked_pose(docked_pose: pybel.Molecule) -> int:
    """
    Remove hydrogen atoms and return the number of atoms in docked pose molecule.
    """
    docked_pose.removeh()
    return len(docked_pose.atoms)


def calculate_rmsd(crystal: pybel.Molecule, docked_pose: pybel.Molecule, fit: bool = False) -> RMSDResult:
    """
    Perform RMSD calculations and update coordinates if required.
    """
    if fit:
        rmsd_result, fitted_result = get_automorphism_rmsd(crystal, docked_pose, fit)
        update_coordinates(docked_pose, fitted_result)
    else:
        rmsd_result = get_automorphism_rmsd(crystal, docked_pose, fit)

    return rmsd_result


def handle_pose_matching(
    out: bool,
    pose_index: int,
    docked_pose: pybel.Molecule,
    rmsd_result: RMSDResult,
    threshold: float | None,
    molecules_dict: dict[int, pybel.Molecule],
    population: dict[int, int],
    out_dict: dict[int, tuple[pybel.Molecule, RMSDResult]],
) -> None:
    """
    Function to handle pose matching and filtering based on 'threshold' parser argument.
    """
    if threshold:
        match_pose, best_match_value = get_best_matching_pose(docked_pose, threshold, molecules_dict)
        if match_pose is not None:
            print_matching_info(pose_index, match_pose, population, best_match_value)
        else:
            save_or_print_info(out, pose_index, docked_pose, rmsd_result, molecules_dict, population, out_dict)
    else:
        save_or_print_info(out, pose_index, docked_pose, rmsd_result, molecules_dict, population, out_dict)


# TODO: Review names and best match rmsd
def get_best_matching_pose(
    docked_pose: pybel.Molecule, threshold: float, molecules_dict: dict[int, pybel.Molecule]
) -> tuple[float | None, float]:
    match_pose = None
    best_match_value = 999999.0

    for did, prevmol in molecules_dict.items():
        tmprmsd = get_automorphism_rmsd(prevmol, docked_pose)

        if tmprmsd < threshold and tmprmsd < best_match_value:
            best_match_value = tmprmsd
            match_pose = did

    return (match_pose, best_match_value)


def print_matching_info(index: int, match_pose: float, population: dict[float, int], best_match_value: float) -> None:
    print(f"Pose {index} matches pose {match_pose} with {best_match_value:.3f} RMSD", file=sys.stderr)
    population[match_pose] += 1


def save_or_print_info(
    out: bool,
    pose_index: int,
    docked_pose: pybel.Molecule,
    rmsd_result: RMSDResult,
    molecules_dict: dict[int, pybel.Molecule],
    population: dict[int, int],
    out_dict: dict[int, tuple[pybel.Molecule, RMSDResult]],
) -> None:
    """
    Function to save and print information based on 'out' parser argument.
    """
    if out:
        out_dict[pose_index] = (docked_pose, rmsd_result)
    print(f"{pose_index}\t{rmsd_result:.2f}")
    molecules_dict[pose_index] = docked_pose
    population[pose_index] = 1


def process_and_save_selected_molecules(
    output_sdf: pybel.Outputfile,
    out_dict: dict[int, tuple[pybel.Molecule, RMSDResult]],
    population: dict[int, int],
):
    for i in sorted(out_dict.keys()):
        molecule, rmsd_result = out_dict[i]
        # Get the number of matches in the thresholding operation
        population_value = population.get(i) or 1
        save_molecule_with_rmsd(output_sdf, molecule, rmsd_result, population_value)


if __name__ == "__main__":
    main()
