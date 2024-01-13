import argparse
import logging
import math
import sys

import numpy
from openbabel import pybel

from rdock_utils.helper import Helper, HelperData, PoseMatchData

logger = logging.getLogger("sdrmsd")

Coordinate = tuple[float, float, float]
SingularValueDecomposition = tuple[numpy.ndarray[float], numpy.ndarray[float], numpy.ndarray[float]]
RMSDResult = float | tuple[float, numpy.ndarray[float]]


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    helper = Helper(args.reference, args.input, args.fit, args.threshold, args.out)

    # Read crystal pose
    crystal_pose = get_crystal_pose(helper.reference_sdf)
    crystal_atoms = len(crystal_pose.atoms)

    # Find the RMSD between the crystal pose and each docked pose
    docked_poses = pybel.readfile("sdf", helper.input_sdf)

    display_fit_message(helper.fit)

    data = HelperData()

    for i, docked_pose in enumerate(docked_poses, start=1):
        atoms_number = process_docked_pose(docked_pose)

        if atoms_number != crystal_atoms:
            data.skipped.append(i)
            continue

        rmsd_result = calculate_rmsd(crystal_pose, docked_pose, helper)
        pose_match_data = PoseMatchData(i, docked_pose, data)
        handle_pose_matching(rmsd_result, pose_match_data, helper.out)

    if helper.out:
        output_sdf = pybel.Outputfile("sdf", helper.out, overwrite=True)
        process_and_save_selected_molecules(output_sdf, data)

    if data.skipped:
        print(f"SKIPPED input molecules due to the number of atom mismatch: {data.skipped}", file=sys.stderr)


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
        type=str,
        help="Path to the SDF file with the reference molecule.",
    )
    parser.add_argument(
        "input",
        type=str,
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


def get_crystal_pose(reference_sdf: str) -> pybel.Molecule:
    """
    Read crystal pose file and remove hydrogen atoms. Returns crystal pose molecule.
    """
    crystal_pose = next(pybel.readfile("sdf", reference_sdf))
    crystal_pose.removeh()
    return crystal_pose


def display_fit_message(fit: bool) -> None:
    message = "FIT" if fit else "NOFIT"
    print(f"POSE\tRMSD_{message}")


def process_docked_pose(docked_pose: pybel.Molecule) -> int:
    """
    Remove hydrogen atoms and return the number of atoms in docked pose molecule.
    """
    docked_pose.removeh()
    return len(docked_pose.atoms)


def calculate_rmsd(crystal: pybel.Molecule, docked_pose: pybel.Molecule, helper: Helper) -> RMSDResult:
    """
    Perform RMSD calculations and update coordinates if required.
    """
    if helper.fit:
        rmsd_result, fitted_result = get_automorphism_rmsd(crystal, docked_pose, helper)
        helper.update_coordinates(docked_pose, fitted_result)
    else:
        rmsd_result = get_automorphism_rmsd(crystal, docked_pose, helper)

    return rmsd_result


def handle_pose_matching(rmsd_result: RMSDResult, pose_match_data: PoseMatchData, helper: Helper) -> None:
    """
    Function to handle pose matching and filtering based on 'threshold' parser argument.
    """
    threshold = helper.threshold
    out = helper.out

    if threshold:
        match_pose, best_match_value = get_best_matching_pose(pose_match_data, threshold)
        if match_pose is not None:
            print_matching_info(pose_match_data, match_pose, best_match_value)
        else:
            save_or_print_info(rmsd_result, pose_match_data, out)
    else:
        save_or_print_info(rmsd_result, pose_match_data, out)


def process_and_save_selected_molecules(output_sdf: pybel.Outputfile, data: HelperData):
    for i in sorted(data.out_dict.keys()):
        molecule, rmsd_result = data.out_dict[i]
        # Get the number of matches in the thresholding operation
        population_value = data.population.get(i) or 1
        save_molecule_with_rmsd(output_sdf, molecule, rmsd_result, population_value)


# TODO: Review names and best match rmsd
def get_best_matching_pose(pose_match_data: PoseMatchData, threshold: float) -> tuple[float | None, float]:
    docked_pose = pose_match_data.docked_pose
    molecules_dict = pose_match_data.helper_data.molecules_dict
    match_pose = None
    best_match_value = 999999.0

    for did, prevmol in molecules_dict.items():
        tmprmsd = get_automorphism_rmsd(prevmol, docked_pose)

        if tmprmsd < threshold and tmprmsd < best_match_value:
            best_match_value = tmprmsd
            match_pose = did

    return (match_pose, best_match_value)


def print_matching_info(pose_match_data: PoseMatchData, match_pose: float, best_match_value: float) -> None:
    index = pose_match_data.pose_index
    population = pose_match_data.helper_data.population

    print(f"Pose {index} matches pose {match_pose} with {best_match_value:.3f} RMSD", file=sys.stderr)
    population[match_pose] += 1


def save_or_print_info(rmsd_result: RMSDResult, pose_match_data: PoseMatchData, out: bool) -> None:
    """
    Function to save and print information based on 'out' parser argument.
    """
    index = pose_match_data.pose_index
    docked_pose = pose_match_data.docked_pose
    out_dict = pose_match_data.helper_data.out_dict
    molecules_dict = pose_match_data.helper_data.molecules_dict
    population = pose_match_data.helper_data.population

    if out:
        out_dict[index] = (docked_pose, rmsd_result)
    print(f"{index}\t{rmsd_result:.2f}")
    molecules_dict[index] = docked_pose
    population[index] = 1


def get_automorphism_rmsd(target: pybel.Molecule, molecule: pybel.Molecule, helper: Helper) -> RMSDResult:
    """
    Use Automorphism to reorder target coordinates to match reference coordinates atom order
    for correct RMSD comparison. Only the lowest RMSD will be returned.

    Returns:
      If fit=False: bestRMSD	(float)					RMSD of the best matching mapping.
      If fit=True:  (bestRMSD, molecCoordinates)	(float, npy.array)	RMSD of best match and its molecule fitted coordinates.
    """
    fit = helper.fit
    mappings = pybel.ob.vvpairUIntUInt()
    bit_vector = pybel.ob.OBBitVec()
    success = pybel.ob.FindAutomorphisms(target.OBMol, mappings)
    target_coordinates = [atom.coords for atom in target]
    index_to_target_coordinates = dict(enumerate(target_coordinates))
    mappose = numpy.array(helper.map_to_crystal(target, molecule))
    sorted_indices = numpy.argsort(mappose[:, 0])
    mappose_result = mappose[sorted_indices][:, 1]
    molecule_coordinates = [atom.coords for atom in molecule]
    pose_coordinates = numpy.array(molecule_coordinates)[mappose_result]
    rmsd_result = math.inf

    # Loop through automorphisms
    for mapping in mappings:
        automorph_coords = [index_to_target_coordinates[i] for i in mapping]
        mapping_rmsd = helper.rmsd(pose_coordinates, automorph_coords)

        # Update result if the current mapping has a lower RMSD
        if mapping_rmsd < rmsd_result:
            rmsd_result = mapping_rmsd

        # Additional fitting if fit=True
        if fit:
            fitted_pose, fitted_rmsd = helper.superpose_3D(
                numpy.array(automorph_coords), numpy.array(pose_coordinates)
            )

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


if __name__ == "__main__":
    main()
