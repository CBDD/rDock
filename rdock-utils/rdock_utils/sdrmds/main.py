import logging

import numpy
from openbabel import pybel

from .parser import get_parser
from .sdrmsd import SDRMSD, PoseMatchData, SDRMSDData

logger = logging.getLogger("SDRMSD")

Coordinate = tuple[float, float, float]
SingularValueDecomposition = tuple[numpy.ndarray[float], numpy.ndarray[float], numpy.ndarray[float]]
RMSDResult = float | tuple[float, numpy.ndarray[float]]


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    sdrmsd = SDRMSD(args.reference, args.input, args.fit, args.threshold, args.out)

    # Find the RMSD between the crystal pose and each docked pose
    docked_poses = pybel.readfile("sdf", sdrmsd.input_sdf)

    display_fit_message(sdrmsd.fit)

    data = SDRMSDData()

    # Read crystal pose
    crystal_pose = get_crystal_pose(sdrmsd.reference_sdf)
    crystal_atoms = len(crystal_pose.atoms)
    for i, docked_pose in enumerate(docked_poses, start=1):
        atoms_number = process_docked_pose(docked_pose)

        if atoms_number != crystal_atoms:
            data.skipped.append(i)
            continue

        rmsd_result = calculate_rmsd(crystal_pose, docked_pose, sdrmsd)
        pose_match_data = PoseMatchData(i, docked_pose, data)
        handle_pose_matching(rmsd_result, pose_match_data, sdrmsd)

    if sdrmsd.out:
        output_sdf = pybel.Outputfile("sdf", sdrmsd.out, overwrite=True)
        process_and_save_selected_molecules(output_sdf, data)

    if data.skipped:
        logger.warning(f"SKIPPED input molecules due to the number of atom mismatch: {data.skipped}")


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


def calculate_rmsd(crystal: pybel.Molecule, docked_pose: pybel.Molecule, sdrmsd: SDRMSD) -> RMSDResult:
    """
    Perform RMSD calculations and update coordinates if required.
    """
    if sdrmsd.fit:
        rmsd_result, fitted_result = sdrmsd.get_automorphism_rmsd(crystal, docked_pose, sdrmsd)
        sdrmsd.update_coordinates(docked_pose, fitted_result)
    else:
        rmsd_result = sdrmsd.get_automorphism_rmsd(crystal, docked_pose, sdrmsd)

    return rmsd_result


def handle_pose_matching(rmsd_result: RMSDResult, pose_match_data: PoseMatchData, sdrmsd: SDRMSD) -> None:
    """
    Function to handle pose matching and filtering based on 'threshold' parser argument.
    """
    threshold = sdrmsd.threshold
    out = sdrmsd.out

    if threshold:
        match_pose, best_match_value = sdrmsd.get_best_matching_pose(pose_match_data, threshold)
        if match_pose is not None:
            print_matching_info(pose_match_data, match_pose, best_match_value)
        else:
            save_and_print_info(rmsd_result, pose_match_data, out)
    else:
        save_and_print_info(rmsd_result, pose_match_data, out)


def process_and_save_selected_molecules(output_sdf: pybel.Outputfile, data: SDRMSDData):
    for i in sorted(data.out_dict.keys()):
        molecule, rmsd_result = data.out_dict[i]
        # Get the number of matches in the thresholding operation
        population_value = data.population.get(i, 1)
        save_molecule_with_rmsd(output_sdf, molecule, rmsd_result, population_value)


def print_matching_info(pose_match_data: PoseMatchData, match_pose: int, best_match_value: float) -> None:
    index = pose_match_data.pose_index
    population = pose_match_data.sdrmsd_data.population
    logger.info(f"Pose {index} matches pose {match_pose} with {best_match_value:.3f} RMSD")
    population[match_pose] += 1


def save_and_print_info(rmsd_result: RMSDResult, pose_match_data: PoseMatchData, out: bool) -> None:
    """
    Function to save and print information based on 'out' parser argument.
    """
    index = pose_match_data.pose_index
    docked_pose = pose_match_data.docked_pose
    out_dict = pose_match_data.sdrmsd_data.out_dict
    molecules_dict = pose_match_data.sdrmsd_data.molecules_dict
    population = pose_match_data.sdrmsd_data.population

    if out:
        out_dict[index] = (docked_pose, rmsd_result)
    print(f"{index}\t{rmsd_result:.2f}")
    molecules_dict[index] = docked_pose
    population[index] = 1


def save_molecule_with_rmsd(
    output_sdf: pybel.Outputfile, molecule: pybel.Molecule, rmsd: float, population_value: int
) -> None:
    new_data = pybel.ob.OBPairData()
    new_data.SetAttribute("RMSD")
    new_data.SetValue(f"{rmsd:.3f}")
    # pasar al helper
    # gestionar que el population sea opcional no booleano
    # if threshold
    pop_data = pybel.ob.OBPairData()
    pop_data.SetAttribute("Population")
    pop_data.SetValue(f"{population_value}")
    molecule.OBMol.CloneData(pop_data)

    molecule.OBMol.CloneData(new_data)  # Add new data
    output_sdf.write(molecule)


if __name__ == "__main__":
    main()
