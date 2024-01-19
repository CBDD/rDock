import logging

import numpy
from openbabel import pybel

from .parser import get_parser
from .sdrmsd import SDRMSD, PoseMatchData, SDRMSDData

logger = logging.getLogger("SDRMSD")

Coordinate = tuple[float, float, float]
SingularValueDecomposition = tuple[numpy.ndarray[float], numpy.ndarray[float], numpy.ndarray[float]]
AutomorphismRMSD = float | tuple[float, numpy.ndarray[float]]


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    sdrmsd = SDRMSD(args.reference, args.input, args.fit, args.threshold, args.out)

    # Find the RMSD between the crystal pose and each docked pose
    docked_poses = pybel.readfile("sdf", sdrmsd.input_sdf)

    sdrmsd.display_fit_message()

    data = SDRMSDData()

    # Read crystal pose
    crystal_pose = sdrmsd.get_crystal_pose()
    crystal_atoms = len(crystal_pose.atoms)
    for i, docked_pose in enumerate(docked_poses, start=1):
        atoms_number = sdrmsd.process_docked_pose(docked_pose)

        if atoms_number != crystal_atoms:
            data.skipped.append(i)
            continue

        rmsd_result = sdrmsd.calculate_rmsd(crystal_pose, docked_pose)
        pose_match_data = PoseMatchData(i, docked_pose, data)
        sdrmsd.handle_pose_matching(rmsd_result, pose_match_data)

    if sdrmsd.out:
        output_sdf = pybel.Outputfile("sdf", sdrmsd.out, overwrite=True)
        sdrmsd.process_and_save_selected_molecules(output_sdf, data)

    if data.skipped:
        logger.warning(f"SKIPPED input molecules due to the number of atom mismatch: {data.skipped}")


if __name__ == "__main__":
    main()
