import argparse
from dataclasses import dataclass


@dataclass
class SDTetherConfig:
    reference_filename: str
    input_filename: str
    output_filename: str
    smarts: str


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="SDTETHER",
        usage="%(prog)s [options] reference.sdf input.sdf output.sdf 'SMARTS'",
        description="sdtether",
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
        "output",
        type=str,
        help="Path to the SDF output file.",
    )
    parser.add_argument(
        "smarts",
        type=str,
        help="SMARTS expression.",
        metavar="SMARTS",
    )
    return parser


def get_config(args: list[str] | None = None) -> SDTetherConfig:
    parser = get_parser()
    parsed_args = parser.parse_args(args)
    return SDTetherConfig(
        reference_filename=parsed_args.reference,
        input_filename=parsed_args.input,
        output_filename=parsed_args.output,
        smarts=parsed_args.smarts,
    )
