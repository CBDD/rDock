import argparse
import sys
from dataclasses import dataclass

from rdock_utils.common import inputs_generator, read_molecules_from_all_inputs


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Set the first title line equal to a given data field")
    parser.add_argument(
        "-f", "--field", type=str, required=True, help="Data field to set the first title line equal to"
    )
    infile_help = "input file[s] to be processed. if not provided, stdin is used."
    parser.add_argument("infiles", type=str, nargs="*", help=infile_help)
    return parser


@dataclass
class SDModifyConfig:
    field: str
    infiles: list[str]


def get_config(argv: list[str] | None = None) -> SDModifyConfig:
    parser = get_parser()
    args = parser.parse_args(argv)
    return SDModifyConfig(field=args.field, infiles=args.infiles)


def main(argv: list[str] | None = None) -> None:
    config = get_config(argv)
    inputs = inputs_generator(config.infiles)
    for index, mol in enumerate(read_molecules_from_all_inputs(inputs), start=1):
        value = mol.get_field(config.field) if config.field != "_REC" else str(index)
        mol.set_title(value)
        mol.write(sys.stdout)
