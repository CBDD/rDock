# Standard Library
import argparse
from logging import getLogger

# Local imports
from .common import inputs_generator, read_molecules_from_all_inputs

logger = getLogger("sdfield")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Adding fields to SD files")
    parser.add_argument("fieldname", type=str, help="name of the field to be added")
    parser.add_argument("value", type=str, help="value of the field to be added")
    infile_help = "input file[s] to be processed. if not provided, stdin is used."
    parser.add_argument("infile", type=str, nargs="*", help=infile_help)
    outfile_help = "output file. if not provided, stdout is used."
    parser.add_argument("-o", "--outfile", default=None, type=str, help=outfile_help)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    inputs = inputs_generator(args.infile)
    for molecule in read_molecules_from_all_inputs(inputs):
        molecule.data[args.fieldname] = args.value
        print(repr(molecule))


if __name__ == "__main__":
    main()
