# Standard Library
import argparse
import sys
from dataclasses import dataclass
from logging import getLogger
from typing import Generator, TextIO

# Local imports
from .common import inputs_generator, read_molecules_from_all_inputs

logger = getLogger("sdfield")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Adding fields to SD files")
    parser.add_argument("fieldname", type=str, help="name of the field to be added")
    parser.add_argument("value", type=str, help="value of the field to be added")
    infiles_help = "input file[s] to be processed. if not provided, stdin is used."
    parser.add_argument("infiles", type=str, nargs="*", help=infiles_help)
    outfile_help = "output file. if not provided, stdout is used."
    parser.add_argument("-o", "--outfile", default=None, type=str, help=outfile_help)

    return parser


@dataclass
class SDFieldConfig:
    fieldname: str
    value: str
    infiles: list[str]
    outfile: str | None

    def get_outfile(self) -> Generator[TextIO, None, None]:
        if self.outfile:
            with open(self.outfile, "w") as f:
                yield f
        else:
            yield sys.stdout


def get_config(argv: list[str] | None = None) -> SDFieldConfig:
    parser = get_parser()
    args = parser.parse_args(argv)
    return SDFieldConfig(fieldname=args.fieldname, value=args.value, infiles=args.infiles, outfile=args.outfile)


def main(argv: list[str] | None = None) -> None:
    config = get_config(argv)
    inputs = inputs_generator(config.infiles)
    output_gen = config.get_outfile()
    output = next(output_gen)
    for molecule in read_molecules_from_all_inputs(inputs):
        molecule.data[config.fieldname] = config.value
        molecule.write(output)


if __name__ == "__main__":
    main()
