# Standard Library
import argparse
import sys
from logging import getLogger
from typing import Iterable, TextIO

# Local imports
from .parser import FastSDMol

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


def inputs_generator(inputs: list[str]) -> Iterable[TextIO]:
    if not inputs:
        yield sys.stdin
    else:
        for infile in inputs:
            yield open(infile, "r")


def get_output(outfile: str | None) -> TextIO:
    if outfile:
        return open(outfile, "w")
    else:
        return sys.stdout


def read_molecules(file: TextIO) -> Iterable[FastSDMol]:
    while mol := FastSDMol.read(file):
        yield mol


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    inputs = inputs_generator(args.infile)
    for source in inputs:
        for molecule in read_molecules(source):
            molecule.data[args.fieldname] = args.value
            print(repr(molecule))


if __name__ == "__main__":
    main()
