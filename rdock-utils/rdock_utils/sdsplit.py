import argparse
import itertools
import logging
import re
import sys
from dataclasses import dataclass
from typing import Generator, TextIO

from .common import inputs_generator, read_molecules_from_all_inputs

logger = logging.getLogger("SDSplit")


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.WARNING)
    config = get_config(argv)
    logging.root.setLevel(config.log_level)
    inputs = inputs_generator(config.infiles)
    batched_molecules = itertools.batched(read_molecules_from_all_inputs(inputs), config.record_size)
    outputs = outputs_generator(config.output_root)
    for molecule_batch, output in zip(batched_molecules, outputs):
        for molecule in molecule_batch:
            molecule.write(output)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Splits SD records into multiple files of equal number of records")
    infile_help = "input file[s] to be processed. if not provided, stdin is used."
    parser.add_argument("infiles", type=str, nargs="*", help=infile_help)
    parser.add_argument(
        "-r",
        "--record-size",
        dest="rec_size",
        default=1000,
        type=int,
        metavar="RecSize",
        help="Record size to split into (default = 1000 records)",
    )
    output_root_help = "Root name for output files (default = tmp)"
    parser.add_argument("-o", default="tmp", type=str, dest="output_root", metavar="OutputRoot", help=output_root_help)
    parser.add_argument("-l", "--log-level", type=str, default="INFO")
    return parser


@dataclass
class SDSplitConfig:
    infiles: list[str]
    record_size: int
    output_root: str
    log_level: str


def get_config(argv: list[str] | None = None) -> SDSplitConfig:
    parser = get_parser()
    argv = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(sanitize_args(argv))
    return SDSplitConfig(
        infiles=args.infiles,
        record_size=args.rec_size,
        output_root=args.output_root,
        log_level=args.log_level,
    )


def outputs_generator(output_root: str) -> Generator[TextIO, None, None]:
    for file_index in itertools.count():
        filename = f"{output_root}{file_index}.sd"
        logger.info(f"Opening {filename}")
        with open(filename, "w") as f:
            yield f


# This function is just for retrocompatibility with the old -<size> argument
def sanitize_args(argv: list[str]) -> list[str]:
    def _replace_invalid_arg(arg: str) -> str:
        if regex.match(arg):
            logger.warning("Record size definition as -<size> is deprecated. Use -r <size> instead.")
            logger.warning(f"Replacing {arg} with -r={arg[1:]}")
            arg = arg.replace("-", "-r=")
        return arg

    regex = re.compile(r"-[0-9]+")
    return [_replace_invalid_arg(arg) for arg in argv]


if __name__ == "__main__":
    main()
