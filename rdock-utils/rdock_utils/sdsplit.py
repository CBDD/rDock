import argparse
import itertools
import logging
from dataclasses import dataclass
from typing import Generator, TextIO

from .common import inputs_generator, read_molecules_from_all_inputs

logger = logging.getLogger("SDSplit")


def main(argv: list[str] | None = None) -> None:
    config = get_config(argv)
    logging.basicConfig(level=config.log_level)
    logger.info(config)
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
    record_size_help = "Record size to split into (default = 1000 records)"
    parser.add_argument("-", dest="rec_size", default=1000, type=int, metavar="RecSize", help=record_size_help)
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
    args = parser.parse_args(argv)
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


if __name__ == "__main__":
    main()
