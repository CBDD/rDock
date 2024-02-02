# Standard Library
import logging
import sys
from typing import Iterable, TextIO

# Local imports
from .parser import FastSDMol

logger = logging.getLogger("rdock-utils:common")


def inputs_generator(inputs: list[str]) -> Iterable[TextIO]:
    if not inputs:
        yield sys.stdin
    else:
        for infile in inputs:
            with open(infile, "r") as f:
                yield f


def read_molecules(file: TextIO) -> Iterable[FastSDMol]:
    while True:
        try:
            mol = FastSDMol.read(file)
            if mol is None:
                break
            yield mol
        except ValueError as e:
            logger.warning(f"error reading molecule: {e}")
