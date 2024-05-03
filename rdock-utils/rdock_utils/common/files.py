import sys
from typing import Generator, TextIO


def inputs_generator(inputs: list[str] | None) -> Generator[TextIO, None, None]:
    if not inputs:
        yield sys.stdin
    else:
        for infile in inputs:
            with open(infile, "r") as f:
                yield f
