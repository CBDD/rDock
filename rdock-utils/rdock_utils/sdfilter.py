# Standard Library
import argparse
import logging
import operator
from pathlib import Path
import shlex
from typing import Any, Callable

# Local imports
from .common import inputs_generator, read_molecules
from .parser import FastSDMol

logger = logging.getLogger("sdfilter")

Filter = Callable[[FastSDMol, int, int], bool]


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filters SD records by data fields")
    filter_help = (
        "Filters can be provided as a string or in a file, one per line. All filters are OR'd together.\n"
        "Filters follow the format:\n"
        "'$<DataField> <Operator> <Value>'\n"
        "where valid operators are: '==', '!=', '<', '>', '<=', and '>=' for general values,\n"
        "'in' and 'not_in' for strings, and 'eq', 'ne', 'lt', 'gt', 'le', and 'ge' \n"
        "for strings for perl version retro-compatibility.\n"
        "_REC (record number), _TITLE1, _TITLE2, and _TITLE3 are provided as a pseudo-data field\n"
        "rdock-utils provides expanded functionality, where two data fields can be compared\n"
        "using the following syntax:\n"
        "'$<DataField1> <Operator> $<DataField2>'\n"
        "also, any combination of literal filters and filter files can be provided\n"
        "filter files including other filters are supported as well, so be careful with recursion\n"
    )
    parser.add_argument("-f", "--filter", type=str, nargs=1, action="extend", help=filter_help, required=True)
    s_help = "If -s option is used, _COUNT (#occurrences of DataField) is provided as a pseudo-data field"
    parser.add_argument("-s", type=str, default=None, help=s_help)
    infile_help = "input file[s] to be processed. if not provided, stdin is used."
    parser.add_argument("infile", type=str, nargs="*", help=infile_help)
    return parser


def get_operand_raw_value(operand: str, molecule: FastSDMol, record: int, field_count: int) -> str | int:
    if operand == "$_REC":
        return record
    elif operand == "$_COUNT":
        return field_count
    elif operand.startswith("$"):
        op = molecule.get_field(operand[1:])
        if op is None:
            logger.warning(f"field {operand} not found in record {record}, assuming literal string")
            return operand
        return op
    else:
        return operand


def get_casted_operands(operand1: str | int, operand2: str | int) -> tuple[str, str] | tuple[float, float]:
    try:
        return (float(operand1), float(operand2))
    except ValueError:
        return (str(operand1), str(operand2))


OPERATION_MAP: dict[str, Callable[[Any, Any], bool]] = {
    "==": operator.eq,
    "!=": operator.ne,
    "<": operator.lt,
    ">": operator.gt,
    "<=": operator.le,
    ">=": operator.ge,
    "eq": operator.eq,
    "ne": operator.ne,
    "lt": operator.lt,
    "gt": operator.gt,
    "le": operator.le,
    "ge": operator.ge,
    "in": lambda x, y: x in y,
    "not_in": lambda x, y: x not in y,
}


def create_filter(filter_str: str) -> list[Filter]:
    tokens = shlex.split(filter_str)
    if len(tokens) == 1 and (path := Path(filter_str)).is_file():
        with open(path, "r") as f:
            return [filter for line in f for filter in create_filter(line.strip())]
    elif len(tokens) != 3:
        raise ValueError(f"invalid filter: {filter_str}")

    if tokens[1] not in OPERATION_MAP:
        raise ValueError(f"invalid operator: {tokens[1]}. expected: {OPERATION_MAP.keys()}")

    def func(molecule: FastSDMol, record: int, field_count: int) -> bool:
        raw_operand1 = get_operand_raw_value(tokens[0], molecule, record, field_count)
        raw_operand2 = get_operand_raw_value(tokens[2], molecule, record, field_count)
        operand1, operand2 = get_casted_operands(raw_operand1, raw_operand2)
        return OPERATION_MAP[tokens[1]](operand1, operand2)

    return [func]


def create_filters(filters: list[str]) -> list[Filter]:
    return [filter for filter_str in filters for filter in create_filter(filter_str)]


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    inputs = inputs_generator(args.infile)
    field_count = 0
    filters = create_filters(args.filter)
    for source in inputs:
        for record, molecule in enumerate(read_molecules(source), start=1):
            if args.s is not None and molecule.get_field(args.s) is not None:
                field_count += 1
            if any(filter(molecule, record, field_count) for filter in filters):
                print(repr(molecule))


if __name__ == "__main__":
    main()
