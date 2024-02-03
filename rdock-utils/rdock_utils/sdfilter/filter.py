# Standard Library
import logging
import operator
import shlex
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Iterable

# Local imports
from rdock_utils.common.SDFParser import FastSDMol

logger = logging.getLogger("sdfilter")


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


class ExpressionContext:
    def __init__(self, summary_field: str | None = None):
        self.summary_field = summary_field
        self.summary_counter: Counter[Any] = Counter()
        self.record = 0

    def _get_symbol_as_field(self, symbol: str, molecule: FastSDMol) -> str:
        value = molecule.get_field(symbol[1:])
        if value is None:
            logger.warning(f"field {symbol} not found in record {self.record}, assuming literal string")
            return symbol
        return value

    def get_operand_raw_value(self, operand: str, molecule: FastSDMol) -> str | int:
        if operand.startswith("$"):
            return self.get_symbol_value(operand, molecule)
        else:
            return operand

    def get_symbol_value(self, symbol: str, molecule: FastSDMol) -> str | int:
        match symbol:
            case "$_REC":
                return self.record
            case "$_COUNT":
                if self.summary_field is None:
                    raise ValueError("summary field not provided")
                summary_field_value = molecule.get_field(self.summary_field)
                return self.summary_counter[summary_field_value]
            case _:
                return self._get_symbol_as_field(symbol, molecule)


class FilterExpression:
    def __init__(self, operand1: str, operator: str, operand2: str, context: ExpressionContext):
        self.operand1 = operand1
        self.operator = operator
        self.operand2 = operand2
        self.context = context

    def evaluate(self, molecule: FastSDMol) -> bool:
        raw_operand1 = self.context.get_operand_raw_value(self.operand1, molecule)
        raw_operand2 = self.context.get_operand_raw_value(self.operand2, molecule)
        operand1, operand2 = get_casted_operands(raw_operand1, raw_operand2)
        return OPERATION_MAP[self.operator](operand1, operand2)


def create_filters(filter_str: str, context: ExpressionContext) -> list[FilterExpression]:
    tokens = shlex.split(filter_str)
    if len(tokens) == 1 and (path := Path(filter_str)).is_file():
        with open(path, "r") as f:
            return [filter for line in f for filter in create_filters(line.strip(), context)]
    elif len(tokens) != 3:
        raise ValueError(f"invalid filter: {filter_str}")

    if tokens[1] not in OPERATION_MAP:
        raise ValueError(f"invalid operator: {tokens[1]}. expected: {OPERATION_MAP.keys()}")

    return [FilterExpression(tokens[0], tokens[1], tokens[2], context)]


def molecules_with_context(molecules: Iterable[FastSDMol], context: ExpressionContext) -> Iterable[FastSDMol]:
    for molecule in molecules:
        context.record += 1
        if context.summary_field is not None:
            context.summary_counter[molecule.get_field(context.summary_field)] += 1
        yield molecule
