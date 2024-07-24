from dataclasses import dataclass
from typing import Any

from rdock_utils.common import Array1DFloat

Filter = dict[str, Any]  # The type for the values is either 'float' or 'int'; 'Any' is used to comply with mypy


@dataclass
class RBHTFinderConfig:
    input: str
    output: str
    threshold: str | None
    name: int
    filters: list[Filter]
    validation: int
    header: bool
    max_time: float
    min_percentage: float

    def __post_init__(self) -> None:
        self.filters = self.get_parsed_filters()

    def get_parsed_filters(self) -> list[Filter]:
        filter_args: list[str] = self.filters  # type: ignore
        parsed_filters = [self._parse_filter(filter) for filter in filter_args]
        # sort filters by step at which they are applied
        parsed_filters.sort(key=lambda n: n["steps"])
        return parsed_filters

    @staticmethod
    def _parse_filter(filter_str: str) -> Filter:
        parsed_filter = {}

        for item in filter_str.split(","):
            key, value = item.split("=")
            parsed_filter[key] = float(value) if key in ("interval", "min", "max") else int(value)
        # User inputs with 1-based numbering whereas python uses 0-based
        parsed_filter["column"] -= 1

        if "interval" not in parsed_filter:
            parsed_filter["interval"] = 1.0

        return parsed_filter


@dataclass
class FilterCombinationResult:
    combination: Array1DFloat
    perc_val: dict[int, float]
    percentages: list[float]
    time: float
