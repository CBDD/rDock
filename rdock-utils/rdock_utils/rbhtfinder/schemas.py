from dataclasses import dataclass, field

from rdock_utils.common import Array1DFloat


@dataclass
class Filter:
    column: int = 0
    steps: int = 0
    min: float = 0.0
    max: float = 0.0
    interval: float = 0.0


@dataclass
class FilterCombinationResult:
    combination: Array1DFloat
    perc_val: dict[int, float]
    percentages: list[float]
    time: float


@dataclass
class MinMax:
    min: float
    max: float


@dataclass
class MinMaxValues:
    values: dict[int | str, MinMax] = field(default_factory=dict)

    def update(self, column: int | str, min_val: float, max_val: float) -> None:
        self.values[column] = MinMax(min=min_val, max=max_val)

    def get(self, column: int | str) -> MinMax:
        return self.values[column]
