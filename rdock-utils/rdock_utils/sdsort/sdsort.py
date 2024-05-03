import itertools
import logging
import math
import sys
from dataclasses import dataclass
from typing import Iterable, TextIO

from rdock_utils.common import FastSDMol, inputs_generator, molecules_with_progress_log, read_molecules_from_all_inputs

logger = logging.getLogger("sdsort")


@dataclass(frozen=True)
class SDSortConfig:
    sorting_field: str
    reverse_sort: bool
    numeric_sort: bool
    fast_mode: bool
    group_key: str
    files: list[str] | None


class SDSort:
    def __init__(self, config: SDSortConfig, output: TextIO = sys.stdout) -> None:
        self.config = config
        self.output = output

    def run(self) -> None:
        inputs = inputs_generator(self.config.files)
        input_molecules = read_molecules_from_all_inputs(inputs)
        molecules = molecules_with_progress_log(input_molecules)
        sort_method = self.sort_records_fast_mode if self.config.fast_mode else self.sorted_records_normal
        sorted_records = sort_method(molecules)

        for molecule in sorted_records:
            molecule.write(self.output)

    def sorted_records_normal(self, molecules: Iterable[FastSDMol]) -> Iterable[FastSDMol]:
        return sorted(molecules, key=self.get_sorting_value, reverse=self.config.reverse_sort)

    def sort_records_fast_mode(self, molecules: Iterable[FastSDMol]) -> Iterable[FastSDMol]:
        grouped_molecules = itertools.groupby(molecules, key=lambda x: x.get_field(self.config.group_key))
        sorted_groups = (self.sorted_records_normal(group_records) for _, group_records in grouped_molecules)
        return itertools.chain.from_iterable(sorted_groups)

    def get_sorting_value(self, molecule: FastSDMol) -> float | str:
        field = molecule.get_field(self.config.sorting_field)

        if not self.config.numeric_sort:
            return field or ""

        try:
            if field is None:
                raise ValueError("Field is missing")
            return float(field)

        except ValueError as e:
            logger.warning(
                f"Field '{self.config.sorting_field}' for molecule {molecule.title}: {e} "
                "Defaulted to to infinity. Consider using sdfilter to remove invalid results"
            )
            return math.inf
