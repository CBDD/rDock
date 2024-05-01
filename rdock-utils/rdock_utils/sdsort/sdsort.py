import itertools
import logging
import math
import sys
from typing import Iterable, TextIO

from ..common import FastSDMol, inputs_generator, molecules_with_progress_log, read_molecules_from_all_inputs
from .parser import SDSortConfig

logger = logging.getLogger("sdsort")


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

        if self.config.numeric_sort:
            try:
                if field is None:
                    raise ValueError
                else:
                    return float(field)

            except ValueError:
                logger_msg = (
                    f"Molecule {molecule.title} has no field '{self.config.sorting_field}'. "
                    "Defaulted to to infinity. Consider using sdfilter to remove invalid results"
                )
                logger.warning(logger_msg)
                return math.inf

        return field or ""
