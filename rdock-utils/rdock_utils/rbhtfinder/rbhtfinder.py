import itertools
import logging
import multiprocessing
import os
from collections import Counter, defaultdict
from functools import partial

import numpy as np
import pandas as pd

from rdock_utils.common import (
    Array1DFloat,
    Array1DInt,
    Array2DFloat,
    Array3DFloat,
    ColumnNamesArray,
    FilterCombination,
    InputData,
    MinMaxValues,
    MinScoreIndices,
    SDReportArray,
)

from .models import Filter, FilterCombinationResult, RBHTFinderConfig

logger = logging.getLogger("RBHTFinder")


class RBHTFinder:
    def __init__(self, config: RBHTFinderConfig) -> None:
        self.config = config

    def run(self) -> None:
        filters_combinations = self.generate_filters_combinations(self.config.filters)
        print(f"{len(filters_combinations)} combinations of filters calculated.")
        distinct_combinations = self.remove_redundant_combinations(filters_combinations, self.config.filters)

        if len(distinct_combinations) == 0:
            raise RuntimeError("No filter combinations could be calculated - check the thresholds specified.")

        print(
            f"{len(distinct_combinations)} combinations of filters remain after removal of redundant combinations. "
            "Starting calculations..."
        )
        sdreport_array, column_names = self.read_data()
        print(f"First few rows of the input array:\n{sdreport_array[:5]}")
        print("Data read in from input file.")
        # Convert to 3D array (molecules x poses x columns)
        molecule_array = self.prepare_array(sdreport_array, self.config.name)
        # Find the top scoring compounds for validation of the filter combinations
        columns = set(filter["column"] for filter in self.config.filters)
        min_score_indices = {
            column: np.argpartition(np.min(molecule_array[:, :, column], axis=1), self.config.validation)[
                : self.config.validation
            ]
            for column in columns
        }
        results = self.process_filter_combinations(molecule_array, min_score_indices, distinct_combinations)
        self.write_output(results, column_names)
        best_filter_combination_index = self.get_best_filter_combination_index(results)
        if self.config.threshold is not None:
            num_poses = molecule_array.shape[1]
            self.handle_threshold(best_filter_combination_index, distinct_combinations, column_names, num_poses)

    def process_filter_combinations(
        self, molecule_array: Array3DFloat, min_score_indices: MinScoreIndices, distinct_combinations: Array2DFloat
    ) -> list[FilterCombinationResult]:
        num_cpus = os.cpu_count() or 1
        with multiprocessing.Pool(num_cpus) as pool:
            function_to_apply = partial(
                self.calculate_results_for_filter_combination,
                molecule_array=molecule_array,
                min_score_indices=min_score_indices,
            )
            results = pool.map(function_to_apply, distinct_combinations)
        return results

    def handle_threshold(
        self,
        combination_index: int,
        distinct_combinations: Array2DFloat,
        column_names: ColumnNamesArray,
        num_poses: int,
    ) -> None:
        if combination_index:
            best_filter_combination: Array1DFloat = distinct_combinations[combination_index]
            self.write_threshold(best_filter_combination, column_names, num_poses)
        else:
            message = "Filter combinations defined are too strict or would take too long to run; no threshold file was written."
            logger.warning(message)

    def read_data(self) -> InputData:
        try:
            data_array, column_names = self.read_data_using_pandas()
        except Exception as e:
            logging.error(f"Error reading data with pandas: {e}")
            data_array, column_names = self.read_data_using_numpy()
        return data_array, column_names

    def read_data_using_pandas(self) -> InputData:
        sdreport_dataframe = pd.read_csv(self.config.input, sep="\t", header=0 if self.config.header else None)

        if self.config.header:
            column_names = sdreport_dataframe.columns.values
        else:
            # Use index names; add 1 to deal with zero-based numbering
            column_names = [f"COL{n+1}" for n in range(len(sdreport_dataframe.columns))]

        sdreport_array = sdreport_dataframe.values
        return sdreport_array, column_names

    def read_data_using_numpy(self) -> InputData:
        np_array = np.loadtxt(self.config.input, dtype=str)

        if self.config.header:
            column_names = np_array[0]
            sdreport_array = np_array[1:]
        else:
            column_names = [f"COL{n+1}" for n in range(np_array.shape[1])]
            sdreport_array = np_array

        return sdreport_array, column_names

    def apply_threshold(self, scored_poses: Array3DFloat, column: int, steps: int, threshold: float) -> Array1DInt:
        """
        Filter out molecules from `scored_poses`, where the minimum score reached (for a specified `column`) after `steps` is more negative than `threshold`.
        """
        # Minimum score after `steps` per molecule
        mins = np.min(scored_poses[:, :steps, column], axis=1)
        # Return those molecules where the minimum score is less than the threshold
        passing_molecules = np.where(mins < threshold)[0]
        return passing_molecules

    def prepare_array(self, data_array: SDReportArray, name_column: int) -> Array3DFloat:
        """
        Convert `sdreport_array` (read directly from the tsv) to 3D array (molecules x poses x columns) and filter out molecules with too few/many poses
        """
        split_indices = (
            np.where(
                data_array[:, name_column] != np.hstack((data_array[1:, name_column], data_array[0, name_column]))
            )[0]
            + 1
        )
        split_array = np.split(data_array, split_indices)
        modal_shape = Counter([array.shape for array in split_array]).most_common(1)[0]
        number_of_poses = modal_shape[0][0]  # Find modal number of poses per molecule in the array
        valid_split_arrays = [
            np.array_split(array, array.shape[0] / number_of_poses)  # type: ignore
            for array in split_array
            if not array.shape[0] % number_of_poses and array.shape[0]
        ]
        flattened_split_array = np.concatenate(valid_split_arrays)

        if len(flattened_split_array) * number_of_poses < data_array.shape[0] * 0.99:
            message = (
                "WARNING: The number of poses provided per molecule is inconsistent. "
                f"Only {len(flattened_split_array)} of {int(data_array.shape[0] / number_of_poses)} molecules have {number_of_poses} poses."
            )
            logger.warning(message)

        molecule_3d_array = np.array(flattened_split_array)
        # Overwrite the name column (should be the only one with dtype=str) so we can force everything to float
        molecule_3d_array[:, :, name_column] = 0
        final_array = molecule_3d_array.astype(float)
        return final_array

    def calculate_results_for_filter_combination(
        self, filter_combination: Array2DFloat, molecule_array: Array3DFloat, min_score_indices: MinScoreIndices
    ) -> FilterCombinationResult:
        """
        For a particular combination of filters, calculate the percentage of molecules that will be filtered, the percentage of top-scoring molecules that will be filtered, and the time taken relative to exhaustive docking
        """
        num_molecules = molecule_array.shape[0]
        num_steps = molecule_array.shape[1]
        # Passing_molecule_indices is a list of indices of molecules which have passed the applied filters. As more filters are applied, it gets smaller. Before any iteration, we initialise with all molecules passing
        passing_molecule_indices = np.arange(num_molecules)
        filter_percentages = []
        number_of_simulated_poses = 0  # Number of poses which we calculate would be generated, we use this to calculate the TIME column in the final output

        for i, threshold in enumerate(filter_combination):
            number_of_simulated_poses += self.calculate_simulated_poses_increment(i, passing_molecule_indices)
            column: int = self.config.filters[i]["column"]
            step: int = self.config.filters[i]["steps"]
            passing_indices = self.apply_threshold(molecule_array, column, step, threshold)
            # All mols which pass the threshold and which were already in passing_molecule_indices, i.e. passed all previous filters
            passing_molecule_indices = np.intersect1d(passing_molecule_indices, passing_indices, assume_unique=True)
            filter_percentages.append(len(passing_molecule_indices) / num_molecules)

        number_of_simulated_poses += len(passing_molecule_indices) * (num_steps - self.config.filters[-1]["steps"])
        perc_val = {
            k: len(np.intersect1d(v, passing_molecule_indices, assume_unique=True)) / self.config.validation
            for k, v in min_score_indices.items()
        }
        time = float(number_of_simulated_poses / np.prod(molecule_array.shape[:2]))
        result = FilterCombinationResult(
            combination=filter_combination,
            perc_val=perc_val,
            percentages=filter_percentages,
            time=time,
        )
        return result

    def calculate_simulated_poses_increment(self, index: int, passing_molecule_indices: Array1DInt) -> int:
        if index:
            # e.g. if there are 5000 mols left after 15 steps and the last filter was at 5 steps, append 5000 * (15 - 5) to number_of_simulated_poses
            increment: int = len(passing_molecule_indices) * (
                self.config.filters[index]["steps"] - self.config.filters[index - 1]["steps"]
            )
        else:
            increment = len(passing_molecule_indices) * self.config.filters[index]["steps"]
        return increment

    def write_output(
        self,
        results: list[FilterCombinationResult],
        column_names: ColumnNamesArray,
        sep: str = "\t",
        end: str = "\n",
    ) -> None:
        """
        Print results as a table. The number of columns varies depending how many columns the user picked.
        """
        with open(self.config.output, "w") as f:
            header = self._get_output_header(results[0], column_names)
            f.write(sep.join(header) + end)
            content_lines = [sep.join(self._get_output_content(result, column_names)) + end for result in results]
            f.writelines(content_lines)

    def _get_output_header(self, result: FilterCombinationResult, column_names: ColumnNamesArray) -> list[str]:
        header = []
        for i in range(len(result.combination)):
            header.extend([f"FILTER{i + 1}", f"NSTEPS{i + 1}", f"THR{i + 1}", f"PERC{i + 1}"])

        for col_index in result.perc_val.keys():
            header.append(f"TOP{self.config.validation}_{column_names[col_index]}")
            header.append(f"ENRICH_{column_names[col_index]}")

        header.append("TIME")

        return header

    def _get_output_content(self, result: FilterCombinationResult, column_names: ColumnNamesArray) -> list[str]:
        content = []

        for i, threshold in enumerate(result.combination):
            column_name = column_names[self.config.filters[i]["column"]]
            steps = self.config.filters[i]["steps"]
            filter_percentage = result.percentages[i] * 100
            content.extend([f"{column_name}", f"{steps}", f"{threshold:.2f}", f"{filter_percentage:.2f}"])

        for value in result.perc_val.values():
            perc_val_percent = value * 100
            enrichment = value / result.percentages[-1] if result.percentages[-1] else float("nan")
            content.append(f"{perc_val_percent:.2f}")
            content.append(f"{enrichment:.2f}")

        content.append(f"{result.time:.4f}")

        return content

    def get_best_filter_combination_index(self, results: list[FilterCombinationResult]) -> int:
        """
        Very debatable how to do this...
        Here we exclude all combinations with TIME < max_time and calculate an "enrichment factor"
        (= percentage of validation compounds / percentage of all compounds); we select the
        threshold with the highest enrichment factor
        """
        min_max_values: MinMaxValues = {}
        # Transpose the `perc_val` data to get columns
        perc_vals = {col: [result.perc_val[col] for result in results] for col in results[0].perc_val}
        min_max_values.update({col: {"min": min(vals), "max": max(vals)} for col, vals in perc_vals.items()})
        time_vals = [result.time for result in results]
        min_max_values["time"] = {"min": min(time_vals), "max": max(time_vals)}
        combination_scores = [self.calculate_combination_score(result, min_max_values) for result in results]
        index = np.argmax(combination_scores)
        return int(index)

    def calculate_combination_score(self, result: FilterCombinationResult, min_max_values: MinMaxValues) -> float:
        if result.time < self.config.max_time and result.percentages[-1] >= self.config.min_percentage / 100:
            col_scores = [
                (result.perc_val[col] - min_max_values[col]["min"])
                / (min_max_values[col]["max"] - min_max_values[col]["min"])
                for col in min_max_values
                if col != "time"
            ]
            time_score = (min_max_values["time"]["max"] - result.time) / (
                min_max_values["time"]["max"] - min_max_values["time"]["min"]
            )
            score = sum(col_scores) + time_score
        else:
            score = 0

        return score

    def write_threshold(
        self,
        best_filter_combination: Array1DFloat,
        column_names: ColumnNamesArray,
        max_number_of_runs: int,
        sep: str = "\n",
        end: str = "\n",
    ) -> None:
        path = self.config.threshold or "default_threshold.txt"
        with open(path, "w") as f:
            content = self._get_threshold_content(best_filter_combination, column_names, max_number_of_runs)
            f.write(sep.join(content) + end)

    def _get_threshold_content(
        self, best_filter_combination: Array1DFloat, column_names: ColumnNamesArray, max_number_of_runs: int
    ) -> list[str]:
        content = []
        # Number of filters to apply
        content.append(f"{len(self.config.filters) + 1}")
        # Get each filter to a separate line
        filter_lines = [
            f'if - {best_filter_combination[i]:.2f} {column_names[filter["column"]]} 1.0 '
            f'if - SCORE.NRUNS  {filter["steps"]} 0.0 -1.0,'
            for i, filter in enumerate(self.config.filters)
        ]
        content.extend(filter_lines)
        # Filter to terminate docking when NRUNS reaches the number of runs used in the input file
        content.append(f"if - SCORE.NRUNS {max_number_of_runs - 1} 0.0 -1.0")
        # Find strictest filters for all columns and apply them again
        filters_by_column = defaultdict(list)
        for i, filter in enumerate(self.config.filters):
            col = filter["column"]
            filters_by_column[col].append(best_filter_combination[i])
        # Number of filters (same as number of columns filtered on)
        content.append(f"{len(filters_by_column)}")
        # Filter
        filter_min_values = [f"- {column_names[col]} {min(values)}," for col, values in filters_by_column.items()]
        content.extend(filter_min_values)
        return content

    def generate_filters_combinations(self, filters: list[Filter]) -> list[FilterCombination]:
        filter_ranges = ((filter["min"], filter["max"] + filter["interval"], filter["interval"]) for filter in filters)
        combinations = (np.arange(*range) for range in filter_ranges)
        filters_combinations = list(itertools.product(*combinations))
        return filters_combinations

    def remove_redundant_combinations(
        self, all_combinations: list[FilterCombination], filters: list[Filter]
    ) -> Array2DFloat:
        all_combinations_array = np.array(all_combinations)
        columns = [filter["column"] for filter in filters]
        indices_per_col = {col: [i for i, c in enumerate(columns) if c == col] for col in set(columns)}
        # Create a mask to keep only valid combinations
        mask = np.ones(len(all_combinations_array), dtype=bool)

        for indices in indices_per_col.values():
            col_data = all_combinations_array[:, indices]
            sorted_data = np.sort(col_data, axis=1)[:, ::-1]  # Sort descending
            is_valid = np.all(col_data == sorted_data, axis=1)  # Check if sorted matches original
            is_unique = np.apply_along_axis(lambda x: len(set(x)) == len(x), 1, col_data)
            mask &= is_valid & is_unique

        valid_combinations: Array2DFloat = all_combinations_array[mask]
        return valid_combinations
