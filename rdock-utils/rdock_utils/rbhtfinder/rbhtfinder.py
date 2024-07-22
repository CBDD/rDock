import itertools
import logging
import multiprocessing
import os
from collections import Counter
from dataclasses import dataclass
from functools import partial

import numpy as np
import pandas as pd

from .parser import Filter, RBHTFinderConfig

logger = logging.getLogger("RBHTFinder")

InputData = tuple[np.ndarray, list[str]]


@dataclass
class RBHTResult:
    filter_combination: np.ndarray
    perc_val: dict[str, int]
    filter_percentages: list[int]
    time: float


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

        # convert to 3D array (molecules x poses x columns)
        molecule_array = self.prepare_array(sdreport_array, self.config.name)

        # find the top scoring compounds for validation of the filter combinations
        min_score_indices: dict[float, np.ndarray] = {}
        for column in set(filter["column"] for filter in self.config.filters):
            min_scores = np.min(molecule_array[:, :, column], axis=1)
            min_score_indices[column] = np.argpartition(min_scores, self.config.validation)[: self.config.validation]

        pool = multiprocessing.Pool(os.cpu_count())
        results = pool.map(
            partial(
                self.calculate_results_for_filter_combination,
                molecule_array=molecule_array,
                min_score_indices=min_score_indices,
            ),
            distinct_combinations,
        )

        self.write_output(results, self.config.filters, self.config.validation, self.config.output, column_names)

        best_filter_combination = self.select_best_filter_combination(
            results, self.config.max_time, self.config.min_percentage
        )
        if self.config.threshold:
            if best_filter_combination:
                self.write_threshold_file(
                    self.config.filters,
                    distinct_combinations[best_filter_combination],
                    self.config.threshold,
                    column_names,
                    molecule_array.shape[1],
                )
            else:
                print(
                    "Filter combinations defined are too strict or would take too long to run; no threshold file was written."
                )
                exit(1)

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
            # use index names; add 1 to deal with zero-based numbering
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

    def apply_threshold(self, scored_poses, column, steps, threshold):
        """
        Filter out molecules from `scored_poses`, where the minimum score reached (for a specified `column`) after `steps` is more negative than `threshold`.
        """
        # minimum score after `steps` per molecule
        mins = np.min(scored_poses[:, :steps, column], axis=1)
        # return those molecules where the minimum score is less than the threshold
        passing_molecules = np.where(mins < threshold)[0]
        return passing_molecules

    def prepare_array(self, data_array: np.ndarray, name_column: int) -> np.ndarray:
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
        modal_shape = Counter([n.shape for n in split_array]).most_common(1)[0]
        number_of_poses = modal_shape[0][0]  # Find modal number of poses per molecule in the array
        split_array_clean = sum(
            [
                np.array_split(n, n.shape[0] / number_of_poses)
                for n in split_array
                if not n.shape[0] % number_of_poses and n.shape[0]
            ],
            [],
        )

        if len(split_array_clean) * number_of_poses < data_array.shape[0] * 0.99:
            message = (
                "WARNING: The number of poses provided per molecule is inconsistent. "
                f"Only {len(split_array_clean)} of {int(data_array.shape[0] / number_of_poses)} moleules have {number_of_poses} poses."
            )
            logger.warning(message)

        molecule_3d_array = np.array(split_array_clean)
        # Overwrite the name column (should be the only one with dtype=str) so we can force everything to float
        molecule_3d_array[:, :, name_column] = 0
        final_array = molecule_3d_array.astype(float)
        return final_array

    def calculate_results_for_filter_combination(
        self, filter_combination: np.ndarray, molecule_array: np.ndarray, min_score_indices: dict[float, np.ndarray]
    ) -> RBHTResult:
        """
        For a particular combination of filters, calculate the percentage of molecules that will be filtered, the percentage of top-scoring molecules that will be filtered, and the time taken relative to exhaustive docking
        """
        # mols_passed_threshold is a list of indices of molecules which have passed the applied filters. As more filters are applied, it gets smaller. Before any iteration, we initialise with all molecules passing
        mols_passed_threshold = list(range(molecule_array.shape[0]))
        filter_percentages = []
        number_of_simulated_poses = 0  # number of poses which we calculate would be generated, we use this to calculate the TIME column in the final output
        for n, threshold in enumerate(filter_combination):
            if n:
                # e.g. if there are 5000 mols left after 15 steps and the last filter was at 5 steps, append 5000 * (15 - 5) to number_of_simulated_poses
                number_of_simulated_poses += len(mols_passed_threshold) * (
                    self.config.filters[n]["steps"] - self.config.filters[n - 1]["steps"]
                )
            else:
                number_of_simulated_poses += len(mols_passed_threshold) * self.config.filters[n]["steps"]
            mols_passed_threshold = [  # all mols which pass the threshold and which were already in mols_passed_threshold, i.e. passed all previous filters
                n
                for n in self.apply_threshold(
                    molecule_array, self.config.filters[n]["column"], self.config.filters[n]["steps"], threshold
                )
                if n in mols_passed_threshold
            ]
            filter_percentages.append(len(mols_passed_threshold) / molecule_array.shape[0])
        number_of_simulated_poses += len(mols_passed_threshold) * (
            molecule_array.shape[1] - self.config.filters[-1]["steps"]
        )
        perc_val = {
            k: len([n for n in v if n in mols_passed_threshold]) / self.config.validation
            for k, v in min_score_indices.items()
        }
        time = number_of_simulated_poses / np.product(molecule_array.shape[:2])
        rbhtresult = RBHTResult(
            filter_combination=filter_combination, perc_val=perc_val, filter_percentages=filter_percentages, time=time
        )
        return rbhtresult

    def write_output(self, results: list[RBHTResult], filters, number_of_validation_mols, output_file, column_names):
        """
        Print results as a table. The number of columns varies depending how many columns the user picked.
        """
        with open(output_file, "w") as f:
            # write header
            for n in range(len(results[0].filter_combination)):
                f.write(f"FILTER{n+1}\tNSTEPS{n+1}\tTHR{n+1}\tPERC{n+1}\t")
            for n in results[0].perc_val:
                f.write(f"TOP{number_of_validation_mols}_{column_names[n]}\t")
                f.write(f"ENRICH_{column_names[n]}\t")
            f.write("TIME\n")

            # write results
            for result in results:
                for n, threshold in enumerate(result.filter_combination):
                    f.write(
                        f"{column_names[filters[n]['column']]}\t{filters[n]['steps']}\t{threshold:.2f}\t{result.filter_percentages[n]*100:.2f}\t"
                    )
                for n in result.perc_val:
                    f.write(f"{result.perc_val[n]*100:.2f}\t")
                    if result.filter_percentages[-1]:
                        f.write(f"{result.perc_val[n]/result.filter_percentages[-1]:.2f}\t")
                    else:
                        f.write("NaN\t")
                f.write(f"{result.time:.4f}\n")
        return

    def select_best_filter_combination(self, results: list[RBHTResult], max_time, min_perc):
        """
        Very debatable how to do this...
        Here we exclude all combinations with TIME < max_time and calculate an "enrichment factor"
        (= percentage of validation compounds / percentage of all compounds); we select the
        threshold with the highest enrichment factor
        """
        min_max_values = {}
        for col in results[0].perc_val.keys():
            vals = [result.perc_val[col] for result in results]
            min_max_values[col] = {"min": min(vals), "max": max(vals)}
        time_vals = [result.time for result in results]
        min_max_values["time"] = {"min": min(time_vals), "max": max(time_vals)}

        combination_scores = [
            sum(
                [
                    (
                        (result.perc_val[col] - min_max_values[col]["min"])
                        / (min_max_values[col]["max"] - min_max_values[col]["min"])
                    )
                    for col in results[0].perc_val.keys()
                ]
                + [
                    (min_max_values["time"]["max"] - result.time)
                    / (min_max_values["time"]["max"] - min_max_values["time"]["min"])
                ]
            )
            if result.time < max_time and result.filter_percentages[-1] >= min_perc / 100
            else 0
            for result in results
        ]
        return np.argmax(combination_scores)

    def write_threshold_file(self, filters, best_filter_combination, threshold_file, column_names, max_number_of_runs):
        with open(threshold_file, "w") as f:
            # write number of filters to apply
            f.write(f"{len(filters) + 1}\n")
            # write each filter to a separate line
            for n, filtr in enumerate(filters):
                f.write(
                    f'if - {best_filter_combination[n]:.2f} {column_names[filtr["column"]]} 1.0 if - SCORE.NRUNS  {filtr["steps"]} 0.0 -1.0,\n'
                )
            # write filter to terminate docking when NRUNS reaches the number of runs used in the input file
            f.write(f"if - SCORE.NRUNS {max_number_of_runs - 1} 0.0 -1.0\n")

            # write final filters - find strictest filters for all columns and apply them again
            filters_by_column = {
                col: [best_filter_combination[n] for n, filtr in enumerate(filters) if filtr["column"] == col]
                for col in set([filtr["column"] for filtr in filters])
            }
            # write number of filters (same as number of columns filtered on)
            f.write(f"{len(filters_by_column)}\n")
            # write filter
            for col, values in filters_by_column.items():
                f.write(f"- {column_names[col]} {min(values)},\n")

    def generate_filters_combinations(self, filters: list[Filter]) -> list[tuple]:
        filter_ranges = ((filter["min"], filter["max"] + filter["interval"], filter["interval"]) for filter in filters)
        combinations = (np.arange(*range) for range in filter_ranges)
        filters_combinations = list(itertools.product(*combinations))
        return filters_combinations

    def remove_redundant_combinations(self, all_combinations: list[tuple], filters: list[Filter]) -> np.ndarray:
        all_combinations_array = np.array(all_combinations)
        columns = [filter["column"] for filter in filters]
        indices_per_col = {col: [i for i, c in enumerate(columns) if c == col] for col in set(columns)}

        # Create a mask to keep only valid combinations
        mask = np.ones(len(all_combinations_array), dtype=bool)

        for _, indices in indices_per_col.items():
            col_data = all_combinations_array[:, indices]
            sorted_data = np.sort(col_data, axis=1)[:, ::-1]  # sort descending
            is_valid = np.all(col_data == sorted_data, axis=1)  # Check if sorted matches original
            is_unique = np.apply_along_axis(lambda x: len(set(x)) == len(x), 1, col_data)
            mask &= is_valid & is_unique

        valid_combinations = all_combinations_array[mask]
        return valid_combinations
