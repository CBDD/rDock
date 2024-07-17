import argparse
from dataclasses import dataclass

Filter = dict[str, float]


@dataclass
class rbhtfinderConfig:
    input: str
    output: str
    threshold: str
    name: int
    filters: list[str]
    validation: int
    header: bool
    max_time: float
    min_percentage: float

    def __post_init__(self):
        self.filters = self.get_parsed_filters()

    def get_parsed_filters(self) -> list[Filter]:
        parsed_filters = [self._parse_filter(filter) for filter in self.filters]
        return parsed_filters

    @staticmethod
    def _parse_filter(filter_str: str) -> Filter:
        parsed_filter = {}

        for item in filter_str.split(","):
            key, value = item.split("=")
            parsed_filter[key] = float(value) if key in ["interval", "min", "max"] else int(value)
        # User inputs with 1-based numbering whereas python uses 0-based
        parsed_filter["column"] -= 1

        return parsed_filter


def get_parser() -> argparse.ArgumentParser:
    description = """
    Estimate the results and computation time of an rDock high-throughput protocol.

    Steps:
    1. Perform exhaustive docking of a small representative part of the entire library.
    2. Store the result of sdreport -t from that exhaustive docking run in a file 
       <sdreport_file>, which will be the input of this script.
    3. Run rbhtfinder, specifying -i <sdreport_file> and an arbitrary number of filters 
       using the -f option, for example, "-f column=6,steps=5,min=0.5,max=1.0,interval=0.1". 
       This example would simulate the effect of applying thresholds on column 6 after 5 poses 
       have been generated, for values between 0.5 and 1.0 (i.e., 0.5, 0.6, 0.7, 0.8, 0.9, 1.0). 
       More than one threshold can be specified, e.g., "-f column=4,steps=5,min=-12,max=-10,
       interval=1 column=4,steps=15,min=-16,max=-15,interval=1" will test the following 
       combinations of thresholds on column 4:
            5   -10     15      -15
            5   -11     15      -15
            5   -12     15      -15
            5   -10     15      -16
            5   -11     15      -16
            5   -12     15      -16
       The number of combinations will increase very rapidly, the more filters are used and the 
       larger the range of values specified for each. It may be sensible to run rbhtfinder several 
       times to explore the effects of various filters independently.

    Output:
    The output of the program consists of the following columns:
            FILTER1 NSTEPS1 THR1    PERC1   TOP500_SCORE.INTER  ENRICH_SCORE.INTER      TIME
            SCORE.INTER     5       -13.00  6.04    72.80   12.05   0.0500
            SCORE.INTER     5       -12.00  9.96    82.80   8.31    0.0500
    The four columns are repeated for each filter specified with the -f option: 
        name of the column on which the filter is applied (FILTER1), 
        number of steps at which the threshold is applied (NSTEPS1), 
        value of the threshold (THR1)   
        and the percentage of poses which pass this filter (PERC1). 
    Additional filters (FILTER2, FILTER3 etc.) are listed in the order that they are applied 
    (i.e., by NSTEPS).

    The final columns provide some overall statistics for the combination of thresholds 
    specified in a row. TOP500_SCORE.INTER gives the percentage of the top-scoring 500 poses, 
    measured by SCORE.INTER, from the whole of <sdreport_file> which are retained after the 
    thresholds are applied. This can be contrasted with the final PERC column. The higher the 
    ratio (the 'enrichment factor'), the better the combination of thresholds. If thresholds are 
    applied on multiple columns, this column will be duplicated for each, e.g. TOP500_SCORE.INTER 
    and TOP500_SCORE.RESTR will give the percentage of the top-scoring poses retained for both of 
    these scoring methods. The exact number of poses used for this validation can be changed from 
    the default 500 using the --validation flag.
    ENRICH_SCORE.INTER gives the enrichment factor as a quick rule-of-thumb to assess the best 
    choice of thresholds. The final column TIME provides an estimate of the time taken to perform 
    docking, as a proportion of the time taken for exhaustive docking. This value should be below 
    0.1.

    After a combination of thresholds has been selected, they need to be encoded into a threshold 
    file which rDock can use as an input. rbhtfinder attempts to help with this task by 
    automatically selecting a combination and writing a threshold file. The combination chosen is 
    that which provides the highest enrichment factor, after all options with a TIME value over 
    0.1 are excluded. This choice should not be blindly followed, so the threshold file should be 
    considered a template that the user modifies as needed.

    Requirements:
    rbhtfinder requires NumPy. Installation of pandas is recommended, but optional; if pandas is 
    not available, loading the input file for calculations will be considerably slower.
    """
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i",
        "--input",
        help="Input from sdreport (tabular separated format).",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file for report on threshold combinations.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        help="Threshold file used by rDock as input.",
        type=str,
    )
    parser.add_argument(
        "-n",
        "--name",
        type=int,
        default=1,  # Actually, index of molecule name in input file is 1 by default
        help="Index of column containing the molecule name. Default is 2.",
    )
    parser.add_argument(
        "-f",
        "--filter",
        nargs="+",
        # action="append" removed in favor of simplicity and to avoid redundancy (unnecessary nested structure),
        type=str,
        help="Filter to apply, e.g. column=4,steps=5,min=-10,max=-15,interval=1 will test applying a filter to column 4 after generation of 5 poses, with threshold values between -10 and -15 tested. The variables column, steps, min and max must all be specified; interval defaults to 1 if not given.",
    )
    parser.add_argument(
        "-v",
        "--validation",
        type=int,
        default=500,
        help="Top-scoring N molecules from input to use for validating threshold combinations. Default is 500.",
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help="Specify if the input file from sdreport contains a header line with column names. If not, output files will describe columns using indices, e.g. COL4, COL5.",
    )
    parser.add_argument(
        "--max-time",
        type=float,
        default=0.1,
        help="Maximum value for time to use when autogenerating a high-throughput protocol - default is 0.1, i.e. 10%% of the time exhaustive docking would take.",
    )
    parser.add_argument(
        "--min-perc",
        type=float,
        default=1.0,
        help="Minimum value for the estimated final percentage of compounds to use when autogenerating a high-throughput protocol - default is 1.",
    )

    args = parser.parse_args()