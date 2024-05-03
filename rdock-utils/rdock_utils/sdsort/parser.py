import argparse
from dataclasses import dataclass


@dataclass(frozen=True)
class SDSortConfig:
    sorting_field: str
    reverse_sort: bool
    numeric_sort: bool
    fast_mode: bool
    group_key: str
    files: list[str] | None


def get_parser() -> argparse.ArgumentParser:
    description = """Sorts SD records by a specified data field.
    Notes:
    - If no SD file list is provided, the script reads from standard input.
    - Output is directed to standard output.
    - Fast mode can be safely used for partial sorting of large SD files of raw docking hits without encountering memory issues. It will sort together consecutive molecules with the same value for a grouping key (the title, by default) instead of sorting all provided molecules together.
    """
    parser = argparse.ArgumentParser(description=description)
    sorting_field_help = "Specify the field for sorting"
    parser.add_argument("--field", "-f", default="SCORE", metavar="DataField", type=str, help=sorting_field_help)
    parser.add_argument("--reverse", "-r", action="store_true", help="Perform a descending sort (default: ascending)")
    parser.add_argument("--numeric", "-n", action="store_true", help="Perform a numeric sort (default: text sort)")
    fast_mode_help = "Enable fast mode: Sort records for each named compound independently (must be consecutive)"
    parser.add_argument("--fast", "-s", action="store_true", help=fast_mode_help)
    name_field_help = "Specify the grouping field for fast sorting mode (default: _TITLE1)"
    parser.add_argument("--group-key", "-id", default="_TITLE1", metavar="NameField", type=str, help=name_field_help)
    infile_help = "input file[s] to be processed. if not provided, stdin is used."
    parser.add_argument("files", nargs="*", type=str, help=infile_help)
    return parser


def get_config(argv: list[str] | None = None) -> SDSortConfig:
    parser = get_parser()
    args = parser.parse_args(argv)
    return SDSortConfig(args.field, args.reverse, args.numeric, args.fast, args.group_key, args.files)
