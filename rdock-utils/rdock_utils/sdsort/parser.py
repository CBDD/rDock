import argparse
import sys
from dataclasses import dataclass


@dataclass(frozen=True)
class SDSortConfig:
    descending_sort: bool
    numeric_sort: bool
    fast_mode: bool
    name_field: str
    data_field: str
    files: list[str]


def get_parser() -> argparse.ArgumentParser:
    description = """Sorts SD records by a specified data field.
    Notes:
    - "_REC" (record #) is provided as a pseudo-data field.
    - If no SD file list is provided, the script reads from standard input.
    - Output is directed to standard output.
    - Fast mode can be safely used for partial sorting of large SD files of raw docking hits without encountering memory issues.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-r", action="store_true", help="Perform a descending sort (default: ascending)")
    parser.add_argument("-n", action="store_true", help="Perform a numeric sort (default: text sort)")
    s_help = "Enable fast mode: Sort records for each named compound independently (must be consecutive)"
    parser.add_argument("-s", action="store_true", help=s_help)
    id_help = "Specify the field for compound names (default: 1st title line)"
    parser.add_argument("-id", metavar="NameField", type=str, help=id_help)
    parser.add_argument("-f", metavar="DataField", type=str, help="Specify the field for sorting")
    parser.add_argument("files", nargs="*", type=str, help="List of SD files to sort")
    return parser

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)


def get_config(argv: list[str] | None = None) -> SDSortConfig:
    parser = get_parser()
    args = parser.parse_args(argv)
    return SDSortConfig(args.r, args.n, args.s, args.id, args.f, args.files)
