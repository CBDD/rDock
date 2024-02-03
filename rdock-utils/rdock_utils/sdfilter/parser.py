# Standard Library
import argparse
import logging
from dataclasses import dataclass

logger = logging.getLogger("sdfilter")


@dataclass
class SDFilterConfig:
    filter: str
    summary_field: str | None
    infile: list[str]


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filters SD records by data fields")
    filter_help = (
        "Filters can be provided as a string or in a file, one per line. All filters are OR'd together.\n"
        "Filters follow the format:\n"
        "'$<DataField> <Operator> <Value>'\n"
        "where valid operators are: '==', '!=', '<', '>', '<=', and '>=' for general values,\n"
        "'in' and 'not_in' for strings, and 'eq', 'ne', 'lt', 'gt', 'le', and 'ge' \n"
        "for strings for perl version retro-compatibility.\n"
        "_REC (record number), _TITLE1, _TITLE2, and _TITLE3 are provided as a pseudo-data field\n"
        "rdock-utils provides expanded functionality, where two data fields can be compared\n"
        "using the following syntax:\n"
        "'$<DataField1> <Operator> $<DataField2>'\n"
        "also, any combination of literal filters and filter files can be provided\n"
        "filter files including other filters are supported as well, so be careful with recursion\n"
    )
    parser.add_argument("-f", "--filter", type=str, help=filter_help, required=True)
    s_help = "If -s option is used, _COUNT (#occurrences of DataField) is provided as a pseudo-data field"
    parser.add_argument("-s", type=str, default=None, help=s_help)
    infile_help = "input file[s] to be processed. if not provided, stdin is used."
    parser.add_argument("infile", type=str, nargs="*", help=infile_help)
    return parser


def get_config(argv: list[str] | None = None) -> SDFilterConfig:
    parser = get_parser()
    args = parser.parse_args(argv)
    return SDFilterConfig(filter=args.filter, summary_field=args.s, infile=args.infile)
