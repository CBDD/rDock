from .parser import get_parser
from .sdrmsd import SDRMSD


def main(argv: list[str] | None = None) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    sdrmsd = SDRMSD(args.reference, args.input, args.fit, args.threshold, args.out)
    sdrmsd.run()


if __name__ == "__main__":
    main()
