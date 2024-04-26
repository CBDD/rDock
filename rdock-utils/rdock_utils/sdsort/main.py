from .parser import get_config
from .sdsort import SDSort


def main(argv: list[str] | None = None) -> None:
    config = get_config(argv)
    sdsort = SDSort(config)
    sdsort.run()


if __name__ == "__main__":
    main()
