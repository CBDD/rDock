from .parser import get_config
from .rbhtfinder import RBHTFinder


def main(argv: list[str] | None = None) -> None:
    config = get_config(argv)
    rbhtfinder = RBHTFinder(config)
    rbhtfinder.run()


if __name__ == "__main__":
    main()
