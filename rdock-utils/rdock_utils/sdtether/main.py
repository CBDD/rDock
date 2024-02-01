from rdock_utils.sdtether.parser import get_config
from rdock_utils.sdtether.sdtether import SDTether


def main(argv: list[str] | None = None) -> None:
    config = get_config(argv)
    sdtether = SDTether(config)
    sdtether.run()


if __name__ == "__main__":
    main()
