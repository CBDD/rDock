from rdock_utils.common import inputs_generator, read_molecules_from_all_inputs

from .filter import ExpressionContext, create_filters, molecules_with_context
from .parser import get_config


def main(argv: list[str] | None = None) -> None:
    config = get_config(argv)
    inputs = inputs_generator(config.infile)
    context = ExpressionContext(config.summary_field)
    filters = create_filters(config.filter, context)
    molecules = molecules_with_context(read_molecules_from_all_inputs(inputs), context)
    for molecule in molecules:
        if any(filter.evaluate(molecule) for filter in filters):
            print(repr(molecule))


if __name__ == "__main__":
    main()
