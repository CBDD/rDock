def generate_result(ids: list[int], data: dict[str, str], fast_mode: bool = False) -> list[tuple[str, str]]:
    if fast_mode:
        return [(data[f"id{id}"][0], data[f"id{id}"][1]) for id in ids]

    return [(f"MOL{id}", data[f"MOL{id}"]) for id in ids]


def get_data(file_path: str, field: str) -> dict[str, str]:
    data = {}

    with open(file_path, "r") as file:
        lines = file.readlines()
        current_molecule = None

        for i, line in enumerate(lines):
            if line.startswith("MOL"):
                current_molecule = line.strip()
                data[current_molecule] = None

            elif line.startswith(f">  <{field}>"):
                data[current_molecule] = lines[i + 1].strip()

    return data


def get_data_fast_mode(file_path: str, field: str) -> dict[str, str]:
    data = {}

    with open(file_path, "r") as file:
        lines = file.readlines()

        for i, line in enumerate(lines):
            if line.startswith("id"):
                current_molecule = line.strip()
                current_title = lines[i - 1].strip()
                data[current_molecule] = None

            elif line.startswith(f">  <{field}>"):
                data[current_molecule] = (current_title, lines[i + 1].strip())

    return data
