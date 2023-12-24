# Standard Library
import logging
from io import StringIO
from typing import Any, TextIO

logger = logging.getLogger("SDParser")


class FastSDMol:
    def __init__(self, lines: list[str], data: dict[str, str]) -> None:
        self.lines = lines
        self.data = data

    @classmethod
    def read(cls, source: TextIO) -> "FastSDMol | None":
        lines: list[str] = []
        data: dict[str, str] = {}

        for line in source:
            if line.startswith("$$$$"):
                break
            if not line.startswith(">"):
                lines.append(line)
                continue

            # dealing with fields
            field_name = cls.parse_field_name(line)
            field_value = source.readline()
            if field_value.startswith("$$$$"):
                logger.warning(
                    f"found end of molecule {lines[0]} while looking for field {field_name} value."
                    " defaulting to empty string."
                )
                data[field_name] = ""
                break
            data[field_name] = field_value.strip("\n")
            discard_line = source.readline()
            if discard_line.startswith("$$$$"):
                logger.warning(f"found end of molecule {lines[0]} while expecting empty line after field {field_name}")
                break

        return cls(lines, data) if len(lines) >= 4 else None

    @staticmethod
    def parse_field_name(field_line: str) -> str:
        field_start = field_line.find("<") + 1
        field_end = field_line.find(">", 1)
        return field_line[field_start:field_end]

    @staticmethod
    def str_field(field_name: str, field_value: Any) -> str:
        return f">  <{field_name}>\n{field_value}\n\n"

    def __repr__(self) -> str:
        str_io = StringIO()
        self.write(str_io)
        return str_io.getvalue()

    def __str__(self) -> str:
        return f"<Molecule {self.lines[0]}>"

    def write(self, dest: TextIO) -> None:
        dest.writelines(self.lines)
        for field_name, field_value in self.data.items():
            dest.write(self.str_field(field_name, field_value))
        dest.write("$$$$")
