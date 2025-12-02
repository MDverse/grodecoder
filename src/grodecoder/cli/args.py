import sys
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from loguru import logger

DEFAULT_OUTPUT_STEM_SUFFIX = "_grodecoder"


def _fatal_error(msg: str, status: int = 1):
    """Prints an error message and exits with status `status`."""
    logger.critical(msg)
    sys.exit(status)


@dataclass
class InputFile:
    path: Path
    valid_extensions: ClassVar[set[str]]

    def __post_init__(self):
        # Ensures paths are pathlib.Path instances.
        self.path = Path(self.path)

        # Ensures paths are valid files.
        path = self.path
        if not path.exists():
            _fatal_error(f"'{path}' does not exist")
        if not path.is_file():
            _fatal_error(f"'{path}' is not a file")
        if path.suffix not in self.valid_extensions:
            _fatal_error(f"'{path}' has an invalid extension (valid extensions are {self.valid_extensions})")
        return path

    @property
    def extension(self) -> str:
        return self.path.suffix

    @property
    def stem(self) -> str:
        return self.path.stem


@dataclass
class StructureFile(InputFile):
    valid_extensions: ClassVar[set[str]] = {".gro", ".pdb", ".tpr", ".psf"}


@dataclass
class CoordinatesFile(InputFile):
    valid_extensions: ClassVar[set[str]] = {".gro", ".pdb", ".tpr", ".psf", ".coor"}


@dataclass
class Arguments:
    """Holds command-line arguments.

    Attrs:
        structure_file (Path): Path to the structure file.
        coordinates_file (Path): Path to the coordinates file.
        bond_threshold (float | None): Threshold for interchain bond detection.
        no_atom_ids (bool): If True, use compact serialization (no atom indices).
        print_to_stdout (bool): Whether to output results to stdout.
    """

    structure_file: StructureFile
    coordinates_file: CoordinatesFile | None = None
    bond_threshold: float | None = None
    no_atom_ids: bool = True
    print_to_stdout: bool = False

    def get_log_filename(self) -> Path:
        return generate_output_log_path(self.structure_file.stem)

    def get_inventory_filename(self) -> Path:
        return generate_output_inventory_path(self.structure_file.stem)


def generate_output_inventory_path(stem: str) -> Path:
    return Path(stem + DEFAULT_OUTPUT_STEM_SUFFIX + ".json")


def generate_output_log_path(stem: str) -> Path:
    return Path(stem + DEFAULT_OUTPUT_STEM_SUFFIX + ".log")
