from datetime import datetime

import MDAnalysis as mda
from loguru import logger

from ._typing import PathLike, UniverseLike
from .identifier import identify
from .models import Decoded
from .toputils import guess_resolution


def _now() -> str:
    """Returns the current date and time formatted string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _read_structure(structure_path: PathLike, coordinates_path: PathLike | None = None) -> mda.Universe:
    """Reads a structure file."""
    if coordinates_path:
        return mda.Universe(structure_path, coordinates_path)
    return mda.Universe(structure_path)


def decode(universe: UniverseLike, bond_threshold: float = 5.0) -> Decoded:
    """Decodes the universe into an inventory of segments."""
    return Decoded(
        inventory=identify(universe, bond_threshold=bond_threshold),
        resolution=guess_resolution(universe),
    )


def decode_structure(
    structure_path: PathLike, coordinates_path: PathLike | None = None,
    bond_threshold: float = 5.0
) -> Decoded:
    """Reads a structure file and decodes it into an inventory of segments."""
    universe = _read_structure(structure_path, coordinates_path)
    assert universe.atoms is not None  # required by type checker for some reason
    logger.info(f"{structure_path}: {len(universe.atoms):,d} atoms")
    return decode(universe, bond_threshold=bond_threshold)

