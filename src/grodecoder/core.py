from datetime import datetime

from loguru import logger

from ._typing import PathLike, UniverseLike
from .identifier import identify
from .io import read_universe
from .models import Decoded
from .toputils import guess_resolution
from .settings import Settings


def _now() -> str:
    """Returns the current date and time formatted string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def decode(universe: UniverseLike, bond_threshold: float = 5.0) -> Decoded:
    """Decodes the universe into an inventory of segments."""
    return Decoded(
        inventory=identify(universe, bond_threshold=bond_threshold),
        resolution=guess_resolution(universe),
    )


def decode_structure(
    structure_path: PathLike,
    settings: Settings,
    coordinates_path: PathLike | None = None,
) -> Decoded:
    """Reads a structure file and decodes it into an inventory of segments."""
    universe = read_universe(structure_path, coordinates_path)
    assert universe.atoms is not None  # required by type checker for some reason
    logger.debug(f"Universe has {len(universe.atoms):,d} atoms")

    cutoff = settings.chain_detection.distance_cutoff or settings.chain_detection.default_distance_cutoff
    return decode(universe, bond_threshold=cutoff)
