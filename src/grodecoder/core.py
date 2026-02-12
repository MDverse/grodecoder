from datetime import datetime

from loguru import logger

from ._typing import PathLike, UniverseLike
from .identifier import identify
from .io import read_universe
from .models import Decoded
from .guesser import guess_resolution
from .settings import get_settings


def _now() -> str:
    """Returns the current date and time formatted string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def decode(universe: UniverseLike) -> Decoded:
    """Decodes the universe into an inventory of segments."""

    settings = get_settings()

    resolution = guess_resolution(universe)
    logger.info(f"Guessed resolution: {resolution}")

    # Guesses the chain dection distance cutoff if not provided by the user.
    chain_detection_settings = get_settings().chain_detection

    if chain_detection_settings.distance_cutoff.is_set():
        value = chain_detection_settings.distance_cutoff.get()
        logger.debug(f"chain detection: using user-defined value: {value:.2f}")
    else:
        logger.debug("chain detection: guessing distance cutoff based on resolution")
        chain_detection_settings.distance_cutoff.guess(resolution)

    distance_cutoff = chain_detection_settings.distance_cutoff.get()

    return Decoded(
        inventory=identify(universe, bond_threshold=distance_cutoff),
        resolution=resolution,
    )


def decode_structure(structure_path: PathLike, coordinates_path: PathLike | None = None) -> Decoded:
    """Reads a structure file and decodes it into an inventory of segments."""
    universe = read_universe(structure_path, coordinates_path)
    assert universe.atoms is not None  # required by type checker for some reason
    logger.debug(f"Universe has {len(universe.atoms):,d} atoms")
    return decode(universe)
