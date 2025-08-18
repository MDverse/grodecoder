import MDAnalysis as mda
from loguru import logger

from . import databases, toputils
from ._typing import AtomGroup, Json, PathLike, Residue, Universe, UniverseLike
from .identifier import identify
from .models import Decoded

__all__ = [
    "databases",
    "identify",
    "toputils",
    "read_topology",
    "AtomGroup",
    "Decoded",
    "Json",
    "PathLike",
    "Residue",
    "Universe",
    "UniverseLike",
]


def read_topology(path: PathLike, psf_path: PathLike | None = None) -> Universe:
    """Reads a topology file."""
    if psf_path:
        return mda.Universe(path, psf_path)
    return mda.Universe(path)


def decode(universe: UniverseLike) -> Decoded:
    """Decodes the universe into an inventory of segments."""
    return Decoded(
        inventory=identify(universe),
        resolution=toputils.guess_resolution(universe),
    )


def decode_topology(path: PathLike, psf_path: PathLike | None = None) -> Decoded:
    """Reads a topology file and decodes it into an inventory of segments."""
    universe = read_topology(path, psf_path)
    logger.debug(f"{path}: {len(universe.atoms):,d} atoms")
    return decode(universe)
