import json
from datetime import datetime

import MDAnalysis as mda
from loguru import logger

from . import databases, toputils
from ._typing import AtomGroup, Json, PathLike, Residue, Universe, UniverseLike
from .identifier import identify
from .models import Decoded, GrodecoderRunOutput, GrodecoderRunOutputRead

__all__ = [
    "databases",
    "identify",
    "toputils",
    "read_structure",
    "AtomGroup",
    "Decoded",
    "GrodecoderRunOutput",
    "GrodecoderRunOutputRead",
    "Json",
    "PathLike",
    "Residue",
    "Universe",
    "UniverseLike",
]

__version__ = "0.0.1"


def _now() -> str:
    """Returns the current date and time formatted string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_structure(path: PathLike, psf_path: PathLike | None = None) -> Universe:
    """Reads a structure file."""
    if psf_path:
        return mda.Universe(path, psf_path)
    return mda.Universe(path)


def read_grodecoder_output(path: PathLike) -> GrodecoderRunOutputRead:
    with open(path) as fileobj:
        return GrodecoderRunOutputRead.model_validate(json.load(fileobj))


def decode(universe: UniverseLike, bond_threshold: float = 5.0) -> Decoded:
    """Decodes the universe into an inventory of segments."""
    return Decoded(
        inventory=identify(universe, bond_threshold=bond_threshold),
        resolution=toputils.guess_resolution(universe),
    )


def decode_structure(
    path: PathLike, psf_path: PathLike | None = None, bond_threshold: float = 5.0
) -> Decoded:
    """Reads a structure file and decodes it into an inventory of segments."""
    universe = read_structure(path, psf_path)
    assert universe.atoms is not None  # required by type checker for some reason
    logger.info(f"{path}: {len(universe.atoms):,d} atoms")
    return decode(universe, bond_threshold=bond_threshold)
