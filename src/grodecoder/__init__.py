import hashlib
from datetime import datetime

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

__version__ = "0.0.1"


def _now()  -> str:
    """Returns the current date and time as an ISO 8601 formatted string."""
    return datetime.now().isoformat()


def read_topology(path: PathLike, psf_path: PathLike | None = None) -> Universe:
    """Reads a topology file."""
    if psf_path:
        return mda.Universe(path, psf_path)
    return mda.Universe(path)


def decode(universe: UniverseLike, topology_checksum: str, bond_threshold: float = 5.0) -> Decoded:
    """Decodes the universe into an inventory of segments."""
    return Decoded(
        inventory=identify(universe, bond_threshold=bond_threshold),
        resolution=toputils.guess_resolution(universe),
        topology_checksum=topology_checksum,
        database_version=databases.__version__,
        grodecoder_version=__version__,
        grodecoder_run_date=_now(),
    )


def _get_checksum(topology_path: PathLike) -> str:
    """Computes a checksum for the topology file."""
    with open(topology_path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def decode_topology(path: PathLike, psf_path: PathLike | None = None, bond_threshold: float = 5.0) -> Decoded:
    """Reads a topology file and decodes it into an inventory of segments."""
    universe = read_topology(path, psf_path)
    logger.debug(f"{path}: {len(universe.atoms):,d} atoms")
    return decode(universe, topology_checksum=_get_checksum(path), bond_threshold=bond_threshold)
