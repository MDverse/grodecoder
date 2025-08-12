import MDAnalysis as mda

from ._typing import Json, PathLike, Universe, UniverseLike, AtomGroup, Residue
from . import databases
from . import toputils
from .identifier import identify


__all__ = [
    "databases",
    "read_topology",
    "number_of_atoms",
    "number_of_residues",
    "toputils",
    "identify",
    "Json",
    "PathLike",
    "Universe",
    "UniverseLike",
    "AtomGroup",
    "Residue",
]


def read_topology(path: PathLike, psf_path: PathLike | None = None) -> Universe:
    """Reads a topology file."""
    if psf_path:
        return mda.Universe(path, psf_path)
    return mda.Universe(path)
