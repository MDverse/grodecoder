import MDAnalysis as mda

from ._typing import Json, PathLike, Universe, UniverseLike, AtomGroup, Residue
from . import databases


__all__ = ["databases", "read_topology", "number_of_atoms", "number_of_residues",
           "Json", "PathLike", "Universe", "UniverseLike", "AtomGroup", "Residue"]


def read_topology(path: PathLike, psf_path: PathLike | None = None) -> Universe:
    """Reads a topology file."""
    if psf_path:
        return mda.Universe(path, psf_path)
    return mda.Universe(path)


def number_of_atoms(universe: UniverseLike) -> int:
    """Returns the number of atoms in the universe."""
    return len(universe.atoms)


def number_of_residues(universe: UniverseLike) -> int:
    """Returns the number of residues in the universe."""
    return len(universe.residues)
