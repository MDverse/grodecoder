import MDAnalysis as mda

from ._typing import PathLike, UniverseLike
from . import databases


__all__ = ["databases", "read_topology", "number_of_atoms", "number_of_residues"]



def read_topology(path: PathLike) -> mda.Universe:
    """Reads a topology file."""
    return mda.Universe(path)


def number_of_atoms(universe: UniverseLike) -> int:
    """Returns the number of atoms in the universe."""
    return len(universe.atoms)


def number_of_residues(universe: UniverseLike) -> int:
    """Returns the number of residues in the universe."""
    return len(universe.residues)
