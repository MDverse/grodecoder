"""Defines type aliases for the project."""

from pathlib import Path
import MDAnalysis as mda

from MDAnalysis.core.groups import AtomGroup, Residue
from MDAnalysis.core.universe import Universe

PathLike = str | Path
UniverseLike = Universe | AtomGroup
Json = dict[str, "Json"] | list["Json"] | str | int | float | bool | None


