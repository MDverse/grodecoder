"""Defines utility functions for working with molecular structures."""

import collections
from typing import Iterable

import numpy as np

from ._typing import Residue, UniverseLike
from .databases import get_amino_acid_name_map, get_nucleotide_name_map
from .models import MolecularResolution


def _first_alpha(s: str) -> str:
    """Returns the first alphabetical character in a string."""
    for char in s:
        if char.isalpha():
            return char
    return ""


def formula(molecule: UniverseLike, elements: Iterable[str] = ("C", "H", "O", "N", "P")) -> str:
    """Returns the formula of the lipid.

    The formula is based on the first residue alone.
    """
    residues = getattr(molecule.atoms, "residues", [])
    assert len(residues) > 0
    atoms = residues[0].atoms
    count = collections.Counter(_first_alpha(atom.name) for atom in atoms)
    formula = "".join(
        f"{element}{count[element]}" if count[element] > 1 else element
        for element in elements
        if count[element] > 0
    )
    return formula


def sequence(atoms: UniverseLike) -> str:
    """Returns the protein sequence from a set of atoms."""
    residue_names = get_amino_acid_name_map()
    residue_names.update(get_nucleotide_name_map())
    return "".join(residue_names.get(residue.resname, "X") for residue in getattr(atoms, "residues", []))


def has_bonds(residue: Residue, cutoff: float = 2.0):
    """Returns True if the residue has any bonds."""
    cutoff_squared = cutoff**2
    distances_squared = (
        np.linalg.norm(residue.atoms.positions[:, None] - residue.atoms.positions, axis=-1) ** 2
    )
    np.fill_diagonal(distances_squared, np.inf)  # Ignore self-pairs
    bonded = distances_squared < cutoff_squared
    return bool(np.any(bonded))


def has_bonds_between(residue1: Residue, residue2: Residue, cutoff: float = 5.0):
    """Returns True if the two residues are bonded."""
    cutoff_squared = cutoff**2
    distances_squared = (
        np.linalg.norm(residue1.atoms.positions[:, None] - residue2.atoms.positions, axis=-1) ** 2
    )
    bonded = distances_squared < cutoff_squared
    return bool(np.any(bonded))


def detect_chains(universe: UniverseLike, cutoff: float = 5.0) -> list[tuple[int, int]]:
    """Detects chains in a set of atoms.

    Iterates over the residues as detected by MDAnalysis and calculates the bonds
    between them. If the distance between two residues is larger than `cutoff`, they are
    considered to be in different chains.

    Parameters
    ----------
    universe : AtomGroup
        The universe to analyze. Typically a protein or a set of residues.

    cutoff : float, optional
        The cutoff distance to determine if two residues are bonded. Default is 5.0.

    Returns
    --------
    list[tuple[int, int]]: pairs of indices of the residues that are in the same chain.

    Example
    -------

    >>> from grodecoder import read_structure
    >>> universe = read_structure("3EAM.pdb")
    >>> protein = universe.select_atoms("protein")
    >>> chains = detect_chains(protein)
    >>> for start, end in chains:
    ...     print(protein.residues[start], protein.residues[end])
    <Residue VAL, 5> <Residue PHE, 315>
    <Residue VAL, 5> <Residue PHE, 315>
    <Residue VAL, 5> <Residue PHE, 315>
    <Residue VAL, 5> <Residue PHE, 315>
    <Residue VAL, 5> <Residue PHE, 315>
    """

    def end_of_chain():
        """The end of a chain is defined as the point where two consecutive residues are not bonded."""
        return not has_bonds_between(current_residue, next_residue, cutoff)

    segments = []

    assert (residues := getattr(universe, "residues", [])) and len(residues) > 0

    start_current_chain = 0
    for i, current_residue in enumerate(residues[:-1]):
        next_residue = residues[i + 1]
        if end_of_chain():
            segments.append((start_current_chain, i))
            start_current_chain = i + 1
    segments.append((start_current_chain, len(residues) - 1))

    return segments


def guess_resolution(universe: UniverseLike) -> MolecularResolution:
    """Guesses the resolution (i.e. all-atom or coarse grain) of the universe.

    The resolution is considered coarse-grained if a residue has at least two atoms within a distance of 2.0 Å.

    Finds the first five residues with at least two atoms and checks if they have bonds.
    If any of them have bonds, the resolution is considered all-atom.
    If none of the first five residues have bonds, the resolution is considered coarse-grained.
    """
    # Select the first five residues with at least two atoms.
    assert (residues := getattr(universe, "residues", [])) and len(residues) > 0
    residues = [residue for residue in residues if len(residue.atoms) >= 2][:5]
    for residue in residues:
        if has_bonds(residue, cutoff=2.0):
            return MolecularResolution.ALL_ATOM
    return MolecularResolution.COARSE_GRAINED
