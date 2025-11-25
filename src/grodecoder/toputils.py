"""Defines utility functions for working with molecular structures."""

import collections
from itertools import islice
from typing import Iterable

import numpy as np
from loguru import logger

from ._typing import Residue, UniverseLike
from .logging import is_logging_debug
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


def has_bonds(residue: Residue, distance_cutoff: float = 2.0):
    """Returns True if the residue has any bonds."""
    return bool(np.any(get_bonds(residue, distance_cutoff)))


def get_bonds(residue: Residue, cutoff: float = 2.0):
    """Returns the bonds between the atoms of a residue."""
    cutoff_squared = cutoff**2
    distances_squared = (
        np.linalg.norm(residue.atoms.positions[:, None] - residue.atoms.positions, axis=-1) ** 2
    )

    # ignore self-pairs and permutations
    distances_squared[np.tril_indices(distances_squared.shape[0])] = np.inf

    return np.argwhere(distances_squared < cutoff_squared)


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

    >>> universe = MDAnalysis.Universe("path/to/structure_file.pdb")
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


def guess_resolution(universe: UniverseLike, cutoff_distance: float) -> MolecularResolution:
    """Guesses the resolution (i.e. all-atom or coarse grain) of the universe.

    The resolution is considered coarse-grained if a residue has at least two atoms within a distance of 
    `cutoff_distance` Å.

    Finds the first five residues with at least two atoms and checks if they have bonds.
    If any of them have bonds, the resolution is considered all-atom.
    If none of the first five residues have bonds, the resolution is considered coarse-grained.
    """

    def debug(msg):
        where = f"{__name__}.guess_resolution"
        logger.debug(f"{where}: {msg}")

    def print_bonds(residue):
        """Print bonds between atoms inside a residue. Used for debug purposes."""

        def distance(atom1, atom2):
            return (np.linalg.norm(atom1.position - atom2.position) ** 2.0) ** 0.5

        bonds = get_bonds(residue, cutoff_distance)
        for bond in bonds:
            left, right = residue.atoms[bond]
            bond_str = f"residue {left.resname}:{left.resid}, atoms {left.name}-{right.name}"
            debug(f"guess_resolution: Found bond: {bond_str} (distance={distance(left, right):.2f})")
        pair_str = f"pair{'s' if len(bonds) > 1 else ''}"
        debug(f"guess_resolution: detected {len(bonds)} {pair_str} with distance < {cutoff_distance=:.2f}")

    debug(f"start ; {cutoff_distance=:.2f}")

    # Makes ty happy.
    assert (residues := getattr(universe, "residues", [])) and len(residues) > 0

    # Selects the first five residues with at least two atoms.
    residues = list(islice((residue for residue in residues if len(residue.atoms) > 1), 10))

    for residue in residues:
        if has_bonds(residue, cutoff_distance):
            if is_logging_debug():
                print_bonds(residue)
            debug("end: detected resolution: ALL_ATOM")
            return MolecularResolution.ALL_ATOM
    debug(
        f"No intra-atomic distance within {cutoff_distance:.2f} Å found in the first {len(residues)} residues"
    )
    debug("end: detected resolution: COARSE_GRAINED")
    return MolecularResolution.COARSE_GRAINED
