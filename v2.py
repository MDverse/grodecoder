from dataclasses import dataclass
from typing import Iterable

import numpy as np

import grodecoder as gd
import grodecoder.v1 as v1
import grodecoder.databases as DB

# DEBUG
import icecream
import time

icecream.install()


# Import the "old" grodecoder for comparison.
# To be removed.
v1.logger.remove()


def select_protein(universe: gd.UniverseLike) -> gd.AtomGroup:
    residue_names = DB.get_amino_acid_names()
    return universe.select_atoms(f"resname {' '.join(residue_names)}")


def has_bonds(residue1: gd.Residue, residue2: gd.Residue, cutoff: float = 5.0):
    """Returns True if the two residues are bonded."""
    # Check if the two residues are bonded
    cutoff_squared = cutoff**2
    distances_squared = (
        np.linalg.norm(residue1.atoms.positions[:, None] - residue2.atoms.positions, axis=-1) ** 2
    )
    bonded = distances_squared < cutoff_squared
    return True if np.any(bonded) else False


def detect_chains(universe: gd.UniverseLike, cutoff: float = 5.0) -> list[tuple[int, int]]:
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

    >>> from grodecoder import read_topology
    >>> universe = read_topology("3EAM.pdb")
    >>> protein = universe.select_atoms("protein")
    >>> chains = detect_chains(protein)
    >>> for start, end in chains:
    ...     print(protein.residues[start], protein.residues[end])
    <Residue VAL, 5> <Residue PHE, 315>
    <Residue VAL, 5> <Residue PHE, 315>
    <Residue VAL, 5> <Residue PHE, 315>
    <Residue VAL, 5> <Residue PHE, 315>
    """

    def is_last_residue():
        return i == len(universe.residues) - 1

    def end_of_chain():
        return is_last_residue() or not has_bonds(universe.residues[i], current_residue)

    segments = []
    start = 0
    for i, current_residue in enumerate(universe.residues[1:]):
        if end_of_chain():
            segments.append((start, i))
            start = i + 1
    return segments


@dataclass
class ResidueCount:
    residue_name: str
    atom_names: list[str]
    description: str
    number_of_atoms: int = 0
    number_of_residues: int = 0

    @classmethod
    def zero(cls, residue_name: str, atom_names: str, description: str):
        """Creates a ResidueCount object with zero atoms and residues."""
        return cls(residue_name, atom_names, description, 0, 0)


def count_residue(universe: gd.UniverseLike, definition: DB.Ion | DB.Solvent) -> ResidueCount:
    """Counts the number of residues and atoms matching a definition in a universe.

    This function is used to count the number of solvent and ion residues, for which the residue name
    and atom names are used to identify the residues.

    Parameters
    ----------
    universe : UniverseLike
        The universe to analyze. Typically the whole universe including ions and solvent residues.

    definition : Ion | Solvent
        The definition of the residue to count.

    Returns
    -------
    ResidueCount: description and actual count of the residue.
    """

    def are_equal_atom_names(lhs: np.ndarray[str], rhs: np.ndarray[str]) -> bool:
        """Checks if two lists of atom names are equal."""
        if len(lhs) != len(rhs):
            return False
        return (lhs == rhs).all()

    residue_name = definition.residue_name
    selection = universe.select_atoms(f"resname {residue_name} and name {' '.join(definition.atom_names)}")

    if len(selection) == 0:
        return ResidueCount.zero(residue_name, definition.atom_names, definition.description)

    # In the database, two definitions can have the same residue name, but different atom names
    # (e.g. TIP3 with atoms [OH2] and TIP3 with atoms [OH2, H1, H2]).
    # In this example, the first definition is a subset of the second one.
    # Therefore, we need to check that the selected residue has the same atom names as the definition.
    if not are_equal_atom_names(selection.residues[0].atoms.names, definition.atom_names):
        return ResidueCount.zero(residue_name, definition.atom_names, definition.description)

    return ResidueCount(
        residue_name,
        definition.atom_names,
        definition.description,
        len(selection.atoms),
        len(selection.residues),
    )


def count_solvents_and_ions(universe: gd.UniverseLike):
    ion_count = [count_residue(universe, ion) for ion in DB.get_ion_definitions()]
    solvent_count = [count_residue(universe, solvent) for solvent in DB.get_solvent_definitions()]
    # Filters out counts with zero atoms.
    return [count for count in ion_count + solvent_count if count.number_of_residues > 0]


def count_lipids(universe: gd.UniverseLike):
    lipid_definitions = DB.get_lipid_definitions()

    sel = universe.select_atoms("resname ADG")
    ic(sel)


    exit()
    counts = []
    for lipid in lipid_definitions:
        selection = universe.select_atoms(f"resname {lipid.residue_name}")
        if len(selection) == 0:
            continue

        counts.append(
            ResidueCount(
                residue_name=lipid.residue_name,
                description=lipid.description,
                number_of_atoms=len(selection.atoms),
                number_of_residues=len(selection.residues),
                atom_names=[],
            )
        )

    ic(counts)
    return counts


def main():
    # topology_path = "grodecoder/data/examples/1BRS.gro"
    # topology_path = "grodecoder/data/examples/2MAT.pdb"
    topology_path = "grodecoder/data/examples/4ZRY.gro"
    universe = gd.read_topology(topology_path)

    count_lipids(universe)

    # start = time.perf_counter()
    # result = count_solvents_and_ions(universe)
    # print(f"solvent and ion count: {time.perf_counter() - start:.2f} seconds")
    # ic(result)

    # protein = select_protein(universe)
    # start = time.perf_counter()
    # detect_chains(protein)
    # print(f"{time.perf_counter() - start:.2f} seconds")
    # exit()
    #
    # start = time.perf_counter()
    # old_fashion(protein)
    # print(f"{time.perf_counter() - start:.2f} seconds")
    #
    # exit()
    #


if __name__ == "__main__":
    main()
