from typing import Iterable, Iterator

from loguru import logger
from MDAnalysis import AtomGroup

from . import databases as DB
from . import toputils
from ._typing import UniverseLike
from .models import Inventory, SmallMolecule, Segment


def _find_methanol(universe: UniverseLike) -> list[int]:
    """Returns the indices of methanol atoms in the universe."""
    met = universe.select_atoms("resname MET")
    methanol = []
    for residue in met.residues:
        carbones = residue.atoms.select_atoms("name C*")
        if len(carbones) == 1:
            methanol += residue.atoms.indices.tolist()
    return methanol


def _select_protein(universe: UniverseLike) -> AtomGroup:
    """Selects the protein atoms from the universe."""
    protein_residue_names = DB.get_amino_acid_names()
    selection_str = f"resname {' '.join(protein_residue_names)}"

    # Exclude methanol residues from the selection.
    methanol = _find_methanol(universe)
    if methanol:
        selection_str += f" and not index {' '.join(map(str, methanol))}"

    return universe.select_atoms(selection_str)


def _select_nucleic(universe: UniverseLike) -> AtomGroup:
    """Selects the nucleic acid atoms from the universe."""
    nucleic_acid_residue_names = DB.get_nucleotide_names()
    selection_str = f"resname {' '.join(nucleic_acid_residue_names)}"
    return universe.select_atoms(selection_str)


def _iter_chains(atoms: AtomGroup) -> Iterator[AtomGroup]:
    """Iterates over the chains of a group of atoms.

    Chains are defined as segments of residues that are connected by bonds.
    """
    if len(atoms) == 0:
        return
    segments = toputils.detect_chains(atoms)
    for start, end in segments:
        yield atoms.residues[start : end + 1].atoms


def _get_protein_segments(atoms: AtomGroup) -> list[Segment]:
    """Returns the protein segments in the universe."""
    protein = _select_protein(atoms)
    return [
        Segment(atoms=atoms, sequence=toputils.sequence(atoms), molecular_type="protein")
        for atoms in _iter_chains(protein)
    ]


def _get_nucleic_segments(atoms: AtomGroup) -> list[Segment]:
    """Returns the nucleic acid segments in the universe."""
    nucleic = _select_nucleic(atoms)
    return [
        Segment(atoms=atoms, sequence=toputils.sequence(atoms), molecular_type="nucleic_acid")
        for atoms in _iter_chains(nucleic)
    ]


def _unique_definitions(definitions: Iterable[DB.ResidueDefinition]) -> dict[str, DB.ResidueDefinition]:
    """Returns a dictionary of unique residue definitions by their residue name."""
    unique_items = {}
    for definition in definitions:
        if definition.residue_name not in unique_items:
            unique_items[definition.residue_name] = definition
    return unique_items


def _remove_identified_atoms(universe: AtomGroup, molecules: list[AtomGroup]) -> AtomGroup:
    """Removes the atoms of the identified molecules from the universe."""
    for molecule in molecules:
        universe -= molecule.atoms
    return universe


def identify_small_molecule(
    universe: UniverseLike,
    definitions: Iterable[DB.ResidueDefinition],
    molecular_type: str = "small_molecule",
) -> list[SmallMolecule]:
    """Counts the residues in the universe based on the definitions.

    Only using the residue name, not the atom names.
    """
    counts = []
    u_definitions = _unique_definitions(definitions)
    for residue_name, definition in u_definitions.items():
        selection = universe.select_atoms(f"resname {definition.residue_name}")
        if len(selection) == 0:
            continue
        residue = SmallMolecule(selection, definition.description, molecular_type=molecular_type)
        counts.append(residue)
    return counts


def identify(universe: UniverseLike) -> Inventory:
    """Identifies the molecules in the universe."""

    # Ensure the universe is an AtomGroup.
    universe = universe.select_atoms("all")

    # Remove identified atoms from the universe along the way to avoid double counting (e.g.
    # 'MET' residues are counted first in the protein, then removed so not counted elsewhere).

    protein = _get_protein_segments(universe)
    universe = _remove_identified_atoms(universe, protein)

    nucleic = _get_nucleic_segments(universe)
    universe = _remove_identified_atoms(universe, nucleic)

    ions = identify_small_molecule(universe, DB.get_ion_definitions(), molecular_type="ion")
    universe = _remove_identified_atoms(universe, ions)

    solvents = identify_small_molecule(universe, DB.get_solvent_definitions(), molecular_type="solvent")
    universe = _remove_identified_atoms(universe, solvents)

    lipids = identify_small_molecule(universe, DB.get_lipid_definitions(), molecular_type="lipid")
    universe = _remove_identified_atoms(universe, lipids)

    others = identify_small_molecule(universe, DB.get_other_definitions(), molecular_type="other")
    universe = _remove_identified_atoms(universe, others)

    # Collecting unknown small molecules that are not defined in the database.
    unknown_molecules = []
    if len(universe.atoms) > 0:
        residue_names = set(universe.residues.resnames)
        for residue_name in residue_names:
            logger.warning(f"residue {residue_name!r} is unidentified, counting as unknown molecule")
            molecule = SmallMolecule(universe.select_atoms(f"resname {residue_name}"))
            unknown_molecules.append(molecule)

    return Inventory(
        segments=protein + nucleic,
        small_molecules=ions + solvents + lipids + others + unknown_molecules,
    )
