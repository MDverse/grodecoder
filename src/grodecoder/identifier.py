import time
from typing import Iterable, Iterator

from loguru import logger
from MDAnalysis import AtomGroup

from . import databases as DB
from . import toputils
from ._typing import UniverseLike
from .models import HasAtoms, Inventory, MolecularType, Segment, SmallMolecule


def identify_small_molecule(
    universe: UniverseLike,
    definitions: Iterable[DB.ResidueDefinition],
    molecular_type: MolecularType = MolecularType.UNKNOWN,
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
        residue = SmallMolecule(
            atoms=selection, description=definition.description, molecular_type=molecular_type
        )
        logger.debug(
            f"identified small molecule {residue.description}: {len(selection.residues)} residues, {len(selection.atoms)} atoms"
        )
        counts.append(residue)
    return counts


def identify(universe: UniverseLike, bond_threshold: float = 5.0) -> Inventory:
    """Identifies the molecules in a structure file."""
    assert universe.atoms is not None  # avoids ty error "Expect Sized"
    timer_start = time.perf_counter()  # do not include structure reading time in the performance measurement
    inventory = _identify(universe, bond_threshold)
    elapsed = time.perf_counter() - timer_start
    logger.info(f"{len(universe.atoms):,d} atoms processed in {elapsed:.2f} seconds")
    return inventory


# ========================================================================================================
#
#   Private functions
#
# ========================================================================================================


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
    logger.debug("excluding possible methanol residues (MET) from protein")
    methanol = _find_methanol(universe)
    if methanol:
        indexes = " ".join(str(i) for i in methanol)  # ty complains when using map
        selection_str += f" and not index {indexes}"

    logger.debug("selecting protein")
    protein = universe.select_atoms(selection_str)
    logger.debug("selecting protein - done")
    return protein


def _select_nucleic(universe: UniverseLike) -> AtomGroup:
    """Selects the nucleic acid atoms from the universe."""
    nucleic_acid_residue_names = DB.get_nucleotide_names()
    selection_str = f"resname {' '.join(nucleic_acid_residue_names)}"
    logger.debug("selecting nucleic")
    nucleic = universe.select_atoms(selection_str)
    logger.debug("selecting nucleic - done")
    return nucleic


def _iter_chains(atoms: AtomGroup, bond_threshold: float = 5.0) -> Iterator[AtomGroup]:
    """Iterates over the chains of a group of atoms.

    Chains are defined as segments of residues that are connected by bonds.
    """
    if len(atoms) == 0:
        return
    logger.debug(f"detecting segments using cutoff distance {bond_threshold:.2f}")
    segments = toputils.detect_chains(atoms, cutoff_distance=bond_threshold)

    n_seg_str = f"{len(segments)} segment" + "s" if len(segments) > 1 else ""
    logger.debug(f"detecting segments - done: found {n_seg_str}")

    for start, end in segments:
        logger.debug(f"yielding segment containing residues {start} to {end}")
        yield atoms.residues[start : end + 1].atoms


def _get_protein_segments(atoms: AtomGroup, bond_threshold: float = 5.0) -> list[Segment]:
    """Returns the protein segments in the universe."""
    protein = _select_protein(atoms)
    return [
        Segment(atoms=atoms, molecular_type=MolecularType.PROTEIN)
        for atoms in _iter_chains(protein, bond_threshold)
    ]


def _get_nucleic_segments(atoms: AtomGroup, bond_threshold: float = 5.0) -> list[Segment]:
    """Returns the nucleic acid segments in the universe."""
    nucleic = _select_nucleic(atoms)
    return [
        Segment(atoms=atoms, molecular_type=MolecularType.NUCLEIC)
        for atoms in _iter_chains(nucleic, bond_threshold)
    ]


def _unique_definitions(definitions: Iterable[DB.ResidueDefinition]) -> dict[str, DB.ResidueDefinition]:
    """Returns a dictionary of unique residue definitions by their residue name."""
    unique_items = {}
    for definition in definitions:
        if definition.residue_name not in unique_items:
            unique_items[definition.residue_name] = definition
    return unique_items


def _remove_identified_atoms(universe: AtomGroup, molecules: list[HasAtoms]) -> AtomGroup:
    """Removes the atoms of the identified molecules from the universe."""
    for molecule in molecules:
        universe -= molecule.atoms
    return universe


def _trim_sequence(sequence: str) -> str:
    """Trims the sequence to a maximum length for logging."""
    max_length = 10
    if len(sequence) > max_length:
        return sequence[: max_length // 2] + "..." + sequence[-max_length // 2 :]
    return sequence


def _log_identified_segments(segments: list[Segment], label: str) -> None:
    """Logs the identified segments."""
    msg = f"Identified {len(segments)} {label} segments{':' if len(segments) > 0 else ''}"
    logger.info(msg)
    for seg in segments:
        logger.info(
            f"  - Segment with {seg.number_of_residues:,d} residues {seg.number_of_atoms:,d} "
            f"atoms and sequence: {_trim_sequence(seg.sequence)}"
        )


def _log_identified_molecules(molecules: list[SmallMolecule], label: str) -> None:
    """Logs the identified small molecules."""
    logger.info(f"Identified {len(molecules)} {label} types{':' if len(molecules) > 0 else ''}")
    for mol in molecules:
        logger.info(
            f"  - {label} {mol.name!r}, {mol.number_of_atoms:,d} atoms, identified as {mol.description!r}"
        )


def _identify(universe: UniverseLike, bond_threshold: float = 5.0) -> Inventory:
    """Identifies the molecules in the universe."""
    logger.debug("Residu identification: start")

    # Ensure the universe is an AtomGroup.
    universe = universe.select_atoms("all")
    total_number_of_atoms = len(universe)

    # Remove identified atoms from the universe along the way to avoid double counting (e.g.
    # 'MET' residues are counted first in the protein, then removed so not counted elsewhere).

    # All ty: ignore[invalid-argument-type] in this block fix ty clear mistake:
    # Expected `list[HasAtoms]`, found `list[Segment]`
    # while `Segment` clearly satisfies the `HasAtoms` Protocol.

    protein = _get_protein_segments(universe, bond_threshold=bond_threshold)
    _log_identified_segments(protein, "protein")
    universe = _remove_identified_atoms(universe, protein)  # ty: ignore[invalid-argument-type]

    nucleic = _get_nucleic_segments(universe, bond_threshold=bond_threshold)
    _log_identified_segments(nucleic, "nucleic acid")
    universe = _remove_identified_atoms(universe, nucleic)  # ty: ignore[invalid-argument-type]

    ions = identify_small_molecule(universe, DB.get_ion_definitions(), molecular_type=MolecularType.ION)
    _log_identified_molecules(ions, "ion")
    universe = _remove_identified_atoms(universe, ions)  # ty: ignore[invalid-argument-type]

    solvents = identify_small_molecule(
        universe, DB.get_solvent_definitions(), molecular_type=MolecularType.SOLVENT
    )
    _log_identified_molecules(solvents, "solvent")
    universe = _remove_identified_atoms(universe, solvents)  # ty: ignore[invalid-argument-type]

    lipids = identify_small_molecule(universe, DB.get_lipid_definitions(), molecular_type=MolecularType.LIPID)
    _log_identified_molecules(lipids, "lipid")
    universe = _remove_identified_atoms(universe, lipids)  # ty: ignore[invalid-argument-type]

    others = identify_small_molecule(
        universe, DB.get_other_definitions(), molecular_type=MolecularType.UNKNOWN
    )
    _log_identified_molecules(others, "other small molecule")
    universe = _remove_identified_atoms(universe, others)  # ty: ignore[invalid-argument-type]

    # Collecting unknown small molecules that are not defined in the database.
    unknown_molecules = []
    if len(universe.atoms) > 0:
        # Get the unique residue names from the remaining atoms in the universe.
        # Not sorting residue names would result in random order in output file.
        residue_names = sorted(set(universe.residues.resnames))
        for residue_name in residue_names:
            logger.warning(f"residue {residue_name!r} is unidentified, counting as unknown molecule")
            molecule = SmallMolecule(
                atoms=universe.select_atoms(f"resname {residue_name}"),
                description="Unknown small molecule",
                molecular_type=MolecularType.UNKNOWN,
            )
            unknown_molecules.append(molecule)

    logger.debug("Creating inventory")
    inventory = Inventory(
        segments=protein + nucleic,
        small_molecules=ions + solvents + lipids + others + unknown_molecules,
        total_number_of_atoms=total_number_of_atoms,
    )

    logger.debug("Residu identification: end")
    return inventory
