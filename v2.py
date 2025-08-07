import time
from pathlib import Path
from typing import Iterable

import numpy as np

import grodecoder as gd
import grodecoder.databases as DB

# DEBUG
import icecream

icecream.install()


class HasAtoms:
    "Interface for objects that have atoms."

    def __init__(self, atoms: gd.AtomGroup):
        """Initializes the SmallMoleculeBase with a set of atoms."""
        self.atoms = atoms

    @property
    def atom_indices(self) -> np.ndarray:
        """Returns the indices of the atoms in the residue."""
        return self.atoms.indices

    @property
    def number_of_atoms(self) -> int:
        """Returns the number of atoms in the residue."""
        return len(self.atoms)

    @property
    def number_of_residues(self) -> int:
        """Returns the number of residues in the residue."""
        return len(self.atoms.residues)

    @property
    def residues(self):
        """Returns the residues in the residue."""
        return self.atoms.residues


class SmallMoleculeBase(HasAtoms):
    """Used to count small molecules in a universe."""


    def __repr__(self) -> str:
        """Returns a string representation of the SmallMolecule."""
        return f"<{self.__class__.__name__} {self.name}({self.number_of_atoms} atoms, {self.number_of_residues} residues)>"

    @property
    def name(self) -> str:
        """Returns the name of the residue."""
        return self.atoms.residues[0].resname

    @property
    def description(self) -> str:
        """Returns the description of the residue."""
        return "Unknown small molecule"

    def dump_json(self) -> dict:
        """Returns a JSON representation of the residue."""
        return {
            "name": self.name,
            "description": self.description,
            "number_of_atoms": self.number_of_atoms,
            "number_of_residues": self.number_of_residues,
            "molecular_type": getattr(self, "molecular_type", "unknown"),
            # "atom_indices": self.atom_indices.tolist(),
        }


class SmallMolecule(SmallMoleculeBase):
    """Used to count small molecules with actual definition in the database."""

    def __init__(self, definition: DB.ResidueDefinition, atoms: gd.AtomGroup):
        """Initializes the SmallMolecule with a definition and a set of atoms."""
        super().__init__(atoms)
        self.definition = definition

    def __repr__(self) -> str:
        """Returns a string representation of the SmallMolecule."""
        return f"<{self.__class__.__name__} {self.name}({self.number_of_atoms} atoms, {self.number_of_residues} residues)>"

    @property
    def name(self) -> str:
        """Returns the name of the residue."""
        return self.definition.residue_name

    @property
    def description(self) -> str:
        """Returns the description of the residue."""
        return self.definition.description


class Ion(SmallMolecule):
    """Ion residue."""

    molecular_type: str = "ion"


class Solvent(SmallMolecule):
    """Solvent residue."""

    molecular_type: str = "solvent"


class Lipid(SmallMolecule):
    """Lipid residue."""

    molecular_type: str = "lipid"

    def dump_json(self) -> dict:
        """Returns a JSON representation of the lipid."""
        data = super().dump_json()
        data["formula"] = self.formula
        return data

    @property
    def formula(self) -> str:
        """Returns the formula of the lipid."""
        atoms = self.atoms.residues[0].atoms  # formula is based on the first residue alone
        elements = [first_alpha(atom.name) for atom in atoms]
        formula = ""
        for element in ("C", "H", "O", "N", "P"):
            count = elements.count(element)
            if count > 1:
                formula += f"{element}{count}"
            elif count > 0:
                formula += f"{element}"
        return formula


class ProteinSegment(HasAtoms):
    """Holds a segment of a protein."""

    def __init__(self, atoms: gd.AtomGroup):
        """Initializes the ProteinSegment with a set of atoms."""
        self.atoms = atoms

    def __repr__(self) -> str:
        return (
            f"<{self.__class__.__name__}({self.number_of_atoms} atoms, {self.number_of_residues} residues)>"
        )

    @property
    def sequence(self) -> str:
        """Returns the sequence of the protein segment."""
        return sequence(self.atoms)


class Protein(HasAtoms):
    """Holds protein information."""

    def __init__(self, atoms: gd.AtomGroup):
        """Initializes the Protein with a set of atoms."""
        self.atoms = atoms
        self._segments = detect_chains(self.atoms)

    def __repr__(self) -> str:
        number_of_segments = len(self._segments)
        return f"<{self.__class__.__name__}({self.number_of_atoms} atoms, {self.number_of_residues} residues, {number_of_segments} segments)>"

    def iter_segments(self):
        """Iterates over the segments of the protein."""
        for start, end in self._segments:
            yield ProteinSegment(self.atoms.atoms.residues[start : end + 1].atoms)


def first_alpha(s: str) -> str:
    """Returns the first alphabetical character in a string."""
    for char in s:
        if char.isalpha():
            return char
    return ""


def sequence(atoms: gd.AtomGroup) -> str:
    """Returns the protein sequence from a set of atoms."""
    seq = ""
    names = {residue.long_name: residue.short_name for residue in DB.get_amino_acid_definitions()}
    for residue in atoms.residues:
        seq += names.get(residue.resname, "X")
    return seq


def find_methanol(universe: gd.UniverseLike) -> list[int]:
    """Returns the indices of methanol atoms in the universe."""
    met = universe.select_atoms("resname MET")
    methanol = []
    for residue in met.residues:
        carbones = residue.atoms.select_atoms("name C*")
        if len(carbones) == 1:
            methanol += residue.atoms.indices.tolist()
    return methanol


def select_protein(universe: gd.UniverseLike) -> Protein:
    """Selects the protein atoms from the universe."""
    protein_residue_names = DB.get_amino_acid_names()
    selection_str = f"resname {' '.join(protein_residue_names)}"

    # Exclude methanol residues from the selection.
    methanol = find_methanol(universe)
    if methanol:
        selection_str += f" and not index {' '.join(map(str, methanol))}"

    return Protein(universe.select_atoms(selection_str))


def has_bonds(residue1: gd.Residue, residue2: gd.Residue, cutoff: float = 5.0):
    """Returns True if the two residues are bonded."""
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
    <Residue VAL, 5> <Residue PHE, 315>
    """

    def end_of_chain():
        return not has_bonds(current_residue, next_residue, cutoff)

    segments = []

    start_current_chain = 0
    for i, current_residue in enumerate(universe.residues[:-1]):
        next_residue = universe.residues[i + 1]
        if end_of_chain():
            segments.append((start_current_chain, i))
            start_current_chain = i + 1
    segments.append((start_current_chain, len(universe.residues) - 1))

    return segments


def unique_definitions(definitions: Iterable[DB.ResidueDefinition]) -> dict[str, DB.ResidueDefinition]:
    """Returns a dictionary of unique residue definitions by their residue name."""
    unique_items = {}
    for definition in definitions:
        if definition.residue_name not in unique_items:
            unique_items[definition.residue_name] = definition
    return unique_items


def count(
    universe: gd.UniverseLike, definitions: Iterable[DB.ResidueDefinition], model: type = SmallMolecule
) -> list[SmallMolecule]:
    """Counts the residues in the universe based on the definitions.

    Only using the residue name, not the atom names.
    """
    counts = []
    u_definitions = unique_definitions(definitions)
    for residue_name, definition in u_definitions.items():
        selection = universe.select_atoms(f"resname {definition.residue_name}")
        if len(selection) == 0:
            continue
        residue = model(definition, selection)
        counts.append(residue)
    return counts


def count_ions(universe: gd.UniverseLike) -> list[Ion]:
    """Counts the ions in the universe."""
    return count(universe, DB.get_ion_definitions(), model=Ion)


def count_solvents(universe: gd.UniverseLike) -> list[Solvent]:
    """Counts the solvents in the universe."""
    return count(universe, DB.get_solvent_definitions(), model=Solvent)


def count_lipids(universe: gd.UniverseLike) -> list[Lipid]:
    """Counts the lipids in the universe."""
    return count(universe, DB.get_lipid_definitions(), model=Lipid)


def count_protein_segments(universe: gd.UniverseLike) -> dict["str", list[ProteinSegment]]:
    """Counts the protein segments in the universe.

    Segments are grouped by their sequence, hence the return type.
    """
    protein = select_protein(universe)

    if len(protein.atoms) == 0:
        return {}

    # Group the protein segments by sequence.
    by_sequence: dict[str, list[ProteinSegment]] = {}
    for segment in protein.iter_segments():
        by_sequence.setdefault(segment.sequence, []).append(segment)

    return by_sequence


def dump_protein_segments_json(by_sequence: dict[str, list[ProteinSegment]]) -> list[dict]:
    """Dumps the protein segments to JSON format.

    Segments are grouped by their sequence.
    """
    segments_json = []
    for sequence, segments in by_sequence.items():
        segments_json.append(
            {
                "sequence": sequence,
                "number_of_atoms": [segment.number_of_atoms for segment in segments],
                "number_of_residues": [segment.number_of_residues for segment in segments],
                "number_of_segments": len(segments),
                # "atom_indices": [segment.atom_indices.tolist() for segment in segments],
                "molecular_type": "protein",
            }
        )
    return segments_json


def actually_count(topology_path: Path):
    universe = gd.read_topology(topology_path).select_atoms(
        "all"
    )  # needs the 'select_atoms' to convert to AtomGroup

    protein_segments = count_protein_segments(universe)

    # Removing previously analyzed parts of the system prevents residues from being counted in several groups
    # (e.g. residue "MET" could be counted as part of a protein segment and as a solvent).
    for segments in protein_segments.values():
        for segment in segments:
            universe -= segment.atoms

    ions = count_ions(universe)
    for ion in ions:
        universe -= ion.atoms

    solvents = count_solvents(universe)
    for solvent in solvents:
        universe -= solvent.atoms

    lipids = count_lipids(universe)
    for lipid in lipids:
        universe -= lipid.atoms

    # Collecting unknown small molecules that are not defined in the database.
    unknown_molecules = []
    if len(universe.atoms) > 0:
        residue_names = set(universe.residues.resnames)
        for residue_name in residue_names:
            molecule = SmallMoleculeBase(universe.select_atoms(f"resname {residue_name}"))
            unknown_molecules.append(molecule)

    json_data = {
        "inventory": (
            [ion.dump_json() for ion in ions]
            + [solvent.dump_json() for solvent in solvents]
            + [lipid.dump_json() for lipid in lipids]
            + [molecule.dump_json() for molecule in unknown_molecules]
            + dump_protein_segments_json(protein_segments)
        )
    }

    return json_data


def main():
    data_root_dir = Path("data/examples")
    topology_files = [
        "1BRS.gro",
        "1QJ8.gro",
        "1QJ8_ETH_ACN_MET_URE_SOL.pdb",
        "1QJ8_membrane.gro",
        "1QJ8_solution.gro",
        "2MAT.pdb",
        "4MQJ_ABCD.gro",
        "4ZRY.gro",
        "5MBA.gro",
        "5ZOA.gro",
        "barstar.gro",
        "DMPC_PI.gro",
        "DNA_start.gro",
        "noriega_AA_CRD_3CAL.gro",
        "noriega_CG_CRD_3CAL.gro",
        "RNA_start.gro",
    ]

    for topology_file in topology_files:
        topology_path = data_root_dir / topology_file
        timer_start = time.perf_counter()
        count = actually_count(topology_path)
        elapsed = time.perf_counter() - timer_start

        n_atoms = len(gd.read_topology(topology_path).atoms)
        print(f"{topology_path.stem}, {n_atoms:,d} atoms, {elapsed:.2f} seconds:")
        for item in count["inventory"]:
            if item["molecular_type"] == "protein":
                seq = item["sequence"][:10]
                segments = item["number_of_segments"]
                atoms = sum(item["number_of_atoms"])
                residues = sum(item["number_of_residues"])
                print(f"  {seq:10s}     {atoms:8,d} atoms  {residues:6,d} residues  {segments:3d} segments")
            else:
                name = item["name"]
                atoms = item["number_of_atoms"]
                residues = item["number_of_residues"]
                print(f"  {name:10s}     {atoms:8,d} atoms  {residues:6,d} residues")
        print()

        
        import json
        reference_path = Path("tests") / "data" / f"{topology_path.stem}.json"
        if reference_path.exists():
            with open(reference_path, "r") as f:
                expected = json.load(f)["inventory"]

            assert len(count["inventory"]) == len(expected), (
                f"{topology_path}: Counted {len(count['inventory'])} items, expected {len(expected)}"
            )





if __name__ == "__main__":
    main()
