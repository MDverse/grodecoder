import itertools
from dataclasses import dataclass, field
from MDAnalysis import AtomGroup


@dataclass(frozen=True)
class SmallMolecule:
    """Small molecules are defined as residues with a single residue name."""

    atoms: AtomGroup
    description: str = field(default="Unknown small molecule")
    molecular_type: str = field(default="small_molecule")

    def __repr__(self) -> str:
        """Returns a string representation of the SmallMolecule."""
        return (
            f"<{self.__class__.__name__}("
            f"name={self.name!r}, "
            f"description={self.description!r}, "
            f"atoms=<AtomGroup with {self.number_of_atoms()} atoms, {self.number_of_residues()} residues>"
            ")>"
        )

    @property
    def name(self) -> str:
        """Returns the name of the residue."""
        return self.atoms.residues[0].resname

    def number_of_atoms(self) -> int:
        """Returns the number of atoms in the small molecule."""
        return len(self.atoms)

    def number_of_residues(self) -> int:
        """Returns the number of residues in the small molecule."""
        return len(self.atoms.residues)

    def dump_json(self) -> dict:
        """Returns a JSON representation of the residue."""
        return {
            "name": self.name,
            "description": self.description,
            "number_of_atoms": self.number_of_atoms(),
            "number_of_residues": self.number_of_residues(),
            "molecular_type": self.molecular_type,
            # "atom_indices": self.atom_indices.tolist(),
        }


@dataclass(frozen=True)
class Segment:
    """A segment is a group of atoms that are connected."""

    atoms: AtomGroup
    sequence: str
    molecular_type: str

    def __repr__(self) -> str:
        """Returns a string representation of the Segment."""
        seq = self.sequence
        pretty_seq = seq[:15] + "..." + seq[-5:] if len(seq) > 20 else seq
        return (
            f"<{self.__class__.__name__}("
            f"molecular_type={self.molecular_type!r}, "
            f"sequence={pretty_seq!r}, "
            f"atoms=<AtomGroup with {self.number_of_atoms()} atoms, {self.number_of_residues()} residues>"
            ")>"
        )

    @property
    def name(self) -> str:
        """Returns the name of the segment."""
        return self.atoms.residues[0].resname

    def number_of_atoms(self) -> int:
        """Returns the number of atoms in the segment."""
        return len(self.atoms)

    def number_of_residues(self) -> int:
        """Returns the number of residues in the segment."""
        return len(self.atoms.residues)


@dataclass(frozen=True)
class Inventory:
    """An inventory of molecules in a universe."""

    small_molecules: list[SmallMolecule] = field(default_factory=list)
    segments: list[Segment] = field(default_factory=list)

    def dump_json(self) -> list[dict]:
        inventory_json = []

        # Group segments by molecular type and sequence.
        by_molecular_type = itertools.groupby(self.segments, key=lambda x: x.molecular_type)
        for molecular_type, segments in by_molecular_type:
            by_sequence = itertools.groupby(segments, key=lambda x: x.sequence)
            for sequence, segments in by_sequence:
                segments_list = list(segments)
                inventory_json.append(
                    {
                        "sequence": sequence,
                        "number_of_segments": len(segments_list),
                        "number_of_atoms": [segment.number_of_atoms() for segment in segments_list],
                        "number_of_residues": [segment.number_of_residues() for segment in segments_list],
                        "molecular_type": molecular_type,
                    }
                )

        # Add small molecules.
        for small_molecule in self.small_molecules:
            inventory_json.append(
                {
                    "name": small_molecule.name,
                    "number_of_atoms": small_molecule.number_of_atoms(),
                    "number_of_residues": small_molecule.number_of_residues(),
                    "description": small_molecule.description,
                    "molecular_type": small_molecule.molecular_type,
                }
            )
        return inventory_json
