from MDAnalysis import AtomGroup
import numpy as np

from .databases import ResidueDefinition
from .toputils import detect_chains, formula, sequence


class HasAtoms:
    "Interface for objects that have atoms."

    def __init__(self, atoms: AtomGroup):
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

    def __init__(self, definition: ResidueDefinition, atoms: AtomGroup):
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
        return formula(self.atoms)


class ProteinSegment(HasAtoms):
    """Holds a segment of a protein."""

    def __init__(self, atoms: AtomGroup):
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

    def __init__(self, atoms: AtomGroup):
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
