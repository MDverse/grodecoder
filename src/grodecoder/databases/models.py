from typing import Protocol
from pydantic import BaseModel


class ResidueDefinition(Protocol):
    """Protocol for a residue definition."""

    residue_name: str
    description: str

    def __hash__(self):
        pass


class Residue(BaseModel):
    """Model for a residue."""

    description: str
    residue_name: str

    def __hash__(self):
        return hash((self.residue_name, self.description))


class Ion(BaseModel):
    """Model for an ion."""

    description: str
    residue_name: str
    atom_names: list[str]

    def __hash__(self):
        return hash((self.residue_name, frozenset(self.atom_names)))


class Solvent(BaseModel):
    """Model for a solvent."""

    description: str
    residue_name: str
    atom_names: list[str]

    def __hash__(self):
        return hash((self.residue_name, frozenset(self.atom_names)))


class Nucleotide(BaseModel):
    """Model for a nucleotide."""

    description: str
    residue_name: str
    short_name: str

    def __hash__(self):
        return hash(self.residue_name)


class AminoAcid(BaseModel):
    """Model for an amino acid."""

    description: str
    long_name: str
    short_name: str

    def __hash__(self):
        return hash(self.long_name)


class Lipid(BaseModel):
    """Model for a lipid."""

    description: str
    residue_name: str

    def __hash__(self):
        return hash(self.residue_name)
