from enum import StrEnum
from pydantic import BaseModel


class ResidueFamily(StrEnum):
    """Types of residues expected in CSML."""

    PROTEIN = "PROTEIN"
    NUCLEIC_ACID = "NUCLEIC_ACID"
    CARBOHYDRATE = "CARBOHYDRATE"
    LIPID = "MEMBRANE"
    SMALL_MOLECULE = "SMALL MOLECULE"
    ION = "ION"
    SOLVENT = "SOLVENT"
    DETERGENT = "DETERGENT"


class Links(BaseModel):
    """Model for CSML links."""

    view: str
    download: str | None = None


class Residue(BaseModel):
    """Model for a CSML residue.

    Attributes:
        name (str): The residue name expected in topology files.
        description (str): The full name of the residue.
        charge (int): The charge of the residue.
        links (CSMLLinks): Links to the CSML residue visualization and download.
        source (str): The name of the topology file from which the residue is taken.
        family (str): The family of the residue (e.g., protein, nucleic acid, etc.).
    """

    name: str
    description: str
    charge: int
    links: Links
    source: str
    family: ResidueFamily

    def __hash__(self):
        return hash((self.name, self.description, self.charge, self.family))
