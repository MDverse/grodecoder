from enum import StrEnum
from pydantic import BaseModel


class ResidueFamily(StrEnum):
    """Types of residues expected in CSML."""

    DETERGENT = "DETERGENT"
    ION = "ION"
    CARBOHYDRATE = "CARBOHYDRATE"
    LIPID = "LIPID"
    POLYMER = "POLYMER"
    PROTEIN = "PROTEIN"
    SMALL_MOLECULE = "SMALL MOLECULE"
    SOLVENT = "SOLVENT"


class Residue(BaseModel):
    """Model for a MAD residue."""

    name: str
    alias: str
    family: ResidueFamily
    link: str

    def __hash__(self):
        return hash((self.name, self.alias, self.family, self.link))
