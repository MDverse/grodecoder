from .api import (
    get_amino_acid_definitions,
    get_amino_acid_names,
    get_amino_acid_name_map,
    get_ion_definitions,
    get_ion_names,
    get_lipid_definitions,
    get_lipid_names,
    get_nucleotide_definitions,
    get_nucleotide_names,
    get_nucleotide_name_map,
    get_other_definitions,
    get_solvent_definitions,
    get_solvent_names,
    DATABASES_DATA_PATH,
)
from .models import (
    ResidueDefinition,
    AminoAcid,
    Ion,
    Lipid,
    Nucleotide,
    Solvent,
)

__all__ = [
    "get_amino_acid_definitions",
    "get_amino_acid_names",
    "get_amino_acid_name_map",
    "get_ion_definitions",
    "get_ion_names",
    "get_lipid_definitions",
    "get_lipid_names",
    "get_nucleotide_definitions",
    "get_nucleotide_names",
    "get_nucleotide_name_map",
    "get_other_definitions",
    "get_solvent_definitions",
    "get_solvent_names",
    "ResidueDefinition",
    "Ion",
    "Lipid",
    "Solvent",
    "AminoAcid",
    "Nucleotide",
]


with open(DATABASES_DATA_PATH / "version.txt") as f:
    __version__ = f.read().strip()
