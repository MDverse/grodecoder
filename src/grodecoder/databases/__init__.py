from functools import lru_cache

from .api import (
    DATABASES_DATA_PATH,
    get_amino_acid_definitions,
    get_amino_acid_name_map,
    get_amino_acid_names,
    get_ion_definitions,
    get_ion_names,
    get_lipid_definitions,
    get_lipid_names,
    get_nucleotide_definitions,
    get_nucleotide_name_map,
    get_nucleotide_names,
    get_other_definitions,
    get_solvent_definitions,
    get_solvent_names,
)
from .models import (
    AminoAcid,
    Ion,
    Lipid,
    Nucleotide,
    ResidueDefinition,
    Solvent,
)

__all__ = [
    "get_database_version",
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


@lru_cache(maxsize=1)
def get_database_version() -> str:
    with open(DATABASES_DATA_PATH / "version.txt") as f:
        return f.read().strip()
