import itertools
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TypeVar

from loguru import logger

from .models import Ion, Solvent, Nucleotide, AminoAcid, Lipid, Residue

from . import mad
from . import csml


DATABASES_DATA_PATH = Path(__file__).parent.parent / "data" / "databases"
assert DATABASES_DATA_PATH.is_dir(), (
    f"Database path {DATABASES_DATA_PATH} does not exist or is not a directory"
)


class NotAFileError(IOError):
    """Raised when a path is not a file."""


def assert_database_exists(path: Path):
    """Asserts that the database exists and is a file."""
    if not path.exists():
        raise FileNotFoundError(f"Database {path!r} does not exist")
    if not path.is_file():
        raise NotAFileError(f"Database {path!r} is not a file")


ION_DB_PATH = DATABASES_DATA_PATH / "ions.json"
assert_database_exists(ION_DB_PATH)

SOLVENT_DB_PATH = DATABASES_DATA_PATH / "solvents.json"
assert_database_exists(SOLVENT_DB_PATH)

AMINO_ACIDS_DB_PATH = DATABASES_DATA_PATH / "amino_acids.json"
assert_database_exists(AMINO_ACIDS_DB_PATH)

NUCLEOTIDES_DB_PATH = DATABASES_DATA_PATH / "nucleotides.json"
assert_database_exists(NUCLEOTIDES_DB_PATH)

MAD_DB_PATH = DATABASES_DATA_PATH / "mad_database.json"
assert_database_exists(MAD_DB_PATH)

CSML_DB_PATH = DATABASES_DATA_PATH / "charmm_csml_database.json"
assert_database_exists(CSML_DB_PATH)


ModelType = TypeVar("ModelType", Ion, Solvent, Nucleotide, AminoAcid)


def _read_database(path: Path, model: type[ModelType]) -> list[ModelType]:
    """Reads a database file and returns a list of entries."""
    with open(path, "r") as f:
        raw_data = [model.model_validate(entry) for entry in json.load(f)]

    # Check for duplicates.
    data = []
    counts = Counter(entry for entry in raw_data)
    for entry, count in counts.items():
        if count > 1:
            logger.warning(f"Duplicate residue {entry!r} found in database")
        else:
            data.append(entry)

    return data


def read_ion_database() -> list[Ion]:
    return _read_database(ION_DB_PATH, Ion)


def read_solvent_database() -> list[Solvent]:
    return _read_database(SOLVENT_DB_PATH, Solvent)


def read_amino_acid_database() -> list[AminoAcid]:
    return _read_database(AMINO_ACIDS_DB_PATH, AminoAcid)


def read_nucleotide_database() -> list[Nucleotide]:
    return _read_database(NUCLEOTIDES_DB_PATH, Nucleotide)


def read_mad_database() -> list[mad.Residue]:
    """Reads the MAD database and returns a list of residues."""
    return _read_database(MAD_DB_PATH, mad.Residue)


def read_csml_database() -> list[csml.Residue]:
    """Reads the CSML database and returns a list of residues."""
    return _read_database(CSML_DB_PATH, csml.Residue)


ION_DB: list[Ion] = read_ion_database()
SOLVENT_DB: list[Solvent] = read_solvent_database()
AMINO_ACIDS_DB: list[AminoAcid] = read_amino_acid_database()
NUCLEOTIDES_DB: list[Nucleotide] = read_nucleotide_database()

MAD_DB: list[mad.Residue] = read_mad_database()
CSML_DB: list[csml.Residue] = read_csml_database()


def _build_lipid_db() -> list[Lipid]:
    _mad_lipid_resnames = {item.alias for item in MAD_DB if item.family == mad.ResidueFamily.LIPID}
    _csml_lipid_resnames = {residue.name for residue in CSML_DB if residue.family == csml.ResidueFamily.LIPID}

    if False:
        # IMPORTANT:
        #   Right now, having duplicates is not an issue has we are only interested in the residue names
        #   and the molecular type (i.e. lipid in this case).
        #   If we ever start checking atom names or other properties, we will need to handle duplicates
        #   properly.
        _duplicates = _mad_lipid_resnames & _csml_lipid_resnames
        if _duplicates:
            logger.warning(f"Duplicate lipid residue names found in MAD and CSML databases: {_duplicates}")

    db = {
        item.alias: Lipid(description=item.name, residue_name=item.alias)
        for item in MAD_DB
        if item.family == mad.ResidueFamily.LIPID
    }
    db.update(
        {
            residue.name: Lipid(description=residue.description, residue_name=residue.name)
            for residue in CSML_DB
            if residue.family == csml.ResidueFamily.LIPID
        }
    )
    return list(db.values())


LIPID_DB: list[Lipid] = _build_lipid_db()


def _build_other_db() -> list[Residue]:
    """Builds a database of other residues that are not ions, solvents, amino acids, or nucleotides."""
    csml_other = {
        residue.name: Residue(residue_name=residue.name, description=residue.description)
        for residue in CSML_DB
        if residue.family
        not in {
            csml.ResidueFamily.PROTEIN,
            csml.ResidueFamily.NUCLEIC_ACID,
            csml.ResidueFamily.LIPID,
            csml.ResidueFamily.SOLVENT,
            csml.ResidueFamily.ION,
        }
    }

    mad_other = {
        residue.alias: Residue(residue_name=residue.alias, description=residue.name)
        for residue in MAD_DB
        if residue.family
        not in {
            mad.ResidueFamily.PROTEIN,
            mad.ResidueFamily.LIPID,
            mad.ResidueFamily.ION,
            mad.ResidueFamily.SOLVENT,
        }
    }

    by_name = {**csml_other, **mad_other}
    return list(by_name.values())


OTHER_DB: list[Residue] = _build_other_db()


def get_other_definitions() -> list[Residue]:
    """Returns the definitions of other residues in the database."""
    return OTHER_DB


def get_ion_definitions() -> list[Ion]:
    """Returns the definitions of the ions in the database."""
    return ION_DB


def get_solvent_definitions() -> list[Solvent]:
    """Returns the definitions of the solvents in the database."""
    return SOLVENT_DB


def get_amino_acid_definitions() -> list[AminoAcid]:
    """Returns the definitions of the amino acids in the database."""
    return AMINO_ACIDS_DB


def get_amino_acid_name_map() -> dict[str, str]:
    """Returns a mapping of amino acid 3-letter names to 1-letter names."""
    return {aa.long_name: aa.short_name for aa in AMINO_ACIDS_DB}


def get_nucleotide_name_map() -> dict[str, str]:
    """Returns a mapping of nucleotide 3-letter names to 1-letter names."""
    return {nucleotide.residue_name: nucleotide.short_name for nucleotide in NUCLEOTIDES_DB}


def get_nucleotide_definitions() -> list[Nucleotide]:
    """Returns the definitions of the nucleotides in the database."""
    return NUCLEOTIDES_DB


def get_lipid_definitions() -> list[Lipid]:
    """Returns the definitions of the lipids in the database."""
    return LIPID_DB


def get_ion_names() -> set[str]:
    """Returns the names of the ions in the database."""
    return set(ion.residue_name for ion in ION_DB)


def get_solvent_names() -> set[str]:
    """Returns the names of the solvents in the database."""
    return set(solvent.residue_name for solvent in SOLVENT_DB)


def get_amino_acid_names() -> set[str]:
    """Returns the names of the amino acids in the database."""
    return set(aa.long_name for aa in AMINO_ACIDS_DB)


def get_nucleotide_names() -> set[str]:
    """Returns the names of the nucleotides in the database."""
    return set(nucleotide.residue_name for nucleotide in NUCLEOTIDES_DB)


def get_lipid_names() -> set[str]:
    """Returns the names of the lipids in the database."""
    return set(lipid.residue_name for lipid in LIPID_DB)


@dataclass(frozen=True)
class ResidueDatabase:
    """Database of residues."""

    ions: list[Ion]
    solvents: list[Solvent]
    amino_acids: list[AminoAcid]
    nucleotides: list[Nucleotide]

    def __post_init__(self):
        names = {
            "ions": {ion.residue_name for ion in self.ions},
            "solvents": {solvent.residue_name for solvent in self.solvents},
            "amino_acids": {aa.long_name for aa in self.amino_acids},
            "nucleotides": {nucleotide.residue_name for nucleotide in self.nucleotides},
        }

        combinations = itertools.combinations(names.keys(), 2)
        for lhs, rhs in combinations:
            duplicates = names[lhs].intersection(names[rhs])
            if duplicates:
                logger.warning(
                    f"Residue names {duplicates} are defined in multiple families: {lhs} and {rhs}"
                )


class ResidueNotFound(Exception):
    """Raised when a residue with a given name and atom names is not found in the database."""


class DuplicateResidue(Exception):
    """Raised when a residue with a given name and atom names is defined multiple times in the database."""


def _find_using_atom_names(
    residue_name: str, atom_names: Iterable[str], database: list[Ion | Solvent]
) -> Ion | Solvent | None:
    candidate_residues = [ion for ion in database if ion.residue_name == residue_name]
    if not candidate_residues:
        return None

    actual_atom_names = set(atom_names)
    matching_residues = [ion for ion in candidate_residues if set(ion.atom_names) == actual_atom_names]

    if len(matching_residues) == 0:
        raise ResidueNotFound(f"No residue '{residue_name}' found with atom names {actual_atom_names}")

    elif len(matching_residues) > 1:
        raise DuplicateResidue(
            f"Multiple residues '{residue_name}' found with atom names {actual_atom_names}"
        )

    return matching_residues[0]
