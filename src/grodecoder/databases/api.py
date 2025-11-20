import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import TypeVar

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


# ION_DB: list[Ion] = read_ion_database()
# SOLVENT_DB: list[Solvent] = read_solvent_database()
# AMINO_ACIDS_DB: list[AminoAcid] = read_amino_acid_database()
# NUCLEOTIDES_DB: list[Nucleotide] = read_nucleotide_database()
# MAD_DB: list[mad.Residue] = read_mad_database()
# CSML_DB: list[csml.Residue] = read_csml_database()
# LIPID_DB: list[Lipid] = _build_lipid_db()
# OTHER_DB: list[Residue] = _build_other_db()


@lru_cache(maxsize=1)
def _get_csml_databse():
    return read_csml_database()


@lru_cache(maxsize=1)
def _get_mad_databse():
    return read_mad_database()


def _build_lipid_db() -> list[Lipid]:
    mad_db = _get_mad_databse()
    csml_db = _get_csml_databse()

    _mad_lipid_resnames = {item.alias for item in mad_db if item.family == mad.ResidueFamily.LIPID}
    _csml_lipid_resnames = {residue.name for residue in csml_db if residue.family == csml.ResidueFamily.LIPID}

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
        for item in mad_db
        if item.family == mad.ResidueFamily.LIPID
    }
    db.update(
        {
            residue.name: Lipid(description=residue.description, residue_name=residue.name)
            for residue in csml_db
            if residue.family == csml.ResidueFamily.LIPID
        }
    )
    return list(db.values())


def _build_other_db() -> list[Residue]:
    """Builds a database of other residues that are not ions, solvents, amino acids, or nucleotides."""
    mad_db = _get_mad_databse()
    csml_db = _get_csml_databse()

    csml_other = {
        residue.name: Residue(residue_name=residue.name, description=residue.description)
        for residue in csml_db
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
        for residue in mad_db
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


@lru_cache(maxsize=1)
def get_ion_definitions() -> list[Ion]:
    """Returns the definitions of the ions in the database."""
    return read_ion_database()


@lru_cache(maxsize=1)
def get_solvent_definitions() -> list[Solvent]:
    """Returns the definitions of the solvents in the database."""
    return read_solvent_database()


@lru_cache(maxsize=1)
def get_amino_acid_definitions() -> list[AminoAcid]:
    """Returns the definitions of the amino acids in the database."""
    return read_amino_acid_database()


@lru_cache(maxsize=1)
def get_nucleotide_definitions() -> list[Nucleotide]:
    """Returns the definitions of the nucleotides in the database."""
    return read_nucleotide_database()


@lru_cache(maxsize=1)
def get_lipid_definitions() -> list[Lipid]:
    """Returns the definitions of the lipids in the database."""
    return _build_lipid_db()


@lru_cache(maxsize=1)
def get_other_definitions() -> list[Residue]:
    """Returns the definitions of other residues in the database."""
    return _build_other_db()


def get_amino_acid_name_map() -> dict[str, str]:
    """Returns a mapping of amino acid 3-letter names to 1-letter names."""
    return {aa.long_name: aa.short_name for aa in get_amino_acid_definitions()}


def get_nucleotide_name_map() -> dict[str, str]:
    """Returns a mapping of nucleotide 3-letter names to 1-letter names."""
    return {nucleotide.residue_name: nucleotide.short_name for nucleotide in get_nucleotide_definitions()}


def get_ion_names() -> set[str]:
    """Returns the names of the ions in the database."""
    return set(ion.residue_name for ion in get_ion_definitions())


def get_solvent_names() -> set[str]:
    """Returns the names of the solvents in the database."""
    return set(solvent.residue_name for solvent in get_solvent_definitions())


def get_amino_acid_names() -> set[str]:
    """Returns the names of the amino acids in the database."""
    return set(aa.long_name for aa in get_amino_acid_definitions())


def get_nucleotide_names() -> set[str]:
    """Returns the names of the nucleotides in the database."""
    return set(nucleotide.residue_name for nucleotide in get_nucleotide_definitions())


def get_lipid_names() -> set[str]:
    """Returns the names of the lipids in the database."""
    return set(lipid.residue_name for lipid in get_lipid_definitions())
