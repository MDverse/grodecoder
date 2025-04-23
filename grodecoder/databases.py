from __future__ import annotations
import json
from pathlib import Path
from typing import Iterable

from ._typing import PathLike, Json


RESIDUE_DATABASE_PATH = Path(__file__).parent / "data" / "databases" / "residues.json"
assert RESIDUE_DATABASE_PATH.is_file()


class ResidueDatabase:
    """Stores the data from the molecule database."""

    _by_residue_name: dict[str, Json]

    def __init__(self, data: dict[str, Json]):
        self._by_residue_name = {}
        for residue in data:
            assert "residue_name" in residue, f"Missing residue_name in {residue}"
            if "residue_name" in residue:
                self._by_residue_name.setdefault(residue["residue_name"], [])
                if self._is_duplicate(residue):
                    raise DuplicateResidue(f"Duplicate residue found in the dataabse: {residue!r}")
                self._by_residue_name[residue["residue_name"]].append(residue)

    def _is_duplicate(self, residue: Json) -> bool:
        """Returns True if the residue is already in the database.

        Two residues are considered duplicates if all attributes except "description"
        are the same.
        """
        existing_residues = self._by_residue_name.get(residue["residue_name"], [])
        for existing_residue in existing_residues:
            # Check if all attributes except "description" are the same
            if all(
                existing_residue.get(attr) == residue.get(attr)
                for attr in existing_residue
                if attr != "description"
            ):
                return True
        return False

    @property
    def by_residue_name(self) -> dict[str, list[Json]]:
        """Returns the molecule database data indexed by residue name."""
        return self._by_residue_name

    @classmethod
    def from_json(cls, path: PathLike) -> ResidueDatabase:
        """Reads the grodecoder residue database.

        The database is a JSON file that contains information about
        different molecules, including their residue names.
        """
        with open(path, "r") as f:
            data = json.load(f)
        return cls(data)

    @property
    def ions(self) -> list[Json]:
        """Returns the ions in the database."""
        return [
            residue
            for residue_list in self._by_residue_name.values()
            for residue in residue_list
            if residue["family"] == "ion"
        ]

    @property
    def solvents(self) -> list[Json]:
        """Returns the solvents in the database."""
        return [
            residue
            for residue_list in self._by_residue_name.values()
            for residue in residue_list
            if residue["family"] == "solvent"
        ]

    @property
    def ion_names(self) -> list[str]:
        """Returns the names of the ions in the database."""
        return [residue["residue_name"] for residue in self.ions]

    @property
    def solvent_names(self) -> list[str]:
        """Returns the names of the solvents in the database."""
        return [residue["residue_name"] for residue in self.solvents]


RESIDUE_DATABASE = ResidueDatabase.from_json(RESIDUE_DATABASE_PATH)


class ResidueNotFound(Exception):
    """Raised when a residue is not found in the database."""


class DuplicateResidue(Exception):
    """Raised when a residue with a given name and atom names is defined multiple times in the database."""


def find_residue(residue_name: str, atom_names: Iterable[str]) -> Json:
    """Find the definitions of a given residue in the database.

    Parameters
    ----------
    residue_name : str
        The residue name to search for.

    atom_names : list[str]
        The atom names to search for.
    """
    candidate_residues = RESIDUE_DATABASE._by_residue_name.get(residue_name, [])
    if not candidate_residues:
        raise ResidueNotFound(f"No residue '{residue_name}' found in the database.")

    actual_atom_names = set(atom_names)

    # Filter the candidate residues to find those with matching atom names if provided.
    matches = []
    for residue in candidate_residues:
        expected_atom_names = set(residue["atom_names"])
        if expected_atom_names == actual_atom_names:
            matches.append(residue)

    if not matches:
        raise ResidueNotFound(f"No definition found for residue '{residue_name}' with atoms {atom_names}")

    if len(matches) > 1:
        raise DuplicateResidue(
            f"Multiple definitions found for residue '{residue_name}' with atoms {atom_names}"
        )
    return matches[0]
