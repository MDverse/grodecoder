import json
from pathlib import Path

from .models import Residue, ResidueFamily


class DuplicateResidueError(Exception):
    """Raised when duplicate residues are found in the database."""


def read_database(path: Path | str) -> list[Residue]:
    """Reads the CSML database from a JSON file."""
    with open(path, "r") as f:
        data = json.load(f)
        db = [Residue.model_validate(entry) for entry in data]
    if len(db) != len(set(residue.name for residue in db)):
        raise DuplicateResidueError("Duplicate residues found in the database")
    return db


__all__ = ["Residue", "ResidueFamily", "read_database"]
