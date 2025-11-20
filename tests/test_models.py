"""Tests for grodecoder.models."""

import json
from pathlib import Path

import pytest

from grodecoder.models import BaseModelWithAtomsRead, InventoryRead

TEST_DATA_ROOT_DIR = Path(__file__).parent / "data" / "regression_data"
TEST_DATA_EXPECTED_RESULTS_DIR = TEST_DATA_ROOT_DIR / "expected_results"
EXPECTED_INVENTORY_FILES = list(TEST_DATA_EXPECTED_RESULTS_DIR.glob("*.json"))


class TestBaseModelWithAtomsRead:
    def test_number_of_atoms_validation_success(self):
        try:
            BaseModelWithAtomsRead(
                atoms=[1, 2, 3],
                number_of_atoms=3,
                number_of_residues=1,
            )
        except Exception as exc:
            assert False, f"exception was raised: {exc}"

    def test_number_of_atoms_validation_fails(self):
        with pytest.raises(ValueError):
            BaseModelWithAtomsRead(
                atoms=[1, 2, 3],
                number_of_atoms=2,  # wrong number of atoms
                number_of_residues=1,
            )


def _read_json(path: str) -> dict:
    with open(path, "rb") as f:
        return json.load(f)


class TestInventoryRead:
    def test_total_number_of_atoms_success(self):
        raw_json = _read_json(EXPECTED_INVENTORY_FILES[0])
        try:
            InventoryRead.model_validate(raw_json["inventory"])
        except Exception as exc:
            assert False, f"exception was raised: {exc}"

    def test_total_number_of_atoms_fail(self):
        raw_json = _read_json(EXPECTED_INVENTORY_FILES[0])
        assert "total_number_of_atoms" in raw_json.get("inventory", {})
        raw_json["inventory"]["total_number_of_atoms"] += 1000  # sets this field to wrong value
        with pytest.raises(ValueError, match=r"field `total_number_of_atoms` \([0-9]+\) does not add up"):
            InventoryRead.model_validate(raw_json["inventory"])
