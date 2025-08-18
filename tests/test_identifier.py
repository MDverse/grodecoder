"""Unit tests for grodecoder.identifier module."""

import warnings
from pathlib import Path
from unittest.mock import Mock, patch


import MDAnalysis as mda
import numpy as np
import pytest

from grodecoder.databases.models import Ion, Solvent
from grodecoder.identifier import (
    _find_methanol,
    _get_nucleic_segments,
    _get_protein_segments,
    _iter_chains,
    _remove_identified_atoms,
    _select_nucleic,
    _select_protein,
    _unique_definitions,
    identify,
    identify_small_molecule,
)
from grodecoder.models import Inventory, Segment, SmallMolecule


def read_pdb_no_warning_missing_elements(path: Path) -> mda.Universe:
    """Reads a PDB file without warning about missing elements."""
    with warnings.catch_warnings() as w:
        warnings.simplefilter("ignore", category=UserWarning)
        return mda.Universe(path)


class TestFindMethanol:
    """Test cases for the _find_methanol function."""

    TEST_DATA_DIR = Path(__file__).parent / "data" / "find_methanol"

    def test_returns_empty_list_when_has_no_methanol(self):
        """Tests that an empty list is returned when system does not contain any "MET" residue."""
        atoms = read_pdb_no_warning_missing_elements(self.TEST_DATA_DIR / "water.pdb")
        result = _find_methanol(atoms)
        assert result == []

    def test_returns_methanol_indices(self):
        """Tests that methanol indices are returned correctly."""
        atoms = mda.Universe(self.TEST_DATA_DIR / "methanol.pdb")  # Contains 200 methanol molecules
        result = _find_methanol(atoms)
        assert isinstance(result, list)
        assert len(result) == 200 * 6  # 6 atoms per methanol molecule
        assert all(isinstance(idx, int) for idx in result)

    def test_does_not_return_protein(self):
        """Tests that methanol indices do not include protein atoms."""
        atoms = read_pdb_no_warning_missing_elements(self.TEST_DATA_DIR / "protein.pdb")
        result = _find_methanol(atoms)
        assert result == []

    def test_methanol_protein_mix(self):
        """Tests that methanol indices are returned correctly when mixed with protein."""
        atoms = read_pdb_no_warning_missing_elements(self.TEST_DATA_DIR / "system.pdb")
        result = _find_methanol(atoms)
        assert isinstance(result, list)
        assert len(result) == 200 * 6


class TestSelectProtein:
    """Test cases for the _select_protein function."""

    TEST_DATA_DIR = Path(__file__).parent / "data" / "find_methanol"

    def test_select_protein_no_methanol(self):
        """Test selecting protein atoms without methanol interference."""
        atoms = read_pdb_no_warning_missing_elements(self.TEST_DATA_DIR / "protein.pdb")

        result = _select_protein(atoms)

        expected_number_of_atoms = 1158  # from VMD
        expected_number_of_residues = 148  # from VMD

        assert len(result.atoms) == expected_number_of_atoms
        assert len(result.residues) == expected_number_of_residues

    def test_select_protein_with_methanol(self):
        """Test selecting protein atoms with methanol interference."""
        atoms = read_pdb_no_warning_missing_elements(self.TEST_DATA_DIR / "system.pdb")

        result = _select_protein(atoms)

        expected_number_of_atoms = 1158
        expected_number_of_residues = 148

        assert len(result.atoms) == expected_number_of_atoms
        assert len(result.residues) == expected_number_of_residues


class TestSelectNucleic:
    """Test cases for the _select_nucleic function."""

    TEST_DATA_DIR = Path(__file__).parent / "data" / "select_nucleic"

    def test_select_dna(self):
        """Test selecting dna atoms."""
        atoms = read_pdb_no_warning_missing_elements(self.TEST_DATA_DIR / "dna_water_ions.gro")

        result = _select_nucleic(atoms)

        expected_number_of_atoms = 1132
        expected_number_of_residues = 36

        assert len(result.atoms) == expected_number_of_atoms
        assert len(result.residues) == expected_number_of_residues

    def test_select_rna(self):
        """Test selecting rna atoms."""
        atoms = read_pdb_no_warning_missing_elements(self.TEST_DATA_DIR / "rna_water_ions.gro")

        result = _select_nucleic(atoms)

        expected_number_of_atoms = 1168
        expected_number_of_residues = 36

        assert len(result.atoms) == expected_number_of_atoms
        assert len(result.residues) == expected_number_of_residues


class TestIterChains:
    """Test cases for the _iter_chains function."""

    TEST_DATA_DIR = Path(__file__).parent / "data"

    def test_iter_chains_empty_atoms(self):
        """Test iterating chains with empty atom group."""
        mock_atoms = Mock()
        mock_atoms.__len__ = Mock(return_value=0)

        result = list(_iter_chains(mock_atoms))

        assert result == []

    def test_iter_chains_single_chain(self):
        """Test iterating chains with a single chain."""
        universe = mda.Universe(self.TEST_DATA_DIR / "barstar.gro")
        protein_atoms = universe.select_atoms("protein")

        result = list(_iter_chains(protein_atoms))

        assert len(result) == 1
        assert len(result[0].atoms) == 1434
        assert len(result[0].residues) == 89

    def test_iter_chains_multiple_chains(self):
        """Test iterating chains with multiple chains."""
        universe = mda.Universe(self.TEST_DATA_DIR / "1BRS.pdb")
        protein_atoms = universe.select_atoms("protein")

        result = list(_iter_chains(protein_atoms))

        # 1BRS.gro is the complex barstar-barnase, 3 times, i.e. total of 6 chains.
        assert len(result) == 6

        # Build the expected result vector.
        # Use MDAnalysis internal selection protocol to get expected number of atoms and residues
        # for each chain.
        # Should be sufficiently accurate for this test, as we made sure that chains are correctly labeled
        # in this particular test file.
        expected = []
        for chain in ("A", "B", "C", "D", "E", "F"):
            chain_atoms = protein_atoms.select_atoms(f"protein and segid {chain}")
            expected.append(
                {
                    "number_of_atoms": len(chain_atoms.atoms),
                    "number_of_residues": len(chain_atoms.residues),
                }
            )

        for i, chain in enumerate(result):
            assert len(chain.atoms) == expected[i]["number_of_atoms"]
            assert len(chain.residues) == expected[i]["number_of_residues"]


class TestGetProteinSegments:
    """Test cases for the _get_protein_segments function."""

    TEST_DATA_DIR = Path(__file__).parent / "data"

    def test_get_protein_segments_single_chain(self):
        universe = mda.Universe(self.TEST_DATA_DIR / "barstar.gro")
        segments = _get_protein_segments(universe)
        assert len(segments) == 1
        assert isinstance(segments[0], Segment)

    def test_get_protein_segments_multiple_chains(self):
        universe = mda.Universe(self.TEST_DATA_DIR / "1BRS.pdb")
        segments = _get_protein_segments(universe)

        assert len(segments) == 6

        # Check that each segment has the correct molecular type
        for segment in segments:
            assert isinstance(segment, Segment)
            assert segment.molecular_type == "protein"


class TestGetNucleicSegments:
    """Test cases for the _get_nucleic_segments function."""

    TEST_DATA_DIR = Path(__file__).parent / "data" / "select_nucleic"

    def test_get_nucleic_segments_dna(self):
        universe = mda.Universe(self.TEST_DATA_DIR / "dna_water_ions.gro")
        segments = _get_nucleic_segments(universe)

        assert len(segments) == 1
        assert isinstance(segments[0], Segment)
        assert segments[0].molecular_type == "nucleic_acid"

    def test_get_nucleic_segments_rna(self):
        universe = mda.Universe(self.TEST_DATA_DIR / "rna_water_ions.gro")
        segments = _get_nucleic_segments(universe)

        assert len(segments) == 1
        assert isinstance(segments[0], Segment)
        assert segments[0].molecular_type == "nucleic_acid"


class TestUniqueDefinitions:
    """Test cases for the _unique_definitions function."""

    def test_unique_definitions_single(self):
        """Test unique definitions with single definition."""
        def1 = Ion(residue_name="ION", description="Ion", atom_names=["I"])
        result = _unique_definitions([def1])

        assert len(result) == 1
        assert result["ION"] == def1

    def test_unique_definitions_duplicates(self):
        """Test unique definitions with duplicate residue names (first one wins)."""
        def1 = Ion(residue_name="ION", description="First Ion", atom_names=["I"])
        def2 = Ion(residue_name="ION", description="Second Ion", atom_names=["I", "I2"])  # unique is atom name insensitive
        result = _unique_definitions([def1, def2])

        assert len(result) == 1
        assert result["ION"] == def1  # First one should be kept
