"""Integration tests for grodecoder.toputils using real MDAnalysis data."""

import pytest
from pathlib import Path

import MDAnalysis as mda

from grodecoder.toputils import (
    formula,
    sequence,
    has_bonds,
    has_bonds_between,
    detect_chains,
    guess_resolution,
    MolecularResolution,
)

TEST_DATA_DIR = Path(__file__).parent.parent / "data"
GRO_SMALL = TEST_DATA_DIR / "barstar_water_ions.gro"
GRO_BIG = TEST_DATA_DIR / "1BRS.pdb"
GRO_CG = TEST_DATA_DIR / "noriega_CG_CRD_3CAL.gro"


class TestTopUtilsSmallUniverse:
    @pytest.fixture
    def small_universe(self):
        """Create a small MDAnalysis Universe for testing."""
        return mda.Universe(GRO_SMALL)

    @pytest.fixture
    def protein_atoms(self, small_universe):
        """Get protein atoms from the small universe."""
        return small_universe.select_atoms("protein")

    def test_formula(self, small_universe):
        # Get the first residue
        first_residue = small_universe.atoms.residues[0]

        result = formula(first_residue)

        assert first_residue.resname == "LYSH"
        assert result == "C6H15ON2"

    def test_sequence(self, protein_atoms):
        """Test sequence function with real protein data."""
        if len(protein_atoms.residues) == 0:
            pytest.skip("No protein residues in test data")

        result = sequence(protein_atoms)

        # Uniprot p11540 (first residue is missing in the structure)
        expected = "MKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGCDITIILS"[1:]

        assert isinstance(result, str)
        assert len(result) == len(protein_atoms.residues)
        assert result == expected

    def test_has_bonds(self, small_universe):
        """Test has_bonds function with real residue data."""
        first_residue = small_universe.residues[0]
        result = has_bonds(first_residue, cutoff=2.0)
        assert result is True

    def test_has_bonds_between_residues(self, small_universe):
        """Test has_bonds_between function with real residue data."""
        residue1 = small_universe.residues[0]
        residue2 = small_universe.residues[1]

        result = has_bonds_between(residue1, residue2, cutoff=2.0)
        assert result is True

    def test_detect_chains(self, protein_atoms):
        """Test detect_chains function with real protein data."""
        result = detect_chains(protein_atoms, cutoff=5.0)

        # Should return a list of tuples of (start, end) indices.
        assert isinstance(result, list)
        assert all(isinstance(segment, tuple) and len(segment) == 2 for segment in result)
        assert all(isinstance(start, int) and isinstance(end, int) for start, end in result)

        # PDB_small has 1 chain
        assert len(result) == 1
        assert result[0] == (0, len(protein_atoms.residues) - 1)


    def test_guess_resolution(self, small_universe):
        """Test guess_resolution function with real data."""
        result = guess_resolution(small_universe)
        assert result == MolecularResolution.ALL_ATOM


def test_guess_resolution_cg():
    """Test guess_resolution with coarse-grained data."""
    universe = mda.Universe(GRO_CG)
    result = guess_resolution(universe)
    assert result == MolecularResolution.COARSE_GRAINED


def test_detect_chains_big():
    """Test detect_chains with a larger universe."""
    universe = mda.Universe(GRO_BIG)
    protein_atoms = universe.select_atoms("protein")
    
    result = detect_chains(protein_atoms)

    # GRO_BIG (1BRS) is the complex barstar-barnase, 3 times, i.e. total of 6 chains.
    lengths = [108, 110, 108, 87, 86, 89]  # from VMD
    expected = []
    for i, length in enumerate(lengths):
        start = sum(lengths[:i])
        end = start + length - 1
        expected.append((start, end))
    # expected: [(0, 107), (108, 217), (218, 325), (326, 412), (413, 498), (499, 587)]

    assert result == expected
