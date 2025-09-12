"""Unit tests for grodecoder.toputils module."""

import numpy as np
from unittest.mock import Mock

from grodecoder.toputils import (
    _first_alpha,
    formula,
    sequence,
    has_bonds,
    has_bonds_between,
    detect_chains,
    guess_resolution,
    MolecularResolution,
)


def dummy_atom(name: str):
    """Create a mock atom with a given name."""
    atom = Mock()
    atom.name = name
    return atom


def dummy_residue(atom_names: list[str] = None, resname: str = "DUM"):
    """Create a mock residue with a list of atom names."""
    residue = Mock()
    if atom_names is not None:
        residue.atoms = [dummy_atom(name) for name in atom_names]
    residue.resname = resname
    return residue


class TestFirstAlpha:
    """Test cases for the _first_alpha function."""

    def test_first_alpha_simple(self):
        assert _first_alpha("abc123") == "a"

    def test_first_alpha_numeric_start(self):
        assert _first_alpha("123abc") == "a"

    def test_first_alpha_mixed(self):
        assert _first_alpha("!@#def456") == "d"

    def test_first_alpha_uppercase(self):
        assert _first_alpha("123XYZ") == "X"

    def test_first_alpha_empty_string(self):
        assert _first_alpha("") == ""

    def test_first_alpha_no_alpha(self):
        assert _first_alpha("12345!@#") == ""


class TestFormula:
    """Test cases for the formula function."""

    def test_formula_basic(self):
        """Test formula calculation with basic elements."""
        mock_molecule = Mock()
        mock_molecule.atoms.residues = [dummy_residue(["C1", "C2", "H1", "O1"])]

        result = formula(mock_molecule)
        assert result == "C2HO"

    def test_formula_custom_elements(self):
        """Test formula with custom element list."""
        mock_molecule = Mock()
        mock_molecule.atoms.residues = [dummy_residue(["C1", "S1", "H1"])]

        result = formula(mock_molecule, elements=("C", "S"))
        assert result == "CS"

    def test_formula_no_matching_elements(self):
        """Test formula when no atoms match requested elements."""
        mock_molecule = Mock()
        mock_molecule.atoms.residues = [dummy_residue(["X1", "Y1"])]

        result = formula(mock_molecule, elements=("C", "H", "O"))
        assert result == ""


def dummy_universe(resnames: list[str]):
    """Create a mock universe with a list of residue names."""
    universe = Mock()
    universe.residues = [dummy_residue(resname=resname) for resname in resnames]
    return universe


class TestSequence:
    """Test cases for the sequence function."""

    def test_sequence_amino_acids(self):
        atoms = dummy_universe(resnames=["ALA", "GLY", "VAL"])
        result = sequence(atoms)
        assert result == "AGV"

    def test_sequence_nucleotides(self):
        atoms = dummy_universe(["DA", "DT", "DG"])
        result = sequence(atoms)
        assert result == "ATG"

    def test_sequence_unknown_residues(self):
        atoms = dummy_universe(resnames=["ALA", "UNK", "GLY"])
        result = sequence(atoms)
        assert result == "AXG"

    def test_sequence_empty(self):
        """Test sequence generation with no residues."""
        atoms = dummy_universe(resnames=[])
        result = sequence(atoms)
        assert result == ""


class TestHasBonds:
    """Test cases for the has_bonds function, which checks if a residue has any bonds."""

    def test_has_bonds_true(self):
        mock_residue = Mock()
        mock_residue.atoms.positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])

        result = has_bonds(mock_residue, cutoff=2.0)
        assert result is True

    def test_has_bonds_false(self):
        mock_residue = Mock()
        mock_residue.atoms.positions = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])

        result = has_bonds(mock_residue, cutoff=2.0)
        assert result is False

    def test_has_bonds_single_atom(self):
        """Test has_bonds with single atom."""
        mock_residue = Mock()
        mock_residue.atoms.positions = np.array([[0.0, 0.0, 0.0]])

        result = has_bonds(mock_residue, cutoff=2.0)
        assert result is False

    def test_has_bonds_custom_cutoff(self):
        """Test has_bonds with custom cutoff."""
        mock_residue = Mock()
        mock_residue.atoms.positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])

        # Should be False with default cutoff
        result1 = has_bonds(mock_residue, cutoff=2.0)
        assert result1 is False

        # Should be True with larger cutoff
        result2 = has_bonds(mock_residue, cutoff=4.0)
        assert result2 is True


class TestHasBondsBetween:
    """Test cases for the has_bonds_between function, which checks if two residues are bonded."""

    def test_has_bonds_between_true(self):
        mock_residue1 = Mock()
        mock_residue1.atoms.positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

        mock_residue2 = Mock()
        mock_residue2.atoms.positions = np.array([[2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])

        result = has_bonds_between(mock_residue1, mock_residue2, cutoff=5.0)
        assert result is True

    def test_has_bonds_between_false(self):
        mock_residue1 = Mock()
        mock_residue1.atoms.positions = np.array([[0.0, 0.0, 0.0]])

        mock_residue2 = Mock()
        mock_residue2.atoms.positions = np.array([[10.0, 0.0, 0.0]])

        result = has_bonds_between(mock_residue1, mock_residue2, cutoff=5.0)
        assert result is False

    def test_has_bonds_between_custom_cutoff(self):
        """Test has_bonds_between with custom cutoff."""
        mock_residue1 = Mock()
        mock_residue1.atoms.positions = np.array([[0.0, 0.0, 0.0]])

        mock_residue2 = Mock()
        mock_residue2.atoms.positions = np.array([[4.0, 0.0, 0.0]])

        # Should be False with small cutoff
        result1 = has_bonds_between(mock_residue1, mock_residue2, cutoff=3.0)
        assert result1 is False

        # Should be True with larger cutoff
        result2 = has_bonds_between(mock_residue1, mock_residue2, cutoff=6.0)
        assert result2 is True


class TestDetectChains:
    """Test cases for the detect_chains function."""

    def test_detect_chains_single_chain(self, monkeypatch):
        monkeypatch.setattr("grodecoder.toputils.has_bonds_between", lambda res1, res2, cutoff: True)

        mock_residue1 = Mock()
        mock_residue2 = Mock()
        mock_residue3 = Mock()

        mock_universe = Mock()
        mock_universe.residues = [mock_residue1, mock_residue2, mock_residue3]

        result = detect_chains(mock_universe, cutoff=5.0)
        assert result == [(0, 2)]

    def test_detect_chains_multiple_chains(self, monkeypatch):
        # First and second residues are bonded, but second and third are not
        def mock_bonds(res1, res2, cutoff):
            if (res1, res2) == (mock_residue1, mock_residue2):
                return True
            return False

        monkeypatch.setattr("grodecoder.toputils.has_bonds_between", mock_bonds)

        mock_residue1 = Mock()
        mock_residue2 = Mock()
        mock_residue3 = Mock()

        mock_universe = Mock()
        mock_universe.residues = [mock_residue1, mock_residue2, mock_residue3]

        result = detect_chains(mock_universe, cutoff=5.0)
        # result: first chain is residue 0 and 1, second chain is residue 2 (starts at 2, ends at 2)
        assert result == [(0, 1), (2, 2)]


class TestGuessResolution:
    """Test cases for the guess_resolution function."""

    def test_guess_resolution_all_atom(self, monkeypatch):
        """Test guess_resolution returns ALL_ATOM when bonds are found."""

        def mock_has_bonds(residue, cutoff):
            return True

        monkeypatch.setattr("grodecoder.toputils.has_bonds", mock_has_bonds)

        mock_residue = Mock()
        mock_residue.atoms = [Mock(), Mock()]

        mock_universe = Mock()
        mock_universe.residues = [mock_residue]

        result = guess_resolution(mock_universe)
        assert result == MolecularResolution.ALL_ATOM

    def test_guess_resolution_coarse_grained(self, monkeypatch):
        def mock_has_bonds(residue, cutoff):
            return False  # simulates no bonds found in any residue

        monkeypatch.setattr("grodecoder.toputils.has_bonds", mock_has_bonds)

        mock_residue = Mock()
        mock_residue.atoms = [Mock(), Mock()]

        mock_universe = Mock()
        mock_universe.residues = [mock_residue] * 5  # Ensure we have enough residues

        result = guess_resolution(mock_universe)
        assert result == MolecularResolution.COARSE_GRAINED

    def test_guess_resolution_mixed_first_has_bonds(self, monkeypatch):
        """Test guess_resolution when at least one residue has bonds."""
        mock_has_bonds = Mock(side_effect=[False, False, True, False, False])
        monkeypatch.setattr("grodecoder.toputils.has_bonds", mock_has_bonds)

        mock_residue = Mock()
        mock_residue.atoms = [Mock(), Mock()]

        mock_universe = Mock()
        mock_universe.residues = [mock_residue] * 5

        result = guess_resolution(mock_universe)
        assert result == MolecularResolution.ALL_ATOM
