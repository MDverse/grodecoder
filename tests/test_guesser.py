import numpy as np
from MDAnalysis import Universe

from grodecoder.guesser import (
    guess_resolution,
    ResolutionValue,
    AllAtomResolutionReason,
    CoarseGrainResolutionReason,
)


def create_mock_universe(n_residues: int = 10, n_atoms_per_residue: int = 10) -> Universe:
    n_atoms = n_residues * n_atoms_per_residue

    resindices = np.repeat(range(n_residues), n_atoms_per_residue)
    assert len(resindices) == n_atoms

    segindices = [0] * n_residues

    universe = Universe.empty(
        n_atoms,
        n_residues=n_residues,
        atom_resindex=resindices,
        residue_segindex=segindices,
        trajectory=True,  # required to add coordinates
    )

    universe.add_TopologyAttr("name", ["C"] * n_atoms)
    universe.add_TopologyAttr("type", ["C"] * n_atoms)
    universe.add_TopologyAttr("resname", ["UNK"] * n_residues)
    universe.add_TopologyAttr("resid", list(range(1, n_residues + 1)))

    return universe


class TestGuessResolutionAllAtom:
    """Test cases for the guess_resolution function."""

    def test_defaulting_to_aa_when_only_one_particle_in_the_system(self):
        universe = create_mock_universe(n_residues=1, n_atoms_per_residue=1)
        result = guess_resolution(universe)

        assert result.value == ResolutionValue.ALL_ATOM
        assert result.reason == AllAtomResolutionReason.DEFAULT

    def test_system_has_hydrogen(self):
        """Ensures a system is detected all-atom if it contains hydrogen atoms."""
        universe = create_mock_universe()
        universe.atoms.types = ["H"] * len(universe.atoms)
        result = guess_resolution(universe)

        assert result.value == ResolutionValue.ALL_ATOM
        assert result.reason == AllAtomResolutionReason.HAS_HYDROGENS

    def test_has_bonds_within_all_atom_cutoff(self):
        """Ensures a system is detected all atoms when atoms are close enough.

        Even if no hydrogen is present in the system, if atoms are within a typical all-atom distance
        to each other, the system should be detected all-atom.
        """
        universe = create_mock_universe()
        n_atoms = len(universe.atoms)

        # Sets distances between atoms to 0.0
        universe.atoms.positions = np.zeros((n_atoms, 3))

        result = guess_resolution(universe)
        assert result.value == ResolutionValue.ALL_ATOM
        assert result.reason == AllAtomResolutionReason.RESIDUES_HAVE_BONDS_WITHIN_CUTOFF

        # Sets distances between atoms to 5.0
        coordinates = np.zeros((n_atoms, 3))
        for i in range(n_atoms):
            coordinates[i] = np.array((5 * i, 0.0, 0.0))
        universe.atoms.positions = coordinates

        result = guess_resolution(universe)
        assert result.value == ResolutionValue.COARSE_GRAIN
        assert result.reason == CoarseGrainResolutionReason.HAS_NO_BOND_WITHIN_ALL_ATOM_CUTOFF

    # def test_guess_resolution_all_atom(self, monkeypatch):
    #     """Test guess_resolution returns ALL_ATOM when bonds are found."""
    #
    #     def mock_has_bonds(residue, cutoff):
    #         return True
    #
    #     monkeypatch.setattr("grodecoder.toputils.has_bonds", mock_has_bonds)
    #
    #     universe = create_mock_universe()
    #     universe.atoms.types = ["C"]
    #
    #
    #     result = guess_resolution(universe)
    #     # ic(result)
    #     assert False

    # def test_guess_resolution_coarse_grained(self, monkeypatch):
    #     def mock_has_bonds(residue, cutoff):
    #         return False  # simulates no bonds found in any residue
    #
    #     monkeypatch.setattr("grodecoder.toputils.has_bonds", mock_has_bonds)
    #
    #     mock_residue = Mock()
    #     mock_residue.atoms = [Mock(), Mock()]
    #
    #     mock_universe = Mock()
    #     mock_universe.residues = [mock_residue] * 5  # Ensure we have enough residues
    #
    #     result = guess_resolution(mock_universe, cutoff_distance=12)
    #     assert result == ResolutionValue.COARSE_GRAINED
    #
    # def test_guess_resolution_mixed_first_has_bonds(self, monkeypatch):
    #     """Test guess_resolution when at least one residue has bonds."""
    #     mock_has_bonds = Mock(side_effect=[False, False, True, False, False])
    #     monkeypatch.setattr("grodecoder.toputils.has_bonds", mock_has_bonds)
    #
    #     mock_residue = Mock()
    #     mock_residue.atoms = [Mock(), Mock()]
    #
    #     mock_universe = Mock()
    #     mock_universe.residues = [mock_residue] * 5
    #
    #     result = guess_resolution(mock_universe, cutoff_distance=12)
    #     assert result == ResolutionValue.ALL_ATOM
