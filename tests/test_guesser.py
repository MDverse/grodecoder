import numpy as np
from MDAnalysis import Universe

from grodecoder.guesser import (
    guess_resolution,
    ResolutionValue,
    AllAtomResolutionReason,
    CoarseGrainResolutionReason,
)


def create_mock_universe(n_residues: int, n_atoms_per_residue: int) -> Universe:
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


def create_mock_all_atom(n_residues: int = 10, n_atoms_per_residue: int = 10) -> Universe:
    return create_mock_universe(n_residues, n_atoms_per_residue)


def create_mock_coarse_grain(n_residues: int = 10, n_atoms_per_residue: int = 3) -> Universe:
    return create_mock_universe(n_residues, n_atoms_per_residue)


class TestGuessResolutionAllAtom:
    """Test cases for the guess_resolution function.

    All test cases in this ensures that all-atom resolution is detected.
    """

    def test_defaulting_to_aa_when_only_one_particle_in_the_system(self):
        universe = create_mock_all_atom(n_residues=1, n_atoms_per_residue=1)
        result = guess_resolution(universe)

        assert result.value == ResolutionValue.ALL_ATOM
        assert result.reason == AllAtomResolutionReason.DEFAULT

    def test_system_has_hydrogen(self):
        """Ensures a system is detected all-atom if it contains hydrogen atoms."""
        universe = create_mock_all_atom()
        universe.atoms.types = ["H"] * len(universe.atoms)
        result = guess_resolution(universe)

        assert result.value == ResolutionValue.ALL_ATOM
        assert result.reason == AllAtomResolutionReason.HAS_HYDROGENS

    def test_has_bonds_within_all_atom_cutoff(self):
        """Ensures a system is detected all atoms when atoms are close enough.

        Even if no hydrogen is present in the system, if atoms are within a typical all-atom distance
        to each other, the system should be detected all-atom.
        """
        universe = create_mock_all_atom()
        n_atoms = len(universe.atoms)

        # Sets distances between atoms to 0.0
        universe.atoms.positions = np.zeros((n_atoms, 3))

        result = guess_resolution(universe)
        assert result.value == ResolutionValue.ALL_ATOM
        assert result.reason == AllAtomResolutionReason.RESIDUES_HAVE_BONDS_WITHIN_CUTOFF

        # Sets distances between atoms to 5.0 (distance > all-atom cutoff)
        coordinates = np.zeros((n_atoms, 3))
        for i in range(n_atoms):
            coordinates[i] = np.array((5 * i, 0.0, 0.0))
        universe.atoms.positions = coordinates

        result = guess_resolution(universe)
        assert result.value == ResolutionValue.COARSE_GRAIN
        assert result.reason == CoarseGrainResolutionReason.HAS_NO_BOND_WITHIN_ALL_ATOM_CUTOFF


class TestGuessResolutionCoarseGrain:
    """Test cases for the guess_resolution function.

    All test cases in this ensures that coarse-grain resolution is detected.
    """

    def test_model_with_one_grain_per_residue(self):
        """Ensures models with a single grain per residue are detected as coarse-grain."""
        universe = create_mock_coarse_grain(n_atoms_per_residue=1)
        result = guess_resolution(universe)

        assert result.value == ResolutionValue.COARSE_GRAIN
        assert result.reason == CoarseGrainResolutionReason.HAS_ONE_GRAIN_PER_RESIDUE

    def test_model_is_martini(self):
        """Ensures Martini models are detected as coarse-grain.

        Martini systems are expected to contain grains named BB*.
        """
        universe = create_mock_coarse_grain()
        universe.atoms[0].name = "BB1"  # a single grain named BB* is supposed to be enough
        result = guess_resolution(universe)

        assert result.value == ResolutionValue.COARSE_GRAIN
        assert result.reason == CoarseGrainResolutionReason.HAS_BB_ATOMS

    def test_protein_has_no_hydrogen(self):
        """Ensures models where protein is detected but no hydrogen atoms are found are detected as coarse-grain."""
        universe = create_mock_coarse_grain()
        universe.residues.resnames = np.repeat("ALA", len(universe.residues))
        result = guess_resolution(universe)

        assert result.value == ResolutionValue.COARSE_GRAIN
        assert result.reason == CoarseGrainResolutionReason.PROTEIN_HAS_NO_HYDROGEN

    def test_no_bond_within_all_atom_distance(self):
        """Ensures models where no bond are found within typical all-atom distance are detected as coarse-grain."""
        universe = create_mock_coarse_grain()

        # Sets distances between atoms to 5.0 (distance > all-atom cutoff)
        n_atoms = len(universe.atoms)
        coordinates = np.zeros((len(universe.atoms), 3))
        for i in range(n_atoms):
            coordinates[i] = np.array((5 * i, 0.0, 0.0))
        universe.atoms.positions = coordinates

        result = guess_resolution(universe)

        assert result.value == ResolutionValue.COARSE_GRAIN
        assert result.reason == CoarseGrainResolutionReason.HAS_NO_BOND_WITHIN_ALL_ATOM_CUTOFF
