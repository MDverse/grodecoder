from enum import StrEnum, auto
from itertools import islice
from typing import Self

import MDAnalysis as MDA
import numpy as np
from loguru import logger
from pydantic import BaseModel, ConfigDict, model_validator

from .toputils import has_bonds


class ResolutionValue(StrEnum):
    ALL_ATOM = auto()
    COARSE_GRAIN = auto()


class AllAtomResolutionReason(StrEnum):
    DEFAULT = auto()
    HAS_HYDROGENS = auto()
    RESIDUES_HAVE_BONDS_WITHIN_CUTOFF = auto()


class CoarseGrainResolutionReason(StrEnum):
    HAS_BB_ATOMS = auto()
    HAS_ONE_GRAIN_PER_RESIDUE = auto()
    HAS_NO_BOND_WITHIN_ALL_ATOM_CUTOFF = auto()
    PROTEIN_HAS_NO_HYDROGEN = auto()


class MolecularResolution(BaseModel):
    value: ResolutionValue
    reason: AllAtomResolutionReason | CoarseGrainResolutionReason | None = None

    model_config = ConfigDict(validate_assignment=True)

    def is_all_atom(self) -> bool:
        return self.value == ResolutionValue.ALL_ATOM

    def is_coarse_grain(self) -> bool:
        return self.value == ResolutionValue.COARSE_GRAIN

    @model_validator(mode="after")
    def check_reason(self) -> Self:
        """Validate that reason is compatible with value."""
        if self.value == ResolutionValue.ALL_ATOM:
            if self.reason is not None and not isinstance(self.reason, AllAtomResolutionReason):
                raise ValueError(
                    f"reason must be AllAtomResolutionReason when value is 'all-atom', "
                    f"got {type(self.reason).__name__}"
                )
        elif self.value == ResolutionValue.COARSE_GRAIN:
            if self.reason is not None and not isinstance(self.reason, CoarseGrainResolutionReason):
                raise ValueError(
                    f"reason must be CoarseGrainResolutionReason when value is 'coarse-grain', "
                    f"got {type(self.reason).__name__}"
                )
        return self

    @classmethod
    def AllAtomWithHydrogen(cls) -> Self:
        """System that contains hydrogen atoms."""
        return cls(value=ResolutionValue.ALL_ATOM, reason=AllAtomResolutionReason.HAS_HYDROGENS)

    @classmethod
    def AllAtomWithStandardBonds(cls) -> Self:
        """System in which atoms within residues have distances that are typical of all-atom models."""
        return cls(
            value=ResolutionValue.ALL_ATOM, reason=AllAtomResolutionReason.RESIDUES_HAVE_BONDS_WITHIN_CUTOFF
        )

    @classmethod
    def CoarseGrainMartini(cls) -> Self:
        """System that contains atoms named BB*."""
        return cls(value=ResolutionValue.COARSE_GRAIN, reason=CoarseGrainResolutionReason.HAS_BB_ATOMS)

    @classmethod
    def CoarseGrainSingleParticle(cls) -> Self:
        """System with residues made of a single particle."""
        return cls(
            value=ResolutionValue.COARSE_GRAIN, reason=CoarseGrainResolutionReason.HAS_ONE_GRAIN_PER_RESIDUE
        )

    @classmethod
    def CoarseGrainOther(cls) -> Self:
        """Other coarse grain systems, typically with bond length between particle greater than standard
        all-atom models."""
        return cls(
            value=ResolutionValue.COARSE_GRAIN,
            reason=CoarseGrainResolutionReason.HAS_NO_BOND_WITHIN_ALL_ATOM_CUTOFF,
        )

    @classmethod
    def CoarseGrainProteinHasNoHydrogen(cls) -> Self:
        """Protein is detected and residues do not have hydrogen atoms."""
        return cls(
            value=ResolutionValue.COARSE_GRAIN, reason=CoarseGrainResolutionReason.PROTEIN_HAS_NO_HYDROGEN
        )


def _has_hydrogen(model: MDA.AtomGroup) -> bool:
    """Returns True if a model contains hydrogen atoms."""
    return "H" in model.atoms.types


def _is_martini(model: MDA.AtomGroup) -> bool:
    """Returns True if a model contains atoms named BB*."""
    return bool(np.any(np.char.startswith(model.atoms.names.astype("U"), "BB")))


def _has_bonds_within_all_atom_cutoff(model: MDA.AtomGroup, cutoff_distance: float) -> bool:
    for residue in model.residues:
        if has_bonds(residue, cutoff_distance):
            return True
    return False


def _has_protein(model: MDA.AtomGroup) -> bool:
    return len(model.select_atoms("protein").atoms) > 0


def _protein_has_hydrogen(model: MDA.AtomGroup) -> bool:
    """Return True if protein residues have hydrogen atoms."""
    return _has_hydrogen(model.select_atoms("protein"))


def guess_resolution(universe, all_atom_cutoff_distance: float = 1.6) -> MolecularResolution:
    """Guesses a system resolution (all-atom or coarse-grain)."""
    # Only one atom in the system: defaulting to all-atom.
    if len(universe.atoms) == 1:
        return MolecularResolution(value=ResolutionValue.ALL_ATOM, reason=AllAtomResolutionReason.DEFAULT)

    # Selects the first five residues with at least two atoms.
    resindexes = list(
        islice((residue.resindex for residue in universe.residues if len(residue.atoms) > 1), 5)
    )

    # If no residue with more than one atom, definitely coarse grain.
    no_residue_with_more_than_1_particle = len(resindexes) == 0
    if no_residue_with_more_than_1_particle:
        logger.debug("No residues with more than one atom: resolution is coarse grain")
        return MolecularResolution.CoarseGrainSingleParticle()

    small_u: MDA.AtomGroup = universe.select_atoms(f"resindex {' '.join(str(i) for i in resindexes)}")

    # If we find any hydrogen atom, it's all-atom.
    if _has_hydrogen(small_u):
        logger.debug("Found hydrogen atoms: resolution is all-atom")
        return MolecularResolution.AllAtomWithHydrogen()

    # If we find any atom named "BB*", it's Martini (coarse grain).
    if _is_martini(small_u):
        logger.debug("Found residues named BB*: resolution is coarse grain")
        return MolecularResolution.CoarseGrainMartini()

    if _has_protein(universe) and not _protein_has_hydrogen(universe):
        logger.debug("Found protein without hydrogen: resolution is coarse grain")
        return MolecularResolution.CoarseGrainProteinHasNoHydrogen()

    # Last chance: if we find any bond within a given distance, it's all-atom.
    # If we reach this point, it means that, for some reason, no hydrogen atom was detected before.
    if _has_bonds_within_all_atom_cutoff(small_u, all_atom_cutoff_distance):
        logger.debug("Found bonds within all-atom distance cutoff: resolution is all-atom")
        return MolecularResolution.AllAtomWithStandardBonds()

    # Coarse grain not detected before, no bonds within cutoff distance, it's coarse grain.
    logger.debug("No bonds found within all-atom distance cutoff: resolution is coarse grain")
    return MolecularResolution.CoarseGrainOther()
