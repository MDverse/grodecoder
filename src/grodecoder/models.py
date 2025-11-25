from __future__ import annotations

from enum import StrEnum
from typing import Protocol

from MDAnalysis import AtomGroup
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    SerializationInfo,
    SerializerFunctionWrapHandler,
    computed_field,
    field_serializer,
    model_serializer,
    model_validator,
)

from . import toputils


class MolecularResolution(StrEnum):
    COARSE_GRAINED = "coarse-grained"
    ALL_ATOM = "all-atom"


class MolecularType(StrEnum):
    PROTEIN = "protein"
    NUCLEIC = "nucleic_acid"
    LIPID = "lipid"
    SOLVENT = "solvent"
    ION = "ion"
    UNKNOWN = "unknown"


class SerializationMode(StrEnum):
    FULL = "full"
    COMPACT = "compact"


class FrozenModel(BaseModel):
    """Base model with frozen configuration to prevent modification.

    Extra attributes are forbidden for more safety.
    """

    model_config = ConfigDict(
        frozen=True,
        extra="forbid",  # forbid extra attributes
    )


class HasAtoms(Protocol):
    """Protocol for models that have an AtomGroup."""

    atoms: AtomGroup

    @property
    def number_of_atoms(self) -> int:
        """Returns the number of atoms."""
        ...

    @property
    def number_of_residues(self) -> int:
        """Returns the number of residues."""
        ...


class FrozenWithAtoms(FrozenModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    atoms: AtomGroup

    @field_serializer("atoms")
    def serialize_atoms(self, value: AtomGroup) -> list[int]:
        """Serializes the AtomGroup to a list of atom indices."""
        return value.indices.tolist()

    @computed_field
    @property
    def number_of_atoms(self) -> int:
        """Returns the number of atoms in the small molecule."""
        return len(self.atoms)

    @computed_field
    @property
    def number_of_residues(self) -> int:
        """Returns the number of residues in the small molecule."""
        return len(self.atoms.residues)

    def _get_serialization_mode(self, info: SerializationInfo) -> SerializationMode:
        """Determines the serialization mode based on the context."""
        if info.context and "serialization_mode" in info.context:
            return SerializationMode(info.context["serialization_mode"])
        return SerializationMode.FULL

    @model_serializer(mode="wrap")
    def serialize(self, handler: SerializerFunctionWrapHandler, info: SerializationInfo) -> dict:
        """Serializes the small molecule, optionally excluding the atoms.

        Example:
            >>> structure_file = Path("tests/data/barstar.gro")
            >>> decoded = decode_structure(structure_file)
            >>> seg = decoded.inventory.segments[0]
            >>> seg.model_dump()
            {"atoms": [0, 1, 2, ...], "sequence": "KKAVI...", "number_of_atoms": 1434, "number_of_residues": 89}
            >>> seg.model_dump(context={"serialization_mode": "compact"})
            {"sequence": "KKAVI...", "number_of_atoms": 1434, "number_of_residues": 89}  # does not include "atoms" key
        """
        self_data = handler(self)
        if self._get_serialization_mode(info) == SerializationMode.COMPACT:
            # In compact mode, we exclude the atoms from the serialized output.
            del self_data["atoms"]
        return self_data


class MolecularTypeMixin(BaseModel):
    molecular_type: MolecularType

    def is_ion(self) -> bool:
        return self.molecular_type == MolecularType.ION

    def is_lipid(self) -> bool:
        return self.molecular_type == MolecularType.LIPID

    def is_solvent(self) -> bool:
        return self.molecular_type == MolecularType.SOLVENT

    def is_unknown(self) -> bool:
        return self.molecular_type == MolecularType.UNKNOWN

    def is_other(self) -> bool:
        """Alias for MolecularTypeMixin.is_unknown()"""
        return self.is_unknown()

    def is_protein(self) -> bool:
        return self.molecular_type == MolecularType.PROTEIN

    def is_nucleic(self) -> bool:
        return self.molecular_type == MolecularType.NUCLEIC


class SmallMolecule(MolecularTypeMixin, FrozenWithAtoms):
    """Small molecules are defined as residues with a single residue name."""

    description: str

    @computed_field
    @property
    def name(self) -> str:
        """Returns the name of the residue."""
        return self.atoms.residues[0].resname


class Segment(MolecularTypeMixin, FrozenWithAtoms):
    """A segment is a group of atoms that are connected."""

    @computed_field
    @property
    def sequence(self) -> str:
        return toputils.sequence(self.atoms)


class Inventory(FrozenModel):
    """An inventory of molecules in a universe."""

    small_molecules: list[SmallMolecule]
    segments: list[Segment]
    total_number_of_atoms: int


class Decoded(FrozenModel):
    """A decoded structure with its inventory."""

    inventory: Inventory
    resolution: MolecularResolution


from .settings import Settings


class GrodecoderRunOutput(BaseModel):
    """Output model for grodecoder results.

    This class is used to store data which will be part of grodecoder json output file.

    `runtime_in_seconds`.
    The program run time in seconds, that appears in grodecoder output files, is added to the json
    model by hand, at the last moment, for better accuracy.

    Comparison.
    To compare two grodecoder runs, only the `decoded` attribute.
    Other attributes, such as database version, grodecoder version, run time, etc. are ignored.
    """

    decoded: Decoded = Field(frozen=True)
    structure_file_checksum: str
    database_version: str
    grodecoder_version: str
    input_settings: Settings


# =========================================================================================================
# Models used for reading back decoded data from JSON files or APIs.


class BaseModelWithAtomsRead(FrozenModel):
    atoms: list[int]
    number_of_atoms: int
    number_of_residues: int

    @model_validator(mode="after")
    def check_number_of_atoms_is_valid(self):
        if len(self.atoms) != self.number_of_atoms:
            raise ValueError(
                f"field `number_of_atoms` ({self.number_of_atoms}) does not match number of atom ids ({len(self.atoms)})"
            )
        return self


class SmallMoleculeRead(MolecularTypeMixin, BaseModelWithAtomsRead):
    name: str
    description: str


class SegmentRead(MolecularTypeMixin, BaseModelWithAtomsRead):
    sequence: str


class InventoryRead(FrozenModel):
    small_molecules: list[SmallMoleculeRead]
    segments: list[SegmentRead]
    total_number_of_atoms: int

    @model_validator(mode="after")
    def check_total_number_of_atoms(self):
        n = sum(item.number_of_atoms for item in self.small_molecules + self.segments)
        if self.total_number_of_atoms != n:
            raise ValueError(
                f"field `total_number_of_atoms` ({self.total_number_of_atoms}) does not add up with the rest of the inventory (found {n} atoms)"
            )
        return self


class DecodedRead(FrozenModel):
    inventory: InventoryRead
    resolution: MolecularResolution

    @classmethod
    def from_decoded(cls, decoded: Decoded) -> DecodedRead:
        """Creates a DecodedRead instance from a Decoded instance.

        Useful to compare a json-loaded DecodedRead with an in-memory Decoded instance
        (e.g. for regression tests).
        """
        small_molecules_read = [
            SmallMoleculeRead(
                name=sm.name,
                description=sm.description,
                molecular_type=sm.molecular_type,
                number_of_atoms=sm.number_of_atoms,
                number_of_residues=sm.number_of_residues,
                atoms=sm.atoms.indices.tolist(),
            )
            for sm in decoded.inventory.small_molecules
        ]

        segments_read = [
            SegmentRead(
                sequence=seg.sequence,
                molecular_type=seg.molecular_type,
                atoms=seg.atoms.indices.tolist(),
                number_of_atoms=seg.number_of_atoms,
                number_of_residues=seg.number_of_residues,
            )
            for seg in decoded.inventory.segments
        ]

        inventory_read = InventoryRead(
            small_molecules=small_molecules_read,
            segments=segments_read,
            total_number_of_atoms=decoded.inventory.total_number_of_atoms,
        )

        return cls(
            inventory=inventory_read,
            resolution=decoded.resolution,
        )


class GrodecoderRunOutputRead(FrozenModel):
    decoded: DecodedRead
    structure_file_checksum: str
    database_version: str
    grodecoder_version: str
    runtime_in_seconds: float
