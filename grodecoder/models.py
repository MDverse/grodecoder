from __future__ import annotations

from enum import StrEnum

from MDAnalysis import AtomGroup
from pydantic import (
    BaseModel,
    ConfigDict,
    SerializationInfo,
    SerializerFunctionWrapHandler,
    computed_field,
    field_serializer,
    model_serializer,
)


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
    """Base model with frozen configuration to prevent modification."""

    model_config = ConfigDict(frozen=True)


class FrozenWithAtoms(FrozenModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    atoms: AtomGroup

    @field_serializer("atoms")
    def serialize_atoms(self, value: AtomGroup) -> list[int]:
        """Serializes the AtomGroup to a list of atom indices."""
        return value.indices.tolist()

    @computed_field
    def number_of_atoms(self) -> int:
        """Returns the number of atoms in the small molecule."""
        return len(self.atoms)

    @computed_field
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
            >>> topology_file = Path("tests/data/barstar.gro")
            >>> decoded = decode_topology(topology_file)
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


class SmallMolecule(FrozenWithAtoms):
    """Small molecules are defined as residues with a single residue name."""

    description: str
    molecular_type: MolecularType

    @computed_field  # type: ignore[prop-decorator]
    @property
    def name(self) -> str:
        """Returns the name of the residue."""
        return self.atoms.residues[0].resname


class Segment(FrozenWithAtoms):
    """A segment is a group of atoms that are connected."""

    sequence: str
    molecular_type: MolecularType  # likely to be protein or nucleic acid


class Inventory(FrozenModel):
    """An inventory of molecules in a universe."""

    small_molecules: list[SmallMolecule]
    segments: list[Segment]


class Decoded(FrozenModel):
    """A decoded structure with its inventory."""

    inventory: Inventory
    resolution: MolecularResolution
    database_version: str


# =========================================================================================================
# Models used for reading back decoded data from JSON files or APIs.


class BaseModelWithAtomsRead(FrozenModel):
    atoms: list[int]
    number_of_atoms: int
    number_of_residues: int


class SmallMoleculeRead(BaseModelWithAtomsRead):
    name: str
    description: str
    molecular_type: MolecularType


class SegmentRead(BaseModelWithAtomsRead):
    sequence: str
    molecular_type: MolecularType


class InventoryRead(FrozenModel):
    small_molecules: list[SmallMoleculeRead]
    segments: list[SegmentRead]


class DecodedRead(FrozenModel):
    inventory: InventoryRead
    resolution: MolecularResolution

    @classmethod
    def from_decoded(cls, decoded: Decoded) -> DecodedRead:
        """Creates a DecodedRead instance from a Decoded instance."""
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

        inventory_read = InventoryRead(small_molecules=small_molecules_read, segments=segments_read)

        return cls(inventory=inventory_read, resolution=decoded.resolution)
