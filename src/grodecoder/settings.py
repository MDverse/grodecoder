from dataclasses import dataclass
from functools import lru_cache
from typing import ClassVar, TYPE_CHECKING

from loguru import logger
from pydantic import ConfigDict, Field, GetJsonSchemaHandler
from pydantic.json_schema import JsonSchemaValue
from pydantic_core import core_schema
from pydantic_settings import BaseSettings
from typing_extensions import Annotated

if TYPE_CHECKING:
    from .models import MolecularResolution


@dataclass(init=False)
class DistanceCutoff:
    default_cutoff_distance_all_atom: ClassVar[float] = 5.0
    default_cutoff_distance_coarse_grain: ClassVar[float] = 6.0
    _user_cutoff_distance: float | None = None
    _guessed_cutoff_distance: float | None = None

    def __init__(self, user_value: float | None = None):
        if user_value is not None:
            self.set(user_value)

    def is_defined(self) -> bool:
        """Returns True if the distance cutoff has been set or guessed."""
        return self._user_cutoff_distance is not None or self._guessed_cutoff_distance is not None

    def is_set(self) -> bool:
        """Returns True if the distance cutoff has been set."""
        return self._user_cutoff_distance is not None

    def is_guessed(self) -> bool:
        """Returns True if the distance cutoff has been guessed."""
        return self._guessed_cutoff_distance is not None

    def get(self) -> float:
        if not self.is_defined():
            raise ValueError("`cutoff_distance` must be set or guessed before it is used.")
        return self._user_cutoff_distance or self._guessed_cutoff_distance  # ty: ignore[invalid-return-type]

    def set(self, value: float):
        if self.is_guessed():
            self._guessed_cutoff_distance = None
        self._user_cutoff_distance = value

    def guess(self, resolution: "MolecularResolution"):
        if resolution.is_all_atom():
            cutoff_distance = self.default_cutoff_distance_all_atom
            logger.debug(
                f"chain detection: using default distance cutoff for all atom structures: {cutoff_distance:.2f}"
            )
        else:
            cutoff_distance = self.default_cutoff_distance_coarse_grain
            logger.debug(
                f"chain detection: using default distance cutoff for coarse grain structures: {cutoff_distance:.2f}"
            )
        self._guessed_cutoff_distance = cutoff_distance


class _DistanceCutoffPydanticAnnotation:
    """Allows to serialize / validate DistanceCutoff using pydantic.


    Examples:
    >>> from grodecoder.settings import ChainDetectionSettings
    >>> cds = ChainDetectionSettings()
    >>> cds.cutoff_distance
    DistanceCutoff(_user_cutoff_distance=None, _guessed_cutoff_distance=None)

    >>> # Float assignement
    >>> cds.cutoff_distance = 12
    >>> cds.cutoff_distance
    DistanceCutoff(_user_cutoff_distance=12.0, _guessed_cutoff_distance=None)

    >>> # None assignment
    >>> cds.cutoff_distance = None
    >>> cds.cutoff_distance
    DistanceCutoff(_user_cutoff_distance=None, _guessed_cutoff_distance=None)

    >>> # Serialization
    >>> cds.cutoff_distance = 12
    >>> cds.model_dump()
    {'cutoff_distance': 12.0}

    >>> # Validation
    >>> as_json = cds.model_dump()
    >>> ChainDetectionSettings.model_validate(as_json)
    ChainDetectionSettings(cutoff_distance=DistanceCutoff(_user_cutoff_distance=12.0, _guessed_cutoff_distance=None))
    """

    @classmethod
    def __get_pydantic_core_schema__(cls, _source_type, _handler) -> core_schema.CoreSchema:
        """
        We return a pydantic_core.CoreSchema that behaves in the following ways:

        * floats will be parsed as `DistanceCutoff` instances with the float as the `_user_cutoff_distance` attribute
        * `DistanceCutoff` instances will be parsed as `DistanceCutoff` instances without any changes
        * Nothing else will pass validation
        * Serialization will always return just a float
        """

        def validate_from_none(value: None) -> DistanceCutoff:
            return DistanceCutoff()

        def validate_from_float(value: float) -> DistanceCutoff:
            result = DistanceCutoff()
            result.set(value)
            return result

        from_none_schema = core_schema.chain_schema(
            [
                core_schema.none_schema(),
                core_schema.no_info_plain_validator_function(validate_from_none),
            ]
        )
        from_float_schema = core_schema.chain_schema(
            [
                core_schema.float_schema(),
                core_schema.no_info_plain_validator_function(validate_from_float),
            ]
        )
        return core_schema.json_or_python_schema(
            json_schema=from_float_schema,
            python_schema=core_schema.union_schema(
                [
                    # check if it's an instance first before doing any further work
                    core_schema.is_instance_schema(DistanceCutoff),
                    from_none_schema,
                    from_float_schema,
                ]
            ),
            serialization=core_schema.plain_serializer_function_ser_schema(lambda instance: instance.get()),
        )

    @classmethod
    def __get_pydantic_json_schema__(
        cls, _core_schema: core_schema.CoreSchema, handler: GetJsonSchemaHandler
    ) -> JsonSchemaValue:
        return handler(core_schema.float_schema())


# We now create an `Annotated` wrapper that we'll use as the annotation for fields on `BaseModel`s, etc.
PydanticDistanceCutoff = Annotated[DistanceCutoff, _DistanceCutoffPydanticAnnotation]


class ChainDetectionSettings(BaseSettings):
    model_config = ConfigDict(validate_assignment=True)
    cutoff_distance: PydanticDistanceCutoff = Field(default_factory=DistanceCutoff)


class ResolutionDetectionSettings(BaseSettings):
    cutoff_distance: float = 1.6


class OutputSettings(BaseSettings):
    # should we output atom ids?
    atom_ids: bool = True


class Settings(BaseSettings):
    resolution_detection: ResolutionDetectionSettings = ResolutionDetectionSettings()
    chain_detection: ChainDetectionSettings = ChainDetectionSettings()
    output: OutputSettings = OutputSettings()

    debug: bool = False


_settings: Settings | None = None


@lru_cache()
def get_settings():
    global _settings
    if _settings is None:
        _settings = Settings()
    return _settings


def get_chain_detection_settings():
    return get_settings().chain_detection
