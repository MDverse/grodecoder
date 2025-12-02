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
    default_distance_cutoff_all_atom: ClassVar[float] = 5.0
    default_distance_cutoff_coarse_grain: ClassVar[float] = 6.0
    _user_distance_cutoff: float | None = None
    _guessed_distance_cutoff: float | None = None

    def __init__(self, user_value: float | None = None):
        if user_value is not None:
            self.set(user_value)

    def is_defined(self) -> bool:
        """Returns True if the distance cutoff has been set or guessed."""
        return any((self._user_distance_cutoff, self._guessed_distance_cutoff))

    def is_set(self) -> bool:
        """Returns True if the distance cutoff has been set."""
        return self._user_distance_cutoff is not None

    def is_guessed(self) -> bool:
        """Returns True if the distance cutoff has been guessed."""
        return self._guessed_distance_cutoff is not None

    def get(self) -> float:
        if not self.is_defined():
            raise ValueError("`distance_cutoff` must be set or guessed before it is used.")
        return self._user_distance_cutoff or self._guessed_distance_cutoff  # ty: ignore[invalid-return-type]

    def set(self, value: float):
        if self.is_guessed():
            self._guessed_distance_cutoff = None
        self._user_distance_cutoff = value

    def guess(self, resolution: "MolecularResolution"):
        if resolution == "ALL_ATOM":
            distance_cutoff = self.default_distance_cutoff_all_atom
            logger.debug(
                f"chain detection: using default distance cutoff for all atom structures: {distance_cutoff:.2f}"
            )
        else:
            distance_cutoff = self.default_distance_cutoff_coarse_grain
            logger.debug(
                f"chain detection: using default distance cutoff for coarse grain structures: {distance_cutoff:.2f}"
            )
        self._guessed_distance_cutoff = distance_cutoff


class _DistanceCutoffPydanticAnnotation:
    """Allows to serialize / validate DistanceCutoff using pydantic.


    Examples:
    >>> from grodecoder.settings import ChainDetectionSettings
    >>> cds = ChainDetectionSettings()
    >>> cds.distance_cutoff
    DistanceCutoff(_user_distance_cutoff=None, _guessed_distance_cutoff=None)

    >>> # Float assignement
    >>> cds.distance_cutoff = 12
    >>> cds.distance_cutoff
    DistanceCutoff(_user_distance_cutoff=12.0, _guessed_distance_cutoff=None)

    >>> # None assignment
    >>> cds.distance_cutoff = None
    >>> cds.distance_cutoff
    DistanceCutoff(_user_distance_cutoff=None, _guessed_distance_cutoff=None)

    >>> # Serialization
    >>> cds.distance_cutoff = 12
    >>> cds.model_dump()
    {'distance_cutoff': 12.0}

    >>> # Validation
    >>> as_json = cds.model_dump()
    >>> ChainDetectionSettings.model_validate(as_json)
    ChainDetectionSettings(distance_cutoff=DistanceCutoff(_user_distance_cutoff=12.0, _guessed_distance_cutoff=None))
    """

    @classmethod
    def __get_pydantic_core_schema__(cls, _source_type, _handler) -> core_schema.CoreSchema:
        """
        We return a pydantic_core.CoreSchema that behaves in the following ways:

        * floats will be parsed as `DistanceCutoff` instances with the float as the `_user_distance_cutoff` attribute
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
    distance_cutoff: PydanticDistanceCutoff = Field(default_factory=DistanceCutoff)


class ResolutionDetectionSettings(BaseSettings):
    distance_cutoff: float = 1.6


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
