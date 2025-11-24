from functools import lru_cache
from typing import ClassVar

from loguru import logger
from pydantic_settings import BaseSettings

from .models import MolecularResolution


class ResolutionDetectionSettings(BaseSettings):
    default_distance_cutoff: ClassVar[float] = 2.0
    distance_cutoff: float | None = None


from dataclasses import dataclass

@dataclass
class DistanceCutoffSetting:
    default_distance_cutoff_all_atom: ClassVar[float] = 5.0
    default_distance_cutoff_coarse_grain: ClassVar[float] = 7.0
    _user_distance_cutoff: float | None = None
    _guessed_distance_cutoff: float | None = None

    def is_defined(self) -> bool:
        """Returns True if the distance cutoff has been set or guessed."""
        return self._user_distance_cutoff or self._guessed_distance_cutoff

    def is_set(self) -> bool:
        """Returns True if the distance cutoff has been set."""
        return self._user_distance_cutoff is not None

    def is_guessed(self) -> bool:
        """Returns True if the distance cutoff has been guessed."""
        return self._guessed_distance_cutoff is not None

    def get(self) -> float:
        if not self.is_defined():
            raise ValueError("`distance_cutoff` must be set or guessed before it is used.")
        return self._user_distance_cutoff or self._guessed_distance_cutoff

    def set(self, value: float):
        if self.is_guessed():
            self._guessed_distance_cutoff = None
        self._user_distance_cutoff = value

    def guess(self, resolution: MolecularResolution):
        if resolution == MolecularResolution.ALL_ATOM:
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


class ChainDetectionSettings(BaseSettings):
    _distance_cutoff: DistanceCutoffSetting = DistanceCutoffSetting()

    @property
    def distance_cutoff(self) -> DistanceCutoffSetting:
        return self._distance_cutoff

    @distance_cutoff.setter
    def distance_cutoff(self, value: float):
        self._distance_cutoff.set(value)


class Settings(BaseSettings):
    resolution_detection: ResolutionDetectionSettings = ResolutionDetectionSettings()
    chain_detection: ChainDetectionSettings = ChainDetectionSettings()

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
