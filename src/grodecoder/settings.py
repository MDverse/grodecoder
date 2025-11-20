from typing import ClassVar

from pydantic_settings import BaseSettings


class ResolutionDetectionSettings(BaseSettings):
    default_distance_cutoff: ClassVar[float] = 2.0
    distance_cutoff: float | None = None


class ChainDetectionSettings(BaseSettings):
    default_distance_cutoff: ClassVar[float] = 5.0
    distance_cutoff: float | None = None


class Settings(BaseSettings):
    resolution_detection: ResolutionDetectionSettings = ResolutionDetectionSettings()
    chain_detection: ChainDetectionSettings = ChainDetectionSettings()

    debug: bool = False
