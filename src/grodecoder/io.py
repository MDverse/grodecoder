"""Grodecoder read/write functions"""

import json
import MDAnalysis as mda

from loguru import logger

from ._typing import PathLike
from .models import GrodecoderRunOutputRead


def read_universe(structure_path: PathLike, coordinates_path: PathLike | None = None) -> mda.Universe:
    """Reads a structure file."""
    source = (structure_path, coordinates_path) if coordinates_path else (structure_path, )
    source_str = ", ".join(str(s) for s in source)
    logger.debug(f"Reading universe from {source_str}")
    universe: mda.Universe | None = None
    try:
        universe = mda.Universe(*source)
    except Exception as e:
        raise IOError("MDAnalysis error while reading universe") from e

    if not hasattr(universe, "trajectory"):
        raise mda.exceptions.NoDataError(f"no coordinates read from {source_str}")

    return universe


def read_grodecoder_output(path: PathLike) -> GrodecoderRunOutputRead:
    with open(path) as fileobj:
        return GrodecoderRunOutputRead.model_validate(json.load(fileobj))
