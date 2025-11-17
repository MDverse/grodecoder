import json

from . import databases, toputils
from ._typing import AtomGroup, Json, PathLike, Residue, Universe, UniverseLike
from .core import decode, decode_structure
from .identifier import identify
from .models import Decoded, GrodecoderRunOutput, GrodecoderRunOutputRead

__all__ = [
    "databases",
    "decode",
    "decode_structure",
    "identify",
    "toputils",
    "read_structure",
    "AtomGroup",
    "Decoded",
    "GrodecoderRunOutput",
    "GrodecoderRunOutputRead",
    "Json",
    "PathLike",
    "Residue",
    "Universe",
    "UniverseLike",
]


def read_grodecoder_output(path: PathLike) -> GrodecoderRunOutputRead:
    with open(path) as fileobj:
        return GrodecoderRunOutputRead.model_validate(json.load(fileobj))
