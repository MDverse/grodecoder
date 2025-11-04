"""Run decoding tests on a batch of structure files to ensure everything goes ok."""

import sys
from pathlib import Path

from icecream import ic
from loguru import logger

from grodecoder import decode_structure

logger.remove()
logger.add(
    sys.stderr,
    level="INFO",
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{message}</level>",
    colorize=True,
    backtrace=True,
    diagnose=True,
)


def test():
    data_root_dir = Path("tests/data/regression_data/input_files")
    structure_files = set(data_root_dir.rglob("*.gro")) | set(data_root_dir.rglob("*.pdb"))

    for structure_path in structure_files:
        decoded = decode_structure(structure_path)
        as_json = decoded.model_dump(context={"serialization_mode": "compact"})
        ic(as_json)


if __name__ == "__main__":
    test()
