import hashlib
import json
import time
from typing import TYPE_CHECKING

from loguru import logger

from ._typing import PathLike
from .core import decode_structure
from .databases import get_database_version
from .models import GrodecoderRunOutput
from .version import get_version
from .settings import Settings

if TYPE_CHECKING:
    from .cli.args import Arguments as CliArgs


def _get_checksum(structure_path: PathLike) -> str:
    """Computes a checksum for the structure file."""
    with open(structure_path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def main(args: "CliArgs"):
    """Main function to process a structure file and count the molecules."""
    start_time = time.perf_counter_ns()
    structure_path = args.structure_file.path
    coordinates_path = args.coordinates_file.path if args.coordinates_file else None

    settings = Settings()

    logger.info(f"Processing structure file: {structure_path}")

    # Decoding.
    decoded = decode_structure(
        structure_path, coordinates_path=coordinates_path, settings = settings
    )

    output = GrodecoderRunOutput(
        decoded=decoded,
        structure_file_checksum=_get_checksum(structure_path),
        database_version=get_database_version(),
        grodecoder_version=get_version(),
    )

    # Serialization.
    serialization_mode = "compact" if args.no_atom_ids else "full"

    # Updates run time as late as possible.
    output_json = output.model_dump(context={"serialization_mode": serialization_mode})
    output_json["runtime_in_seconds"] = (time.perf_counter_ns() - start_time) / 1e9

    # Output results: to stdout or writes to a file.
    if args.print_to_stdout:
        print(json.dumps(output_json, indent=2))
    else:
        inventory_filename = args.get_inventory_filename()
        with open(inventory_filename, "w") as f:
            f.write(json.dumps(output_json, indent=2))
        logger.info(f"Results written to {inventory_filename}")
