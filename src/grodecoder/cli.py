import sys
from pathlib import Path

import click
from loguru import logger

import grodecoder as gd

GRODECODER_STEM_SUFFIX = "_grodecoder"


class PathToStructureFile(click.ParamType):
    """Custom click parameter type for validating structure files."""

    name = "structure file"

    def convert(self, value, param, ctx):
        """Convert the input value to a Path object."""
        path = Path(value)
        if not path.exists():
            self.fail(f"'{path}' does not exist", param, ctx)
        if not path.is_file():
            self.fail(f"'{path}' is not a file", param, ctx)
        valid_extensions = (".gro", ".pdb", ".coor", ".crd")
        if path.suffix not in valid_extensions:
            self.fail(
                f"'{path}' has an invalid extension (valid extensions are {valid_extensions})", param, ctx
            )
        return path


def get_output_stem(structure_path: Path) -> str:
    """Generates the output stem for log and result files based on the structure file name."""
    return structure_path.stem + GRODECODER_STEM_SUFFIX


def generate_output_filename(structure_path: Path, extension: str) -> Path:
    """Generates the output filename with the given extension."""
    return Path(get_output_stem(structure_path) + extension)


def generate_output_inventory_filename(structure_path: Path) -> Path:
    """Generates the output JSON filename."""
    return generate_output_filename(structure_path, ".json")


def generate_log_filename(structure_path: Path) -> Path:
    """Generates the output log filename."""
    return generate_output_filename(structure_path, ".log")


def setup_logging(logfile: Path, debug: bool = False):
    """Sets up logging configuration."""
    fmt = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> <level>{level}: {message}</level>"
    level = "DEBUG" if debug else "INFO"
    logger.remove()
    logger.add(sys.stderr, level=level, format=fmt, colorize=True)
    logger.add(logfile, level=level, format=fmt, colorize=False, mode="w")


def main(structure_path: Path, bond_threshold: float, compact_serialization: bool, output_to_stdout: bool):
    """Main function to process a structure file and count the molecules.

    Args:
        structure_path (Path): Path to the structure file.
        bond_threshold (float): Threshold for interchain bond detection.
        compact_serialization (bool): If True, use compact serialization (no atom indices).
        output_to_stdout (bool): Whether to output results to stdout.
    """
    logger.info(f"Processing structure file: {structure_path}")

    # Decoding.
    decoded = gd.decode_structure(structure_path, bond_threshold=bond_threshold)

    # Serialization.
    serialization_mode = "compact" if compact_serialization else "full"
    json_string = decoded.model_dump_json(indent=2, context={"serialization_mode": serialization_mode})

    # Output results: to stdout or writes to a file.
    if output_to_stdout:
        print(json_string)
    else:
        output_file = generate_output_inventory_filename(structure_path)
        with open(output_file, "w") as f:
            f.write(json_string)
        logger.info(f"Results written to {output_file}")


@click.command()
@click.argument("structure_path", type=PathToStructureFile())
@click.option(
    "--bond-threshold",
    default=5.0,
    type=float,
    help="Threshold for interchain bond detection (default: 5 Å)",
)
@click.option("--no-atom-ids", is_flag=True, help="do not output the atom indice array")
@click.option("-s", "--stdout", is_flag=True, help="Output the results to stdout in JSON format")
def cli(structure_path, bond_threshold, no_atom_ids, stdout):
    """Command-line interface for processing structure files."""
    logfile = generate_log_filename(structure_path)
    setup_logging(logfile)
    main(structure_path, bond_threshold, no_atom_ids, stdout)
