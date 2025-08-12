import sys
import time
from pathlib import Path

import click
from loguru import logger

import grodecoder as gd
from grodecoder.models import Inventory


class PathToTopologyFile(click.ParamType):
    """Custom click parameter type for validating topology files."""

    name = "topology file"

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


def setup_logging(debug: bool = False):
    """Sets up logging configuration."""
    fmt = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{message}</level>"
    level = "DEBUG" if debug else "INFO"
    logger.remove()
    logger.add(sys.stderr, level=level, format=fmt, colorize=True, backtrace=True, diagnose=True)


def actually_count(universe: gd.UniverseLike) -> Inventory:
    """Identifies and counts the molecules in a topology file."""
    timer_start = time.perf_counter()  # do not include topology reading time in the count
    inventory = gd.identify(universe)
    elapsed = time.perf_counter() - timer_start
    logger.debug(f"{len(universe.atoms):,d} atoms processed in {elapsed:.2f} seconds")
    return inventory


def main(topology_path: Path) -> dict:
    """Main function to process a topology file and count the molecules."""
    logger.info(f"Processing topology file: {topology_path}")

    universe = gd.read_topology(topology_path)
    logger.debug(f"{topology_path}: {len(universe.atoms):,d} atoms")

    inventory = actually_count(universe)

    # Debug mode: print the inventory
    for molecule in inventory.small_molecules:
        logger.debug(f"{topology_path}: {molecule}")
    for segment in inventory.segments:
        logger.debug(f"{topology_path}: {segment}")

    # Create JSON data
    json_data = {
        "resolution": gd.toputils.guess_resolution(universe),
        "inventory": inventory,
    }
    return json_data


@click.command()
@click.argument("topology_path", type=PathToTopologyFile())
@click.option("--debug", is_flag=True, help="Enable debug mode for detailed logging")
def cli(topology_path, debug):
    """Command-line interface for processing topology files."""
    setup_logging(debug)
    main(topology_path)


if __name__ == "__main__":
    main()
