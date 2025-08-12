import sys
import time
from pathlib import Path

import click
from loguru import logger

import grodecoder as gd
from grodecoder.models import Inventory

# DEBUG
import icecream

icecream.install()



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


def test():
    data_root_dir = Path("data/examples")
    topology_files = [
        "1BRS.gro",
        "1QJ8.gro",
        "1QJ8_ETH_ACN_MET_URE_SOL.pdb",
        "1QJ8_membrane.gro",
        "1QJ8_solution.gro",
        "2MAT.pdb",
        "4MQJ_ABCD.gro",
        "4ZRY.gro",
        "5MBA.gro",
        "5ZOA.gro",
        "barstar.gro",
        "DMPC_PI.gro",
        "DNA_start.gro",
        "noriega_AA_CRD_3CAL.gro",
        "noriega_CG_CRD_3CAL.gro",
        "RNA_start.gro",
    ]

    for topology_file in topology_files:
        topology_path = data_root_dir / topology_file
        main(topology_path, debug=False)


def main(topology_path: Path, debug: bool = False) -> dict:
    """Main function to process a topology file and count the molecules."""
    setup_logging(debug)
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
@click.argument("topology_path", type=click.Path(exists=True, path_type=Path))
@click.option("--debug", is_flag=True, help="Enable debug mode for detailed logging")
def cli(topology_path, debug):
    """Command-line interface for processing topology files."""
    main(topology_path)


if __name__ == "__main__":
    # main()
    test()
