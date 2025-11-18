import sys
import warnings
from pathlib import Path

import click
from loguru import logger

from ..main import main as grodecoder_main
from .args import Arguments as CliArgs
from .args import CoordinatesFile, StructureFile


def setup_logging(logfile: Path, debug: bool = False):
    """Sets up logging configuration."""
    fmt = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> <level>{level}: {message}</level>"
    level = "DEBUG" if debug else "INFO"
    logger.remove()
    logger.add(sys.stderr, level=level, format=fmt, colorize=True)
    logger.add(logfile, level=level, format=fmt, colorize=False, mode="w")

    # Setup loguru to capture warnings (typically MDAnalysis warnings)
    def showwarning(message, *args, **kwargs):
        logger.opt(depth=2).warning(message)

    warnings.showwarning = showwarning

@click.command()
@click.argument("structure_file", type=StructureFile)
@click.argument("coordinates_file", type=CoordinatesFile, required=False)
@click.option(
    "--bond-threshold",
    default=5.0,
    type=float,
    help="Threshold for interchain bond detection (default: 5 Å)",
)
@click.option("--no-atom-ids", is_flag=True, help="do not output the atom indice array")
@click.option(
    "-s",
    "--stdout",
    metavar="print_to_stdout",
    is_flag=True,
    help="Output the results to stdout in JSON format",
)
def cli(**kwargs):
    """Command-line interface for processing structure files."""
    args = CliArgs(
        structure_file=kwargs["structure_file"],
        coordinates_file=kwargs["coordinates_file"],
        no_atom_ids=kwargs["no_atom_ids"],
        print_to_stdout=kwargs["stdout"],
    )

    logfile = args.get_log_filename()
    setup_logging(logfile)
    grodecoder_main(args)


__all__ = ["cli"]
