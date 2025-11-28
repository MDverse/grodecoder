import click

from ..main import main as grodecoder_main
from .args import Arguments as CliArgs
from .args import CoordinatesFile, StructureFile
from ..logging import setup_logging


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
@click.option("-v", "--verbose", is_flag=True, help="show debug messages")
def cli(**kwargs):
    """Command-line interface for processing structure files."""
    args = CliArgs(
        structure_file=kwargs["structure_file"],
        coordinates_file=kwargs["coordinates_file"],
        no_atom_ids=kwargs["no_atom_ids"],
        print_to_stdout=kwargs["stdout"],
    )

    logfile = args.get_log_filename()
    setup_logging(logfile, kwargs["verbose"])
    grodecoder_main(args)


__all__ = ["cli"]
