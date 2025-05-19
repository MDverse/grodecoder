"""Creates the database of molecules based on CHARMM data."""

import json
import sys
from pathlib import Path

import click

import grodecoder


@click.command()
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default="charmm_csml_database.json",
    help="Output file for the CHARMM database.",
)
def main(output: Path):
    """Creates the CHARMM database for grodecoder by scrapping the source available at charmm-gui.org."""

    if output.exists():
        click.secho(
            f"Output file {output} already exists. Please remove it or choose a different name.",
            err=True,
            fg="red",
        )
        sys.exit(1)

    db = grodecoder.databases.csml.fetch()

    # Save the data to a JSON file.
    with open(output, "w") as f:
        json.dump([residue.model_dump() for residue in db], f, indent=2)


if __name__ == "__main__":
    main()
