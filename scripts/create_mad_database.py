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
    default="mad_database.json",
    help="Output file for the MAD database.",
)
def main(output: Path):
    """Creates the MAD database for grodecoder by scrapping the source available at mad.ibcp.fr."""

    if output.exists():
        click.secho(
            f"Output file {output} already exists. Please remove it or choose a different name.",
            err=True,
            fg="red",
        )
        sys.exit(1)

    db = grodecoder.databases.mad.fetch()

    # Save the data to a JSON file.
    with open(output, "w") as f:
        json.dump(db, f, indent=2)


if __name__ == "__main__":
    main()
