"""Creates the MAD and CSLM databases for grodecoder.

At the moment these databases are used for identifying lipids.
"""

import json
import sys
from pathlib import Path

import click
from loguru import logger

from grodecoder.databases.mad.web_crawler import fetch as fetch_mad
from grodecoder.databases.csml.web_crawler import fetch as fetch_csml


def build_mad_database(output_path: Path):
    """Creates the MAD database for grodecoder by scrapping the source available at mad.ibcp.fr."""
    logger.info("Building MAD database...")
    db = fetch_mad()
    with open(output_path, "w") as f:
        json.dump([residue.model_dump() for residue in db], f, indent=2)
    logger.info(f"Wrote MAD database with {len(db)} elements to {output_path}")


def build_csml_database(output_path: Path):
    """Creates the CHARMM database for grodecoder by scrapping the source available at charmm-gui.org."""
    logger.info("Building CHARMM database...")
    db = fetch_csml()
    with open(output_path, "w") as f:
        json.dump([residue.model_dump() for residue in db], f, indent=2)
    logger.info(f"Wrote CHARMM database with {len(db)} elements to {output_path}")


def fatal_file_exists(path: Path):
    msg = f"Output file {path} already exists. Please remove it or choose a different name."
    if path.exists():
        click.secho(msg, err=True, fg="red")
        sys.exit(1)


@click.command()
def main():
    """Creates the CHARMM database for grodecoder by scrapping the source available at charmm-gui.org."""

    mad_output = Path("mad_database.json")
    csml_output = Path("charmm_csml_database.json")

    if mad_output.exists():
        fatal_file_exists(mad_output)

    if csml_output.exists():
        fatal_file_exists(csml_output)

    build_mad_database(mad_output)
    build_csml_database(csml_output)


if __name__ == "__main__":
    main()
