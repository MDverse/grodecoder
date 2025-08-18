import json
import sys
from pathlib import Path

import click
from loguru import logger

import grodecoder as gd


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


def main(topology_path: Path, output_to_stdout: bool) -> dict:
    """Main function to process a topology file and count the molecules."""
    logger.info(f"Processing topology file: {topology_path}")

    # Decoding.
    json_data = gd.decode_topology(topology_path).dump_json() 

    # Output results: to stdout or writes to a file.
    if output_to_stdout:
        print(json.dumps(json_data, indent=2))
    else:
        output_file = topology_path.with_suffix(".json").name
        with open(output_file, "w") as f:
            json.dump(json_data, f, indent=2)
        logger.info(f"Results written to {output_file}")


@click.command()
@click.argument("topology_path", type=PathToTopologyFile())
@click.option("--stdout", is_flag=True, help="Output the results to stdout in JSON format")
@click.option("--debug", is_flag=True, help="Enable debug mode for detailed logging")
def cli(topology_path, stdout, debug):
    """Command-line interface for processing topology files."""
    setup_logging(debug)
    main(topology_path, stdout)


if __name__ == "__main__":
    main()
