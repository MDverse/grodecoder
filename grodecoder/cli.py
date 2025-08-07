import argparse
from pathlib import Path

import grodecoder


class NotATopologyFileError(argparse.ArgumentTypeError):
    """Raised when the given file is not a topology (based on its extension)."""


def is_a_structure_file(filepath: str) -> Path:
    """Check if the given filepath points to an existing structure file (.gro, .pdb).

    Parameters
    ----------
        filepath : str
            Path of the file.

    Raises
    ------
        argparse.ArgumentTypeError
            If the given filepath is not an existing file,
            or if it does not have a '.gro' or '.pdb' extension.

    Returns
    -------
        str
            The validated path.
    """
    path = Path(filepath)

    if not path.exists():
        raise argparse.ArgumentTypeError(f"'{path}' does not exist")

    if not path.is_file():
        raise argparse.ArgumentTypeError(f"'{path}' is not a file")

    valid_extensions = (".gro", ".pdb", ".coor", ".crd")
    if path.suffix not in valid_extensions:
        raise NotATopologyFileError(f"'{path!s}': invalid extension '{path.suffix}' (valid extensions are {valid_extensions})")

    return path


def is_a_valid_threshold(threshold: str) -> str | float:
    """Check if the given threshold is permitted).

    Parameters
    ----------
        threshold : str | float
            "auto" or a positive float number.

    Raises
    ------
        argparse.ArgumentTypeError: If the given threshold is not "auto" or a positive number.

    Returns
    -------
        str | float
            "auto" or the validated threshold as a float number.
    """
    if threshold == "auto":
        return threshold
    try:
        threshold_as_float = float(threshold)
    except ValueError as _e:
        raise argparse.ArgumentTypeError("Argument should 'auto' or a positive number")
    if not threshold_as_float > 0.0:
        raise argparse.ArgumentTypeError("Argument should be > 0")
    return threshold_as_float


def parse_arg() -> argparse.Namespace:
    """Parse command-line arguments.

    This function uses the argparse module to parse the command-line arguments
    provided by the user. It sets up the argument parser with information about
    the program, its usage, and the available options.

    Ressources
    ----------
    - https://docs.python.org/3/library/argparse.html

    Return
    ------
        argparse.Namespace
            An object containing the parsed arguments as attributes.

    """
    parser = argparse.ArgumentParser(
        prog="grodecoder",
        description="Extract molecules from a structure file (.gro, .pdb).",
    )
    parser.add_argument(
        "--topology",
        type=is_a_structure_file,
        required=True,
        help="path to topology file",
    )
    parser.add_argument(
        "--psf",
        type=Path,
        help="topology file path (.psf)",
    )
    parser.add_argument(
        "--checkconnectivity",
        help="Add edges and degre in the fingerprint. Default: False.",
        action="store_true",
    )
    parser.add_argument(
        "--bondthreshold",
        type=is_a_valid_threshold,
        help="Choose the method to calculate the atom pairs. If we know the resolution of the system is coarse grain enter a threshold (a positiv float number) or we don't know so choose 'auto'",
        default="auto",
    )
    parser.add_argument(
        "--querypdb",
        help="Add PDB id and their putative name in the JSON file for the protein. Default: False.",
        action="store_true",
    )
    args = parser.parse_args()

    if args.topology.suffix == ".coor" and not args.psf:
        parser.error("The '--psf' argument is required when the topology file is a .coor file.")

    if args.psf and args.topology.suffix != ".coor":
        parser.error("The '--psf' argument is valid only when the topology file is a .coor file.")

    return args


def main():
    raise NotImplementedError
    # args = parse_arg()
    # grodecoder.v1.main(
    #     args.topology,
    #     args.psf,
    #     check_connectivity=args.checkconnectivity,
    #     bond_threshold=args.bondthreshold,
    #     query_pdb=args.querypdb,
    # )



if __name__ == "__main__":
    main()
