"""Retrieves information about a specific protein sequence (given by a fasta file) from the RCSB PDB API.

Usage:
    python search_into_PDB.py --fasta [fasta input file]
"""

__authors__ = ("Karine DUONG", "Pierre POULAIN")
__contact__ = "pierre.poulain@u-paris.fr"


import argparse
from Bio import SeqIO
import json
from pathlib import Path
import requests

from loguru import logger


def API_PDB_search_based_sequence(sequence: str, max_element: int = 10) -> list[str]:
    """This function searches the Protein Data Bank (PDB) database based on a given sequence and returns a list of PDB IDs.

    Ressources
    ----------
        https://education.molssi.org/python-scripting-biochemistry/chapters/rcsb_api.html
        https://search.rcsb.org/index.html#return-count

    Parameters
    ----------
        sequence : str
            The amino acid sequence to search for in the PDB database.
        max_element : int, optional
            The maximum number of PDB IDs to retrieve. Defaults to 10.

    Returns
    -------
        list
            A list of PDB IDs corresponding to the sequences found in the PDB database.
    """
    logger.info("Searching PDB IDs based on a given sequence")
    my_query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.1,
                "identity_cutoff": 0,
                "sequence_type": "protein",
                "value": sequence,
            },
        },
        "return_type": "polymer_entity",
        # "return_type": "entry",
        "request_options": {
            "results_verbosity": "compact",
            "results_content_type": ["experimental"],
            "scoring_strategy": "combined",
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }

    try:
        my_query = json.dumps(my_query)
        data = requests.get(
            f"https://search.rcsb.org/rcsbsearch/v2/query?json={my_query}"
        )
    except requests.exceptions.HTTPError as http_err:
        logger.error(f"HTTP error occurred: {http_err}")
        return []
    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request error occurred: {req_err}")
        return []

    try:
        results = data.json()
        return results["result_set"][:max_element]
    except requests.exceptions.JSONDecodeError as json_err:
        logger.error(f"JSON decode error occurred: {json_err}")
        return []
    except KeyError as key_err:
        logger.error(f"Key error occurred: {key_err}")
        return []


def get_macromolecular_names(PDB_ID, polymer_entity_id, set_names) -> set[str]:
    """Retrieve macromolecular names associated with a given PDB ID and polymer entity ID.

    References
    ----------
    https://www.rcsb.org/docs/programmatic-access/web-services-overview#data-api
    https://data.rcsb.org/index.html#data-api
    https://data.rcsb.org/redoc/index.html#tag/Entity-Service/operation/getPolymerEntityById

    Parameters
    ----------
        PDB_ID : str
            The PDB ID of the structure.
        polymer_entity_id : str
            The polymer entity ID within the PDB structure.
        set_names : set
            A set containing existing macromolecular names. Retrieved names will be added to this set.

    Returns
    -------
        set
            A set containing the updated list of macromolecular names.

    Raises
    ------
        requests.exceptions.RequestException
            If an error occurs while making the API request.
    """
    logger.info(
        f"Retrieving macromolecular name for PDB ID {PDB_ID} polymer {polymer_entity_id}"
    )
    try:
        my_query = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{PDB_ID}/{polymer_entity_id}"
        )
        my_query.raise_for_status()

        results = my_query.json()
        rcsb_macromolecular_names_combined = results["rcsb_polymer_entity"][
            "rcsb_macromolecular_names_combined"
        ]
        for index in range(len(rcsb_macromolecular_names_combined)):
            result = results["rcsb_polymer_entity"][
                "rcsb_macromolecular_names_combined"
            ][index]["name"].upper()
            set_names.add(result)
        return set_names

    except requests.exceptions.RequestException as err:
        print(f"Error: {err}")
        return set_names


def treat_PDB_ID_to_macromolecular_names(PDB_ID_polymer_entity_id: str) -> set[str]:
    """Retrieve macromolecular names for multiple PDB IDs and polymer entity IDs.

    Parameters
    ----------
        PDB_ID_polymer_entity_id : list of str
            A list containing strings in the format '(PDB_ID)_(polymer_entity_id)'.

    Returns
    -------
        set
            A set containing unique macromolecular names extracted from the specified PDB IDs and polymer entity IDs.
    """
    set_macromolecular_names: set[str] = set()

    for ID_entity in PDB_ID_polymer_entity_id:
        pdb_ID, polymer_entity_id = ID_entity.split("_")
        get_macromolecular_names(pdb_ID, polymer_entity_id, set_macromolecular_names)
    return set_macromolecular_names


def get_macromolecular_names_bis(my_query_json: dict[str, str]) -> list[str]:
    """Retrieve all the macromolecular names for this PDB id and this polymer id

    Args:
        my_query_json: dict
            The request from the PDB API

    Returns:
        list[str]
            All the macromolecular names for this PDB id and this polymer id
    """
    list_macromolecular_names = []
    rcsb_macromolecular_names_combined = my_query_json["rcsb_polymer_entity"][
        "rcsb_macromolecular_names_combined"
    ]
    for index in range(len(rcsb_macromolecular_names_combined)):
        list_macromolecular_names.append(
            my_query_json["rcsb_polymer_entity"]["rcsb_macromolecular_names_combined"][
                index
            ]["name"].upper()
        )
    return list_macromolecular_names


def get_organism_names(my_query_json: dict) -> list:
    """Retrieve the organism names for this PDB id and this polymer id

    Args:
        my_query_json: dict
            The request from the PDB API

    Returns:
        list[str]
            The organism names for this PDB id and this polymer id
    """
    return my_query_json["rcsb_entity_source_organism"][0]["ncbi_scientific_name"]


def get_info_one_pdb_id(PDB_ID_polymer_entity_id: str) -> dict[str, list[str]]:
    """Retrieves information about a specific PDB ID and its associated polymer entity from the PDB API.

    Args:
        PDB_ID_polymer_entity_id: str
            One PDB id retrieve from the PDB API

    Returns:
        dict[str, list[str]]
            Dictionnary that have information about this PDB id
    """
    pdb_ID, polymer_entity_id = PDB_ID_polymer_entity_id.split("_")
    dict_pdb_info = {}

    try:
        my_query = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_ID}/{polymer_entity_id}"
        )
        my_query.raise_for_status()
        results = my_query.json()

        dict_pdb_info["pdb_ID"] = pdb_ID
        dict_pdb_info["polymer_id"] = polymer_entity_id
        dict_pdb_info["names"] = get_macromolecular_names_bis(results)
        dict_pdb_info["organism"] = get_organism_names(results)
        return dict_pdb_info

    except requests.exceptions.RequestException as err:
        print(f"Error: {err}")
        return dict_pdb_info


def main(filepath: str):
    """Based on sequence from the input FASTA file, search and print the first 10 PDB id (by default) and their macromolecular names.

    Parameters
    ----------
        filepath: str
            The path to the input FASTA file.
    """
    records = list(SeqIO.parse(filepath, "fasta"))
    sequences = [str(seq.seq) for seq in records]

    list_dict_info_pdb = []
    # dict_IdPDB_seq = {}
    for index, sequence in enumerate(sequences):
        results = API_PDB_search_based_sequence(sequence)
        for pdb_id in results:
            list_dict_info_pdb.append(get_info_one_pdb_id(pdb_id))
        #     set_macromolecular_names = treat_PDB_ID_to_macromolecular_names(results)
        #     dict_IdPDB_seq[index] = {"PDB_ID": results,
        #                             "sequence": sequence,
        #                             "macromolecular_names": set_macromolecular_names}

        for i in list_dict_info_pdb:
            print(i)
    # for index_sequence, values in enumerate(dict_IdPDB_seq.values(), start=1):
    #     PDB_ID, sequence, macromolecular_names = values.values()
    #     logger.success(f"FASTA sequence {index_sequence} (first 20 AA): {sequence[:20]}")
    #     logger.success(f"Putative macromolecular names: {macromolecular_names}")
    #     logger.success(f"Putative PDB IDs (first 10 match): {PDB_ID}")


def is_a_fasta_file(filepath):
    """Check if the given filepath points to an existing fasta file.

    Parameters
    ----------
        filepath : str
            Path of the file.

    Raises
    ------
        argparse.ArgumentTypeError
            If the given filepath is not an existing file,
            or if it does not have a fasta extension.

    Returns
    -------
        str
            The validated path.
    """
    filename = Path(filepath)
    if not Path.is_file(filename):
        raise argparse.ArgumentTypeError(f"{filepath} does not exist")

    if filename.suffix != ".fasta":
        raise argparse.ArgumentTypeError(f"{filepath} is not a .fasta file")
    return filepath


def parse_arg():
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
        prog="search_into_PDB",
        description="Seach PDB IDs and macromolecular names in PDB database, based on an extracted sequence from a fasta file.",
        usage="search_into_PDB.py [-h] --fasta fasta_file",
    )
    parser.add_argument(
        "--fasta",
        type=is_a_fasta_file,
        help="fasta file path",
        required=True,
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arg()
    main(args.fasta)
