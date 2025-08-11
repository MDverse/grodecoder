import sys
import time
from pathlib import Path
from typing import Iterable

import click
from loguru import logger

import grodecoder as gd
import grodecoder.databases as DB
from grodecoder.models import (
    HasAtoms,
    Ion,
    Lipid,
    Nucleotide,
    Protein,
    ProteinSegment,
    SmallMolecule,
    SmallMoleculeBase,
    Solvent,
)
from grodecoder import Json

# DEBUG
import icecream

icecream.install()

logger.remove()
logger.add(
    sys.stderr, level="INFO", format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{message}</level>"
)


def select_protein(universe: gd.UniverseLike) -> Protein:
    """Selects the protein atoms from the universe."""
    protein_residue_names = DB.get_amino_acid_names()
    selection_str = f"resname {' '.join(protein_residue_names)}"

    # Exclude methanol residues from the selection.
    methanol = gd.toputils.find_methanol(universe)
    if methanol:
        selection_str += f" and not index {' '.join(map(str, methanol))}"

    return Protein(universe.select_atoms(selection_str))


def unique_definitions(definitions: Iterable[DB.ResidueDefinition]) -> dict[str, DB.ResidueDefinition]:
    """Returns a dictionary of unique residue definitions by their residue name."""
    unique_items = {}
    for definition in definitions:
        if definition.residue_name not in unique_items:
            unique_items[definition.residue_name] = definition
    return unique_items


def count(
    universe: gd.UniverseLike, definitions: Iterable[DB.ResidueDefinition], model: type = SmallMolecule
) -> list[SmallMolecule]:
    """Counts the residues in the universe based on the definitions.

    Only using the residue name, not the atom names.
    """
    counts = []
    u_definitions = unique_definitions(definitions)
    for residue_name, definition in u_definitions.items():
        selection = universe.select_atoms(f"resname {definition.residue_name}")
        if len(selection) == 0:
            continue
        residue = model(definition, selection)
        counts.append(residue)
    return counts


def count_ions(universe: gd.UniverseLike) -> list[Ion]:
    """Counts the ions in the universe."""
    return count(universe, DB.get_ion_definitions(), model=Ion)


def count_solvents(universe: gd.UniverseLike) -> list[Solvent]:
    """Counts the solvents in the universe."""
    return count(universe, DB.get_solvent_definitions(), model=Solvent)


def count_nucleotides(universe: gd.UniverseLike) -> list[Nucleotide]:
    """Counts the nucleotides in the universe."""
    return count(universe, DB.get_nucleotide_definitions(), model=Nucleotide)


def count_lipids(universe: gd.UniverseLike) -> list[Lipid]:
    """Counts the lipids in the universe."""
    return count(universe, DB.get_lipid_definitions(), model=Lipid)


def count_protein_segments(universe: gd.UniverseLike) -> dict["str", list[ProteinSegment]]:
    """Counts the protein segments in the universe.

    Segments are grouped by their sequence, hence the return type.
    """
    protein = select_protein(universe)

    if len(protein.atoms) == 0:
        return {}

    # Group the protein segments by sequence.
    by_sequence: dict[str, list[ProteinSegment]] = {}
    for segment in protein.iter_segments():
        by_sequence.setdefault(segment.sequence, []).append(segment)

    return by_sequence


def dump_protein_segments_json(by_sequence: dict[str, list[ProteinSegment]]) -> list[dict]:
    """Dumps the protein segments to JSON format.

    Segments are grouped by their sequence.
    """
    segments_json = []
    for sequence, segments in by_sequence.items():
        segments_json.append(
            {
                "sequence": sequence,
                "number_of_atoms": [segment.number_of_atoms for segment in segments],
                "number_of_residues": [segment.number_of_residues for segment in segments],
                "number_of_segments": len(segments),
                # "atom_indices": [segment.atom_indices.tolist() for segment in segments],
                "molecular_type": "protein",
            }
        )
    return segments_json


def remove_identified_atoms(universe: gd.UniverseLike, molecules: list[HasAtoms]) -> gd.UniverseLike:
    """Removes the atoms of the identified molecules from the universe."""
    for molecule in molecules:
        universe -= molecule.atoms
    return universe


def actually_count(topology_path: Path):
    """Identifies and counts the molecules in a topology file."""

    # Along the way, we will remove the counted atoms from the universe to avoid counting them multiple times.
    # This is necessary because some residues can be counted as part of a protein segment and as a solvent or ion.
    # For example, the residue "MET" could be counted as part of a protein segment and as a solvent.
    #
    # This is why the `universe` is created with `select_atoms("all")`: it has to be converted to AtomGroup
    # to allow removing atoms from it later-on.

    universe = gd.read_topology(topology_path).select_atoms("all")

    n_atoms = len(universe.atoms)  # used for logging later-on
    log_prefix = f"{topology_path}"
    logger.debug(f"{log_prefix}: Read {n_atoms:,d} atoms")

    timer_start = time.perf_counter()  # do not include topology reading time in the count

    protein_segments = count_protein_segments(universe)
    for segments in protein_segments.values():
        universe = remove_identified_atoms(universe, segments)
    logger.debug(f"{log_prefix}: Found {len(protein_segments)} protein segments")

    ions = count_ions(universe)
    universe = remove_identified_atoms(universe, ions)
    logger.debug(f"{log_prefix}: Found {len(ions)} ions")

    solvents = count_solvents(universe)
    universe = remove_identified_atoms(universe, solvents)
    logger.debug(f"{log_prefix}: Found {len(solvents)} solvents")

    lipids = count_lipids(universe)
    universe = remove_identified_atoms(universe, lipids)
    logger.debug(f"{log_prefix}: Found {len(lipids)} lipids")

    nucleotides = count_nucleotides(universe)
    universe = remove_identified_atoms(universe, nucleotides)
    logger.debug(f"{log_prefix}: Found {len(nucleotides)} nucleotides")

    # Collecting unknown small molecules that are not defined in the database.
    unknown_molecules = []
    if len(universe.atoms) > 0:
        residue_names = set(universe.residues.resnames)
        for residue_name in residue_names:
            logger.warning(
                f"{log_prefix}: residue {residue_name!r} is unidentified, counting as unknown molecule"
            )
            molecule = SmallMoleculeBase(universe.select_atoms(f"resname {residue_name}"))
            unknown_molecules.append(molecule)

    json_data = {
        "inventory": (
            [ion.dump_json() for ion in ions]
            + [solvent.dump_json() for solvent in solvents]
            + [lipid.dump_json() for lipid in lipids]
            + [nucleotide.dump_json() for nucleotide in nucleotides]
            + [molecule.dump_json() for molecule in unknown_molecules]
            + dump_protein_segments_json(protein_segments)
        )
    }

    elapsed = time.perf_counter() - timer_start
    logger.debug(f"{log_prefix}: {n_atoms:,d} atoms processed in {elapsed:.2f} seconds")
    return json_data


def print_inventory(inventory: dict):
    """Prints the inventory of molecules in a human-readable format."""
    logger.info("Inventory:")
    for molecule in inventory:
        if molecule["molecular_type"] == "protein":
            seq = molecule["sequence"]
            name = f"{seq[0:10]}..."

            if len(name) > 10:
                name += "..."
                name += molecule["sequence"][-min(10, (len(seq) - 10)) :]

            segments = molecule["number_of_segments"]
            atoms = molecule["number_of_atoms"][0]
            residues = molecule["number_of_residues"][0]
            if segments > 1:
                logger.info(f"  {name}: {segments} segments, {residues} residues each, {atoms} atoms each")
            else:
                logger.info(f"  {name}: {segments} segment, {residues} residues, {atoms} atoms")
        else:
            name = molecule["name"]
            atoms = molecule["number_of_atoms"]
            residues = molecule["number_of_residues"]
            logger.info(f"  {name:>5s}: {residues} residues, {atoms} atoms")


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
        top = gd.read_topology(topology_path)

        description: Json = {
            "resolution": gd.toputils.guess_resolution(top),
        }

        description["inventory"] = actually_count(topology_path)["inventory"]
        ic(topology_path, description)


@click.command()
@click.argument("topology_path", type=click.Path(exists=True, path_type=Path))
def main(topology_path):
    logger.info(f"Processing topology file: {topology_path}")
    inventory = actually_count(topology_path)["inventory"]
    print_inventory(inventory)


if __name__ == "__main__":
    # main()
    test()
