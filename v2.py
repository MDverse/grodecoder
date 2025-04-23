from dataclasses import dataclass
from loguru import logger
from typing import Iterable

import grodecoder as gd
import grodecoder.v1 as v1
import grodecoder.databases as DB

# DEBUG
import icecream
import time

icecream.install()


# Import the "old" grodecoder for comparison.
# To be removed.
v1.logger.remove()


@dataclass
class SmallMoleculeCount:
    description: str
    residue_name: str
    number_of_atoms: int = 0
    number_of_residues: int = 0


class NoDefinitionFoundError(Exception):
    """Base class for all definition not found errors."""

    def __init__(self, residue_name: str, atom_names: Iterable[str]):
        self.residue_name = residue_name
        self.atom_names = list(atom_names)


class NoDefinitionWithMatchingAtoms(NoDefinitionFoundError):
    """Exception raised when no definition is found for a given residue name."""

    def __str__(self):
        return f"No definition found for residue '{self.residue_name}' with atoms {self.atom_names}"


class MultipleDefinitionsWithMatchingAtoms(NoDefinitionFoundError):
    """Exception raised when multiple definitions are found for a given residue name."""

    def __str__(self):
        return f"Multiple definitions found for residue '{self.residue_name}' with atoms {self.atom_names}"


def count_solvents_and_ions(universe) -> SmallMoleculeCount:
    ion_or_solvent = set(DB.RESIDUE_DATABASE.ion_names) | set(DB.RESIDUE_DATABASE.solvent_names)

    counts = {}
    for residue_name in ion_or_solvent:
        selection = universe.select_atoms(f"resname {residue_name}")
        if len(selection) == 0:
            continue

        for residue in selection.residues:
            match = DB.find_residue(residue_name, residue.atoms.names)
            key = str(match)
            counts.setdefault(key, SmallMoleculeCount(match["description"], residue_name))
            counts[key].number_of_atoms += len(residue.atoms)
            counts[key].number_of_residues += 1

    return list(counts.values())



def main():
    topology_path = "grodecoder/data/examples/1BRS.gro"
    universe = gd.read_topology(topology_path)
    c = count_solvents_and_ions(universe)
    ic(c)

    exit()

    # water = universe.select_atoms("resname TIP3")
    # expected_number_of_atoms = len(water)
    # expected_number_of_residues = len(water.residues)
    #
    # start = time.perf_counter()
    # result = count_solvents_and_ions(universe)
    # print(f"method 2: {time.perf_counter() - start:.2f} seconds")
    # ic(result)


if __name__ == "__main__":
    main()
