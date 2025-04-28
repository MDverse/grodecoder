from dataclasses import dataclass
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


def count_solvents_and_ions(universe) -> SmallMoleculeCount:
    """Counts the number of solvent and ion residues and atoms in a universe.

    Solvent and ion residues are identified using both the residue name and the atom names.
    """
    ion_or_solvent = DB.get_ion_names() | DB.get_solvent_names()
    counts = {}
    for residue_name in ion_or_solvent:
        selection = universe.select_atoms(f"resname {residue_name}")
        if len(selection) == 0:
            continue

        for residue in selection.residues:
            match = DB.api.find_ion_or_solvent(residue_name, residue.atoms.names)
            key = str(match)
            counts.setdefault(key, SmallMoleculeCount(match.description, residue_name))
            counts[key].number_of_atoms += len(residue.atoms)
            counts[key].number_of_residues += 1

    return list(counts.values())


def select_protein(universe):
    residue_names = DB.get_amino_acid_names()
    return universe.select_atoms(f"resname {' '.join(residue_names)}")


def main():
    # topology_path = "grodecoder/data/examples/1BRS.gro"
    topology_path = "grodecoder/data/examples/2MAT.pdb"
    universe = gd.read_topology(topology_path)
    
    protein = select_protein(universe)


    ic(protein.atoms.segids)
    ic(protein.atoms.chainIDs)
    exit()

    # If "SYSTEM" is in the segids, segment ids have to be guessed.
    guess_chains = "SYSTEM" in protein.atoms.segids
    if guess_chains:
        print("Guessing chains...")
        protein.guess_chains()
    


    exit()

    start = time.perf_counter()
    c = count_solvents_and_ions(universe)
    ic(c)
    print(f"method 1: {time.perf_counter() - start:.2f} seconds")

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
