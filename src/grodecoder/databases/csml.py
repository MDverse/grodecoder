"""CHARMM Small Molecule Library (CSML) residues."""

import json
from enum import StrEnum
from pathlib import Path

import bs4
import requests
from pydantic import BaseModel

CHARMM_GUI_BASE_URL = "https://www.charmm-gui.org"
CHARMM_CSML_URL = CHARMM_GUI_BASE_URL + "/?doc=archive&lib=csml"
CHARMM_CSML_VIEW_URL_FORMAT = CHARMM_GUI_BASE_URL + "/?doc=visualization.ngl.archive&pdb_id={}&arg=csml"
CHARMM_CSML_DOWNLOAD_URL_FORMAT = CHARMM_GUI_BASE_URL + "/{}"


class ResidueFamily(StrEnum):
    """Types of residues expected in CSML."""

    PROTEIN = "PROTEIN"
    NUCLEIC_ACID = "NUCLEIC_ACID"
    CARBOHYDRATE = "CARBOHYDRATE"
    LIPID = "MEMBRANE"
    SMALL_MOLECULE = "SMALL MOLECULE"
    ION = "ION"
    SOLVENT = "SOLVENT"
    DETERGENT = "DETERGENT"


class Links(BaseModel):
    """Model for CSML links."""

    view: str
    download: str | None = None


class Residue(BaseModel):
    """Model for a CSML residue.

    Attributes:
        name (str): The residue name expected in topology files.
        description (str): The full name of the residue.
        charge (int): The charge of the residue.
        links (CSMLLinks): Links to the CSML residue visualization and download.
        source (str): The name of the topology file from which the residue is taken.
        family (str): The family of the residue (e.g., protein, nucleic acid, etc.).
    """

    name: str
    description: str
    charge: int
    links: Links
    source: str
    family: ResidueFamily

    def __hash__(self):
        return hash((self.name, self.description, self.charge, self.family))


# Maps topology file names to molecule families.
TOPOLOGY_FILES = {
    "top_all36_prot.rtf": ResidueFamily.PROTEIN,
    "top_all36_na.rtf": ResidueFamily.NUCLEIC_ACID,
    "top_all36_carb.rtf": ResidueFamily.CARBOHYDRATE,
    "top_all36_lipid.rtf": ResidueFamily.LIPID,
    "top_all36_cgenff.rtf": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_moreions.str": ResidueFamily.ION,
    "toppar_all36_nano_lig.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_nano_lig_patch.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_synthetic_polymer.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_synthetic_polymer_patch.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_polymer_solvent.str": ResidueFamily.SOLVENT,
    "toppar_water_ions.str": ResidueFamily.ION,  # !! To be hacked to get the water molecules.
    "toppar_all36_prot_arg0.str": ResidueFamily.PROTEIN,
    "toppar_all36_prot_c36m_d_aminoacids.str": ResidueFamily.PROTEIN,
    "toppar_all36_prot_fluoro_alkanes.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_prot_heme.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_prot_na_combined.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_prot_retinol.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_prot_model.str": ResidueFamily.PROTEIN,
    "toppar_all36_prot_modify_res.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_na_nad_ppi.str": ResidueFamily.NUCLEIC_ACID,
    "toppar_all36_na_rna_modified.str": ResidueFamily.NUCLEIC_ACID,
    "toppar_all36_lipid_sphingo.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_archaeal.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_bacterial.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_cardiolipin.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_cholesterol.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_dag.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_inositol.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_lnp.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_lps.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_mycobacterial.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_miscellaneous.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_model.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_prot.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_tag.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_yeast.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_hmmm.str": ResidueFamily.LIPID,
    "toppar_all36_lipid_detergent.str": ResidueFamily.LIPID,  # !! NOPE: contains detergents as well.
    "toppar_all36_lipid_ether.str": ResidueFamily.LIPID,
    "toppar_all36_carb_glycolipid.str": ResidueFamily.LIPID,
    "toppar_all36_carb_glycopeptide.str": ResidueFamily.LIPID,
    "toppar_all36_carb_imlab.str": ResidueFamily.LIPID,
    "toppar_all36_label_spin.str": ResidueFamily.SMALL_MOLECULE,
    "toppar_all36_label_fluorophore.str": ResidueFamily.SMALL_MOLECULE,
}


def _get_page_source(url: str) -> str:
    """Fetches the HTML source code of a webpage."""
    response = requests.get(url)
    if not response.ok:
        raise Exception(f"Failed to fetch {url}: {response.status_code}")
    return response.text


def _parse_residue(row: bs4.Tag, source: str, family: str) -> Residue:
    """Parses a residue entry from a table row.

    Args:
        row (bs4.Tag): The table row containing the residue information.
        source (str): The name of the topology file from which the residue is taken.
        family (str): The family of the residue (e.g., protein, nucleic acid, etc.).
    """

    def get_download_link() -> str:
        """Returns the download link for the residue (empty string if not found)."""
        cursor = cells[4].find("a")
        if cursor:
            href = cursor["href"]
            return CHARMM_CSML_DOWNLOAD_URL_FORMAT.format(href)
        return ""

    def get_view_link() -> str:
        """Returns the view link for the residue."""
        return CHARMM_CSML_VIEW_URL_FORMAT.format(csmid.lower())

    cells = row.find_all("td")

    # Expects 6 cells, last one is empty.
    assert len(cells) == 6, f"Expected 6 cells, got {len(cells)}"

    csmid = cells[0].text.strip()
    name = cells[1].text.strip()
    charge = int(cells[2].text)
    view_link = get_view_link()
    download_link = get_download_link()

    return Residue(
        name=csmid,
        description=name,
        charge=charge,
        links=Links(view=view_link, download=download_link),
        source=source,
        family=family,
    )


def _parse_topology_file(cursor: bs4.Tag) -> list[Residue]:
    """Parses a topology file entry from a table row.

    Args:
        cursor (bs4.Tag): The "th" tag containing the topology file name.
    """
    topology_file_name = cursor.find("b").text.strip()
    if topology_file_name not in TOPOLOGY_FILES:
        raise ValueError(f"Unknown topology file: {topology_file_name}")
    family = TOPOLOGY_FILES[topology_file_name]
    residues = [
        _parse_residue(row, topology_file_name, family) for row in cursor.find_next("tbody").find_all("tr")
    ]
    return residues


def _fix_database(residues: list[Residue]) -> list[Residue]:
    """Fixes the family of the residues known to beeing assigned the wrong one in the first place."""
    assert len(residues) == len(set(residue.name for residue in residues)), "Duplicate residues found"

    by_id = {residue.name: residue for residue in residues}

    # From toppar_water_ions.str
    by_id["TIP3"].family = ResidueFamily.SOLVENT
    by_id["TP3M"].family = ResidueFamily.SOLVENT

    # From toppar_all36_lipid_detergent.str

    # propanesulfonate
    by_id["CHAPS"].family = ResidueFamily.SMALL_MOLECULE
    by_id["CHAPSO"].family = ResidueFamily.SMALL_MOLECULE

    # cetrimonium bromide
    by_id["CTB10"].family = ResidueFamily.DETERGENT
    by_id["CTB11"].family = ResidueFamily.DETERGENT
    by_id["CTB12"].family = ResidueFamily.DETERGENT
    by_id["CTB13"].family = ResidueFamily.DETERGENT
    by_id["CTB14"].family = ResidueFamily.DETERGENT
    by_id["CTB15"].family = ResidueFamily.DETERGENT
    by_id["CTB16"].family = ResidueFamily.DETERGENT

    #  phosphocholine
    by_id["C6DHPC"].family = ResidueFamily.DETERGENT
    by_id["C7DHPC"].family = ResidueFamily.DETERGENT

    by_id["CYFOS3"].family = ResidueFamily.DETERGENT
    by_id["CYFOS4"].family = ResidueFamily.DETERGENT
    by_id["CYFOS5"].family = ResidueFamily.DETERGENT
    by_id["CYFOS6"].family = ResidueFamily.DETERGENT
    by_id["CYFOS7"].family = ResidueFamily.DETERGENT

    by_id["FOIS11"].family = ResidueFamily.DETERGENT
    by_id["FOIS9"].family = ResidueFamily.DETERGENT

    by_id["FOS10"].family = ResidueFamily.DETERGENT
    by_id["FOS12"].family = ResidueFamily.DETERGENT
    by_id["FOS13"].family = ResidueFamily.DETERGENT
    by_id["FOS14"].family = ResidueFamily.DETERGENT
    by_id["FOS15"].family = ResidueFamily.DETERGENT
    by_id["FOS16"].family = ResidueFamily.DETERGENT

    by_id["UFOS10"].family = ResidueFamily.DETERGENT

    # dimethylamine-n-oxide
    by_id["DDAO"].family = ResidueFamily.DETERGENT
    by_id["DDAOP"].family = ResidueFamily.DETERGENT

    return list(by_id.values())


def _parse_csml_html(source: str) -> list[Residue]:
    """Parses the HTML source code and extracts CSML residues."""
    soup = bs4.BeautifulSoup(source, "html.parser")

    # Residues are listed in a single table with id "topo_table".
    # For each topology file, there is a separate th and tbody.
    table = soup.find("table", {"id": "topo_table"})

    # List of all topology files from which residues are taken.
    topology_files = table.find_next().find_next_siblings("th")

    # Flattens the list of residues.
    residues = [
        residue for topology_file in topology_files for residue in _parse_topology_file(topology_file)
    ]

    return residues


def fetch() -> list[Residue]:
    """Fetches the CSML residues from the CHARMM-GUI website."""
    source = _get_page_source(CHARMM_CSML_URL)
    residues = _parse_csml_html(source)
    residues = _fix_database(residues)
    return residues


def read_database(path: Path | str) -> list[Residue]:
    """Reads the CSML database from a JSON file."""
    with open(path, "r") as f:
        data = json.load(f)
        db = [Residue.model_validate(entry) for entry in data]
    assert len(db) == len(set(residue.name for residue in db)), f"{path}: Duplicate residues found in the database"
    return db
