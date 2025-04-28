"""CHARMM Small Molecule Library (CSML) residues."""

import bs4
import requests
from pydantic import BaseModel


CHARMM_GUI_BASE_URL = "https://www.charmm-gui.org"
CHARMM_CSML_URL = CHARMM_GUI_BASE_URL + "/?doc=archive&lib=csml"
CHARMM_CSML_VIEW_URL_FORMAT = CHARMM_GUI_BASE_URL + "/?doc=visualization.ngl.archive&pdb_id={}&arg=csml"
CHARMM_CSML_DOWNLOAD_URL_FORMAT = CHARMM_GUI_BASE_URL + "/{}"


class CSMLLinks(BaseModel):
    """Model for CSML links."""

    view: str
    download: str | None = None


class CSMLResidue(BaseModel):
    """Model for a CSML residue.

    Attributes
        source: str
            The name of the topology file from which the residue is taken.
    """

    csml_id: str
    description: str
    charge: int
    links: CSMLLinks
    source: str


def _get_page_source(url: str) -> str:
    """Fetches the HTML source code of a webpage."""
    response = requests.get(url)
    return response.text


def _parse_residue(row: bs4.Tag, source: str) -> CSMLResidue:
    """Parses a residue entry from a table row."""
    cells = row.find_all("td")

    # Expects 6 cells, last one is empty.
    assert len(cells) == 6, f"Expected 6 cells, got {len(cells)}"

    csmid = cells[0].text.strip()
    name = cells[1].text.strip()
    charge = int(cells[2].text)
    view_link = CHARMM_CSML_VIEW_URL_FORMAT.format(csmid.lower())

    # Download link may not exist.
    download_link = ""
    if cells[4].find("a"):
        download_link = CHARMM_CSML_DOWNLOAD_URL_FORMAT.format(cells[4].find("a")["href"])

    return CSMLResidue(
        csml_id=csmid,
        description=name,
        charge=charge,
        links=CSMLLinks(view=view_link, download=download_link),
        source=source,
    )


def _parse_topology_file(cursor: bs4.Tag) -> list[dict]:
    """Parses a topology file entry from a table row.

    Args:
        cursor (bs4.Tag): The "th" tag containing the topology file name.
    """
    topology_file_name = cursor.find("b").text.strip()
    residues = [
        _parse_residue(row, topology_file_name).model_dump()
        for row in cursor.find_next("tbody").find_all("tr")
    ]
    return residues


def fetch() -> list[CSMLResidue]:
    """Fetches the CSML residues from the CHARMM-GUI website."""
    source = _get_page_source(CHARMM_CSML_URL)
    soup = bs4.BeautifulSoup(source, "html.parser")

    # Residues are listed in a single table with id "topo_table".
    # For each topology file, there is a separate th and tbody.
    table = soup.find("table", {"id": "topo_table"})

    topology_files = table.find_next().find_next_siblings("th")

    return [residue for topology_file in topology_files for residue in _parse_topology_file(topology_file)]
