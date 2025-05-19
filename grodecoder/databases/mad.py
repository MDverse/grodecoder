"""MAD database."""

import json
from enum import StrEnum
from pathlib import Path

import bs4
from playwright.sync_api import sync_playwright
from pydantic import BaseModel

MAD_BASE_URL = "https://mad.ibcp.fr"
MAD_EXPLORE_URL = f"{MAD_BASE_URL}/explore"


class ResidueFamily(StrEnum):
    """Types of residues expected in CSML."""

    DETERGENT = "DETERGENT"
    ION = "ION"
    CARBOHYDRATE = "CARBOHYDRATE"
    LIPID = "LIPID"
    POLYMER = "POLYMER"
    PROTEIN = "PROTEIN"
    SMALL_MOLECULE = "SMALL MOLECULE"
    SOLVENT = "SOLVENT"


# Family to actual name as read on MAD's website
FAMILY_TO_MAD = {
    ResidueFamily.DETERGENT: ("Surfactants",),
    ResidueFamily.CARBOHYDRATE: ("Sugars",),
    ResidueFamily.PROTEIN: ("Amino acids",),
    ResidueFamily.LIPID: (
        "Lipids",
        "Lipids",
    ),
    ResidueFamily.SMALL_MOLECULE: (
        "Small molecule",
        "Synthetic nanoparticles",
        "Phytochemicals",
        "Solvents, Small molecule",
    ),
    ResidueFamily.ION: ("Ions",),
    ResidueFamily.SOLVENT: ("Solvents",),
    ResidueFamily.POLYMER: ("Polymers",),
}

MAD_TO_FAMILY = {name: key for key, value in FAMILY_TO_MAD.items() for name in value}


class Residue(BaseModel):
    """Model for a MAD residue."""

    name: str
    alias: str
    family: ResidueFamily
    link: str

    def __hash__(self):
        return hash((self.name, self.alias, self.family, self.link))


class Progress:
    """Tracks the progress of the parsing."""

    def __init__(self):
        self.current_page = 1
        self.last_entry = 0
        self.number_of_entries = 0

    def update(self, page: bs4.Tag):
        """Updates the progress with the current page."""
        footer = page.find("tfoot").find_all("tr")[0]
        text = footer.find("p").text  # Format is "1-20 of 100"
        tokens = text.split()
        self.last_entry = int(tokens[0].split("-")[-1])
        self.number_of_entries = int(tokens[-1])
        self.current_page += 1

    def done(self) -> bool:
        """Returns if the parsing is done, i.e. last entry is the number of entries."""
        return self.number_of_entries > 0 and self.last_entry == self.number_of_entries


class MADParser:
    """Parses the MAD website."""

    def __init__(self):
        self.residues = []

    def _load_page(self, browser, page_number) -> bs4.BeautifulSoup:
        url = f"{MAD_EXPLORE_URL}?page={page_number}"
        page = browser.new_page()
        page.goto(url)
        page.wait_for_selector("tbody")
        return bs4.BeautifulSoup(page.content(), "html.parser")

    def _parse_current_page(self, soup: bs4.BeautifulSoup):
        """Parses the current page and updates the list of residues."""
        tables = soup.find_all("tbody")
        assert len(tables) == 1, f"Expected 1 table, got {len(tables)}"
        rows = tables[0].find_all("tr")
        self.residues += [_parse_residue(row) for row in rows]

    def parse(self):
        """Parses the MAD website and returns a list of residues."""
        progress = Progress()
        with sync_playwright() as pw:
            browser = pw.chromium.launch(headless=True)

            while not progress.done():
                soup = self._load_page(browser, progress.current_page)
                self._parse_current_page(soup)
                progress.update(soup)

            browser.close()

        return self.residues


def _parse_residue(row: bs4.Tag) -> Residue:
    cells = row.find_all("td")

    name = cells[0].text.strip()
    alias = cells[1].text.strip()

    family_text = cells[2].text.strip()
    if family_text not in MAD_TO_FAMILY:
        raise KeyError(f"Unknown residue family: '{family_text}' ({name=}, {alias=})")
    family = MAD_TO_FAMILY[family_text]

    href = cells[0].find("a")["href"]
    link = f"{MAD_BASE_URL}{href}"

    return Residue(
        name=name,
        alias=alias,
        family=family,
        link=link,
    )

def _fix_database(residues: list[Residue]) -> list[Residue]:
    """Fixes the family of the residues known to beeing assigned the wrong one in the first place."""
    assert len(residues) == len(set(residue.name for residue in residues)), "Duplicate residues found"
    return residues


def fetch() -> list[Residue]:
    """Parses the MAD website and returns a list of residues."""
    parser = MADParser()
    residues = parser.parse()
    residues = _fix_database(residues)
    return residues


def read_database(path: Path | str) -> dict[str, Residue]:
    """Reads the CSML database from a JSON file."""
    with open(path, "r") as f:
        data = json.load(f)
        db = [Residue.model_validate(entry) for entry in data]
    return db


if __name__ == "__main__":
    fetch()
