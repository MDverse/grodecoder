from pathlib import Path

import pytest

import grodecoder as gd
from grodecoder.models import DecodedRead

TEST_DATA_ROOT_DIR = Path(__file__).parent / "data" / "regression_data"
TEST_DATA_INPUT_FILES_DIR = TEST_DATA_ROOT_DIR / "input_files"
TEST_DATA_EXPECTED_RESULTS_DIR = TEST_DATA_ROOT_DIR / "expected_results"


assert TEST_DATA_ROOT_DIR.exists()
assert TEST_DATA_INPUT_FILES_DIR.exists()
assert TEST_DATA_EXPECTED_RESULTS_DIR.exists()

TEST_STRUCTURE_FILES = list(TEST_DATA_INPUT_FILES_DIR.glob("*.gro")) + list(
    TEST_DATA_INPUT_FILES_DIR.glob("*.pdb")
)


def expected_results_file(structure_file: Path) -> Path:
    """Returns the expected results file path for a given structure file."""
    return TEST_DATA_EXPECTED_RESULTS_DIR / structure_file.with_suffix(".json").name


def read_expected_json_file(file_path: Path) -> DecodedRead:
    """Reads a JSON file and returns its content."""
    with open(file_path, "r") as file:
        return DecodedRead.model_validate_json(file.read())


def assert_equal_molecules(result: dict, expected: dict) -> None:
    expected_by_name = {
        molecule.get("name", molecule.get("sequence")): molecule for molecule in expected["inventory"]
    }

    assert result["resolution"] == expected["resolution"]

    # Check that all expected molecules are present in the result.
    # This is order independent.
    for molecule in result["inventory"]:
        molecule_name = molecule.get("name", molecule.get("sequence"))
        assert molecule_name in expected_by_name, f"Missing expected molecule: {molecule_name}"
        expected_molecule = expected_by_name[molecule_name]
        assert molecule == expected_molecule, f"Molecule {molecule_name} does not match expected results."


@pytest.mark.parametrize("structure_file", TEST_STRUCTURE_FILES)
def test_regression(structure_file: Path):
    expected_results_json = expected_results_file(structure_file)

    assert structure_file.exists()
    assert expected_results_json.exists()

    result = gd.decode_structure(structure_file)

    actual = gd.models.DecodedRead.from_decoded(result)
    expected = read_expected_json_file(expected_results_json)

    assert actual == expected, f"Decoded result does not match expected for {structure_file.name}"
