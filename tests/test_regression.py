from pathlib import Path

import json
import pytest

import grodecoder as gd

TEST_DATA_ROOT_DIR = Path(__file__).parent / "data" / "regression_data"
TEST_DATA_INPUT_FILES_DIR = TEST_DATA_ROOT_DIR / "input_files"
TEST_DATA_EXPECTED_RESULTS_DIR = TEST_DATA_ROOT_DIR / "expected_results"

assert TEST_DATA_ROOT_DIR.exists()
assert TEST_DATA_INPUT_FILES_DIR.exists()
assert TEST_DATA_EXPECTED_RESULTS_DIR.exists()


def expected_results_file(topology_file: Path) -> Path:
    """Returns the expected results file path for a given topology file."""
    return TEST_DATA_EXPECTED_RESULTS_DIR / topology_file.with_suffix(".json").name


def all_test_data_files() -> list[Path]:
    """Returns a list of all test data files in the input files directory."""
    foo = list(TEST_DATA_INPUT_FILES_DIR.glob("*.gro")) + list(TEST_DATA_INPUT_FILES_DIR.glob("*.pdb"))
    return foo


def read_json_file(file_path: Path) -> dict:
    """Reads a JSON file and returns its content as a dictionary."""
    with open(file_path, "r") as file:
        return json.load(file)


def compare_molecules(result: dict, expected: dict) -> None:
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


@pytest.mark.parametrize("topology_file", all_test_data_files())
def test_regression(topology_file: Path):
    expected_results_json = expected_results_file(topology_file)

    assert topology_file.exists()
    assert expected_results_json.exists()

    result = gd.decode_topology(topology_file).dump_json()
    expected = read_json_file(expected_results_json)

    compare_molecules(result, expected)


def test_coucou():
    topology_file = TEST_DATA_INPUT_FILES_DIR / "5ZOA.gro"
    expected_results_json = expected_results_file(topology_file)

    result = gd.decode_topology(topology_file).dump_json()
    expected = read_json_file(expected_results_json)

    compare_molecules(result, expected)
