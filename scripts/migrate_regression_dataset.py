from pathlib import Path

import grodecoder as gd


TEST_DATA_ROOT_DIR = Path("tests") / "data" / "regression_data"
TEST_DATA_INPUT_FILES_DIR = TEST_DATA_ROOT_DIR / "input_files"
OUTPUT_DIR = Path("expected_results")

assert TEST_DATA_ROOT_DIR.exists()
assert TEST_DATA_INPUT_FILES_DIR.exists()
assert OUTPUT_DIR.exists(), (
    f"Output directory '{OUTPUT_DIR}' does not exist. Please create it before running the script."
)


def all_test_data_files() -> list[Path]:
    """Returns a list of all test data files in the input files directory."""
    return list(TEST_DATA_INPUT_FILES_DIR.glob("*.gro")) + list(TEST_DATA_INPUT_FILES_DIR.glob("*.pdb"))


def main():
    """Main function to print all test data files."""
    test_data_files = all_test_data_files()
    for file in test_data_files:
        output_file = OUTPUT_DIR / file.with_suffix(".json").name
        print(f"Processing {file} -> {output_file}")
        decoded = gd.decode_topology(file)
        with open(output_file, "w") as f:
            f.write(decoded.model_dump_json(indent=2))


if __name__ == "__main__":
    main()
