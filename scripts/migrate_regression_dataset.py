from pathlib import Path
from loguru import logger

import grodecoder as gd


TEST_DATA_ROOT_DIR = Path("tests") / "data" / "regression_data"
TEST_DATA_INPUT_FILES_DIR = TEST_DATA_ROOT_DIR / "input_files"
OUTPUT_DIR = TEST_DATA_ROOT_DIR / "expected_results"

assert TEST_DATA_ROOT_DIR.exists()
assert TEST_DATA_INPUT_FILES_DIR.exists()


def backup_output_dir():
    """Backs up the existing output directory by renaming it."""
    def new_backup_dir(index) -> Path:
        return OUTPUT_DIR.parent / f"{OUTPUT_DIR.name}#{index}"

    index = 1
    backup_dir = new_backup_dir(index)
    while backup_dir.exists():
        index += 1
        backup_dir = new_backup_dir(index)

    logger.info(f"Backing up expected results to {backup_dir}")
    OUTPUT_DIR.rename(backup_dir)


if OUTPUT_DIR.exists():
    backup_output_dir()
OUTPUT_DIR.mkdir(parents=True, exist_ok=False)


def all_test_data_files() -> list[Path]:
    """Returns a list of all test data files in the input files directory."""
    return list(TEST_DATA_INPUT_FILES_DIR.glob("*.gro")) + list(TEST_DATA_INPUT_FILES_DIR.glob("*.pdb"))


def main():
    """Main function to print all test data files."""
    test_data_files = all_test_data_files()
    for file in test_data_files:
        output_file = OUTPUT_DIR / file.with_suffix(".json").name
        logger.info(f"Processing {file} -> {output_file}")
        decoded = gd.decode_topology(file)
        with open(output_file, "w") as f:
            f.write(decoded.model_dump_json(indent=2))
    logger.info("Migration completed.")


if __name__ == "__main__":
    main()
