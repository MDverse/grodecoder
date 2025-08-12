import sys
from pathlib import Path

import icecream
from loguru import logger

from grodecoder.cli import main

icecream.install()

logger.remove()
logger.add(
    sys.stderr,
    level="INFO",
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{message}</level>",
    colorize=True,
    backtrace=True,
    diagnose=True,
)


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
        main(topology_path)


if __name__ == "__main__":
    test()
