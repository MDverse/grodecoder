"""grodecoder logging utilities."""

import sys
import warnings
from pathlib import Path

from loguru import logger

from .settings import get_settings


def setup_logging(logfile: Path):
    """Sets up logging configuration."""
    debug = get_settings().debug

    fmt = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> <level>{level}: {message}</level>"
    level = "DEBUG" if debug else "INFO"

    logger.remove()

    # Screen logger.
    logger.add(sys.stderr, level=level, format=fmt, colorize=True)

    # File logger
    logger.add(logfile, level=level, format=fmt, colorize=False, mode="w")

    # Sets up loguru to capture warnings (typically MDAnalysis warnings)
    def showwarning(message, *args, **kwargs):
        logger.opt(depth=2).warning(message)

    warnings.showwarning = showwarning  # ty: ignore invalid-assignment


def is_logging_debug() -> bool:
    """Returns True if at least one logging handler is set to level DEBUG."""
    print("COUCOU", get_logging_level())
    return "DEBUG" in get_logging_level()


def get_logging_level() -> list[str]:
    """Returns the list of logging level names (one value per handler)."""
    core_logger = logger._core  # ty: ignore unresolved-attribute
    level_dict = {level.no: level.name for level in core_logger.levels.values()}
    return [level_dict[h.levelno] for h in core_logger.handlers.values()]
