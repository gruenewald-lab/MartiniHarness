"""Martini Harness - A framework for dispatching MD simulation campaigns."""

import logging
from pathlib import Path

from .config import SimulationConfig, CampaignConfig
from .simulation import SimulationCampaign

__version__ = "0.1.0"
__author__ = "Fgruenewald"

# Configure package-level logger
logger = logging.getLogger(__name__)


def get_package_data_path(relative_path: str) -> Path:
    """Get the absolute path to a file in the package data directory.

    Args:
        relative_path: Relative path from the package's data directory.

    Returns:
        Absolute Path to the requested file.
    """
    return Path(__file__).parent / "data" / relative_path


def setup_logging(level: int = logging.INFO) -> None:
    """Setup logging for the package.

    Args:
        level: Logging level (e.g., logging.DEBUG, logging.INFO).
    """
    handler = logging.StreamHandler()
    handler.setLevel(level)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)

    logger.addHandler(handler)
    logger.setLevel(level)
