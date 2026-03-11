#!/usr/bin/env python3
"""
logging_setup.py
----------------
Centralised logging configuration for the GPCR B1-family pipeline.

All pipeline scripts can call ``setup_logger()`` to obtain a consistently
formatted logger instead of using bare ``print()`` calls.

Usage::

    from logging_setup import setup_logger

    log = setup_logger("extract_ncbi_cache")
    log.info("Starting extraction …")
    log.warning("Accession %s not found", acc)
"""

import logging
import os
import sys
from datetime import datetime
from typing import Optional


_LOG_FORMAT = "%(asctime)s [%(levelname)-8s] %(name)s: %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logger(
    name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    file_level: Optional[int] = None,
) -> logging.Logger:
    """Create and return a logger with console (and optionally file) handlers.

    Parameters
    ----------
    name:
        Logger name — usually the calling module's ``__name__`` or script stem.
    level:
        Log level for the console handler (default: INFO).
    log_file:
        If given, also write logs to this file path.
    file_level:
        Log level for the file handler (default: same as *level*).

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(name)
    if logger.handlers:
        # Already configured — return as-is to avoid duplicate handlers.
        return logger

    logger.setLevel(min(level, file_level or level))
    formatter = logging.Formatter(_LOG_FORMAT, datefmt=_DATE_FORMAT)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if log_file:
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(file_level or level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_run_log_path(outputs_root: str, script_name: str) -> str:
    """Return a timestamped log file path inside *outputs_root*/logs/.

    Example::

        get_run_log_path("/data/outputs", "extract_ncbi_cache")
        # → "/data/outputs/logs/extract_ncbi_cache_20260311_143022.log"
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = os.path.join(outputs_root, "logs")
    return os.path.join(log_dir, f"{script_name}_{timestamp}.log")
