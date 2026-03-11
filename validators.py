#!/usr/bin/env python3
"""
validators.py
-------------
Input validation helpers for the GPCR B1-family pipeline.

These functions raise ``ValueError`` (or ``SystemExit`` for fatal CLI errors)
when inputs do not meet the expected constraints.  Import them at the top of
any pipeline script to fail fast with clear error messages::

    from validators import validate_email, validate_input_dir, validate_cache_json

    validate_email(args.email)
    validate_input_dir(args.input_dir)
    validate_cache_json(args.cache)
"""

import json
import os
import re
from typing import Optional


# ---------------------------------------------------------------------------
# Entrez / NCBI
# ---------------------------------------------------------------------------

_EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


def validate_email(email: str) -> str:
    """Raise ``ValueError`` if *email* does not look like a valid e-mail address.

    Returns *email* unchanged so it can be used inline::

        Entrez.email = validate_email(args.email)
    """
    if not _EMAIL_RE.match(email):
        raise ValueError(
            f"Invalid Entrez e-mail address: {email!r}\n"
            "NCBI requires a valid e-mail for Entrez requests."
        )
    return email


def validate_sleep(sleep: float) -> float:
    """Raise ``ValueError`` if *sleep* is negative or too short for NCBI limits.

    NCBI allows at most 3 requests per second without an API key, so the
    minimum safe sleep is ~0.34 s.
    """
    if sleep < 0:
        raise ValueError(f"Sleep interval must be non-negative, got {sleep}")
    if sleep < 0.1:
        raise ValueError(
            f"Sleep interval {sleep} s is too short — NCBI enforces ≤3 req/s. "
            "Use at least 0.34 s."
        )
    return sleep


# ---------------------------------------------------------------------------
# Filesystem
# ---------------------------------------------------------------------------

def validate_input_dir(path: Optional[str], default_name: str = "phmmer_inputs") -> str:
    """Return the resolved input directory, raising ``SystemExit`` if missing.

    If *path* is ``None``, falls back to *default_name* next to this file.
    """
    resolved = path or os.path.join(os.path.dirname(__file__), default_name)
    if not os.path.isdir(resolved):
        raise SystemExit(
            f"Input directory not found: {resolved}\n"
            "Use --input-dir to specify the folder containing PHMMER result files."
        )
    return resolved


def validate_output_dir(path: str, create: bool = True) -> str:
    """Return *path*, optionally creating it.  Raise ``ValueError`` if it is a file."""
    if os.path.isfile(path):
        raise ValueError(f"Output path exists as a file, not a directory: {path}")
    if create:
        os.makedirs(path, exist_ok=True)
    return path


def validate_file_exists(path: str, label: str = "File") -> str:
    """Raise ``SystemExit`` if *path* does not exist or is not a regular file."""
    if not os.path.isfile(path):
        raise SystemExit(f"{label} not found: {path}")
    return path


# ---------------------------------------------------------------------------
# JSON cache
# ---------------------------------------------------------------------------

_REQUIRED_CACHE_KEYS = {"generated_at", "file_records", "accession_data", "gene_summary"}


def validate_cache_json(path: str) -> dict:
    """Load and validate the extraction cache JSON.

    Raises
    ------
    SystemExit
        If the file is missing, is not valid JSON, or is missing required keys.

    Returns
    -------
    dict
        The parsed cache dictionary.
    """
    validate_file_exists(path, label="Cache JSON")
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Cache file is not valid JSON: {path}\n{exc}") from exc

    missing = _REQUIRED_CACHE_KEYS - data.keys()
    if missing:
        raise SystemExit(
            f"Cache JSON is missing required keys: {sorted(missing)}\n"
            f"File: {path}"
        )
    return data


# ---------------------------------------------------------------------------
# MAFFT / tree
# ---------------------------------------------------------------------------

_VALID_MAFFT_MODES = {"quick", "auto", "linsi", "ginsi", "einsi"}


def validate_mafft_mode(mode: str) -> str:
    """Raise ``ValueError`` if *mode* is not a recognised MAFFT preset."""
    if mode not in _VALID_MAFFT_MODES:
        raise ValueError(
            f"Unknown MAFFT mode: {mode!r}. "
            f"Choose one of: {sorted(_VALID_MAFFT_MODES)}"
        )
    return mode
