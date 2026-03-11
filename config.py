#!/usr/bin/env python3
"""
config.py
---------
Central configuration constants for the GPCR B1-family pipeline.

Import from any pipeline module to avoid magic numbers and hardcoded strings::

    from config import ENTREZ_EMAIL, SLEEP_BETWEEN_REQUESTS, GROUP_ORDER
"""

# ---------------------------------------------------------------------------
# NCBI / Entrez
# ---------------------------------------------------------------------------

ENTREZ_EMAIL: str = "a71364@ualg.pt"
"""E-mail passed to the NCBI Entrez API (required by NCBI policy)."""

SLEEP_BETWEEN_REQUESTS: float = 0.34
"""Seconds to pause between consecutive Entrez requests (≤3 req/s without API key)."""

PROGRESS_LOG_STEP: int = 10
"""Print a progress line every *n* percent of a loop."""

# ---------------------------------------------------------------------------
# Pipeline defaults
# ---------------------------------------------------------------------------

DEFAULT_PHMMER_INPUT_DIR: str = "phmmer_inputs"
"""Relative path (from repo root) to the directory with PHMMER result files."""

DEFAULT_MAFFT_MODE: str = "auto"
"""MAFFT alignment strategy used when none is specified on the CLI."""

DEFAULT_TREE_PREFIX: str = "ALL_CLASSIFIED_AUTOHQ"
"""Prefix for phylogeny output files (alignment, tree, PDFs)."""

# ---------------------------------------------------------------------------
# Taxonomy / classification
# ---------------------------------------------------------------------------

GROUP_ORDER: list[str] = ["CALCR", "CRHR", "GCGR", "SCTR", "PTHR", "PDFR"]
"""Display order of receptor super-groups in Excel workbooks."""

CLASSIFICATION_ORDER: list[str] = [
    "VIPR2",
    "VIPR1",
    "GLP2R",
    "GLP1R",
    "PTH2R",
    "PTH1R",
    "CRHR2",
    "CRHR1",
    "SCTR",
    "GCGR",
    "GIPR",
    "GHRHR",
    "ADCYAP1R1",
    "CALCR",
    "CALCRL",
    "PDFR",
]
"""Priority order for subfamily classification (first match wins)."""

GROUP_SUBFAMILIES: dict[str, list[str]] = {
    "CALCR": ["CALCRL", "CALCR"],
    "CRHR": ["CRHR1", "CRHR2"],
    "GCGR": ["GCGR", "GLP1R", "GLP2R", "GIPR"],
    "SCTR": ["GHRHR", "ADCYAP1R1", "VIPR1", "VIPR2", "SCTR"],
    "PTHR": ["PTH1R", "PTH2R"],
    "PDFR": ["PDFR"],
}
"""Mapping from super-group name to its constituent subfamily names."""

SUBFAMILY_TO_GROUP: dict[str, str] = {
    subfamily: group
    for group, subfamilies in GROUP_SUBFAMILIES.items()
    for subfamily in subfamilies
}
"""Reverse mapping: subfamily → super-group (auto-generated from GROUP_SUBFAMILIES)."""

# ---------------------------------------------------------------------------
# Known NCBI misannotations
# ---------------------------------------------------------------------------

ACCESSION_SUBFAMILY_OVERRIDES: dict[str, str] = {
    "XP_052822567": "PDFR",
    "XP_052822567.1": "PDFR",
}
"""Accessions whose subfamily classification is forced regardless of description."""

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

PIPELINE_VERSION: str = "1.1.0"
"""Semantic version string embedded in generated outputs and the --version flag."""

OUTPUTS_ROOT_NAME: str = "outputs"
"""Top-level directory (relative to repo root) for all generated files."""
