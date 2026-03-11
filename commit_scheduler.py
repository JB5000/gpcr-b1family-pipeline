#!/usr/bin/env python3
"""
commit_scheduler.py
-------------------
Generates a commit development plan for a project:
- Scans the repository structure
- Suggests which modules/folders to improve
- Proposes a list of meaningful commit messages with priorities

Usage:
    python commit_scheduler.py [--n-commits N] [--seed SEED]
"""

import argparse
import os
import random
import datetime


COMMIT_TEMPLATES = [
    ("feat", "Add {module} utility module for shared helpers"),
    ("feat", "Add {module} configuration constants"),
    ("feat", "Add input validation for {module}"),
    ("feat", "Add unit tests for {module}"),
    ("feat", "Add progress logging to {module}"),
    ("feat", "Add --version flag and CLI improvements to {module}"),
    ("feat", "Add summary statistics export to {module}"),
    ("refactor", "Add type hints to {module}"),
    ("refactor", "Refactor {module} for better readability"),
    ("refactor", "Extract helper functions from {module}"),
    ("fix", "Fix edge case in {module} for empty inputs"),
    ("fix", "Fix path handling in {module} for cross-platform support"),
    ("docs", "Update README with usage examples and badges"),
    ("docs", "Add docstrings to all public functions in {module}"),
    ("docs", "Add CHANGELOG.md with release notes"),
    ("chore", "Update requirements.txt with dev dependencies"),
    ("chore", "Add .gitignore entries for outputs and caches"),
    ("chore", "Add pre-commit hooks configuration"),
    ("test", "Add integration test for full pipeline run"),
    ("test", "Add test fixtures and conftest.py"),
]

FOLDERS = {
    "src": ["pipeline_gpcr_staged.py", "master_pipeline_gpcr.py",
            "extract_ncbi_cache.py", "generate_outputs.py",
            "build_excels_from_cache.py", "build_phylogeny.py",
            "render_tree_pdf.py", "export_species_fastas_4letter.py"],
    "new_modules": ["utils", "config", "validators", "logging_setup",
                    "summary_stats", "ncbi_utils", "phylo_utils"],
    "tests": ["test_pipeline", "test_validators", "test_utils",
              "test_extract_cache", "test_build_excels"],
}


def scan_repo(repo_root: str) -> list[str]:
    """Return a list of Python files in the repo root."""
    modules = []
    for f in os.listdir(repo_root):
        if f.endswith(".py") and not f.startswith("_"):
            modules.append(os.path.splitext(f)[0])
    return modules


def generate_plan(n_commits: int, seed: int, repo_root: str) -> list[dict]:
    random.seed(seed)
    modules = scan_repo(repo_root)
    all_modules = modules + FOLDERS["new_modules"] + FOLDERS["tests"]

    plan = []
    used = set()
    attempts = 0
    while len(plan) < n_commits and attempts < 200:
        attempts += 1
        commit_type, template = random.choice(COMMIT_TEMPLATES)
        module = random.choice(all_modules)
        message = template.format(module=module)
        key = (commit_type, message)
        if key in used:
            continue
        used.add(key)
        plan.append({
            "index": len(plan) + 1,
            "type": commit_type,
            "message": message,
            "target": module,
        })
    return plan


def print_plan(plan: list[dict], n_commits: int, today: str) -> None:
    print(f"\n{'='*60}")
    print(f"  COMMIT DEVELOPMENT PLAN  —  {today}")
    print(f"  Total commits planned: {n_commits}")
    print(f"{'='*60}\n")
    for item in plan:
        print(f"  [{item['index']:02d}] {item['type'].upper():8s} | {item['message']}")
        print(f"        Target module: {item['target']}\n")
    print("="*60)
    print(f"\nRun these commits in order. Each should correspond to a")
    print(f"real code change in the target module.")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(description="Generate a commit development plan.")
    parser.add_argument("--n-commits", type=int, default=9,
                        help="Number of commits to plan (default: 9)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility (default: 42)")
    parser.add_argument("--repo-root", type=str,
                        default=os.path.dirname(os.path.abspath(__file__)),
                        help="Path to the repository root")
    args = parser.parse_args()

    today = datetime.date.today().isoformat()
    plan = generate_plan(args.n_commits, args.seed, args.repo_root)
    print_plan(plan, args.n_commits, today)


if __name__ == "__main__":
    main()
