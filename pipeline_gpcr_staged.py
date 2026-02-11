#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
from datetime import datetime


def script_path(name):
    return os.path.join(os.path.dirname(__file__), name)


def run_cmd(cmd, env=None):
    print("\n$", " ".join(cmd))
    subprocess.run(cmd, check=True, env=env)


def default_run_output_root(workspace_root):
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join(workspace_root, "outputs", f"run_{run_id}")


def add_retrieve_args(parser):
    parser.add_argument(
        "--input-dir",
        default=os.path.join(os.path.dirname(__file__), "phmmer_inputs"),
        help="Folder with PHMMER inputs",
    )
    parser.add_argument(
        "--cache-out",
        required=True,
        help="Output JSON cache path (slow step output)",
    )
    parser.add_argument(
        "--email",
        default="a71364@ualg.pt",
        help="Entrez email",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.34,
        help="Delay between Entrez calls",
    )
    parser.add_argument(
        "--progress-step",
        type=int,
        default=10,
        help="Progress reporting step percent",
    )


def add_classify_args(parser, output_root_required=False):
    parser.add_argument("--cache", required=True, help="JSON cache from retrieve step")
    parser.add_argument(
        "--output-root",
        required=output_root_required,
        default=None,
        help="Output run folder for generated excels",
    )


def add_tree_args(parser):
    parser.add_argument("--run-root", required=True, help="Run folder with per-species excels")
    parser.add_argument(
        "--mafft-mode",
        default="auto",
        choices=["quick", "auto", "linsi", "ginsi", "einsi"],
        help="MAFFT alignment strategy",
    )
    parser.add_argument(
        "--prefix",
        default="ALL_CLASSIFIED_AUTOHQ",
        help="Prefix for tree/alignment outputs",
    )
    parser.add_argument(
        "--skip-render",
        action="store_true",
        help="Build tree only (skip PDF rendering)",
    )


def run_retrieve(args):
    cmd = [
        sys.executable,
        script_path("extract_ncbi_cache.py"),
        "--input-dir",
        args.input_dir,
        "--cache-out",
        args.cache_out,
        "--email",
        args.email,
        "--sleep",
        str(args.sleep),
        "--progress-step",
        str(args.progress_step),
    ]
    run_cmd(cmd)


def run_classify(args):
    cmd = [
        sys.executable,
        script_path("build_excels_from_cache.py"),
        "--cache",
        args.cache,
    ]
    if args.output_root:
        cmd.extend(["--output-root", args.output_root])
    run_cmd(cmd)


def run_tree(args):
    cmd = [
        sys.executable,
        script_path("build_phylogeny.py"),
        "--run-root",
        args.run_root,
        "--mafft-mode",
        args.mafft_mode,
        "--prefix",
        args.prefix,
    ]
    run_cmd(cmd)

    if args.skip_render:
        return

    tree_file = os.path.join(args.run_root, "phylogeny_all", f"{args.prefix}.fasttree.nwk")
    family_pdf = os.path.join(
        args.run_root,
        "phylogeny_all",
        f"{args.prefix}.family_colors.pdf",
    )
    subfamily_pdf = os.path.join(
        args.run_root,
        "phylogeny_all",
        f"{args.prefix}.subfamily_colors.pdf",
    )

    render_env = os.environ.copy()
    mpl_dir = os.path.join(args.run_root, ".mplconfig")
    os.makedirs(mpl_dir, exist_ok=True)
    render_env["MPLCONFIGDIR"] = mpl_dir

    run_cmd([
        sys.executable,
        script_path("render_tree_pdf.py"),
        "-i",
        tree_file,
        "-o",
        family_pdf,
        "--color-mode",
        "family",
    ], env=render_env)

    run_cmd([
        sys.executable,
        script_path("render_tree_pdf.py"),
        "-i",
        tree_file,
        "-o",
        subfamily_pdf,
        "--color-mode",
        "subfamily",
    ], env=render_env)


def main():
    parser = argparse.ArgumentParser(
        description="Staged GPCR pipeline: retrieval -> classification/excel -> tree"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_retrieve = sub.add_parser("retrieve", help="Run only NCBI retrieval and save cache")
    add_retrieve_args(p_retrieve)

    p_classify = sub.add_parser("classify", help="Run only classification/excel from cache")
    add_classify_args(p_classify)

    p_tree = sub.add_parser("tree", help="Run only tree + render from run folder")
    add_tree_args(p_tree)

    p_from_cache = sub.add_parser(
        "from-cache",
        help="Run classify + tree starting from existing cache",
    )
    add_classify_args(p_from_cache, output_root_required=True)
    p_from_cache.add_argument(
        "--mafft-mode",
        default="auto",
        choices=["quick", "auto", "linsi", "ginsi", "einsi"],
        help="MAFFT alignment strategy",
    )
    p_from_cache.add_argument(
        "--prefix",
        default="ALL_CLASSIFIED_AUTOHQ",
        help="Prefix for tree/alignment outputs",
    )
    p_from_cache.add_argument(
        "--skip-render",
        action="store_true",
        help="Build tree only (skip PDF rendering)",
    )

    p_full = sub.add_parser("full", help="Run full pipeline from NCBI retrieval")
    add_retrieve_args(p_full)
    p_full.add_argument(
        "--output-root",
        default=None,
        help="Output run folder for excels/tree",
    )
    p_full.add_argument(
        "--mafft-mode",
        default="auto",
        choices=["quick", "auto", "linsi", "ginsi", "einsi"],
        help="MAFFT alignment strategy",
    )
    p_full.add_argument(
        "--prefix",
        default="ALL_CLASSIFIED_AUTOHQ",
        help="Prefix for tree/alignment outputs",
    )
    p_full.add_argument(
        "--skip-render",
        action="store_true",
        help="Build tree only (skip PDF rendering)",
    )

    args = parser.parse_args()

    if args.command == "retrieve":
        run_retrieve(args)
        return

    if args.command == "classify":
        run_classify(args)
        return

    if args.command == "tree":
        run_tree(args)
        return

    if args.command == "from-cache":
        run_classify(args)
        tree_args = argparse.Namespace(
            run_root=args.output_root,
            mafft_mode=args.mafft_mode,
            prefix=args.prefix,
            skip_render=args.skip_render,
        )
        run_tree(tree_args)
        return

    if args.command == "full":
        run_retrieve(args)
        workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        output_root = args.output_root or default_run_output_root(workspace_root)
        classify_args = argparse.Namespace(cache=args.cache_out, output_root=output_root)
        run_classify(classify_args)
        tree_args = argparse.Namespace(
            run_root=output_root,
            mafft_mode=args.mafft_mode,
            prefix=args.prefix,
            skip_render=args.skip_render,
        )
        run_tree(tree_args)
        return


if __name__ == "__main__":
    main()
