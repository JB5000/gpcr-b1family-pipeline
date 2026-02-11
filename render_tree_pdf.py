#!/usr/bin/env python3
import argparse
import os
import re
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from Bio import Phylo

FAMILY_COLORS = {
    "CALCR": "#d97706",
    "CRHR": "#7c3aed",
    "GCGR": "#059669",
    "SCTR": "#dc2626",
    "PTHR": "#2563eb",
}

SUBFAMILY_COLORS = {
    "CALCRL": "#f59e0b",
    "CALCR": "#b45309",
    "CRHR1": "#8b5cf6",
    "CRHR2": "#6d28d9",
    "GCGR": "#10b981",
    "GLP1R": "#34d399",
    "GLP2R": "#059669",
    "GIPR": "#047857",
    "GHRHR": "#ef4444",
    "ADCYAP1R1": "#f97316",
    "VIPR1": "#dc2626",
    "VIPR2": "#b91c1c",
    "SCTR": "#991b1b",
    "PTH1R": "#3b82f6",
    "PTH2R": "#1d4ed8",
}

LABEL_RE = re.compile(r"^([^_]+)_(.+)_(CALCR|CRHR|GCGR|SCTR|PTHR)_([A-Za-z0-9]+)$")


def parse_args():
    parser = argparse.ArgumentParser(description="Render a Newick tree to PDF")
    parser.add_argument("-i", "--input", required=True, help="Input Newick tree file")
    parser.add_argument("-o", "--output", required=True, help="Output PDF file path")
    parser.add_argument(
        "--no-midpoint-root",
        action="store_true",
        help="Disable midpoint rooting before rendering",
    )
    parser.add_argument(
        "--no-ladderize",
        action="store_true",
        help="Disable ladderize ordering before rendering",
    )
    parser.add_argument(
        "--hide-legend",
        action="store_true",
        help="Hide family/subfamily legend block",
    )
    parser.add_argument(
        "--color-mode",
        choices=["family", "subfamily"],
        default="subfamily",
        help="Color tips by major family or by subfamily",
    )
    return parser.parse_args()


def parse_label_fields(label):
    match = LABEL_RE.match(str(label or ""))
    if not match:
        return None
    species, accession, family, subfamily = match.groups()
    return species, accession, family, subfamily


def pretty_tip_label(clade):
    # Hide internal-node default labels like "Clade".
    name = getattr(clade, "name", None)
    if not name:
        return None
    fields = parse_label_fields(name)
    if not fields:
        return None
    species, accession, family, subfamily = fields
    return f"{species} {accession} | {family} > {subfamily}"


def add_legend(ax, color_mode):
    if color_mode == "family":
        lines = [
            "Color mode: family",
            "CALCR",
            "CRHR",
            "GCGR",
            "SCTR",
            "PTHR",
        ]
    else:
        lines = [
            "Color mode: subfamily",
            "CALCR: CALCRL, CALCR",
            "CRHR: CRHR1, CRHR2",
            "GCGR: GCGR, GLP1R, GLP2R, GIPR",
            "SCTR: GHRHR, ADCYAP1R1, VIPR1, VIPR2, SCTR",
            "PTHR: PTH1R, PTH2R",
        ]
    ax.text(
        1.01,
        1.0,
        "\n".join(lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7,
        family="monospace",
        bbox={
            "facecolor": "#ffffff",
            "edgecolor": "#334155",
            "linewidth": 0.6,
            "pad": 0.4,
        },
    )


def main():
    args = parse_args()
    tree = Phylo.read(args.input, "newick")
    if not args.no_midpoint_root:
        try:
            tree.root_at_midpoint()
        except Exception:
            pass
    if not args.no_ladderize:
        tree.ladderize()
    fig = plt.figure(figsize=(19, 50))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax, label_func=pretty_tip_label)

    for text in ax.texts:
        label = text.get_text()
        if " | " not in label or " > " not in label:
            continue
        meta = label.split(" | ", 1)[1]
        family, subfamily = [p.strip() for p in meta.split(">", 1)]
        if args.color_mode == "family":
            edge = FAMILY_COLORS.get(family, "#111827")
        else:
            edge = SUBFAMILY_COLORS.get(subfamily, FAMILY_COLORS.get(family, "#111827"))
        text.set_color(edge)
        text.set_fontsize(6)
        text.set_bbox({
            "facecolor": "#f8fafc",
            "edgecolor": edge,
            "linewidth": 0.4,
            "pad": 0.2,
        })

    if not args.hide_legend:
        add_legend(ax, args.color_mode)

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 0.88, 1))
    fig.savefig(args.output, format="pdf")


if __name__ == "__main__":
    main()
