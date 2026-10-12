"""Render three scientific SVG/PDF figures from actual compact case outputs."""

import os

for key in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[key] = "1"

import argparse
import csv
import json
import math
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
matplotlib.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "svg.fonttype": "none",
        "svg.hashsalt": "mechanical-lie-v0.1.0",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    }
)
BLUE, ORANGE, GRAY = "#236b8e", "#cc7431", "#5a6472"


def save(figure, output, name, preview=None):
    for suffix in ("svg", "pdf"):
        figure.savefig(
            output / f"{name}.{suffix}",
            metadata={
                "Creator": "mechanical-lie plot_cases.py",
                **({"Date": None} if suffix == "svg" else {"CreationDate": None}),
            },
        )
    if preview is not None:
        preview.mkdir(parents=True, exist_ok=True)
        figure.savefig(preview / f"{name}.png", dpi=160)
    plt.close(figure)


def main(data, output, preview=None):
    output.mkdir(parents=True, exist_ok=True)
    examples = json.loads((data / "tree_examples.json").read_text())["examples"]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4))
    for ax, item in zip(axes, examples):
        tree = item["trees"][0]["tree"]
        nodes = []
        edges = []

        def walk(node, parent=None):
            i = len(nodes)
            nodes.append(node[0])
            if parent is not None:
                edges.append((parent, i))
            for child in node[1]:
                walk(child, i)

        walk(tree)
        positions = [(i - (len(nodes) - 1) / 2, 0) for i in range(len(nodes))]
        for i, j in edges:
            ax.plot([positions[i][0], positions[j][0]], [0, 0], color=BLUE, lw=2.4)
        for i, count in enumerate(nodes):
            x, y = positions[i]
            ax.scatter([x], [y], s=280, color=BLUE, zorder=3)
            for k in range(count):
                angle = math.pi / 2 + (k - (count - 1) / 2) * 0.6
                xp, yp = x + 0.85 * math.cos(angle), y + 0.85 * math.sin(angle)
                ax.plot([x, xp], [y, yp], ls="--", color=ORANGE, lw=1.8)
                ax.scatter([xp], [yp], s=95, facecolors="white", edgecolors=ORANGE, zorder=3)
        ax.set_title(item["name"], fontsize=14, pad=13)
        ax.text(0, -0.8, item["polynomial"], ha="center", va="top", fontsize=11, color=GRAY)
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.45, 1.35)
        ax.axis("off")
    fig.subplots_adjust(top=0.72, bottom=0.16, wspace=0.28)
    fig.suptitle("Tree contractions for V(q) = q²/2 + q⁴/4", fontsize=15, y=0.99)
    fig.text(
        0.5,
        0.015,
        "● potential derivative vertex      solid: contracted derivative edge      dashed: momentum leaf",
        ha="center",
        fontsize=9,
        color=GRAY,
    )
    save(fig, output, "bracket_trees", preview)
    rows = list(csv.DictReader((ROOT / "numerical_case/summaries/graded_dimensions.csv").open()))
    n = [int(r["order"]) for r in rows]
    dim = [int(r["dimension"]) for r in rows]
    positive = [int(r["positive_dimension"]) for r in rows]
    pot = [d - g for d, g in zip(dim, positive)]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.7))
    axes[0].plot(n, dim, "o-", color=BLUE, label="Total source dimension")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Dimension (log scale)")
    axes[0].set_title("Exact graded dimension growth")
    axes[0].legend(frameon=False)
    axes[1].plot(n, positive, "o-", color=BLUE, label="Positive momentum degree")
    odd = [(order, count) for order, count in zip(n, pot) if count]
    axes[1].plot(
        [x[0] for x in odd],
        [x[1] for x in odd],
        "s-",
        color=ORANGE,
        label="Potential part (odd orders only)",
    )
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Component dimension (log scale)")
    axes[1].set_title("Positive / potential decomposition")
    axes[1].legend(frameon=False, fontsize=9, loc="upper left")
    for ax in axes:
        ax.set_xlabel("Original-letter order")
        ax.set_xticks([1, 4, 8, 12, 16, 18])
        ax.grid(axis="y", alpha=0.18)
    fig.tight_layout()
    save(fig, output, "graded_dimensions", preview)
    pattern = json.loads((data / "matrix_pattern.json").read_text())
    stats = list(csv.DictReader((data / "sparsity.csv").open()))
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.8))
    points = pattern["nonzero_positions"]
    axes[0].scatter([x[1] for x in points], [x[0] for x in points], s=11, color=BLUE, marker="s")
    axes[0].set_xlim(-0.7, pattern["columns"] - 0.3)
    axes[0].set_ylim(pattern["rows"] - 0.3, -0.7)
    axes[0].set_xlabel("Source basis column")
    axes[0].set_ylabel("Canonical tree row")
    density = len(points) / (pattern["rows"] * pattern["columns"])
    axes[0].set_title(
        f"Order {pattern['order']}, momentum degree {pattern['momentum_degree']}\n{pattern['rows']} × {pattern['columns']}; {len(points)} nonzeros ({density:.1%})",
        fontsize=11,
    )
    orders = sorted({int(r["order"]) for r in stats})
    support = [
        max(int(r["largest_column_support"]) for r in stats if int(r["order"]) == k) for k in orders
    ]
    axes[1].plot(orders, support, "o-", color=ORANGE)
    axes[1].set_xlabel("Original-letter order")
    axes[1].set_ylabel("Nonzero tree coefficients")
    axes[1].set_title("Largest actual column support")
    axes[1].set_xticks(orders)
    axes[1].grid(axis="y", alpha=0.18)
    fig.tight_layout()
    save(fig, output, "sparsity", preview)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--preview-dir", type=Path, help="Optional PNG previews for visual inspection"
    )
    args = parser.parse_args()
    main(args.data, args.output, args.preview_dir)
