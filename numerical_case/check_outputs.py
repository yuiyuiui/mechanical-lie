"""Compare exact finite outputs with published summaries and check figure files."""

import argparse
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def main(data: Path, figures: Path):
    for name in ["checks.json", "sparsity.csv", "matrix_pattern.json", "tree_examples.json"]:
        actual = (data / name).read_bytes()
        expected = (ROOT / "numerical_case/summaries" / name).read_bytes()
        if actual != expected:
            raise RuntimeError(f"Computed summary differs: {name}")
    for name in ["bracket_trees", "graded_dimensions", "sparsity"]:
        svg = figures / (name + ".svg")
        root = ET.fromstring(svg.read_bytes())
        if not root.tag.endswith("svg") or not (figures / (name + ".pdf")).read_bytes().startswith(
            b"%PDF"
        ):
            raise RuntimeError(f"Invalid figure pair: {name}")
        print(f"{svg.name}: {hashlib.sha256(svg.read_bytes()).hexdigest()}")
    print("All four exact compact summaries match; all three SVG/PDF pairs parse.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()
    main(args.data, args.figures)
