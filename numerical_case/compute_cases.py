"""Compute fixed exact cases; write new outputs and compact reproducible summaries."""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import shutil
import subprocess

from mechanical_lie import Polynomial, bracket, kinetic, phase_space, potential
from mechanical_lie._engine import core as e
from mechanical_lie.validation import rebuild_basis, verify_certificate_block

ROOT = Path(__file__).resolve().parents[1]


def polynomial_record(p):
    return {
        "dimension": p.dimension,
        "terms": [
            {"exponents": list(a), "coefficient": str(c)} for a, c in sorted(p.terms.items())
        ],
    }


def full_coefficient_checks():
    results = []
    for d in (1, 2):
        q, _ = phase_space(d)
        v = q[0] ** 2 / 2 + q[0] ** 4 / 4 if d == 1 else q[0] ** 2 / 2 + q[0] * q[1] + q[1] ** 4 / 4
        t = Polynomial.kinetic(d)
        from mechanical_lie import Hamiltonian, poisson

        def px(x):
            return v if not x else poisson(px(x[1]), poisson(px(x[0]), t))

        def pw(w):
            if len(w) == 1:
                return t if w[0] == e.A else poisson(px(w[0][1]), t)
            a, b = e.standard_split(w)
            return poisson(pw(a), pw(b))

        alphabet = tuple(sorted([e.A] + [e.Y(x) for n in range(1, 4) for x in e.magmatic(n)]))
        comparisons = 0
        for order in range(1, 7):
            pairs = [
                (Hamiltonian(e.theta_lyndon(w)), pw(w))
                for w in e.words_of_weight(order, alphabet)
                if e.lyndon(w)
            ]
            if order % 2:
                pairs += [(Hamiltonian(e.theta_X(x)), px(x)) for x in e.magmatic((order + 1) // 2)]
            for expression, reference in pairs:
                if expression.realize(v) != reference:
                    raise RuntimeError(f"Full polynomial mismatch at dimension {d}, order {order}")
                comparisons += 1
        results.append(
            {
                "dimension": d,
                "max_order": 6,
                "comparison_kind": "complete coefficient equality",
                "comparisons": comparisons,
                "potential": polynomial_record(v),
            }
        )
    return results


def compute(output: Path):
    output.mkdir(parents=True, exist_ok=False)
    with contextlib.redirect_stdout(io.StringIO()):
        summary = rebuild_basis(10, output / "matrices")
    with gzip.open(output / "matrices/rank_certificates.json.gz", "rt") as stream:
        blocks = json.load(stream)
    sparsity = []
    for block in blocks:
        rank = verify_certificate_block(block)
        rows = sorted({row for column in block["columns"] for row in column})
        nnz = sum(len(column) for column in block["columns"])
        sparsity.append(
            {
                "order": block["order"],
                "momentum_degree": block["degree"],
                "rows": len(rows),
                "columns": block["dimension"],
                "rank": rank,
                "nonzeros": nnz,
                "dense_entries": len(rows) * block["dimension"],
                "largest_column_support": max(map(len, block["columns"])),
            }
        )
    with (output / "sparsity.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(sparsity[0]))
        writer.writeheader()
        writer.writerows(sparsity)
    largest = max(blocks, key=lambda b: len({r for c in b["columns"] for r in c}) * b["dimension"])
    rows = sorted({row for column in largest["columns"] for row in column})
    positions = {row: i for i, row in enumerate(rows)}
    pattern = {
        "order": largest["order"],
        "momentum_degree": largest["degree"],
        "rows": len(rows),
        "columns": largest["dimension"],
        "nonzero_positions": [
            [positions[row], j] for j, c in enumerate(largest["columns"]) for row in c
        ],
    }
    (output / "matrix_pattern.json").write_text(json.dumps(pattern, indent=2) + "\n")
    q, _ = phase_space(1)
    v = q[0] ** 2 / 2 + q[0] ** 4 / 4
    a, b = kinetic(), potential()
    examples = []
    for name, expression in [
        ("V", b),
        ("{V,T}", bracket(b, a)),
        ("{V,{V,T}}", bracket(b, bracket(b, a))),
    ]:
        examples.append(
            {
                "name": name,
                "trees": [
                    {"tree": tree, "coefficient": str(c)} for tree, c in expression.terms.items()
                ],
                "polynomial": str(expression.realize(v)),
                "polynomial_terms": polynomial_record(expression.realize(v)),
            }
        )
    (output / "tree_examples.json").write_text(
        json.dumps({"potential": polynomial_record(v), "examples": examples}, indent=2) + "\n"
    )
    checks = {
        "fixed_case_order": 10,
        "prime": 1000003,
        "total_columns": summary["total_dimension"],
        "blocks_replayed": len(blocks),
        "checked_stages": summary["checked_stages"],
        "skipped_stages": summary["skipped_stages"],
        "validation_status": summary["validation_status"],
        "full_polynomial_checks": full_coefficient_checks(),
    }
    (output / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    historical = list(
        csv.DictReader((ROOT / "numerical_case/summaries/graded_dimensions.csv").open())
    )
    rebuilt = list(csv.DictReader((output / "matrices/rank_summary.csv").open()))
    if [r["dimension"] for r in rebuilt] != [r["dimension"] for r in historical[:10]]:
        raise RuntimeError("Fresh lower-order dimensions disagree with the historical scalar table")
    git = shutil.which("git")
    revision = None
    if git:
        status = subprocess.run(
            [git, "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True
        )
        if status.returncode == 0:
            revision = status.stdout.strip()
    source_files = sorted((ROOT / "src").rglob("*.py")) + [Path(__file__).resolve()]
    manifest = {
        "case": "fixed symbolic demonstration",
        "python": platform.python_version(),
        "actual_parameters": {
            "reconstruction_max_order": 10,
            "prime": 1000003,
            "coefficient_checks_max_order": 6,
            "dimensions": [1, 2],
            "potential_degree": 4,
            "workers": 1,
        },
        "commit_observed": revision,
        "source_sha256": {
            f.relative_to(ROOT).as_posix(): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in source_files
        },
        "compact_output_sha256": {
            name: hashlib.sha256((output / name).read_bytes()).hexdigest()
            for name in ["checks.json", "sparsity.csv", "matrix_pattern.json", "tree_examples.json"]
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "columns": checks["total_columns"],
                "blocks": len(blocks),
                "full_polynomial_comparisons": sum(
                    x["comparisons"] for x in checks["full_polynomial_checks"]
                ),
                "status": checks["validation_status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    compute(parser.parse_args().output)
