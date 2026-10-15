"""Bounded matched exact-output comparisons of resource and finite-field policies."""

from __future__ import annotations

import argparse
import contextlib
from fractions import Fraction
import gzip
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import tempfile
import time
import tracemalloc

THREAD_CAPS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
for _key in THREAD_CAPS:
    os.environ[_key] = "1"

import sympy as sp

from mechanical_lie import (
    Polynomial,
    __version__,
    bracket,
    kinetic,
    modular_rank_report,
    phase_space,
    poisson,
    potential,
)
from mechanical_lie.validation import rebuild_basis, verify_certificate_block

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPEATS = 3
MAX_REPEATS = 20


def output_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def polynomial_terms(value):
    return [[list(powers), str(c)] for powers, c in sorted(value.terms.items())]


def measure_variants(variants, repeats, validate):
    """Warm each choice, rotate trial order, and time the whole API/pipeline call."""
    samples = [[] for _ in variants]
    peaks = [[] for _ in variants]
    for _, operation, _ in variants:
        validate(operation())
    for trial in range(repeats):
        for offset in range(len(variants)):
            index = (offset + trial) % len(variants)
            operation = variants[index][1]
            tracemalloc.start()
            start = time.perf_counter()
            value = operation()
            elapsed = time.perf_counter() - start
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            validate(value)
            samples[index].append(elapsed)
            peaks[index].append(peak)
    return [
        dict(
            configuration=name,
            effective_configuration=config,
            elapsed_seconds=samples[i],
            median_elapsed_seconds=statistics.median(samples[i]),
            peak_traced_bytes=max(peaks[i]),
            validation="PASS",
        )
        for i, (name, _, config) in enumerate(variants)
    ]


def realization_cases(repeats):
    records = []
    for dimension in (1, 2, 3):
        setup_start = time.perf_counter()
        q, _ = phase_space(dimension)
        v = sum((x**2 / 2 + x**4 / 4 for x in q), Polynomial.constant(dimension, 0))
        v += sum(
            (q[i] * q[i + 1] / 3 for i in range(dimension - 1)), Polynomial.constant(dimension, 0)
        )
        h, direct = potential(), v
        for generator in (kinetic(), kinetic(), potential(), kinetic()):
            h = bracket(generator, h)
        input_setup = time.perf_counter() - setup_start
        t = Polynomial.kinetic(dimension)
        for generator in (t, t, v, t):
            direct = poisson(generator, direct)
        # A second independent implementation checks every coefficient.
        variables = sp.symbols(
            " ".join([f"q{i}" for i in range(dimension)] + [f"p{i}" for i in range(dimension)])
        )
        sv = sum(x * x / 2 + x**4 / 4 for x in variables[:dimension])
        sv += sum(variables[i] * variables[i + 1] / 3 for i in range(dimension - 1))
        st = sum(x * x / 2 for x in variables[dimension:])
        reference = sv
        for generator in (st, st, sv, st):
            reference = sp.expand(
                sum(
                    sp.diff(generator, variables[i]) * sp.diff(reference, variables[dimension + i])
                    - sp.diff(generator, variables[dimension + i])
                    * sp.diff(reference, variables[i])
                    for i in range(dimension)
                )
            )
        expected = Polynomial(
            dimension,
            {
                powers: c
                for powers, c in (
                    (powers, Fraction(int(c.p), int(c.q)))
                    for powers, c in sp.Poly(reference, *variables).terms()
                )
            },
        )
        if direct != expected:
            raise RuntimeError("Independent polynomial references disagree")
        plan = h.realization_plan(v)
        configurations = [
            ("default", {}),
            ("minimum-admitting-cap", {"max_assignments": plan.largest_assignments}),
            ("larger-explicit-cap", {"max_assignments": 200000}),
        ]
        variants = [
            (name, lambda kw=kwargs: h.realize(v, **kw), h.realization_plan(v, **kwargs).to_dict())
            for name, kwargs in configurations
        ]

        def validate(value):
            if value != expected:
                raise RuntimeError("Realization differs from complete independent coefficients")

        rows = measure_variants(variants, repeats, validate)
        for row in rows:
            row.update(
                category="realization",
                case=f"coupled-quartic-d{dimension}",
                input_setup_seconds=input_setup,
                dimension=dimension,
                input_potential=polynomial_terms(v),
                source_tree_count=len(h.terms),
                output_terms=len(expected.terms),
                exact_output_sha256=output_hash(polynomial_terms(expected)),
                exact_reference="complete direct Poisson and independent SymPy coefficients",
            )
        records.extend(rows)
    return records


def rank_cases(repeats):
    records = []
    for size in (6, 12):
        columns = [{f"r{j}": 1} | ({f"r{j - 1}": 2} if j else {}) for j in range(size)]
        reference = int(
            sp.Matrix([[column.get(f"r{i}", 0) for column in columns] for i in range(size)]).rank()
        )
        variants = []
        for name, kwargs in [
            ("default", {}),
            ("explicit-small-field", {"prime": 2}),
            ("explicit-prime-101", {"prime": 101}),
        ]:
            report = modular_rank_report(columns, **kwargs)
            config = {k: report.to_dict()[k] for k in ("prime", "backend", "selection_reason")}
            variants.append((name, lambda kw=kwargs: modular_rank_report(columns, **kw), config))

        def validate(value):
            if value.rank != reference or not value.rational_independence_established:
                raise RuntimeError("Full-rank comparison failed on a unimodular integer matrix")

        rows = measure_variants(variants, repeats, validate)
        for row in rows:
            row.update(
                category="modular-rank",
                case=f"unit-triangular-{size}",
                columns=columns,
                rational_reference_rank=reference,
                observed_rank=reference,
                integer_nonzeros=sum(len(c) for c in columns),
                comparison_scope="Same full-column-independence question on unimodular inputs; fields differ",
            )
        records.extend(rows)
    bad_default = modular_rank_report([{"x": 1000003}])
    alternative = modular_rank_report([{"x": 1000003}], prime=101)
    if (
        bad_default.rank != 0
        or bad_default.rational_independence_established
        or alternative.rank != 1
    ):
        raise RuntimeError("Bad-prime interpretation regression")
    return records, {
        "columns": [{"x": 1000003}],
        "rational_reference_rank": 1,
        "default": bad_default.to_dict(),
        "explicit_prime_101": alternative.to_dict(),
        "interpretation": "Deficiency at the default prime is not rational dependence.",
    }


def reconstruction_cases(output, repeats):
    raw = output / "raw"
    raw.mkdir()
    gold = []

    def operation(kwargs):
        with tempfile.TemporaryDirectory(dir=raw) as temp:
            target = Path(temp) / "case"
            with contextlib.redirect_stdout(io.StringIO()):
                result = rebuild_basis(6, target, **kwargs)
            with gzip.open(target / "rank_certificates.json.gz", "rt", encoding="utf8") as stream:
                blocks = json.load(stream)
            for block in blocks:
                verify_certificate_block(block)
            matrices = [[b["order"], b["degree"], b["columns"]] for b in blocks]
            digest = output_hash(matrices)
            if not gold:
                gold.append(digest)
            if (
                digest != gold[0]
                or result["total_dimension"] != 16
                or not result["all_checks_passed"]
            ):
                raise RuntimeError(
                    "Prime changed the integer construction or a finite check failed"
                )
            return result

    variants = []
    for name, kwargs in [("default", {}), ("explicit-prime-101", {"prime": 101})]:
        from mechanical_lie.validation import reconstruction_plan

        variants.append(
            (name, lambda kw=kwargs: operation(kw), reconstruction_plan(6, **kwargs).to_dict())
        )
    rows = measure_variants(variants, repeats, lambda r: None)
    for row in rows:
        row.update(
            category="reconstruction",
            case="finite-order-6",
            max_order=6,
            total_columns=16,
            integer_matrix_sha256=gold[0],
            exact_reference="unchanged integer matrices, known dimensions and separate certificate replay",
        )
    return rows


def benchmark(output: Path, repeats: int = DEFAULT_REPEATS):
    if type(repeats) is not int or not 1 <= repeats <= MAX_REPEATS:
        raise ValueError(f"repeats must be an integer between 1 and {MAX_REPEATS}")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    rank_records, bad_prime = rank_cases(repeats)
    records = realization_cases(repeats) + rank_records + reconstruction_cases(output, repeats)
    git = shutil.which("git")
    revision = (
        subprocess.run([git, "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True)
        if git
        else None
    )
    sources = sorted((ROOT / "src/mechanical_lie").rglob("*.py")) + [Path(__file__).resolve()]
    result = dict(
        schema="mechanical-lie-policy-benchmark-v1",
        version=__version__,
        command=f"python bench/benchmark_policies.py --output <fresh-directory> --repeats {repeats}",
        environment={
            "python": platform.python_version(),
            "platform": platform.platform(),
            "sympy": importlib.metadata.version("sympy"),
            "workers": 1,
            "random_seed": None,
            "randomness": "none",
            "native_thread_caps": {key: os.environ[key] for key in THREAD_CAPS},
        },
        git_commit_observed=revision.stdout.strip()
        if revision and revision.returncode == 0
        else None,
        source_sha256={
            f.relative_to(ROOT).as_posix(): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in sources
        },
        measurement={
            "repeats": repeats,
            "warmups_per_configuration": 1,
            "trial_order": "rotating round-robin on the same process and inputs",
            "time": "instrumented warm API/pipeline calls include policy and input validation; imports, independent references and final result comparison excluded",
            "memory": "peak Python allocations during the call, not total process RSS or a universal memory bound",
            "reconstruction_time": "includes temporary output setup, full rebuild, replay and cleanup",
        },
        records=records,
        bad_prime_control=bad_prime,
        conclusion="All admitted cap choices preserve complete exact outputs. Caps control admission, not precision. The finite-field comparisons do not establish a generally optimal prime or rational dependence.",
    )
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                status="PASS",
                configurations=len(records),
                repeats=repeats,
                bad_prime_control="PASS",
                version=__version__,
            )
        )
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument(
        "--repeats",
        type=int,
        default=DEFAULT_REPEATS,
        help=f"Timed samples per choice (1–{MAX_REPEATS})",
    )
    args = parser.parse_args()
    benchmark(args.output, args.repeats)
