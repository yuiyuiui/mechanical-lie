"""Paper-scale exact algebra audits; finite checks, separate from the all-order proof."""

from __future__ import annotations
import csv
import gzip
import hashlib
import json
import time
from pathlib import Path
from collections import defaultdict
from .core import (
    A,
    B,
    Y,
    add,
    adjoint_x,
    clean,
    comb,
    decode_marker,
    derivative_A,
    is_prime,
    left_action_word,
    lie_expansion,
    lyndon,
    magmatic,
    mechanical_bracket,
    project_large_labels,
    retract_marker,
    sparse_rank,
    star,
    theta_X,
    theta_lyndon,
    tree_degree,
    tree_order,
    words_of_weight,
)
from .errors import ValidationError
from .._policy import DEFAULT_PRIME, marker_for_order


def run_audit(
    max_order: int, output: Path, prime: int = DEFAULT_PRIME, skip_tree_rank: bool = False
) -> dict:
    """Rebuild every basis column through max_order and save replayable certificates."""
    M = marker_for_order(max_order)
    if not is_prime(prime):
        raise ValueError("prime modulus required")
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    N = max_order
    alphabet = tuple(sorted([A] + [Y(t) for n in range(1, N // 2 + 1) for t in magmatic(n)]))
    marker = comb(M)
    b = {(0, ()): 1}
    a = {None: 1}
    ba = mechanical_bracket(b, a)
    if not ba == {(1, ()): 1}:
        raise ValidationError("Invariant failed: ba == {(1, ()): 1}")
    bba = mechanical_bracket(b, ba)
    if not bba == theta_X(star(B, B)):
        raise ValidationError("Invariant failed: bba == theta_X(star(B, B))")
    if not mechanical_bracket(b, bba) == {}:
        raise ValidationError("Invariant failed: mechanical_bracket(b, bba) == {}")
    calibrations = [b, a, ba, bba, mechanical_bracket(ba, a)]
    for p in calibrations:
        for q in calibrations:
            if not mechanical_bracket(p, q) == {t: -c for t, c in mechanical_bracket(q, p).items()}:
                raise ValidationError(
                    "Invariant failed: mechanical_bracket(p, q) == {t: -c for t, c in mechanical_bracket(q, p).items()}"
                )
            for r in calibrations:
                jac = add(
                    add(
                        mechanical_bracket(p, mechanical_bracket(q, r)),
                        mechanical_bracket(q, mechanical_bracket(r, p)),
                    ),
                    mechanical_bracket(r, mechanical_bracket(p, q)),
                )
                if not not jac:
                    raise ValidationError("Invariant failed: not jac")
    counts = {
        "degree_one_action": 0,
        "higher_degree_projection": 0,
        "retraction": 0,
        "pbw_purity": 0,
        "marker_word_decoding": 0,
        "jacobi_triples": 125,
    }
    y_alphabet = tuple((g for g in alphabet if g != A))
    coded = {}
    for n in range(0, N + 1, 2):
        for w in words_of_weight(n, y_alphabet):
            t = left_action_word(w, marker)
            if not decode_marker(t, M) == w:
                raise ValidationError("Invariant failed: decode_marker(t, M) == w")
            if not (t not in coded or coded[t] == w):
                raise ValidationError("Invariant failed: t not in coded or coded[t] == w")
            coded[t] = w
            counts["marker_word_decoding"] += 1
    table = []
    rank_certificates = []
    for n in range(1, N + 1):
        basis = [w for w in words_of_weight(n, alphabet) if lyndon(w)]
        blocks = defaultdict(list)
        block_names = defaultdict(list)
        for w in basis:
            f = lie_expansion(w)
            if not (f and f.get(w) == 1):
                raise ValidationError("Invariant failed: f and f.get(w) == 1")
            k = sum((g == A for g in w))
            ax = adjoint_x(marker, f)
            if not ax:
                raise ValidationError("Invariant failed: ax")
            if k == 0:
                if not all((not word and z is not None for word, z in ax)):
                    raise ValidationError(
                        "Invariant failed: all((not word and z is not None for word, z in ax))"
                    )
                expected = defaultdict(int)
                for v, c in f.items():
                    expected[(), left_action_word(v, marker)] -= c * (-1) ** len(v)
                if not ax == clean(expected):
                    raise ValidationError("Invariant failed: ax == clean(expected)")
                counts["degree_one_action"] += 1
            else:
                if not all((z is None for word, z in ax)):
                    raise ValidationError("Invariant failed: all((z is None for word, z in ax))")
                pure = {word: c for (word, z), c in ax.items()}
                derivative = derivative_A(f, marker)
                if not project_large_labels(pure, M) == derivative:
                    raise ValidationError(
                        "Invariant failed: project_large_labels(pure, M) == derivative"
                    )
                if not retract_marker(derivative, marker) == {v: k * c for v, c in f.items()}:
                    raise ValidationError(
                        "Invariant failed: retract_marker(derivative, marker) == {v: k * c for v, c in f.items()}"
                    )
                counts["higher_degree_projection"] += 1
                counts["retraction"] += 1
            counts["pbw_purity"] += 1
            if not skip_tree_rank:
                image = theta_lyndon(w)
                if not image:
                    raise ValidationError("Invariant failed: image")
                if not all((tree_degree(t) == k + 1 and tree_order(t) == n for t in image)):
                    raise ValidationError(
                        "Invariant failed: all((tree_degree(t) == k + 1 and tree_order(t) == n for t in image))"
                    )
                blocks[k + 1].append(image)
                block_names[k + 1].append(repr(w))
        if n % 2:
            for x in magmatic((n + 1) // 2):
                if not skip_tree_rank:
                    blocks[0].append(theta_X(x))
                    block_names[0].append(repr(x))
        dimension = len(basis) + (len(magmatic((n + 1) // 2)) if n % 2 else 0)
        ranks = {}
        if not skip_tree_rank:
            for d, cols in sorted(blocks.items()):
                rank, cert = sparse_rank(cols, prime)
                if not rank == len(cols):
                    raise ValidationError((n, d, rank, len(cols)))
                ranks[d] = rank
                serialized = [
                    {repr(t): c for t, c in sorted(col.items(), key=lambda z: repr(z[0]))}
                    for col in cols
                ]
                payload = {
                    "order": n,
                    "degree": d,
                    "prime": prime,
                    "dimension": len(cols),
                    "rank": rank,
                    "columns": serialized,
                    "pivots": cert,
                    "source_basis": block_names[d],
                    "basis_kind": "ComMag binary trees" if d == 0 else "weighted Lyndon words",
                }
                canonical_json = json.dumps(payload, sort_keys=True, separators=(",", ":"))
                payload["sha256_without_hash"] = hashlib.sha256(canonical_json.encode()).hexdigest()
                rank_certificates.append(payload)
        row = {
            "order": n,
            "dimension": dimension,
            "positive_dimension": len(basis),
            "rank": sum(ranks.values()) if ranks else None,
            "largest_column_support": max(
                (len(col) for cols in blocks.values() for col in cols), default=0
            ),
        }
        table.append(row)
        print(json.dumps(row), flush=True)
    with (out / "rank_summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    if rank_certificates:
        with gzip.open(out / "rank_certificates.json.gz", "wt", encoding="utf8") as fh:
            json.dump(rank_certificates, fh, ensure_ascii=False)
    result = {
        "max_order": N,
        "marker_leaves": M,
        "prime": prime,
        "checks": counts,
        "total_dimension": sum((r["dimension"] for r in table)),
        "total_positive_dimension": sum((r["positive_dimension"] for r in table)),
        "exact_integer_arithmetic": True,
        "all_checks_passed": None if skip_tree_rank else True,
        "all_executed_checks_passed": True,
        "validation_status": "PARTIAL" if skip_tree_rank else "COMPLETE",
        "checked_stages": [
            "generator signs",
            "antisymmetry",
            "Jacobi",
            "PBW lowering",
            "marker decoding",
        ]
        + (
            []
            if skip_tree_rank
            else ["tree reconstruction", "modular rank", "certificate generation"]
        ),
        "skipped_stages": ["tree reconstruction", "modular rank", "certificate generation"]
        if skip_tree_rank
        else [],
        "elapsed_seconds": round(time.time() - start, 3),
        "scope": "Finite audits only. The all-order argument is separate.",
    }
    (out / "audit_summary.json").write_text(json.dumps(result, indent=2), encoding="utf8")
    print(json.dumps(result, indent=2), flush=True)
    return result
