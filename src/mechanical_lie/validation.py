"""Rebuild finite basis matrices and replay integer certificates."""

from pathlib import Path
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from ._engine.audit import run_audit
from ._engine.certificates import digest, verify_block
from ._engine.errors import ValidationError
from ._policy import DEFAULT_PRIME, marker_for_order, prime_choice


@dataclass(frozen=True)
class ReconstructionPlan:
    """Chosen finite target and the deterministic marker/prime policy."""

    max_order: int
    marker_leaves: int
    prime: int
    prime_reason: str
    marker_reason: str
    verification: str = "all-stages"

    def to_dict(self) -> dict:
        return asdict(self)


def reconstruction_plan(max_order: int, *, prime: int = DEFAULT_PRIME) -> ReconstructionPlan:
    """Inspect configuration without enumerating a basis or creating files."""
    marker = marker_for_order(max_order)
    prime, reason = prime_choice(prime)
    return ReconstructionPlan(
        max_order,
        marker,
        prime,
        reason,
        "floor(max_order/2)+1 strictly exceeds the source side-label leaf budget; "
        "the requested order is never reduced to make computation cheaper.",
    )


def rebuild_basis(max_order: int, output: str | Path, *, prime: int = DEFAULT_PRIME) -> dict:
    """Reconstruct and validate every basis column through order 2--18.

    The output directory must be new. Higher orders grow rapidly; the quick
    examples use order 10 or below. This is finite validation, not a theorem.
    """
    plan = reconstruction_plan(max_order, prime=prime)
    result = run_audit(plan.max_order, Path(output), plan.prime)
    result["effective_configuration"] = plan.to_dict()
    # Keep the saved summary and the returned diagnostics identical.
    import json

    (Path(output) / "audit_summary.json").write_text(json.dumps(result, indent=2), encoding="utf8")
    return result


def verify_certificate_block(data: Mapping) -> int:
    """Check a supplied full-column-rank integer matrix and all saved pivots.

    This does not reconstruct the brackets that produced the matrix. Missing
    fields and malformed records raise ValidationError, not silent success.
    """
    required = {
        "order",
        "degree",
        "prime",
        "dimension",
        "rank",
        "columns",
        "pivots",
        "source_basis",
        "sha256_without_hash",
    }
    if not isinstance(data, Mapping) or not required <= data.keys():
        raise ValidationError("Missing certificate fields")
    if any(type(data[k]) is not int for k in ("order", "degree", "prime", "dimension", "rank")):
        raise ValidationError("Certificate metadata must contain integers")
    if data["order"] < 1 or data["degree"] < 0 or data["dimension"] < 1:
        raise ValidationError("Invalid certificate order, degree or dimension")
    if any(not isinstance(data[k], list) for k in ("columns", "pivots", "source_basis")):
        raise ValidationError("Certificate columns, pivots and basis must be lists")
    for column in data["columns"]:
        if not isinstance(column, dict) or any(
            not isinstance(row, str) or type(value) is not int for row, value in column.items()
        ):
            raise ValidationError("Matrix columns must map string rows to integer coefficients")
    try:
        return verify_block(dict(data))
    except (KeyError, TypeError, ValueError) as error:
        if isinstance(error, ValidationError):
            raise
        raise ValidationError(f"Invalid certificate: {error}") from error


__all__ = [
    "ReconstructionPlan",
    "digest",
    "rebuild_basis",
    "reconstruction_plan",
    "verify_certificate_block",
]
