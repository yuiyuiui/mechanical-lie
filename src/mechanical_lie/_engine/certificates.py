"""Replay saved matrices without using the mechanical-tree construction."""

from __future__ import annotations

import hashlib
import json
from .errors import ValidationError
from .linear import is_prime


def digest(payload: dict) -> str:
    """Hash canonical JSON, excluding the hash field itself."""
    data = {k: v for k, v in payload.items() if k != "sha256_without_hash"}
    return hashlib.sha256(
        json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def verify_block(data: dict) -> int:
    """Check metadata, integrity and every elimination pivot; return the rank.

    This certifies a supplied integer matrix. Its construction from brackets
    is checked separately by rebuilding all columns with run_audit.
    """
    if digest(data) != data["sha256_without_hash"]:
        raise ValidationError("Certificate digest mismatch")
    p = data["prime"]
    if not is_prime(p):
        raise ValidationError("Certificate modulus must be prime")
    count = data["dimension"]
    if not count == len(data["columns"]) == len(data["pivots"]) == len(data["source_basis"]):
        raise ValidationError("Certificate dimensions disagree")
    pivots = {}
    for column, certificate in zip(data["columns"], data["pivots"]):
        v = {row: c % p for row, c in column.items() if c % p}
        while v:
            row = min(v)
            value = v[row]
            if row not in pivots:
                if certificate != {"pivot": row, "value": value}:
                    raise ValidationError("Elimination pivot mismatch")
                inverse = pow(value, -1, p)
                pivots[row] = {r: c * inverse % p for r, c in v.items()}
                break
            for r, c in pivots[row].items():
                new = (v.get(r, 0) - value * c) % p
                if new:
                    v[r] = new
                else:
                    v.pop(r, None)
        else:
            raise ValidationError("A column is dependent")
    if not len(pivots) == count == data["rank"]:
        raise ValidationError("Reported full rank is incorrect")
    return count
