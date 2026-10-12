"""Exact finite-field rank for integer sparse columns."""

from collections.abc import Mapping
from ._engine.linear import sparse_rank
from .polynomial import positive_int


def modular_rank(columns, prime: int = 1000003) -> int:
    """Rank over F_prime of integer columns whose row labels are strings.

    Full modular column rank implies rational independence. Modular deficiency
    alone does not imply rational dependence. No rational rank is returned.
    """
    positive_int(prime, "prime")
    checked = []
    for column in columns:
        if not isinstance(column, Mapping):
            raise TypeError("Each column must be a mapping")
        if any(not isinstance(row, str) or type(c) is not int for row, c in column.items()):
            raise TypeError("Row labels must be strings and matrix coefficients must be integers")
        checked.append(dict(column))
    return sparse_rank(checked, prime)[0]
