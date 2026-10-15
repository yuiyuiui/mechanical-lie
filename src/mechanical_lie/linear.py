"""Exact finite-field rank for integer sparse columns."""

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from ._engine.linear import sparse_rank
from ._policy import DEFAULT_PRIME, prime_choice


@dataclass(frozen=True)
class ModularRankReport:
    """A finite-field result; deficiency does not establish rational dependence."""

    rank: int
    column_count: int
    row_count: int
    prime: int
    selection_reason: str
    backend: str = "sparse-modular"

    @property
    def rational_independence_established(self) -> bool:
        return self.rank == self.column_count

    def to_dict(self) -> dict:
        return asdict(self) | {
            "rational_independence_established": self.rational_independence_established
        }


def modular_rank_report(columns, prime: int = DEFAULT_PRIME) -> ModularRankReport:
    """Return the same exact rank plus effective field, dimensions and reason."""
    prime, reason = prime_choice(prime)  # Validate before consuming input iterators.
    checked = []
    rows = set()
    for column in columns:
        if not isinstance(column, Mapping):
            raise TypeError("Each column must be a mapping")
        if any(not isinstance(row, str) or type(c) is not int for row, c in column.items()):
            raise TypeError("Row labels must be strings and matrix coefficients must be integers")
        checked.append(dict(column))
        rows.update(row for row, c in column.items() if c)
    rank = sparse_rank(checked, prime)[0]
    return ModularRankReport(rank, len(checked), len(rows), prime, reason)


def modular_rank(columns, prime: int = DEFAULT_PRIME) -> int:
    """Rank over F_prime of integer columns whose row labels are strings.

    Full modular column rank implies rational independence. Modular deficiency
    alone does not imply rational dependence. No rational rank is returned.
    """
    return modular_rank_report(columns, prime).rank
