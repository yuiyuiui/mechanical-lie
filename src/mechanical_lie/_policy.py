"""Shared numerical choices; exact outputs are never truncated to fit a policy."""

from ._engine.linear import is_prime

DEFAULT_MAX_ASSIGNMENTS = 100_000
DEFAULT_PRIME = 1_000_003
MAX_RECONSTRUCTION_ORDER = 18


def assignment_limit(value: int) -> tuple[int, str]:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("max_assignments must be an integer; no accuracy/auto mode is supported")
    if value < 1:
        raise ValueError("max_assignments must be positive")
    if value == DEFAULT_MAX_ASSIGNMENTS:
        return value, (
            "Standard per-tree safeguard of 100000 assignments; an engineering limit, "
            "not an accuracy setting or a total memory/time guarantee."
        )
    return value, "Explicit per-tree safeguard; accepted calculations remain complete and exact."


def prime_choice(value: int) -> tuple[int, str]:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("prime must be an integer; a finite field cannot be inferred from data")
    if not is_prime(value):
        raise ValueError("prime must be prime (at least 2)")
    if value == DEFAULT_PRIME:
        return value, (
            "Fixed standard prime 1000003, balancing modest integer arithmetic with a "
            "non-tiny field; no matrix-dependent prime search or rational-rank claim."
        )
    return value, "Explicit prime defines the finite field; rank may differ from other fields."


def marker_for_order(max_order: int) -> int:
    if isinstance(max_order, bool) or not isinstance(max_order, int):
        raise TypeError("max_order must be an integer")
    if not 2 <= max_order <= MAX_RECONSTRUCTION_ORDER:
        raise ValueError(f"max_order must be between 2 and {MAX_RECONSTRUCTION_ORDER}")
    return max_order // 2 + 1
