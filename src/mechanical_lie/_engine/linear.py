"""Finite-field column elimination. Full rank modulo a prime implies full rank over the rationals."""

from __future__ import annotations


def sparse_rank(columns: list[dict], prime: int) -> tuple:
    if not is_prime(prime):
        raise ValueError("A prime modulus is required")
    pivots = {}
    cert = []
    for column in columns:
        v = {repr(t): c % prime for t, c in column.items() if c % prime}
        while v:
            k = min(v)
            a = v[k]
            if k not in pivots:
                inv = pow(a, -1, prime)
                pivots[k] = {j: c * inv % prime for j, c in v.items()}
                cert.append({"pivot": k, "value": a})
                break
            for j, c in pivots[k].items():
                value = (v.get(j, 0) - a * c) % prime
                if value:
                    v[j] = value
                elif j in v:
                    del v[j]
        else:
            cert.append(None)
    return (len(pivots), cert)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True
