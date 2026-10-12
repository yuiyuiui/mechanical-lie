"""Canonical nonplanar binary trees for the free nonunital commutative magmatic algebra."""

from __future__ import annotations
from functools import lru_cache

B = ()
A = ("A",)


@lru_cache(None)
def leaves(t: tuple) -> int:
    return 1 if not t else leaves(t[0]) + leaves(t[1])


def star(a: tuple, b: tuple) -> tuple:
    return (a, b) if a <= b else (b, a)


@lru_cache(None)
def magmatic(n: int) -> tuple:
    if n < 1:
        return ()
    if n == 1:
        return (B,)
    ans = set()
    for i in range(1, n // 2 + 1):
        for a in magmatic(i):
            for b in magmatic(n - i):
                ans.add(star(a, b))
    return tuple(sorted(ans))


def comb(n: int) -> tuple:
    if n < 1:
        raise ValueError("positive number of leaves required")
    t = B
    for _ in range(n - 1):
        t = star(B, t)
    return t


def Y(t: tuple) -> tuple:
    return ("Y", t)


def weight(g: tuple) -> int:
    return 1 if g == A else 2 * leaves(g[1])
