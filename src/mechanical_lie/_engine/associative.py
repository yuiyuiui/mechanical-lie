"""Sparse integer associative words, Lyndon bracketing and exact degree lowering."""

from __future__ import annotations
from collections import defaultdict
from functools import lru_cache
from .errors import ValidationError, cached_mapping
from typing import Iterator
from .magmatic import A, Y, leaves, star, weight


def clean(p: dict) -> dict:
    return {k: v for k, v in p.items() if v}


def add(a: dict, b: dict, scale: int = 1) -> dict:
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + scale * v
        if not out[k]:
            del out[k]
    return out


def multiply(a: dict, b: dict) -> dict:
    out = defaultdict(int)
    for w, c in a.items():
        for v, d in b.items():
            out[w + v] += c * d
    return clean(out)


def commutator(a: dict, b: dict) -> dict:
    return add(multiply(a, b), multiply(b, a), -1)


@lru_cache(None)
def lyndon(w: tuple) -> bool:
    return bool(w) and all((w < w[i:] + w[:i] for i in range(1, len(w))))


@lru_cache(None)
def standard_split(w: tuple) -> tuple:
    if len(w) < 2 or not lyndon(w):
        raise ValueError("nontrivial Lyndon word required")
    for i in range(1, len(w)):
        if lyndon(w[i:]):
            if not lyndon(w[:i]):
                raise ValidationError("Invariant failed: lyndon(w[:i])")
            return (w[:i], w[i:])
    raise AssertionError("missing Lyndon factorization")


@cached_mapping
def lie_expansion(w: tuple) -> dict:
    if len(w) == 1:
        return {w: 1}
    u, v = standard_split(w)
    return commutator(lie_expansion(u), lie_expansion(v))


def words_of_weight(n: int, alphabet: tuple) -> Iterator[tuple]:
    if n == 0:
        yield ()
        return
    for g in alphabet:
        k = weight(g)
        if k <= n:
            for v in words_of_weight(n - k, alphabet):
                yield ((g,) + v)


@cached_mapping
def normal_x_word(x: tuple, w: tuple) -> dict:
    if not w:
        return {((), x): 1}
    a, tail = (w[0], w[1:])
    out = defaultdict(int)
    for (v, z), c in normal_x_word(x, tail).items():
        out[(a,) + v, z] += c
    if a == A:
        out[(Y(x),) + tail, None] += 1
    else:
        for k, c in normal_x_word(star(x, a[1]), tail).items():
            out[k] += c
    return clean(out)


def adjoint_x(x: tuple, f: dict) -> dict:
    out = defaultdict(int)
    for w, c in f.items():
        for k, d in normal_x_word(x, w).items():
            out[k] += c * d
        out[w, x] -= c
    return clean(out)


def derivative_A(f: dict, x: tuple) -> dict:
    out = defaultdict(int)
    for w, c in f.items():
        for i, g in enumerate(w):
            if g == A:
                out[w[:i] + (Y(x),) + w[i + 1 :]] += c
    return clean(out)


def project_large_labels(f: dict, M: int) -> dict:
    return {w: c for w, c in f.items() if all((g == A or leaves(g[1]) <= M for g in w))}


def retract_marker(f: dict, x: tuple) -> dict:
    out = defaultdict(int)
    for w, c in f.items():
        out[tuple((A if g == Y(x) else g for g in w))] += c
    return clean(out)


def left_action_word(w: tuple, marker: tuple) -> tuple:
    t = marker
    for g in reversed(w):
        if not g[0] == "Y":
            raise ValidationError("Invariant failed: g[0] == 'Y'")
        t = star(g[1], t)
    return t


def decode_marker(t: tuple, M: int) -> tuple:
    out = []
    while leaves(t) > M:
        if not (t and len(t) == 2):
            raise ValidationError("Invariant failed: t and len(t) == 2")
        a, b = t
        large = [u for u in (a, b) if leaves(u) >= M]
        if not len(large) == 1:
            raise ValidationError("Invariant failed: len(large) == 1")
        main = large[0]
        side = b if main == a else a
        if not leaves(side) < M:
            raise ValidationError("Invariant failed: leaves(side) < M")
        out.append(Y(side))
        t = main
    if not leaves(t) == M:
        raise ValidationError("Invariant failed: leaves(t) == M")
    return tuple(out)
