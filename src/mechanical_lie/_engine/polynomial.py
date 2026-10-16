"""Independent sparse-polynomial Poisson calibration and iterated lowering checks.

Sparse polynomials map exponent tuples in (q_1,...,q_d,p_1,...,p_d)
to exact rational coefficients. Only the Python standard library is required.
"""

from __future__ import annotations
import itertools
import math
import random
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from . import core as e
from .errors import ValidationError


def pmul(p, q):
    out = defaultdict(Fraction)
    for a, c in p.items():
        for b, d in q.items():
            out[tuple((x + y for x, y in zip(a, b)))] += c * d
    return e.clean(out)


def pdiff(p, i):
    out = {}
    for a, c in p.items():
        if a[i]:
            aa = list(a)
            aa[i] -= 1
            out[tuple(aa)] = c * a[i]
    return out


def poisson(p, q, d):
    """Contract monomial pairs directly, retaining exact coefficients.

    Integer-valued rational coefficients stay integers during contraction;
    the public Polynomial constructor restores its uniform Fraction interface.
    No derivative or product polynomial is materialized for each coordinate.
    """

    def scalar(c):
        return c.numerator if isinstance(c, Fraction) and c.denominator == 1 else c

    left = [(a, scalar(c)) for a, c in p.items()]
    right = [(b, scalar(c)) for b, c in q.items()]
    out = {}
    for a, c in left:
        for b, coefficient in right:
            product = None
            powers = None
            for i in range(d):
                factor = a[i] * b[i + d] - a[i + d] * b[i]
                if not factor:
                    continue
                if product is None:
                    product = c * coefficient
                    powers = tuple(x + y for x, y in zip(a, b))
                alpha = list(powers)
                alpha[i] -= 1
                alpha[i + d] -= 1
                key = tuple(alpha)
                out[key] = out.get(key, 0) + product * factor
    return e.clean(out)


def peval(p, point):
    return sum((c * math.prod((x**a for x, a in zip(point, alpha))) for alpha, c in p.items()))


def exponents(d, max_degree):
    for alpha in itertools.product(range(max_degree + 1), repeat=d):
        if sum(alpha) <= max_degree:
            yield alpha


def polynomial_calibration(max_order=8, d=2, seed=519, *, details=None):
    """Compare every input column with direct Poisson differentiation at a jet.

    If supplied, ``details`` is populated with the full potential, evaluation
    point and per-column values, using strings for exact rational numbers.
    """
    rng = random.Random(seed)
    v = {a + (0,) * d: Fraction(rng.choice([-3, -2, -1, 1, 2, 3])) for a in exponents(d, 5)}
    kinetic = {}
    for i in range(d):
        a = [0] * (2 * d)
        a[d + i] = 2
        kinetic[tuple(a)] = Fraction(1, 2)
    p0 = tuple(range(1, d + 1))
    point = (0,) * d + p0
    jet = {a[:d]: c * math.prod((math.factorial(x) for x in a[:d])) for a, c in v.items()}
    if details is not None:
        details.update(
            variables=[f"q{i + 1}" for i in range(d)] + [f"p{i + 1}" for i in range(d)],
            potential=[{"exponents": a, "coefficient": str(c)} for a, c in sorted(v.items())],
            evaluation_point=[str(x) for x in point],
            comparisons=[],
        )

    @lru_cache(maxsize=4096)
    def px(t):
        if t == ():
            return v
        a, b = t
        return poisson(px(b), poisson(px(a), kinetic, d), d)

    @lru_cache(maxsize=4096)
    def pg(g):
        return kinetic if g == e.A else poisson(px(g[1]), kinetic, d)

    @lru_cache(maxsize=4096)
    def pw(w):
        if len(w) == 1:
            return pg(w[0])
        a, b = e.standard_split(w)
        return poisson(pw(a), pw(b), d)

    @lru_cache(None)
    def eval_tree(t):
        if t is None:
            return Fraction(sum((p * p for p in p0)), 2)
        colors, adj = e.decode_mechanical(t)
        k = len(colors)
        core_edges = [(i, j) for i in range(k) for j in adj[i] if i < j]
        free_edges = [i for i in range(k) for _ in range(colors[i])]
        total = 0
        for labels in itertools.product(range(d), repeat=len(core_edges) + len(free_edges)):
            inds = [[0] * d for _ in range(k)]
            for label, (i, j) in zip(labels, core_edges):
                inds[i][label] += 1
                inds[j][label] += 1
            momenta = 1
            for label, i in zip(labels[len(core_edges) :], free_edges):
                inds[i][label] += 1
                momenta *= p0[label]
            total += momenta * math.prod((jet.get(tuple(a), 0) for a in inds))
        return total

    alphabet = tuple(
        sorted([e.A] + [e.Y(t) for n in range(1, max_order // 2 + 1) for t in e.magmatic(n)])
    )
    checks = 0
    for n in range(1, max_order + 1):
        for w in e.words_of_weight(n, alphabet):
            if not e.lyndon(w):
                continue
            direct = peval(pw(w), point)
            tree = sum((c * eval_tree(t) for t, c in e.theta_lyndon(w).items()))
            if not direct == tree:
                raise ValidationError((n, w, direct, tree))
            if details is not None:
                details["comparisons"].append(
                    {
                        "order": n,
                        "kind": "Lyndon",
                        "input": repr(w),
                        "direct": str(direct),
                        "tree": str(tree),
                    }
                )
            checks += 1
        if n % 2:
            for x in e.magmatic((n + 1) // 2):
                direct = peval(px(x), point)
                tree = sum((c * eval_tree(t) for t, c in e.theta_X(x).items()))
                if not direct == tree:
                    raise ValidationError((n, x, direct, tree))
                if details is not None:
                    details["comparisons"].append(
                        {
                            "order": n,
                            "kind": "potential",
                            "input": repr(x),
                            "direct": str(direct),
                            "tree": str(tree),
                        }
                    )
                checks += 1
    return {"dimension": d, "potential_degree": 5, "seed": seed, "checks": checks}


@lru_cache(None)
def balanced(n):
    if n == 1:
        return ()
    return e.star(balanced(n // 2), balanced(n - n // 2))


def iterated_lowering(max_order=8):
    alphabet = tuple(
        sorted([e.A] + [e.Y(t) for n in range(1, max_order // 2 + 1) for t in e.magmatic(n)])
    )
    traces = []
    for n in range(1, max_order + 1):
        for w in e.words_of_weight(n, alphabet):
            if not e.lyndon(w):
                continue
            f = e.lie_expansion(w)
            degree = 1 + sum((g == e.A for g in w))
            initial_degree = degree
            steps = []
            while degree > 0:
                max_weight = max((sum((e.leaves(g[1]) for g in word if g != e.A)) for word in f))
                marker = balanced(max_weight + 1)
                value = e.adjoint_x(marker, f)
                if not value:
                    raise ValidationError("Invariant failed: value")
                steps.append(
                    {
                        "degree_before": degree,
                        "marker_leaves": max_weight + 1,
                        "pbw_input_support": len(f),
                        "pbw_output_support": len(value),
                    }
                )
                if degree == 1:
                    if not all((not word and z is not None for word, z in value)):
                        raise ValidationError(
                            "Invariant failed: all((not word and z is not None for word, z in value))"
                        )
                    final = {z: c for (word, z), c in value.items()}
                    if not final:
                        raise ValidationError("Invariant failed: final")
                else:
                    if not all((z is None for word, z in value)):
                        raise ValidationError(
                            "Invariant failed: all((z is None for word, z in value))"
                        )
                    f = {word: c for (word, z), c in value.items()}
                degree -= 1
            traces.append(
                {
                    "order": n,
                    "initial_degree": initial_degree,
                    "steps": steps,
                    "nonzero_final_potential_terms": len(final),
                }
            )
    return traces
