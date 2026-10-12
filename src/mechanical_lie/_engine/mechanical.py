"""Exact Euclidean mechanical-tree contractions with the Poisson sign f_q g_p - f_p g_q."""

from __future__ import annotations
from collections import defaultdict
from functools import lru_cache
from .errors import cached_mapping
from .magmatic import A
from .associative import clean, standard_split


@lru_cache(None)
def decode_mechanical(code: tuple) -> tuple:
    colors = []
    adj = []

    def walk(node, parent=None):
        i = len(colors)
        colors.append(node[0])
        adj.append([])
        if parent is not None:
            adj[i].append(parent)
            adj[parent].append(i)
        for child in node[1]:
            walk(child, i)

    walk(code)
    return (tuple(colors), tuple((tuple(a) for a in adj)))


def canonical(colors: tuple, adj: tuple) -> tuple:
    validate_graph(colors, adj)

    def root(i, parent):
        return (colors[i], tuple(sorted((root(j, i) for j in adj[i] if j != parent))))

    return min((root(i, -1) for i in range(len(colors))))


@cached_mapping
def add_momentum(code: tuple) -> dict:
    colors, adj = decode_mechanical(code)
    out = defaultdict(int)
    for i in range(len(colors)):
        cc = list(colors)
        cc[i] += 1
        out[canonical(tuple(cc), adj)] += 1
    return dict(out)


@cached_mapping
def tree_bracket(u, v) -> dict:
    if u == v:
        return {}
    if u is None:
        return {t: -c for t, c in add_momentum(v).items()}
    if v is None:
        return add_momentum(u)
    cu, au = decode_mechanical(u)
    cv, av = decode_mechanical(v)
    if not any(cu) and (not any(cv)):
        return {}
    k = len(cu)
    base_c = cu + cv
    base_a = [list(z) for z in au] + [[j + k for j in z] for z in av]
    out = defaultdict(int)
    for i in range(k):
        for j in range(len(cv)):
            if cu[i] == 0 and cv[j] == 0:
                continue
            aa = [z[:] for z in base_a]
            aa[i].append(k + j)
            aa[k + j].append(i)
            adj = tuple((tuple(z) for z in aa))
            if cv[j]:
                cc = list(base_c)
                cc[k + j] -= 1
                out[canonical(tuple(cc), adj)] += cv[j]
            if cu[i]:
                cc = list(base_c)
                cc[i] -= 1
                out[canonical(tuple(cc), adj)] -= cu[i]
    return clean(out)


def mechanical_bracket(p: dict, q: dict) -> dict:
    out = defaultdict(int)
    for u, c in p.items():
        for v, d in q.items():
            for t, e in tree_bracket(u, v).items():
                out[t] += c * d * e
    return clean(out)


@cached_mapping
def theta_X(t: tuple) -> dict:
    if not t:
        return {(0, ()): 1}
    x, z = t
    return mechanical_bracket(theta_X(z), mechanical_bracket(theta_X(x), {None: 1}))


@cached_mapping
def theta_generator(g: tuple) -> dict:
    return {None: 1} if g == A else mechanical_bracket(theta_X(g[1]), {None: 1})


@cached_mapping
def theta_lyndon(w: tuple) -> dict:
    if len(w) == 1:
        return theta_generator(w[0])
    u, v = standard_split(w)
    return mechanical_bracket(theta_lyndon(u), theta_lyndon(v))


def tree_degree(t) -> int:
    return 2 if t is None else sum(decode_mechanical(t)[0])


def tree_order(t) -> int:
    return 1 if t is None else 2 * len(decode_mechanical(t)[0]) + tree_degree(t) - 1


def validate_graph(colors: tuple, adj: tuple) -> None:
    """Reject graphs that are not finite nonempty simple trees."""
    n = len(colors)
    if not n or len(adj) != n or any((not isinstance(c, int) or c < 0 for c in colors)):
        raise ValueError("Nonnegative integer colors and a nonempty graph required")
    for i, row in enumerate(adj):
        if len(set(row)) != len(row) or any(
            (not isinstance(j, int) or not 0 <= j < n or j == i for j in row)
        ):
            raise ValueError("Invalid adjacency")
        if any((i not in adj[j] for j in row)):
            raise ValueError("Adjacency must be symmetric")
    reached = {0}
    stack = [0]
    while stack:
        for j in adj[stack.pop()]:
            if j not in reached:
                reached.add(j)
                stack.append(j)
    if len(reached) != n or sum(map(len, adj)) != 2 * (n - 1):
        raise ValueError("Graph must be a connected tree")
