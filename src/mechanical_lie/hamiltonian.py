"""Canonical tree expressions for Euclidean mechanical Hamiltonians."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from itertools import product

from ._engine import core as engine
from .polynomial import Polynomial, exact, phase_space, positive_int


def _tree(code):
    if code is None:
        return None

    def check(node):
        if not isinstance(node, tuple) or len(node) != 2:
            raise ValueError("A mechanical tree node must be (momentum_count, children_tuple)")
        color, children = node
        if type(color) is not int or color < 0 or not isinstance(children, tuple):
            raise ValueError(
                "Momentum counts must be nonnegative integers; children must be a tuple"
            )
        for child in children:
            check(child)

    check(code)
    colors, adjacency = engine.decode_mechanical(code)
    return engine.canonical(colors, adjacency)


@dataclass(frozen=True, init=False)
class Hamiltonian:
    """Finite rational linear combination of canonical mechanical trees.

    A node is ``(momentum_count, children_tuple)``; ``None`` denotes kinetic
    energy. Arithmetic is linear; ordinary multiplication is not a supported
    tree operation. Use ``bracket`` or ``commutative_product`` instead.
    """

    _items: tuple

    def __init__(self, terms: Mapping[tuple | None, int | Fraction]):
        if not isinstance(terms, Mapping):
            raise TypeError("terms must be a mapping from mechanical trees to exact coefficients")
        clean = {}
        for tree, coefficient in terms.items():
            tree = _tree(tree)
            clean[tree] = clean.get(tree, 0) + exact(coefficient)
        object.__setattr__(
            self,
            "_items",
            tuple(sorted(((t, c) for t, c in clean.items() if c), key=lambda item: repr(item[0]))),
        )

    @property
    def terms(self) -> dict:
        """Return a copy; modifying it cannot alter this Hamiltonian."""
        return dict(self._items)

    @property
    def momentum_degrees(self) -> tuple[int, ...]:
        return tuple(sorted({engine.tree_degree(t) for t, _ in self._items}))

    @property
    def orders(self) -> tuple[int, ...]:
        """Original-letter orders present in the expression; empty for zero."""
        return tuple(sorted({engine.tree_order(t) for t, _ in self._items}))

    def __bool__(self) -> bool:
        return bool(self._items)

    def __add__(self, other) -> Hamiltonian:
        if not isinstance(other, Hamiltonian):
            return NotImplemented
        return Hamiltonian(engine.add(self.terms, other.terms))

    def __neg__(self) -> Hamiltonian:
        return Hamiltonian({t: -c for t, c in self._items})

    def __sub__(self, other) -> Hamiltonian:
        if not isinstance(other, Hamiltonian):
            return NotImplemented
        return self + -other

    def __mul__(self, scalar: int | Fraction) -> Hamiltonian:
        scalar = exact(scalar)
        return Hamiltonian({t: c * scalar for t, c in self._items})

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | Fraction) -> Hamiltonian:
        scalar = exact(scalar)
        if not scalar:
            raise ZeroDivisionError("Division by zero")
        return self * (1 / scalar)

    def realize(self, potential: Polynomial, *, max_assignments: int = 100000) -> Polynomial:
        """Contract every tree to a complete exact polynomial for V(q).

        Reject momentum-dependent potentials and a tree whose index-assignment
        count exceeds the explicit cap. This cap does not bound polynomial
        term growth or the total runtime across an arbitrary expression.
        """
        if not isinstance(potential, Polynomial):
            raise TypeError("potential must be a Polynomial")
        d = potential.dimension
        if any(any(a[d:]) for a in potential.terms):
            raise ValueError("The potential must depend only on position variables")
        positive_int(max_assignments, "max_assignments")
        result = Polynomial.constant(d, 0)
        _, momenta = phase_space(d)
        derivatives = {}
        for tree, coefficient in self._items:
            if tree is None:
                result += coefficient * Polynomial.kinetic(d)
                continue
            colors, adjacency = engine.decode_mechanical(tree)
            edges = [(i, j) for i, row in enumerate(adjacency) for j in row if i < j]
            free = [i for i, count in enumerate(colors) for _ in range(count)]
            assignments = d ** (len(edges) + len(free))
            if assignments > max_assignments:
                raise ValueError(
                    f"Tree requires {assignments} index assignments; cap is {max_assignments}"
                )
            value = Polynomial.constant(d, 0)
            for labels in product(range(d), repeat=len(edges) + len(free)):
                counts = [[0] * d for _ in colors]
                for label, (i, j) in zip(labels, edges):
                    counts[i][label] += 1
                    counts[j][label] += 1
                monomial = Polynomial.constant(d, 1)
                for label, vertex in zip(labels[len(edges) :], free):
                    counts[vertex][label] += 1
                    monomial *= momenta[label]
                for alpha in map(tuple, counts):
                    if alpha not in derivatives:
                        derivative = potential
                        for i, count in enumerate(alpha):
                            for _ in range(count):
                                derivative = derivative.derivative(i)
                        derivatives[alpha] = derivative
                    monomial *= derivatives[alpha]
                value += monomial
            result += coefficient * value
        return result


def kinetic() -> Hamiltonian:
    """Universal T = sum(p_i**2)/2."""
    return Hamiltonian({None: 1})


def potential() -> Hamiltonian:
    """Universal position potential V(q)."""
    return Hamiltonian({(0, ()): 1})


def bracket(left: Hamiltonian, right: Hamiltonian) -> Hamiltonian:
    """Euclidean mechanical Poisson bracket with sign f_q*g_p - f_p*g_q."""
    if not isinstance(left, Hamiltonian) or not isinstance(right, Hamiltonian):
        raise TypeError("bracket expects two Hamiltonian objects")
    return Hamiltonian(engine.mechanical_bracket(left.terms, right.terms))


def commutative_product(left: Hamiltonian, right: Hamiltonian) -> Hamiltonian:
    """Gradient product left*right = {right,{left,T}}, for degree-zero inputs.

    This product is commutative and nonassociative, not pointwise multiplication.
    """
    if not isinstance(left, Hamiltonian) or not isinstance(right, Hamiltonian):
        raise TypeError("commutative_product expects two Hamiltonian objects")
    if set(left.momentum_degrees + right.momentum_degrees) - {0}:
        raise ValueError("Gradient products require degree-zero potentials")
    return bracket(right, bracket(left, kinetic()))
