"""Immutable sparse rational polynomials in ordered position/momentum variables."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction

from ._engine.polynomial import pdiff, peval, pmul, poisson as _poisson


def exact(value: int | Fraction) -> Fraction:
    """Accept explicit exact coefficients, never silently convert floats."""
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("An int or fractions.Fraction is required; floats are not exact inputs")
    return Fraction(value)


def positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 1:
        raise ValueError(f"{name} must be positive")
    return value


@dataclass(frozen=True, init=False)
class Polynomial:
    """Exact polynomial with exponents ``(q1,...,qd,p1,...,pd)``.

    ``terms`` maps exponent tuples to integers or Fractions. Zero coefficients
    are removed. Arithmetic preserves dimension and returns new objects.
    """

    dimension: int
    _items: tuple

    def __init__(self, dimension: int, terms: Mapping[tuple[int, ...], int | Fraction]):
        positive_int(dimension, "dimension")
        if not isinstance(terms, Mapping):
            raise TypeError("terms must be a mapping of exponent tuples to exact coefficients")
        clean = {}
        for powers, coefficient in terms.items():
            if not isinstance(powers, tuple) or len(powers) != 2 * dimension:
                raise ValueError("Each exponent tuple must have exactly 2 * dimension entries")
            if any(type(n) is not int or n < 0 for n in powers):
                raise ValueError("Exponents must be nonnegative integers")
            coefficient = exact(coefficient)
            if coefficient:
                clean[powers] = coefficient
        object.__setattr__(self, "dimension", dimension)
        object.__setattr__(self, "_items", tuple(sorted(clean.items())))

    @property
    def terms(self) -> dict[tuple[int, ...], Fraction]:
        """Return a fresh sparse coefficient dictionary."""
        return dict(self._items)

    @classmethod
    def constant(cls, dimension: int, value: int | Fraction) -> Polynomial:
        positive_int(dimension, "dimension")
        return cls(dimension, {(0,) * (2 * dimension): exact(value)})

    @classmethod
    def variable(cls, dimension: int, index: int) -> Polynomial:
        positive_int(dimension, "dimension")
        if type(index) is not int or not 0 <= index < 2 * dimension:
            raise ValueError("Variable index must be in range(2 * dimension)")
        powers = [0] * (2 * dimension)
        powers[index] = 1
        return cls(dimension, {tuple(powers): 1})

    @classmethod
    def kinetic(cls, dimension: int) -> Polynomial:
        """Return the Euclidean kinetic energy sum(p_i**2)/2."""
        _, momenta = phase_space(dimension)
        return sum((p**2 / 2 for p in momenta), cls.constant(dimension, 0))

    def _coerce(self, other) -> Polynomial:
        if not isinstance(other, Polynomial):
            other = Polynomial.constant(self.dimension, exact(other))
        if self.dimension != other.dimension:
            raise ValueError("Polynomial dimensions must agree")
        return other

    def __bool__(self) -> bool:
        return bool(self._items)

    def __add__(self, other) -> Polynomial:
        other = self._coerce(other)
        terms = self.terms
        for powers, c in other._items:
            terms[powers] = terms.get(powers, 0) + c
        return Polynomial(self.dimension, terms)

    __radd__ = __add__

    def __neg__(self) -> Polynomial:
        return Polynomial(self.dimension, {a: -c for a, c in self._items})

    def __sub__(self, other) -> Polynomial:
        return self + -self._coerce(other)

    def __rsub__(self, other) -> Polynomial:
        return self._coerce(other) - self

    def __mul__(self, other) -> Polynomial:
        other = self._coerce(other)
        return Polynomial(self.dimension, pmul(self.terms, other.terms))

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | Fraction) -> Polynomial:
        scalar = exact(scalar)
        if not scalar:
            raise ZeroDivisionError("Division by zero")
        return self * (1 / scalar)

    def __pow__(self, exponent: int) -> Polynomial:
        if type(exponent) is not int or exponent < 0:
            raise ValueError("Polynomial powers must be nonnegative integers")
        value = Polynomial.constant(self.dimension, 1)
        base = self
        while exponent:
            if exponent & 1:
                value = value * base
            exponent >>= 1
            if exponent:
                base = base * base
        return value

    def derivative(self, index: int) -> Polynomial:
        if type(index) is not int or not 0 <= index < 2 * self.dimension:
            raise ValueError("Derivative index must be in range(2 * dimension)")
        return Polynomial(self.dimension, pdiff(self.terms, index))

    def evaluate(self, point: Sequence[int | Fraction]) -> Fraction:
        if len(point) != 2 * self.dimension:
            raise ValueError("Evaluation point must have exactly 2 * dimension entries")
        return Fraction(peval(self.terms, tuple(exact(x) for x in point)))

    def __str__(self) -> str:
        if not self:
            return "0"
        names = [f"q{i + 1}" for i in range(self.dimension)] + [
            f"p{i + 1}" for i in range(self.dimension)
        ]
        result = []
        for powers, coefficient in sorted(self._items, reverse=True):
            factors = [name if n == 1 else f"{name}^{n}" for name, n in zip(names, powers) if n]
            monomial = "*".join(factors)
            magnitude = abs(coefficient)
            term = monomial if magnitude == 1 and monomial else str(magnitude)
            if monomial and magnitude != 1:
                term += "*" + monomial
            result.append(("-" if coefficient < 0 else "+", term))
        first_sign, first = result[0]
        return (
            ("-" if first_sign == "-" else "")
            + first
            + "".join(f" {sign} {term}" for sign, term in result[1:])
        )


def phase_space(dimension: int) -> tuple[tuple[Polynomial, ...], tuple[Polynomial, ...]]:
    """Return position and momentum variables, each as a tuple of length d."""
    positive_int(dimension, "dimension")
    variables = tuple(Polynomial.variable(dimension, i) for i in range(2 * dimension))
    return variables[:dimension], variables[dimension:]


def poisson(left: Polynomial, right: Polynomial) -> Polynomial:
    """Compute sum(left_q * right_p - left_p * right_q) exactly."""
    if not isinstance(left, Polynomial) or not isinstance(right, Polynomial):
        raise TypeError("poisson expects two Polynomial objects")
    if left.dimension != right.dimension:
        raise ValueError("Polynomial dimensions must agree")
    return Polynomial(left.dimension, _poisson(left.terms, right.terms, left.dimension))
