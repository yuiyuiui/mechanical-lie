"""Full coefficient comparisons with independent symbolic differentiation."""

from fractions import Fraction
import unittest

import sympy as sp

from mechanical_lie import Polynomial, phase_space, poisson


class OptimizedPoissonTests(unittest.TestCase):
    def test_dense_integer_and_fractional_polynomials_against_sympy(self):
        for dimension in (1, 2, 3):
            variables = sp.symbols(f"x:{2 * dimension}")
            q, p = phase_space(dimension)
            f = (sum(q) + sum(p) + 1) ** 3 / 3
            g = (sum((i + 1) * x for i, x in enumerate(q + p)) - 1) ** 2 / 7

            def expression(poly):
                return sum(
                    sp.Rational(c.numerator, c.denominator)
                    * sp.prod(x**n for x, n in zip(variables, powers))
                    for powers, c in poly.terms.items()
                )

            a, b = expression(f), expression(g)
            expected = sum(
                sp.diff(a, variables[i]) * sp.diff(b, variables[i + dimension])
                - sp.diff(a, variables[i + dimension]) * sp.diff(b, variables[i])
                for i in range(dimension)
            )
            self.assertEqual(sp.expand(expression(poisson(f, g)) - expected), 0)
            self.assertEqual(poisson(f, g), -poisson(g, f))
            self.assertFalse(poisson(f, f))

    def test_constants_empty_and_cancellation_remain_exact(self):
        q, p = phase_space(1)
        self.assertFalse(poisson(q[0], Polynomial.constant(1, Fraction(2, 3))))
        self.assertEqual(poisson(q[0] * p[0], q[0] * p[0] ** 2), q[0] * p[0] ** 2)
