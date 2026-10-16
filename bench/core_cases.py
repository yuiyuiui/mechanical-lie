"""Dense exact polynomial brackets with independent symbolic differentiation."""

import hashlib
import sympy as sp
from mechanical_lie import phase_space, poisson


def make(case):
    dimension, a_degree, b_degree = {"small": (1, 3, 3), "medium": (2, 4, 3), "large": (2, 5, 4)}[
        case
    ]
    q, p = phase_space(dimension)
    f = (sum(q) + sum(p) + 1) ** a_degree
    g = (sum(q) - sum(p) + 2) ** b_degree
    variables = sp.symbols(f"x:{2 * dimension}")

    def symbolic(poly):
        return sum(
            sp.Rational(c.numerator, c.denominator) * sp.prod(x**n for x, n in zip(variables, a))
            for a, c in poly.terms.items()
        )

    sf, sg = symbolic(f), symbolic(g)
    expected = sp.Poly(
        sp.expand(
            sum(
                sp.diff(sf, variables[i]) * sp.diff(sg, variables[i + dimension])
                - sp.diff(sf, variables[i + dimension]) * sp.diff(sg, variables[i])
                for i in range(dimension)
            )
        ),
        *variables,
    )
    reference = {a: (int(c.p), int(c.q)) for a, c in expected.terms()}

    def call():
        return poisson(f, g)

    def validate(result):
        actual = {a: (c.numerator, c.denominator) for a, c in result.terms.items()}
        if actual != reference:
            raise ArithmeticError("Full exact SymPy coefficient comparison failed")
        digest = hashlib.sha256(repr(sorted(actual.items())).encode()).hexdigest()
        return {
            "passed": True,
            "full_polynomial_reference": "SymPy derivatives",
            "terms": len(actual),
            "output_sha256": digest,
        }

    return (
        call,
        validate,
        {
            "dimension": dimension,
            "left_terms": len(f.terms),
            "right_terms": len(g.terms),
            "arithmetic": "exact rational",
        },
    )
