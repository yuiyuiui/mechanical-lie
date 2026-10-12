"""Public contracts, hand calculations and independent exact references."""

from fractions import Fraction
import tempfile
from pathlib import Path
import unittest
import sympy as sp

from mechanical_lie import (
    Hamiltonian,
    Polynomial,
    ValidationError,
    binary_trees,
    bracket,
    commutative_product,
    kinetic,
    modular_rank,
    phase_space,
    poisson,
    potential,
    potential_tree,
)
from mechanical_lie._engine.audit import run_audit
from mechanical_lie.validation import verify_certificate_block


def sympy_expression(poly):
    variables = sp.symbols(
        " ".join(
            [f"q{i + 1}" for i in range(poly.dimension)]
            + [f"p{i + 1}" for i in range(poly.dimension)]
        )
    )
    return sp.expand(
        sum(
            sp.Rational(c.numerator, c.denominator)
            * sp.prod(x**a for x, a in zip(variables, powers))
            for powers, c in poly.terms.items()
        )
    )


class PolynomialContracts(unittest.TestCase):
    def test_hand_calculation(self):
        q, p = phase_space(1)
        v, t = q[0] ** 3 / 3, p[0] ** 2 / 2
        self.assertEqual(poisson(v, t), q[0] ** 2 * p[0])
        self.assertEqual(poisson(v, poisson(v, t)), q[0] ** 4)

    def test_exact_arithmetic_and_evaluation(self):
        q, p = phase_space(1)
        f = (q[0] + p[0] / 3) ** 2
        self.assertEqual(f.evaluate((Fraction(1, 2), 3)), Fraction(9, 4))
        self.assertEqual(f.derivative(0), 2 * q[0] + 2 * p[0] / 3)

    def test_dimensions_and_exponents(self):
        for d in (0, -1):
            with self.assertRaises(ValueError):
                Polynomial(d, {})
        with self.assertRaises(TypeError):
            Polynomial(True, {})
        for a in ((1,), (1, 0, 0), (True, 0), (-1, 0), (1.0, 0)):
            with self.assertRaises(ValueError):
                Polynomial(1, {a: 1})
        with self.assertRaises(ValueError):
            poisson(Polynomial(1, {}), Polynomial(2, {}))
        with self.assertRaises(ValueError):
            Polynomial(1, {}).evaluate((1,))
        with self.assertRaises(ValueError):
            Polynomial(1, {}).evaluate((1, 2, 3))

    def test_reject_float_and_boolean_coefficients(self):
        for c in (0.5, True, "1/2"):
            with self.assertRaises(TypeError):
                Polynomial(1, {(0, 0): c})
        q, _ = phase_space(1)
        with self.assertRaises(TypeError):
            q[0].evaluate((1.0, 0))
        with self.assertRaises(TypeError):
            q[0] + 0.5

    def test_power_and_division_contracts(self):
        q, _ = phase_space(1)
        self.assertEqual(q[0] ** 0, Polynomial.constant(1, 1))
        with self.assertRaises(ZeroDivisionError):
            q[0] / 0
        for power in (-1, 0.5, True):
            with self.assertRaises(ValueError):
                q[0] ** power
        with self.assertRaises(ValueError):
            q[0].derivative(2)

    def test_immutable_terms_and_hash(self):
        q, _ = phase_space(1)
        before = hash(q[0])
        q[0].terms.clear()
        self.assertEqual(hash(q[0]), before)
        self.assertEqual(str(q[0]), "q1")

    def test_sympy_poisson_nonseparable_reference(self):
        q, p = phase_space(2)
        f = (q[0] + q[1]) ** 3 / 3 + p[0] * q[1]
        g = p[0] ** 2 / 2 + p[1] ** 3 / 3 + q[0] * p[1]
        sf, sg = sympy_expression(f), sympy_expression(g)
        qs = sp.symbols("q1 q2")
        ps = sp.symbols("p1 p2")
        reference = sum(
            sp.diff(sf, x) * sp.diff(sg, y) - sp.diff(sf, y) * sp.diff(sg, x)
            for x, y in zip(qs, ps)
        )
        self.assertEqual(sp.expand(reference - sympy_expression(poisson(f, g))), 0)


class HamiltonianContracts(unittest.TestCase):
    def test_quartic_full_polynomials(self):
        q, p = phase_space(1)
        v = q[0] ** 2 / 2 + q[0] ** 4 / 4
        a, b = kinetic(), potential()
        self.assertEqual(bracket(b, a).realize(v), (q[0] + q[0] ** 3) * p[0])
        self.assertEqual(bracket(b, bracket(b, a)).realize(v), (q[0] + q[0] ** 3) ** 2)
        self.assertEqual(bracket(a, bracket(a, b)).realize(v), (1 + 3 * q[0] ** 2) * p[0] ** 2)

    def test_sympy_tree_contraction_reference(self):
        q, p = phase_space(2)
        v = q[0] ** 2 / 2 + q[0] * q[1] + q[1] ** 4 / 4
        a, b = kinetic(), potential()
        qs = sp.symbols("q1 q2")
        ps = sp.symbols("p1 p2")
        sv = sympy_expression(v)
        st = sum(x**2 for x in ps) / 2

        def pb(f, g):
            return sp.expand(
                sum(
                    sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x)
                    for x, y in zip(qs, ps)
                )
            )

        cases = [(a, st), (b, sv)]
        for i, j in [(1, 0), (0, 1), (1, 2), (0, 2), (0, 4), (1, 5)]:
            cases.append((bracket(cases[i][0], cases[j][0]), pb(cases[i][1], cases[j][1])))
        for h, reference in cases:
            self.assertEqual(sp.expand(sympy_expression(h.realize(v)) - reference), 0)

    def test_jacobi_with_rational_coefficients(self):
        a, b = kinetic() / 3, potential() * Fraction(2, 5)
        c = bracket(b, a)
        self.assertFalse(
            bracket(a, bracket(b, c)) + bracket(b, bracket(c, a)) + bracket(c, bracket(a, b))
        )

    def test_canonicalization_and_cancellation(self):
        x = Hamiltonian({(1, ((0, ()),)): 1})
        y = Hamiltonian({(0, ((1, ()),)): 1})
        self.assertEqual(x, y)
        self.assertFalse(x - y)
        self.assertEqual((x - y).orders, ())
        x.terms.clear()
        self.assertTrue(x)

    def test_grading_and_gradient_product(self):
        a, b = kinetic(), potential()
        self.assertEqual(bracket(b, a).momentum_degrees, (1,))
        self.assertEqual(bracket(b, a).orders, (2,))
        self.assertEqual(commutative_product(b, b), bracket(b, bracket(b, a)))
        self.assertEqual(potential_tree(binary_trees(2)[0]), commutative_product(b, b))
        with self.assertRaises(ValueError):
            commutative_product(a, b)

    def test_tree_errors(self):
        for tree in ((-1, ()), (True, ()), (0, None), (0,), (0, ((0,),))):
            with self.assertRaises(ValueError):
                Hamiltonian({tree: 1})
        with self.assertRaises(ValueError):
            potential_tree(((),))
        with self.assertRaises(ValueError):
            binary_trees(0)
        with self.assertRaises(TypeError):
            binary_trees(True)

    def test_realization_preconditions_and_guard(self):
        q, p = phase_space(2)
        with self.assertRaises(ValueError):
            potential().realize(q[0] + p[0])
        with self.assertRaises(ValueError):
            Hamiltonian({(3, ()): 1}).realize(q[0] ** 4, max_assignments=1)
        with self.assertRaises(ValueError):
            potential().realize(q[0], max_assignments=0)
        with self.assertRaises(TypeError):
            bracket(potential(), q[0])
        with self.assertRaises(TypeError):
            kinetic() * 0.5


class ValidationContracts(unittest.TestCase):
    def test_fixed_case_manifest_without_git(self):
        import contextlib
        import importlib.util
        import io
        import json
        from unittest.mock import patch

        script = Path(__file__).resolve().parents[1] / "numerical_case/compute_cases.py"
        spec = importlib.util.spec_from_file_location("fixed_cases", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with (
            tempfile.TemporaryDirectory() as temp,
            patch.object(module.shutil, "which", return_value=None),
        ):
            output = Path(temp) / "case"
            with contextlib.redirect_stdout(io.StringIO()):
                module.compute(output)
            record = json.loads((output / "manifest.json").read_text())
            self.assertIsNone(record["commit_observed"])
            self.assertEqual(record["actual_parameters"]["reconstruction_max_order"], 10)
            self.assertEqual(record["actual_parameters"]["coefficient_checks_max_order"], 6)
            self.assertEqual(record["actual_parameters"]["potential_degree"], 4)

    def test_modular_rank_hand_cases(self):
        self.assertEqual(modular_rank([{"x": 2}], 2), 0)
        self.assertEqual(modular_rank([{"x": 2}], 3), 1)
        self.assertEqual(modular_rank([{"x": 2, "y": 1}, {"x": 1, "y": 1}], 101), 2)
        self.assertEqual(modular_rank([], 101), 0)
        with self.assertRaises(ValueError):
            modular_rank([{"x": 1}], 15)
        with self.assertRaises(TypeError):
            modular_rank([{"x": 0.5}], 101)
        with self.assertRaises(TypeError):
            modular_rank([{1: 1}], 101)

    def test_missing_certificate_fields(self):
        with self.assertRaises(ValidationError):
            verify_certificate_block({})

    def test_partial_status_is_explicit(self):
        with tempfile.TemporaryDirectory() as temp:
            result = run_audit(3, Path(temp) / "partial", skip_tree_rank=True)
            self.assertIsNone(result["all_checks_passed"])
            self.assertEqual(result["validation_status"], "PARTIAL")
            self.assertTrue(result["skipped_stages"])


if __name__ == "__main__":
    unittest.main()
