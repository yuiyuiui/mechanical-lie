"""Output and failure contracts for defaults, explicit choices and automatic plans."""

import contextlib
from dataclasses import FrozenInstanceError
from fractions import Fraction
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import sympy as sp

from mechanical_lie import (
    Hamiltonian,
    Polynomial,
    bracket,
    kinetic,
    modular_rank,
    modular_rank_report,
    phase_space,
    potential,
)
from mechanical_lie import algebra as a
from mechanical_lie.validation import rebuild_basis, reconstruction_plan


class RealizationPolicies(unittest.TestCase):
    def setUp(self):
        q, _ = phase_space(2)
        self.v = q[0] ** 2 / 2 + q[0] * q[1] + q[1] ** 4 / 4
        self.h = bracket(potential(), bracket(potential(), kinetic()))

    def test_default_and_exact_boundary_against_independent_gradient(self):
        x, y = sp.symbols("x y")
        reference = sp.Poly(
            sp.diff(x * x / 2 + x * y + y**4 / 4, x) ** 2
            + sp.diff(x * x / 2 + x * y + y**4 / 4, y) ** 2,
            x,
            y,
        )
        expected = Polynomial(
            2, {powers + (0, 0): Fraction(int(c.p), int(c.q)) for powers, c in reference.terms()}
        )
        plan = self.h.realization_plan(self.v)
        self.assertEqual(plan.tree_assignments, (2,))
        self.assertEqual(self.h.realize(self.v), expected)
        self.assertEqual(self.h.realize(self.v, max_assignments=2), expected)
        self.assertEqual(self.h.realize(self.v, max_assignments=200000), expected)
        with self.assertRaisesRegex(ValueError, "requires 2.*cap is 1"):
            self.h.realize(self.v, max_assignments=1)

    def test_input_dependent_assignment_counts(self):
        expression = Hamiltonian({(2, ((0, ()),)): 1})  # E=1, F=2.
        for dimension in (1, 2, 3):
            q, _ = phase_space(dimension)
            plan = expression.realization_plan(q[0] ** 4, max_assignments=8)
            self.assertEqual(plan.total_assignments, dimension**3)
            self.assertEqual(plan.within_limit, dimension <= 2)

    def test_plan_is_pure_deterministic_and_immutable(self):
        with patch("mechanical_lie.hamiltonian.phase_space", side_effect=AssertionError):
            first = self.h.realization_plan(self.v)
            second = self.h.realization_plan(self.v)
        self.assertEqual(first, second)
        self.assertIn("engineering limit", first.selection_reason)
        self.assertEqual(json.loads(json.dumps(first.to_dict()))["largest_assignments"], 2)
        with self.assertRaises(FrozenInstanceError):
            first.max_assignments = 1

    def test_all_guards_before_any_polynomial_contraction(self):
        expression = potential() + Hamiltonian({(4, ()): 1})
        with patch("mechanical_lie.hamiltonian.phase_space", side_effect=AssertionError):
            with self.assertRaisesRegex(ValueError, "requires 16"):
                expression.realize(self.v, max_assignments=1)

    def test_zero_kinetic_and_constant_potentials(self):
        zero = Hamiltonian({})
        self.assertEqual(zero.realization_plan(self.v).total_assignments, 0)
        self.assertEqual(kinetic().realization_plan(self.v).tree_assignments, (0,))
        self.assertEqual(kinetic().realize(self.v, max_assignments=1), Polynomial.kinetic(2))
        self.assertFalse(self.h.realize(Polynomial.constant(2, 3)))

    def test_invalid_caps_and_unavailable_backend(self):
        for value in (True, 1.5, "auto", None):
            with self.assertRaises(TypeError):
                self.h.realize(self.v, max_assignments=value)
        for value in (0, -1):
            with self.assertRaises(ValueError):
                self.h.realization_plan(self.v, max_assignments=value)
        with self.assertRaises(TypeError):
            self.h.realize(self.v, backend="float")


class RankPolicies(unittest.TestCase):
    def test_legacy_integer_subclass_tuning_remains_supported(self):
        from enum import IntEnum

        class Choice(IntEnum):
            SMALL = 3

        q, _ = phase_space(1)
        self.assertEqual(potential().realize(q[0], max_assignments=Choice.SMALL), q[0])
        self.assertEqual(modular_rank([{"x": 1}], prime=Choice.SMALL), 1)
        self.assertEqual(reconstruction_plan(Choice.SMALL).max_order, 3)

    def test_default_override_and_rational_reference(self):
        columns = [{"x": 2, "y": 1}, {"x": 1, "y": 1}]
        reference = int(sp.Matrix([[2, 1], [1, 1]]).rank())
        for kwargs in ({}, {"prime": 2}, {"prime": 101}):
            report = modular_rank_report(iter(columns), **kwargs)
            self.assertEqual(report.rank, reference)
            self.assertEqual(modular_rank(columns, **kwargs), reference)
            self.assertTrue(report.rational_independence_established)
            self.assertEqual((report.column_count, report.row_count), (2, 2))

    def test_bad_default_prime_is_not_rational_dependence(self):
        columns = [{"x": 1000003}]
        report = modular_rank_report(columns)
        self.assertEqual(report.rank, 0)
        self.assertFalse(report.rational_independence_established)
        self.assertEqual(sp.Matrix([[1000003]]).rank(), 1)
        self.assertEqual(modular_rank_report(columns, prime=101).rank, 1)

    def test_zero_rows_zero_columns_and_dependent_columns(self):
        empty = modular_rank_report([])
        self.assertEqual((empty.rank, empty.row_count, empty.column_count), (0, 0, 0))
        self.assertTrue(empty.rational_independence_established)  # Empty family.
        zero = modular_rank_report([{"zero": 0}, {}])
        self.assertEqual((zero.rank, zero.row_count, zero.column_count), (0, 0, 2))
        self.assertFalse(
            modular_rank_report([{"x": 1}, {"x": 2}]).rational_independence_established
        )

    def test_reject_prime_before_consuming_columns(self):
        def unavailable():
            raise AssertionError("invalid policy must not consume data")
            yield {}

        for value in (0, 1, 15):
            with self.assertRaises(ValueError):
                modular_rank(unavailable(), prime=value)
        for value in (True, 101.0, "auto", None):
            with self.assertRaises(TypeError):
                modular_rank(unavailable(), prime=value)

    def test_rank_report_immutable_and_reproducible(self):
        first = modular_rank_report([{"x": 1}])
        self.assertEqual(first, modular_rank_report([{"x": 1}]))
        self.assertIn("no matrix-dependent", first.selection_reason)
        self.assertEqual(first.to_dict()["prime"], first.prime)
        with self.assertRaises(FrozenInstanceError):
            first.rank = 0


class ReconstructionPolicies(unittest.TestCase):
    def test_automatic_marker_on_even_odd_orders_and_upper_boundary(self):
        for order, marker in ((2, 2), (3, 2), (6, 4), (7, 4), (10, 6), (18, 10)):
            plan = reconstruction_plan(order)
            self.assertEqual((plan.max_order, plan.marker_leaves), (order, marker))
            self.assertGreater(2 * marker, order)
            self.assertEqual(plan.verification, "all-stages")
        self.assertEqual(reconstruction_plan(6), reconstruction_plan(6))

    def test_explicit_prime_preserves_target_and_saved_configuration(self):
        with tempfile.TemporaryDirectory() as temp:
            for prime in (101, 1000003):
                output = Path(temp) / str(prime)
                with contextlib.redirect_stdout(io.StringIO()):
                    result = rebuild_basis(3, output, prime=prime)
                self.assertEqual(result["total_dimension"], 5)
                self.assertEqual(result["validation_status"], "COMPLETE")
                self.assertEqual(result["skipped_stages"], [])
                plan = reconstruction_plan(3, prime=prime).to_dict()
                self.assertEqual(result["effective_configuration"], plan)
                self.assertEqual(json.loads((output / "audit_summary.json").read_text()), result)

    def test_invalid_targets_do_not_create_output(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "new"
            for value in (1, 19, 0):
                with self.assertRaises(ValueError):
                    rebuild_basis(value, output)
            for value in (True, 3.0, "auto", None):
                with self.assertRaises(TypeError):
                    rebuild_basis(value, output)
            with self.assertRaises(ValueError):
                rebuild_basis(3, output, prime=15)
            self.assertFalse(output.exists())


class AdvancedParameterContracts(unittest.TestCase):
    def test_exact_word_scalars_and_marker_identity(self):
        word = (a.A, a.Y(a.B))
        f = {word: Fraction(1, 3)}
        marker = a.comb(4)
        self.assertEqual(a.retract_marker(a.derivative_A(f, marker), marker), f)
        self.assertEqual(a.commutator({(): Fraction(1, 2)}, f), {})
        self.assertEqual(a.leaves(a.star(a.B, marker)), 5)

    def test_projection_zero_and_nonzero_thresholds(self):
        f = {(a.A,): 1, (a.Y(a.B),): 2, (a.Y(a.comb(2)),): 3}
        self.assertEqual(a.project_large_labels(f, 0), {(a.A,): 1})
        self.assertEqual(a.project_large_labels(f, 1), {(a.A,): 1, (a.Y(a.B),): 2})
        self.assertEqual(a.project_large_labels(f, 2), f)

    def test_invalid_trees_words_and_thresholds(self):
        for invalid in (((),), (1, ()), [(), ()]):
            for function in (a.Y, a.leaves, a.potential_tree):
                with self.assertRaises(ValueError):
                    function(invalid)
        for invalid in ((a.A, a.A), (), (("unknown",),)):
            with self.assertRaises(ValueError):
                a.lie_expansion(invalid)
        with self.assertRaises(TypeError):
            a.lyndon([a.A])
        with self.assertRaises(TypeError):
            a.comb(True)
        for invalid in (-1, 1.0, True):
            with self.assertRaises(ValueError):
                a.project_large_labels({}, invalid)
        with self.assertRaises(TypeError):
            a.commutator({(a.A,): 0.5}, {})


if __name__ == "__main__":
    unittest.main()
