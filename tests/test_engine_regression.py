"""Independent identities, reference cases and adversarial certificate tests."""

from fractions import Fraction
import itertools
from pathlib import Path
import tempfile
import unittest

from mechanical_lie._engine import ValidationError
from mechanical_lie._engine import core as e
from mechanical_lie._engine.audit import run_audit
from mechanical_lie._engine.certificates import digest, verify_block
from mechanical_lie._engine.polynomial import poisson, polynomial_calibration, iterated_lowering


class MagmaticTests(unittest.TestCase):
    def test_independent_tree_counts(self):
        counts = [0, 1]
        for n in range(2, 10):
            total = sum(counts[i] * counts[n - i] for i in range(1, (n + 1) // 2))
            if n % 2 == 0:
                k = counts[n // 2]
                total += k * (k + 1) // 2
            counts.append(total)
        self.assertEqual(counts[1:], [1, 1, 1, 2, 3, 6, 11, 23, 46])
        self.assertEqual([len(e.magmatic(n)) for n in range(1, 10)], counts[1:])

    def test_commutativity_not_associativity(self):
        b, s = e.B, e.star(e.B, e.B)
        self.assertEqual(e.star(b, s), e.star(s, b))
        self.assertNotEqual(e.star(e.star(b, b), s), e.star(b, e.star(b, s)))

    def test_marker_word_separation_and_empty_word(self):
        marker = e.comb(5)
        for length in range(4):
            for word in itertools.product((e.Y(e.B), e.Y(e.star(e.B, e.B))), repeat=length):
                tree = e.left_action_word(word, marker)
                self.assertEqual(e.decode_marker(tree, 5), word)

    def test_ambiguous_marker_rejected(self):
        with self.assertRaises(ValidationError):
            e.decode_marker(e.star(e.comb(3), e.comb(3)), 3)

    def test_invalid_comb(self):
        with self.assertRaises(ValueError):
            e.comb(0)


class AssociativeTests(unittest.TestCase):
    def test_lie_sign(self):
        a, b = e.A, e.Y(e.B)
        self.assertEqual(e.lie_expansion((a, b)), {(a, b): 1, (b, a): -1})

    def test_associative_jacobi(self):
        a, b, c = ({(e.A,): 1}, {(e.Y(e.B),): 1}, {(e.Y(e.star(e.B, e.B)),): 1})
        value = e.add(
            e.add(e.commutator(a, e.commutator(b, c)), e.commutator(b, e.commutator(c, a))),
            e.commutator(c, e.commutator(a, b)),
        )
        self.assertEqual(value, {})

    def test_pure_A_euler_retraction(self):
        marker = e.comb(4)
        f = {(e.A, e.A, e.Y(e.B)): 3}
        self.assertEqual(
            e.retract_marker(e.derivative_A(f, marker), marker), {w: 2 * c for w, c in f.items()}
        )

    def test_lowering_A(self):
        self.assertEqual(e.adjoint_x(e.B, {(e.A,): 1}), {((e.Y(e.B),), None): 1})

    def test_lowering_Y_sign(self):
        self.assertEqual(e.adjoint_x(e.B, {(e.Y(e.B),): 1}), {((), e.star(e.B, e.B)): 1})

    def test_projection_then_retraction(self):
        marker = e.comb(4)
        w = (e.A, e.A, e.Y(e.B))
        f = e.lie_expansion(w)
        ax = e.adjoint_x(marker, f)
        self.assertTrue(all(z is None for _, z in ax))
        projected = e.project_large_labels({v: c for (v, _), c in ax.items()}, 4)
        self.assertEqual(projected, e.derivative_A(f, marker))
        self.assertEqual(e.retract_marker(projected, marker), {v: 2 * c for v, c in f.items()})

    def test_invalid_lyndon_split(self):
        with self.assertRaises(ValueError):
            e.standard_split((e.A, e.A))

    def test_cached_outputs_cannot_be_poisoned(self):
        w = (e.A, e.Y(e.B))
        before = e.lie_expansion(w)
        e.lie_expansion(w).clear()
        self.assertEqual(e.lie_expansion(w), before)


class MechanicalTests(unittest.TestCase):
    def test_poisson_generator_sign(self):
        self.assertEqual(e.mechanical_bracket({(0, ()): 1}, {None: 1}), {(1, ()): 1})

    def test_potential_brackets_zero(self):
        for a in e.magmatic(3):
            self.assertEqual(e.mechanical_bracket(e.theta_X(a), e.theta_X(e.B)), {})

    def test_gradient_square(self):
        self.assertEqual(e.theta_X(e.star(e.B, e.B)), {(0, ((0, ()),)): 1})

    def test_canonical_relabelling(self):
        colors = (0, 1, 2, 0)
        adj = ((1,), (0, 2), (1, 3), (2,))
        reference = e.canonical(colors, adj)
        for perm in itertools.permutations(range(4)):
            pos = {v: i for i, v in enumerate(perm)}
            self.assertEqual(
                e.canonical(
                    tuple(colors[v] for v in perm),
                    tuple(tuple(pos[j] for j in adj[v]) for v in perm),
                ),
                reference,
            )

    def test_invalid_graph(self):
        for colors, adj in [
            ((0,), ((0,),)),
            ((0, 0), ((), ())),
            ((0, 0, 0), ((1, 2), (0, 2), (0, 1))),
        ]:
            with self.assertRaises(ValueError):
                e.canonical(colors, adj)

    def test_jacobi_nontrivial(self):
        a, b, y = {None: 1}, {(0, ()): 1}, {(1, ()): 1}
        value = e.add(
            e.add(
                e.mechanical_bracket(a, e.mechanical_bracket(b, y)),
                e.mechanical_bracket(b, e.mechanical_bracket(y, a)),
            ),
            e.mechanical_bracket(y, e.mechanical_bracket(a, b)),
        )
        self.assertEqual(value, {})


class PolynomialTests(unittest.TestCase):
    def test_fractional_poisson_hand_calculation(self):
        # In one dimension {q^3/3,p^2/2}=q^2 p.
        v = {(3, 0): Fraction(1, 3)}
        t = {(0, 2): Fraction(1, 2)}
        self.assertEqual(poisson(v, t, 1), {(2, 1): Fraction(1)})
        self.assertEqual(poisson(v, poisson(v, t, 1), 1), {(4, 0): Fraction(1)})

    def test_two_dimensional_sign(self):
        self.assertEqual(poisson({(0, 1, 0, 0): 1}, {(0, 0, 0, 1): 1}, 2), {(0, 0, 0, 0): 1})

    def test_direct_polynomial_crosscheck(self):
        self.assertEqual(polynomial_calibration(6, 2, 519)["checks"], 16)

    def test_iterated_lowering(self):
        traces = iterated_lowering(4)
        self.assertEqual(len(traces), 5)
        self.assertTrue(all(t["nonzero_final_potential_terms"] > 0 for t in traces))


class CertificateTests(unittest.TestCase):
    def block(self):
        columns = [{"a": 2, "b": 1}, {"a": 1, "b": 1}]
        rank, pivots = e.sparse_rank(columns, 101)
        data = {
            "order": 1,
            "degree": 0,
            "prime": 101,
            "dimension": 2,
            "rank": rank,
            "columns": [{repr(row): c for row, c in col.items()} for col in columns],
            "pivots": pivots,
            "source_basis": ["u", "v"],
        }
        data["sha256_without_hash"] = digest(data)
        return data

    def test_replay(self):
        self.assertEqual(verify_block(self.block()), 2)

    def test_corrupt_data(self):
        data = self.block()
        data["columns"][0]["a"] = 3
        with self.assertRaises(ValidationError):
            verify_block(data)

    def test_corrupt_pivot_with_new_hash(self):
        data = self.block()
        data["pivots"][0]["value"] = 3
        data["sha256_without_hash"] = digest(data)
        with self.assertRaises(ValidationError):
            verify_block(data)

    def test_nonprime_modulus(self):
        with self.assertRaises(ValueError):
            e.sparse_rank([{"a": 1}], 15)
        data = self.block()
        data["prime"] = 15
        data["sha256_without_hash"] = digest(data)
        with self.assertRaises(ValidationError):
            verify_block(data)

    def test_dependent_matrix(self):
        self.assertEqual(e.sparse_rank([{"a": 1}, {"a": 2}], 101)[0], 1)

    def test_finite_audit_known_dimensions(self):
        with tempfile.TemporaryDirectory() as temp:
            data = run_audit(6, Path(temp) / "data")
            self.assertEqual(data["total_dimension"], 16)
            self.assertTrue(data["all_checks_passed"])

    def test_output_collision_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileExistsError):
                run_audit(2, Path(temp))


if __name__ == "__main__":
    unittest.main()
