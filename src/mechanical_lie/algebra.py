"""Commutative binary trees, associative Lie words and degree lowering.

Advanced word operations use tuples of A=('A',) and Y(tree)=('Y',tree).
Their precise PBW and marker preconditions are documented in docs/api.md.
"""

from collections.abc import Mapping
from ._engine import associative as _words
from ._engine import magmatic as _trees
from ._engine.magmatic import A, B, magmatic
from ._engine.mechanical import theta_X
from .hamiltonian import Hamiltonian
from .polynomial import exact, positive_int


def _check_tree(tree):
    if not isinstance(tree, tuple) or len(tree) not in (0, 2):
        raise ValueError("A binary tree must be () or a pair of binary trees")
    if tree:
        _check_tree(tree[0])
        _check_tree(tree[1])
        if tree[0] > tree[1]:
            raise ValueError("Binary tree pairs must be in canonical sorted order")


def _check_word(word):
    if not isinstance(word, tuple):
        raise TypeError("A word must be a tuple of A or Y(tree) letters")
    for letter in word:
        if letter == A:
            continue
        if not isinstance(letter, tuple) or len(letter) != 2 or letter[0] != "Y":
            raise ValueError("A letter must be A=('A',) or Y(tree)=('Y',tree)")
        _check_tree(letter[1])


def _checked_polynomial(value):
    if not isinstance(value, Mapping):
        raise TypeError("An associative polynomial must map words to exact coefficients")
    for word, coefficient in value.items():
        _check_word(word)
        exact(coefficient)
    return dict(value)


def leaves(t: tuple) -> int:
    """Count leaves of a canonical commutative binary tree."""
    _check_tree(t)
    return _trees.leaves(t)


def star(a: tuple, b: tuple) -> tuple:
    """Form a canonical sorted pair from two canonical binary trees."""
    _check_tree(a)
    _check_tree(b)
    return _trees.star(a, b)


def comb(n: int) -> tuple:
    """Construct the n-leaf comb; n is a positive integer."""
    positive_int(n, "n")
    return _trees.comb(n)


def Y(t: tuple) -> tuple:
    """Construct the potential-derived letter for a canonical binary tree."""
    _check_tree(t)
    return _trees.Y(t)


def lyndon(w: tuple) -> bool:
    """Test the fixed lexicographic Lyndon condition; the empty word is false."""
    _check_word(w)
    return _words.lyndon(w)


def lie_expansion(w: tuple) -> dict:
    """Standard Lyndon bracketing; require a nonempty Lyndon word."""
    if not lyndon(w):
        raise ValueError("A nonempty Lyndon word is required")
    return _words.lie_expansion(w)


def commutator(a: Mapping, b: Mapping) -> dict:
    """Associative a*b-b*a on exact word polynomials."""
    return _words.commutator(_checked_polynomial(a), _checked_polynomial(b))


def adjoint_x(x: tuple, f: Mapping) -> dict:
    """Compute the normal PBW form of [x,f]; no purity promise is inferred."""
    _check_tree(x)
    return _words.adjoint_x(x, _checked_polynomial(f))


def derivative_A(f: Mapping, x: tuple) -> dict:
    """Replace each occurrence of A in turn by Y(x)."""
    _check_tree(x)
    return _words.derivative_A(_checked_polynomial(f), x)


def project_large_labels(f: Mapping, M: int) -> dict:
    """Delete words with a Y label exceeding nonnegative leaf threshold M."""
    if isinstance(M, bool) or not isinstance(M, int) or M < 0:
        raise ValueError("M must be a nonnegative integer leaf threshold")
    return _words.project_large_labels(_checked_polynomial(f), M)


def retract_marker(f: Mapping, x: tuple) -> dict:
    """Replace exactly the letter Y(x) by A."""
    _check_tree(x)
    return _words.retract_marker(_checked_polynomial(f), x)


def binary_trees(leaves_count: int) -> tuple:
    """Enumerate canonical nonplanar commutative binary trees with n leaves."""
    positive_int(leaves_count, "leaves_count")
    return magmatic(leaves_count)


def potential_tree(tree: tuple) -> Hamiltonian:
    """Map a commutative binary tree to its gradient-product potential."""

    _check_tree(tree)
    return Hamiltonian(theta_X(tree))


__all__ = [
    "A",
    "B",
    "Y",
    "adjoint_x",
    "binary_trees",
    "comb",
    "commutator",
    "derivative_A",
    "leaves",
    "lie_expansion",
    "lyndon",
    "potential_tree",
    "project_large_labels",
    "retract_marker",
    "star",
]
