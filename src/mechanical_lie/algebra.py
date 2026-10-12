"""Commutative binary trees, associative Lie words and degree lowering.

Advanced word operations use tuples of A=('A',) and Y(tree)=('Y',tree).
Their precise PBW and marker preconditions are documented in docs/api.md.
"""

from ._engine.associative import (
    adjoint_x,
    commutator,
    derivative_A,
    lie_expansion,
    lyndon,
    project_large_labels,
    retract_marker,
)
from ._engine.magmatic import A, B, Y, comb, leaves, magmatic, star
from ._engine.mechanical import theta_X
from .hamiltonian import Hamiltonian
from .polynomial import positive_int


def binary_trees(leaves_count: int) -> tuple:
    """Enumerate canonical nonplanar commutative binary trees with n leaves."""
    positive_int(leaves_count, "leaves_count")
    return magmatic(leaves_count)


def potential_tree(tree: tuple) -> Hamiltonian:
    """Map a commutative binary tree to its gradient-product potential."""

    def check(node):
        if not isinstance(node, tuple) or len(node) not in (0, 2):
            raise ValueError("A binary tree must be () or a pair of binary trees")
        if node:
            check(node[0])
            check(node[1])
            if node[0] > node[1]:
                raise ValueError("Binary tree pairs must be in canonical sorted order")

    check(tree)
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
