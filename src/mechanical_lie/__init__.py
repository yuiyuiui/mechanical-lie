"""Exact symbolic Lie/Poisson operations for Euclidean mechanical Hamiltonians."""

from ._engine.errors import ValidationError
from .algebra import binary_trees, potential_tree
from .hamiltonian import Hamiltonian, bracket, commutative_product, kinetic, potential
from .linear import modular_rank
from .polynomial import Polynomial, phase_space, poisson

__version__ = "0.1.0"
__all__ = [
    "Hamiltonian",
    "Polynomial",
    "ValidationError",
    "binary_trees",
    "bracket",
    "commutative_product",
    "kinetic",
    "modular_rank",
    "phase_space",
    "poisson",
    "potential",
    "potential_tree",
]
