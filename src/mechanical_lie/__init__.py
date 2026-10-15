"""Exact symbolic Lie/Poisson operations for Euclidean mechanical Hamiltonians."""

from ._engine.errors import ValidationError
from .algebra import binary_trees, potential_tree
from .hamiltonian import (
    Hamiltonian,
    RealizationPlan,
    bracket,
    commutative_product,
    kinetic,
    potential,
)
from .linear import ModularRankReport, modular_rank, modular_rank_report
from .polynomial import Polynomial, phase_space, poisson

__version__ = "0.2.0"
__all__ = [
    "Hamiltonian",
    "Polynomial",
    "RealizationPlan",
    "ModularRankReport",
    "ValidationError",
    "binary_trees",
    "bracket",
    "commutative_product",
    "kinetic",
    "modular_rank",
    "modular_rank_report",
    "phase_space",
    "poisson",
    "potential",
    "potential_tree",
]
