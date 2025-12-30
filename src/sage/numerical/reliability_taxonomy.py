"""
Numerical reliability failure taxonomy.

This module defines common failure categories encountered
in numerical computations.
"""

from enum import Enum, auto


class NumericalFailure(Enum):
    OVERFLOW = auto()
    UNDERFLOW = auto()
    CANCELLATION = auto()
    ILL_CONDITIONED = auto()
    NON_CONVERGENCE = auto()
    LOSS_OF_SIGNIFICANCE = auto()
    DOMAIN_ERROR = auto()

