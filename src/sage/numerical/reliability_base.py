"""
Base classes for numerical reliability checks.

This module defines a common interface for checking the numerical
reliability of approximate solutions produced by numerical algorithms.
"""

class NumericalReliabilityResult:
    """
    Stores the result of a numerical reliability check.
    """

    def __init__(self, value, residual, reliable, tolerance):
        self.value = value
        self.residual = residual
        self.reliable = reliable
        self.tolerance = tolerance

    def __repr__(self):
        status = "reliable" if self.reliable else "unreliable"
        return (
            f"NumericalReliabilityResult("
            f"value={self.value}, residual={self.residual}, "
            f"tolerance={self.tolerance}, status={status})"
        )


class NumericalReliabilityChecker:
    """
    Abstract base class for numerical reliability checks.
    """

    def compute_residual(self, *args, **kwargs):
        raise NotImplementedError("compute_residual must be implemented")

    def check(self, value, tolerance=1e-12, *args, **kwargs):
        residual = self.compute_residual(value, *args, **kwargs)
        reliable = abs(residual) <= tolerance
        return NumericalReliabilityResult(
            value=value,
            residual=residual,
            reliable=reliable,
            tolerance=tolerance,
        )

