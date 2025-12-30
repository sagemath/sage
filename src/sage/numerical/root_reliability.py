from sage.numerical.reliability_base import NumericalReliabilityChecker
from sage.numerical.reliability_diagnostics import NumericalReliabilityResult


class RootResidualReliability(NumericalReliabilityChecker):
    """
    Residual-based reliability check for numerical roots.
    """

    name = "root_residual"

    def check(self, f, root, tol=1e-10):
        try:
            residual = abs(f(root))
            return NumericalReliabilityResult(
                reliable=residual <= tol,
                residual=residual,
                message="Residual check passed"
                if residual <= tol
                else "Residual too large",
            )
        except Exception as e:
            return NumericalReliabilityResult(
                reliable=False,
                residual=None,
                message=str(e),
            )


def residual_check(f, root, tol=1e-10):
    """
    Convenience wrapper for residual-based root reliability.
    """
    checker = RootResidualReliability()
    return checker.check(f, root, tol)

