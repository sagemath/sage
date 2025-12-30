from sage.numerical.reliability_base import NumericalReliabilityChecker


class RootReliabilityChecker(NumericalReliabilityChecker):
    """
    Reliability checker for numerical roots.
    """

    def __init__(self, function):
        self.function = function

    def compute_residual(self, value):
        return self.function(value)


def check_root(function, root, tolerance=1e-12):
    """
    Check numerical reliability of a computed root.

    EXAMPLE::

        sage: f = lambda x: x^2 - 2
        sage: r = sqrt(2).n()
        sage: check_root(f, r)
    """
    checker = RootReliabilityChecker(function)
    return checker.check(root, tolerance=tolerance)

