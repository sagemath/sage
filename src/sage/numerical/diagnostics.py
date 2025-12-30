# sage/numerical/diagnostics.py

class NumericalDiagnostics:
    """
    Container for numerical reliability diagnostics.

    This object stores quantitative indicators that help decide
    whether a numerical result can be trusted.

    EXAMPLES::

        sage: from sage.numerical.diagnostics import NumericalDiagnostics
        sage: d = NumericalDiagnostics(residual=1e-5, precision=53)
        sage: d.is_reliable(tol=1e-6)
        False
    """

    def __init__(self, residual=None, condition=None, precision=None, message=None):
        self.residual = residual
        self.condition = condition
        self.precision = precision
        self.message = message

    def is_reliable(self, tol=1e-10):
        """
        Decide whether the numerical result is reliable.

        INPUT:
        - ``tol`` -- tolerance for residual-based reliability
        """
        if self.residual is None:
            return True
        return self.residual <= tol

    def summary(self):
        """
        Return a short human-readable summary.
        """
        parts = []
        if self.residual is not None:
            parts.append(f"residual={self.residual}")
        if self.condition is not None:
            parts.append(f"condition={self.condition}")
        if self.precision is not None:
            parts.append(f"precision={self.precision}")
        if self.message:
            parts.append(self.message)
        return ", ".join(parts)

