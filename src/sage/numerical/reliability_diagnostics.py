class NumericalReliabilityResult:
    """
    Container for numerical reliability diagnostics.
    """

    def __init__(self, reliable, residual=None, message=""):
        self.reliable = bool(reliable)
        self.residual = residual
        self.message = message

    def __bool__(self):
        return self.reliable

    def __repr__(self):
        return (
            f"NumericalReliabilityResult("
            f"reliable={self.reliable}, "
            f"residual={self.residual}, "
            f"message={self.message!r})"
        )

