class NumericalReliabilityChecker:
    """
    Base class for numerical reliability checks.
    """

    name = "base"

    def check(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement check()")

