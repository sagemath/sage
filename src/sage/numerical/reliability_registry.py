class NumericalReliabilityRegistry:
    """
    Registry for numerical reliability checkers.
    """

    def __init__(self):
        self._registry = {}

    def register(self, checker):
        self._registry[checker.name] = checker

    def get(self, name):
        return self._registry.get(name)

    def available(self):
        return sorted(self._registry.keys())

