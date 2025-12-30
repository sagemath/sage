"""
Numerical Reliability Registry

Central registry for numerical reliability checkers.
"""

class NumericalReliabilityRegistry:
    def __init__(self):
        self._checkers = {}

    def register(self, name, checker):
        if name in self._checkers:
            raise ValueError(f"Checker '{name}' already registered")
        self._checkers[name] = checker

    def check(self, name, *args, **kwargs):
        if name not in self._checkers:
            raise KeyError(f"No checker registered under '{name}'")
"""
Numerical Reliability Registry

Central registry for numerical reliability checkers.
"""

class NumericalReliabilityRegistry:
    def __init__(self):
        self._checkers = {}

    def register(self, name, checker):
        if name in self._checkers:
            raise ValueError(f"Checker '{name}' already registered")
        self._checkers[name] = checker

    def check(self, name, *args, **kwargs):
        if name not in self._checkers:
            raise KeyError(f"No checker registered under '{name}'")
        return self._checkers[name].check(*args, **kwargs)


# Global default registry
default_registry = NumericalReliabilityRegistry()

