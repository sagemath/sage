"""
Numerical Reliability Registry

Central registry for numerical reliability checkers.
"""

class NumericalReliabilityRegistry:
    def __init__(self):
        self._checkers = {}

    def register(self, name, checker):
        if name in self._checkers:
"""
Diagnostics container for numerical reliability checks.
"""

class NumericalReliabilityResult:
    def __init__(self, reliable, residual=None, tolerance=None, message=None, suggestion=None):
        self.reliable = reliable
        self.residual = residual
        self.tolerance = tolerance
        self.message = message
        self.suggestion = suggestion

    def __repr__(self):
        status = "reliable" if self.reliable else "unreliable"
        return f"<NumericalReliabilityResult {status}>"

