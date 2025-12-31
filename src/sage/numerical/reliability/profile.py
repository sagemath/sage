class ReliabilityProfile:
    """
    Semantic description of numerical result reliability.

    This object intentionally carries qualitative information only.
    It does not attempt to prove correctness.
    """

    def __init__(self):
        self.trust_level = "UNKNOWN"
        self.failure_signals = []
        self.assumptions = []

    def add_failure_signal(self, signal):
        if signal not in self.failure_signals:
            self.failure_signals.append(signal)

    def add_assumption(self, assumption):
        if assumption not in self.assumptions:
            self.assumptions.append(assumption)

    def set_trust_level(self, level):
        self.trust_level = level

    def summary(self):
        return {
            "trust_level": self.trust_level,
            "failure_signals": self.failure_signals,
            "assumptions": self.assumptions,
        }

