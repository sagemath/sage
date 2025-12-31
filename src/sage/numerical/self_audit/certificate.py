class AuditCertificate:
    """
    Experimental audit record for a numerical computation.
    """

    def __init__(self):
        self.invariants = {}
        self.residuals = {}
        self.assumptions = []
        self.conditioning = "unknown"
        self.failure_probability = "unknown"
        self.warnings = []

    def add_invariant(self, name, value):
        self.invariants[name] = value

    def add_residual(self, name, value):
        self.residuals[name] = value

    def add_assumption(self, text):
        self.assumptions.append(text)

    def add_warning(self, text):
        self.warnings.append(text)

    def set_conditioning(self, value):
        self.conditioning = value

    def set_failure_probability(self, value):
        self.failure_probability = value

    def summary(self):
        return {
            "invariants": self.invariants,
            "residuals": self.residuals,
            "assumptions": self.assumptions,
            "conditioning": self.conditioning,
            "failure_probability": self.failure_probability,
            "warnings": self.warnings,
        }

