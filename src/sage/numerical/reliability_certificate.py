class NumericalReliabilityCertificate:
    """
    Certificate explaining why a numerical result is reliable.
    """

    def __init__(self, method, metrics, tolerance, passed):
        self.method = method
        self.metrics = metrics
        self.tolerance = tolerance
        self.passed = passed

    def explain(self):
        lines = [f"Method: {self.method}"]
        for k, v in self.metrics.items():
            lines.append(f"{k}: {v}")
        lines.append(f"Tolerance: {self.tolerance}")
        lines.append("Result: PASSED" if self.passed else "Result: FAILED")
        return "\n".join(lines)

