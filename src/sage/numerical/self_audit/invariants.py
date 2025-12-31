from sage.numerical.self_audit.invariant import ExecutableInvariant

class SquareInvariant(ExecutableInvariant):
    def __init__(self, target):
        self.target = target

    def evaluate(self, x):
        return abs(x*x - self.target)

