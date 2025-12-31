import math
from sage.numerical.reliability.profile import ReliabilityProfile

def reliable_sqrt(x):
    value = math.sqrt(x)

    profile = ReliabilityProfile()
    profile.add_assumption("x > 0")

    residual = abs(value * value - x)

    if residual < 1e-12:
        profile.set_trust_level("HIGH")
    else:
        profile.set_trust_level("LOW")
        profile.add_failure_signal("LOSS_OF_SIGNIFICANCE")

    return value, profile

