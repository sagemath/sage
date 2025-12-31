from sage.numerical.self_audit.certificate import AuditCertificate
from sage.numerical.self_audit.invariants import SquareInvariant

def audit_sqrt(input_value, result):
    cert = AuditCertificate()

    cert.add_assumption("input >= 0")

    invariant = SquareInvariant(input_value)
    residual = invariant.evaluate(result)

    cert.add_invariant("x^2 ≈ {}".format(input_value), residual)
    cert.add_residual("sqrt_residual", residual)

    if residual < 1e-12:
        cert.set_conditioning("well-conditioned")
        cert.set_failure_probability("negligible")
    else:
        cert.set_conditioning("ill-conditioned")
        cert.set_failure_probability("elevated")
        cert.add_warning("High residual detected")

    return cert

