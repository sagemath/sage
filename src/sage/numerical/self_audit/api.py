from sage.numerical.self_audit.engine import audit_sqrt

def n_audit(expr, digits=None):
    """
    Return numerical value with audit certificate.
    """
    value = expr.n(digits)
    cert = audit_sqrt(expr, value)
    return value, cert

