"""
Residual-based reliability checks for numerical roots.

This module provides certificate-based reliability checks
for numerical root approximations.
"""

from sage.numerical.reliability_certify import certify_residual


def certify_root(f, root, tol=1e-10):
    """
    Certify the numerical reliability of a root approximation.

    INPUT:

    - ``f`` -- callable
    - ``root`` -- numeric value
    - ``tol`` -- tolerance (default: 1e-10)

    OUTPUT:

    - ``NumericalReliabilityCertificate``

    EXAMPLES::

        sage: from sage.numerical.root_reliability import certify_root
        sage: f = lambda x: x^2 - 2
        sage: r = sqrt(2).n()
        sage: cert = certify_root(f, r)
        sage: cert.passed
        True
    """
    return certify_residual(f, root, tol)


def residual_check(f, root, tol=1e-10):
    """
    Boolean wrapper for backward compatibility.

    This preserves the old API while internally using
    certificate-based reliability checking.
    """
    cert = certify_root(f, root, tol)
    return cert.passed

