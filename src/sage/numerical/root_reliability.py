"""
Residual-based reliability checks for numerical roots.

EXAMPLES::

    sage: from sage.numerical.root_reliability import residual_check
    sage: f = lambda x: x^2 - 2
    sage: r = sqrt(2).n()
    sage: residual_check(f, r)
    True
"""

def residual_check(f, root, tol=1e-10):
    """
    Check if a numerical root is reliable by residual size.

    INPUT:

    - ``f`` -- callable
    - ``root`` -- numeric value
    - ``tol`` -- tolerance (default: 1e-10)

    OUTPUT:

    - ``True`` or ``False``
    """
    try:
        return abs(f(root)) <= tol
    except Exception:
        return False

