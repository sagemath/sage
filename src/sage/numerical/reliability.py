"""
Utilities for checking numerical reliability of computed results.
"""

def residual_check(f, root, tol=1e-10):
    """
    Check numerical reliability of a computed root.

    INPUT:
        - f: symbolic expression
        - root: numeric root
        - tol: tolerance

    OUTPUT:
        (is_reliable, residual)
    """
    try:
        res = abs(f(x=root))
    except TypeError:
        res = abs(f.subs(x=root))

    return res < tol, res

