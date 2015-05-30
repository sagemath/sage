"""
Solving quadratic equations.

Interface to the PARI/GP quadratic forms code of Denis Simon.

AUTHORS:

- Denis Simon (GP code)

- Nick Alexander (Sage interface)

- Jeroen Demeyer (2014-09-23): use PARI instead of GP scripts,
  return vectors instead of tuples (:trac:`16997`).
"""

#*****************************************************************************
#       Copyright (C) 2008 Nick Alexander
#       Copyright (C) 2014 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, QQ
from sage.modules.free_module_element import vector

def qfsolve(G):
    r"""
    Find a solution `x = (x_0,...,x_n)` to `x G x^t = 0` for an
    `n \times n`-matrix ``G`` over `\QQ`.

    OUTPUT:

    If a solution exists, return a vector of rational numbers `x`.
    Otherwise, returns `-1` if no solution exists over the reals or a
    prime `p` if no solution exists over the `p`-adic field `\QQ_p`.

    ALGORITHM:

    Uses PARI/GP function ``qfsolve``.

    EXAMPLES::

        sage: from sage.quadratic_forms.qfsolve import qfsolve
        sage: M = Matrix(QQ, [[0, 0, -12], [0, -12, 0], [-12, 0, -1]]); M
        [  0   0 -12]
        [  0 -12   0]
        [-12   0  -1]
        sage: sol = qfsolve(M); sol
        (1, 0, 0)
        sage: sol.parent()
        Vector space of dimension 3 over Rational Field

        sage: M = Matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sage: ret = qfsolve(M); ret
        -1
        sage: ret.parent()
        Integer Ring

        sage: M = Matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, -7]])
        sage: qfsolve(M)
        7

        sage: M = Matrix(QQ, [[3, 0, 0, 0], [0, 5, 0, 0], [0, 0, -7, 0], [0, 0, 0, -11]])
        sage: qfsolve(M)
        (3, -4, -3, -2)
    """
    ret = G._pari_().qfsolve()
    if ret.type() == 't_COL':
        return vector(QQ, ret)
    return ZZ(ret)

def qfparam(G, sol):
    r"""
    Parametrizes the conic defined by the matrix ``G``.

    INPUT:

     - ``G`` -- a `3 \times 3`-matrix over `\QQ`.

     - ``sol`` -- a triple of rational numbers providing a solution
       to sol*G*sol^t = 0.

    OUTPUT:

    A triple of polynomials that parametrizes all solutions of
    x*G*x^t = 0 up to scaling.

    ALGORITHM:

    Uses PARI/GP function ``qfparam``.

    EXAMPLES::

        sage: from sage.quadratic_forms.qfsolve import qfsolve, qfparam
        sage: M = Matrix(QQ, [[0, 0, -12], [0, -12, 0], [-12, 0, -1]]); M
        [  0   0 -12]
        [  0 -12   0]
        [-12   0  -1]
        sage: sol = qfsolve(M)
        sage: ret = qfparam(M, sol); ret
        (-12*t^2 - 1, 24*t, 24)
        sage: ret.parent()
        Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in t over Rational Field
    """
    R = QQ['t']
    t = R.gen()
    mat = G._pari_().qfparam(sol)
    # Interpret the rows of mat as coefficients of polynomials
    return vector(R, mat.Col())
