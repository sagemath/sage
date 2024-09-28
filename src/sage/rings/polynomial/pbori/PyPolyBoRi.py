# sage_setup: distribution = sagemath-brial
# sage.doctest: needs sage.rings.polynomial.pbori
r"""
PolyBoRi's interface to libpolybori/BRiAL

This file makes interfaces to PolyBoRi's runtime libraries available in Python via sage.


AUTHOR:

- The PolyBoRi Team, 2007-2012

    EXAMPLES::

        sage: from sage.rings.polynomial.pbori.pbori import *
        sage: from sage.rings.polynomial.pbori.blocks import declare_ring
        sage: r=declare_ring(["x0","x1","x2","y0","y1","y2"], globals())
        sage: x0>x1
        True
        sage: x0>x1*x2
        True
        sage: y0>y1
        True
        sage: y0>y1*y2
        True

        sage: r = r.clone(ordering=dlex)
        sage: r(x0) > r(x1)
        True
        sage: r(x0) > r(x1*x2)
        False

        sage: r = r.clone(ordering=dp_asc)
        sage: r(x0) > r(x1)
        False
        sage: r(x0) > r(x1*x2)
        False

        sage: r = r.clone(ordering=block_dlex, blocks=[3])
        sage: r(x0) > r(x1)
        True
        sage: r(x0) > r(x1*x2)
        False
        sage: r(x0) > r(y0*y1*y2)
        True

        sage: r = r.clone(ordering=block_dp_asc)
        sage: r(x0) > r(x1)
        False
        sage: r(x0) > r(y0)
        False
        sage: r(x0) > r(x1*x2)
        False

        sage: r = r.clone(ordering=block_dp_asc, blocks=[3])
        sage: r(x0) > r(y0)
        True

        sage: r(x0) > r(y0*y1)
        True

        sage: r = r.clone(names=["z17", "z7"])
        sage: [r.variable(idx) for idx in range(3)]
        [z17, z7, x2]
        sage: r = r.clone(names='abcde')
        sage: [r.variable(idx) for idx in range(6)]
        [a, b, c, d, e, y2]
"""

from .pbori import (order_dict, TermOrder_from_pb_order, BooleanPolynomialRing,
                    BooleanPolynomialVector, MonomialFactory,
                    PolynomialFactory, VariableFactory, add_up_polynomials)
from .pbori import gauss_on_polys as _gauss_on_polys

import weakref


OrderCode = type('OrderCode', (object,), order_dict)


def Ring(n, order='lp', names=None, blocks=None):
    if blocks is None:
        blocks = []
    pbnames = names
    if pbnames is None:
        pbnames = ['x(' + str(idx) + ')' for idx in range(n)]
    order = TermOrder_from_pb_order(n, order, blocks)
    return BooleanPolynomialRing(n, names=pbnames, order=order)


BoolePolynomialVector = BooleanPolynomialVector


# todo: PolyBoRi's original interface uses its WeakRingPtr here
def WeakRingRef(ring):
    return weakref.weakref(ring)


Monomial = MonomialFactory()
Polynomial = PolynomialFactory()
Variable = VariableFactory()


_add_up_polynomials = add_up_polynomials


def add_up_polynomials(polys, init):
    r"""
    Add up the polynomials in polys (which should be a
    ``BoolePolynomialVector`` or a sequence of ???
    """
    if not isinstance(polys, BoolePolynomialVector):
        vec = BoolePolynomialVector
        for p in polys:
            vec.append(p)
        polys = vec

    return _add_up_polynomials(polys, init)


def gauss_on_polys(l):
    vec = BoolePolynomialVector(l)
    return list(_gauss_on_polys(vec))
