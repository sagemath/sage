# sage_setup: distribution = sagemath-brial
# sage.doctest: needs sage.rings.polynomial.pbori
r"""
parallel.py
PolyBoRi

Created by Michael Brickenstein on 2008-10-31.
Copyright 2008 The PolyBoRi Team
"""
from zlib import compress, decompress
import copyreg

from .pbori import if_then_else, BooleSet, CCuddNavigator
from .PyPolyBoRi import (Polynomial, Ring, WeakRingRef, Monomial, Variable)
from .gbcore import groebner_basis


def to_fast_pickable(l):
    r"""
    Convert a list of polynomials into a builtin Python value, which is fast pickable and compact.

    INPUT:

    - ``l`` -- a list of Boolean polynomials

    OUTPUT:

    It is converted to a tuple consisting of
    - codes referring to the polynomials
    - list of conversions of nodes.
    The nodes are sorted, so that n occurs before n.else_branch(), n.then_branch()
    nodes are only listed, if they are not constant.

    A node is converted in this way:
    0 -> 0
    1 -> 1
    if_then_else(v,t,e) -> (v, index of then branch +2, index of else branch +2)
    the shift of +2 is for the constant values implicitly contained in the list.
    Each code c refers to the c-2-th position in the conversion list, if c >=2, else to
    the corresponding Boolean constant if c in {0, 1}

    EXAMPLES::

        sage: from sage.rings.polynomial.pbori import Ring, Polynomial
        sage: from sage.rings.polynomial.pbori.parallel import to_fast_pickable, from_fast_pickable
        sage: r = Ring(1000)
        sage: x = r.variable
        sage: to_fast_pickable([Polynomial(1, r)])
        [[1], []]
        sage: to_fast_pickable([Polynomial(0, r)])
        [[0], []]
        sage: to_fast_pickable([x(0)])
        [[2], [(0, 1, 0)]]
        sage: to_fast_pickable([x(0)*x(1)+x(1)])
        [[2], [(0, 3, 3), (1, 1, 0)]]
        sage: to_fast_pickable([x(1)])
        [[2], [(1, 1, 0)]]
        sage: to_fast_pickable([x(0)+1])
        [[2], [(0, 1, 1)]]
        sage: to_fast_pickable([x(0)*x(1)])
        [[2], [(0, 3, 0), (1, 1, 0)]]
        sage: to_fast_pickable([x(0)*x(1)+x(1)])
        [[2], [(0, 3, 3), (1, 1, 0)]]
        sage: to_fast_pickable([x(0)*x(1)+x(2)])
        [[2], [(0, 3, 4), (1, 1, 0), (2, 1, 0)]]
        sage: p=x(5)*x(23) + x(5)*x(24)*x(59) + x(5) + x(6)*x(23)*x(89) + x(6)*x(60)*x(89) + x(23) + x(24)*x(89) + x(24) + x(60)*x(89) + x(89) + 1
        sage: from_fast_pickable(to_fast_pickable([p]), r)==[p]
        True
        sage: to_fast_pickable([x(0)*x(1), Polynomial(0, r), Polynomial(1, r), x(3)])
        [[2, 0, 1, 4], [(0, 3, 0), (1, 1, 0), (3, 1, 0)]]
    """
    if not l:
        return [[], []]

    f = l[0]
    f = f.set()
    r = f.ring()
    one = r.one().navigation()
    zero = r.zero().navigation()
    nodes = set()

    def find_navs(nav):
        if nav not in nodes and not nav.constant():
            nodes.add(nav)
            find_navs(nav.then_branch())
            find_navs(nav.else_branch())
    for f in l:
        f_nav = f.set().navigation()
        find_navs(f_nav)

    nodes_sorted = sorted(nodes, key=CCuddNavigator.value)
    nodes2i = {one: 1, zero: 0}
    for (i, n) in enumerate(nodes_sorted):
        nodes2i[n] = i + 2

    for i in range(len(nodes_sorted)):
        n = nodes_sorted[i]
        t = nodes2i[n.then_branch()]
        e = nodes2i[n.else_branch()]
        nodes_sorted[i] = (n.value(), t, e)

    return [[nodes2i[f.set().navigation()] for f in l], nodes_sorted]


def from_fast_pickable(l, r):
    r"""
    Undo the operation :func:`to_fast_pickable`.

    The first argument is an object created by :func:`to_fast_pickable`.

    For the specified format, see the documentation of :func:`to_fast_pickable`.
    The second argument is ring, in which this polynomial should be created.

    INPUT:

    See OUTPUT of :func:`to_fast_pickable`.

    OUTPUT: list of Boolean polynomials

    EXAMPLES::

        sage: from sage.rings.polynomial.pbori import Ring
        sage: from sage.rings.polynomial.pbori.parallel import from_fast_pickable
        sage: r = Ring(1000)
        sage: x = r.variable
        sage: from_fast_pickable([[1], []], r)
        [1]
        sage: from_fast_pickable([[0], []], r)
        [0]
        sage: from_fast_pickable([[2], [(0, 1, 0)]], r)
        [x(0)]
        sage: from_fast_pickable([[2], [(1, 1, 0)]], r)
        [x(1)]
        sage: from_fast_pickable([[2], [(0, 1, 1)]], r)
        [x(0) + 1]
        sage: from_fast_pickable([[2], [(0, 3, 0), (1, 1, 0)]], r)
        [x(0)*x(1)]
        sage: from_fast_pickable([[2], [(0, 3, 3), (1, 1, 0)]], r)
        [x(0)*x(1) + x(1)]
        sage: from_fast_pickable([[2], [(0, 3, 4), (1, 1, 0), (2, 1, 0)]], r)
        [x(0)*x(1) + x(2)]
        sage: from_fast_pickable([[2, 0, 1, 4], [(0, 3, 0), (1, 1, 0), (3, 1, 0)]], r)
        [x(0)*x(1), 0, 1, x(3)]
    """
    i2poly = {0: r.zero(), 1: r.one()}
    (indices, terms) = l

    for i in reversed(range(len(terms))):
        (v, t, e) = terms[i]
        t = i2poly[t]
        e = i2poly[e]
        terms[i] = if_then_else(v, t, e)
        i2poly[i + 2] = terms[i]
    return [Polynomial(i2poly[i]) for i in indices]


def _calculate_gb_with_keywords(args):
    (I, kwds_as_single_arg) = args
    return groebner_basis(I, **kwds_as_single_arg)


def _decode_polynomial(code):
    return from_fast_pickable(*code)[0]


def _encode_polynomial(poly):
    return (to_fast_pickable([poly]), poly.ring())


def pickle_polynomial(self):
    return (_decode_polynomial, (_encode_polynomial(self), ))


copyreg.pickle(Polynomial, pickle_polynomial)


def pickle_bset(self):
    return (BooleSet, (Polynomial(self), ))


copyreg.pickle(BooleSet, pickle_bset)


def pickle_monom(self):
    return (Monomial, (list(self.variables()),))


copyreg.pickle(Monomial, pickle_monom)


def pickle_var(self):
    return (Variable, (self.index(), self.ring()))


copyreg.pickle(Variable, pickle_var)


def _decode_ring(code):
    import os
    (identifier, data, varnames, blocks) = code

    global _polybori_parallel_rings
    try:
        _polybori_parallel_rings
    except NameError:
        _polybori_parallel_rings = {}

    for key in [key for key in _polybori_parallel_rings
                if not _polybori_parallel_rings[key][0]()]:
        del _polybori_parallel_rings[key]

    if identifier in _polybori_parallel_rings:
        ring = _polybori_parallel_rings[identifier][0]()
    else:
        ring = None

    if not ring:
        varnames = decompress(varnames).split('\n')
        (nvars, ordercode) = data
        ring = Ring(nvars, ordercode, names=varnames, blocks=blocks)
        storage_data = (WeakRingRef(ring), code)
        _polybori_parallel_rings[identifier] = storage_data
        _polybori_parallel_rings[(ring.id(), os.getpid())] = storage_data

    return ring


def _encode_ring(ring):
    import os
    identifier = (ring.id(), os.getpid())

    global _polybori_parallel_rings
    try:
        _polybori_parallel_rings
    except NameError:
        _polybori_parallel_rings = {}

    for key in [key for key in _polybori_parallel_rings
                if not _polybori_parallel_rings[key][0]()]:
        del _polybori_parallel_rings[key]

    if identifier in _polybori_parallel_rings:
        code = _polybori_parallel_rings[identifier][1]
    else:
        nvars = ring.n_variables()
        data = (nvars, ring.get_order_code())
        varnames = '\n'.join(str(ring.variable(idx))
                             for idx in range(nvars))
        blocks = list(ring.blocks())
        code = (identifier, data, compress(varnames), blocks[:-1])
        _polybori_parallel_rings[identifier] = (WeakRingRef(ring), code)

    return code


def pickle_ring(self):
    return (_decode_ring, (_encode_ring(self), ))


copyreg.pickle(Ring, pickle_ring)


def groebner_basis_first_finished(I, *l):
    r"""

    INPUT:

    - ``I`` -- ideal
    - ``l`` -- keyword dictionaries, which will be keyword arguments to
      ``groebner_basis``

    OUTPUT:

    - tries to compute ``groebner_basis(I, **kwd)`` for kwd in l
    - returns the result of the first terminated computation

    EXAMPLES::

        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import Ring
        sage: r = Ring(1000)
        sage: ideal = [r.variable(1)*r.variable(2)+r.variable(2)+r.variable(1)]
        sage: from sage.rings.polynomial.pbori.parallel import groebner_basis_first_finished
        sage: groebner_basis_first_finished(ideal, dict(heuristic=True), dict(heuristic=False))
        [x1, x2]
    """
    if not I:
        return []

    from multiprocessing import Pool

    pool = Pool(processes=len(l))
    it = pool.imap_unordered(_calculate_gb_with_keywords,
                             [(I, kwds) for kwds in l])
    res = next(it)

    pool.terminate()

    return res
