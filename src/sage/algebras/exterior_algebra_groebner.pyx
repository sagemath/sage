# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules
r"""
Exterior algebras Gröbner bases

This contains the backend implementations in Cython for the Gröbner bases
of exterior algebra.

AUTHORS:

- Trevor K. Karn, Travis Scrimshaw (July 2022): Initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2022 Trevor K. Karn <karnx018 at umn.edu>
#                 (C) 2022 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_check
from sage.libs.gmp.mpz cimport mpz_sizeinbase, mpz_setbit, mpz_tstbit, mpz_sgn
from sage.data_structures.bitset_base cimport (bitset_init, bitset_first,
                                               bitset_next, bitset_set_to, bitset_len)
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport richcmp, rich_to_bool
from sage.data_structures.blas_dict cimport iaxpy
from copy import copy

cdef inline long degree(FrozenBitset X) noexcept:
    """
    Compute the degree of ``X``.
    """
    return bitset_len(X._bitset)


cdef inline CliffordAlgebraElement build_monomial(Parent E, FrozenBitset supp):
    """
    Helper function for the fastest way to build a monomial.
    """
    return <CliffordAlgebraElement> E.element_class(E, {supp: E._base.one()})

cdef class GBElement:
    """
    Helper class for storing an element with its leading support both as
    a :class:`FrozenBitset` and an :class:`Integer`.
    """
    def __init__(self, CliffordAlgebraElement x, FrozenBitset ls, Integer n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import GBElement
            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: X = GBElement(a, a.leading_support(), 1)
            sage: TestSuite(X).run(skip='_test_pickling')
        """
        self.elt = x
        self.ls = ls
        self.lsi = n

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import GBElement
            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: X = GBElement(a, a.leading_support(), 1)
            sage: hash(X) == 1
            True
        """
        return int(self.lsi)

    def __richcmp__(self, other, int op):
        """
        Rich compare ``self`` with ``other`` by ``op``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import GBElement
            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: X = GBElement(a, a.leading_support(), 1)
            sage: Y = GBElement(a*b, (a*b).leading_support(), 3)
            sage: X == X
            True
            sage: X == Y
            False
            sage: X != Y
            True
        """
        if self is other:
            return rich_to_bool(op, 0)
        return richcmp(self.elt, (<GBElement> other).elt, op)

cdef class GroebnerStrategy:
    """
    A strategy for computing a Gröbner basis.

    INPUT:

    - ``I`` -- the ideal
    """
    def __init__(self, I):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import GroebnerStrategy
            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([a + 1], side='left')
            sage: G = GroebnerStrategy(I)
            sage: TestSuite(G).run(skip='_test_pickling')
        """
        self.ideal = I
        self.groebner_basis = (None,)
        self.E = <Parent> I.ring()
        self.homogeneous = I._homogeneous
        self.rank = Integer(self.E.ngens())
        if self.homogeneous or I.side() == "left":
            self.side = 0
        elif I.side() == "right":
            self.side = 1
        else:
            self.side = 2

    cdef inline FrozenBitset leading_support(self, CliffordAlgebraElement f):
        """
        Return the leading support of the exterior algebra element ``f``.
        """
        cdef dict mc = <dict> f._monomial_coefficients
        return self.int_to_bitset(max(self.bitset_to_int(k) for k in mc))

    cdef inline partial_S_poly_left(self, GBElement f, GBElement g):
        r"""
        Compute one half of the `S`-polynomial for ``f`` and ``g``.

        This computes:

        .. MATH::

            LCM(LM(f), LM(g)) / LT(f) \cdot f.
        """
        cdef FrozenBitset D = <FrozenBitset> g.ls.difference(f.ls)
        cdef GBElement ret = self.prod_term_GB(D, f)
        inv = ~ret.elt[ret.ls]
        for k in ret.elt._monomial_coefficients:
            ret.elt._monomial_coefficients[k] *= inv
        return ret

    cdef inline partial_S_poly_right(self, GBElement f, GBElement g):
        r"""
        Compute one half of the `S`-polynomial for ``f`` and ``g``.

        This computes:

        .. MATH::

            f \cdot LCM(LM(f), LM(g)) / LT(f).
        """
        cdef FrozenBitset D = <FrozenBitset> g.ls.difference(f.ls)
        cdef GBElement ret = self.prod_GB_term(f, D)
        inv = ~ret.elt[ret.ls]
        for k in ret.elt._monomial_coefficients:
            ret.elt._monomial_coefficients[k] *= inv
        return ret

    cdef inline GBElement build_elt(self, CliffordAlgebraElement f):
        """
        Convert ``f`` into a ``GBElement``.
        """
        cdef dict mc = <dict> f._monomial_coefficients
        # if not mc:
        #     return GBElement(f, FrozenBitset(), -1)
        cdef Integer r = <Integer> max(self.bitset_to_int(k) for k in mc)
        return GBElement(f, self.int_to_bitset(r), r)

    cdef inline GBElement prod_GB_term(self, GBElement f, FrozenBitset t):
        """
        Return the GBElement corresponding to ``f * t``.

        .. WARNING::

            This assumes the leading support is ``f.ls._union(t)``.
        """
        ret = f.elt._mul_self_term(self.E, self.E._base.one())
        cdef FrozenBitset ls = <FrozenBitset> f.ls._union(t)
        return GBElement(<CliffordAlgebraElement> ret, ls, self.bitset_to_int(ls))

    cdef inline GBElement prod_term_GB(self, FrozenBitset t, GBElement f):
        """
        Return the GBElement corresponding to ``t * f``.

        .. WARNING::

            This assumes the leading support is ``f.ls._union(t)``.
        """
        ret = f.elt._mul_term_self(t, self.E._base.one())
        cdef FrozenBitset ls = <FrozenBitset> f.ls._union(t)
        return GBElement(<CliffordAlgebraElement> ret, ls, self.bitset_to_int(ls))

    cdef inline bint build_S_poly(self, GBElement f, GBElement g) noexcept:
        r"""
        Check to see if we should build the `S`-polynomial.

        For homogeneous ideals, we throw out all those pairs `(f, g)` such that

        .. MATH::

            LCM(LM(f), LM(g)) == LM(f) \cdot LM(g).
        """
        if not self.homogeneous:
            return True

        return (<FrozenBitset> f.ls.intersection(g.ls)).isempty()

    cdef inline set preprocessing(self, list P, list G):
        """
        Perform the preprocessing step.
        """
        cdef GBElement f, g, f0, f1
        cdef set additions

        cdef set L = set()
        if self.side == 1:
            for f0, f1 in P:
                if self.build_S_poly(f0, f1):
                    L.add(self.partial_S_poly_right(f0, f1))
                    L.add(self.partial_S_poly_right(f1, f0))
        else:  # We compute a left Gröbner basis for two-sided ideals
            for f0, f1 in P:
                if self.build_S_poly(f0, f1):
                    L.add(self.partial_S_poly_left(f0, f1))
                    L.add(self.partial_S_poly_left(f1, f0))

        if self.side == 2 and not self.homogeneous:
            # Add in all S-poly times positive degree monomials
            one = self.E._base.one()
            additions = set((<GBElement> f).elt._mul_self_term(t, one) for t in self.E._indices for f in L)
            L.update(self.build_elt(f) for f in additions if f)

        cdef set done = set((<GBElement> f).ls for f in L)
        cdef set monL = set()
        for f in L:
            monL.update(f.elt._monomial_coefficients)
        monL.difference_update(done)

        while monL:
            m = self.int_to_bitset(max(self.bitset_to_int(k) for k in monL))
            done.add(m)
            monL.remove(m)
            for g in G:
                lm = g.ls
                if lm <= m:
                    f = self.prod_term_GB(<FrozenBitset> m.difference(lm), g)
                    if f in L:
                        break
                    monL.update(set(f.elt._monomial_coefficients) - done)
                    L.add(f)
                    break
        return L

    cdef inline list reduction(self, list P, list G):
        """
        Perform the reduction of ``P`` mod ``G`` in ``E``.
        """
        cdef set L = self.preprocessing(P, G)
        cdef Py_ssize_t i
        from sage.matrix.constructor import matrix
        cdef Integer r = Integer(2) ** self.rank - Integer(1)  # r for "rank" or "reverso"
        M = matrix({(i, r - self.bitset_to_int(<FrozenBitset> m)): c
                    for i, f in enumerate(L)
                    for m, c in (<GBElement> f).elt._monomial_coefficients.items()},
                   sparse=True)
        M.echelonize()  # Do this in place
        lead_supports = set((<GBElement> f).lsi for f in L)
        return [GBElement(self.E.element_class(self.E, {self.int_to_bitset(r - Integer(j)): c for j, c in M[i].iteritems()}),
                          self.int_to_bitset(Integer(r - p)),
                          Integer(r - p))
                for i, p in enumerate(M.pivots())
                if r - Integer(p) not in lead_supports]

    def compute_groebner(self, reduced=True):
        r"""
        Compute the (``reduced``) left/right Gröbner basis for the ideal ``I``.

        EXAMPLES::

            sage: E.<y, x> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x*y - x, x*y - 1], side='left')
            sage: I.groebner_basis()  # indirect doctest
            (1,)
            sage: J = E.ideal([x*y - x, 2*x*y - 2], side='left')
            sage: J.groebner_basis()  # indirect doctest
            (1,)

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([a+b*c], side='left')
            sage: I.groebner_basis()  # indirect doctest
            (b*c + a,)
            sage: I = E.ideal([a+b*c], side='twosided')
            sage: I.groebner_basis()  # indirect doctest
            (a*b, a*c, b*c + a, a*d)
        """
        cdef FrozenBitset p0, p1
        cdef long deg
        cdef Py_ssize_t i, j
        cdef GBElement f0, f1
        cdef list G = [], Gp
        cdef dict constructed = {}
        cdef CliffordAlgebraElement f

        for f in self.ideal.gens():
            if not f:  # Remove 0s
                continue
            f0 = self.build_elt(f)
            if f0.lsi in constructed:
                if f0 in constructed[f0.lsi]:  # Already there
                    continue
                constructed[f0.lsi].add(f0)
            else:
                constructed[f0.lsi] = set([f0])
            G.append(f0)

        if self.side == 2 and not self.homogeneous:
            # Add in all S-poly times positive degree monomials
            one = self.E._base.one()
            for t in self.E._indices:
                for f0 in G:
                    f = f0.elt._mul_self_term(t, one)
                    if not f:
                        continue
                    f1 = self.build_elt(f)
                    if f1.lsi in constructed:
                        if f1 in constructed[f1.lsi]:  # Already there
                            continue
                        constructed[f1.lsi].add(f1)
                    else:
                        constructed[f1.lsi] = set([f1])
                    G.append(f1)

        cdef Py_ssize_t n = len(G)
        cdef dict P = {}

        for i in range(n):
            f0 = <GBElement> G[i]
            p0 = f0.ls
            for j in range(i+1, n):
                f1 = <GBElement> G[j]
                p1 = f1.ls
                deg = degree(<FrozenBitset> (p0._union(p1)))
                if deg in P:
                    P[deg].append((f0, f1))
                else:
                    P[deg] = [(f0, f1)]

        while P:
            sig_check()
            Pp = P.pop(min(P))  # The selection: lowest lcm degree
            Gp = self.reduction(Pp, G)
            # Add the elements Gp to G when a new element is found
            for f0 in Gp:
                if f0.lsi in constructed:
                    if f0 in constructed[f0.lsi]:  # Already there
                        continue
                    constructed[f0.lsi].add(f0)
                else:
                    constructed[f0.lsi] = set([f0])
                G.append(f0)
            # Find the degrees of the new pairs
            for j in range(n, len(G)):
                f1 = G[j]
                p1 = f1.ls
                for i in range(j):
                    f0 = G[i]
                    p0 = f0.ls
                    deg = degree(<FrozenBitset> (p0._union(p1)))
                    if deg in P:
                        P[deg].append((f0, f1))
                    else:
                        P[deg] = [(f0, f1)]
            n = len(G)

        G.sort(key=lambda x: (<GBElement> x).lsi)

        if not reduced:
            self.groebner_basis = tuple([(<GBElement> f0).elt for f0 in G if (<GBElement> f0).elt])
            return
        self.reduced_gb(G)

    cdef int reduced_gb(self, list G) except -1:
        """
        Convert the Gröbner basis ``G`` into a reduced Gröbner basis.
        """
        cdef Py_ssize_t i, j, k
        cdef GBElement f0, f1

        # Now that we have a Gröbner basis, we make this into a reduced Gröbner basis
        cdef Py_ssize_t num_zeros = 0
        cdef Py_ssize_t n = len(G)
        cdef set pairs = set((i, j) for i in range(n) for j in range(n) if i != j)

        while pairs:
            sig_check()
            i, j = pairs.pop()
            f0 = <GBElement> G[i]
            f1 = <GBElement> G[j]
            assert f0.elt._monomial_coefficients is not f1.elt._monomial_coefficients, (i, j)
            # We perform the classical reduction algorithm here on each pair
            # TODO: Make this faster by using the previous technique?
            if self.reduce_single(f0.elt, f1.elt):
                if f0.elt:
                    G[i] = self.build_elt(f0.elt)
                    pairs.update((k, i) for k in range(n) if k != i)
                else:
                    G[i] = GBElement(f0.elt, FrozenBitset(), Integer(2)**self.rank + 1)
                    num_zeros += 1
                    pairs.difference_update((k, i) for k in range(n) if k != i)
                    pairs.difference_update((i, k) for k in range(n) if k != i)

        G.sort(key=lambda x: (<GBElement> x).lsi)
        for i in range(len(G)-num_zeros):
            f0 = <GBElement> G[i]
            if f0.elt:
                inv = ~f0.elt[f0.ls]
                for key in f0.elt._monomial_coefficients:
                    f0.elt._monomial_coefficients[key] *= inv
        self.groebner_basis = tuple([f0.elt for f0 in G[:len(G)-num_zeros]])
        return 0

    def reduce_computed_gb(self):
        """
        Convert the computed Gröbner basis to a reduced Gröbner basis.

            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x+y*z])
            sage: I.groebner_basis(reduced=False)
            (x*y, x*z, y*z + x, x*y*z)
            sage: I._groebner_strategy.reduce_computed_gb()
            sage: I._groebner_strategy.groebner_basis
            (x*y, x*z, y*z + x)
        """
        if self.groebner_basis == [(None,)]:
            raise ValueError("Gröbner basis not yet computed")
        cdef list G = [self.build_elt(f) for f in self.groebner_basis]
        self.reduced_gb(G)

    cpdef CliffordAlgebraElement reduce(self, CliffordAlgebraElement f):
        """
        Reduce ``f`` modulo the ideal with Gröbner basis ``G``.

        EXAMPLES::

            sage: E.<a,b,c,d,e> = ExteriorAlgebra(QQ)
            sage: rels = [c*d*e - b*d*e + b*c*e - b*c*d,
            ....:         c*d*e - a*d*e + a*c*e - a*c*d,
            ....:         b*d*e - a*d*e + a*b*e - a*b*d,
            ....:         b*c*e - a*c*e + a*b*e - a*b*c,
            ....:         b*c*d - a*c*d + a*b*d - a*b*c]
            sage: I = E.ideal(rels)
            sage: GB = I.groebner_basis()
            sage: I._groebner_strategy.reduce(a*b*e)
            a*b*e
            sage: I._groebner_strategy.reduce(b*d*e)
            a*b*d - a*b*e + a*d*e
            sage: I._groebner_strategy.reduce(c*d*e)
            a*c*d - a*c*e + a*d*e
            sage: I._groebner_strategy.reduce(a*b*c*d*e)
            0
            sage: I._groebner_strategy.reduce(a*b*c*d)
            0
            sage: I._groebner_strategy.reduce(E.zero())
            0

        Check :issue:`37108` is fixed::

            sage: E = ExteriorAlgebra(QQ, 6)
            sage: E.inject_variables(verbose=False)
            sage: gens = [-e0*e1*e2 + e0*e1*e5 - e0*e2*e3 - e0*e3*e5 + e1*e2*e3 + e1*e3*e5,
            ....:         e1*e2 - e1*e5 + e2*e5, e0*e2 - e0*e4 + e2*e4,
            ....:         e3*e4 - e3*e5 + e4*e5, e0*e1 - e0*e3 + e1*e3]
            sage: I = E.ideal(gens)
            sage: S = E.quo(I)
            sage: I.reduce(e1*e3*e4*e5)
            0
        """
        if not f:
            return f
        # Make a copy to mutate
        f = type(f)(f._parent, copy(f._monomial_coefficients))
        was_reduced = True
        while was_reduced:
            was_reduced = False
            for g in self.groebner_basis:
                was_reduced = self.reduce_single(f, g)
                if was_reduced:
                    break
        return f

    cdef bint reduce_single(self, CliffordAlgebraElement f, CliffordAlgebraElement g) except -1:
        r"""
        Reduce ``f`` by ``g``.

        .. WARNING::

            This modifies the element ``f``.
        """
        cdef FrozenBitset lm = self.leading_support(g), s, t
        cdef bint did_reduction = True, was_reduced=False
        cdef CliffordAlgebraElement gp

        one = self.E._base.one()
        while did_reduction:
            did_reduction = False
            for s in f._monomial_coefficients:
                if lm.issubset(s):
                    t = s
                    did_reduction = True
                    was_reduced = True
                    break
            if did_reduction:
                if self.side == 0:
                    gp = g._mul_term_self(t - lm, one)
                else:
                    gp = g._mul_self_term(t - lm, one)
                coeff = f[t] / gp._monomial_coefficients[t]
                iaxpy(-coeff, gp._monomial_coefficients, f._monomial_coefficients)
        return was_reduced

    cdef Integer bitset_to_int(self, FrozenBitset X):
        raise NotImplementedError

    cdef FrozenBitset int_to_bitset(self, Integer n):
        raise NotImplementedError

    def sorted_monomials(self, as_dict=False):
        """
        Helper function to display the monomials in their term order
        used by ``self``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import *
            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x, y])
            sage: GroebnerStrategyNegLex(I).sorted_monomials()
            [1, x, y, x*y, z, x*z, y*z, x*y*z]
            sage: GroebnerStrategyDegLex(I).sorted_monomials()
            [1, x, y, z, x*y, x*z, y*z, x*y*z]
            sage: GroebnerStrategyDegRevLex(I).sorted_monomials()
            [1, z, y, x, y*z, x*z, x*y, x*y*z]

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([a, b])
            sage: GroebnerStrategyNegLex(I).sorted_monomials()
            [1,
             a,
             b, a*b,
             c, a*c, b*c, a*b*c,
             d, a*d, b*d, a*b*d, c*d, a*c*d, b*c*d, a*b*c*d]
            sage: GroebnerStrategyDegLex(I).sorted_monomials()
            [1,
             a, b, c, d,
             a*b, a*c, a*d, b*c, b*d, c*d,
             a*b*c, a*b*d, a*c*d, b*c*d,
             a*b*c*d]
            sage: GroebnerStrategyDegRevLex(I).sorted_monomials()
            [1,
             d, c, b, a,
             c*d, b*d, a*d, b*c, a*c, a*b,
             b*c*d, a*c*d, a*b*d, a*b*c,
             a*b*c*d]
        """
        cdef FrozenBitset X
        cdef Integer i
        cdef list D = [self.bitset_to_int(X) for X in self.E._indices]
        D.sort()
        if as_dict:
            return {i: build_monomial(self.E, self.int_to_bitset(i)) for i in D}
        return [build_monomial(self.E, self.int_to_bitset(i)) for i in D]

    def monomial_to_int(self):
        """
        Helper function to display the monomials in their term order
        used by ``self``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import *
            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([a, b])
            sage: GroebnerStrategyDegLex(I).monomial_to_int()
            {1: 0,
             a: 1, b: 2, c: 3, d: 4,
             a*b: 5, a*c: 6, a*d: 7, b*c: 8, b*d: 9, c*d: 10,
             a*b*c: 11, a*b*d: 12, a*c*d: 13, b*c*d: 14,
             a*b*c*d: 15}
            sage: GroebnerStrategyDegRevLex(I).monomial_to_int()
            {1: 0,
             a: 4, b: 3, c: 2, d: 1,
             a*b: 10, a*c: 9, a*d: 7, b*c: 8, b*d: 6, c*d: 5,
             a*b*c: 14, a*b*d: 13, a*c*d: 12, b*c*d: 11,
             a*b*c*d: 15}
        """
        B = self.E.basis()
        return {B[X]: self.bitset_to_int(X) for X in self.E._indices}

cdef class GroebnerStrategyNegLex(GroebnerStrategy):
    """
    Gröbner basis strategy implementing neglex ordering.
    """
    cdef inline Integer bitset_to_int(self, FrozenBitset X):
        """
        Convert ``X`` to an :class:`Integer`.
        """
        cdef Integer ret = Integer(0)
        cdef long elt = bitset_first(X._bitset)
        while elt >= 0:
            mpz_setbit(ret.value, elt)
            elt = bitset_next(X._bitset, elt + 1)
        return ret

    cdef inline FrozenBitset int_to_bitset(self, Integer n):
        """
        Convert a nonnegative integer ``n`` to a :class:`FrozenBitset`.
        """
        cdef size_t i

        if mpz_sgn(n.value) == 0:
            return FrozenBitset()

        cdef FrozenBitset ret = <FrozenBitset> FrozenBitset()
        cdef size_t s = mpz_sizeinbase(n.value, 2)
        bitset_init(ret._bitset, s)
        for i in range(s):
            bitset_set_to(ret._bitset, i, mpz_tstbit(n.value, i))
        return ret

cdef class GroebnerStrategyDegRevLex(GroebnerStrategy):
    """
    Gröbner basis strategy implementing degree revlex ordering.
    """
    cdef inline Integer bitset_to_int(self, FrozenBitset X):
        """
        Convert ``X`` to an :class:`Integer`.
        """
        if X.isempty():
            return Integer(0)

        cdef Integer n = self.rank
        cdef long i, deg = degree(X)
        cdef long r = 1
        cdef Integer t = Integer(0)

        cdef long elt = bitset_first(X._bitset)
        while elt >= 0:
            t += Integer(elt).binomial(r)
            r += 1
            elt = bitset_next(X._bitset, elt + 1)
        return Integer(sum(n.binomial(i) for i in range(deg+1)) - t - 1)

    cdef inline FrozenBitset int_to_bitset(self, Integer n):
        """
        Convert a nonnegative integer ``n`` to a :class:`FrozenBitset`.
        """
        if mpz_sgn(n.value) == 0:
            return FrozenBitset()

        cdef Py_ssize_t deg = 0
        cdef Integer binom = Integer(1)
        while n >= binom:
            n -= binom
            deg += 1
            binom = self.rank.binomial(deg)

        # TODO: Cythonize the from_rank
        from sage.combinat.combination import from_rank
        return FrozenBitset([self.rank - val - 1 for val in from_rank(n, self.rank, deg)])

cdef class GroebnerStrategyDegLex(GroebnerStrategy):
    """
    Gröbner basis strategy implementing degree lex ordering.
    """
    cdef inline Integer bitset_to_int(self, FrozenBitset X):
        """
        Convert ``X`` to an :class:`Integer`.
        """
        if X.isempty():
            return Integer(0)

        cdef Integer n = self.rank
        cdef long i, deg = degree(X)
        cdef long r = deg
        cdef Integer t = Integer(0)

        cdef long elt = bitset_first(X._bitset)
        while elt >= 0:
            t += Integer(n - 1 - elt).binomial(r)
            r -= 1
            elt = bitset_next(X._bitset, elt + 1)
        return Integer(sum(n.binomial(i) for i in range(deg+1)) - t - 1)

    cdef inline FrozenBitset int_to_bitset(self, Integer n):
        """
        Convert a nonnegative integer ``n`` to a :class:`FrozenBitset`.
        """
        if mpz_sgn(n.value) == 0:
            return FrozenBitset()

        cdef Py_ssize_t deg = 0
        cdef Integer binom = Integer(1)
        while n >= binom:
            n -= binom
            deg += 1
            binom = self.rank.binomial(deg)

        from sage.combinat.combination import from_rank
        return FrozenBitset(from_rank(n, self.rank, deg))
