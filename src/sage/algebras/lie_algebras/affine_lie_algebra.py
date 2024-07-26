"""
Affine Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.element import parent
from sage.categories.kac_moody_algebras import KacMoodyAlgebras

from sage.algebras.lie_algebras.lie_algebra import LieAlgebra, FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import UntwistedAffineLieAlgebraElement
from sage.combinat.root_system.cartan_type import CartanType
from sage.categories.cartesian_product import cartesian_product
from sage.rings.integer_ring import ZZ
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Set_generic


class AffineLieAlgebra(FinitelyGeneratedLieAlgebra):
    r"""
    An (untwisted) affine Lie algebra.

    Note that the derived subalgebra of the Kac-Moody algebra is the
    affine Lie algebra.

    INPUT:

    Can be one of the following:

    - a base ring and an affine Cartan type: constructs the affine
      (Kac-Moody) Lie algebra of the classical Lie algebra in the
      bracket representation over the base ring

    - a classical Lie algebra: constructs the corresponding affine
      (Kac-Moody) Lie algebra

    There is the optional argument ``kac_moody``, which can be set
    to ``False`` to obtain the affine Lie algebra instead of the affine
    Kac-Moody algebra.

    .. SEEALSO::

        - :class:`UntwistedAffineLieAlgebra`
        - :class:`TwistedAffineLieAlgebra`

    REFERENCES:

    - [Ka1990]_
    """
    @staticmethod
    def __classcall_private__(cls, arg0, cartan_type=None, kac_moody=True):
        """
        Parse input to ensure a unique representation.

        INPUT:

        - ``arg0`` -- a simple Lie algebra or a base ring
        - ``cartan_type`` -- a Cartan type

        EXAMPLES::

            sage: L1 = lie_algebras.Affine(QQ, ['A', 4, 1])
            sage: cl = lie_algebras.sl(QQ, 5)
            sage: L2 = lie_algebras.Affine(cl)
            sage: L1 is L2
            True
            sage: cl.affine() is L1
            True

            sage: L1 = LieAlgebra(QQ, cartan_type=['A', 5, 2])
            sage: L2 = lie_algebras.Affine(QQ, ['A', 5, 2])
            sage: L1 is L2
            True
        """
        if isinstance(arg0, LieAlgebra):
            ct = arg0.cartan_type()
            if not ct.is_finite():
                raise ValueError("the base Lie algebra is not simple")
            cartan_type = ct.affine()
            g = arg0
        else:
            # arg0 is the base ring
            cartan_type = CartanType(cartan_type)
            if not cartan_type.is_affine():
                raise ValueError("the Cartan type must be affine")
            if cartan_type.is_untwisted_affine():
                g = LieAlgebra(arg0, cartan_type=cartan_type.classical())

        if cartan_type.is_untwisted_affine():
            return UntwistedAffineLieAlgebra(g, kac_moody=kac_moody)
        return TwistedAffineLieAlgebra(arg0, cartan_type, kac_moody=kac_moody)

    def __init__(self, g, cartan_type, names, kac_moody):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: asl = lie_algebras.Affine(QQ, ['D', 4, 1])
            sage: TestSuite(asl).run()
        """
        self._g = g
        self._cartan_type = cartan_type

        if kac_moody:
            names += ['d']
        self._kac_moody = kac_moody

        names = tuple(names)
        self._ordered_indices = names
        R = g.base_ring()
        cat = KacMoodyAlgebras(R).WithBasis()
        if not self._cartan_type.is_untwisted_affine():
            cat = cat.Subobjects()
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, names, category=cat)

    @cached_method
    def basis(self):
        r"""
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['D', 4, 1])
            sage: B = g.basis()
            sage: al = RootSystem(['D',4]).root_lattice().simple_roots()
            sage: B[al[1]+al[2]+al[4],4]
            (E[alpha[1] + alpha[2] + alpha[4]])#t^4
            sage: B[-al[1]-2*al[2]-al[3]-al[4],2]
            (E[-alpha[1] - 2*alpha[2] - alpha[3] - alpha[4]])#t^2
            sage: B[al[4],-2]
            (E[alpha[4]])#t^-2
            sage: B['c']
            c
            sage: B['d']
            d

            sage: g = LieAlgebra(QQ, cartan_type=['D', 4, 2], kac_moody=False)
            sage: B = g.basis()
            sage: it = iter(B)
            sage: [next(it) for _ in range(3)]
            [c, (E[alpha[1]])#t^0, (E[alpha[2]])#t^0]
            sage: B['c']
            c
            sage: B['d']
            0
        """
        if self._cartan_type.is_untwisted_affine():
            K = cartesian_product([self._g.basis().keys(), ZZ])
        else:
            K = TwistedAffineIndices(self._cartan_type)
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        c = FiniteEnumeratedSet(['c'])
        if self._kac_moody:
            d = FiniteEnumeratedSet(['d'])
            keys = DisjointUnionEnumeratedSets([c, d, K])
        else:
            keys = DisjointUnionEnumeratedSets([c, K])
        return Family(keys, self.monomial)

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A',1])
            sage: A = g.affine()
            sage: D = A.derived_subalgebra()
            sage: A(D.an_element())
            (E[alpha[1]] + h1 + E[-alpha[1]])#t^0
             + (E[-alpha[1]])#t^1 + (E[alpha[1]])#t^-1 + c
            sage: A(g.an_element())
            (E[alpha[1]] + h1 + E[-alpha[1]])#t^0
        """
        P = parent(x)
        if P is self.derived_subalgebra():
            return self.element_class(self, x.t_dict(), x.c_coefficient(),
                                      x.d_coefficient())
        if P == self._g:
            zero = self.base_ring().zero()
            return self.element_class(self, {0: x}, zero, zero)
        return super()._element_constructor_(x)

    def _coerce_map_from_(self, R):
        """
        Return the coerce map from ``R`` to ``self`` or ``True`` if
        a coerce map exists.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G',2])
            sage: A = g.affine()
            sage: A.has_coerce_map_from(g)
            True
            sage: D = A.derived_subalgebra()
            sage: A.has_coerce_map_from(D)
            True
        """
        if R is self.derived_subalgebra() or R is self._g:
            return True
        return super()._coerce_map_from_(R)

    def derived_series(self):
        """
        Return the derived series of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B',3,1])
            sage: g.derived_series()
            [Affine Kac-Moody algebra of ['B', 3] in the Chevalley basis,
             Affine Lie algebra of ['B', 3] in the Chevalley basis]
            sage: g.lower_central_series()
            [Affine Kac-Moody algebra of ['B', 3] in the Chevalley basis,
             Affine Lie algebra of ['B', 3] in the Chevalley basis]

            sage: D = g.derived_subalgebra()
            sage: D.derived_series()
            [Affine Lie algebra of ['B', 3] in the Chevalley basis]
        """
        if self._kac_moody:
            return [self, self.derived_subalgebra()]
        return [self]

    lower_central_series = derived_series

    def is_nilpotent(self):
        """
        Return ``False`` as ``self`` is semisimple.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B',3,1])
            sage: g.is_nilpotent()
            False
            sage: g.is_solvable()
            False
        """
        return False

    is_solvable = is_nilpotent

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['C',3,1])
            sage: g.cartan_type()
            ['C', 3, 1]
        """
        return self._cartan_type

    def classical(self):
        r"""
        Return the classical Lie algebra of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['F',4,1])
            sage: g.classical()
            Lie algebra of ['F', 4] in the Chevalley basis

            sage: so5 = lie_algebras.so(QQ, 5, 'matrix')
            sage: A = so5.affine()
            sage: A.classical() == so5
            True
        """
        return self._g

    @cached_method
    def zero(self):
        r"""
        Return the element `0`.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['F',4,1])
            sage: g.zero()
            0
        """
        zero = self.base_ring().zero()
        return self.element_class(self, {}, zero, zero)

    @cached_method
    def c(self):
        r"""
        Return the canonical central element `c` of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A',3,1])
            sage: g.c()
            c
        """
        R = self.base_ring()
        return self.element_class(self, {}, R.one(), R.zero())

    @cached_method
    def d(self):
        r"""
        Return the canonical derivation `d` of ``self``.

        If ``self`` is the affine Lie algebra, then this returns `0`.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A',3,1])
            sage: g.d()
            d
            sage: D = g.derived_subalgebra()
            sage: D.d()
            0
        """
        if not self._kac_moody:
            return self.zero()
        R = self.base_ring()
        return self.element_class(self, {}, R.zero(), R.one())

    @cached_method
    def lie_algebra_generators(self):
        r"""
        Return the Lie algebra generators of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A',1,1])
            sage: list(g.lie_algebra_generators())
            [(E[alpha[1]])#t^0,
             (E[-alpha[1]])#t^0,
             (h1)#t^0,
             (E[-alpha[1]])#t^1,
             (E[alpha[1]])#t^-1,
             c,
             d]

            sage: L = LieAlgebra(QQ, cartan_type=['A',5,2])
            sage: list(L.lie_algebra_generators())
            [(E[alpha[1]])#t^0,
             (E[alpha[2]])#t^0,
             (E[alpha[3]])#t^0,
             (E[-alpha[1]])#t^0,
             (E[-alpha[2]])#t^0,
             (E[-alpha[3]])#t^0,
             (h1)#t^0,
             (h2)#t^0,
             (h3)#t^0,
             (E[-alpha[1] - 2*alpha[2] - alpha[3]])#t^1,
             (E[alpha[1] + 2*alpha[2] + alpha[3]])#t^-1,
             c,
             d]
        """
        zero = self.base_ring().zero()
        d = {}
        if self._kac_moody:
            d['d'] = self.d()
        d['c'] = self.c()

        try:
            finite_gens = dict(self._g.lie_algebra_generators(True))
        except TypeError:
            finite_gens = dict(self._g.lie_algebra_generators())
        for k, g in finite_gens.items():
            d[k] = self.element_class(self, {0: g}, zero, zero)

        if self._cartan_type.is_untwisted_affine():
            # e_0 = f_{\theta} t
            d['e0'] = self.element_class(self, {1: self._g.highest_root_basis_elt(False)},
                                         zero, zero)
            # f_0 = e_{\theta} t^-1
            d['f0'] = self.element_class(self, {-1: self._g.highest_root_basis_elt(True)},
                                         zero, zero)
        elif self._cartan_type.type() != 'BC':
            a = self._cartan_type.a()
            Q = self._g._Q
            theta = Q._from_dict({i: a[i] for i in Q.index_set()}, remove_zeros=False)
            # e_0 = f_{\theta} t
            d['e0'] = self.element_class(self, {1: self._g.basis()[-theta]},
                                         zero, zero)
            # f_0 = e_{\theta} t^-1
            d['f0'] = self.element_class(self, {-1: self._g.basis()[theta]},
                                         zero, zero)
        else:
            n = self._g.cartan_type().rank()
            a = self._cartan_type.a()
            Q = self._g._Q
            theta = Q._from_dict({i: ZZ(2) for i in Q.index_set()}, remove_zeros=False)
            # e_0 = f_{\theta} t
            d[f'e{n}'] = self.element_class(self, {1: self._g1.basis()[-theta]},
                                         zero, zero)
            # f_0 = e_{\theta} t^-1
            d[f'f{n}'] = self.element_class(self, {-1: self._g1.basis()[theta]},
                                         zero, zero)

        return Family(self.variable_names(), d.__getitem__)

    def e(self, i=None):
        """
        Return the generators `e` of ``self``.

        INPUT:

        - ``i`` -- (optional) if specified, return just the
          generator `e_i`

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 3, 1])
            sage: list(g.e())
            [(E[-alpha[1] - 2*alpha[2] - 2*alpha[3]])#t^1,
             (E[alpha[1]])#t^0, (E[alpha[2]])#t^0, (E[alpha[3]])#t^0]
            sage: g.e(2)
            (E[alpha[2]])#t^0
        """
        gens = self.lie_algebra_generators()
        if i is None:
            I = self._cartan_type.index_set()
            d = {j: gens[f'e{j}'] for j in I}
            return Family(I, d.__getitem__)
        return gens[f'e{i}']

    def f(self, i=None):
        """
        Return the generators `f` of ``self``.

        INPUT:

        - ``i`` -- (optional) if specified, return just the
          generator `f_i`

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 5, 2])
            sage: list(g.f())
            [(E[alpha[1] + 2*alpha[2] + alpha[3]])#t^-1,
             (E[-alpha[1]])#t^0, (E[-alpha[2]])#t^0, (E[-alpha[3]])#t^0]
            sage: g.f(2)
            (E[-alpha[2]])#t^0
        """
        gens = self.lie_algebra_generators()
        if i is None:
            I = self._cartan_type.index_set()
            d = {j: gens[f'f{j}'] for j in I}
            return Family(I, d.__getitem__)
        return gens[f'f{i}']

    def monomial(self, m):
        r"""
        Construct the monomial indexed by ``m``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B',4,1])
            sage: al = RootSystem(['B',4]).root_lattice().simple_roots()
            sage: g.monomial((al[1]+al[2]+al[3],4))
            (E[alpha[1] + alpha[2] + alpha[3]])#t^4
            sage: g.monomial((-al[1]-al[2]-2*al[3]-2*al[4],2))
            (E[-alpha[1] - alpha[2] - 2*alpha[3] - 2*alpha[4]])#t^2
            sage: g.monomial((al[4],-2))
            (E[alpha[4]])#t^-2
            sage: g.monomial('c')
            c
            sage: g.monomial('d')
            d
        """
        if m == 'c':
            return self.c()
        if m == 'd':
            return self.d()
        G = self._g.basis()
        zero = self.base_ring().zero()
        return self.element_class(self, {m[1]: G[m[0]]}, zero, zero)


class UntwistedAffineLieAlgebra(AffineLieAlgebra):
    r"""
    An untwisted affine Lie algebra.

    Let `R` be a ring.  Given a finite-dimensional simple Lie algebra
    `\mathfrak{g}` over `R`, the affine Lie algebra
    `\widehat{\mathfrak{g}}^{\prime}` associated to `\mathfrak{g}` is
    defined as

    .. MATH::

        \widehat{\mathfrak{g}}' = \bigl( \mathfrak{g} \otimes
        R[t, t^{-1}] \bigr) \oplus R c,

    where `c` is the canonical central element and `R[t, t^{-1}]` is the
    Laurent polynomial ring over `R`. The Lie bracket is defined as

    .. MATH::

        [x \otimes t^m + \lambda c, y \otimes t^n + \mu c] =
        [x, y] \otimes t^{m+n} + m \delta_{m,-n} ( x | y ) c,

    where `( x | y )` is the Killing form on `\mathfrak{g}`.

    There is a canonical derivation `d` on `\widehat{\mathfrak{g}}'`
    that is defined by

    .. MATH::

        d(x \otimes t^m + \lambda c) = a \otimes m t^m,

    or equivalently by `d = t \frac{d}{dt}`.

    The affine Kac-Moody algebra `\widehat{\mathfrak{g}}` is formed by
    adjoining the derivation `d` such that

    .. MATH::

        \widehat{\mathfrak{g}} = \bigl( \mathfrak{g} \otimes R[t,t^{-1}]
        \bigr) \oplus R c \oplus R d.

    Specifically, the bracket on `\widehat{\mathfrak{g}}` is defined as

    .. MATH::

        [t^m \otimes x \oplus \lambda c \oplus \mu d, t^n \otimes y \oplus
        \lambda_1 c \oplus \mu_1 d] = \bigl( t^{m+n} [x,y] + \mu n t^n \otimes
        y - \mu_1 m t^m \otimes x\bigr) \oplus m \delta_{m,-n} (x|y) c .

    EXAMPLES:

    We begin by constructing an affine Kac-Moody algebra of type `G_2^{(1)}`
    from the classical Lie algebra of type `G_2`::

        sage: g = LieAlgebra(QQ, cartan_type=['G',2])
        sage: A = g.affine()
        sage: A
        Affine Kac-Moody algebra of ['G', 2] in the Chevalley basis

    Next, we construct the generators and perform some computations::

        sage: A.inject_variables()
        Defining e1, e2, f1, f2, h1, h2, e0, f0, c, d
        sage: e1.bracket(f1)
        (h1)#t^0
        sage: e0.bracket(f0)
        (-h1 - 2*h2)#t^0 + 8*c
        sage: e0.bracket(f1)
        0
        sage: A[d, f0]
        (-E[3*alpha[1] + 2*alpha[2]])#t^-1
        sage: A([[e0, e2], [[[e1, e2], [e0, [e1, e2]]], e1]])
        (-6*E[-3*alpha[1] - alpha[2]])#t^2
        sage: f0.bracket(f1)
        0
        sage: f0.bracket(f2)
        (E[3*alpha[1] + alpha[2]])#t^-1
        sage: A[h1+3*h2, A[[[f0, f2], f1], [f1,f2]] + f1] - f1
        (2*E[alpha[1]])#t^-1

    We can construct its derived subalgebra, the affine Lie algebra
    of type `G_2^{(1)}`. In this case, there is no canonical derivation,
    so the generator `d` is `0`::

        sage: D = A.derived_subalgebra()
        sage: D.d()
        0
    """
    def __init__(self, g, kac_moody):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: asl = lie_algebras.Affine(QQ, ['A',4,1])
            sage: TestSuite(asl).run()
        """
        cartan_type = g.cartan_type().affine()
        names = list(g.variable_names()) + ['e0', 'f0', 'c']
        super().__init__(g, cartan_type, names, kac_moody)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['D',4,1])
            sage: g
            Affine Kac-Moody algebra of ['D', 4] in the Chevalley basis
            sage: g.derived_subalgebra()
            Affine Lie algebra of ['D', 4] in the Chevalley basis
        """
        base = "Affine "
        rep = repr(self._g)
        if self._kac_moody:
            old_len = len(rep)
            rep = rep.replace("Lie", "Kac-Moody")
            if len(rep) == old_len: # We did not replace anything
                base += "Kac-Moody "
        return base + rep

    def derived_subalgebra(self):
        """
        Return the derived subalgebra of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B',3,1])
            sage: g
            Affine Kac-Moody algebra of ['B', 3] in the Chevalley basis
            sage: D = g.derived_subalgebra(); D
            Affine Lie algebra of ['B', 3] in the Chevalley basis
            sage: D.derived_subalgebra() == D
            True
        """
        if self._kac_moody:
            return UntwistedAffineLieAlgebra(self._g, kac_moody=False)
        return self

    Element = UntwistedAffineLieAlgebraElement


class TwistedAffineLieAlgebra(AffineLieAlgebra):
    r"""
    A twisted affine Lie algebra.

    A twisted affine Lie algebra is an affine Lie algebra for
    type `X_N^{(r)}` with `r > 1`. We realize this inside an
    untwisted affine Kac--Moody Lie algebra following Chapter 8
    of [Ka1990]_.

    Let `\overline{\mathfrak{g}}` be the classical Lie algebra by
    taking the index set `I \setminus \{\epsilon\}`, where
    `\epsilon = 0` unless `\epsilon = n` for `X_N^{(r)} = A_{2n}^{(2)}`,
    for the twisted affine Lie algebra `\widetilde{\mathfrak{g}}`.
    Let `\mathfrak{g}` be the basic Lie algebra of type `X_N`.
    We realize `\overline{\mathfrak{g}}` as the fixed-point subalgebra
    `\mathfrak{g}^{(0)}` of `\mathfrak{g}` under the order `r` diagram
    automorphism `\mu`. This naturally acts on the `\zeta_r` (a primitive
    `r`-th root of unity) eigenspace `\mathfrak{g}^{(1)}` of `\mu`,
    which is the highest weight representation corresponding to
    the small adjoint (where the weight spaces are the short roots
    of `\overline{\mathfrak{g}}`). The *twisted affine (Kac-Moody)
    Lie algebra* `\widehat{\mathfrak{g}}` is constructed as the
    subalgebra of `X_N^{(1)}` given by

    .. MATH::

        \sum_{i \in \ZZ} \mathfrak{g}^{(i \mod 2)} \otimes t^i
        \oplus R c \oplus R d,

    where `R` is the base ring.

    We encode our basis by using the classical Lie algebra except
    for type `A_{2n}^{(2)}`. For type `A_{2n}^{(2)}`, the fixed-point
    algebra `\mathfrak{g}^{(0)}` is of type `B_n` using the index set
    `\{0, \ldots, n-1\}`. For `\mathfrak{g}^{(1)}`, we identify the
    weights in this representation with the roots of type `B_n` and
    the double all of its short roots.
    """
    def __init__(self, R, cartan_type, kac_moody):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.Affine(QQ, ['A', 5, 2])
            sage: TestSuite(g).run()

            sage: g = lie_algebras.Affine(QQ, ['D', 4, 2])
            sage: TestSuite(g).run()

            sage: g = lie_algebras.Affine(QQ, ['D', 5, 2])
            sage: TestSuite(g).run()

            sage: g = lie_algebras.Affine(QQ, ['A', 6, 2])
            sage: TestSuite(g).run(skip=['_test_elements'])  # _test_monomial_coefficients fails

            sage: g = lie_algebras.Affine(QQ, ['A', 2, 2])
            sage: TestSuite(g).run(skip=['_test_elements'])  # _test_monomial_coefficients fails

            sage: g = lie_algebras.Affine(QQ, ['E', 6, 2])
            sage: TestSuite(g).run()  # long time

            sage: g = lie_algebras.Affine(QQ, ['D', 4, 3])
            sage: TestSuite(g).run()  # long time
        """
        # basic setup for AffineLieAlgebra
        if cartan_type.type() == 'BC':
            classical = cartan_type.classical().dual()
            n = classical.rank()
            classical = classical.relabel({n-i: i for i in range(n)})
        else:
            classical = cartan_type.classical()
        g = LieAlgebra(R, cartan_type=classical)
        n = classical.rank()
        names = ['e%s' % i for i in range(1, n+1)]
        names.extend('f%s' % i for i in range(1, n+1))
        if cartan_type.type() == 'BC':
            names.extend('h%s' % i for i in range(n))
        else:
            names.extend('h%s' % i for i in range(1, n+1))
        names += ['e0', 'f0', 'c']
        super().__init__(g, cartan_type, names, kac_moody)

        # setup the ambient simply-laced algebra
        basic_ct = cartan_type.basic_untwisted()
        if cartan_type.dual().type() == 'B':
            ep = [(i, i+1) for i in range(1, n)]
            ep.extend((i+1, i) for i in range(n, 2*n-1))
        elif cartan_type.dual().type() == 'F':
            ep = [(1, 3), (3, 4), (5, 4), (6, 5), (4, 2)]
        elif cartan_type.dual().type() == 'G':
            ep = [(1, 2), (3, 2), (4, 2)]
        else:
            ep = basic_ct.dynkin_diagram().to_undirected().edges(labels=False, sort=False)

        if self._cartan_type.dual().type() == 'G':
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            RP = PolynomialRing(R, 'x')
            Rext = RP.quotient(RP.gen(0)**3 - 1)
            self._basic = LieAlgebra(Rext, cartan_type=basic_ct, epsilon=ep)
        else:
            self._basic = LieAlgebra(R, cartan_type=basic_ct, epsilon=ep)
        self._ambient = self._basic.affine(kac_moody=self._kac_moody)

        # setup the embeddings
        basic_ct = self._basic.cartan_type()
        if self._cartan_type.dual().type() == 'G':
            gens = basic_ct.dynkin_diagram().automorphism_group().gens()
            auto = gens[0] * gens[1]
        else:
            auto = basic_ct.dynkin_diagram().automorphism_group().gen(0)
        basic_Q = basic_ct.root_system().root_lattice()
        basic_p_roots = basic_Q.positive_roots()
        visited = set()
        orbits = []
        for al in basic_p_roots:
            if al in visited:
                continue
            visited.add(al)
            O = [al]
            cur = basic_Q._from_dict({auto(i): c for i, c in al}, remove_zeros=False)
            while cur != al:
                O.append(cur)
                visited.add(cur)
                cur = basic_Q._from_dict({auto(i): c for i, c in cur}, remove_zeros=False)
            orbits.append(O)

        finite_ct = self._g.cartan_type()
        Q = finite_ct.root_system().root_lattice()
        I = finite_ct.index_set()
        a = finite_ct.symmetrizer()
        ord = auto.order()

        if self._cartan_type.dual().type() == 'F':
            # For E_6^(2), we need to take into account the lack of index sets matching
            reindex = {2: 4, 4: 3, 3: 2, 1: 1}

            def build_root(O):
                return Q._from_dict({reindex[i]: c * (ord // a[reindex[i]]) / len(O) for i, c in sum(O) if i in reindex},
                                    remove_zeros=False)
        elif self._cartan_type.type() == 'BC':
            reindex = {n-i: i for i in range(finite_ct.rank())}

            def build_root(O):
                return Q._from_dict({reindex[i]: c * (ord // len(O)) for i, c in sum(O) if i in reindex},
                                    remove_zeros=False)
        else:

            def build_root(O):
                return Q._from_dict({i: c * (ord // a[i]) / len(O) for i, c in sum(O) if i in I},
                                    remove_zeros=False)

        self._root_mapping = {build_root(O): O for O in orbits}
        for r in list(self._root_mapping.keys()):
            self._root_mapping[-r] = [-s for s in self._root_mapping[r]]
        if self._cartan_type.type() == 'BC':
            assert {r for r in self._root_mapping if len(self._root_mapping[r]) > 1} == set(Q.roots())
            if self._cartan_type.rank() == 2:
                # Special case since sl_2 has only 1 root length
                assert {r / 2 for r in self._root_mapping if len(self._root_mapping[r]) == 1} == set(Q.roots())
            else:
                assert {r / 2 for r in self._root_mapping if len(self._root_mapping[r]) == 1} == set(Q.short_roots())
            from sage.combinat.free_module import CombinatorialFreeModule
            X = sorted(self._root_mapping, key=str)
            self._g1 = CombinatorialFreeModule(R, X, prefix='E')
        else:
            assert set(self._root_mapping) == set(Q.roots())
        al = Q.simple_roots()
        ac = Q.simple_coroots()
        for i in I:
            self._root_mapping[ac[i]] = [r.associated_coroot() for r in self._root_mapping[al[i]]]
        self._inverse_root_map = {O[0]: r for r, O in self._root_mapping.items()}

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['D', 4, 2])
            sage: g
            Twisted affine Kac-Moody algebra of type ['C', 3, 1]^* over Rational Field
            sage: g.derived_subalgebra()
            Twisted affine Lie algebra of type ['C', 3, 1]^* over Rational Field
        """
        rep = "Twisted affine "
        rep += "Kac-Moody " if self._kac_moody else "Lie "
        rep += f"algebra of type {self._cartan_type} over {self.base_ring()}"
        return rep

    def _test_classical_subalgebra(self, **options):
        r"""
        Test the Chevalley basis properties for the classical subalgebra
        of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 5, 2])
            sage: L._test_classical_subalgebra()
            sage: L = LieAlgebra(QQ, cartan_type=['D', 4, 2])
            sage: L._test_classical_subalgebra()
        """
        tester = self._tester(**options)
        B = self.basis()
        roots = set(self._g._Q.roots())
        ac = list(self._g._Q.simple_coroots())
        from sage.misc.misc import some_tuples
        for r, s in some_tuples(roots, 2, tester._max_runs):
            ret = B[r,0].bracket(B[s,0])
            if r + s in roots:
                tester.assertEqual(list(ret.support()), [(r+s, 0)], f"obtained [{r}, {s}] == {ret}")
            elif r == -s:
                supp = {(ac, 0) for ac in r.associated_coroot().monomials()}
                tester.assertEqual(set(ret.support()), supp, f"obtained [{r}, {s}] == {ret}")
            else:
                tester.assertEqual(ret, self.zero(), f"nonzero for [{r}, {s}]")

    def derived_subalgebra(self):
        r"""
        Return the derived subalgebra of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 5, 2])
            sage: g
            Twisted affine Kac-Moody algebra of type ['B', 3, 1]^* over Rational Field
            sage: D = g.derived_subalgebra(); D
            Twisted affine Lie algebra of type ['B', 3, 1]^* over Rational Field
            sage: D.derived_subalgebra() == D
            True
        """
        if self._kac_moody:
            return TwistedAffineLieAlgebra(self.base_ring(), self._cartan_type, kac_moody=False)
        return self

    def ambient(self):
        r"""
        Return the ambient untwisted affine Lie algebra of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 5, 2])
            sage: g.ambient()
            Affine Kac-Moody algebra of ['A', 5] in the Chevalley basis
        """
        return self._ambient

    def retract(self, x):
        r"""
        Retract the element ``x`` from the ambient untwisted affine Lie
        algebra into ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 5, 2])
            sage: it = iter(g.basis())
            sage: elts = [next(it) for _ in range(20)]
            sage: elts
            [c,
             d,
             (E[alpha[1]])#t^0,
             (E[alpha[2]])#t^0,
             (E[alpha[3]])#t^0,
             (E[alpha[1] + alpha[2]])#t^0,
             (E[alpha[2] + alpha[3]])#t^0,
             (E[2*alpha[2] + alpha[3]])#t^0,
             (E[alpha[1] + alpha[2] + alpha[3]])#t^0,
             (E[2*alpha[1] + 2*alpha[2] + alpha[3]])#t^0,
             (E[alpha[1] + 2*alpha[2] + alpha[3]])#t^0,
             (E[-alpha[1]])#t^0,
             (E[-alpha[2]])#t^0,
             (E[-alpha[3]])#t^0,
             (E[-alpha[1] - alpha[2]])#t^0,
             (E[-alpha[2] - alpha[3]])#t^0,
             (E[-2*alpha[2] - alpha[3]])#t^0,
             (E[-alpha[1] - alpha[2] - alpha[3]])#t^0,
             (E[-2*alpha[1] - 2*alpha[2] - alpha[3]])#t^0,
             (E[-alpha[1] - 2*alpha[2] - alpha[3]])#t^0]
            sage: all(g.retract(g.to_ambient(x)) == x for x in elts)
            True
        """
        t_dict = x.t_dict()
        c_coeff = x.c_coefficient()
        d_coeff = x.d_coefficient()
        if self._cartan_type.dual().type() == 'G':
            R = self.base_ring()
            for i in t_dict:
                t_dict[i] = self._g._from_dict({self._inverse_root_map[r]: R(c.lift())
                                                for r, c in t_dict[i] if r in self._inverse_root_map},
                                               remove_zeros=False)
        elif self._cartan_type.type() == 'BC':
            for i in t_dict:
                if i % 2:
                    t_dict[i] = self._g1._from_dict({self._inverse_root_map[r]: c for r, c in t_dict[i]
                                                     if r in self._inverse_root_map},
                                                    remove_zeros=False)
                else:
                    t_dict[i] = self._g._from_dict({self._inverse_root_map[r]: c for r, c in t_dict[i]
                                                    if r in self._inverse_root_map},
                                                   remove_zeros=False)
        else:
            for i in t_dict:
                t_dict[i] = self._g._from_dict({self._inverse_root_map[r]: c for r, c in t_dict[i]
                                                if r in self._inverse_root_map},
                                               remove_zeros=False)
        return self.element_class(self, t_dict, c_coeff, d_coeff)

    @lazy_attribute
    def to_ambient(self):
        r"""
        Lift the element ``x`` from the ambient untwisted affine Lie
        algebra into ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 5, 2])
            sage: g.to_ambient
            Generic morphism:
              From: Twisted affine Kac-Moody algebra of type ['B', 3, 1]^* over Rational Field
              To:   Affine Kac-Moody algebra of ['A', 5] in the Chevalley basis
        """
        one = self.base_ring().one()

        if self._cartan_type.type() == 'BC':
            mone = -one

            def basis_map(r):
                O = self._root_mapping[r]
                return self._basic._from_dict({O[0]: one, O[1]: mone**(1+O[1].height())},
                                              remove_zeros=False)
        else:

            def basis_map(r):
                return self._basic._from_dict({s: one for s in self._root_mapping[r]}, remove_zeros=False)

        if self._cartan_type.dual().type() == 'G':
            zeta3 = self._basic.base_ring().gen()

            def basis_alt(r):
                return self._basic._from_dict({s: zeta3**ind for ind, s in enumerate(self._root_mapping[r])},
                                              remove_zeros=False)
        elif self._cartan_type.type() == 'BC':

            def basis_alt(r):
                O = self._root_mapping[r]
                if len(O) == 1:
                    return self._basic.monomial(O[0])
                return self._basic._from_dict({O[0]: one, O[1]: mone**O[1].height()},
                                              remove_zeros=False)
        else:
            mone = -one

            def basis_alt(r):
                return self._basic._from_dict({s: mone**ind for ind, s in enumerate(self._root_mapping[r])},
                                              remove_zeros=False)

        def lift_map(elt):
            t_dict = elt.t_dict()
            c_coeff = elt.c_coefficient()
            d_coeff = elt.d_coefficient()
            for i in t_dict:
                if i % 2:
                    t_dict[i] = self._basic.linear_combination((basis_alt(r), c)
                                                               for r, c in t_dict[i])
                else:
                    t_dict[i] = self._basic.linear_combination((basis_map(r), c)
                                                               for r, c in t_dict[i])
            return self._ambient.element_class(self._ambient, t_dict, c_coeff, d_coeff)

        return self.module_morphism(function=lift_map, codomain=self._ambient)

    class Element(UntwistedAffineLieAlgebraElement):
        def _bracket_(self, y):
            r"""
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::

                sage: g = LieAlgebra(QQ, cartan_type=['D', 5, 2])
                sage: e0, e1, e2, e3, e4 = g.e()
                sage: f0, f1, f2, f3, f4 = g.f()
                sage: B = g.basis()
                sage: Q = g.classical().cartan_type().root_system().root_lattice()
                sage: h1, h2, h3, h4 = [B[ac, 0] for ac in Q.simple_coroots()]
                sage: e1._bracket_(e2)
                (-E[alpha[1] + alpha[2]])#t^0
                sage: e1._bracket_(e3)
                0
                sage: e0._bracket_(e1)
                (-E[-alpha[2] - alpha[3] - alpha[4]])#t^1
                sage: e0._bracket_(e2)
                0
                sage: f1._bracket_(f2)
                (E[-alpha[1] - alpha[2]])#t^0
                sage: f1._bracket_(f3)
                0
                sage: f0._bracket_(f1)
                (E[alpha[2] + alpha[3] + alpha[4]])#t^-1
                sage: f0._bracket_(f2)
                0
                sage: g[f0, e0]
                (2*h1 + 2*h2 + 2*h3 + h4)#t^0 + -32*c
                sage: g([f1, [e1, e2]])
                (E[alpha[2]])#t^0
                sage: g[h1, e0]
                (-E[-alpha[1] - alpha[2] - alpha[3] - alpha[4]])#t^1
                sage: g[h2, f0]
                0
            """
            P = parent(self)
            ax = P.to_ambient(self)
            ay = P.to_ambient(y)
            return P.retract(ax.bracket(ay))


class TwistedAffineIndices(UniqueRepresentation, Set_generic):
    r"""
    The indices for the basis of a twisted affine Lie algebra.

    INPUT:

    - ``cartan_type`` -- the Cartan type of twisted affine type Lie algebra

    EXAMPLES::

        sage: from sage.algebras.lie_algebras.affine_lie_algebra import TwistedAffineIndices
        sage: I = TwistedAffineIndices(['A', 3, 2])
        sage: it = iter(I)
        sage: [next(it) for _ in range(20)]
        [(alpha[1], 0), (alpha[2], 0), (alpha[1] + alpha[2], 0),
         (2*alpha[1] + alpha[2], 0), (-alpha[1], 0), (-alpha[2], 0),
         (-alpha[1] - alpha[2], 0), (-2*alpha[1] - alpha[2], 0),
         (alphacheck[1], 0), (alphacheck[2], 0), (alpha[1], 1),
         (alpha[1] + alpha[2], 1), (-alpha[1], 1), (-alpha[1] - alpha[2], 1),
         (alphacheck[1], 1), (alpha[1], -1), (alpha[1] + alpha[2], -1),
         (-alpha[1], -1), (-alpha[1] - alpha[2], -1), (alphacheck[1], -1)]

        sage: I = TwistedAffineIndices(['A', 4, 2])
        sage: it = iter(I)
        sage: [next(it) for _ in range(20)]
        [(alpha[0], 0), (alpha[1], 0), (alpha[0] + alpha[1], 0),
         (2*alpha[0] + alpha[1], 0), (-alpha[0], 0), (-alpha[1], 0),
         (-alpha[0] - alpha[1], 0), (-2*alpha[0] - alpha[1], 0),
         (alphacheck[0], 0), (alphacheck[1], 0), (alpha[0], 1), (alpha[1], 1),
         (alpha[0] + alpha[1], 1), (2*alpha[0] + alpha[1], 1), (-alpha[0], 1),
         (-alpha[1], 1), (-alpha[0] - alpha[1], 1), (-2*alpha[0] - alpha[1], 1),
         (2*alpha[0], 1), (2*alpha[0] + 2*alpha[1], 1)]

        sage: I = TwistedAffineIndices(['A', 2, 2])
        sage: it = iter(I)
        sage: [next(it) for _ in range(10)]
        [(alpha[0], 0), (-alpha[0], 0), (alphacheck[0], 0), (alpha[0], 1),
         (-alpha[0], 1), (2*alpha[0], 1), (-2*alpha[0], 1),
         (alphacheck[0], 1), (alpha[0], -1), (-alpha[0], -1)]
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

            sage: from sage.algebras.lie_algebras.affine_lie_algebra import TwistedAffineIndices
            sage: I1 = TwistedAffineIndices(CartanType(['C', 4, 1]).dual())
            sage: I2 = TwistedAffineIndices(['D', 5, 2])
            sage: I1 is I2
            True
            sage: I = TwistedAffineIndices(['C', 4, 1])
            Traceback (most recent call last):
            ...
            ValueError: the Cartan type must be a twisted affine type
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_affine() or cartan_type.is_untwisted_affine():
            raise ValueError("the Cartan type must be a twisted affine type")
        return super().__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.affine_lie_algebra import TwistedAffineIndices
            sage: I = TwistedAffineIndices(['D', 4, 2])
            sage: TestSuite(I).run()
        """
        self._cartan_type = cartan_type
        if cartan_type.type() == 'BC':
            finite_ct = cartan_type.classical().dual()
            n = finite_ct.rank()
            Q = finite_ct.relabel({n-i: i for i in range(n)}).root_system().root_lattice()
            self._roots = tuple(Q.roots())
            self._ac = tuple(Q.simple_coroots())
            CP = cartesian_product([range(3)] * n)
            if cartan_type.rank() == 2:
                self._short_roots = self._roots + tuple(2*r for r in Q.roots())
            else:
                self._short_roots = self._roots + tuple(2*r for r in Q.short_roots())
            self._short_roots += self._ac
            facade = cartesian_product([self._short_roots, ZZ])
        else:
            Q = cartan_type.classical().root_system().root_lattice()
            self._roots = tuple(Q.roots())
            self._ac = tuple(Q.simple_coroots())
            self._short_roots = tuple(Q.short_roots())
            ac = Q.simple_coroots()
            self._short_roots += tuple([ac[i] for i in Q.index_set() if Q.simple_root(i).is_short_root()])
            facade = cartesian_product([self._roots + self._ac, ZZ])
        from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
        super().__init__(facade=facade, category=InfiniteEnumeratedSets())

    def __contains__(self, x):
        """
        Return if ``x`` is contained in ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.affine_lie_algebra import TwistedAffineIndices
            sage: I = TwistedAffineIndices(['D', 4, 2])
            sage: Q = RootSystem(['B', 3]).root_lattice()
            sage: all((r, 4) in I for r in Q.roots())
            True
            sage: all((r, 3) in I for r in Q.short_roots())
            True
            sage: all((r, 3) not in I for r in Q.long_roots())
            True
            sage: list(I.an_element()) in I  # lists are not included
            False
            sage: (5, Q) in I
            False
            sage: (5, 5) in I
            False
            sage: (Q.simple_root(1), Q.simple_root(1)) in I
            False
            sage: (Q.simple_coroot(2), 1) in I
            False
            sage: (Q.simple_coroot(3), 1) in I
            True
        """
        if x not in self._facade_for[0]:
            return False
        x = self._facade_for[0](x)
        # self._short_roots also contains the corresponding coroots
        return (x[1] % 2 == 0) or x[0] in self._short_roots

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.affine_lie_algebra import TwistedAffineIndices
            sage: I = TwistedAffineIndices(['D', 3, 2])
            sage: it = iter(I)
            sage: [next(it) for _ in range(22)]
            [(alpha[1], 0), (alpha[2], 0), (alpha[1] + 2*alpha[2], 0), (alpha[1] + alpha[2], 0),
             (-alpha[1], 0), (-alpha[2], 0), (-alpha[1] - 2*alpha[2], 0), (-alpha[1] - alpha[2], 0),
             (alphacheck[1], 0), (alphacheck[2], 0), (alpha[2], 1), (alpha[1] + alpha[2], 1),
             (-alpha[2], 1), (-alpha[1] - alpha[2], 1), (alphacheck[2], 1), (alpha[2], -1),
             (alpha[1] + alpha[2], -1), (-alpha[2], -1), (-alpha[1] - alpha[2], -1),
             (alphacheck[2], -1), (alpha[1], 2), (alpha[2], 2)]
        """
        if self._cartan_type.type() == 'BC':
            finite_ct = self._cartan_type.classical().dual()
            n = finite_ct.rank()
            finite_ct = finite_ct.relabel({n-i: i for i in range(n)})
        else:
            finite_ct = self._cartan_type.classical()
        Q = finite_ct.root_system().root_lattice()
        P = self._facade_for[0]
        for i in ZZ:
            if i % 2:
                # self._short_roots also contains the corresponding coroots
                for r in self._short_roots:
                    yield P((r, i))
            else:
                for r in self._roots:
                    yield P((r, i))
                for r in self._ac:
                    yield P((r, i))
