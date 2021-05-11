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
from sage.structure.element import parent
from sage.categories.lie_algebras import LieAlgebras

from sage.algebras.lie_algebras.lie_algebra import LieAlgebra, FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import UntwistedAffineLieAlgebraElement
from sage.combinat.root_system.cartan_type import CartanType
from sage.categories.cartesian_product import cartesian_product
from sage.rings.integer_ring import ZZ
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family


class AffineLieAlgebra(FinitelyGeneratedLieAlgebra):
    r"""
    An (untwisted) affine Lie algebra.

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

            sage: L1 = lie_algebras.Affine(QQ, ['A',4,1])
            sage: cl = lie_algebras.sl(QQ, 5)
            sage: L2 = lie_algebras.Affine(cl)
            sage: L1 is L2
            True
            sage: cl.affine() is L1
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
            g = LieAlgebra(arg0, cartan_type=cartan_type.classical())

        if not cartan_type.is_untwisted_affine():
            raise NotImplementedError("only currently implemented for untwisted affine types")
        return super(AffineLieAlgebra, cls).__classcall__(cls, g, kac_moody)

    def __init__(self, g, kac_moody):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: asl = lie_algebras.Affine(QQ, ['A',4,1])
            sage: TestSuite(asl).run()
        """
        self._g = g
        self._cartan_type = g.cartan_type().affine()
        R = g.base_ring()
        names = list(g.variable_names()) + ['e0', 'f0', 'c']

        if kac_moody:
            names += ['d']
        self._kac_moody = kac_moody

        names = tuple(names)
        self._ordered_indices = names
        cat = LieAlgebras(R).WithBasis()
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, names, category=cat)

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

    @cached_method
    def basis(self):
        r"""
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['D',4,1])
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
        """
        K = cartesian_product([self._g.basis().keys(), ZZ])
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
        return super(AffineLieAlgebra, self)._element_constructor_(x)

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
        return super(AffineLieAlgebra, self)._coerce_map_from_(R)

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
            return AffineLieAlgebra(self._g, kac_moody=False)
        return self

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
        for k,g in finite_gens.items():
            d[k] = self.element_class(self, {0: g}, zero, zero)
        # e_0 = f_{\theta} t
        d['e0'] = self.element_class(self, {1: self._g.highest_root_basis_elt(False)},
                                     zero, zero)
        # f_0 = e_{\theta} t^-1
        d['f0'] = self.element_class(self, {-1: self._g.highest_root_basis_elt(True)},
                                     zero, zero)
        return Family(self.variable_names(), d.__getitem__)

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

    Element = UntwistedAffineLieAlgebraElement

