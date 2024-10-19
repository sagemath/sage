r"""
Verma Modules

AUTHORS:

- Travis Scrimshaw (2017-06-30): Initial version

.. TODO::

    Implement a :class:`sage.categories.pushout.ConstructionFunctor`
    and return as the ``construction()``.
"""

# ****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.categories.modules import Modules
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom, Homset
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.combinat.free_module import CombinatorialFreeModule
from sage.modules.free_module_element import vector
from sage.sets.family import Family
from sage.structure.richcmp import richcmp
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


class ModulePrinting:
    """
    Helper mixin class for printing the module vectors.
    """
    def __init__(self, vector_name='v'):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.verma_module import ModulePrinting
            sage: MP = ModulePrinting()
            sage: TestSuite(MP).run(skip="_test_pickling")
        """
        self.__vector_name = vector_name

    def _repr_generator(self, m):
        r"""
        Return a string representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 4)
            sage: La = L.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: M = L.verma_module(-1/2*La[1] + 3/7*La[2])
            sage: f1, f2 = L.f()
            sage: x = M.pbw_basis()(L([f1, [f1, f2]]))
            sage: v = x * M.highest_weight_vector()
            sage: M._repr_generator(v.leading_support())
            'f[-2*alpha[1] - alpha[2]]*v[(-1/14, 3/7)]'

            sage: M.highest_weight_vector()
            v[(-1/14, 3/7)]
            sage: 2 * M.highest_weight_vector()
            2*v[(-1/14, 3/7)]
        """
        ret = super()._repr_generator(m)
        if ret == '1':
            ret = ''
        else:
            ret += '*'
        return ret + self.__vector_name + "[{}]".format(self._weight)

    def _latex_generator(self, m):
        r"""
        Return a latex representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 4)
            sage: La = L.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: M = L.verma_module(-1/2*La[1] + 3/7*La[2])
            sage: f1, f2 = L.f()
            sage: x = M.pbw_basis()(L([f1, [f1, f2]]))
            sage: v = x * M.highest_weight_vector()
            sage: M._latex_generator(v.leading_support())
            f_{-2 \alpha_{1} - \alpha_{2}} v_{-\frac{1}{14} e_{0} + \frac{3}{7} e_{1}}

            sage: latex(2 * M.highest_weight_vector())
            2  v_{-\frac{1}{14} e_{0} + \frac{3}{7} e_{1}}
            sage: latex(M.highest_weight_vector())
             v_{-\frac{1}{14} e_{0} + \frac{3}{7} e_{1}}
        """
        ret = super()._latex_generator(m)
        if ret == '1':
            ret = ''
        from sage.misc.latex import latex
        return ret + " {}_{{{}}}".format(self.__vector_name, latex(self._weight))

    _repr_term = _repr_generator
    _latex_term = _latex_generator


class VermaModule(ModulePrinting, CombinatorialFreeModule):
    r"""
    A Verma module.

    Let `\lambda` be a weight and `\mathfrak{g}` be a Kac--Moody Lie
    algebra with a fixed Borel subalgebra `\mathfrak{b} = \mathfrak{h}
    \oplus \mathfrak{g}^+`. The *Verma module* `M_{\lambda}` is a
    `U(\mathfrak{g})`-module given by

    .. MATH::

        M_{\lambda} := U(\mathfrak{g}) \otimes_{U(\mathfrak{b})} F_{\lambda},

    where `F_{\lambda}` is the `U(\mathfrak{b})` module such that
    `h \in U(\mathfrak{h})` acts as multiplication by
    `\langle \lambda, h \rangle` and `U(\mathfrak{g}^+) F_{\lambda} = 0`.

    INPUT:

    - ``g`` -- a Lie algebra
    - ``weight`` -- a weight

    EXAMPLES::

        sage: L = lie_algebras.sl(QQ, 3)
        sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
        sage: M = L.verma_module(2*La[1] + 3*La[2])
        sage: pbw = M.pbw_basis()
        sage: E1,E2,F1,F2,H1,H2 = [pbw(g) for g in L.gens()]
        sage: v = M.highest_weight_vector()
        sage: x = F2^3 * F1 * v
        sage: x
        f[-alpha[2]]^3*f[-alpha[1]]*v[2*Lambda[1] + 3*Lambda[2]]
        sage: F1 * x
        f[-alpha[2]]^3*f[-alpha[1]]^2*v[2*Lambda[1] + 3*Lambda[2]]
         + 3*f[-alpha[2]]^2*f[-alpha[1]]*f[-alpha[1] - alpha[2]]*v[2*Lambda[1] + 3*Lambda[2]]
        sage: E1 * x
        2*f[-alpha[2]]^3*v[2*Lambda[1] + 3*Lambda[2]]
        sage: H1 * x
        3*f[-alpha[2]]^3*f[-alpha[1]]*v[2*Lambda[1] + 3*Lambda[2]]
        sage: H2 * x
        -2*f[-alpha[2]]^3*f[-alpha[1]]*v[2*Lambda[1] + 3*Lambda[2]]

    REFERENCES:

    - :wikipedia:`Verma_module`
    """
    def __init__(self, g, weight, basis_key=None, prefix='f', **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + 4*La[2])
            sage: TestSuite(M).run()
            sage: M = L.verma_module(La[1] - 2*La[2])
            sage: TestSuite(M).run()

            sage: L = lie_algebras.sp(QQ, 4)
            sage: La = L.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: M = L.verma_module(-1/2*La[1] + 3/7*La[2])
            sage: TestSuite(M).run()
        """
        if basis_key is not None:
            self._basis_key = basis_key
        else:
            self._basis_key = g._basis_key

        self._weight = weight

        R = g.base_ring()
        self._g = g
        self._pbw = g.pbw_basis(basis_key=self._triangular_key)
        monomials = IndexedFreeAbelianMonoid(g._negative_half_index_set(),
                                             prefix,
                                             sorting_key=self._monoid_key,
                                             **kwds)
        CombinatorialFreeModule.__init__(self, R, monomials,
                                         prefix='', bracket=False, latex_bracket=False,
                                         sorting_key=self._monomial_key,
                                         category=Modules(R).WithBasis().Graded())
        ModulePrinting.__init__(self)

    def _triangular_key(self, x):
        r"""
        Return a key for sorting for the index ``x`` that respects
        the triangular decomposition by `U^-, U^0, U^+`.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1])
            sage: sorted(L.basis().keys(), key=L._basis_key)
            [alpha[2], alpha[1], alpha[1] + alpha[2],
             alphacheck[1], alphacheck[2],
             -alpha[2], -alpha[1], -alpha[1] - alpha[2]]
            sage: sorted(L.basis().keys(), key=M._triangular_key)
            [-alpha[2], -alpha[1], -alpha[1] - alpha[2],
             alphacheck[1], alphacheck[2],
             alpha[2], alpha[1], alpha[1] + alpha[2]]

            sage: def neg_key(x):
            ....:     return -L.basis().keys().index(x)
            sage: sorted(L.basis().keys(), key=neg_key)
            [-alpha[1] - alpha[2], -alpha[1], -alpha[2],
             alphacheck[2], alphacheck[1],
             alpha[1] + alpha[2], alpha[1], alpha[2]]
            sage: N = L.verma_module(La[1], basis_key=neg_key)
            sage: sorted(L.basis().keys(), key=N._triangular_key)
            [-alpha[1] - alpha[2], -alpha[1], -alpha[2],
             alphacheck[2], alphacheck[1],
             alpha[1] + alpha[2], alpha[1], alpha[2]]
        """
        return (self._g._part_on_basis(x), self._basis_key(x))

    def _monoid_key(self, x):
        """
        Return a key for comparison in the underlying monoid of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1])
            sage: monoid = M.basis().keys()
            sage: prod(monoid.gens())  # indirect doctest
            f[-alpha[2]]*f[-alpha[1]]*f[-alpha[1] - alpha[2]]
            sage: [M._monoid_key(x) for x in monoid.an_element()._sorted_items()]
            [5, 6, 7]

            sage: def neg_key(x):
            ....:     return -L.basis().keys().index(x)
            sage: M = L.verma_module(La[1], basis_key=neg_key)
            sage: monoid = M.basis().keys()
            sage: prod(monoid.gens())  # indirect doctest
            f[-alpha[1] - alpha[2]]*f[-alpha[1]]*f[-alpha[2]]
            sage: [M._monoid_key(x) for x in monoid.an_element()._sorted_items()]
            [-7, -6, -5]
        """
        return self._basis_key(x[0])

    def _monomial_key(self, x):
        """
        Compute the key for ``x`` so that the comparison is done by
        triangular decomposition and then reverse degree lexicographic order.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1])
            sage: pbw = M.pbw_basis()
            sage: f1,f2 = pbw(L.f(1)), pbw(L.f(2))
            sage: f1 * f2 * f1 * M.highest_weight_vector()  # indirect doctest
            f[-alpha[2]]*f[-alpha[1]]^2*v[Lambda[1]]
             + f[-alpha[1]]*f[-alpha[1] - alpha[2]]*v[Lambda[1]]

            sage: def neg_key(x):
            ....:     return -L.basis().keys().index(x)
            sage: M = L.verma_module(La[1], basis_key=neg_key)
            sage: f1 * f2 * f1 * M.highest_weight_vector()  # indirect doctest
            f[-alpha[1]]^2*f[-alpha[2]]*v[Lambda[1]]
             - f[-alpha[1] - alpha[2]]*f[-alpha[1]]*v[Lambda[1]]
        """
        return (-len(x), [self._triangular_key(l) for l in x.to_word_list()])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['E',6])
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(2*La[1] + 3*La[2] - 5*La[5])
            sage: M
            Verma module with highest weight 2*Lambda[1] + 3*Lambda[2] - 5*Lambda[5]
             of Lie algebra of ['E', 6] in the Chevalley basis
        """
        return "Verma module with highest weight {} of {}".format(self._weight, self._g)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['E',7])
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = L.verma_module(2*La[1] + 7*La[4] - 3/4*La[7])
            sage: latex(M)
            M_{2 \Lambda_{1} + 7 \Lambda_{4} - \frac{3}{4} \Lambda_{7}}
        """
        from sage.misc.latex import latex
        return "M_{{{}}}".format(latex(self._weight))

    def lie_algebra(self):
        r"""
        Return the underlying Lie algebra of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.so(QQ, 9)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = L.verma_module(La[3] - 1/2*La[1])
            sage: M.lie_algebra()
            Lie algebra of ['B', 4] in the Chevalley basis
        """
        return self._g

    def pbw_basis(self):
        r"""
        Return the PBW basis of the underlying Lie algebra
        used to define ``self``.

        EXAMPLES::

            sage: L = lie_algebras.so(QQ, 8)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[2] - 2*La[3])
            sage: M.pbw_basis()
            Universal enveloping algebra of Lie algebra of ['D', 4] in the Chevalley basis
             in the Poincare-Birkhoff-Witt basis
        """
        return self._pbw

    poincare_birkhoff_witt_basis = pbw_basis

    @cached_method
    def highest_weight_vector(self):
        """
        Return the highest weight vector of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] - 3*La[2])
            sage: M.highest_weight_vector()
            v[Lambda[1] - 3*Lambda[2]]
        """
        one = self.base_ring().one()
        return self._from_dict({self._indices.one(): one},
                               remove_zeros=False, coerce=False)

    def gens(self):
        r"""
        Return the generators of ``self`` as a `U(\mathfrak{g})`-module.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] - 3*La[2])
            sage: M.gens()
            (v[Lambda[1] - 3*Lambda[2]],)
        """
        return (self.highest_weight_vector(),)

    def highest_weight(self):
        r"""
        Return the highest weight of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.so(QQ, 7)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = L.verma_module(4*La[1] - 3/2*La[2])
            sage: M.highest_weight()
            4*Lambda[1] - 3/2*Lambda[2]
        """
        return self._weight

    def dual(self):
        r"""
        Return the dual module `M(\lambda)^{\vee}` in Category `\mathcal{O}`.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = L.verma_module(2*La[1])
            sage: Mc = M.dual()

            sage: Mp = L.verma_module(-2*La[1])
            sage: Mp.dual() is Mp
            True
        """
        if self.is_simple():
            return self
        from sage.algebras.lie_algebras.bgg_dual_module import BGGDualModule
        return BGGDualModule(self)

    def degree_on_basis(self, m):
        r"""
        Return the degree (or weight) of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(2*La[1] + 3*La[2])
            sage: v = M.highest_weight_vector()
            sage: M.degree_on_basis(v.leading_support())
            2*Lambda[1] + 3*Lambda[2]

            sage: pbw = M.pbw_basis()
            sage: G = list(pbw.gens())
            sage: f1, f2 = L.f()
            sage: x = pbw(f1.bracket(f2)) * pbw(f1) * v
            sage: x.degree()
            -Lambda[1] + 3*Lambda[2]
        """
        P = self._weight.parent()
        return self._weight + P.sum(P(e * self._g.degree_on_basis(k))
                                    for k,e in m.dict().items())

    def _coerce_map_from_(self, R):
        r"""
        Return if there is a coercion map from ``R`` to ``self``.

        There is a coercion map from ``R`` if and only if

        - there is a coercion from ``R`` into the base ring;
        - ``R`` is a Verma module over the same Lie algebra and
          there is a nonzero Verma module morphism from ``R``
          into ``self``.

        EXAMPLES::

            sage: L = lie_algebras.so(QQ, 8)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: Mpp = L.verma_module(M.highest_weight().dot_action([1,2]) + La[1])
            sage: M._coerce_map_from_(Mp) is not None
            True
            sage: Mp._coerce_map_from_(M)
            sage: M._coerce_map_from_(Mpp)
            sage: M._coerce_map_from_(ZZ)
            True
        """
        if self.base_ring().has_coerce_map_from(R):
            return True
        if isinstance(R, VermaModule) and R._g is self._g:
            H = Hom(R, self)
            if H.dimension() == 1:
                return H.natural_map()
        return super()._coerce_map_from_(R)

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + 2*La[2])
            sage: M(3)
            3*v[Lambda[1] + 2*Lambda[2]]
            sage: pbw = M.pbw_basis()
            sage: [M(g) for g in pbw.gens()]
            [0,
             0,
             0,
             v[Lambda[1] + 2*Lambda[2]],
             2*v[Lambda[1] + 2*Lambda[2]],
             f[-alpha[2]]*v[Lambda[1] + 2*Lambda[2]],
             f[-alpha[1]]*v[Lambda[1] + 2*Lambda[2]],
             f[-alpha[1] - alpha[2]]*v[Lambda[1] + 2*Lambda[2]]]
        """
        if x in self.base_ring():
            return self._from_dict({self._indices.one(): x})
        if isinstance(x, self._pbw.element_class):
            return self.highest_weight_vector()._acted_upon_(x, False)
        return super()._element_constructor_(self, x)

    def contravariant_form(self, x, y):
        r"""
        Return the contravariant form of ``x`` and ``y``.

        Let `C(x, y)` denote the (universal) contravariant form on
        `U(\mathfrak{g})`. Then the contravariant form on `M(\lambda)` is
        given by evaluating `C(x, y) \in U(\mathfrak{h})` at `\lambda`.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(2*La[1])
            sage: U = M.pbw_basis()
            sage: v = M.highest_weight_vector()
            sage: e, h, f = U.algebra_generators()
            sage: elts = [f^k * v for k in range(8)]; elts
            [v[2*Lambda[1]], f[-alpha[1]]*v[2*Lambda[1]],
             f[-alpha[1]]^2*v[2*Lambda[1]], f[-alpha[1]]^3*v[2*Lambda[1]],
             f[-alpha[1]]^4*v[2*Lambda[1]], f[-alpha[1]]^5*v[2*Lambda[1]],
             f[-alpha[1]]^6*v[2*Lambda[1]], f[-alpha[1]]^7*v[2*Lambda[1]]]
            sage: matrix([[M.contravariant_form(x, y) for x in elts] for y in elts])
            [1 0 0 0 0 0 0 0]
            [0 2 0 0 0 0 0 0]
            [0 0 4 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
        """
        pbw = self._pbw
        I = pbw._indices
        xlift = pbw.element_class(pbw, {I(m._monomial): c for m, c in x._monomial_coefficients.items()})
        ylift = pbw.element_class(pbw, {I(m._monomial): c for m, c in y._monomial_coefficients.items()})
        univ = pbw.contravariant_form(xlift, ylift)
        la = self._weight
        R = self.base_ring()
        return R.sum(c * R.prod(la.scalar(k) ** e for k, e in m._monomial.items())
                     for m, c in univ._monomial_coefficients.items())

    @lazy_attribute
    def _dominant_data(self):
        r"""
        Return the closest to dominant weight in the dot orbit of
        the highest weight of ``self`` and the corresponding reduced word.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: M._dominant_data
            (Lambda[1] + Lambda[2], [])
            sage: M = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: M._dominant_data
            (Lambda[1] + Lambda[2], [1, 2])
            sage: M = L.verma_module(-4*La[1] - La[2])
            sage: M._dominant_data
            (-Lambda[1] + 2*Lambda[2], [1, 2])
        """
        P = self._weight.parent()
        wt, w = (self._weight + P.rho()).to_dominant_chamber(reduced_word=True)
        return (wt - P.rho(), w)

    def is_singular(self):
        r"""
        Return if ``self`` is a singular Verma module.

        A Verma module `M_{\lambda}` is *singular* if there does not
        exist a dominant weight `\tilde{\lambda}` that is in the dot
        orbit of `\lambda`. We call a Verma module *regular* otherwise.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: M.is_singular()
            False
            sage: M = L.verma_module(La[1] - La[2])
            sage: M.is_singular()
            True
            sage: M = L.verma_module(2*La[1] - 10*La[2])
            sage: M.is_singular()
            False
            sage: M = L.verma_module(-2*La[1] - 2*La[2])
            sage: M.is_singular()
            False
            sage: M = L.verma_module(-4*La[1] - La[2])
            sage: M.is_singular()
            True
        """
        return not self._dominant_data[0].is_dominant()

    def is_simple(self):
        r"""
        Return if ``self`` is a simple module.

        A Verma module `M_{\lambda}` is simple if and only if `\lambda`
        is *Verma antidominant* in the sense

        .. MATH::

            \langle \lambda + \rho, \alpha^{\vee} \rangle \notin \ZZ_{>0}

        for all positive roots `\alpha`.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: L.verma_module(La[1] + La[2]).is_simple()
            False
            sage: L.verma_module(-La[1] - La[2]).is_simple()
            True
            sage: L.verma_module(3/2*La[1] + 1/2*La[2]).is_simple()
            False
            sage: L.verma_module(3/2*La[1] + 1/3*La[2]).is_simple()
            True
            sage: L.verma_module(-3*La[1] + 2/3*La[2]).is_simple()
            True
        """
        return self._weight.is_verma_dominant(positive=False)

    def is_projective(self):
        r"""
        Return if ``self`` is a projective module in Category `\mathcal{O}`.

        A Verma module `M_{\lambda}` is projective (in Category `\mathcal{O}`
        if and only if `\lambda` is *Verma dominant* in the sense

        .. MATH::

            \langle \lambda + \rho, \alpha^{\vee} \rangle \notin \ZZ_{<0}

        for all positive roots `\alpha`.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: L.verma_module(La[1] + La[2]).is_projective()
            True
            sage: L.verma_module(-La[1] - La[2]).is_projective()
            True
            sage: L.verma_module(3/2*La[1] + 1/2*La[2]).is_projective()
            True
            sage: L.verma_module(3/2*La[1] + 1/3*La[2]).is_projective()
            True
            sage: L.verma_module(-3*La[1] + 2/3*La[2]).is_projective()
            False
        """
        return self._weight.is_verma_dominant()

    def homogeneous_component_basis(self, d):
        r"""
        Return a basis for the ``d``-th homogeneous component of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: P = L.cartan_type().root_system().weight_lattice()
            sage: La = P.fundamental_weights()
            sage: al = P.simple_roots()
            sage: mu = 2*La[1] + 3*La[2]
            sage: M = L.verma_module(mu)
            sage: M.homogeneous_component_basis(mu - al[2])
            [f[-alpha[2]]*v[2*Lambda[1] + 3*Lambda[2]]]
            sage: M.homogeneous_component_basis(mu - 3*al[2])
            [f[-alpha[2]]^3*v[2*Lambda[1] + 3*Lambda[2]]]
            sage: M.homogeneous_component_basis(mu - 3*al[2] - 2*al[1])
            [f[-alpha[2]]*f[-alpha[1] - alpha[2]]^2*v[2*Lambda[1] + 3*Lambda[2]],
             f[-alpha[2]]^2*f[-alpha[1]]*f[-alpha[1] - alpha[2]]*v[2*Lambda[1] + 3*Lambda[2]],
             f[-alpha[2]]^3*f[-alpha[1]]^2*v[2*Lambda[1] + 3*Lambda[2]]]
            sage: M.homogeneous_component_basis(mu - La[1])
            Family ()
        """
        diff = (d - self._weight)._to_root_vector()
        if diff is None or not all(coeff <= 0 and coeff in ZZ for coeff in diff):
            return Family([])
        return sorted(self._homogeneous_component_f(diff))

    weight_space_basis = homogeneous_component_basis

    @cached_method
    def _homogeneous_component_f(self, d):
        r"""
        Return a basis of the PBW given by ``d`` expressed in the
        root lattice in terms of the simple roots.

        INPUT:

        - ``d`` -- the coefficients of the simple roots as a vector

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: sorted(M._homogeneous_component_f(vector([-1,-2])), key=str)
            [f[-alpha[2]]*f[-alpha[1] - alpha[2]]*v[Lambda[1] + Lambda[2]],
             f[-alpha[2]]^2*f[-alpha[1]]*v[Lambda[1] + Lambda[2]]]
            sage: sorted(M._homogeneous_component_f(vector([-5,-4])), key=str)
            [f[-alpha[1]]*f[-alpha[1] - alpha[2]]^4*v[Lambda[1] + Lambda[2]],
             f[-alpha[2]]*f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]^3*v[Lambda[1] + Lambda[2]],
             f[-alpha[2]]^2*f[-alpha[1]]^3*f[-alpha[1] - alpha[2]]^2*v[Lambda[1] + Lambda[2]],
             f[-alpha[2]]^3*f[-alpha[1]]^4*f[-alpha[1] - alpha[2]]*v[Lambda[1] + Lambda[2]],
             f[-alpha[2]]^4*f[-alpha[1]]^5*v[Lambda[1] + Lambda[2]]]
        """
        if not d:
            return frozenset([self.highest_weight_vector()])
        f = {i: self._pbw(g) for i,g in enumerate(self._g.f())}
        basis = d.parent().basis() # Standard basis vectors
        ret = set()

        def degree(m):
            m = m.dict()
            if not m:
                return d.parent().zero()
            return sum(e * self._g.degree_on_basis(k) for k,e in m.items()).to_vector()
        for i in f:
            if d[i] == 0:
                continue
            for b in self._homogeneous_component_f(d + basis[i]):
                temp = f[i] * b
                ret.update([self.monomial(m) for m in temp.support() if degree(m) == d])
        return frozenset(ret)

    def _Hom_(self, Y, category=None, **options):
        r"""
        Return the homset from ``self`` to ``Y`` in the
        category ``category``.

        INPUT:

        - ``Y`` -- an object
        - ``category`` -- a subcategory of :class:`Crystals`() or ``None``

        The sole purpose of this method is to construct the homset as a
        :class:`~sage.algebras.lie_algebras.verma_module.VermaModuleHomset`.
        If ``category`` is specified and is not a subcategory of
        ``self.category()``, a :exc:`TypeError` is raised instead.

        This method is not meant to be called directly. Please use
        :func:`sage.categories.homset.Hom` instead.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(3*La[1] - 3*La[2])
            sage: H = Hom(M, Mp)
            sage: type(H)
            <...VermaModuleHomset_with_category_with_equality_by_id'>
        """
        from sage.algebras.lie_algebras.bgg_dual_module import BGGDualModule, SimpleModule
        if not ((isinstance(Y, (VermaModule, SimpleModule))
                 or (isinstance(Y, BGGDualModule) and Y._module is self))
                and self._g is Y.lie_algebra()):
            raise TypeError("{} must be an object in Category O of {}".format(Y, self._g))
        if category is not None and not category.is_subcategory(self.category()):
            raise TypeError("{} is not a subcategory of {}".format(category, self.category()))
        return VermaModuleHomset(self, Y)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES:

            Check that other PBW algebras have an action::

                sage: L = lie_algebras.sp(QQ, 6)
                sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
                sage: M = L.verma_module(La[1] - 3*La[2])
                sage: PBW = L.pbw_basis()
                sage: F1 = PBW(L.f(1))
                sage: F1 * M.highest_weight_vector()
                f[-alpha[1]]*v[Lambda[1] - 3*Lambda[2]]
                sage: F1.parent() is M.pbw_basis()
                False
                sage: F1 * M.highest_weight_vector()
                f[-alpha[1]]*v[Lambda[1] - 3*Lambda[2]]
                sage: E1 = PBW(L.e(1))
                sage: E1 * F1
                PBW[alpha[1]]*PBW[-alpha[1]]
                sage: E1 * F1 * M.highest_weight_vector()
                v[Lambda[1] - 3*Lambda[2]]
                sage: M.pbw_basis()(E1 * F1)
                PBW[-alpha[1]]*PBW[alpha[1]] + PBW[alphacheck[1]]
            """
            P = self.parent()
            # Check for scalars first
            # TODO: Pass by these checks if a PBW basis element of the Lie algebra
            if scalar in P.base_ring():
                # Don't have this be a super call
                return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)

            # Check for Lie algebra elements
            try:
                scalar = P._g(scalar)
            except (ValueError, TypeError):
                pass

            # Check for PBW elements
            try:
                scalar = P._pbw(scalar)
            except (ValueError, TypeError):
                # Cannot be made into a PBW element, so propagate it up
                return CombinatorialFreeModule.Element._acted_upon_(self,
                        scalar, self_on_left)

            # We only implement x * self, i.e., as a left module
            if self_on_left:
                return None

            # Lift ``self`` to the PBW basis and do multiplication there
            mc = self._monomial_coefficients
            d = {P._pbw._indices(x.dict()): mc[x] for x in mc} # Lift the index set
            ret = scalar * P._pbw._from_dict(d, remove_zeros=False, coerce=False)

            # Now have ``ret`` act on the highest weight vector
            d = {}
            for m in ret._monomial_coefficients:
                c = ret._monomial_coefficients[m]
                mp = {}
                for k,e in reversed(m._sorted_items()):
                    part = P._g._part_on_basis(k)
                    if part > 0:
                        mp = None
                        break
                    elif part == 0:
                        c *= P._g._weight_action(k, P._weight)**e
                    else:
                        mp[k] = e
                # This term is 0, so nothing to do
                if mp is None:
                    continue
                # Convert back to an element of the indexing set
                mp = P._indices(mp)
                if mp in d:
                    d[mp] += c
                else:
                    d[mp] = c
            return P._from_dict(d)

        _lmul_ = _acted_upon_
        _rmul_ = _acted_upon_


#####################################################################
# Morphisms and Homset


class VermaModuleMorphism(Morphism):
    r"""
    A morphism of a Verma module to another module in Category `\mathcal{O}`.
    """
    def __init__(self, parent, scalar):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: phi = Hom(Mp, M).natural_map()
            sage: TestSuite(phi).run()
        """
        self._scalar = scalar
        Morphism.__init__(self, parent)

    def _repr_type(self):
        """
        Return a string describing the specific type of this map,
        to be used when printing ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: phi = Hom(Mp, M).natural_map()
            sage: phi._repr_type()
            'Verma module'
        """
        return "Verma module"

    def _repr_defn(self):
        r"""
        Return a string describing the definition of ``self``,
        to be used when printing ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: phi = Hom(Mp, M).natural_map()
            sage: phi._repr_defn()
            'v[-5*Lambda[1] + Lambda[2]] |--> f[-alpha[2]]^2*f[-alpha[1]]^4*v[Lambda[1]
              + Lambda[2]] + 8*f[-alpha[2]]*f[-alpha[1]]^3*f[-alpha[1] - alpha[2]]*v[Lambda[1]
              + Lambda[2]] + 12*f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]^2*v[Lambda[1] + Lambda[2]]'

            alpha[1]]^2*f[-alpha[1] - alpha[2]]^2*v[Lambda[1] + Lambda[2]]'
            sage: psi = Hom(M, Mp).natural_map()
            sage: psi
            Verma module morphism:
              From: Verma module with highest weight Lambda[1] + Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              To:   Verma module with highest weight -5*Lambda[1] + Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              Defn: v[Lambda[1] + Lambda[2]] |--> 0
        """
        v = self.domain().highest_weight_vector()
        if not self._scalar:
            return "{} |--> {}".format(v, self.codomain().zero())
        return "{} |--> {}".format(v, self._scalar * self.parent().highest_weight_image())

    def _richcmp_(self, other, op):
        r"""
        Return whether this morphism and ``other`` satisfy ``op``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: H = Hom(Mp, M)
            sage: H(1) < H(2)
            True
            sage: H(2) < H(1)
            False
            sage: H.zero() == H(0)
            True
            sage: H(3) <= H(3)
            True
        """
        return richcmp(self._scalar, other._scalar, op)

    def _call_(self, x):
        r"""
        Apply this morphism to ``x``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: pbw = M.pbw_basis()
            sage: f1, f2 = pbw.f()
            sage: v = Mp.highest_weight_vector()
            sage: phi = Hom(Mp, M).natural_map()
            sage: phi(f1 * v) == f1 * phi(v)
            True
            sage: phi(f2 * f1 * v) == f2 * f1 * phi(v)
            True
            sage: phi(f1 * f2 * f1 * v) == f1 * f2 * f1 * phi(v)
            True

            sage: Mpp = L.verma_module(M.highest_weight().dot_action([1,2]) + La[1])
            sage: psi = Hom(Mpp, M).natural_map()
            sage: v = Mpp.highest_weight_vector()
            sage: psi(v)
            0
        """
        if not self._scalar or not self.parent().highest_weight_image():
            return self.codomain().zero()
        mc = x.monomial_coefficients(copy=False)
        return self.codomain().linear_combination((self._on_basis(m), self._scalar * c)
                                                  for m,c in mc.items())

    def _on_basis(self, m):
        r"""
        Return the image of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: pbw = M.pbw_basis()
            sage: f1, f2 = pbw.f()
            sage: v = Mp.highest_weight_vector()
            sage: phi = Hom(Mp, M).natural_map()
            sage: phi._on_basis((f1 * v).leading_support()) == f1 * phi(v)
            True
        """
        vec = self.parent().highest_weight_image()
        if not vec:
            return vec
        pbw = self.codomain()._pbw
        return pbw.monomial(pbw._indices(m.dict())) * vec

    def _add_(self, other):
        """
        Add ``self`` and ``other``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: phi = Hom(Mp, M).natural_map()
            sage: (phi + 3/2 * phi)._scalar
            5/2
        """
        return type(self)(self.parent(), self._scalar + other._scalar)

    def _sub_(self, other):
        """
        Subtract ``self`` and ``other``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: phi = Hom(Mp, M).natural_map()
            sage: (phi - 3/2 * phi)._scalar
            -1/2
        """
        return type(self)(self.parent(), self._scalar - other._scalar)

    def _acted_upon_(self, other, self_on_left):
        """
        Return the action of ``other`` on ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: phi = Hom(Mp, M).natural_map()
            sage: phi._scalar
            1
            sage: (0 * phi)._scalar
            0
            sage: R.<x> = QQ[]
            sage: x * phi
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: ...
        """
        R = self.parent().base_ring()
        if other not in R:
            return None
        return type(self)(self.parent(), R(other) * self._scalar)

    def _composition_(self, right, homset):
        r"""
        Return the composition of ``self`` and ``right``.

        INPUT:

        - ``self``, ``right`` -- maps
        - ``homset`` -- a homset

        ASSUMPTION:

        The codomain of ``right`` is contained in the domain of ``self``.
        This assumption is not verified.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: Mpp = L.verma_module(M.highest_weight().dot_action([1,2]) + La[1])
            sage: phi = Hom(Mp, M).natural_map()
            sage: psi = Hom(Mpp, Mp).natural_map()
            sage: xi = phi * psi
            sage: xi._scalar
            0
        """
        if (isinstance(right, VermaModuleMorphism)
            and right.domain()._g is self.codomain()._g):
            return homset.element_class(homset, right._scalar * self._scalar)
        return super()._composition_(right, homset)

    def is_injective(self):
        r"""
        Return if ``self`` is injective or not.

        A morphism `\phi : M \to M'` from a Verma module `M` to another
        Verma module `M'` is injective if and only if `\dim \hom(M, M') = 1`
        and `\phi \neq 0`. If `M'` is a dual Verma or simple module, then
        the result is not injective.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: Mpp = L.verma_module(M.highest_weight().dot_action([1,2]) + La[1])
            sage: phi = Hom(Mp, M).natural_map()
            sage: phi.is_injective()
            True
            sage: (0 * phi).is_injective()
            False
            sage: psi = Hom(Mpp, Mp).natural_map()
            sage: psi.is_injective()
            False
        """
        if not isinstance(self.codomain(), VermaModule):
            return False
        return bool(self._scalar)

    def is_surjective(self):
        r"""
        Return if ``self`` is surjective or not.

        A morphism `\phi : M \to M'` from a Verma module `M` to another
        Verma module `M'` is surjective if and only if the domain is
        equal to the codomain and it is not the zero morphism.

        If `M'` is a simple module, then this surjective if and only if
        `\dim \hom(M, M') = 1` and `\phi \neq 0`.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: phi = Hom(M, M).natural_map()
            sage: phi.is_surjective()
            True
            sage: (0 * phi).is_surjective()
            False
            sage: psi = Hom(Mp, M).natural_map()
            sage: psi.is_surjective()
            False
        """
        if not bool(self._scalar):
            return False

        if isinstance(self.codomain(), VermaModule):
            return self.domain() == self.codomain()

        from sage.algebras.lie_algebras.bgg_dual_module import SimpleModule
        if isinstance(self.codomain(), SimpleModule):
            return self.domain().highest_weight() == self.codomain().highest_weight()

        return False

    def image(self):
        r"""
        Return the image of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(La[1] + 2*La[2])
            sage: Mp = g.verma_module(La[1] + 3*La[2])
            sage: phi = Hom(M, Mp).natural_map()
            sage: phi.image()
            Free module generated by {} over Rational Field
            sage: Mc = M.dual()
            sage: phi = Hom(M, Mc).natural_map()
            sage: L = phi.image(); L
            Simple module with highest weight Lambda[1] + 2*Lambda[2] of
             Lie algebra of ['B', 2] in the Chevalley basis
            sage: psi = Hom(M, L).natural_map()
            sage: psi.image()
            Simple module with highest weight Lambda[1] + 2*Lambda[2] of
             Lie algebra of ['B', 2] in the Chevalley basis
        """
        C = self.codomain()
        if not bool(self._scalar):
            return C.submodule([])

        if isinstance(C, VermaModule):
            if self.domain() == C:
                return C
            raise NotImplementedError("submodules of Verma modules not yet implemented")

        from sage.algebras.lie_algebras.bgg_dual_module import BGGDualModule, SimpleModule
        if isinstance(C, BGGDualModule) and isinstance(C._module, VermaModule):
            return SimpleModule(C.lie_algebra(), C.highest_weight(), prefix=C._indices.prefix(),
                                basis_key=C._module._basis_key)

        if isinstance(self.codomain(), SimpleModule):
            return self.codomain()


class VermaModuleHomset(Homset):
    r"""
    The set of morphisms from a Verma module to another module in
    Category `\mathcal{O}` considered as `U(\mathfrak{g})`-representations.

    This currently assumes the codomain is a Verma module, its dual,
    or a simple module.

    Let `M_{w \cdot \lambda}` and `M_{w' \cdot \lambda'}` be
    Verma modules, `\cdot` is the dot action, and `\lambda + \rho`,
    `\lambda' + \rho` are dominant weights. Then we have

    .. MATH::

        \dim \hom(M_{w \cdot \lambda}, M_{w' \cdot \lambda'}) = 1

    if and only if `\lambda = \lambda'` and `w' \leq w` in Bruhat
    order. Otherwise the homset is 0 dimensional.

    If the codomain is a dual Verma module `M_{\mu}^{\vee}`, then the
    homset is `\delta_{\lambda\mu}` dimensional. When `\mu = \lambda`,
    the image is the simple module `L_{\lambda}`.
    """
    def __call__(self, x, **options):
        r"""
        Construct a morphism in this homset from ``x`` if possible.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: Mpp = L.verma_module(M.highest_weight().dot_action([1,2,1]))
            sage: phi = Hom(Mp, M).natural_map()
            sage: Hom(Mpp, M)(phi)
            Verma module morphism:
              From: Verma module with highest weight -3*Lambda[1] - 3*Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              To:   Verma module with highest weight Lambda[1] + Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              Defn: v[-3*Lambda[1] - 3*Lambda[2]] |-->
                     f[-alpha[2]]^4*f[-alpha[1]]^4*v[Lambda[1] + Lambda[2]]
                       + 8*f[-alpha[2]]^3*f[-alpha[1]]^3*f[-alpha[1] - alpha[2]]*v[Lambda[1] + Lambda[2]]
                       + 12*f[-alpha[2]]^2*f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]^2*v[Lambda[1] + Lambda[2]]

            sage: psi = Hom(Mpp, Mp).natural_map()
            sage: Hom(Mpp, M)(psi)
            Verma module morphism:
              From: Verma module with highest weight -3*Lambda[1] - 3*Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              To:   Verma module with highest weight Lambda[1] + Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              Defn: v[-3*Lambda[1] - 3*Lambda[2]] |-->
                     f[-alpha[2]]^4*f[-alpha[1]]^4*v[Lambda[1] + Lambda[2]]
                      + 8*f[-alpha[2]]^3*f[-alpha[1]]^3*f[-alpha[1] - alpha[2]]*v[Lambda[1] + Lambda[2]]
                      + 12*f[-alpha[2]]^2*f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]^2*v[Lambda[1] + Lambda[2]]
        """
        if isinstance(x, VermaModuleMorphism):
            if x.parent() is self:
                return x
            if x.parent() == self:
                x._set_parent(self) # needed due to non-uniqueness of homsets
                return x

            if x.domain() != self.domain():
                x = x * Hom(self.domain(), x.domain()).natural_map()
            if x.codomain() != self.codomain():
                x = Hom(x.codomain(), self.codomain()).natural_map() * x

            return x

        if x in self.base_ring():
            if self.singular_vector() is None:
                return self.zero()
            return self.element_class(self, self.base_ring()(x))

        return super().__call__(x, **options)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([2]))
            sage: H = Hom(Mp, M)
            sage: H._an_element_()
            Verma module morphism:
              From: Verma module with highest weight 3*Lambda[1] - 3*Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              To:   Verma module with highest weight Lambda[1] + Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              Defn: v[3*Lambda[1] - 3*Lambda[2]] |-->
                     f[-alpha[2]]^2*v[Lambda[1] + Lambda[2]]
        """
        return self.natural_map()

    def highest_weight_image(self):
        r"""
        Return the image of the highest weight vector of the domain
        in the codomain.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['C', 3])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(La[1] + 2*La[3])
            sage: Mc = M.dual()
            sage: H = Hom(M, Mc)
            sage: H.highest_weight_image()
            v[Lambda[1] + 2*Lambda[3]]^*
            sage: L = H.natural_map().image()
            sage: Hp = Hom(M, L)
            sage: Hp.highest_weight_image()
            u[Lambda[1] + 2*Lambda[3]]
        """
        C = self.codomain()
        if isinstance(C, VermaModule):
            # singular_vector() is cached, so we can safely call it twice
            if self.singular_vector() is None:
                return C.zero()
            return self.singular_vector()
        # Otherwise, it is a dual Verma or a simple, so the image
        #   must be the highest weight vector.
        if self.domain().highest_weight() == C.highest_weight():
            return C.highest_weight_vector()
        return C.zero()

    @cached_method
    def singular_vector(self):
        r"""
        Return the singular vector in the codomain corresponding
        to the domain's highest weight element or ``None`` if no
        such element exists.

        ALGORITHM:

        We essentially follow the algorithm laid out in [deG2005]_.
        We split the main computation into two cases. If there exists
        an `i` such that `\langle \lambda + \rho, \alpha_i^{\vee}
        \rangle = m > 0` (i.e., the weight `\lambda` is `i`-dominant
        with respect to the dot action), then we use the `\mathfrak{sl}_2`
        relation on `M_{s_i \cdot \lambda} \to M_{\lambda}` to
        construct the singular vector `f_i^m v_{\lambda}`. Otherwise
        we find the shortest root `\alpha` such that `\langle \lambda
        + \rho, \alpha^{\vee} \rangle > 0` and explicitly compute the
        kernel with respect to the weight basis elements. We iterate
        this until we reach `\mu`.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: la = La[1] - La[3]
            sage: mu = la.dot_action([1,2])
            sage: M = L.verma_module(la)
            sage: Mp = L.verma_module(mu)
            sage: H = Hom(Mp, M)
            sage: v = H.singular_vector(); v
            f[-alpha[2]]*f[-alpha[1]]^3*v[Lambda[1] - Lambda[3]]
             + 3*f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]*v[Lambda[1] - Lambda[3]]
            sage: v.degree() == Mp.highest_weight()
            True

        ::

            sage: L = LieAlgebra(QQ, cartan_type=['F', 4])
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: la = La[1] + La[2] - La[3]
            sage: mu = la.dot_action([1,2,3,2])
            sage: M = L.verma_module(la)
            sage: Mp = L.verma_module(mu)
            sage: H = Hom(Mp, M)
            sage: v = H.singular_vector()
            sage: pbw = M.pbw_basis()
            sage: E = [pbw(e) for e in L.e()]
            sage: all(e * v == M.zero() for e in E)  # long time
            True
            sage: v.degree() == Mp.highest_weight()
            True

        When `w \cdot \lambda \notin \lambda + Q^-`, there does not
        exist a singular vector::

            sage: L = lie_algebras.sl(QQ, 4)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: la = 3/7*La[1] - 1/2*La[3]
            sage: mu = la.dot_action([1,2])
            sage: M = L.verma_module(la)
            sage: Mp = L.verma_module(mu)
            sage: H = Hom(Mp, M)
            sage: H.singular_vector() is None
            True

        When we need to apply a non-simple reflection, we can compute
        the singular vector (see :issue:`36793`)::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module((0*La[1]).dot_action([1]))
            sage: Mp = g.verma_module((0*La[1]).dot_action([1,2]))
            sage: H = Hom(Mp, M)
            sage: v = H.singular_vector(); v
            1/2*f[-alpha[2]]*f[-alpha[1]]*v[-2*Lambda[1] + Lambda[2]]
             + f[-alpha[1] - alpha[2]]*v[-2*Lambda[1] + Lambda[2]]
            sage: pbw = M.pbw_basis()
            sage: E = [pbw(e) for e in g.e()]
            sage: all(e * v == M.zero() for e in E)
            True
            sage: v.degree() == Mp.highest_weight()
            True

        TESTS::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: al = L.cartan_type().root_system().root_lattice().simple_roots()
            sage: M = L.verma_module(La[1] + La[2])
            sage: pbw = M.pbw_basis()
            sage: E = {i: pbw(L.e(i)) for i in L.cartan_type().index_set()}
            sage: all(not E[i] * Hom(L.verma_module(mu), M).singular_vector()
            ....:     for i in L.cartan_type().index_set()
            ....:     for mu in M.highest_weight().dot_orbit())
            True
        """
        if self.is_endomorphism_set():
            return self.codomain().highest_weight_vector()
        if self.domain()._dominant_data[0] != self.codomain()._dominant_data[0]:
            return None

        from sage.combinat.root_system.coxeter_group import CoxeterGroup
        from sage.matrix.constructor import matrix
        W = CoxeterGroup(self.domain()._g._cartan_type)
        # We take the inverse to account for the left versus right action
        wp = W.from_reduced_word(reversed(self.domain()._dominant_data[1]))
        w = W.from_reduced_word(reversed(self.codomain()._dominant_data[1]))
        if not w.bruhat_le(wp):
            return None
        C = self.codomain()
        pbw = C._pbw
        F = pbw.f()
        E = pbw.e()
        index_set = F.keys()
        cur_w = w
        rho = C._weight.parent().rho()
        ac = C._weight.parent().simple_coroots()
        elt = pbw.one()
        wt = C._weight
        pos_roots_by_ht = C._g._cartan_type.root_system().root_lattice().positive_roots_by_height()
        assert all(sum(rt.coefficients()) == 1 for rt in pos_roots_by_ht[:len(index_set)])
        # for this, we don't need to check the simple roots
        pos_roots_by_ht = pos_roots_by_ht[len(index_set):]

        while cur_w != wp:
            ind = None
            for i in cur_w.descents(side='right', positive=True):
                exp = (wt + rho).scalar(ac[i])
                if exp not in ZZ or exp <= 0:
                    continue
                # We need to check that the result is still smaller in Bruhat order
                next_w = cur_w.apply_simple_reflection_right(i)
                if not next_w.bruhat_le(wp):
                    continue
                ind = i
                # favor a path in weak order so we only do sl_2 relations
                if not next_w.weak_le(wp, side="right"):
                    continue
                break
            if ind is None:  # no simple root; need a more general approach
                # We search for the shortest root that can be applied to minimize
                #   the size of the basis needed to compute the kernel.
                for rt in pos_roots_by_ht:
                    exp = (wt + rho).scalar(rt.associated_coroot())
                    # We need to check that the result is still smaller in Bruhat order
                    i, wd = rt.to_simple_root(reduced_word=True)
                    refl = wd + (i,) + tuple(reversed(wd))
                    next_w = cur_w.apply_reflections(refl, side='right', word_type="simple")
                    if exp not in ZZ or exp <= 0:
                        continue
                    if not next_w.bruhat_le(wp):
                        continue
                    # We construct the Verma module of the appropriate weight in
                    #   order to reduce the dimension and number of multiplications.
                    Mp = C._g.verma_module(wt)
                    basis = sorted(Mp._homogeneous_component_f(-rt.to_vector()), key=str)
                    for i in index_set:
                        image = [E[i] * b for b in basis]
                        supp = set()
                        for vec in image:
                            supp.update(vec._monomial_coefficients)
                        supp = sorted(supp, key=pbw._monomial_key)
                        if not supp:  # everything is in the kernel
                            continue
                        M = matrix(pbw.base_ring(), [[v[s] for v in image] for s in supp])
                        ker = M.right_kernel_matrix()
                        basis = [C.linear_combination((basis[j], c) for j, c in kv.iteritems())
                                 for kv in ker.rows()]

                    assert len(basis) == 1
                    if Mp is C:  # We've constructed the element in the codomain
                        assert next_w == wp
                        assert basis[0].degree() == self.domain().highest_weight()
                        return basis[0]
                    pbw_elt = pbw.element_class(pbw, {pbw._indices(m._monomial): c
                                                      for m, c in basis[0]._monomial_coefficients.items()})
                    elt = pbw_elt * elt
                    wt = wt.dot_action(refl)
                    cur_w = next_w
                    break
                else:
                    #assert False, "unable to find root"
                    # Have a more explicit check at the beginning using the integral
                    #   orbit action for the correct version of dominance; see, e.g.,
                    #   Humphreys "Representations of Semisimple Lie Algebras in the BGG Category O".
                    return None
            else:
                # Construct the singular vector by iterated embeddings of Verma
                #   modules from the sl_2 relations (without constructing
                #   the modules themselves)
                elt = F[ind]**ZZ(exp) * elt
                wt = wt.dot_action([ind])
                cur_w = cur_w.apply_simple_reflection_right(ind)
        ret = C.highest_weight_vector()._acted_upon_(elt, False)
        assert ret.degree() == self.domain().highest_weight()
        return ret

    @cached_method
    def natural_map(self):
        """
        Return the "natural map" of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([2]))
            sage: H = Hom(Mp, M)
            sage: H.natural_map()
            Verma module morphism:
              From: Verma module with highest weight 3*Lambda[1] - 3*Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              To:   Verma module with highest weight Lambda[1] + Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              Defn: v[3*Lambda[1] - 3*Lambda[2]] |-->
                     f[-alpha[2]]^2*v[Lambda[1] + Lambda[2]]

            sage: Mp = L.verma_module(La[1] + 2*La[2])
            sage: H = Hom(Mp, M)
            sage: H.natural_map()
            Verma module morphism:
              From: Verma module with highest weight Lambda[1] + 2*Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              To:   Verma module with highest weight Lambda[1] + Lambda[2]
                     of Lie algebra of ['A', 2] in the Chevalley basis
              Defn: v[Lambda[1] + 2*Lambda[2]] |--> 0
        """
        if not self.highest_weight_image():
            return self.zero()
        return self.element_class(self, self.base_ring().one())

    @cached_method
    def zero(self):
        """
        Return the zero morphism of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = L.verma_module(La[1] + 2/3*La[2])
            sage: Mp = L.verma_module(La[2] - La[3])
            sage: H = Hom(Mp, M)
            sage: H.zero()
            Verma module morphism:
              From: Verma module with highest weight Lambda[2] - Lambda[3]
                     of Lie algebra of ['C', 3] in the Chevalley basis
              To:   Verma module with highest weight Lambda[1] + 2/3*Lambda[2]
                     of Lie algebra of ['C', 3] in the Chevalley basis
              Defn: v[Lambda[2] - Lambda[3]] |--> 0
        """
        return self.element_class(self, self.base_ring().zero())

    def dimension(self):
        r"""
        Return the dimension of ``self`` (as a vector space over
        the base ring).

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([2]))
            sage: H = Hom(Mp, M)
            sage: H.dimension()
            1

            sage: Mp = L.verma_module(La[1] + 2*La[2])
            sage: H = Hom(Mp, M)
            sage: H.dimension()
            0
        """
        if not self.highest_weight_image():
            return ZZ.zero()
        return ZZ.one()

    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([2]))
            sage: H = Hom(Mp, M)
            sage: list(H.basis()) == [H.natural_map()]
            True

            sage: Mp = L.verma_module(La[1] + 2*La[2])
            sage: H = Hom(Mp, M)
            sage: H.basis()
            Family ()
        """
        if not self.highest_weight_image():
            return Family([])
        return Family([self.natural_map()])

    Element = VermaModuleMorphism
