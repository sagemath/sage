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


class VermaModule(CombinatorialFreeModule):
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
    `\langle \lambda, h \rangle` and `U\mathfrak{g}^+) F_{\lambda} = 0`.

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

    def _triangular_key(self, x):
        """
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

    def _repr_generator(self, m):
        r"""
        Return a string representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 4)
            sage: La = L.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: M = L.verma_module(-1/2*La[1] + 3/7*La[2])
            sage: f1, f2 = L.f(1), L.f(2)
            sage: x = M.pbw_basis()(L([f1, [f1, f2]]))
            sage: v = x * M.highest_weight_vector()
            sage: M._repr_generator(v.leading_support())
            'f[-2*alpha[1] - alpha[2]]*v[(-1/14, 3/7)]'

            sage: M.highest_weight_vector()
            v[(-1/14, 3/7)]
            sage: 2 * M.highest_weight_vector()
            2*v[(-1/14, 3/7)]
        """
        ret = super(VermaModule, self)._repr_generator(m)
        if ret == '1':
            ret = ''
        else:
            ret += '*'
        return ret + "v[{}]".format(self._weight)

    def _latex_generator(self, m):
        r"""
        Return a latex representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 4)
            sage: La = L.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: M = L.verma_module(-1/2*La[1] + 3/7*La[2])
            sage: f1, f2 = L.f(1), L.f(2)
            sage: x = M.pbw_basis()(L([f1, [f1, f2]]))
            sage: v = x * M.highest_weight_vector()
            sage: M._latex_generator(v.leading_support())
            f_{-2 \alpha_{1} - \alpha_{2}} v_{-\frac{1}{14} e_{0} + \frac{3}{7} e_{1}}

            sage: latex(2 * M.highest_weight_vector())
            2  v_{-\frac{1}{14} e_{0} + \frac{3}{7} e_{1}}
            sage: latex(M.highest_weight_vector())
             v_{-\frac{1}{14} e_{0} + \frac{3}{7} e_{1}}
        """
        ret = super(VermaModule, self)._latex_generator(m)
        if ret == '1':
            ret = ''
        from sage.misc.latex import latex
        return ret + " v_{{{}}}".format(latex(self._weight))

    _repr_term = _repr_generator
    _latex_term = _latex_generator

    def lie_algebra(self):
        """
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
        """
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
          there is a non-zero Verma module morphism from ``R``
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
        return super(VermaModule, self)._coerce_map_from_(R)

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
        return super(VermaModule, self)._element_constructor_(self, x)

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
        diff = _convert_wt_to_root(d - self._weight)
        if diff is None or not all(coeff <= 0 and coeff in ZZ for coeff in diff):
            return Family([])
        return sorted(self._homogeneous_component_f(diff))

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
        ``self.category()``, a ``TypeError`` is raised instead.

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
        if not (isinstance(Y, VermaModule) and self._g is Y._g):
            raise TypeError("{} must be a Verma module of {}".format(Y, self._g))
        if category is not None and not category.is_subcategory(self.category()):
            raise TypeError("{} is not a subcategory of {}".format(category, self.category()))
        return VermaModuleHomset(self, Y)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

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
## Morphisms and Homset

class VermaModuleMorphism(Morphism):
    """
    A morphism of Verma modules.
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
        return "{} |--> {}".format(v, self._scalar * self.parent().singular_vector())

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
            sage: f1, f2 = pbw(L.f(1)), pbw(L.f(2))
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
        if not self._scalar or self.parent().singular_vector() is None:
            return self.codomain().zero()
        mc = x.monomial_coefficients(copy=False)
        return self.codomain().linear_combination((self._on_basis(m), self._scalar * c)
                                                  for m,c in mc.items())

    def _on_basis(self, m):
        """
        Return the image of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 3)
            sage: La = L.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = L.verma_module(La[1] + La[2])
            sage: Mp = L.verma_module(M.highest_weight().dot_action([1,2]))
            sage: pbw = M.pbw_basis()
            sage: f1, f2 = pbw(L.f(1)), pbw(L.f(2))
            sage: v = Mp.highest_weight_vector()
            sage: phi = Hom(Mp, M).natural_map()
            sage: phi._on_basis((f1 * v).leading_support()) == f1 * phi(v)
            True
        """
        pbw = self.codomain()._pbw
        return pbw.monomial(pbw._indices(m.dict())) * self.parent().singular_vector()

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
        - homset -- a homset

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
        return super(VermaModuleMorphism, self)._composition_(right, homset)

    def is_injective(self):
        r"""
        Return if ``self`` is injective or not.

        A Verma module morphism `\phi : M \to M'` is injective if
        and only if `\dim \hom(M, M') = 1` and `\phi \neq 0`.

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
        return self.parent().singular_vector() is not None and bool(self._scalar)

    def is_surjective(self):
        """
        Return if ``self`` is surjective or not.

        A Verma module morphism is surjective if and only if the
        domain is equal to the codomain and it is not the zero
        morphism.

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
        return self.domain() == self.codomain() and bool(self._scalar)

class VermaModuleHomset(Homset):
    r"""
    The set of morphisms from one Verma module to another
    considered as `U(\mathfrak{g})`-representations.

    Let `M_{w \cdot \lambda}` and `M_{w' \cdot \lambda'}` be
    Verma modules, `\cdot` is the dot action, and `\lambda + \rho`,
    `\lambda' + \rho` are dominant weights. Then we have

    .. MATH::

        \dim \hom(M_{w \cdot \lambda}, M_{w' \cdot \lambda'}) = 1

    if and only if `\lambda = \lambda'` and `w' \leq w` in Bruhat
    order. Otherwise the homset is 0 dimensional.
    """
    def __call__(self, x, **options):
        """
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

        return super(VermaModuleHomset, self).__call__(x, **options)

    def _an_element_(self):
        """
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

    @cached_method
    def singular_vector(self):
        r"""
        Return the singular vector in the codomain corresponding
        to the domain's highest weight element or ``None`` if no
        such element exists.

        ALGORITHM:

        We essentially follow the algorithm laid out in [deG2005]_.
        We use the `\mathfrak{sl}_2` relation on
        `M_{s_i \cdot \lambda} \to M_{\lambda}`, where
        `\langle \lambda + \delta, \alpha_i^{\vee} \rangle = m > 0`,
        i.e., the weight `\lambda` is `i`-dominant with respect to
        the dot action. From here, we construct the singular vector
        `f_i^m v_{\lambda}`. We iterate this until we reach `\mu`.

        EXAMPLES::

            sage: L = lie_algebras.sp(QQ, 6)
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: la = La[1] - La[3]
            sage: mu = la.dot_action([1,2])
            sage: M = L.verma_module(la)
            sage: Mp = L.verma_module(mu)
            sage: H = Hom(Mp, M)
            sage: H.singular_vector()
            f[-alpha[2]]*f[-alpha[1]]^3*v[Lambda[1] - Lambda[3]]
             + 3*f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]*v[Lambda[1] - Lambda[3]]

        ::

            sage: L = LieAlgebra(QQ, cartan_type=['F',4])
            sage: La = L.cartan_type().root_system().weight_space().fundamental_weights()
            sage: la = La[1] + La[2] - La[3]
            sage: mu = la.dot_action([1,2,3,2])
            sage: M = L.verma_module(la)
            sage: Mp = L.verma_module(mu)
            sage: H = Hom(Mp, M)
            sage: v = H.singular_vector()
            sage: pbw = M.pbw_basis()
            sage: E = [pbw(e) for e in L.e()]
            sage: all(e * v == M.zero() for e in E)
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
        W = CoxeterGroup(self.domain()._g._cartan_type)
        wp = W.from_reduced_word(self.domain()._dominant_data[1])
        w = W.from_reduced_word(self.codomain()._dominant_data[1])
        if not w.bruhat_le(wp):
            return None
        C = self.codomain()
        pbw = C._pbw
        f = C._g.f()
        F = {i: pbw(f[i]) for i in f.keys()}
        red_word = (wp * ~w).reduced_word()
        rho = C._weight.parent().rho()
        ac = C._weight.parent().simple_coroots()
        elt = pbw.one()
        wt = C._weight
        # Construct the singular vector by iterated embeddings of Verma
        #   modules (without constructing the modules themselves)
        for i in reversed(red_word):
            exp = (wt + rho).scalar(ac[i])
            if exp not in ZZ or exp < 0:
                return None
            elt = F[i]**ZZ(exp) * elt
            wt = wt.dot_action([i])
        return C.highest_weight_vector()._acted_upon_(elt, False)

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
        if self.singular_vector() is None:
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
        """
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
        if self.singular_vector() is None:
            return ZZ.zero()
        return ZZ.one()

    def basis(self):
        """
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
        if self.singular_vector() is None:
            return Family([])
        return Family([self.natural_map()])

    Element = VermaModuleMorphism


def _convert_wt_to_root(wt):
    r"""
    Helper function to express ``wt`` as a linear combination
    of simple roots.

    INPUT:

    - ``wt`` -- an element of a weight lattice realization

    OUTPUT:

    A vector over `\QQ` representing ``wt`` as a linear combination
    of simple roots.

    EXAMPLES::

        sage: from sage.algebras.lie_algebras.verma_module import _convert_wt_to_root
        sage: P = RootSystem(['A',3]).weight_lattice()
        sage: La = P.fundamental_weights()
        sage: [_convert_wt_to_root(al) for al in P.simple_roots()]
        [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        sage: _convert_wt_to_root(La[1] + La[2])
        (5/4, 3/2, 3/4)

        sage: L = RootSystem(['A',3]).ambient_space()
        sage: e = L.basis()
        sage: _convert_wt_to_root(e[0] + 3*e[3])
        sage: _convert_wt_to_root(e[0] - e[1])
        (1, 0, 0)
        sage: _convert_wt_to_root(e[0] + 2*e[1] - 3*e[2])
        (1, 3, 0)
    """
    v = wt.to_vector().change_ring(QQ)
    al = [a.to_vector() for a in wt.parent().simple_roots()]
    b = v.parent().linear_dependence([v] + al)
    if len(b) != 1 or b[0] == 0:
        return None
    b = b[0]  # Get the actual vector that gives the linear dependency
    # Get v as a linear combination of the simple roots
    return vector(QQ, [-x / b[0] for x in b[1:]])
