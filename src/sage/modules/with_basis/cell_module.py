# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.combinat
r"""
Cell modules
"""
# ****************************************************************************
#       Copyright (C) 2015-2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.structure.element import Element
from sage.combinat.free_module import CombinatorialFreeModule
from sage.data_structures.blas_dict import linear_combination
from sage.modules.with_basis.subquotient import QuotientModuleWithBasis


class CellModule(CombinatorialFreeModule):
    r"""
    A cell module.

    Let `R` be a commutative ring. Let `A` be a cellular `R`-algebra
    with cell datum `(\Lambda, i, M, C)`. A *cell module* `W(\lambda)`
    is the `R`-module given by `R\{C_s \mid s \in M(\lambda)\}` with
    an action of `a \in A` given by `a C_s = \sum_{u \in M(\lambda)}
    r_a(u, s) C_u`, where `r_a(u, s)` is the same as those
    given by the cellular condition:

    .. MATH::

        a C^{\lambda}_{st} = \sum_{u \in M(\lambda)} r_a(u, s)
        C_{ut}^{\lambda} +
        \sum_{\substack{\mu < \lambda \\ x,y \in M(\mu)}} R C^{\mu}_{xy}.

    INPUT:

    - ``A`` -- a cellular algebra
    - ``mu`` -- an element of the cellular poset of ``A``

    .. SEEALSO::

        :class:`~sage.algebras.cellular_basis.CellularBasis`

    AUTHORS:

    - Travis Scrimshaw (2015-11-5): Initial version

    REFERENCES:

    - [GrLe1996]_
    - [KX1998]_
    - [Mat1999]_
    - :wikipedia:`Cellular_algebra`
    - http://webusers.imj-prg.fr/~bernhard.keller/ictp2006/lecturenotes/xi.pdf
    """
    @staticmethod
    def __classcall_private__(cls, A, mu, **kwds):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W1 = S.cell_module([2,1])
            sage: W2 = S.cell_module([2,1], prefix='W')
            sage: W1 is W2
            True
        """
        mu = A.cell_poset()(mu)
        kwds['prefix'] = kwds.get('prefix', 'W')
        return super().__classcall__(cls, A, mu, **kwds)

    def __init__(self, A, mu, **kwds):
        r"""
        Initialize ``self``.

        TESTS::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: TestSuite(W).run()
        """
        self._algebra = A
        self._la = mu
        cat = ModulesWithBasis(A.base_ring()).FiniteDimensional()
        CombinatorialFreeModule.__init__(self, A.base_ring(),
                                         A.cell_module_indices(mu),
                                         category=cat, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: S.cell_module([2,1])
            Cell module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field
        """
        return "Cell module indexed by {} of {}".format(self._la, self._algebra)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: latex(W)
            W_{...}\left(...\right)
        """
        from sage.misc.latex import latex
        return "W_{{{}}}\\left({}\\right)".format(latex(self._algebra), latex(self._la))

    def cellular_algebra(self):
        r"""
        Return the cellular algebra of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: W.cellular_algebra() is S.cellular_basis()
            True
            sage: S.has_coerce_map_from(W.cellular_algebra())
            True
        """
        return self._algebra

    @cached_method
    def _action_basis(self, x, s):
        """
        Return the action of the basis element of the cellular algebra
        indexed by ``x`` on the basis element of ``self`` indexed by ``s``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: W = S.cell_module([2,1])
            sage: s = StandardTableaux([2,1])([[1,2],[3]])
            sage: mu = Partition([2,1])
            sage: W._action_basis((mu, s, s), s)
            W[[1, 2], [3]]
        """
        t = self.basis().keys()[0]
        B = self._algebra.basis()
        elt = B[x] * B[(self._la, s, t)]
        mc = elt.monomial_coefficients(copy=False)
        d = {k[1]: mc[k] for k in mc if k[0] == self._la}
        return self._from_dict(d, remove_zeros=False)

    @cached_method
    def _bilinear_form_on_basis(self, s, t):
        """
        Return the result of the bilinear form on the basis elements
        indexed by ``s`` and ``t``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: K = W.basis().keys()
            sage: matrix([[W._bilinear_form_on_basis(s, t) for t in K] for s in K])
            [1 0]
            [0 1]
        """
        B = self._algebra.basis()
        elt = B[(self._la, s, s)] * B[(self._la, t, t)]
        return elt[(self._la, s, t)]

    def bilinear_form(self, x, y):
        r"""
        Return the bilinear form on ``x`` and ``y``.

        The cell module `W(\lambda)` has a canonical bilinear form
        `\Phi_{\lambda} : W(\lambda) \times W(\lambda) \to W(\lambda)`
        given by

        .. MATH::

            C_{ss}^{\lambda} C_{tt}^{\lambda} = \Phi_{\lambda}(C_s, C_t)
            C_{st}^{\lambda} +
            \sum_{\substack{\mu < \lambda \\ x,y \in M(\mu)}} R C^{\mu}_{xy}.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: elt = W.an_element(); elt
            2*W[[1, 2], [3]] + 2*W[[1, 3], [2]]
            sage: W.bilinear_form(elt, elt)
            8
        """
        R = self.base_ring()
        return R.sum(self._bilinear_form_on_basis(s, t) * cx * cy
                     for s, cx in x for t, cy in y)

    def bilinear_form_matrix(self, ordering=None):
        """
        Return the matrix corresponding to the bilinear form
        of ``self``.

        INPUT:

        - ``ordering`` -- (optional) an ordering of the indices

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: W.bilinear_form_matrix()
            [1 0]
            [0 1]
        """
        if ordering is None:
            ordering = sorted(self.basis().keys())
        else:
            sordering = set(ordering)
            if sordering != set(self.basis().keys()) or len(sordering) != len(ordering):
                raise ValueError("not an ordering of the basis indices")
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(self.base_ring(), len(ordering))
        return MS([[self._bilinear_form_on_basis(s, t) for t in ordering]
                   for s in ordering])

    @cached_method
    def nonzero_bilinear_form(self):
        """
        Return ``True`` if the bilinear form of ``self`` is nonzero.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: W.nonzero_bilinear_form()
            True
        """
        C = list(self.basis().keys())
        # Since the bilinear form is symmetric, it is sufficient
        #   to check on the upper triangular part
        return any(self._bilinear_form_on_basis(s, t)
                   for i, s in enumerate(C) for t in C[i:])

    @cached_method
    def radical_basis(self):
        """
        Return a basis of the radical of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: W.radical_basis()
            ()
        """
        mat = self.bilinear_form_matrix(self.get_order())
        return tuple([self.from_vector(v) for v in mat.left_kernel().basis()])

    def radical(self):
        r"""
        Return the radical of ``self``.

        Let `W(\lambda)` denote a cell module. The *radical* of `W(\lambda)`
        is defined as

        .. MATH::

            \operatorname{rad}(\lambda) := \{x \in W(\lambda) \mid
            \Phi_{\lambda}(x, y)\},

        and note that it is a submodule of `W(\lambda)`.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: R = W.radical(); R
            Radical of Cell module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field
            sage: R.basis()
            Finite family {}
        """
        radical = self.submodule(self.radical_basis(),
                                 category=self.category().Subobjects(),
                                 already_echelonized=True)
        radical.rename("Radical of {}".format(self))
        return radical

    def simple_module(self):
        r"""
        Return the corresponding simple module of ``self``.

        Let `W(\lambda)` denote a cell module. The simple module `L(\lambda)`
        is defined as `W(\lambda) / \operatorname{rad}(\lambda)`,
        where `\operatorname{rad}(\lambda)` is the radical of the
        bilinear form `\Phi_{\lambda}`.

        .. SEEALSO::

            :meth:`radical`

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: L = W.simple_module(); L
            Simple module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field
            sage: L.has_coerce_map_from(W)
            True
        """
        if not self.nonzero_bilinear_form():
            raise ValueError("bilinear form is 0; no corresponding simple module")
        return SimpleModule(self.radical())

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: S = SymmetricGroupAlgebra(QQ, 3)
                sage: C = S.cellular_basis()
                sage: W = S.cell_module([2,1])
                sage: elt = W.an_element(); elt
                2*W[[1, 2], [3]] + 2*W[[1, 3], [2]]
                sage: sc = C.an_element(); sc
                2*C([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]])
                 + 3*C([2, 1], [[1, 3], [2]], [[1, 2], [3]])
                 + 2*C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
                sage: sc * elt
                10*W[[1, 3], [2]]

                sage: s = S.an_element(); s
                [1, 2, 3] + 2*[1, 3, 2] + 3*[2, 1, 3] + [3, 1, 2]
                sage: s * elt
                11*W[[1, 2], [3]] - 3/2*W[[1, 3], [2]]

                sage: 1/2 * elt
                W[[1, 2], [3]] + W[[1, 3], [2]]
            """
            # Check for elements coercable to the base ring first
            ret = CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)
            if ret is not None:
                return ret
            # See message in CombinatorialFreeModuleElement._acted_upon_
            P = self.parent()
            if isinstance(scalar, Element) and scalar.parent() is not P._algebra:
                # Temporary needed by coercion (see Polynomial/FractionField tests).
                if not P._algebra.has_coerce_map_from(scalar.parent()):
                    return None
                scalar = P._algebra(scalar)

            if self_on_left:
                raise NotImplementedError
                # scalar = scalar.cellular_involution()
            mc = self._monomial_coefficients
            scalar_mc = scalar.monomial_coefficients(copy=False)
            D = linear_combination([(P._action_basis(x, k)._monomial_coefficients,
                                     scalar_mc[x] * mc[k])
                                    for k in mc for x in scalar_mc],
                                   factor_on_left=False)

            return P._from_dict(D, remove_zeros=False)

        # For backward compatibility
        _lmul_ = _acted_upon_
        _rmul_ = _acted_upon_


class SimpleModule(QuotientModuleWithBasis):
    r"""
    A simple module of a cellular algebra.

    Let `W(\lambda)` denote a cell module. The simple module `L(\lambda)`
    is defined as `W(\lambda) / \operatorname{rad}(\lambda)`,
    where `\operatorname{rad}(\lambda)` is the radical of the
    bilinear form `\Phi_{\lambda}`.
    """
    def __init__(self, submodule):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: L = S.cell_module([2,1]).simple_module()
            sage: TestSuite(L).run()
        """
        cat = ModulesWithBasis(submodule.category().base_ring()).Quotients()
        QuotientModuleWithBasis.__init__(self, submodule, cat)
        # We set some print options since QuotientModuleWithBasis doesn't take them
        self.print_options(prefix='L', bracket=False)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: S.cell_module([2,1]).simple_module()
            Simple module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field
        """
        W = self.ambient()
        return "Simple module indexed by {} of {}".format(W._la, W._algebra)

    def _coerce_map_from_(self, A):
        """
        Return a coercion map from ``A`` if one exists.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: L = W.simple_module()
            sage: L._coerce_map_from_(W)
            Generic morphism:
              From: Cell module indexed by [2, 1] of Cellular basis
                 of Symmetric group algebra of order 3 over Rational Field
              To:   Simple module indexed by [2, 1] of Cellular basis
                 of Symmetric group algebra of order 3 over Rational Field
        """
        if A == self._ambient:
            return A.module_morphism(self.retract, codomain=self)
        return super()._coerce_map_from_(A)

    class Element(QuotientModuleWithBasis.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: S = SymmetricGroupAlgebra(QQ, 3)
                sage: C = S.cellular_basis()
                sage: L = S.cell_module([2,1]).simple_module()
                sage: elt = L.an_element(); elt
                2*L[[1, 2], [3]] + 2*L[[1, 3], [2]]
                sage: sc = C.an_element(); sc
                2*C([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]])
                 + 3*C([2, 1], [[1, 3], [2]], [[1, 2], [3]])
                 + 2*C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
                sage: sc * elt
                10*L[[1, 3], [2]]

                sage: s = S.an_element(); s
                [1, 2, 3] + 2*[1, 3, 2] + 3*[2, 1, 3] + [3, 1, 2]
                sage: s * elt
                11*L[[1, 2], [3]] - 3/2*L[[1, 3], [2]]

                sage: 1/2 * elt
                L[[1, 2], [3]] + L[[1, 3], [2]]
            """
            ret = CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)
            if ret is not None:
                return ret
            P = self.parent()
            ret = P.lift(self)._acted_upon_(scalar, self_on_left)
            if ret is None:
                return None
            return P.retract(ret)

        # For backward compatibility
        _lmul_ = _acted_upon_
        _rmul_ = _acted_upon_
