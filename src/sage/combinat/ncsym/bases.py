r"""
Bases for `NCSym`

AUTHORS:

- Travis Scrimshaw (08-04-2013): Initial version
"""
# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.all import ModulesWithBasis, tensor, Hom
from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.combinat.free_module import CombinatorialFreeModule


class NCSymBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for a basis of `NCSym` or its dual.
    """
    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        INPUT:

        - ``x`` -- a set partition or list of lists of integers

        EXAMPLES::

            sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
            sage: m([[1,3],[2]])
            m{{1, 3}, {2}}
            sage: m(SetPartition([[1,3],[2]]))
            m{{1, 3}, {2}}
        """
        if isinstance(x, (list, tuple)):
            x = SetPartition(x)
        return super(NCSymBasis_abstract, self)._element_constructor_(x)


class NCSymOrNCSymDualBases(Category_realization_of_parent):
    r"""
    Base category for the category of bases of symmetric functions
    in non-commuting variables or its Hopf dual for the common code.
    """
    def super_categories(self):
        r"""
        Return the super categories of bases of (the Hopf dual of) the
        symmetric functions in non-commuting variables.

        OUTPUT:

        - a list of categories

        TESTS::

            sage: from sage.combinat.ncsym.bases import NCSymOrNCSymDualBases
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: NCSymOrNCSymDualBases(NCSym).super_categories()
            [Category of realizations of Symmetric functions in
              non-commuting variables over the Rational Field,
             Category of graded hopf algebras with basis over Rational Field,
             Join of Category of realizations of hopf algebras over Rational Field
              and Category of graded algebras over Rational Field
              and Category of graded coalgebras over Rational Field]
        """
        R = self.base().base_ring()
        from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
        return [self.base().Realizations(),
                GradedHopfAlgebrasWithBasis(R),
                GradedHopfAlgebras(R).Realizations()]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.ncsym.bases import NCSymOrNCSymDualBases
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: NCSymOrNCSymDualBases(NCSym)
            Category of bases of NCSym or NCSym^* over the Rational Field
        """
        return "Category of bases of NCSym or NCSym^* over the {}".format(self.base().base_ring())

    class ParentMethods:
        def _repr_(self):
            """
            TESTS::

                sage: SymmetricFunctionsNonCommutingVariables(QQ).m()
                Symmetric functions in non-commuting variables over the Rational Field in the monomial basis
                sage: SymmetricFunctionsNonCommutingVariables(QQ).m().dual_basis()
                Dual symmetric functions in non-commuting variables over the Rational Field in the w basis
                sage: SymmetricFunctionsNonCommutingVariables(QQ).chi()
                Symmetric functions in non-commuting variables over the Rational Field in the
                 supercharacter basis with parameter q=2
                sage: SymmetricFunctionsNonCommutingVariables(QQ['q'].fraction_field()).rho('q')
                Symmetric functions in non-commuting variables over the Fraction Field
                 of Univariate Polynomial Ring in q over Rational Field in the
                 deformed_coarse_powersum basis with parameter q
            """
            str = "{} in the {} basis".format(self.realization_of(), self._realization_name())
            if hasattr(self, '_q'):
                str += " with parameter q"
                if repr(self._q) != 'q':
                    str += "=" + repr(self._q)
            return str

        def __getitem__(self, i):
            """
            Return the basis element indexed by ``i``.

            INPUT:

            - ``i`` -- a set partition or a list of list of integers

            EXAMPLES::

                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w[[1], [2,3]]
                w{{1}, {2, 3}}
                sage: w[{1}, (2,3)]
                w{{1}, {2, 3}}
                sage: w[[]]
                w{}
            """
            if isinstance(i, SetPartition):
                return self.monomial(i)
            if i == []:
                return self.one()
            if not isinstance(i, tuple):
                i = (i,)
            return self.monomial(SetPartition(i))

        @cached_method
        def one_basis(self):
            r"""
            Return the index of the basis element containing `1`.

            OUTPUT:

            - The empty set partition

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m.one_basis()
                {}
                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w.one_basis()
                {}
            """
            return SetPartition([])

        def counit_on_basis(self, A):
            r"""
            The counit is defined by sending all elements of positive degree
            to zero.

            INPUT:

            - ``A`` -- a set partition

            OUTPUT:

            - either the ``0`` or the ``1`` of the base ring of ``self``

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m.counit_on_basis(SetPartition([[1,3], [2]]))
                0
                sage: m.counit_on_basis(SetPartition([]))
                1
                sage: w = SymmetricFunctionsNonCommutingVariables(QQ).dual().w()
                sage: w.counit_on_basis(SetPartition([[1,3], [2]]))
                0
                sage: w.counit_on_basis(SetPartition([]))
                1
            """
            if len(A) != 0:
                return self.base_ring().zero()
            return self.base_ring().one()

        def duality_pairing(self, x, y):
            r"""
            Compute the pairing between an element of ``self`` and an element
            of the dual.

            Carry out this computation by converting ``x`` to the `\mathbf{m}`
            basis and ``y`` to the `\mathbf{w}` basis.

            INPUT:

            - ``x`` -- an element of symmetric functions in non-commuting
              variables
            - ``y`` -- an element of the dual of symmetric functions in
              non-commuting variables

            OUTPUT:

            - an element of the base ring of ``self``

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: h = NCSym.h()
                sage: w = NCSym.m().dual_basis()
                sage: matrix([[h(A).duality_pairing(w(B)) for A in SetPartitions(3)] for B in SetPartitions(3)])
                [6 2 2 2 1]
                [2 2 1 1 1]
                [2 1 2 1 1]
                [2 1 1 2 1]
                [1 1 1 1 1]
                sage: (h[[1,2],[3]] + 3*h[[1,3],[2]]).duality_pairing(2*w[[1,3],[2]] + w[[1,2,3]] + 2*w[[1,2],[3]])
                32
            """
            m = self.realization_of().m()
            x = m(x)
            y = m.dual_basis()(y)
            return sum(coeff * y[I] for (I, coeff) in x)

        def duality_pairing_matrix(self, basis, degree):
            r"""
            The matrix of scalar products between elements of `NCSym` and
            elements of `NCSym^*`.

            INPUT:

            - ``basis`` -- a basis of the dual Hopf algebra
            - ``degree`` -- a non-negative integer

            OUTPUT:

            - the matrix of scalar products between the basis ``self`` and the
              basis ``basis`` in the dual Hopf algebra of degree ``degree``

            EXAMPLES:

            The matrix between the `\mathbf{m}` basis and the
            `\mathbf{w}` basis::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: m = NCSym.m()
                sage: w = NCSym.dual().w()
                sage: m.duality_pairing_matrix(w, 3)
                [1 0 0 0 0]
                [0 1 0 0 0]
                [0 0 1 0 0]
                [0 0 0 1 0]
                [0 0 0 0 1]

            Similarly for some of the other basis of `NCSym` and the `\mathbf{w}`
            basis::

                sage: e = NCSym.e()
                sage: e.duality_pairing_matrix(w, 3)
                [0 0 0 0 1]
                [0 0 1 1 1]
                [0 1 0 1 1]
                [0 1 1 0 1]
                [1 1 1 1 1]
                sage: p = NCSym.p()
                sage: p.duality_pairing_matrix(w, 3)
                [1 0 0 0 0]
                [1 1 0 0 0]
                [1 0 1 0 0]
                [1 0 0 1 0]
                [1 1 1 1 1]
                sage: cp = NCSym.cp()
                sage: cp.duality_pairing_matrix(w, 3)
                [1 0 0 0 0]
                [1 1 0 0 0]
                [0 0 1 0 0]
                [1 0 0 1 0]
                [1 1 1 1 1]
                sage: x = NCSym.x()
                sage: w.duality_pairing_matrix(x, 3)
                [ 0  0  0  0  1]
                [ 1  0 -1 -1  1]
                [ 1 -1  0 -1  1]
                [ 1 -1 -1  0  1]
                [ 2 -1 -1 -1  1]

            A base case test::

                sage: m.duality_pairing_matrix(w, 0)
                [1]
            """
            from sage.matrix.constructor import matrix
            # TODO: generalize to keys indexing the basis of the graded component
            return matrix(self.base_ring(),
                    [[self.duality_pairing(self[I], basis[J])
                      for J in SetPartitions(degree)]
                     for I in SetPartitions(degree)])

    class ElementMethods:
        def duality_pairing(self, other):
            r"""
            Compute the pairing between ``self`` and an element ``other`` of the dual.

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: m = NCSym.m()
                sage: w = m.dual_basis()
                sage: elt = m[[1,3],[2]] - 3*m[[1,2],[3]]
                sage: elt.duality_pairing(w[[1,3],[2]])
                1
                sage: elt.duality_pairing(w[[1,2],[3]])
                -3
                sage: elt.duality_pairing(w[[1,2]])
                0
                sage: e = NCSym.e()
                sage: w[[1,3],[2]].duality_pairing(e[[1,3],[2]])
                0
            """
            return self.parent().duality_pairing(self, other)


class NCSymBases(Category_realization_of_parent):
    r"""
    Category of bases of symmetric functions in non-commuting variables.

    EXAMPLES::

        sage: from sage.combinat.ncsym.bases import NCSymBases
        sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
        sage: NCSymBases(NCSym)
        Category of bases of symmetric functions in non-commuting variables over the Rational Field
    """
    def super_categories(self):
        r"""
        Return the super categories of bases of the Hopf dual of the
        symmetric functions in non-commuting variables.

        OUTPUT:

        - a list of categories

        TESTS::

            sage: from sage.combinat.ncsym.bases import NCSymBases
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: NCSymBases(NCSym).super_categories()
            [Category of bases of NCSym or NCSym^* over the Rational Field]
        """
        return [NCSymOrNCSymDualBases(self.base())]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.ncsym.bases import NCSymBases
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: NCSymBases(NCSym)
            Category of bases of symmetric functions in non-commuting variables over the Rational Field
        """
        return "Category of bases of symmetric functions in non-commuting"\
               " variables over the {}".format(self.base().base_ring())

    class ParentMethods:
        def from_symmetric_function(self, f):
            """
            Return the image of the symmetric function ``f`` in ``self``.

            This is performed by converting to the monomial basis and
            extending the method :meth:`sum_of_partitions` linearly.  This is a
            linear map from the symmetric functions to the symmetric functions
            in non-commuting variables that does not preserve the product or
            coproduct structure of the Hopf algebra.

            .. SEEALSO:: :meth:`to_symmetric_function`

            INPUT:

            - ``f`` -- a symmetric function

            OUTPUT:

            - an element of ``self``

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: Sym = SymmetricFunctions(QQ)
                sage: e = NCSym.e()
                sage: elem = Sym.e()
                sage: elt = e.from_symmetric_function(elem[2,1,1]); elt
                1/12*e{{1}, {2}, {3, 4}} + 1/12*e{{1}, {2, 3}, {4}} + 1/12*e{{1}, {2, 4}, {3}}
                 + 1/12*e{{1, 2}, {3}, {4}} + 1/12*e{{1, 3}, {2}, {4}} + 1/12*e{{1, 4}, {2}, {3}}
                sage: elem(elt.to_symmetric_function())
                e[2, 1, 1]
                sage: e.from_symmetric_function(elem[4])
                1/24*e{{1, 2, 3, 4}}
                sage: p = NCSym.p()
                sage: pow = Sym.p()
                sage: elt = p.from_symmetric_function(pow[2,1,1]); elt
                1/6*p{{1}, {2}, {3, 4}} + 1/6*p{{1}, {2, 3}, {4}} + 1/6*p{{1}, {2, 4}, {3}}
                 + 1/6*p{{1, 2}, {3}, {4}} + 1/6*p{{1, 3}, {2}, {4}} + 1/6*p{{1, 4}, {2}, {3}}
                sage: pow(elt.to_symmetric_function())
                p[2, 1, 1]
                sage: p.from_symmetric_function(pow[4])
                p{{1, 2, 3, 4}}
                sage: h = NCSym.h()
                sage: comp = Sym.complete()
                sage: elt = h.from_symmetric_function(comp[2,1,1]); elt
                1/12*h{{1}, {2}, {3, 4}} + 1/12*h{{1}, {2, 3}, {4}} + 1/12*h{{1}, {2, 4}, {3}}
                 + 1/12*h{{1, 2}, {3}, {4}} + 1/12*h{{1, 3}, {2}, {4}} + 1/12*h{{1, 4}, {2}, {3}}
                sage: comp(elt.to_symmetric_function())
                h[2, 1, 1]
                sage: h.from_symmetric_function(comp[4])
                1/24*h{{1, 2, 3, 4}}
            """
            m = self.realization_of().m()
            return self(m.from_symmetric_function(f))

        def primitive(self, A, i=1):
            r"""
            Return the primitive associated to ``A`` in ``self``.

            .. SEEALSO::

                :meth:`~sage.combinat.ncsym.ncsym.SymmetricFunctionsNonCommutingVariables.powersum.primitive`

            INPUT:

            - ``A`` -- a set partition
            - ``i`` -- a positive integer

            OUTPUT:

            - an element of ``self``

            EXAMPLES::

                sage: e = SymmetricFunctionsNonCommutingVariables(QQ).e()
                sage: elt = e.primitive(SetPartition([[1,3],[2]])); elt
                e{{1, 2}, {3}} - e{{1, 3}, {2}}
                sage: elt.coproduct()
                e{} # e{{1, 2}, {3}} - e{} # e{{1, 3}, {2}} + e{{1, 2}, {3}} # e{} - e{{1, 3}, {2}} # e{}
            """
            p = self.realization_of().p()
            return self(p.primitive(A, i))

        @abstract_method(optional=True)
        def internal_coproduct_on_basis(self, i):
            """
            The internal coproduct of the algebra on the basis (optional).

            INPUT:

            - ``i`` -- the indices of an element of the basis of ``self``

            OUTPUT:

            - an element of the tensor squared of ``self``

            EXAMPLES::

                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).m()
                sage: m.internal_coproduct_on_basis(SetPartition([[1,2]]))
                m{{1, 2}} # m{{1, 2}}
            """

        @lazy_attribute
        def internal_coproduct(self):
            r"""
            Compute the internal coproduct of ``self``.

            If :meth:`internal_coproduct_on_basis()` is available, construct
            the internal coproduct morphism from ``self`` to ``self``
            `\otimes` ``self`` by extending it by linearity. Otherwise, this uses
            :meth:`internal_coproduct_by_coercion()`, if available.

            OUTPUT:

            - an element of the tensor squared of ``self``

            EXAMPLES::

                sage: cp = SymmetricFunctionsNonCommutingVariables(QQ).cp()
                sage: cp.internal_coproduct(cp[[1,3],[2]] - 2*cp[[1]])
                -2*cp{{1}} # cp{{1}} + cp{{1, 2, 3}} # cp{{1, 3}, {2}} + cp{{1, 3}, {2}} # cp{{1, 2, 3}}
                 + cp{{1, 3}, {2}} # cp{{1, 3}, {2}}
            """
            if self.internal_coproduct_on_basis is not NotImplemented:
                return Hom(self, tensor([self, self]),
                           ModulesWithBasis(self.base_ring()))(on_basis=self.internal_coproduct_on_basis)
            elif hasattr(self, "internal_coproduct_by_coercion"):
                return self.internal_coproduct_by_coercion

        def internal_coproduct_by_coercion(self, x):
            r"""
            Return the internal coproduct by coercing the element to the powersum basis.

            INPUT:

            - ``x`` -- an element of ``self``

            OUTPUT:

            - an element of the tensor squared of ``self``

            EXAMPLES::

                sage: h = SymmetricFunctionsNonCommutingVariables(QQ).h()
                sage: h[[1,3],[2]].internal_coproduct() # indirect doctest
                2*h{{1}, {2}, {3}} # h{{1}, {2}, {3}} - h{{1}, {2}, {3}} # h{{1, 3}, {2}}
                 - h{{1, 3}, {2}} # h{{1}, {2}, {3}} + h{{1, 3}, {2}} # h{{1, 3}, {2}}
            """
            R = self.realization_of().a_realization()
            return self.tensor_square().sum(coeff * tensor([self(R[A]), self(R[B])])
                                            for ((A, B), coeff) in R(x).internal_coproduct())

    class ElementMethods:
        def expand(self, n, alphabet='x'):
            r"""
            Expand the symmetric function into ``n`` non-commuting
            variables in an alphabet, which by default is ``'x'``.

            This computation is completed by coercing the element ``self``
            into the monomial basis and computing the expansion in
            the ``alphabet`` there.

            INPUT:

            - ``n`` -- the number of variables in the expansion
            - ``alphabet`` -- (default: ``'x'``) the alphabet in which
              ``self`` is to be expanded

            OUTPUT:

            - an expansion of ``self`` into the ``n`` non-commuting
              variables specified by ``alphabet``

            EXAMPLES::

                sage: h = SymmetricFunctionsNonCommutingVariables(QQ).h()
                sage: h[[1,3],[2]].expand(3)
                2*x0^3 + x0^2*x1 + x0^2*x2 + 2*x0*x1*x0 + x0*x1^2 + x0*x1*x2 + 2*x0*x2*x0
                 + x0*x2*x1 + x0*x2^2 + x1*x0^2 + 2*x1*x0*x1 + x1*x0*x2 + x1^2*x0 + 2*x1^3
                 + x1^2*x2 + x1*x2*x0 + 2*x1*x2*x1 + x1*x2^2 + x2*x0^2 + x2*x0*x1 + 2*x2*x0*x2
                 + x2*x1*x0 + x2*x1^2 + 2*x2*x1*x2 + x2^2*x0 + x2^2*x1 + 2*x2^3
                sage: x = SymmetricFunctionsNonCommutingVariables(QQ).x()
                sage: x[[1,3],[2]].expand(3)
                -x0^2*x1 - x0^2*x2 - x0*x1^2 - x0*x1*x2 - x0*x2*x1 - x0*x2^2 - x1*x0^2
                 - x1*x0*x2 - x1^2*x0 - x1^2*x2 - x1*x2*x0 - x1*x2^2 - x2*x0^2 - x2*x0*x1
                 - x2*x1*x0 - x2*x1^2 - x2^2*x0 - x2^2*x1
            """
            m = self.parent().realization_of().monomial()
            return m(self).expand(n, alphabet)

        def to_symmetric_function(self):
            r"""
            Compute the projection of an element of symmetric function in
            non-commuting variables to the symmetric functions.

            The projection of a monomial symmetric function in non-commuting
            variables indexed by the set partition ``A`` is defined as

            .. MATH::

                \mathbf{m}_A \mapsto m_{\lambda(A)} \prod_i n_i(\lambda(A))!

            where `\lambda(A)` is the partition associated with `A` by
            taking the sizes of the parts and `n_i(\mu)` is the
            multiplicity of `i` in `\mu`.  For other bases this map is extended
            linearly.

            OUTPUT:

            - an element of the symmetric functions in the monomial basis

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: e = NCSym.e()
                sage: h = NCSym.h()
                sage: p = NCSym.p()
                sage: cp = NCSym.cp()
                sage: x = NCSym.x()
                sage: cp[[1,3],[2]].to_symmetric_function()
                m[2, 1]
                sage: x[[1,3],[2]].to_symmetric_function()
                -6*m[1, 1, 1] - 2*m[2, 1]
                sage: e[[1,3],[2]].to_symmetric_function()
                2*e[2, 1]
                sage: h[[1,3],[2]].to_symmetric_function()
                2*h[2, 1]
                sage: p[[1,3],[2]].to_symmetric_function()
                p[2, 1]
            """
            m = self.parent().realization_of().monomial()
            return m(self).to_symmetric_function()

        def to_wqsym(self):
            r"""
            Return the image of ``self`` under the canonical
            inclusion map `NCSym \to WQSym`.

            The canonical inclusion map `NCSym \to WQSym` is
            an injective homomorphism of algebras. It sends a
            basis element `\mathbf{m}_A` of `NCSym` to the sum of
            basis elements `\mathbf{M}_P` of `WQSym`, where `P`
            ranges over all ordered set partitions that become
            `A` when the ordering is forgotten.
            This map is denoted by `\theta` in [BZ05]_ (17).

            .. SEEALSO::

                :class:`WordQuasiSymmetricFunctions` for a
                definition of `WQSym`.

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: e = NCSym.e()
                sage: h = NCSym.h()
                sage: p = NCSym.p()
                sage: cp = NCSym.cp()
                sage: x = NCSym.x()
                sage: m = NCSym.m()
                sage: m[[1,3],[2]].to_wqsym()
                M[{1, 3}, {2}] + M[{2}, {1, 3}]
                sage: x[[1,3],[2]].to_wqsym()
                -M[{1}, {2}, {3}] - M[{1}, {2, 3}] - M[{1}, {3}, {2}]
                 - M[{1, 2}, {3}] - M[{2}, {1}, {3}] - M[{2}, {3}, {1}]
                 - M[{2, 3}, {1}] - M[{3}, {1}, {2}] - M[{3}, {1, 2}]
                 - M[{3}, {2}, {1}]
                sage: (4*p[[1,3],[2]]-p[[1]]).to_wqsym()
                -M[{1}] + 4*M[{1, 2, 3}] + 4*M[{1, 3}, {2}] + 4*M[{2}, {1, 3}]
            """
            parent = self.parent()
            NCSym = parent.realization_of()
            R = parent.base_ring()
            m = NCSym.monomial()
            from sage.combinat.chas.wqsym import WordQuasiSymmetricFunctions
            M = WordQuasiSymmetricFunctions(R).M()
            from itertools import permutations
            OSP = M.basis().keys()

            def to_wqsym_on_m_basis(A):
                # Return the image of `\mathbf{m}_A` under the inclusion
                # map `NCSym \to WQSym`.
                l = len(A)
                return M.sum_of_terms(((OSP([A[ui] for ui in u]), 1)
                                       for u in permutations(range(l))),
                                      distinct=True)
            return M.linear_combination((to_wqsym_on_m_basis(A), coeff)
                                        for A, coeff in m(self))

        def internal_coproduct(self):
            r"""
            Return the internal coproduct of ``self``.

            The internal coproduct is defined on the power sum basis as

            .. MATH::

                \mathbf{p}_A \mapsto \mathbf{p}_A \otimes \mathbf{p}_A

            and the map is extended linearly.

            OUTPUT:

            - an element of the tensor square of the basis of ``self``

            EXAMPLES::

                sage: x = SymmetricFunctionsNonCommutingVariables(QQ).x()
                sage: x[[1,3],[2]].internal_coproduct()
                x{{1}, {2}, {3}} # x{{1, 3}, {2}} + x{{1, 3}, {2}} # x{{1}, {2}, {3}}
                 + x{{1, 3}, {2}} # x{{1, 3}, {2}}
            """
            return self.parent().internal_coproduct(self)

        def omega(self):
            r"""
            Return the involution `\omega` applied to ``self``.

            The involution `\omega` is defined by

            .. MATH::

                \mathbf{e}_A \mapsto \mathbf{h}_A

            and the result is extended linearly.

            OUTPUT:

            - an element in the same basis as ``self``

            EXAMPLES::

                sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
                sage: m = NCSym.m()
                sage: m[[1,3],[2]].omega()
                -2*m{{1, 2, 3}} - m{{1, 3}, {2}}
                sage: p = NCSym.p()
                sage: p[[1,3],[2]].omega()
                -p{{1, 3}, {2}}
                sage: cp = NCSym.cp()
                sage: cp[[1,3],[2]].omega()
                -2*cp{{1, 2, 3}} - cp{{1, 3}, {2}}
                sage: x = NCSym.x()
                sage: x[[1,3],[2]].omega()
                -2*x{{1}, {2}, {3}} - x{{1, 3}, {2}}
            """
            P = self.parent()
            e = P.realization_of().e()
            h = P.realization_of().h()
            return P(h.sum_of_terms(e(self)))


class MultiplicativeNCSymBases(Category_realization_of_parent):
    r"""
    Category of multiplicative bases of symmetric functions in non-commuting variables.

    A multiplicative basis is one for which `\mathbf{b}_A \mathbf{b}_B = \mathbf{b}_{A|B}`
    where `A|B` is the :meth:`~sage.combinat.set_partition.SetPartition.pipe` operation
    on set partitions.

    EXAMPLES::

        sage: from sage.combinat.ncsym.bases import MultiplicativeNCSymBases
        sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
        sage: MultiplicativeNCSymBases(NCSym)
        Category of multiplicative bases of symmetric functions in non-commuting variables over the Rational Field
    """
    def super_categories(self):
        r"""
        Return the super categories of bases of the Hopf dual of the
        symmetric functions in non-commuting variables.

        OUTPUT:

        - a list of categories

        TESTS::

            sage: from sage.combinat.ncsym.bases import MultiplicativeNCSymBases
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: MultiplicativeNCSymBases(NCSym).super_categories()
            [Category of bases of symmetric functions in non-commuting variables over the Rational Field]
        """
        return [NCSymBases(self.base())]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.ncsym.bases import MultiplicativeNCSymBases
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: MultiplicativeNCSymBases(NCSym)
            Category of multiplicative bases of symmetric functions in non-commuting variables over the Rational Field
        """
        return "Category of multiplicative bases of symmetric functions in non-commuting"\
               " variables over the {}".format(self.base().base_ring())

    class ParentMethods:
        def product_on_basis(self, A, B):
            r"""
            The product on basis elements.

            The product on a multiplicative basis is given by
            `\mathbf{b}_A \cdot \mathbf{b}_B = \mathbf{b}_{A | B}`.

            The bases `\{ \mathbf{e}, \mathbf{h}, \mathbf{x}, \mathbf{cp}, \mathbf{p},
            \mathbf{chi}, \mathbf{rho} \}` are all multiplicative.

            INPUT:

            - ``A``, ``B`` -- set partitions

            OUTPUT:

            - an element in the basis ``self``

            EXAMPLES::

                sage: e = SymmetricFunctionsNonCommutingVariables(QQ).e()
                sage: h = SymmetricFunctionsNonCommutingVariables(QQ).h()
                sage: x = SymmetricFunctionsNonCommutingVariables(QQ).x()
                sage: cp = SymmetricFunctionsNonCommutingVariables(QQ).cp()
                sage: p = SymmetricFunctionsNonCommutingVariables(QQ).p()
                sage: chi = SymmetricFunctionsNonCommutingVariables(QQ).chi()
                sage: rho = SymmetricFunctionsNonCommutingVariables(QQ).rho()
                sage: A = SetPartition([[1], [2, 3]])
                sage: B = SetPartition([[1], [3], [2,4]])
                sage: e.product_on_basis(A, B)
                e{{1}, {2, 3}, {4}, {5, 7}, {6}}
                sage: h.product_on_basis(A, B)
                h{{1}, {2, 3}, {4}, {5, 7}, {6}}
                sage: x.product_on_basis(A, B)
                x{{1}, {2, 3}, {4}, {5, 7}, {6}}
                sage: cp.product_on_basis(A, B)
                cp{{1}, {2, 3}, {4}, {5, 7}, {6}}
                sage: p.product_on_basis(A, B)
                p{{1}, {2, 3}, {4}, {5, 7}, {6}}
                sage: chi.product_on_basis(A, B)
                chi{{1}, {2, 3}, {4}, {5, 7}, {6}}
                sage: rho.product_on_basis(A, B)
                rho{{1}, {2, 3}, {4}, {5, 7}, {6}}
                sage: e.product_on_basis(A,B)==e(h(e(A))*h(e(B)))
                True
                sage: h.product_on_basis(A,B)==h(x(h(A))*x(h(B)))
                True
                sage: x.product_on_basis(A,B)==x(h(x(A))*h(x(B)))
                True
                sage: cp.product_on_basis(A,B)==cp(p(cp(A))*p(cp(B)))
                True
                sage: p.product_on_basis(A,B)==p(e(p(A))*e(p(B)))
                True
            """
            return self.monomial(A.pipe(B))

    class ElementMethods:
        pass


class NCSymDualBases(Category_realization_of_parent):
    r"""
    Category of bases of dual symmetric functions in non-commuting variables.

    EXAMPLES::

        sage: from sage.combinat.ncsym.bases import NCSymDualBases
        sage: DNCSym = SymmetricFunctionsNonCommutingVariables(QQ).dual()
        sage: NCSymDualBases(DNCSym)
        Category of bases of dual symmetric functions in non-commuting variables over the Rational Field
    """
    def super_categories(self):
        r"""
        Return the super categories of bases of the Hopf dual of the
        symmetric functions in non-commuting variables.

        OUTPUT:

        - a list of categories

        TESTS::

            sage: from sage.combinat.ncsym.bases import NCSymBases
            sage: NCSym = SymmetricFunctionsNonCommutingVariables(QQ)
            sage: NCSymBases(NCSym).super_categories()
            [Category of bases of NCSym or NCSym^* over the Rational Field]
        """
        return [NCSymOrNCSymDualBases(self.base())]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.ncsym.bases import NCSymDualBases
            sage: DNCSym = SymmetricFunctionsNonCommutingVariables(QQ).dual()
            sage: NCSymDualBases(DNCSym)
            Category of bases of dual symmetric functions in non-commuting variables over the Rational Field
        """
        return "Category of bases of dual symmetric functions in non-commuting"\
               " variables over the {}".format(self.base().base_ring())
