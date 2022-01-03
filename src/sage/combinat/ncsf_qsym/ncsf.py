# -*- coding: utf-8 -*-
"""
Non-Commutative Symmetric Functions
"""
#*****************************************************************************
#       Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#                     2012 Franco Saliola <saliola@gmail.com>,
#                     2012 Chris Berg <chrisjamesberg@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

########################################
# TODO:
# 1. Make Coercion run faster between multiple bases.
########################################

from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.arith.misc import factorial
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.rings import Rings
from sage.categories.fields import Fields
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.combinat.composition import Compositions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.ncsf_qsym.generic_basis_code import BasesOfQSymOrNCSF
from sage.combinat.ncsf_qsym.combinatorics import (coeff_pi, coeff_lp,
        coeff_sp, coeff_ell, m_to_s_stat, number_of_fCT, number_of_SSRCT, compositions_order)
from sage.combinat.partition import Partition
from sage.combinat.permutation import Permutations
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.combinat.sf.sf import SymmetricFunctions


class NonCommutativeSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The abstract algebra of non-commutative symmetric functions.

    We construct the abstract algebra of non-commutative symmetric
    functions over the rational numbers::

        sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
        sage: NCSF
        Non-Commutative Symmetric Functions over the Rational Field
        sage: S = NCSF.complete()
        sage: R = NCSF.ribbon()
        sage: S[2,1]*R[1,2]
        S[2, 1, 1, 2] - S[2, 1, 3]

    NCSF is the unique free (non-commutative!) graded connected algebra with
    one generator in each degree::

        sage: NCSF.category()
        Join of Category of hopf algebras over Rational Field
            and Category of graded algebras over Rational Field
            and Category of monoids with realizations
            and Category of graded coalgebras over Rational Field
            and Category of coalgebras over Rational Field with realizations
            and Category of cocommutative coalgebras over Rational Field

        sage: [S[i].degree() for i in range(10)]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    We use the Sage standard renaming idiom to get shorter outputs::

        sage: NCSF.rename("NCSF")
        sage: NCSF
        NCSF

    NCSF has many representations as a concrete algebra. Each of them
    has a distinguished basis, and its elements are expanded in this
    basis. Here is the `\Psi`
    (:class:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Psi`)
    representation::

        sage: Psi = NCSF.Psi()
        sage: Psi
        NCSF in the Psi basis

    Elements of ``Psi`` are linear combinations of basis elements indexed
    by compositions::

      sage: Psi.an_element()
      2*Psi[] + 2*Psi[1] + 3*Psi[1, 1]

    The basis itself is accessible through::

        sage: Psi.basis()
        Lazy family (Term map from Compositions of non-negative integers...
        sage: Psi.basis().keys()
        Compositions of non-negative integers

    To construct an element one can therefore do::

        sage: Psi.basis()[Composition([2,1,3])]
        Psi[2, 1, 3]

    As this is rather cumbersome, the following abuses of notation are
    allowed::

        sage: Psi[Composition([2, 1, 3])]
        Psi[2, 1, 3]
        sage: Psi[[2, 1, 3]]
        Psi[2, 1, 3]
        sage: Psi[2, 1, 3]
        Psi[2, 1, 3]

    or even::

        sage: Psi[(i for i in [2, 1, 3])]
        Psi[2, 1, 3]

    Unfortunately, due to a limitation in Python syntax, one cannot use::

        sage: Psi[]       # not implemented

    Instead, you can use::

        sage: Psi[[]]
        Psi[]

    Now, we can construct linear combinations of basis elements::

        sage: Psi[2,1,3] + 2 * (Psi[4] + Psi[2,1])
        2*Psi[2, 1] + Psi[2, 1, 3] + 2*Psi[4]

    .. rubric:: Algebra structure

    To start with, ``Psi`` is a graded algebra, the grading being induced by
    the size of compositions. The one is the basis element indexed by the empty
    composition::

        sage: Psi.one()
        Psi[]
        sage: S.one()
        S[]
        sage: R.one()
        R[]

    As we have seen above, the ``Psi`` basis is multiplicative; that is
    multiplication is induced by linearity from the concatenation of
    compositions::

        sage: Psi[1,3] * Psi[2,1]
        Psi[1, 3, 2, 1]
        sage: (Psi.one() + 2 * Psi[1,3]) * Psi[2, 4]
        2*Psi[1, 3, 2, 4] + Psi[2, 4]

    .. rubric:: Hopf algebra structure

    ``Psi`` is further endowed with a coalgebra structure. The coproduct
    is an algebra morphism, and therefore determined by its values on
    the generators; those are primitive::

        sage: Psi[1].coproduct()
        Psi[] # Psi[1] + Psi[1] # Psi[]
        sage: Psi[2].coproduct()
        Psi[] # Psi[2] + Psi[2] # Psi[]

    The coproduct, being cocommutative on the generators, is
    cocommutative everywhere::

        sage: Psi[1,2].coproduct()
        Psi[] # Psi[1, 2] + Psi[1] # Psi[2] + Psi[1, 2] # Psi[] + Psi[2] # Psi[1]

    The algebra and coalgebra structures on ``Psi`` combine to form a
    bialgebra structure, which cooperates with the grading to form a
    connected graded bialgebra. Thus, as any connected graded bialgebra,
    ``Psi`` is a Hopf algebra. Over ``QQ`` (or any other `\QQ`-algebra),
    this Hopf algebra ``Psi`` is isomorphic to the tensor algebra of
    its space of primitive elements.

    The antipode is an anti-algebra morphism; in the ``Psi`` basis, it
    sends the generators to their opposites and changes their sign if
    they are of odd degree::

        sage: Psi[3].antipode()
        -Psi[3]
        sage: Psi[1,3,2].antipode()
        -Psi[2, 3, 1]
        sage: Psi[1,3,2].coproduct().apply_multilinear_morphism(lambda be,ga: Psi(be)*Psi(ga).antipode())
        0

    The counit is defined by sending all elements of positive degree to
    zero::

        sage: S[3].degree(), S[3,1,2].degree(), S.one().degree()
        (3, 6, 0)
        sage: S[3].counit()
        0
        sage: S[3,1,2].counit()
        0
        sage: S.one().counit()
        1
        sage: (S[3] - 2*S[3,1,2] + 7).counit()
        7
        sage: (R[3] - 2*R[3,1,2] + 7).counit()
        7

    It is possible to change the prefix used to display the basis
    elements using the method
    :meth:`~sage.structure.indexed_generators.IndexedGenerators.print_options`.
    Say that for instance one wanted to display the
    :class:`~NonCommutativeSymmetricFunctions.Complete` basis as having
    a prefix ``H`` instead of the default ``S``::

        sage: H = NCSF.complete()
        sage: H.an_element()
        2*S[] + 2*S[1] + 3*S[1, 1]
        sage: H.print_options(prefix='H')
        sage: H.an_element()
        2*H[] + 2*H[1] + 3*H[1, 1]
        sage: H.print_options(prefix='S') #restore to 'S'

    .. rubric:: Concrete representations

    NCSF admits the concrete realizations defined in [NCSF1]_::

        sage: Phi        = NCSF.Phi()
        sage: Psi        = NCSF.Psi()
        sage: ribbon     = NCSF.ribbon()
        sage: complete   = NCSF.complete()
        sage: elementary = NCSF.elementary()

    To change from one basis to another, one simply does::

        sage: Phi(Psi[1])
        Phi[1]
        sage: Phi(Psi[3])
        -1/4*Phi[1, 2] + 1/4*Phi[2, 1] + Phi[3]

    In general, one can mix up different bases in computations::

        sage: Phi[1] * Psi[1]
        Phi[1, 1]

    Some of the changes of basis are easy to guess::

        sage: ribbon(complete[1,3,2])
        R[1, 3, 2] + R[1, 5] + R[4, 2] + R[6]

    This is the sum of all fatter compositions. Using the usual
    Möbius function for the boolean lattice, the inverse change of
    basis is given by the alternating sum of all fatter compositions::

        sage: complete(ribbon[1,3,2])
        S[1, 3, 2] - S[1, 5] - S[4, 2] + S[6]

    The analogue of the elementary basis is the sum over
    all finer compositions than the 'complement' of the composition
    in the ribbon basis::

        sage: Composition([1,3,2]).complement()
        [2, 1, 2, 1]
        sage: ribbon(elementary([1,3,2]))
        R[1, 1, 1, 1, 1, 1] + R[1, 1, 1, 2, 1] + R[2, 1, 1, 1, 1] + R[2, 1, 2, 1]

    By Möbius inversion on the composition poset, the ribbon
    basis element corresponding to a composition `I` is then the
    alternating sum over all compositions fatter than the
    complement composition of `I` in the elementary basis::

        sage: elementary(ribbon[2,1,2,1])
        L[1, 3, 2] - L[1, 5] - L[4, 2] + L[6]

    The `\Phi`
    (:class:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Phi`)
    and `\Psi` bases are computed by changing to and from the
    :class:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Complete`
    basis. The expansion of `\Psi` basis is given in Proposition 4.5
    of [NCSF1]_ by the formulae

    .. MATH::

        S^I = \sum_{J \geq I} \frac{1}{\pi_u(J,I)} \Psi^J

    and

    .. MATH::

        \Psi^I = \sum_{J \geq I} (-1)^{\ell(J)-\ell(I)} lp(J,I) S^J

    where the coefficients `\pi_u(J,I)` and `lp(J,I)` are coefficients in the
    methods :meth:`~sage.combinat.ncsf_qsym.combinatorics.coeff_pi` and
    :meth:`~sage.combinat.ncsf_qsym.combinatorics.coeff_lp` respectively.  For
    example::

        sage: Psi(complete[3])
        1/6*Psi[1, 1, 1] + 1/3*Psi[1, 2] + 1/6*Psi[2, 1] + 1/3*Psi[3]
        sage: complete(Psi[3])
        S[1, 1, 1] - 2*S[1, 2] - S[2, 1] + 3*S[3]

    The
    :class:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Phi`
    basis is another analogue of the power sum basis from the algebra of
    symmetric functions and the expansion in the Complete basis is given in
    Proposition 4.9 of [NCSF1]_ by the formulae

    .. MATH::

        S^I = \sum_{J \geq I} \frac{1}{sp(J,I)} \Phi^J

    and

    .. MATH::

        \Phi^I = \sum_{J \geq I} (-1)^{\ell(J)-\ell(I)}
        \frac{\prod_i I_i}{\ell(J,I)} S^J

    where the coefficients `sp(J,I)` and `\ell(J,I)` are coefficients in the
    methods :meth:`~sage.combinat.ncsf_qsym.combinatorics.coeff_sp` and
    :meth:`~sage.combinat.ncsf_qsym.combinatorics.coeff_ell` respectively.
    For example::

        sage: Phi(complete[3])
        1/6*Phi[1, 1, 1] + 1/4*Phi[1, 2] + 1/4*Phi[2, 1] + 1/3*Phi[3]
        sage: complete(Phi[3])
        S[1, 1, 1] - 3/2*S[1, 2] - 3/2*S[2, 1] + 3*S[3]

    Here is how to fetch the conversion morphisms::

        sage: f = complete.coerce_map_from(elementary); f
        Generic morphism:
          From: NCSF in the Elementary basis
          To:   NCSF in the Complete basis
        sage: g = elementary.coerce_map_from(complete); g
        Generic morphism:
          From: NCSF in the Complete basis
          To:   NCSF in the Elementary basis
        sage: f.category()
        Category of homsets of unital magmas and right modules over Rational Field and
          left modules over Rational Field
        sage: f(elementary[1,2,2])
        S[1, 1, 1, 1, 1] - S[1, 1, 1, 2] - S[1, 2, 1, 1] + S[1, 2, 2]
        sage: g(complete[1,2,2])
        L[1, 1, 1, 1, 1] - L[1, 1, 1, 2] - L[1, 2, 1, 1] + L[1, 2, 2]
        sage: h = f*g; h
        Composite map:
          From: NCSF in the Complete basis
          To:   NCSF in the Complete basis
          Defn:   Generic morphism:
                  From: NCSF in the Complete basis
                  To:   NCSF in the Elementary basis
                then
                  Generic morphism:
                  From: NCSF in the Elementary basis
                  To:   NCSF in the Complete basis
        sage: h(complete[1,3,2])
        S[1, 3, 2]

    .. rubric:: Additional concrete representations

    NCSF has some additional bases which appear in the literature::

        sage: Monomial                 = NCSF.Monomial()
        sage: Immaculate               = NCSF.Immaculate()
        sage: dualQuasisymmetric_Schur = NCSF.dualQuasisymmetric_Schur()

    The :class:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Monomial`
    basis was introduced in [Tev2007]_ and the
    :class:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Immaculate`
    basis was introduced in [BBSSZ2012]_.  The
    :class:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Quasisymmetric_Schur`
    were defined in [QSCHUR]_ and the dual basis is implemented here as
    :class:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.dualQuasisymmetric_Schur`.
    Refer to the documentation for the use and definition of these bases.

    .. TODO::

        - implement fundamental, forgotten, and simple (coming
          from the simple modules of HS_n) bases.

    We revert back to the original name from our custom short name NCSF::

        sage: NCSF
        NCSF
        sage: NCSF.rename()
        sage: NCSF
        Non-Commutative Symmetric Functions over the Rational Field

    TESTS::

        sage: TestSuite(Phi).run()
        sage: TestSuite(Psi).run()
        sage: TestSuite(complete).run()

    """
    def __init__(self, R):
        r"""
        TESTS::

            sage: NCSF1 = NonCommutativeSymmetricFunctions(FiniteField(23))
            sage: NCSF2 = NonCommutativeSymmetricFunctions(Integers(23))
            sage: TestSuite(NonCommutativeSymmetricFunctions(QQ)).run()
        """
        # change the line below to assert(R in Rings()) once MRO issues from #15536, #15475 are resolved
        assert(R in Fields() or R in Rings())  # side effect of this statement assures MRO exists for R
        self._base = R  # Won't be needed once CategoryObject won't override base_ring
        cat = GradedHopfAlgebras(R).WithRealizations().Cocommutative()
        Parent.__init__(self, category=cat)

        # COERCION METHODS
        Psi = self.Psi()
        Phi = self.Phi()
        complete = self.complete()
        elementary = self.elementary()
        ribbon = self.ribbon()

        # complete to ribbon, and back
        complete.module_morphism(ribbon.sum_of_fatter_compositions,
                                 codomain=ribbon).register_as_coercion()
        ribbon.module_morphism(complete.alternating_sum_of_fatter_compositions,
                               codomain=complete).register_as_coercion()

        complete.algebra_morphism(elementary.alternating_sum_of_compositions,
                                  codomain=elementary).register_as_coercion()
        elementary.algebra_morphism(complete.alternating_sum_of_compositions,
                                    codomain=complete).register_as_coercion()

        complete.algebra_morphism(Psi._from_complete_on_generators,
                                  codomain=Psi).register_as_coercion()
        Psi.algebra_morphism(Psi._to_complete_on_generators,
                             codomain=complete).register_as_coercion()

        complete.algebra_morphism(Phi._from_complete_on_generators,
                                  codomain=Phi).register_as_coercion()
        Phi.algebra_morphism(Phi._to_complete_on_generators,
                             codomain=complete).register_as_coercion()

    def _repr_(self): # could be taken care of by the category
        r"""
        EXAMPLES::

            sage: N = NonCommutativeSymmetricFunctions(ZZ)
            sage: N._repr_()
            'Non-Commutative Symmetric Functions over the Integer Ring'
        """
        return "Non-Commutative Symmetric Functions over the %s" % self.base_ring()

    def a_realization(self):
        r"""
        Gives a realization of the algebra of non-commutative symmetric functions. This
        particular realization is the complete basis of non-commutative symmetric functions.

        OUTPUT:

        - The realization of the non-commutative symmetric functions in the
          complete basis.

        EXAMPLES::

            sage: NonCommutativeSymmetricFunctions(ZZ).a_realization()
            Non-Commutative Symmetric Functions over the Integer Ring in the Complete basis
        """
        return self.complete()

    _shorthands = tuple(['S', 'R', 'L', 'Phi', 'Psi', 'nM', 'I', 'dQS', 'dYQS', 'ZL', 'ZR'])

    def dual(self):
        r"""
        Return the dual to the non-commutative symmetric functions.

        OUTPUT:

        - The dual of the non-commutative symmetric functions over a ring. This
          is the algebra of quasi-symmetric functions over the ring.

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: NCSF.dual()
            Quasisymmetric functions over the Rational Field
        """
        from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
        return QuasiSymmetricFunctions(self.base_ring())

    class Bases(Category_realization_of_parent):
        """
        Category of bases of non-commutative symmetric functions.

        EXAMPLES::

            sage: N = NonCommutativeSymmetricFunctions(QQ)
            sage: N.Bases()
            Category of bases of Non-Commutative Symmetric Functions over the Rational Field
            sage: R = N.Ribbon()
            sage: R in N.Bases()
            True
        """
        def super_categories(self):
            r"""
            Return the super categories of the category of bases of the
            non-commutative symmetric functions.

            OUTPUT:

            - list

            TESTS::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: N.Bases().super_categories()
                [Category of bases of Non-Commutative Symmetric Functions or Quasisymmetric functions over the Rational Field,
                 Category of realizations of graded modules with internal product over Rational Field]

            """
            R = self.base().base_ring()
            from .generic_basis_code import GradedModulesWithInternalProduct
            return [BasesOfQSymOrNCSF(self.base()),
                    GradedModulesWithInternalProduct(R).Realizations()]

        class ParentMethods:

            def to_symmetric_function_on_basis(self, I):
                r"""
                The image of the basis element indexed by ``I`` under the map
                to the symmetric functions.

                This default implementation does a change of basis and
                computes the image in the complete basis.

                INPUT:

                - ``I`` -- a composition

                OUTPUT:

                - The image of the non-commutative basis element of
                  ``self`` indexed by the composition ``I`` under the map from
                  non-commutative symmetric functions to the symmetric
                  functions. This will be a symmetric function.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: I = N.Immaculate()
                    sage: I.to_symmetric_function(I[1,3])
                    -h[2, 2] + h[3, 1]
                    sage: I.to_symmetric_function(I[1,2])
                    0
                    sage: Phi = N.Phi()
                    sage: Phi.to_symmetric_function_on_basis([3,1,2])==Phi.to_symmetric_function(Phi[3,1,2])
                    True
                    sage: Phi.to_symmetric_function_on_basis([])
                    h[]
                """
                S = self.realization_of().complete()
                return S.to_symmetric_function(S(self[I]))

            @lazy_attribute
            def to_symmetric_function(self):
                r"""
                Morphism to the algebra of symmetric functions.

                This is constructed by extending the computation on the basis
                or by coercion to the complete basis.

                OUTPUT:

                - The module morphism from the basis ``self`` to the symmetric
                  functions which corresponds to taking a commutative image.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: x = R.an_element(); x
                    2*R[] + 2*R[1] + 3*R[1, 1]
                    sage: R.to_symmetric_function(x)
                    2*s[] + 2*s[1] + 3*s[1, 1]
                    sage: nM = N.Monomial()
                    sage: nM.to_symmetric_function(nM[3,1])
                    h[1, 1, 1, 1] - 7/2*h[2, 1, 1] + h[2, 2] + 7/2*h[3, 1] - 2*h[4]
                """
                on_basis = self.to_symmetric_function_on_basis
                codom = on_basis([]).parent()
                return self.module_morphism(on_basis, codomain=codom)

            def immaculate_function(self, xs):
                r"""
                Return the immaculate function corresponding to the
                integer vector ``xs``, written in the basis ``self``.

                If `\alpha` is any integer vector -- i.e., an element of
                `\ZZ^m` for some `m \in \NN` --, the *immaculate function
                corresponding to* `\alpha` is a non-commutative symmetric
                function denoted by `\mathfrak{S}_{\alpha}`. One way to
                define this function is by setting

                .. MATH::

                    \mathfrak{S}_{\alpha}
                    = \sum_{\sigma \in S_m} (-1)^{\sigma}
                    S_{\alpha_1 + \sigma(1) - 1}
                    S_{\alpha_2 + \sigma(2) - 2}
                    \cdots
                    S_{\alpha_m + \sigma(m) - m},

                where `\alpha` is written in the form
                `(\alpha_1, \alpha_2, \ldots, \alpha_m)`, and where `S`
                stands for the complete basis
                (:class:`~NonCommutativeSymmetricFunctions.Complete`).

                The immaculate function `\mathfrak{S}_{\alpha}` first
                appeared in [BBSSZ2012]_ (where it was defined
                differently, but the definition we gave above appeared as
                Theorem 3.27).

                The immaculate functions `\mathfrak{S}_{\alpha}` for
                `\alpha` running over all compositions (i.e., finite
                sequences of positive integers) form a basis of `NCSF`.
                This is the *immaculate basis*
                (:class:`~NonCommutativeSymmetricFunctions.Immaculate`).

                INPUT:

                - ``xs`` -- list (or tuple or any iterable -- possibly a
                  composition) of integers

                OUTPUT:

                The immaculate function `\mathfrak{S}_{xs}`
                written in the basis ``self``.

                EXAMPLES:

                Let us first check that, for ``xs`` a composition, we get
                the same as the result of
                ``self.realization_of().I()[xs]``::

                    sage: def test_comp(xs):
                    ....:     NSym = NonCommutativeSymmetricFunctions(QQ)
                    ....:     I = NSym.I()
                    ....:     return I[xs] == I.immaculate_function(xs)
                    sage: def test_allcomp(n):
                    ....:     return all( test_comp(c) for c in Compositions(n) )
                    sage: test_allcomp(1)
                    True
                    sage: test_allcomp(2)
                    True
                    sage: test_allcomp(3)
                    True

                Now some examples of non-composition immaculate
                functions::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: I = NSym.I()
                    sage: I.immaculate_function([0, 1])
                    0
                    sage: I.immaculate_function([0, 2])
                    -I[1, 1]
                    sage: I.immaculate_function([-1, 1])
                    -I[]
                    sage: I.immaculate_function([2, -1])
                    0
                    sage: I.immaculate_function([2, 0])
                    I[2]
                    sage: I.immaculate_function([2, 0, 1])
                    0
                    sage: I.immaculate_function([1, 0, 2])
                    -I[1, 1, 1]
                    sage: I.immaculate_function([2, 0, 2])
                    -I[2, 1, 1]
                    sage: I.immaculate_function([0, 2, 0, 2])
                    I[1, 1, 1, 1] + I[1, 2, 1]
                    sage: I.immaculate_function([2, 0, 2, 0, 2])
                    I[2, 1, 1, 1, 1] + I[2, 1, 2, 1]

                TESTS:

                Basis-independence::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: L = NSym.L()
                    sage: S = NSym.S()
                    sage: L(S.immaculate_function([0, 2, 0, 2])) == L.immaculate_function([0, 2, 0, 2])
                    True
                """
                S = self.realization_of().S()
                res = S.zero()
                m = len(xs)
                ys = [xs_i - i - 1 for i, xs_i in enumerate(xs)]
                for s in Permutations(m):
                    psco = [ys[i] + s_i for i, s_i in enumerate(s)]
                    if not all(j >= 0 for j in psco):
                        continue
                    psco2 = [j for j in psco if j != 0]
                    pr = s.sign() * S[psco2]
                    res += pr
                return self(res)

        class ElementMethods:

            def verschiebung(self, n):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the `n`-th Verschiebung operator.

                The `n`-th Verschiebung operator `\mathbf{V}_n` is defined
                to be the map from the `\mathbf{k}`-algebra of noncommutative
                symmetric functions to itself that sends the complete function
                `S^I` indexed by a composition `I = (i_1, i_2, \ldots , i_k)`
                to `S^{(i_1/n, i_2/n, \ldots , i_k/n)}` if all of the numbers
                `i_1, i_2, \ldots, i_k` are divisible by `n`, and to `0`
                otherwise. This operator `\mathbf{V}_n` is a Hopf algebra
                endomorphism. For every positive integer `r` with `n \mid r`,
                it satisfies

                .. MATH::

                    \mathbf{V}_n(S_r) = S_{r/n},
                    \quad \mathbf{V}_n(\Lambda_r) = (-1)^{r - r/n} \Lambda_{r/n},
                    \quad \mathbf{V}_n(\Psi_r) = n \Psi_{r/n},
                    \quad \mathbf{V}_n(\Phi_r) = n \Phi_{r/n}

                (where `S_r` denotes the `r`-th complete non-commutative
                symmetric function, `\Lambda_r` denotes the `r`-th elementary
                non-commutative symmetric function, `\Psi_r` denotes the `r`-th
                power-sum non-commutative symmetric function of the first kind,
                and `\Phi_r` denotes the `r`-th power-sum non-commutative
                symmetric function of the second kind). For every positive
                integer `r` with `n \nmid r`, it satisfes

                .. MATH::

                    \mathbf{V}_n(S_r) = \mathbf{V}_n(\Lambda_r)
                    = \mathbf{V}_n(\Psi_r) = \mathbf{V}_n(\Phi_r) = 0.

                The `n`-th Verschiebung operator is also called the `n`-th
                Verschiebung endomorphism.

                It is a lift of the `n`-th Verschiebung operator on the ring
                of symmetric functions
                (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung`)
                to the ring of noncommutative symmetric functions.

                The action of the `n`-th Verschiebung operator can also be
                described on the ribbon Schur functions. Namely, every
                composition `I` of size `n \ell` satisfies

                .. MATH::

                    \mathbf{V}_n ( R_I )
                    = (-1)^{\ell(I) - \ell(J)} \cdot R_{J / n},

                where `J` denotes the meet of the compositions `I` and
                `(\underbrace{n, n, \ldots, n}_{|I|/n \mbox{ times}})`,
                where `\ell(I)` is the length of `I`, and
                where `J / n` denotes the composition obtained by dividing
                every entry of `J` by `n`.
                For a composition `I` of size not divisible by `n`, we have
                `\mathbf{V}_n( R_I ) = 0`.

                .. SEEALSO::

                    :meth:`frobenius method of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.frobenius>`,
                    :meth:`verschiebung method of Sym
                    <sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung>`

                INPUT:

                - ``n`` -- a positive integer

                OUTPUT:

                The result of applying the `n`-th Verschiebung operator (on the
                ring of noncommutative symmetric functions) to ``self``.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: S = NSym.S()
                    sage: S[3,2].verschiebung(2)
                    0
                    sage: S[6,4].verschiebung(2)
                    S[3, 2]
                    sage: (S[9,1] - S[8,2] + 2*S[6,4] - 3*S[3] + 4*S[[]]).verschiebung(2)
                    4*S[] + 2*S[3, 2] - S[4, 1]
                    sage: (S[3,3] - 2*S[2]).verschiebung(3)
                    S[1, 1]
                    sage: S([4,2]).verschiebung(1)
                    S[4, 2]
                    sage: R = NSym.R()
                    sage: R([4,2]).verschiebung(2)
                    R[2, 1]

                Being Hopf algebra endomorphisms, the Verschiebung operators
                commute with the antipode::

                    sage: all( R(I).verschiebung(2).antipode()
                    ....:      == R(I).antipode().verschiebung(2)
                    ....:      for I in Compositions(4) )
                    True

                They lift the Verschiebung operators of the ring of symmetric
                functions::

                    sage: all( S(I).verschiebung(2).to_symmetric_function()
                    ....:      == S(I).to_symmetric_function().verschiebung(2)
                    ....:      for I in Compositions(4) )
                    True

                The Frobenius operators on `QSym` are adjoint to the
                Verschiebung operators on `NSym` with respect to the duality
                pairing::

                    sage: QSym = QuasiSymmetricFunctions(ZZ)
                    sage: M = QSym.M()
                    sage: all( all( M(I).frobenius(3).duality_pairing(S(J))
                    ....:           == M(I).duality_pairing(S(J).verschiebung(3))
                    ....:           for I in Compositions(2) )
                    ....:      for J in Compositions(3) )
                    True
                """
                # Convert to the homogeneous basis, there apply Verschiebung
                # componentwise, then convert back.
                parent = self.parent()
                S = parent.realization_of().S()
                C = parent._indices
                dct = {C([i // n for i in I]): coeff
                       for (I, coeff) in S(self) if all(i % n == 0 for i in I)}
                return parent(S._from_dict(dct))

            def bernstein_creation_operator(self, n):
                r"""
                Return the image of ``self`` under the `n`-th Bernstein
                creation operator.

                Let `n` be an integer. The `n`-th Bernstein creation
                operator `\mathbb{B}_n` is defined as the endomorphism of
                the space `NSym` of noncommutative symmetric functions
                which sends every `f` to

                .. MATH::

                    \sum_{i \geq 0} (-1)^i H_{n+i} F_{1^i}^\perp,

                where usual notations are in place (the letter `H` stands
                for the complete basis of `NSym`, the letter `F` stands
                for the fundamental basis of the algebra `QSym` of
                quasisymmetric functions, and `F_{1^i}^\perp` means
                skewing (:meth:`~sage.combinat.ncsf_qsym.generic_basis_code.BasesOfQSymOrNCSF.ElementMethods.skew_by`)
                by `F_{1^i}`). Notice that `F_{1^i}` is nothing other than the
                elementary symmetric function `e_i`.

                This has been introduced in [BBSSZ2012]_, section 3.1, in
                analogy to the Bernstein creation operators on the
                symmetric functions
                (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.bernstein_creation_operator`),
                and studied further in [BBSSZ2012]_, mainly in the context
                of immaculate functions
                (:class:`~NonCommutativeSymmetricFunctions.Immaculate`).
                In fact, if `(\alpha_1, \alpha_2, \ldots, \alpha_m)` is
                an `m`-tuple of integers, then

                .. MATH::

                    \mathbb{B}_n I_{(\alpha_1, \alpha_2, \ldots, \alpha_m)}
                    = I_{(n, \alpha_1, \alpha_2, \ldots, \alpha_m)},

                where `I_{(\alpha_1, \alpha_2, \ldots, \alpha_m)}` is the
                immaculate function associated to the `m`-tuple
                `(\alpha_1, \alpha_2, \ldots, \alpha_m)` (see
                :meth:`~NonCommutativeSymmetricFunctions.Bases.ParentMethods.immaculate_function`).

                EXAMPLES:

                We get the immaculate functions by repeated application of
                Bernstein creation operators::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: I = NSym.I()
                    sage: S = NSym.S()
                    sage: def immaculate_by_bernstein(xs):
                    ....:     # immaculate function corresponding to integer
                    ....:     # tuple ``xs``, computed by iterated application
                    ....:     # of Bernstein creation operators.
                    ....:     res = S.one()
                    ....:     for i in reversed(xs):
                    ....:         res = res.bernstein_creation_operator(i)
                    ....:     return res
                    sage: import itertools
                    sage: all( immaculate_by_bernstein(p) == I.immaculate_function(p)
                    ....:      for p in itertools.product(range(-1, 3), repeat=3))
                    True

                Some examples::

                    sage: S[3,2].bernstein_creation_operator(-2)
                    S[2, 1]
                    sage: S[3,2].bernstein_creation_operator(-1)
                    S[1, 2, 1] - S[2, 2] - S[3, 1]
                    sage: S[3,2].bernstein_creation_operator(0)
                    -S[1, 2, 2] - S[1, 3, 1] + S[2, 2, 1] + S[3, 2]
                    sage: S[3,2].bernstein_creation_operator(1)
                    S[1, 3, 2] - S[2, 2, 2] - S[2, 3, 1] + S[3, 2, 1]
                    sage: S[3,2].bernstein_creation_operator(2)
                    S[2, 3, 2] - S[3, 2, 2] - S[3, 3, 1] + S[4, 2, 1]
                """
                # We use the definition of this operator.
                parent = self.parent()
                res = parent.zero()
                if not self:
                    return res
                max_degree = max(sum(m) for m, c in self)
                # ``max_degree`` is now the maximum degree in which ``self``
                # has a nonzero coefficient.
                NSym = parent.realization_of()
                S = NSym.S()
                F = NSym.dual().F()
                for i in range(max_degree + 1):
                    if n + i > 0:
                        res += (-1) ** i * S[n + i] * self.skew_by(F[[1] * i])
                    elif n + i == 0:
                        res += (-1) ** i * self.skew_by(F[[1] * i])
                return res

            def star_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the star involution.

                The star involution is defined as the algebra antihomomorphism
                `NCSF \to NCSF` which, for every positive integer `n`, sends
                the `n`-th complete non-commutative symmetric function `S_n` to
                `S_n`. Denoting by `f^{\ast}` the image of an element
                `f \in NCSF` under this star involution, it can be shown that
                every composition `I` satisfies

                .. MATH::

                    (S^I)^{\ast} = S^{I^r}, \quad
                    (\Lambda^I)^{\ast} = \Lambda^{I^r}, \quad
                    R_I^{\ast} = R_{I^r}, \quad
                    (\Phi^I)^{\ast} = \Phi^{I^r},

                where `I^r` denotes the reversed composition of `I`, and
                standard notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary basis,
                `R` for the ribbon basis, and `\Phi` for that of the power-sums
                of the second kind). The star involution is an involution and a
                coalgebra automorphism of `NCSF`. It is an automorphism of the
                graded vector space `NCSF`. Under the canonical isomorphism
                between the `n`-th graded component of `NCSF` and the descent
                algebra of the symmetric group `S_n` (see
                :meth:`to_descent_algebra`), the star involution (restricted to
                the `n`-th graded component) corresponds to the automorphism
                of the descent algebra given by
                `x \mapsto \omega_n x \omega_n`, where `\omega_n` is the
                permutation `(n, n-1, \ldots, 1) \in S_n` (written here in
                one-line notation). If `\pi` denotes the projection from `NCSF`
                to the ring of symmetric functions
                (:meth:`to_symmetric_function`), then `\pi(f^{\ast}) = \pi(f)`
                for every `f \in NCSF`.

                The star involution on `NCSF` is adjoint to the star involution
                on `QSym` by the standard adjunction between `NCSF` and `QSym`.

                The star involution has been denoted by `\rho` in [LMvW13]_,
                section 3.6.
                See [NCSF2]_, section 2.3 for the properties of this map.

                .. SEEALSO::

                    :meth:`star involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`,
                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: S = NSym.S()
                    sage: S[3,2].star_involution()
                    S[2, 3]
                    sage: S[6,3].star_involution()
                    S[3, 6]
                    sage: (S[9,1] - S[8,2] + 2*S[6,4] - 3*S[3] + 4*S[[]]).star_involution()
                    4*S[] + S[1, 9] - S[2, 8] - 3*S[3] + 2*S[4, 6]
                    sage: (S[3,3] - 2*S[2]).star_involution()
                    -2*S[2] + S[3, 3]
                    sage: S([4,2]).star_involution()
                    S[2, 4]
                    sage: R = NSym.R()
                    sage: R([4,2]).star_involution()
                    R[2, 4]
                    sage: R.zero().star_involution()
                    0
                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: Phi = NSym.Phi()
                    sage: Phi([2,1]).star_involution()
                    Phi[1, 2]

                The Psi basis doesn't behave as nicely::

                    sage: Psi = NSym.Psi()
                    sage: Psi([2,1]).star_involution()
                    Psi[1, 2]
                    sage: Psi([3,1]).star_involution()
                    1/2*Psi[1, 1, 2] - 1/2*Psi[1, 2, 1] + Psi[1, 3]

                The star involution commutes with the antipode::

                    sage: all( R(I).star_involution().antipode()
                    ....:      == R(I).antipode().star_involution()
                    ....:      for I in Compositions(4) )
                    True

                Checking the relation with the descent algebra described
                above::

                    sage: def descent_test(n):
                    ....:     DA = DescentAlgebra(QQ, n)
                    ....:     NSym = NonCommutativeSymmetricFunctions(QQ)
                    ....:     S = NSym.S()
                    ....:     DAD = DA.D()
                    ....:     w_n = DAD(set(range(1, n)))
                    ....:     for I in Compositions(n):
                    ....:         if not (S[I].star_involution()
                    ....:                 == w_n * S[I].to_descent_algebra(n) * w_n):
                    ....:             return False
                    ....:         return True
                    sage: all( descent_test(i) for i in range(4) )
                    True
                    sage: all( descent_test(i) for i in range(6) ) # long time
                    True

                Testing the `\pi(f^{\ast}) = \pi(f)` relation noticed above::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = NSym.R()
                    sage: all( R(I).star_involution().to_symmetric_function()
                    ....:      == R(I).to_symmetric_function()
                    ....:      for I in Compositions(4) )
                    True

                The star involution on `QSym` is adjoint to the star involution
                on `NSym` with respect to the duality pairing::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.M()
                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = NSym.S()
                    sage: all( all( M(I).star_involution().duality_pairing(S(J))
                    ....:           == M(I).duality_pairing(S(J).star_involution())
                    ....:           for I in Compositions(2) )
                    ....:      for J in Compositions(3) )
                    True
                """
                # Convert to the homogeneous basis, there apply the star
                # involution componentwise, then convert back.
                parent = self.parent()
                S = parent.realization_of().S()
                dct = {I.reversed(): coeff for (I, coeff) in S(self)}
                return parent(S._from_dict(dct))

            def omega_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the omega involution.

                The omega involution is defined as the algebra antihomomorphism
                `NCSF \to NCSF` which, for every positive integer `n`, sends
                the `n`-th complete non-commutative symmetric function `S_n` to
                the `n`-th elementary non-commutative symmetric function
                `\Lambda_n`. This omega involution is denoted by `\omega`. It
                can be shown that every composition `I` satisfies

                .. MATH::

                    \omega(S^I) = \Lambda^{I^r}, \quad
                    \omega(\Lambda^I) = S^{I^r}, \quad
                    \omega(R_I) = R_{I^t}, \quad
                    \omega(\Phi^I) = (-1)^{|I|-\ell(I)} \Phi^{I^r},
                    \omega(\Psi^I) = (-1)^{|I|-\ell(I)} \Psi^{I^r},

                where `I^r` denotes the reversed composition of `I`, and
                `I^t` denotes the conjugate composition of `I`, and `\ell(I)`
                denotes the length of the
                composition `I`, and standard notations for classical bases
                of `NCSF` are being used (`S` for the complete basis,
                `\Lambda` for the elementary basis, `R` for the ribbon
                basis, `\Phi` for that of the power-sums of the second
                kind, and `\Psi` for that of the power-sums of the first
                kind). More generally, if `f` is a homogeneous element of
                `NCSF` of degree `n`, then

                .. MATH::

                    \omega(f) = (-1)^n S(f),

                where `S` denotes the antipode of `NCSF`.

                The omega involution `\omega` is an involution and a
                coalgebra automorphism of `NCSF`. It is an automorphism of the
                graded vector space `NCSF`. If `\pi` denotes the projection
                from `NCSF` to the ring of symmetric functions
                (:meth:`to_symmetric_function`), then
                `\pi(\omega(f)) = \omega(\pi(f))` for every `f \in NCSF`,
                where the `\omega` on the right hand side denotes the omega
                automorphism of `Sym`.

                The omega involution on `NCSF` is adjoint to the omega
                involution on `QSym` by the standard adjunction between `NCSF`
                and `QSym`.

                The omega involution has been denoted by `\omega` in [LMvW13]_,
                section 3.6.
                See [NCSF1]_, section 3.1 for the properties of this map.

                .. SEEALSO::

                    :meth:`omega involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.omega_involution>`,
                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: S = NSym.S()
                    sage: L = NSym.L()
                    sage: L(S[3,2].omega_involution())
                    L[2, 3]
                    sage: L(S[6,3].omega_involution())
                    L[3, 6]
                    sage: L(S[1,3].omega_involution())
                    L[3, 1]
                    sage: L((S[9,1] - S[8,2] + 2*S[6,4] - 3*S[3] + 4*S[[]]).omega_involution()) # long time
                    4*L[] + L[1, 9] - L[2, 8] - 3*L[3] + 2*L[4, 6]
                    sage: L((S[3,3] - 2*S[2]).omega_involution())
                    -2*L[2] + L[3, 3]
                    sage: L(S([4,2]).omega_involution())
                    L[2, 4]
                    sage: R = NSym.R()
                    sage: R([4,2]).omega_involution()
                    R[1, 2, 1, 1, 1]
                    sage: R.zero().omega_involution()
                    0
                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: Phi = NSym.Phi()
                    sage: Phi([2,1]).omega_involution()
                    -Phi[1, 2]
                    sage: Psi = NSym.Psi()
                    sage: Psi([2,1]).omega_involution()
                    -Psi[1, 2]
                    sage: Psi([3,1]).omega_involution()
                    Psi[1, 3]

                Testing the `\pi(\omega(f)) = \omega(\pi(f))` relation noticed
                above::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = NSym.R()
                    sage: all( R(I).omega_involution().to_symmetric_function()
                    ....:      == R(I).to_symmetric_function().omega_involution()
                    ....:      for I in Compositions(4) )
                    True

                The omega involution on `QSym` is adjoint to the omega
                involution on `NSym` with respect to the duality pairing::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.M()
                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = NSym.S()
                    sage: all( all( M(I).omega_involution().duality_pairing(S(J))
                    ....:           == M(I).duality_pairing(S(J).omega_involution())
                    ....:           for I in Compositions(2) )
                    ....:      for J in Compositions(3) )
                    True
                """
                # Use the `\omega(f) = (-1)^n S(f)` formula.
                return self.antipode().degree_negation()

            def psi_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the involution `\psi`.

                The involution `\psi` is defined as the linear map
                `NCSF \to NCSF` which, for every composition `I`, sends the
                complete noncommutative symmetric function `S^I` to the
                elementary noncommutative symmetric function `\Lambda^I`.
                It can be shown that every composition `I` satisfies

                .. MATH::

                    \psi(R_I) = R_{I^c}, \quad \psi(S^I) = \Lambda^I, \quad
                    \psi(\Lambda^I) = S^I, \quad
                    \psi(\Phi^I) = (-1)^{|I| - \ell(I)} \Phi^I

                where `I^c` denotes the complement of the composition `I`, and
                `\ell(I)` denotes the length of `I`, and where standard
                notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary basis,
                `\Phi` for the basis of the power sums of the second kind,
                and `R` for the ribbon basis). The map `\psi` is an involution
                and a graded Hopf algebra automorphism of `NCSF`. If `\pi`
                denotes the projection from `NCSF` to the ring of symmetric
                functions (:meth:`to_symmetric_function`), then
                `\pi(\psi(f)) = \omega(\pi(f))` for every `f \in NCSF`, where
                the `\omega` on the right hand side denotes the omega
                automorphism of `Sym`.

                The involution `\psi` of `NCSF` is adjoint to the involution
                `\psi` of `QSym` by the standard adjunction between `NCSF` and
                `QSym`.

                The involution `\psi` has been denoted by `\psi` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`psi involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: R = NSym.R()
                    sage: R[3,2].psi_involution()
                    R[1, 1, 2, 1]
                    sage: R[6,3].psi_involution()
                    R[1, 1, 1, 1, 1, 2, 1, 1]
                    sage: (R[9,1] - R[8,2] + 2*R[2,4] - 3*R[3] + 4*R[[]]).psi_involution()
                    4*R[] - 3*R[1, 1, 1] + R[1, 1, 1, 1, 1, 1, 1, 1, 2] - R[1, 1, 1, 1, 1, 1, 1, 2, 1] + 2*R[1, 2, 1, 1, 1]
                    sage: (R[3,3] - 2*R[2]).psi_involution()
                    -2*R[1, 1] + R[1, 1, 2, 1, 1]
                    sage: R([2,1,1]).psi_involution()
                    R[1, 3]
                    sage: S = NSym.S()
                    sage: S([2,1]).psi_involution()
                    S[1, 1, 1] - S[2, 1]
                    sage: S.zero().psi_involution()
                    0
                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: Phi = NSym.Phi()
                    sage: Phi([2,1]).psi_involution()
                    -Phi[2, 1]
                    sage: Phi([3,1]).psi_involution()
                    Phi[3, 1]

                The Psi basis doesn't behave as nicely::

                    sage: Psi = NSym.Psi()
                    sage: Psi([2,1]).psi_involution()
                    -Psi[2, 1]
                    sage: Psi([3,1]).psi_involution()
                    1/2*Psi[1, 2, 1] - 1/2*Psi[2, 1, 1] + Psi[3, 1]

                The involution `\psi` commutes with the antipode::

                    sage: all( R(I).psi_involution().antipode()
                    ....:      == R(I).antipode().psi_involution()
                    ....:      for I in Compositions(4) )
                    True

                Testing the `\pi(\psi(f)) = \omega(\pi(f))` relation noticed
                above::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = NSym.R()
                    sage: all( R(I).psi_involution().to_symmetric_function()
                    ....:      == R(I).to_symmetric_function().omega()
                    ....:      for I in Compositions(4) )
                    True

                The involution `\psi` of `QSym` is adjoint to the involution
                `\psi` of `NSym` with respect to the duality pairing::

                    sage: QSym = QuasiSymmetricFunctions(QQ)
                    sage: M = QSym.M()
                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = NSym.S()
                    sage: all( all( M(I).psi_involution().duality_pairing(S(J))
                    ....:           == M(I).duality_pairing(S(J).psi_involution())
                    ....:           for I in Compositions(2) )
                    ....:      for J in Compositions(3) )
                    True
                """
                # Convert to the ribbon basis, there apply the psi
                # involution componentwise, then convert back.
                parent = self.parent()
                R = parent.realization_of().R()
                dct = {I.complement(): coeff for (I, coeff) in R(self)}
                return parent(R._from_dict(dct))

            def left_padded_kronecker_product(self, x):
                r"""
                Return the left-padded Kronecker product of ``self`` and
                ``x`` in the basis of ``self``.

                The left-padded Kronecker product is a bilinear map
                mapping two non-commutative symmetric functions to
                another, not necessarily preserving degree.
                It can be defined as follows: Let `*` denote the internal
                product (:meth:`~sage.combinat.ncsf_qsym.generic_basis_code.GradedModulesWithInternalProduct.ElementMethods.internal_product`)
                on the space of non-commutative symmetric functions. For any
                composition `I`, let `S^I` denote the complete homogeneous
                symmetric function indexed by `I`. For any compositions
                `\alpha`, `\beta`, `\gamma`, let
                `g^{\gamma}_{\alpha, \beta}` denote the coefficient of
                `S^{\gamma}` in the internal product
                `S^{\alpha} * S^{\beta}`.
                For every composition `I = (i_1, i_2, \ldots, i_k)`
                and every integer `n > \left\lvert I \right\rvert`, define the
                *`n`-completion of `I`* to be the composition
                `(n - \left\lvert I \right\rvert, i_1, i_2, \ldots, i_k)`;
                this `n`-completion is denoted by `I[n]`.
                Then, for any compositions `\alpha` and `\beta` and every
                integer `n > \left\lvert \alpha \right\rvert
                + \left\lvert\beta\right\rvert`, we can write the
                internal product `S^{\alpha[n]} * S^{\beta[n]}` in the form

                .. MATH::

                    S^{\alpha[n]} * S^{\beta[n]} =
                    \sum_{\gamma} g^{\gamma[n]}_{\alpha[n], \beta[n]}
                    S^{\gamma[n]}

                with `\gamma` ranging over all compositions. The
                coefficients `g^{\gamma[n]}_{\alpha[n], \beta[n]}`
                are independent on `n`. These coefficients
                `g^{\gamma[n]}_{\alpha[n], \beta[n]}` are denoted by
                `\widetilde{g}^{\gamma}_{\alpha, \beta}`, and the
                non-commutative symmetric function

                .. MATH::

                    \sum_{\gamma} \widetilde{g}^{\gamma}_{\alpha, \beta}
                    S^{\gamma}

                is said to be the *left-padded Kronecker product* of
                `S^{\alpha}` and `S^{\beta}`. By bilinearity, this
                extends to a definition of a left-padded Kronecker product
                of any two non-commutative symmetric functions.

                The left-padded Kronecker product on the non-commutative
                symmetric functions lifts the left-padded Kronecker
                product on the symmetric functions. More precisely: Let
                `\pi` denote the canonical projection
                (:meth:`to_symmetric_function`) from the non-commutative
                symmetric functions to the symmetric functions. Then, any
                two non-commutative symmetric functions `f` and `g`
                satisfy

                .. MATH::

                    \pi(f \overline{*} g) = \pi(f) \overline{*} \pi(g),

                where the `\overline{*}` on the left-hand side denotes the
                left-padded Kronecker product on the non-commutative
                symmetric functions, and the `\overline{*}` on the
                right-hand side denotes the left-padded Kronecker product
                on the symmetric functions.

                INPUT:

                - ``x`` -- element of the ring of non-commutative
                  symmetric functions over the same base ring as ``self``

                OUTPUT:

                - the left-padded Kronecker product of ``self`` with ``x``
                  (an element of the ring of non-commutative symmetric
                  functions in the same basis as ``self``)

                AUTHORS:

                - Darij Grinberg (15 Mar 2014)

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = NSym.S()
                    sage: S[2,1].left_padded_kronecker_product(S[3])
                    S[1, 1, 1, 1] + S[1, 2, 1] + S[2, 1] + S[2, 1, 1, 1] + S[2, 2, 1] + S[3, 2, 1]
                    sage: S[2,1].left_padded_kronecker_product(S[1])
                    S[1, 1, 1] + S[1, 2, 1] + S[2, 1]
                    sage: S[1].left_padded_kronecker_product(S[2,1])
                    S[1, 1, 1] + S[2, 1] + S[2, 1, 1]
                    sage: S[1,1].left_padded_kronecker_product(S[2])
                    S[1, 1] + 2*S[1, 1, 1] + S[2, 1, 1]
                    sage: S[1].left_padded_kronecker_product(S[1,2,1])
                    S[1, 1, 1, 1] + S[1, 2, 1] + S[1, 2, 1, 1] + S[2, 1, 1]
                    sage: S[2].left_padded_kronecker_product(S[3])
                    S[1, 2] + S[2, 1, 1] + S[3, 2]

                Taking the left-padded Kronecker product with
                `1 = S^{\empty}` is the identity map on the ring of
                non-commutative symmetric functions::

                    sage: all( S[Composition([])].left_padded_kronecker_product(S[lam])
                    ....:      == S[lam].left_padded_kronecker_product(S[Composition([])])
                    ....:      == S[lam] for i in range(4)
                    ....:      for lam in Compositions(i) )
                    True

                Here is a rule for the left-padded Kronecker product of
                `S_1` (this is the same as `S^{(1)}`) with any complete
                homogeneous function: Let `I` be a composition.
                Then, the left-padded Kronecker product of `S_1` and
                `S^I` is `\sum_K a_K S^K`, where the sum runs
                over all compositions `K`, and the coefficient `a_K` is
                defined as the number of ways to obtain `K` from `I` by
                one of the following two operations:

                - Insert a `1` at the end of `I`.
                - Subtract `1` from one of the entries of `I` (and remove
                  the entry if it thus becomes `0`), and insert a `1` at
                  the end of `I`.

                We check this for compositions of size `\leq 4`::

                    sage: def mults1(I):
                    ....:     # Left left-padded Kronecker multiplication by S[1].
                    ....:     res = S[I[:] + [1]]
                    ....:     for k in range(len(I)):
                    ....:         I2 = I[:]
                    ....:         if I2[k] == 1:
                    ....:             I2 = I2[:k] + I2[k+1:]
                    ....:         else:
                    ....:             I2[k] -= 1
                    ....:         res += S[I2 + [1]]
                    ....:     return res
                    sage: all( mults1(I) == S[1].left_padded_kronecker_product(S[I])
                    ....:      for i in range(5) for I in Compositions(i) )
                    True

                A similar rule can be made for the left-padded Kronecker
                product of any complete homogeneous function with `S_1`:
                Let `I` be a composition. Then, the left-padded Kronecker
                product of `S^I` and `S_1` is `\sum_K b_K S^K`, where the
                sum runs over all compositions `K`, and the coefficient
                `b_K` is defined as the number of ways to obtain `K` from
                `I` by one of the following two operations:

                - Insert a `1` at the front of `I`.
                - Subtract `1` from one of the entries of `I` (and remove
                  the entry if it thus becomes `0`), and insert a `1`
                  right after this entry.

                We check this for compositions of size `\leq 4`::

                    sage: def mults2(I):
                    ....:     # Left left-padded Kronecker multiplication by S[1].
                    ....:     res = S[[1] + I[:]]
                    ....:     for k in range(len(I)):
                    ....:         I2 = I[:]
                    ....:         i2k = I2[k]
                    ....:         if i2k != 1:
                    ....:             I2 = I2[:k] + [i2k-1, 1] + I2[k+1:]
                    ....:         res += S[I2]
                    ....:     return res
                    sage: all( mults2(I) == S[I].left_padded_kronecker_product(S[1])
                    ....:      for i in range(5) for I in Compositions(i) )
                    True

                Checking the
                `\pi(f \overline{*} g) = \pi(f) \overline{*} \pi(g)`
                equality::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: R = NSym.R()
                    sage: def testpi(n):
                    ....:     for I in Compositions(n):
                    ....:         for J in Compositions(n):
                    ....:             a = R[I].to_symmetric_function()
                    ....:             b = R[J].to_symmetric_function()
                    ....:             x = a.left_padded_kronecker_product(b)
                    ....:             y = R[I].left_padded_kronecker_product(R[J])
                    ....:             y = y.to_symmetric_function()
                    ....:             if x != y:
                    ....:                 return False
                    ....:     return True
                    sage: testpi(3)
                    True

                TESTS::

                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: (2*S([])).left_padded_kronecker_product(3*S([]))
                    6*S[]

                Different bases and base rings::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: S = NSym.S()
                    sage: L = NSym.L()
                    sage: L(S[2].left_padded_kronecker_product(L[2]))
                    L[1, 1, 1] + L[2] + L[2, 1, 1] - L[2, 2]
                    sage: S(L[2].left_padded_kronecker_product(S[1,1]))
                    S[1, 1] + 2*S[1, 1, 1] + S[1, 1, 1, 1] - S[1, 1, 2]

                    sage: NSym = NonCommutativeSymmetricFunctions(CyclotomicField(12))
                    sage: S = NSym.S()
                    sage: L = NSym.L()
                    sage: v = L[2].left_padded_kronecker_product(L[2]); v
                    L[1, 1] + L[1, 1, 1] + (-1)*L[2] + L[2, 2]
                    sage: parent(v)
                    Non-Commutative Symmetric Functions over the Cyclotomic Field of order 12 and degree 4 in the Elementary basis

                    sage: NSym = NonCommutativeSymmetricFunctions(Zmod(14))
                    sage: S = NSym.S()
                    sage: L = NSym.L()
                    sage: v = L[2].left_padded_kronecker_product(L[2]); v
                    L[1, 1] + L[1, 1, 1] + 13*L[2] + L[2, 2]
                    sage: parent(v)
                    Non-Commutative Symmetric Functions over the Ring of integers modulo 14 in the Elementary basis

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: R = NSym.R()
                    sage: v = R[1].left_padded_kronecker_product(R[1]); parent(v)
                    Non-Commutative Symmetric Functions over the Integer Ring in the Ribbon basis
                """
                _Compositions = Compositions()
                parent = self.parent()
                comp_parent = parent.realization_of().S()
                comp_self = comp_parent(self)
                comp_x = comp_parent(x)
                # Now, comp_self and comp_x are the same as self and x, but in the
                # S (=complete homogeneous) basis, which we call comp_parent.
                result = comp_parent.zero()
                for lam, a in comp_self:
                    # lam is a composition, a is an element of the base ring.
                    lam_list = lam._list
                    if not lam._list:
                        # Special handling for the empty composition. The left-padded
                        # Kronecker product of 1 with any non-commutative symmetric
                        # function f is f.
                        result += a * comp_x
                        continue
                    sum_lam = sum(lam_list)
                    for mu, b in comp_x:
                        # mu is a composition, b is an element of the base ring.
                        mu_list = mu._list
                        if not mu_list:
                            # Special handling for the empty composition.
                            result += a * b * comp_parent(lam)
                            continue
                        # Now, both lam and mu are nonempty.
                        sum_mu = sum(mu_list)
                        stab = 1 + sum_lam + sum_mu
                        S_lam_stabilized = comp_parent(_Compositions([stab - sum_lam] + lam_list))
                        S_mu_stabilized = comp_parent(_Compositions([stab - sum_mu] + mu_list))
                        lam_star_mu = S_lam_stabilized.internal_product(S_mu_stabilized)
                        # lam_star_mu is now a non-commutative symmetric function
                        # in the S-basis.
                        for nu, c in lam_star_mu:
                            # nu is a composition of the integer stab, c is an element
                            # of the base ring.
                            nu_unstabilized = _Compositions(nu[1:])
                            result += a * b * c * comp_parent(nu_unstabilized)
                return parent(result)

            def to_descent_algebra(self, n=None):
                r"""
                Return the image of the ``n``-th degree homogeneous component
                of ``self`` in the descent algebra of `S_n` over the same
                base ring as ``self``.

                This is based upon the canonical isomorphism from the
                `n`-th degree homogeneous component of the algebra of
                noncommutative symmetric functions to the descent algebra
                of `S_n`. This isomorphism maps the inner product of
                noncommutative symmetric functions either to the product
                in the descent algebra of `S_n` or to its opposite
                (depending on how the latter is defined).

                If ``n`` is not specified, it will be taken to be the highest
                homogeneous component of ``self``.

                OUTPUT:

                - The image of the ``n``-th homogeneous component of ``self``
                  under the isomorphism into the descent algebra of `S_n`
                  over the same base ring as ``self``.

                EXAMPLES::

                    sage: S = NonCommutativeSymmetricFunctions(ZZ).S()
                    sage: S[2,1].to_descent_algebra(3)
                    B[2, 1]
                    sage: (S[1,2,1] - 3 * S[1,1,2]).to_descent_algebra(4)
                    -3*B[1, 1, 2] + B[1, 2, 1]
                    sage: S[2,1].to_descent_algebra(2)
                    0
                    sage: S[2,1].to_descent_algebra()
                    B[2, 1]
                    sage: S.zero().to_descent_algebra().parent()
                    Descent algebra of 0 over Integer Ring in the subset basis
                    sage: (S[1,2,1] - 3 * S[1,1,2]).to_descent_algebra(1)
                    0
                """
                if n is None:
                    if self.is_zero():
                        n = 0
                    else:
                        n = self.degree()
                from sage.combinat.descent_algebra import DescentAlgebra
                S = NonCommutativeSymmetricFunctions(self.base_ring()).S()
                S_expansion = S(self)
                B = DescentAlgebra(self.base_ring(), n).B()
                return B.sum(coeff * B[I] for I, coeff in S_expansion if sum(I) == n)

            def to_symmetric_group_algebra(self):
                r"""
                Return the image of a non-commutative symmetric function into the
                symmetric group algebra where the ribbon basis element indexed by a
                composition is associated with the sum of all permutations which have
                descent set equal to said composition. In compliance with the anti-
                isomorphism between the descent algebra and the non-commutative
                symmetric functions, we index descent positions by the reversed composition.

                OUTPUT:

                - The image of ``self`` under the embedding of the `n`-th
                  degree homogeneous component of the non-commutative
                  symmetric functions in the symmetric group algebra of
                  `S_n`. This can behave unexpectedly when ``self`` is
                  not homogeneous.

                EXAMPLES::

                    sage: R=NonCommutativeSymmetricFunctions(QQ).R()
                    sage: R[2,1].to_symmetric_group_algebra()
                    [1, 3, 2] + [2, 3, 1]
                    sage: R([]).to_symmetric_group_algebra()
                    []

                TESTS:

                Sending a noncommutative symmetric function to the symmetric
                group algebra directly has the same result as going through
                the descent algebra::

                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: SGA4 = SymmetricGroupAlgebra(QQ, 4)
                    sage: D4 = DescentAlgebra(QQ, 4).D()
                    sage: all( S[C].to_symmetric_group_algebra()
                    ....:      == SGA4(D4(S[C].to_descent_algebra(4)))
                    ....:      for C in Compositions(4) )
                    True
                """
                S = NonCommutativeSymmetricFunctions(self.base_ring()).S()
                S_expansion = S(self)
                return sum(S_expansion.coefficient(I)*S._to_symmetric_group_algebra_on_basis(I) for I in S_expansion.support())
                # TODO:
                # This is ugly (uses global sum function) and undefined if self
                # is not homogeneous. Improve?

            def to_symmetric_function(self):
                r"""
                Return the commutative image of a non-commutative symmetric function.

                OUTPUT:

                - The commutative image of ``self``. This will be a symmetric function.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: x = R.an_element(); x
                    2*R[] + 2*R[1] + 3*R[1, 1]
                    sage: x.to_symmetric_function()
                    2*s[] + 2*s[1] + 3*s[1, 1]
                    sage: y = N.Phi()[1,3]
                    sage: y.to_symmetric_function()
                    h[1, 1, 1, 1] - 3*h[2, 1, 1] + 3*h[3, 1]
                """
                return self.parent().to_symmetric_function(self)

            chi = to_symmetric_function

            def to_ncsym(self):
                r"""
                Return the image of ``self`` under the injective
                algebra homomorphism `\kappa : NSym \to NCSym`
                that fixes the symmetric functions.

                As usual, `NCSym` denotes the ring of symmetric
                functions in non-commuting variables.
                Let `S_n` denote a generator of the complete basis.
                The algebra homomorphism `\kappa : NSym \to NCSym`
                is defined by

                .. MATH::

                    S_n \mapsto \sum_{A \vdash [n]}
                    \frac{\lambda(A)! \lambda(A)^!}{n!} \mathbf{m}_A .

                It has the property that the canonical maps
                `\chi : NCSym \to Sym` and `\rho : NSym \to Sym`
                satisfy `\chi \circ \kappa = \rho`.

                .. NOTE::

                    A remark in [BRRZ08]_  makes it clear that the embedding
                    of `NSym` into `NCSym` that preserves the projection into
                    the symmetric functions is not unique.  While this seems
                    to be a natural embedding, any free set of algebraic
                    generators of `NSym` can be sent to a set of free elements
                    in `NCSym` to form another embedding.

                .. SEEALSO::

                    :class:`NonCommutativeSymmetricFunctions` for a definition
                    of `NCSym`.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = N.complete()
                    sage: S[2].to_ncsym()
                    1/2*m{{1}, {2}} + m{{1, 2}}
                    sage: S[1,2,1].to_ncsym()
                    1/2*m{{1}, {2}, {3}, {4}} + 1/2*m{{1}, {2}, {3, 4}} + m{{1}, {2, 3}, {4}}
                     + m{{1}, {2, 3, 4}} + 1/2*m{{1}, {2, 4}, {3}} + 1/2*m{{1, 2}, {3}, {4}}
                     + 1/2*m{{1, 2}, {3, 4}} + m{{1, 2, 3}, {4}} + m{{1, 2, 3, 4}}
                     + 1/2*m{{1, 2, 4}, {3}} + 1/2*m{{1, 3}, {2}, {4}} + 1/2*m{{1, 3}, {2, 4}}
                     + 1/2*m{{1, 3, 4}, {2}} + 1/2*m{{1, 4}, {2}, {3}} + m{{1, 4}, {2, 3}}
                    sage: S[1,2].to_ncsym()
                    1/2*m{{1}, {2}, {3}} + m{{1}, {2, 3}} + 1/2*m{{1, 2}, {3}}
                     + m{{1, 2, 3}} + 1/2*m{{1, 3}, {2}}
                    sage: S[[]].to_ncsym()
                    m{}

                    sage: R = N.ribbon()
                    sage: x = R.an_element(); x
                    2*R[] + 2*R[1] + 3*R[1, 1]
                    sage: x.to_ncsym()
                    2*m{} + 2*m{{1}} + 3/2*m{{1}, {2}}
                    sage: R[2,1].to_ncsym()
                    1/3*m{{1}, {2}, {3}} + 1/6*m{{1}, {2, 3}}
                     + 2/3*m{{1, 2}, {3}} + 1/6*m{{1, 3}, {2}}

                    sage: Phi = N.Phi()
                    sage: Phi[1,2].to_ncsym()
                    m{{1}, {2, 3}} + m{{1, 2, 3}}
                    sage: Phi[1,3].to_ncsym()
                    -1/4*m{{1}, {2}, {3, 4}} - 1/4*m{{1}, {2, 3}, {4}} + m{{1}, {2, 3, 4}}
                     + 1/2*m{{1}, {2, 4}, {3}} - 1/4*m{{1, 2}, {3, 4}} - 1/4*m{{1, 2, 3}, {4}}
                     + m{{1, 2, 3, 4}} + 1/2*m{{1, 2, 4}, {3}} + 1/2*m{{1, 3}, {2, 4}}
                     - 1/4*m{{1, 3, 4}, {2}} - 1/4*m{{1, 4}, {2, 3}}

                TESTS:

                Check that `\chi \circ \kappa = \rho`::

                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                    sage: mon = SymmetricFunctions(QQ).monomial()
                    sage: all(S[c].to_ncsym().to_symmetric_function()
                    ....:     == S[c].to_symmetric_function()
                    ....:     for i in range(5) for c in Compositions(i))
                    True

                We also check that the `NCSym` monomials agree on the homogeneous
                and complete basis::

                    sage: h = SymmetricFunctions(QQ).h()
                    sage: all(m.from_symmetric_function(h[i])
                    ....:     == S[i].to_ncsym() for i in range(6))
                    True
                """
                from sage.combinat.ncsym.ncsym import SymmetricFunctionsNonCommutingVariables
                from sage.combinat.set_partition import SetPartitions
                P = self.parent()
                S = P.realization_of().complete()
                R = P.base_ring()
                m = SymmetricFunctionsNonCommutingVariables(R).m()
                SP = SetPartitions()

                def on_basis(I):
                    if not I._list:
                        return m.one()

                    def c_num(A):
                        return R(prod(factorial(i) for i in A.shape()))
                    return prod(m.sum_of_terms([(SP(A), R(c_num(A) / factorial(n)))
                                                for A in SetPartitions(n)], distinct=True)
                                for n in I)

                return m.linear_combination((on_basis(I), coeff)
                                            for I, coeff in S(self))

            def to_fqsym(self):
                r"""
                Return the image of the non-commutative symmetric
                function ``self`` under the morphism
                `\iota : NSym \to FQSym`.

                This morphism is the injective algebra homomorphism
                `NSym \to FQSym` that sends each Complete generator
                `S_n` to `F_{[1, 2, \ldots, n]}`. It is the inclusion
                map, if we regard both `NSym` and `FQSym` as rings of
                noncommutative power series.

                .. SEEALSO::

                    :class:`FreeQuasisymmetricFunctions` for a definition
                    of `FQSym`.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: x = 2*R[[]] + 2*R[1] + 3*R[2]
                    sage: x.to_fqsym()
                    2*F[] + 2*F[1] + 3*F[1, 2]
                    sage: R[2,1].to_fqsym()
                    F[1, 3, 2] + F[3, 1, 2]
                    sage: x = R.an_element(); x
                    2*R[] + 2*R[1] + 3*R[1, 1]
                    sage: x.to_fqsym()
                    2*F[] + 2*F[1] + 3*F[2, 1]

                    sage: y = N.Phi()[1,2]
                    sage: y.to_fqsym()
                    F[1, 2, 3] - F[1, 3, 2] + F[2, 1, 3] + F[2, 3, 1]
                     - F[3, 1, 2] - F[3, 2, 1]

                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S[2].to_fqsym()
                    F[1, 2]
                    sage: S[1,2].to_fqsym()
                    F[1, 2, 3] + F[2, 1, 3] + F[2, 3, 1]
                    sage: S[2,1].to_fqsym()
                    F[1, 2, 3] + F[1, 3, 2] + F[3, 1, 2]
                    sage: S[1,2,1].to_fqsym()
                    F[1, 2, 3, 4] + F[1, 2, 4, 3] + F[1, 4, 2, 3]
                     + F[2, 1, 3, 4] + F[2, 1, 4, 3] + F[2, 3, 1, 4]
                     + F[2, 3, 4, 1] + F[2, 4, 1, 3] + F[2, 4, 3, 1]
                     + F[4, 1, 2, 3] + F[4, 2, 1, 3] + F[4, 2, 3, 1]
                """
                from sage.combinat.fqsym import FreeQuasisymmetricFunctions
                P = self.parent()
                S = P.realization_of().complete()
                F = FreeQuasisymmetricFunctions(P.base_ring()).F()

                def on_basis(I):
                    return F.prod(F[Permutations(i)(range(1, i+1))] for i in I)
                return F.linear_combination((on_basis(I), coeff)
                                            for I, coeff in S(self))

            def to_fsym(self):
                r"""
                Return the image of ``self`` under the natural map to
                `FSym`.

                There is an injective Hopf algebra morphism from `NSym` to
                `FSym` (see
                :class:`~sage.combinat.chas.fsym.FreeSymmetricFunctions`),
                which maps the ribbon `R_\alpha` indexed by a composition
                `\alpha` to the sum of all tableaux whose descent
                composition is `\alpha`.
                If we regard `NSym` as a Hopf subalgebra of `FQSym` via
                the morphism `\iota : NSym \to FQSym` (implemented as
                :meth:`to_fqsym`), then this injective morphism is just
                the inclusion map.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: x = 2*R[[]] + 2*R[1] + 3*R[2]
                    sage: x.to_fsym()
                    2*G[] + 2*G[1] + 3*G[12]
                    sage: R[2,1].to_fsym()
                    G[12|3]
                    sage: R[1,2].to_fsym()
                    G[13|2]
                    sage: R[2,1,2].to_fsym()
                    G[12|35|4] + G[125|3|4]
                    sage: x = R.an_element(); x
                    2*R[] + 2*R[1] + 3*R[1, 1]
                    sage: x.to_fsym()
                    2*G[] + 2*G[1] + 3*G[1|2]

                    sage: y = N.Phi()[1,2]
                    sage: y.to_fsym()
                    -G[1|2|3] - G[12|3] + G[123] + G[13|2]

                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S[2].to_fsym()
                    G[12]
                    sage: S[2,1].to_fsym()
                    G[12|3] + G[123]
                """
                from sage.combinat.chas.fsym import FreeSymmetricFunctions
                G = FreeSymmetricFunctions(self.base_ring()).G()
                return G(self)

            def expand(self, n, alphabet='x'):
                r"""
                Expand the noncommutative symmetric function into an
                element of a free algebra in ``n`` indeterminates of
                an alphabet, which by default is ``'x'``.

                INPUT:

                - ``n`` -- a nonnegative integer; the number of variables
                  in the expansion
                - ``alphabet`` -- (default: ``'x'``); the alphabet in
                  which ``self`` is to be expanded

                OUTPUT:

                - An expansion of ``self`` into the ``n`` variables
                  specified by ``alphabet``.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = NSym.S()
                    sage: S[3].expand(3)
                    x0^3 + x0^2*x1 + x0^2*x2 + x0*x1^2 + x0*x1*x2
                     + x0*x2^2 + x1^3 + x1^2*x2 + x1*x2^2 + x2^3
                    sage: L = NSym.L()
                    sage: L[3].expand(3)
                    x2*x1*x0
                    sage: L[2].expand(3)
                    x1*x0 + x2*x0 + x2*x1
                    sage: L[3].expand(4)
                    x2*x1*x0 + x3*x1*x0 + x3*x2*x0 + x3*x2*x1
                    sage: Psi = NSym.Psi()
                    sage: Psi[2, 1].expand(3)
                    x0^3 + x0^2*x1 + x0^2*x2 + x0*x1*x0 + x0*x1^2 + x0*x1*x2
                     + x0*x2*x0 + x0*x2*x1 + x0*x2^2 - x1*x0^2 - x1*x0*x1
                     - x1*x0*x2 + x1^2*x0 + x1^3 + x1^2*x2 + x1*x2*x0
                     + x1*x2*x1 + x1*x2^2 - x2*x0^2 - x2*x0*x1 - x2*x0*x2
                     - x2*x1*x0 - x2*x1^2 - x2*x1*x2 + x2^2*x0 + x2^2*x1 + x2^3

                One can use a different set of variables by adding an optional
                argument ``alphabet=...``::

                    sage: L[3].expand(4, alphabet="y")
                    y2*y1*y0 + y3*y1*y0 + y3*y2*y0 + y3*y2*y1

                TESTS::

                    sage: (3*S([])).expand(2)
                    3
                    sage: L[4,2].expand(0)
                    0
                    sage: S([]).expand(0)
                    1
                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: S = NSym.S()
                    sage: S[3].expand(3)
                    x0^3 + x0^2*x1 + x0^2*x2 + x0*x1^2 + x0*x1*x2
                     + x0*x2^2 + x1^3 + x1^2*x2 + x1*x2^2 + x2^3

                .. TODO::

                    So far this is only implemented on the elementary
                    basis, and everything else goes through coercion.
                    Maybe it is worth shortcutting some of the other
                    bases?
                """
                NSym = self.parent().realization_of()
                L = NSym.L()
                from sage.algebras.free_algebra import FreeAlgebra
                P = FreeAlgebra(NSym.base_ring(), n, alphabet)
                x = P.gens()

                def image_of_L_k(k, i):
                    # Return the expansion of `L_k` (for `k` nonnegative
                    # integer) in the first `i` of the variables.
                    if k == 0:
                        return P.one()
                    if k > i:
                        return P.zero()
                    return x[i-1] * image_of_L_k(k - 1, i - 1) + image_of_L_k(k, i - 1)

                def on_basis(comp):
                    return P.prod((image_of_L_k(k, n) for k in comp))
                return L._apply_module_morphism(L(self), on_basis, codomain=P)

    class MultiplicativeBases(Category_realization_of_parent):
        """
        Category of multiplicative bases of non-commutative symmetric functions.

        EXAMPLES::

            sage: N = NonCommutativeSymmetricFunctions(QQ)
            sage: N.MultiplicativeBases()
            Category of multiplicative bases of Non-Commutative Symmetric Functions over the Rational Field

        The complete basis is a multiplicative basis, but the ribbon basis is not::

            sage: N.Complete() in N.MultiplicativeBases()
            True
            sage: N.Ribbon() in N.MultiplicativeBases()
            False
        """
        def super_categories(self):
            r"""
            Return the super categories of the category of multiplicative
            bases of the non-commutative symmetric functions.

            OUTPUT:

            - list

            TESTS::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: N.MultiplicativeBases().super_categories()
                [Category of bases of Non-Commutative Symmetric Functions over the Rational Field]

            """
            return [self.base().Bases()]

        class ParentMethods:

            @cached_method
            def algebra_generators(self):
                """
                Return the algebra generators of a given multiplicative basis of
                non-commutative symmetric functions.

                OUTPUT:

                - The family of generators of the multiplicative basis ``self``.

                EXAMPLES::

                    sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                    sage: f = Psi.algebra_generators()
                    sage: f
                    Lazy family (<lambda>(i))_{i in Positive integers}
                    sage: f[1], f[2], f[3]
                    (Psi[1], Psi[2], Psi[3])
                """
                from sage.sets.family import Family
                from sage.sets.positive_integers import PositiveIntegers
                return Family(PositiveIntegers(), lambda i: self.monomial(self._indices([i])))

            def product_on_basis(self, composition1, composition2):
                """
                Return the product of two basis elements from the multiplicative basis.
                Multiplication is just concatenation on compositions.

                INPUT:

                - ``composition1``, ``composition2`` -- integer compositions

                OUTPUT:

                - The product of the two non-commutative symmetric functions
                  indexed by ``composition1`` and ``composition2`` in the
                  multiplicative basis ``self``. This will be again
                  a non-commutative symmetric function.

                EXAMPLES::

                    sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                    sage: Psi[3,1,2] * Psi[4,2] == Psi[3,1,2,4,2]
                    True
                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S.product_on_basis(Composition([2,1]), Composition([1,2]))
                    S[2, 1, 1, 2]
                """
                return self.monomial(composition1 + composition2)

            def algebra_morphism(self, on_generators, **keywords):
                """
                Given a map defined on the generators of the multiplicative
                basis ``self``, return the algebra morphism that extends
                this map to the whole algebra of non-commutative symmetric
                functions.

                INPUT:

                - ``on_generators`` -- a function defined on the index set of
                  the generators (that is, on the positive integers)
                - ``anti`` -- a boolean; defaults to ``False``
                - ``category`` -- a category; defaults to ``None``

                OUTPUT:

                - The algebra morphism of ``self`` which is defined by
                  ``on_generators`` in the basis ``self``. When ``anti``
                  is set to ``True``, an algebra anti-morphism is
                  computed instead of an algebra morphism.

                EXAMPLES::

                    sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                    sage: Psi = NCSF.Psi()
                    sage: double = lambda i: Psi[i,i]
                    sage: f = Psi.algebra_morphism(double, codomain = Psi)
                    sage: f
                    Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
                    sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
                    2*Psi[] + 3*Psi[1, 1, 3, 3, 2, 2] + Psi[2, 2, 4, 4]
                    sage: f.category()
                    Category of endsets of unital magmas and right modules over Rational Field and left modules over Rational Field

                When extra properties about the morphism are known, one
                can specify the category of which it is a morphism::

                    sage: negate = lambda i: -Psi[i]
                    sage: f = Psi.algebra_morphism(negate, codomain = Psi, category = GradedHopfAlgebrasWithBasis(QQ))
                    sage: f
                    Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
                    sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
                    2*Psi[] - 3*Psi[1, 3, 2] + Psi[2, 4]
                    sage: f.category()
                    Category of endsets of hopf algebras over Rational Field and graded modules over Rational Field

                If ``anti`` is true, this returns an anti-algebra morphism::

                    sage: f = Psi.algebra_morphism(double, codomain = Psi, anti=True)
                    sage: f
                    Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
                    sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
                    2*Psi[] + 3*Psi[2, 2, 3, 3, 1, 1] + Psi[4, 4, 2, 2]
                    sage: f.category()
                    Category of endsets of modules with basis over Rational Field
                """
                from sage.combinat.ncsf_qsym.generic_basis_code import AlgebraMorphism
                return AlgebraMorphism(self, on_generators, **keywords)

            def to_symmetric_function_on_generators( self, i ):
                r"""
                Morphism of the generators to symmetric functions.

                This is constructed by coercion to the complete basis
                and applying the morphism.

                OUTPUT:

                - The module morphism from the basis ``self`` to the symmetric
                  functions which corresponds to taking a commutative image.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: Phi = N.Phi()
                    sage: Phi.to_symmetric_function_on_generators(3)
                    h[1, 1, 1] - 3*h[2, 1] + 3*h[3]
                    sage: Phi.to_symmetric_function_on_generators(0)
                    h[]
                    sage: Psi = N.Psi()
                    sage: Psi.to_symmetric_function_on_generators(3)
                    h[1, 1, 1] - 3*h[2, 1] + 3*h[3]
                    sage: L = N.elementary()
                    sage: L.to_symmetric_function_on_generators(3)
                    h[1, 1, 1] - 2*h[2, 1] + h[3]
                """
                S = self.realization_of().a_realization()
                if not i:
                    return S.to_symmetric_function_on_basis([])
                return S.to_symmetric_function(S(self([i])))

            @lazy_attribute
            def to_symmetric_function(self):
                r"""
                Morphism to the algebra of symmetric functions.

                This is constructed by extending the algebra morphism
                by the image of the generators.

                OUTPUT:

                - The module morphism from the basis ``self`` to the symmetric
                  functions which corresponds to taking a commutative image.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = N.complete()
                    sage: S.to_symmetric_function(S[1,3])
                    h[3, 1]
                    sage: Phi = N.Phi()
                    sage: Phi.to_symmetric_function(Phi[1,3])
                    h[1, 1, 1, 1] - 3*h[2, 1, 1] + 3*h[3, 1]
                    sage: Psi = N.Psi()
                    sage: Psi.to_symmetric_function(Psi[1,3])
                    h[1, 1, 1, 1] - 3*h[2, 1, 1] + 3*h[3, 1]
                """
                codom = self.to_symmetric_function_on_generators(1).parent()
                return self.algebra_morphism(self.to_symmetric_function_on_generators, codomain = codom)

            @lazy_attribute
            def antipode(self):
                r"""
                Return the antipode morphism on the basis ``self``.

                The antipode of `NSym` is closely related to the omega
                involution; see
                :meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.omega_involution`
                for the latter.

                OUTPUT:

                - The antipode module map from non-commutative symmetric
                  functions on basis ``self``.

                EXAMPLES::

                    sage: S=NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S.antipode
                    Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
                """
                if hasattr(self, "antipode_on_generators"):
                    return self.algebra_morphism(self.antipode_on_generators, codomain = self, anti = True)
                else:
                    return NotImplemented

            @lazy_attribute
            def coproduct(self):
                r"""
                Return the coproduct morphism in the basis ``self``.

                OUTPUT:

                - The coproduct module map from non-commutative symmetric
                  functions on basis ``self``.

                EXAMPLES::

                    sage: S=NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S.coproduct
                    Generic morphism:
                      From: Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
                      To:   Non-Commutative Symmetric Functions over the Rational Field in the Complete basis # Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
                """
                from sage.categories.all import tensor
                if hasattr(self, "coproduct_on_generators"):
                    return self.algebra_morphism(self.coproduct_on_generators, codomain = tensor([self, self]))
                else:
                    return NotImplemented

    class MultiplicativeBasesOnGroupLikeElements(Category_realization_of_parent):
        r"""
        Category of multiplicative bases on grouplike elements of
        non-commutative symmetric functions.

        Here, a "multiplicative basis on grouplike elements" means
        a multiplicative basis whose generators `(f_1, f_2, f_3, \ldots )`
        satisfy

        .. MATH::

            \Delta(f_i) = \sum_{j=0}^{i} f_j \otimes f_{i-j}

        with `f_0 = 1`. (In other words, the generators are to form a
        divided power sequence in the sense of a coalgebra.) This
        does not mean that the generators are grouplike, but means that
        the element `1 + f_1 + f_2 + f_3 + \cdots` in the completion of
        the ring of non-commutative symmetric functions with respect
        to the grading is grouplike.

        EXAMPLES::

            sage: N = NonCommutativeSymmetricFunctions(QQ)
            sage: N.MultiplicativeBasesOnGroupLikeElements()
            Category of multiplicative bases on group like elements of Non-Commutative Symmetric Functions over the Rational Field

        The complete basis is a multiplicative basis, but the ribbon basis is not::

            sage: N.Complete() in N.MultiplicativeBasesOnGroupLikeElements()
            True
            sage: N.Ribbon() in N.MultiplicativeBasesOnGroupLikeElements()
            False
        """
        def super_categories(self):
            r"""
            Return the super categories of the category of multiplicative
            bases of group-like elements of the non-commutative symmetric
            functions.

            OUTPUT:

            - list

            TESTS::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: N.MultiplicativeBasesOnGroupLikeElements().super_categories()
                [Category of multiplicative bases of Non-Commutative Symmetric Functions over the Rational Field]
            """
            return [self.base().MultiplicativeBases()]

        class ParentMethods:

            def antipode_on_basis(self, composition):
                """
                Return the application of the antipode to a basis element.

                INPUT:

                - ``composition`` -- a composition

                OUTPUT:

                - The image of the basis element indexed by ``composition``
                  under the antipode map.

                EXAMPLES::

                    sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                    sage: S.antipode_on_basis(Composition([2,1]))
                    -S[1, 1, 1] + S[1, 2]
                    sage: S[1].antipode() # indirect doctest
                    -S[1]
                    sage: S[2].antipode() # indirect doctest
                    S[1, 1] - S[2]
                    sage: S[3].antipode() # indirect doctest
                    -S[1, 1, 1] + S[1, 2] + S[2, 1] - S[3]
                    sage: S[2,3].coproduct().apply_multilinear_morphism(lambda be,ga: S(be)*S(ga).antipode())
                    0
                    sage: S[2,3].coproduct().apply_multilinear_morphism(lambda be,ga: S(be).antipode()*S(ga))
                    0
                """
                return (-1)**len(composition) * self.alternating_sum_of_finer_compositions(composition.reversed())

            # @cached_method?
            def coproduct_on_generators(self, i):
                r"""
                Return the image of the `i^{th}` generator of the algebra under
                the coproduct.

                INPUT:

                - ``i`` -- a positive integer

                OUTPUT:

                - The result of applying the coproduct to the `i^{th}`
                  generator of ``self``.

                EXAMPLES::

                    sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                    sage: S.coproduct_on_generators(3)
                    S[] # S[3] + S[1] # S[2] + S[2] # S[1] + S[3] # S[]

                TESTS::

                    sage: S.coproduct_on_generators(0)
                    Traceback (most recent call last):
                    ...
                    ValueError: Not a positive integer: 0
                """
                if i < 1:
                    raise ValueError("Not a positive integer: {}".format(i))

                def C(i):
                    return self._indices([i]) if i else self._indices([])
                T = self.tensor_square()
                return T.sum_of_monomials( (C(j), C(i-j)) for j in range(i+1) )

    class MultiplicativeBasesOnPrimitiveElements(Category_realization_of_parent):
        r"""
        Category of multiplicative bases of the non-commutative symmetric
        functions whose generators are primitive elements.

        An element `x` of a bialgebra is *primitive* if
        `\Delta(x) = x \otimes 1 + 1 \otimes x`, where
        `\Delta` is the coproduct of the bialgebra.

        Given a multiplicative basis and knowing the coproducts and antipodes
        of its generators, one can compute the coproduct and the antipode of
        any element, since they are respectively algebra morphisms and
        anti-morphisms. See :meth:`~ParentMethods.antipode_on_generators` and
        :meth:`~ParentMethods.coproduct_on_generators`.

        .. TODO:: this could be generalized to any free algebra.

        EXAMPLES::

            sage: N = NonCommutativeSymmetricFunctions(QQ)
            sage: N.MultiplicativeBasesOnPrimitiveElements()
            Category of multiplicative bases on primitive elements of Non-Commutative Symmetric Functions over the Rational Field

        The Phi and Psi bases are multiplicative bases whose generators
        are primitive elements, but the complete and ribbon bases are not::

            sage: N.Phi() in N.MultiplicativeBasesOnPrimitiveElements()
            True
            sage: N.Psi() in N.MultiplicativeBasesOnPrimitiveElements()
            True
            sage: N.Complete() in N.MultiplicativeBasesOnPrimitiveElements()
            False
            sage: N.Ribbon() in N.MultiplicativeBasesOnPrimitiveElements()
            False
        """

        def super_categories(self):
            r"""
            Return the super categories of the category of multiplicative
            bases of primitive elements of the non-commutative symmetric
            functions.

            OUTPUT:

            - list

            TESTS::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: N.MultiplicativeBasesOnPrimitiveElements().super_categories()
                [Category of multiplicative bases of Non-Commutative Symmetric Functions over the Rational Field]
            """
            return [self.base().MultiplicativeBases()]

        class ParentMethods:

            def antipode_on_generators(self, i):
                r"""
                Return the image of a generator of a primitive basis of
                the non-commutative symmetric functions under the antipode
                map.

                INPUT:

                - ``i`` -- a positive integer

                OUTPUT:

                - The image of the `i`-th generator of the multiplicative
                  basis ``self`` under the antipode of the algebra of
                  non-commutative symmetric functions.

                EXAMPLES::

                    sage: Psi=NonCommutativeSymmetricFunctions(QQ).Psi()
                    sage: Psi.antipode_on_generators(2)
                    -Psi[2]

                TESTS::

                    sage: Psi.antipode_on_generators(0)
                    Traceback (most recent call last):
                    ...
                    ValueError: Not a positive integer: 0
                """
                if i < 1:
                    raise ValueError("Not a positive integer: {}".format(i))
                return - self.algebra_generators()[i]

            def coproduct_on_generators(self, i):
                r"""
                Return the image of the `i^{th}` generator of the
                multiplicative basis ``self`` under the coproduct.

                INPUT:

                - ``i`` -- a positive integer

                OUTPUT:

                - The result of applying the coproduct to the
                  `i^{th}` generator of ``self``.

                EXAMPLES::

                    sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                    sage: Psi.coproduct_on_generators(3)
                    Psi[] # Psi[3] + Psi[3] # Psi[]

                TESTS::

                    sage: Psi.coproduct_on_generators(0)
                    Traceback (most recent call last):
                    ...
                    ValueError: Not a positive integer: 0
                """
                if i < 1:
                    raise ValueError("Not a positive integer: {}".format(i))
                x = self.algebra_generators()[i]
                from sage.categories.all import tensor
                return tensor([self.one(), x]) + tensor([x, self.one()])

    class Ribbon(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of non-commutative symmetric functions in the
        Ribbon basis.

        The Ribbon basis is defined in Definition 3.12 of [NCSF1]_, where
        it is denoted by `(R_I)_I`. It is connected to the complete
        basis of the ring of non-commutative symmetric functions by the
        following relation: For every composition `I`, we have

        .. MATH::

            R_I = \sum_J (-1)^{\ell(I) - \ell(J)} S^J,

        where the sum is over all compositions `J` which are coarser than
        `I` and `\ell(I)` denotes the length of `I`. (See the proof of
        Proposition 4.13 in [NCSF1]_.)

        The elements of the Ribbon basis are commonly referred to as the
        ribbon Schur functions.

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: R = NCSF.Ribbon(); R
            Non-Commutative Symmetric Functions over the Rational Field in the Ribbon basis
            sage: R.an_element()
            2*R[] + 2*R[1] + 3*R[1, 1]

        The following are aliases for this basis::

            sage: NCSF.ribbon()
            Non-Commutative Symmetric Functions over the Rational Field in the Ribbon basis
            sage: NCSF.R()
            Non-Commutative Symmetric Functions over the Rational Field in the Ribbon basis
        """
        def __init__(self, NCSF):
            r"""
            EXAMPLES::

                sage: R = NonCommutativeSymmetricFunctions(ZZ).ribbon(); R
                Non-Commutative Symmetric Functions over the Integer Ring in the Ribbon basis
                sage: TestSuite(R).run()

            TESTS::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: R.one().coproduct()
                R[] # R[]
                sage: R.one().antipode()
                R[]

            We include some sanity tests to verify that conversions between the
            ribbon basis and other bases work the way they should::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: all(S(R(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(R(S(R[comp])) == R[comp] for comp in Compositions(5))
                True

            ::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: L = NonCommutativeSymmetricFunctions(QQ).elementary()
                sage: all(L(R(L[comp])) == L[comp] for comp in Compositions(5))
                True
                sage: all(R(L(R[comp])) == R[comp] for comp in Compositions(5))
                True
            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='R', bracket=False,
                                             category=NCSF.Bases())

        def dual(self):
            r"""
            Return the dual basis to the ribbon basis of the non-commutative symmetric
            functions. This is the Fundamental basis of the quasi-symmetric functions.

            OUTPUT:

            - The fundamental basis of the quasi-symmetric functions.

            EXAMPLES::

                sage: R=NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: R.dual()
                Quasisymmetric functions over the Rational Field in the Fundamental basis
            """
            return self.realization_of().dual().Fundamental()

        def product_on_basis(self, I, J):
            r"""
            Return the product of two ribbon basis elements of the non-commutative
            symmetric functions.

            INPUT:

            - ``I``, ``J`` -- compositions

            OUTPUT:

            - The product of the ribbon functions indexed by ``I`` and ``J``.

            EXAMPLES::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: R[1,2,1] * R[3,1]
                R[1, 2, 1, 3, 1] + R[1, 2, 4, 1]
                sage: ( R[1,2] + R[3] ) * ( R[3,1] + R[1,2,1] )
                R[1, 2, 1, 2, 1] + R[1, 2, 3, 1] + R[1, 3, 2, 1] + R[1, 5, 1] + R[3, 1, 2, 1] + R[3, 3, 1] + R[4, 2, 1] + R[6, 1]

            TESTS::

                sage: R[[]] * R[3,1]
                R[3, 1]
                sage: R[1,2,1] * R[[]]
                R[1, 2, 1]
                sage: R.product_on_basis(Composition([2,1]), Composition([1]))
                R[2, 1, 1] + R[2, 2]
            """
            if not I._list:
                return self.monomial(J)
            elif not J._list:
                return self.monomial(I)
            else:
                return self.monomial(self._indices(I[:] + J[:])) + \
                       self.monomial(self._indices(I[:-1] + [I[-1]+J[0]] + J[1:]))

        def antipode_on_basis(self, composition):
            """
            Return the application of the antipode to a basis element
            of the ribbon basis ``self``.

            INPUT:

            - ``composition`` -- a composition

            OUTPUT:

            - The image of the basis element indexed by ``composition``
              under the antipode map.

            EXAMPLES::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: R.antipode_on_basis(Composition([2,1]))
                -R[2, 1]
                sage: R[3,1].antipode() # indirect doctest
                R[2, 1, 1]
                sage: R[[]].antipode() # indirect doctest
                R[]

            We check that the implementation of the antipode at hand does
            not contradict the generic one::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: all( S(R[I].antipode()) == S(R[I]).antipode()
                ....:      for I in Compositions(4) )
                True
            """
            if composition.size() % 2 == 0:
                return self[composition.conjugate()]
            else:
                return - self[composition.conjugate()]

        def to_symmetric_function_on_basis(self, I):
            r"""
            Return the commutative image of a ribbon basis element of the
            non-commutative symmetric functions.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The commutative image of the ribbon basis element indexed by
              ``I``. This will be expressed as a symmetric function in the
              Schur basis.

            EXAMPLES::

                sage: R=NonCommutativeSymmetricFunctions(QQ).R()
                sage: R.to_symmetric_function_on_basis(Composition([3,1,1]))
                s[3, 1, 1]
                sage: R.to_symmetric_function_on_basis(Composition([4,2,1]))
                s[4, 2, 1] + s[5, 1, 1] + s[5, 2]
                sage: R.to_symmetric_function_on_basis(Composition([]))
                s[]
            """
            s = SymmetricFunctions(self.base_ring()).schur()
            if not I:
                return s([])
            return s(I.to_skew_partition())

        class Element(CombinatorialFreeModule.Element):

            def verschiebung(self, n):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the `n`-th Verschiebung operator.

                The `n`-th Verschiebung operator `\mathbf{V}_n` is defined
                to be the map from the `\mathbf{k}`-algebra of noncommutative
                symmetric functions to itself that sends the complete function
                `S^I` indexed by a composition `I = (i_1, i_2, \ldots , i_k)`
                to `S^{(i_1/n, i_2/n, \ldots , i_k/n)}` if all of the numbers
                `i_1, i_2, \ldots, i_k` are divisible by `n`, and to `0`
                otherwise. This operator `\mathbf{V}_n` is a Hopf algebra
                endomorphism. For every positive integer `r` with `n \mid r`,
                it satisfies

                .. MATH::

                    \mathbf{V}_n(S_r) = S_{r/n},
                    \quad \mathbf{V}_n(\Lambda_r) = (-1)^{r - r/n} \Lambda_{r/n},
                    \quad \mathbf{V}_n(\Psi_r) = n \Psi_{r/n},
                    \quad \mathbf{V}_n(\Phi_r) = n \Phi_{r/n}

                (where `S_r` denotes the `r`-th complete non-commutative
                symmetric function, `\Lambda_r` denotes the `r`-th elementary
                non-commutative symmetric function, `\Psi_r` denotes the `r`-th
                power-sum non-commutative symmetric function of the first kind,
                and `\Phi_r` denotes the `r`-th power-sum non-commutative
                symmetric function of the second kind). For every positive
                integer `r` with `n \nmid r`, it satisfes

                .. MATH::

                    \mathbf{V}_n(S_r) = \mathbf{V}_n(\Lambda_r)
                    = \mathbf{V}_n(\Psi_r) = \mathbf{V}_n(\Phi_r) = 0.

                The `n`-th Verschiebung operator is also called the `n`-th
                Verschiebung endomorphism.

                It is a lift of the `n`-th Verschiebung operator on the ring
                of symmetric functions
                (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung`)
                to the ring of noncommutative symmetric functions.

                The action of the `n`-th Verschiebung operator can also be
                described on the ribbon Schur functions. Namely, every
                composition `I` of size `n \ell` satisfies

                .. MATH::

                    \mathbf{V}_n ( R_I )
                    = (-1)^{\ell(I) - \ell(J)} \cdot R_{J / n},

                where `J` denotes the meet of the compositions `I` and
                `(\underbrace{n, n, \ldots, n}_{|I|/n \mbox{ times}})`,
                where `\ell(I)` is the length of `I`, and
                where `J / n` denotes the composition obtained by
                dividing every entry of `J` by `n`.
                For a composition `I` of size not divisible by `n`, we
                have `\mathbf{V}_n ( R_I ) = 0`.

                .. SEEALSO::

                    :meth:`verschiebung method of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.verschiebung>`,
                    :meth:`frobenius method of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.frobenius>`,
                    :meth:`verschiebung method of Sym
                    <sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung>`

                INPUT:

                - ``n`` -- a positive integer

                OUTPUT:

                The result of applying the `n`-th Verschiebung operator (on the
                ring of noncommutative symmetric functions) to ``self``.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: R = NSym.R()
                    sage: R([4,2]).verschiebung(2)
                    R[2, 1]
                    sage: R([2,1]).verschiebung(3)
                    -R[1]
                    sage: R([3]).verschiebung(2)
                    0
                    sage: R([]).verschiebung(2)
                    R[]
                    sage: R([5, 1]).verschiebung(3)
                    -R[2]
                    sage: R([5, 1]).verschiebung(6)
                    -R[1]
                    sage: R([5, 1]).verschiebung(2)
                    -R[3]
                    sage: R([1, 2, 3, 1]).verschiebung(7)
                    -R[1]
                    sage: R([1, 2, 3, 1]).verschiebung(5)
                    0
                    sage: (R[1] - R[2] + 2*R[3]).verschiebung(1)
                    R[1] - R[2] + 2*R[3]

                TESTS:

                The current implementation on the ribbon basis gives the
                same results as the default implementation::

                    sage: S = NSym.S()
                    sage: def test_ribbon(N, n):
                    ....:     for I in Compositions(N):
                    ....:         if S(R[I].verschiebung(n)) != S(R[I]).verschiebung(n):
                    ....:             return False
                    ....:     return True
                    sage: test_ribbon(4, 2)
                    True
                    sage: test_ribbon(6, 2)
                    True
                    sage: test_ribbon(6, 3)
                    True
                    sage: test_ribbon(8, 4)     # long time
                    True
                """
                parent = self.parent()
                C = parent._indices

                def ribbon_mapper(I, coeff):
                    # return \mathbf{V}_n ( coeff * R_I ) as pair
                    # (composition, coefficient)
                    M = sum(I)
                    m = M // n
                    J = I.meet([n] * m)
                    Jn = C([j // n for j in J])
                    if (len(I) - len(J)) % 2:
                        return (Jn, - coeff)
                    else:
                        return (Jn, coeff)
                return parent.sum_of_terms([ribbon_mapper(I, coeff)
                                            for (I, coeff) in self
                                            if sum(I) % n == 0])

            def star_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the star involution.

                The star involution is defined as the algebra antihomomorphism
                `NCSF \to NCSF` which, for every positive integer `n`, sends
                the `n`-th complete non-commutative symmetric function `S_n` to
                `S_n`. Denoting by `f^{\ast}` the image of an element
                `f \in NCSF` under this star involution, it can be shown that
                every composition `I` satisfies

                .. MATH::

                    (S^I)^{\ast} = S^{I^r}, \quad
                    (\Lambda^I)^{\ast} = \Lambda^{I^r}, \quad
                    R_I^{\ast} = R_{I^r}, \quad
                    (\Phi^I)^{\ast} = \Phi^{I^r},

                where `I^r` denotes the reversed composition of `I`, and
                standard notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary basis,
                `R` for the ribbon basis, and `\Phi` for that of the power-sums
                of the second kind). The star involution is an involution and a
                coalgebra automorphism of `NCSF`. It is an automorphism of the
                graded vector space `NCSF`. Under the canonical isomorphism
                between the `n`-th graded component of `NCSF` and the descent
                algebra of the symmetric group `S_n` (see
                :meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_descent_algebra`),
                the star involution (restricted to
                the `n`-th graded component) corresponds to the automorphism
                of the descent algebra given by
                `x \mapsto \omega_n x \omega_n`, where `\omega_n` is the
                permutation `(n, n-1, \ldots, 1) \in S_n` (written here in
                one-line notation). If `\pi` denotes the projection from `NCSF`
                to the ring of symmetric functions
                (:meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_symmetric_function`),
                then `\pi(f^{\ast}) = \pi(f)` for every `f \in NCSF`.

                The star involution on `NCSF` is adjoint to the star involution
                on `QSym` by the standard adjunction between `NCSF` and `QSym`.

                The star involution has been denoted by `\rho` in [LMvW13]_,
                section 3.6.
                See [NCSF2]_, section 2.3 for the properties of this map.

                .. SEEALSO::

                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`,
                    :meth:`star involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`,
                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: R = NSym.R()
                    sage: R[3,1,4,2].star_involution()
                    R[2, 4, 1, 3]
                    sage: R[4,1,2].star_involution()
                    R[2, 1, 4]
                    sage: (R[1] - R[2] + 2*R[5,4] - 3*R[3] + 4*R[[]]).star_involution()
                    4*R[] + R[1] - R[2] - 3*R[3] + 2*R[4, 5]
                    sage: (R[3,3] - 21*R[1]).star_involution()
                    -21*R[1] + R[3, 3]
                    sage: R([14,1]).star_involution()
                    R[1, 14]

                The implementation at hand is tailored to the ribbon basis.
                It is equivalent to the generic implementation via the
                complete basis::

                    sage: S = NSym.S()
                    sage: all( S(R[I].star_involution()) == S(R[I]).star_involution()
                    ....:      for I in Compositions(4) )
                    True
                """
                parent = self.parent()
                dct = {I.reversed(): coeff for (I, coeff) in self}
                return parent._from_dict(dct)

    R = ribbon = Ribbon

    class Complete(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of non-commutative symmetric functions in the
        Complete basis.

        The Complete basis is defined in Definition 3.4 of [NCSF1]_, where
        it is denoted by `(S^I)_I`. It is a multiplicative basis, and is
        connected to the elementary generators `\Lambda_i` of the ring of
        non-commutative symmetric functions by the following relation:
        Define a non-commutative symmetric function `S_n` for every
        nonnegative integer `n` by the power series identity

        .. MATH::

            \sum_{k \geq 0} t^k S_k
            = \left( \sum_{k \geq 0} (-t)^k \Lambda_k \right)^{-1},

        with `\Lambda_0` denoting `1`. For every composition
        `(i_1, i_2, \ldots, i_k)`, we have
        `S^{(i_1, i_2, \ldots, i_k)} = S_{i_1} S_{i_2} \cdots S_{i_k}`.

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: S = NCSF.Complete(); S
            Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
            sage: S.an_element()
            2*S[] + 2*S[1] + 3*S[1, 1]

        The following are aliases for this basis::

            sage: NCSF.complete()
            Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
            sage: NCSF.S()
            Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
        """
        def __init__(self, NCSF):
            r"""
            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: TestSuite(S).run()
            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='S', bracket=False,
                                             category=NCSF.MultiplicativeBasesOnGroupLikeElements())

        def dual(self):
            r"""
            Return the dual basis to the complete basis of non-commutative symmetric
            functions. This is the Monomial basis of quasi-symmetric functions.

            OUTPUT:

            - The Monomial basis of quasi-symmetric functions.

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.dual()
                Quasisymmetric functions over the Rational Field in the Monomial basis
            """
            return self.realization_of().dual().Monomial()

        def internal_product_on_basis(self, I, J):
            r"""
            The internal product of two non-commutative symmetric
            complete functions.

            See :meth:`~sage.combinat.ncsf_qsym.generic_basis_code.GradedModulesWithInternalProduct.ElementMethods.internal_product`
            for a thorough documentation of this operation.

            INPUT:

            - ``I``, ``J`` -- compositions

            OUTPUT:

            - The internal product of the complete non-commutative symmetric
              function basis elements indexed by ``I`` and ``J``, expressed in
              the complete basis.

            EXAMPLES::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: S = N.complete()
                sage: S.internal_product_on_basis([2,2],[1,2,1])
                2*S[1, 1, 1, 1] + S[1, 1, 2] + S[2, 1, 1]
                sage: S.internal_product_on_basis([2,2],[1,2])
                0
            """
            from sage.combinat.integer_matrices import IntegerMatrices
            IM = IntegerMatrices(I, J)
            return self.sum_of_monomials(IM.to_composition(m) for m in IM)

        def to_symmetric_function_on_basis(self, I):
            r"""
            The commutative image of a complete element

            The commutative image of a basis element is obtained by sorting
            the indexing composition of the basis element and the output
            is in the complete basis of the symmetric functions.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The commutative image of the complete basis element
              indexed by ``I``. The result is the complete symmetric function
              indexed by the partition obtained by sorting ``I``.

            EXAMPLES::

                sage: S=NonCommutativeSymmetricFunctions(QQ).complete()
                sage: S.to_symmetric_function_on_basis([2,1,3])
                h[3, 2, 1]
                sage: S.to_symmetric_function_on_basis([])
                h[]
            """
            h = SymmetricFunctions(self.base_ring()).complete()
            return h[Partition(sorted(I,reverse=True))]

        @lazy_attribute
        def to_symmetric_function(self):
            r"""
            Morphism to the algebra of symmetric functions.

            This is constructed by extending the computation on the
            complete basis.

            OUTPUT:

            - The module morphism from the basis ``self`` to the symmetric
              functions which corresponds to taking a commutative image.

            EXAMPLES::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: S = N.complete()
                sage: S.to_symmetric_function(S[3,1,2])
                h[3, 2, 1]
                sage: S.to_symmetric_function(S[[]])
                h[]
            """
            on_basis = self.to_symmetric_function_on_basis
            codom = SymmetricFunctions(self.base_ring()).complete()
            return self.module_morphism(on_basis, codomain=codom)

        def _to_symmetric_group_algebra_on_basis(self, I):
            r"""
            Return the image of the complete non-commutative symmetric function indexed
            by the composition ``I`` in the symmetric group algebra under the canonical
            embedding of the degree-`|I|` homogeneous non-commutative symmetric
            functions into the `|I|`-th symmetric group algebra.

            This embedding sends the complete basis element indexed by the composition
            ``I`` to the sum of all permutations whose descent composition is fatter
            than ``I`` (that is, all permutations whose right descent set is contained
            in the subset corresponding to ``I``).

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The sum of all permutations of `\{ 1, 2, \ldots, |I| \}` with right
              descent set contained in ``I``.

            EXAMPLES::

                sage: S=NonCommutativeSymmetricFunctions(QQ).S()
                sage: S._to_symmetric_group_algebra_on_basis(Composition([1,2]))
                [1, 2, 3] + [2, 1, 3] + [3, 1, 2]
                sage: S._to_symmetric_group_algebra_on_basis(Composition([]))
                []
                sage: S._to_symmetric_group_algebra_on_basis(Composition([1]))
                [1]
            """
            n = sum(I)
            from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
            from sage.sets.set import Set
            if n == 0:
                return SymmetricGroupAlgebra(self.base_ring(),n).one()
            sga = SymmetricGroupAlgebra(self.base_ring(),n)
            J = [j-1 for j in I.to_subset()]
            return sga.sum_of_monomials( p for K in Set(J).subsets()
                                         for p in Permutations(descents=(K,n)) )

        class Element(CombinatorialFreeModule.Element):
            """
            An element in the Complete basis.
            """
            def psi_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the involution `\psi`.

                The involution `\psi` is defined as the linear map
                `NCSF \to NCSF` which, for every composition `I`, sends the
                complete noncommutative symmetric function `S^I` to the
                elementary noncommutative symmetric function `\Lambda^I`.
                It can be shown that every composition `I` satisfies

                .. MATH::

                    \psi(R_I) = R_{I^c}, \quad \psi(S^I) = \Lambda^I, \quad
                    \psi(\Lambda^I) = S^I, \quad
                    \psi(\Phi^I) = (-1)^{|I| - \ell(I)} \Phi^I

                where `I^c` denotes the complement of the composition `I`, and
                `\ell(I)` denotes the length of `I`, and where standard
                notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary
                basis, `\Phi` for the basis of the power sums of the second
                kind, and `R` for the ribbon basis). The map `\psi` is an
                involution and a graded Hopf algebra automorphism of `NCSF`.
                If `\pi` denotes the projection from `NCSF` to the ring of
                symmetric functions
                (:meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_symmetric_function`),
                then `\pi(\psi(f)) = \omega(\pi(f))` for every `f \in NCSF`, where
                the `\omega` on the right hand side denotes the omega
                automorphism of `Sym`.

                The involution `\psi` of `NCSF` is adjoint to the involution
                `\psi` of `QSym` by the standard adjunction between `NCSF` and
                `QSym`.

                The involution `\psi` has been denoted by `\psi` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`psi involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: S = NSym.S()
                    sage: L = NSym.L()
                    sage: S[3,1].psi_involution()
                    S[1, 1, 1, 1] - S[1, 2, 1] - S[2, 1, 1] + S[3, 1]
                    sage: L(S[3,1].psi_involution())
                    L[3, 1]
                    sage: S[[]].psi_involution()
                    S[]
                    sage: S[1,1].psi_involution()
                    S[1, 1]
                    sage: (S[2,1] - 2*S[2]).psi_involution()
                    -2*S[1, 1] + S[1, 1, 1] + 2*S[2] - S[2, 1]

                The implementation at hand is tailored to the complete basis.
                It is equivalent to the generic implementation via the
                ribbon basis::

                    sage: R = NSym.R()
                    sage: all( R(S[I].psi_involution()) == R(S[I]).psi_involution()
                    ....:      for I in Compositions(4) )
                    True
                """
                parent = self.parent()
                return parent.sum( (-1) ** (I.size() - len(I)) * coeff
                                   * parent.alternating_sum_of_finer_compositions(I)
                                   for I, coeff in self._monomial_coefficients.items() )

    S = complete = Complete

    class Elementary(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of non-commutative symmetric functions in the
        Elementary basis.

        The Elementary basis is defined in Definition 3.4 of [NCSF1]_,
        where it is denoted by `(\Lambda^I)_I`. It is a multiplicative
        basis, and is obtained from the elementary generators
        `\Lambda_i` of the ring of non-commutative symmetric functions
        through the formula
        `\Lambda^{(i_1, i_2, \ldots, i_k)}
        = \Lambda_{i_1} \Lambda_{i_2} \cdots \Lambda_{i_k}`
        for every composition `(i_1, i_2, \ldots, i_k)`.

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: L = NCSF.Elementary(); L
            Non-Commutative Symmetric Functions over the Rational Field in the Elementary basis
            sage: L.an_element()
            2*L[] + 2*L[1] + 3*L[1, 1]

        The following are aliases for this basis::

            sage: NCSF.elementary()
            Non-Commutative Symmetric Functions over the Rational Field in the Elementary basis
            sage: NCSF.L()
            Non-Commutative Symmetric Functions over the Rational Field in the Elementary basis
        """
        def __init__(self, NCSF):
            r"""
            EXAMPLES::

                sage: NonCommutativeSymmetricFunctions(ZZ).elementary()
                Non-Commutative Symmetric Functions over the Integer Ring in the Elementary basis

            TESTS::

                sage: L = NonCommutativeSymmetricFunctions(ZZ).elementary()
                sage: TestSuite(L).run()

            We include a sanity test to verify the conversion between the
            elementary and complete basis works the way it should::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: L = NonCommutativeSymmetricFunctions(QQ).elementary()
                sage: L(S[3])
                L[1, 1, 1] - L[1, 2] - L[2, 1] + L[3]
                sage: L(S[2,1])
                L[1, 1, 1] - L[2, 1]
                sage: L(S[1,2])
                L[1, 1, 1] - L[1, 2]
                sage: L(S[1,1,1])
                L[1, 1, 1]
                sage: S(L[3])
                S[1, 1, 1] - S[1, 2] - S[2, 1] + S[3]
                sage: S(L[2,1])
                S[1, 1, 1] - S[2, 1]
                sage: S(L[1,2])
                S[1, 1, 1] - S[1, 2]
                sage: S(L[1,1,1])
                S[1, 1, 1]
                sage: all(S(L(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(L(S(L[comp])) == L[comp] for comp in Compositions(5))
                True
            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='L', bracket=False,
                                             category=NCSF.MultiplicativeBasesOnGroupLikeElements())

        class Element(CombinatorialFreeModule.Element):

            def verschiebung(self, n):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the `n`-th Verschiebung operator.

                The `n`-th Verschiebung operator `\mathbf{V}_n` is defined
                to be the map from the `\mathbf{k}`-algebra of noncommutative
                symmetric functions to itself that sends the complete function
                `S^I` indexed by a composition `I = (i_1, i_2, \ldots , i_k)`
                to `S^{(i_1/n, i_2/n, \ldots , i_k/n)}` if all of the numbers
                `i_1, i_2, \ldots, i_k` are divisible by `n`, and to `0`
                otherwise. This operator `\mathbf{V}_n` is a Hopf algebra
                endomorphism. For every positive integer `r` with `n \mid r`,
                it satisfies

                .. MATH::

                    \mathbf{V}_n(S_r) = S_{r/n},
                    \quad \mathbf{V}_n(\Lambda_r) = (-1)^{r - r/n} \Lambda_{r/n},
                    \quad \mathbf{V}_n(\Psi_r) = n \Psi_{r/n},
                    \quad \mathbf{V}_n(\Phi_r) = n \Phi_{r/n}

                (where `S_r` denotes the `r`-th complete non-commutative
                symmetric function, `\Lambda_r` denotes the `r`-th elementary
                non-commutative symmetric function, `\Psi_r` denotes the `r`-th
                power-sum non-commutative symmetric function of the first kind,
                and `\Phi_r` denotes the `r`-th power-sum non-commutative
                symmetric function of the second kind). For every positive
                integer `r` with `n \nmid r`, it satisfes

                .. MATH::

                    \mathbf{V}_n(S_r) = \mathbf{V}_n(\Lambda_r)
                    = \mathbf{V}_n(\Psi_r) = \mathbf{V}_n(\Phi_r) = 0.

                The `n`-th Verschiebung operator is also called the `n`-th
                Verschiebung endomorphism.

                It is a lift of the `n`-th Verschiebung operator on the ring
                of symmetric functions
                (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung`)
                to the ring of noncommutative symmetric functions.

                The action of the `n`-th Verschiebung operator can also be
                described on the ribbon Schur functions. Namely, every
                composition `I` of size `n \ell` satisfies

                .. MATH::

                    \mathbf{V}_n ( R_I )
                    = (-1)^{\ell(I) - \ell(J)} \cdot R_{J / n},

                where `J` denotes the meet of the compositions `I` and
                `(\underbrace{n, n, \ldots, n}_{|I|/n \mbox{ times}})`,
                where `\ell(I)` is the length of `I`, and
                where `J / n` denotes the composition obtained by
                dividing every entry of `J` by `n`.
                For a composition `I` of size not divisible by `n`, we
                have `\mathbf{V}_n ( R_I ) = 0`.

                .. SEEALSO::

                    :meth:`verschiebung method of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.verschiebung>`,
                    :meth:`frobenius method of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.frobenius>`,
                    :meth:`verschiebung method of Sym
                    <sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung>`

                INPUT:

                - ``n`` -- a positive integer

                OUTPUT:

                The result of applying the `n`-th Verschiebung operator (on the
                ring of noncommutative symmetric functions) to ``self``.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: L = NSym.L()
                    sage: L([4,2]).verschiebung(2)
                    -L[2, 1]
                    sage: L([2,4]).verschiebung(2)
                    -L[1, 2]
                    sage: L([6]).verschiebung(2)
                    -L[3]
                    sage: L([2,1]).verschiebung(3)
                    0
                    sage: L([3]).verschiebung(2)
                    0
                    sage: L([]).verschiebung(2)
                    L[]
                    sage: L([5, 1]).verschiebung(3)
                    0
                    sage: L([5, 1]).verschiebung(6)
                    0
                    sage: L([5, 1]).verschiebung(2)
                    0
                    sage: L([1, 2, 3, 1]).verschiebung(7)
                    0
                    sage: L([7]).verschiebung(7)
                    L[1]
                    sage: L([1, 2, 3, 1]).verschiebung(5)
                    0
                    sage: (L[1] - L[2] + 2*L[3]).verschiebung(1)
                    L[1] - L[2] + 2*L[3]

                TESTS:

                The current implementation on the Elementary basis gives
                the same results as the default implementation::

                    sage: S = NSym.S()
                    sage: def test_L(N, n):
                    ....:     for I in Compositions(N):
                    ....:         if S(L[I].verschiebung(n)) != S(L[I]).verschiebung(n):
                    ....:             return False
                    ....:     return True
                    sage: test_L(4, 2)
                    True
                    sage: test_L(6, 2)
                    True
                    sage: test_L(6, 3)
                    True
                    sage: test_L(8, 4)     # long time
                    True
                """
                parent = self.parent()
                C = parent._indices
                return parent.sum_of_terms([(C([i // n for i in I]),
                                            coeff * (-1) ** (sum(I) * (n-1) // n))
                                            for (I, coeff) in self
                                            if all(i % n == 0 for i in I)],
                                           distinct=True)

            def star_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the star involution.

                The star involution is defined as the algebra antihomomorphism
                `NCSF \to NCSF` which, for every positive integer `n`, sends
                the `n`-th complete non-commutative symmetric function `S_n` to
                `S_n`. Denoting by `f^{\ast}` the image of an element
                `f \in NCSF` under this star involution, it can be shown that
                every composition `I` satisfies

                .. MATH::

                    (S^I)^{\ast} = S^{I^r}, \quad
                    (\Lambda^I)^{\ast} = \Lambda^{I^r}, \quad
                    R_I^{\ast} = R_{I^r}, \quad
                    (\Phi^I)^{\ast} = \Phi^{I^r},

                where `I^r` denotes the reversed composition of `I`, and
                standard notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary basis,
                `R` for the ribbon basis, and `\Phi` for that of the power-sums
                of the second kind). The star involution is an involution and a
                coalgebra automorphism of `NCSF`. It is an automorphism of the
                graded vector space `NCSF`. Under the canonical isomorphism
                between the `n`-th graded component of `NCSF` and the descent
                algebra of the symmetric group `S_n` (see
                :meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_descent_algebra`),
                the star involution (restricted to
                the `n`-th graded component) corresponds to the automorphism
                of the descent algebra given by
                `x \mapsto \omega_n x \omega_n`, where `\omega_n` is the
                permutation `(n, n-1, \ldots, 1) \in S_n` (written here in
                one-line notation). If `\pi` denotes the projection from `NCSF`
                to the ring of symmetric functions
                (:meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_symmetric_function`),
                then `\pi(f^{\ast}) = \pi(f)` for every `f \in NCSF`.

                The star involution on `NCSF` is adjoint to the star involution
                on `QSym` by the standard adjunction between `NCSF` and `QSym`.

                The star involution has been denoted by `\rho` in [LMvW13]_,
                section 3.6.
                See [NCSF2]_, section 2.3 for the properties of this map.

                .. SEEALSO::

                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`,
                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: L = NSym.L()
                    sage: L[3,3,2,3].star_involution()
                    L[3, 2, 3, 3]
                    sage: L[6,3,3].star_involution()
                    L[3, 3, 6]
                    sage: (L[1,9,1] - L[8,2] + 2*L[6,4] - 3*L[3] + 4*L[[]]).star_involution()
                    4*L[] + L[1, 9, 1] - L[2, 8] - 3*L[3] + 2*L[4, 6]
                    sage: (L[3,3] - 2*L[2]).star_involution()
                    -2*L[2] + L[3, 3]
                    sage: L([4,1]).star_involution()
                    L[1, 4]

                The implementation at hand is tailored to the elementary basis.
                It is equivalent to the generic implementation via the
                complete basis::

                    sage: S = NSym.S()
                    sage: all( S(L[I].star_involution()) == S(L[I]).star_involution()
                    ....:      for I in Compositions(4) )
                    True
                """
                parent = self.parent()
                dct = {I.reversed(): coeff for (I, coeff) in self}
                return parent._from_dict(dct)

            def psi_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the involution `\psi`.

                The involution `\psi` is defined as the linear map
                `NCSF \to NCSF` which, for every composition `I`, sends the
                complete noncommutative symmetric function `S^I` to the
                elementary noncommutative symmetric function `\Lambda^I`.
                It can be shown that every composition `I` satisfies

                .. MATH::

                    \psi(R_I) = R_{I^c}, \quad \psi(S^I) = \Lambda^I, \quad
                    \psi(\Lambda^I) = S^I, \quad
                    \psi(\Phi^I) = (-1)^{|I| - \ell(I)} \Phi^I

                where `I^c` denotes the complement of the composition `I`, and
                `\ell(I)` denotes the length of `I`, and where standard
                notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary basis,
                `\Phi` for the basis of the power sums of the second kind,
                and `R` for the ribbon basis). The map `\psi` is an involution
                and a graded Hopf algebra automorphism of `NCSF`. If `\pi`
                denotes the projection from `NCSF` to the ring of symmetric functions
                (:meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_symmetric_function`),
                then `\pi(\psi(f)) = \omega(\pi(f))` for every `f \in NCSF`, where
                the `\omega` on the right hand side denotes the omega
                automorphism of `Sym`.

                The involution `\psi` of `NCSF` is adjoint to the involution
                `\psi` of `QSym` by the standard adjunction between `NCSF` and
                `QSym`.

                The involution `\psi` has been denoted by `\psi` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`psi involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = NSym.S()
                    sage: L = NSym.L()
                    sage: L[3,1].psi_involution()
                    L[1, 1, 1, 1] - L[1, 2, 1] - L[2, 1, 1] + L[3, 1]
                    sage: S(L[3,1].psi_involution())
                    S[3, 1]
                    sage: L[[]].psi_involution()
                    L[]
                    sage: L[1,1].psi_involution()
                    L[1, 1]
                    sage: (L[2,1] - 2*L[2]).psi_involution()
                    -2*L[1, 1] + L[1, 1, 1] + 2*L[2] - L[2, 1]

                The implementation at hand is tailored to the elementary basis.
                It is equivalent to the generic implementation via the
                ribbon basis::

                    sage: R = NSym.R()
                    sage: all( R(L[I].psi_involution()) == R(L[I]).psi_involution()
                    ....:      for I in Compositions(3) )
                    True
                    sage: all( R(L[I].psi_involution()) == R(L[I]).psi_involution()
                    ....:      for I in Compositions(4) )
                    True
                """
                parent = self.parent()
                return parent.sum( (-1) ** (I.size() - len(I)) * coeff
                                   * parent.alternating_sum_of_finer_compositions(I)
                                   for I, coeff in
                                   self._monomial_coefficients.items() )

    L = elementary = Elementary

    class Psi(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of non-commutative symmetric functions in the
        Psi basis.

        The Psi basis is defined in Definition 3.4 of [NCSF1]_, where
        it is denoted by `(\Psi^I)_I`. It is a multiplicative basis, and
        is connected to the elementary generators `\Lambda_i` of the ring
        of non-commutative symmetric functions by the following relation:
        Define a non-commutative symmetric function `\Psi_n` for every
        positive integer `n` by the power series identity

        .. MATH::

            \frac{d}{dt} \sigma(t)
            = \sigma(t) \cdot \left( \sum_{k \geq 1} t^{k-1} \Psi_k \right),

        where

        .. MATH::

            \sigma(t) = \left( \sum_{k \geq 0} (-t)^k \Lambda_k \right)^{-1}

        and where `\Lambda_0` denotes `1`. For every composition
        `(i_1, i_2, \ldots, i_k)`, we have
        `\Psi^{(i_1, i_2, \ldots, i_k)}
        = \Psi_{i_1} \Psi_{i_2} \cdots \Psi_{i_k}`.

        The `\Psi`-basis is a basis only when the base ring is a
        `\QQ`-algebra (although the `\Psi^I` can be defined over any base
        ring). The elements of the `\Psi`-basis are known as the
        "power-sum non-commutative symmetric functions of the first kind".
        The generators `\Psi_n` correspond to the Dynkin
        (quasi-)idempotents in the descent algebras of the symmetric
        groups (see [NCSF1]_, 5.2 for details).

        Another (equivalent) definition of `\Psi_n` is

        .. MATH::

            \Psi_n = \sum_{i=0}^{n-1} (-1)^i R_{1^i, n-i},

        where `R` denotes the ribbon basis of `NCSF`, and where `1^i`
        stands for `i` repetitions of the integer `1`.

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: Psi = NCSF.Psi(); Psi
            Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
            sage: Psi.an_element()
            2*Psi[] + 2*Psi[1] + 3*Psi[1, 1]

        Checking the equivalent definition of `\Psi_n`::

            sage: def test_psi(n):
            ....:     NCSF = NonCommutativeSymmetricFunctions(ZZ)
            ....:     R = NCSF.R()
            ....:     Psi = NCSF.Psi()
            ....:     a = R.sum([(-1) ** i * R[[1]*i + [n-i]]
            ....:                for i in range(n)])
            ....:     return a == R(Psi[n])
            sage: test_psi(2)
            True
            sage: test_psi(3)
            True
            sage: test_psi(4)
            True
        """
        def __init__(self, NCSF):
            r"""
            TESTS:

            We include a sanity test to verify the conversion to
            and from the complete basis works the way it should::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi(); Psi
                Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
                sage: all(S(Psi(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(Psi(S(Psi[comp])) == Psi[comp] for comp in Compositions(5))
                True
            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='Psi', bracket=False,
                                             category=NCSF.MultiplicativeBasesOnPrimitiveElements())

        def _from_complete_on_generators(self, n):
            r"""
            Expand a complete generator of non-commutative symmetric
            functions in the Psi basis.

            INPUT:

            - ``n`` -- a positive integer

            OUTPUT:

            - The expansion of the complete generator indexed by ``n`` into the
              Psi basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi._from_complete_on_generators(1)
                Psi[1]
                sage: Psi._from_complete_on_generators(2)
                1/2*Psi[1, 1] + 1/2*Psi[2]
                sage: Psi._from_complete_on_generators(3)
                1/6*Psi[1, 1, 1] + 1/3*Psi[1, 2] + 1/6*Psi[2, 1] + 1/3*Psi[3]
            """
            # Equation (58) of NCSF I article
            one = self.base_ring().one()
            I = self._indices([n])
            return self.sum_of_terms( ( (J, one/coeff_pi(J,I)) for J in Compositions(n) ),
                                      distinct=True )

        def _to_complete_on_generators(self, n):
            r"""
            Expand a `\Psi` basis element of non-commutative symmetric
            functions in the complete basis.

            This formula is given in Proposition 4.5 of [NCSF1]_ which states

            .. MATH::

                \Psi_n = \sum_{J \models n} (-1)^{\ell(J)-1} lp(J,I) S^J.

            The coefficient `lp(J,I)` is given in the function
            :meth:`sage.combinat.ncsf_qsym.combinatorics.coeff_lp`.

            INPUT:

            - ``n`` -- a positive integer

            OUTPUT:

            - The expansion of the `\Psi` function indexed by ``n`` in the
              complete basis.

            TESTS::

                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi._to_complete_on_generators(1)
                S[1]
                sage: Psi._to_complete_on_generators(2)
                -S[1, 1] + 2*S[2]
                sage: Psi._to_complete_on_generators(3)
                S[1, 1, 1] - 2*S[1, 2] - S[2, 1] + 3*S[3]
            """
            minus_one = -self.base_ring().one()
            complete = self.realization_of().complete()
            return complete.sum_of_terms( ((J, minus_one**(len(J)+1)*coeff_lp(J,[n]))
                        for J in Compositions(n)), distinct=True )

        def internal_product_on_basis_by_bracketing(self, I, J):
            r"""
            The internal product of two elements of the Psi basis.

            See :meth:`~sage.combinat.ncsf_qsym.generic_basis_code.GradedModulesWithInternalProduct.ElementMethods.internal_product`
            for a thorough documentation of this operation.

            This is an implementation using [NCSF2]_ Lemma 3.10.
            It is fast when the length of `I` is small, but can get
            very slow otherwise. Therefore it is not being used by
            default for internally multiplying Psi functions.

            INPUT:

            - ``I``, ``J`` -- compositions

            OUTPUT:

            - The internal product of the elements of the Psi basis of
              `NSym` indexed by ``I`` and ``J``, expressed in the Psi
              basis.

            AUTHORS:

            - Travis Scrimshaw, 29 Mar 2014

            EXAMPLES::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: Psi = N.Psi()
                sage: Psi.internal_product_on_basis_by_bracketing([2,2],[1,2,1])
                0
                sage: Psi.internal_product_on_basis_by_bracketing([1,2,1],[2,1,1])
                4*Psi[1, 2, 1]
                sage: Psi.internal_product_on_basis_by_bracketing([2,1,1],[1,2,1])
                4*Psi[2, 1, 1]
                sage: Psi.internal_product_on_basis_by_bracketing([1,2,1], [1,1,1,1])
                0
                sage: Psi.internal_product_on_basis_by_bracketing([3,1], [1,2,1])
                -Psi[1, 2, 1] + Psi[2, 1, 1]
                sage: Psi.internal_product_on_basis_by_bracketing([1,2,1], [3,1])
                0
                sage: Psi.internal_product_on_basis_by_bracketing([2,2],[1,2])
                0
                sage: Psi.internal_product_on_basis_by_bracketing([4], [1,2,1])
                -Psi[1, 1, 2] + 2*Psi[1, 2, 1] - Psi[2, 1, 1]

            TESTS:

            The internal product computed by this method is identical with
            the one obtained by coercion to the complete basis::

                sage: S = N.S()
                sage: def psi_int_test(n):
                ....:     for I in Compositions(n):
                ....:         for J in Compositions(n):
                ....:             a = S(Psi.internal_product_on_basis_by_bracketing(I, J))
                ....:             b = S(Psi[I]).internal_product(S(Psi[J]))
                ....:             if a != b:
                ....:                 return False
                ....:     return True
                sage: all( psi_int_test(i) for i in range(4) )
                True
                sage: psi_int_test(4)   # long time
                True
            """
            # The algorithm used here is described in
            # :meth:`generic_basis_code.GradedModulesWithInternalProduct.ElementMethods.internal_product`.
            if sum(I) != sum(J):
                return self.zero()
            p = len(I)
            q = len(J)
            if p > q:
                return self.zero()
            if p == q:
                Is = sorted(I, reverse=True)
                Js = sorted(J, reverse=True)
                if Is != Js:
                    return 0
                return Partition(Is).centralizer_size() * self[I]

            # If we're still here, we must have p < q.
            def Gamma(K):
                r"""
                Compute `\Gamma_K` for a nonempty composition `K` (which
                can be encoded as a list). See the doc of
                :meth:`~sage.combinat.ncsf_qsym.generic_basis_code.GradedModulesWithInternalProduct.ElementMethods.internal_product`
                for a definition of this.
                """
                k1 = K[0]
                res = k1 * self[k1]
                for k in K[1:]:
                    Psik = self[k]
                    res = res * Psik - Psik * res
                return res

            # Special case when I = [n], there is exactly one ordered set
            #   partition and letting this case through would mean another
            #   case check during the backtracking algorithm
            if p == 1:
                return Gamma(J)

            # We now need to sum over all ordered set partitions
            # `(K_1, K_2, \ldots, K_p)` of `\{ 1, 2, \ldots, q \}`
            # into `p` parts such that
            # each `0 \leq k < p` satisfies `|J_{K_k}| = I_k`.
            # To do so, we will encode such partitions as lists
            # of subsets (which themselves are encoded as lists,
            # in increasing order, with every entry decremented by
            # 1 so as to simplify indexing).
            # We create a variable K which traverses
            # (among other things) these ordered set partitions in
            # lexicographic order (on lists of lists of integers,
            # NOT flattened). It follows a backtracking
            # algorithm; when not backtracking, the last entry
            # of its last part will be "exploring" different
            # possible values.
            K = [[-1]]
            cur_sum = 0
            # base will be the set of elements that are currently
            # not in K (again, all decremented by 1). Here, the
            # last entry of the last part of K does not count as
            # being in K when we are between ordered set partitions.
            base = set(range(q))

            result = self.zero()
            while True:
                # If K is too long or there is nothing more to add:
                # backtrack by removing the last part of K.
                if len(K) > p or not base:
                    # Remove the last part from K.
                    base.union(K.pop()[:-1])
                    # We don't need checks here since p > 0 and all parts
                    #   have size > 0 or we couldn't have added everything to the first
                    part = K[-1]
                    base.add(part[-1])
                    # Similarly, we can just continue on
                else:
                    part = K[-1]
                # part is now the last part of K.

                # Find a part `K_k` such that `|J_{K_k}| = I_k`
                Ik = I[len(K) - 1] # -1 for indexing
                cur_sum = sum(J[j] for j in part[:-1]) # The last entry hasn't been added yet

                while cur_sum != Ik:
                    part[-1] += 1
                    # If we can't add the value (because it is too large
                    #    for the base): backtrack
                    if part[-1] >= q:
                        part.pop()
                        if not part:
                            break
                        base.add(part[-1])
                        cur_sum -= J[part[-1]]
                    elif part[-1] in base and cur_sum + J[part[-1]] <= Ik:
                        cur_sum += J[part[-1]]
                        base.remove(part[-1])
                        if cur_sum < Ik: # Still more work to do
                            part.append(part[-1])

                # If the last part is empty (i.e. we didn't find a part): backtrack
                if not part:
                    K.pop()
                    if not K:
                        break
                    base.add(K[-1][-1])
                    continue

                if not base and len(K) == p:
                    # We've found such a set partition
                    result += self.prod(Gamma(tuple(J[j] for j in S)) for S in K)
                    base.add(part[-1])
                else:
                    # Otherwise create a new part
                    K.append([-1])

            return result

        class Element(CombinatorialFreeModule.Element):

            def verschiebung(self, n):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the `n`-th Verschiebung operator.

                The `n`-th Verschiebung operator `\mathbf{V}_n` is defined
                to be the map from the `\mathbf{k}`-algebra of noncommutative
                symmetric functions to itself that sends the complete function
                `S^I` indexed by a composition `I = (i_1, i_2, \ldots , i_k)`
                to `S^{(i_1/n, i_2/n, \ldots , i_k/n)}` if all of the numbers
                `i_1, i_2, \ldots, i_k` are divisible by `n`, and to `0`
                otherwise. This operator `\mathbf{V}_n` is a Hopf algebra
                endomorphism. For every positive integer `r` with `n \mid r`,
                it satisfies

                .. MATH::

                    \mathbf{V}_n(S_r) = S_{r/n},
                    \quad \mathbf{V}_n(\Lambda_r) = (-1)^{r - r/n} \Lambda_{r/n},
                    \quad \mathbf{V}_n(\Psi_r) = n \Psi_{r/n},
                    \quad \mathbf{V}_n(\Phi_r) = n \Phi_{r/n}

                (where `S_r` denotes the `r`-th complete non-commutative
                symmetric function, `\Lambda_r` denotes the `r`-th elementary
                non-commutative symmetric function, `\Psi_r` denotes the `r`-th
                power-sum non-commutative symmetric function of the first kind,
                and `\Phi_r` denotes the `r`-th power-sum non-commutative
                symmetric function of the second kind). For every positive
                integer `r` with `n \nmid r`, it satisfes

                .. MATH::

                    \mathbf{V}_n(S_r) = \mathbf{V}_n(\Lambda_r)
                    = \mathbf{V}_n(\Psi_r) = \mathbf{V}_n(\Phi_r) = 0.

                The `n`-th Verschiebung operator is also called the `n`-th
                Verschiebung endomorphism.

                It is a lift of the `n`-th Verschiebung operator on the ring
                of symmetric functions
                (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung`)
                to the ring of noncommutative symmetric functions.

                The action of the `n`-th Verschiebung operator can also be
                described on the ribbon Schur functions. Namely, every
                composition `I` of size `n \ell` satisfies

                .. MATH::

                    \mathbf{V}_n ( R_I )
                    = (-1)^{\ell(I) - \ell(J)} \cdot R_{J / n},

                where `J` denotes the meet of the compositions `I` and
                `(\underbrace{n, n, \ldots, n}_{|I|/n \mbox{ times}})`,
                where `\ell(I)` is the length of `I`, and
                where `J / n` denotes the composition obtained by
                dividing every entry of `J` by `n`.
                For a composition `I` of size not divisible by `n`, we
                have `\mathbf{V}_n ( R_I ) = 0`.

                .. SEEALSO::

                    :meth:`verschiebung method of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.verschiebung>`,
                    :meth:`frobenius method of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.frobenius>`,
                    :meth:`verschiebung method of Sym
                    <sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung>`

                INPUT:

                - ``n`` -- a positive integer

                OUTPUT:

                The result of applying the `n`-th Verschiebung operator (on the
                ring of noncommutative symmetric functions) to ``self``.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: Psi = NSym.Psi()
                    sage: Psi([4,2]).verschiebung(2)
                    4*Psi[2, 1]
                    sage: Psi([2,4]).verschiebung(2)
                    4*Psi[1, 2]
                    sage: Psi([6]).verschiebung(2)
                    2*Psi[3]
                    sage: Psi([2,1]).verschiebung(3)
                    0
                    sage: Psi([3]).verschiebung(2)
                    0
                    sage: Psi([]).verschiebung(2)
                    Psi[]
                    sage: Psi([5, 1]).verschiebung(3)
                    0
                    sage: Psi([5, 1]).verschiebung(6)
                    0
                    sage: Psi([5, 1]).verschiebung(2)
                    0
                    sage: Psi([1, 2, 3, 1]).verschiebung(7)
                    0
                    sage: Psi([7]).verschiebung(7)
                    7*Psi[1]
                    sage: Psi([1, 2, 3, 1]).verschiebung(5)
                    0
                    sage: (Psi[1] - Psi[2] + 2*Psi[3]).verschiebung(1)
                    Psi[1] - Psi[2] + 2*Psi[3]

                TESTS:

                The current implementation on the Psi basis gives the
                same results as the default implementation::

                    sage: S = NSym.S()
                    sage: def test_psi(N, n):
                    ....:     for I in Compositions(N):
                    ....:         if S(Psi[I].verschiebung(n)) != S(Psi[I]).verschiebung(n):
                    ....:             return False
                    ....:     return True
                    sage: test_psi(4, 2)
                    True
                    sage: test_psi(6, 2)
                    True
                    sage: test_psi(6, 3)
                    True
                    sage: test_psi(8, 4)     # long time
                    True
                """
                parent = self.parent()
                C = parent._indices
                return parent.sum_of_terms([(C([i // n for i in I]),
                                            coeff * (n ** len(I)))
                                            for (I, coeff) in self
                                            if all(i % n == 0 for i in I)],
                                           distinct=True)

    class Phi(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of non-commutative symmetric functions in the
        Phi basis.

        The Phi basis is defined in Definition 3.4 of [NCSF1]_, where
        it is denoted by `(\Phi^I)_I`. It is a multiplicative basis, and
        is connected to the elementary generators `\Lambda_i` of the ring
        of non-commutative symmetric functions by the following relation:
        Define a non-commutative symmetric function `\Phi_n` for every
        positive integer `n` by the power series identity

        .. MATH::

            \sum_{k\geq 1} t^k \frac{1}{k} \Phi_k
            = -\log \left( \sum_{k \geq 0} (-t)^k \Lambda_k \right),

        with `\Lambda_0` denoting `1`. For every composition
        `(i_1, i_2, \ldots, i_k)`, we have
        `\Phi^{(i_1, i_2, \ldots, i_k)}
        = \Phi_{i_1} \Phi_{i_2} \cdots \Phi_{i_k}`.

        The `\Phi`-basis is well-defined only when the base ring is a
        `\QQ`-algebra. The elements of the `\Phi`-basis are known as the
        "power-sum non-commutative symmetric functions of the second
        kind".

        The generators `\Phi_n` are related to the (first) Eulerian
        idempotents in the descent algebras of the symmetric groups (see
        [NCSF1]_, 5.4 for details).

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: Phi = NCSF.Phi(); Phi
            Non-Commutative Symmetric Functions over the Rational Field in the Phi basis
            sage: Phi.an_element()
            2*Phi[] + 2*Phi[1] + 3*Phi[1, 1]
        """
        def __init__(self, NCSF):
            r"""
            TESTS:

            We include a sanity test to verify the conversion to
            and from the complete basis works the way it should::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Phi = NonCommutativeSymmetricFunctions(QQ).Phi(); Phi
                Non-Commutative Symmetric Functions over the Rational Field in the Phi basis
                sage: all(S(Phi(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(Phi(S(Phi[comp])) == Phi[comp] for comp in Compositions(5))
                True

            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='Phi', bracket=False,
                                             category=NCSF.MultiplicativeBasesOnPrimitiveElements())

        def _from_complete_on_generators(self, n):
            r"""
            Expand a complete basis element of non-commutative symmetric
            functions in the `\Phi` basis.

            INPUT:

            - ``n`` -- a positive integer

            OUTPUT:

            - The expansion of the complete function indexed by ``n`` in the
              `\Phi` basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Phi = NonCommutativeSymmetricFunctions(QQ).Phi()
                sage: Phi._from_complete_on_generators(1)
                Phi[1]
                sage: Phi._from_complete_on_generators(2)
                1/2*Phi[1, 1] + 1/2*Phi[2]
                sage: Phi._from_complete_on_generators(3)
                1/6*Phi[1, 1, 1] + 1/4*Phi[1, 2] + 1/4*Phi[2, 1] + 1/3*Phi[3]
            """
            # Proposition 4.9 of NCSF I article
            one = self.base_ring().one()
            return self.sum_of_terms( ( (J, one / coeff_sp(J,[n])) for J in Compositions(n) ),
                                      distinct=True )

        def _to_complete_on_generators(self, n):
            r"""
            Expand a `\Phi` basis element of non-commutative symmetric
            functions in the complete basis.

            INPUT:

            - ``n`` -- a positive integer

            OUTPUT:

            - The expansion of the `\Phi` function indexed by ``n`` in the
              complete basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Phi = NonCommutativeSymmetricFunctions(QQ).Phi()
                sage: Phi._to_complete_on_generators(1)
                S[1]
                sage: Phi._to_complete_on_generators(2)
                -S[1, 1] + 2*S[2]
                sage: Phi._to_complete_on_generators(3)
                S[1, 1, 1] - 3/2*S[1, 2] - 3/2*S[2, 1] + 3*S[3]
            """
            # Proposition 4.9 of NCSF I article
            minus_one = -self.base_ring().one()
            complete = self.realization_of().complete()
            return complete.sum_of_terms( ( (J, minus_one**(len(J)+1) * n / coeff_ell(J,[n]))
                                            for J in Compositions(n) ),
                                          distinct=True )

        class Element(CombinatorialFreeModule.Element):

            def verschiebung(self, n):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the `n`-th Verschiebung operator.

                The `n`-th Verschiebung operator `\mathbf{V}_n` is defined
                to be the map from the `\mathbf{k}`-algebra of noncommutative
                symmetric functions to itself that sends the complete function
                `S^I` indexed by a composition `I = (i_1, i_2, \ldots , i_k)`
                to `S^{(i_1/n, i_2/n, \ldots , i_k/n)}` if all of the numbers
                `i_1, i_2, \ldots, i_k` are divisible by `n`, and to `0`
                otherwise. This operator `\mathbf{V}_n` is a Hopf algebra
                endomorphism. For every positive integer `r` with `n \mid r`,
                it satisfies

                .. MATH::

                    \mathbf{V}_n(S_r) = S_{r/n},
                    \quad \mathbf{V}_n(\Lambda_r) = (-1)^{r - r/n} \Lambda_{r/n},
                    \quad \mathbf{V}_n(\Psi_r) = n \Psi_{r/n},
                    \quad \mathbf{V}_n(\Phi_r) = n \Phi_{r/n}

                (where `S_r` denotes the `r`-th complete non-commutative
                symmetric function, `\Lambda_r` denotes the `r`-th elementary
                non-commutative symmetric function, `\Psi_r` denotes the `r`-th
                power-sum non-commutative symmetric function of the first kind,
                and `\Phi_r` denotes the `r`-th power-sum non-commutative
                symmetric function of the second kind). For every positive
                integer `r` with `n \nmid r`, it satisfes

                .. MATH::

                    \mathbf{V}_n(S_r) = \mathbf{V}_n(\Lambda_r)
                    = \mathbf{V}_n(\Psi_r) = \mathbf{V}_n(\Phi_r) = 0.

                The `n`-th Verschiebung operator is also called the `n`-th
                Verschiebung endomorphism.

                It is a lift of the `n`-th Verschiebung operator on the ring
                of symmetric functions
                (:meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung`)
                to the ring of noncommutative symmetric functions.

                The action of the `n`-th Verschiebung operator can also be
                described on the ribbon Schur functions. Namely, every
                composition `I` of size `n \ell` satisfies

                .. MATH::

                    \mathbf{V}_n ( R_I )
                    = (-1)^{\ell(I) - \ell(J)} \cdot R_{J / n},

                where `J` denotes the meet of the compositions `I` and
                `(\underbrace{n, n, \ldots, n}_{|I|/n \mbox{ times}})`,
                where `\ell(I)` is the length of `I`, and
                where `J / n` denotes the composition obtained by
                dividing every entry of `J` by `n`.
                For a composition `I` of size not divisible by `n`, we
                have `\mathbf{V}_n ( R_I ) = 0`.

                .. SEEALSO::

                    :meth:`verschiebung method of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.verschiebung>`,
                    :meth:`frobenius method of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.frobenius>`,
                    :meth:`verschiebung method of Sym
                    <sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung>`

                INPUT:

                - ``n`` -- a positive integer

                OUTPUT:

                The result of applying the `n`-th Verschiebung operator (on the
                ring of noncommutative symmetric functions) to ``self``.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(ZZ)
                    sage: Phi = NSym.Phi()
                    sage: Phi([4,2]).verschiebung(2)
                    4*Phi[2, 1]
                    sage: Phi([2,4]).verschiebung(2)
                    4*Phi[1, 2]
                    sage: Phi([6]).verschiebung(2)
                    2*Phi[3]
                    sage: Phi([2,1]).verschiebung(3)
                    0
                    sage: Phi([3]).verschiebung(2)
                    0
                    sage: Phi([]).verschiebung(2)
                    Phi[]
                    sage: Phi([5, 1]).verschiebung(3)
                    0
                    sage: Phi([5, 1]).verschiebung(6)
                    0
                    sage: Phi([5, 1]).verschiebung(2)
                    0
                    sage: Phi([1, 2, 3, 1]).verschiebung(7)
                    0
                    sage: Phi([7]).verschiebung(7)
                    7*Phi[1]
                    sage: Phi([1, 2, 3, 1]).verschiebung(5)
                    0
                    sage: (Phi[1] - Phi[2] + 2*Phi[3]).verschiebung(1)
                    Phi[1] - Phi[2] + 2*Phi[3]

                TESTS:

                The current implementation on the Phi basis gives the
                same results as the default implementation::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = NSym.S()
                    sage: Phi = NSym.Phi()
                    sage: def test_phi(N, n):
                    ....:     for I in Compositions(N):
                    ....:         if S(Phi[I].verschiebung(n)) != S(Phi[I]).verschiebung(n):
                    ....:             return False
                    ....:     return True
                    sage: test_phi(4, 2)
                    True
                    sage: test_phi(6, 2)
                    True
                    sage: test_phi(6, 3)
                    True
                    sage: test_phi(8, 4)     # long time
                    True
                """
                parent = self.parent()
                C = parent._indices
                return parent.sum_of_terms([(C([i // n for i in I]),
                                            coeff * (n ** len(I)))
                                            for (I, coeff) in self
                                            if all(i % n == 0 for i in I)],
                                           distinct=True)

            def star_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the star involution.

                The star involution is defined as the algebra antihomomorphism
                `NCSF \to NCSF` which, for every positive integer `n`, sends
                the `n`-th complete non-commutative symmetric function `S_n` to
                `S_n`. Denoting by `f^{\ast}` the image of an element
                `f \in NCSF` under this star involution, it can be shown that
                every composition `I` satisfies

                .. MATH::

                    (S^I)^{\ast} = S^{I^r}, \quad
                    (\Lambda^I)^{\ast} = \Lambda^{I^r}, \quad
                    R_I^{\ast} = R_{I^r}, \quad
                    (\Phi^I)^{\ast} = \Phi^{I^r},

                where `I^r` denotes the reversed composition of `I`, and
                standard notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary basis,
                `R` for the ribbon basis, and `\Phi` for that of the power-sums
                of the second kind). The star involution is an involution and a
                coalgebra automorphism of `NCSF`. It is an automorphism of the
                graded vector space `NCSF`. Under the canonical isomorphism
                between the `n`-th graded component of `NCSF` and the descent
                algebra of the symmetric group `S_n` (see
                :meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_descent_algebra`),
                the star involution (restricted to
                the `n`-th graded component) corresponds to the automorphism
                of the descent algebra given by
                `x \mapsto \omega_n x \omega_n`, where `\omega_n` is the
                permutation `(n, n-1, \ldots, 1) \in S_n` (written here in
                one-line notation). If `\pi` denotes the projection from `NCSF`
                to the ring of symmetric functions
                (:meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_symmetric_function`),
                then `\pi(f^{\ast}) = \pi(f)` for every `f \in NCSF`.

                The star involution on `NCSF` is adjoint to the star involution
                on `QSym` by the standard adjunction between `NCSF` and `QSym`.

                The star involution has been denoted by `\rho` in [LMvW13]_,
                section 3.6.
                See [NCSF2]_, section 2.3 for the properties of this map.

                .. SEEALSO::

                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`,
                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: Phi = NSym.Phi()
                    sage: Phi[3,1,1,4].star_involution()
                    Phi[4, 1, 1, 3]
                    sage: Phi[4,2,1].star_involution()
                    Phi[1, 2, 4]
                    sage: (Phi[1,4] - Phi[2,3] + 2*Phi[5,4] - 3*Phi[3] + 4*Phi[[]]).star_involution()
                    4*Phi[] - 3*Phi[3] - Phi[3, 2] + Phi[4, 1] + 2*Phi[4, 5]
                    sage: (Phi[3,3] + 3*Phi[1]).star_involution()
                    3*Phi[1] + Phi[3, 3]
                    sage: Phi([2,1]).star_involution()
                    Phi[1, 2]

                The implementation at hand is tailored to the Phi basis.
                It is equivalent to the generic implementation via the
                complete basis::

                    sage: S = NSym.S()
                    sage: all( S(Phi[I].star_involution()) == S(Phi[I]).star_involution()
                    ....:      for I in Compositions(4) )
                    True
                """
                parent = self.parent()
                dct = {I.reversed(): coeff for (I, coeff) in self}
                return parent._from_dict(dct)

            def psi_involution(self):
                r"""
                Return the image of the noncommutative symmetric function
                ``self`` under the involution `\psi`.

                The involution `\psi` is defined as the linear map
                `NCSF \to NCSF` which, for every composition `I`, sends the
                complete noncommutative symmetric function `S^I` to the
                elementary noncommutative symmetric function `\Lambda^I`.
                It can be shown that every composition `I` satisfies

                .. MATH::

                    \psi(R_I) = R_{I^c}, \quad \psi(S^I) = \Lambda^I, \quad
                    \psi(\Lambda^I) = S^I, \quad
                    \psi(\Phi^I) = (-1)^{|I| - \ell(I)} \Phi^I

                where `I^c` denotes the complement of the composition `I`, and
                `\ell(I)` denotes the length of `I`, and where standard
                notations for classical bases of `NCSF` are being used
                (`S` for the complete basis, `\Lambda` for the elementary basis,
                `\Phi` for the basis of the power sums of the second kind,
                and `R` for the ribbon basis). The map `\psi` is an involution
                and a graded Hopf algebra automorphism of `NCSF`. If `\pi`
                denotes the projection from `NCSF` to the ring of symmetric functions
                (:meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_symmetric_function`),
                then `\pi(\psi(f)) = \omega(\pi(f))` for every `f \in NCSF`, where
                the `\omega` on the right hand side denotes the omega
                automorphism of `Sym`.

                The involution `\psi` of `NCSF` is adjoint to the involution
                `\psi` of `QSym` by the standard adjunction between `NCSF` and
                `QSym`.

                The involution `\psi` has been denoted by `\psi` in [LMvW13]_,
                section 3.6.

                .. SEEALSO::

                    :meth:`psi involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`psi involution of QSym
                    <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Bases.ElementMethods.psi_involution>`,
                    :meth:`star involution of NCSF
                    <sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.star_involution>`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: Phi = NSym.Phi()
                    sage: Phi[3,2].psi_involution()
                    -Phi[3, 2]
                    sage: Phi[2,2].psi_involution()
                    Phi[2, 2]
                    sage: Phi[[]].psi_involution()
                    Phi[]
                    sage: (Phi[2,1] - 2*Phi[2]).psi_involution()
                    2*Phi[2] - Phi[2, 1]
                    sage: Phi(0).psi_involution()
                    0

                The implementation at hand is tailored to the Phi basis.
                It is equivalent to the generic implementation via the
                ribbon basis::

                    sage: R = NSym.R()
                    sage: all( R(Phi[I].psi_involution()) == R(Phi[I]).psi_involution()
                    ....:      for I in Compositions(4) )
                    True
                """
                parent = self.parent()
                dct = {I: (-1) ** (I.size() - len(I)) * coeff for (I, coeff) in self}
                return parent._from_dict(dct)

    class Monomial(CombinatorialFreeModule, BindableClass):
        r"""
        The monomial basis defined in Tevlin's paper [Tev2007]_.

        The monomial basis is well-defined only when the base ring is a
        `\QQ`-algebra. It is the basis denoted by `(M^I)_I` in [Tev2007]_.

        TESTS::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: nM = NCSF.monomial(); nM
            Non-Commutative Symmetric Functions over the Rational Field in the Monomial basis
            sage: nM([1,1])*nM([2])
            3*nM[1, 1, 2] + nM[1, 3] + nM[2, 2]
            sage: R = NCSF.ribbon()
            sage: nM(R[1,3,1])
            11*nM[1, 1, 1, 1, 1] + 8*nM[1, 1, 2, 1] + 8*nM[1, 2, 1, 1] + 5*nM[1, 3, 1] + 8*nM[2, 1, 1, 1] + 5*nM[2, 2, 1] + 5*nM[3, 1, 1] + 2*nM[4, 1]
        """

        def __init__(self, NCSF):
            r"""
            TESTS:

            We include a sanity test to verify the conversion to
            and from the complete basis works the way it should::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: nM = NonCommutativeSymmetricFunctions(QQ).Monomial(); nM
                Non-Commutative Symmetric Functions over the Rational Field in the Monomial basis
                sage: all(S(nM(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(nM(S(nM[comp])) == nM[comp] for comp in Compositions(5))
                True

            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='nM', bracket=False,
                                             category=NCSF.Bases())
            category = self.category()

            NCSF = NonCommutativeSymmetricFunctions(self.base_ring())
            S = NCSF.complete()
            Psi = NCSF.Psi()
            to_S = self.module_morphism(
                    on_basis = self._to_complete_on_basis,
                    codomain = S,
                    category = category)
            to_S.register_as_coercion()

            from_psi = Psi.module_morphism(
                        on_basis = self._from_psi_on_basis,
                        codomain = self,
                        category = category)
            from_psi.register_as_coercion()

        def _to_complete_on_basis(self, I):
            r"""
            Expand a Monomial basis element of non-commutative symmetric functions
            in the complete basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The expansion of the Monomial function indexed by ``I`` in
              the complete basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: nM = NonCommutativeSymmetricFunctions(QQ).nM()
                sage: nM._to_complete_on_basis(Composition([1,1,1]))
                S[1, 1, 1] - S[1, 2] - S[2, 1] + S[3]
                sage: nM._to_complete_on_basis(Composition([1,2]))
                -S[1, 1, 1] + 2*S[1, 2] + 1/2*S[2, 1] - 3/2*S[3]
                sage: nM._to_complete_on_basis(Composition([2,1]))
                -S[1, 1, 1] + S[1, 2] + 3/2*S[2, 1] - 3/2*S[3]
                sage: nM._to_complete_on_basis(Composition([3]))
                S[1, 1, 1] - 2*S[1, 2] - S[2, 1] + 3*S[3]
            """
            S = NonCommutativeSymmetricFunctions(self.base_ring()).S()
            return S.sum_of_terms( ( (K, m_to_s_stat(self.base_ring(),I,K))
                                     for K in Compositions(sum(I)) ),
                                   distinct=True )
            # Note: sum(I) works both if I is a list and if I is a composition
            # (although the latter case doesn't work in IPython, cf.
            # trac #15163).

        def _from_psi_on_basis(self, I):
            r"""
            Expand a Psi basis element of non-commutative symmetric functions
            in the Monomial basis.

            INPUT:

            - ``self`` - the Monomial basis of non-commutative symmetric functions
            - ``I`` - a composition

            OUTPUT:

            - The expansion of the Psi function indexed by ``I`` into the Monomial
              basis.

            TESTS::

                sage: nM=NonCommutativeSymmetricFunctions(QQ).nM()
                sage: nM._from_psi_on_basis(Composition([3]))
                nM[3]
                sage: nM._from_psi_on_basis(Composition([1,2]))
                2*nM[1, 2] + nM[3]
                sage: nM._from_psi_on_basis(Composition([2,1]))
                2*nM[2, 1] + nM[3]
                sage: nM._from_psi_on_basis(Composition([1,1,1]))
                6*nM[1, 1, 1] + 2*nM[1, 2] + 4*nM[2, 1] + nM[3]
            """
            M = NonCommutativeSymmetricFunctions(self.base_ring()).nM()
            sum_of_elements = M.zero()
            for J in Compositions(I.size()):
                if I.is_finer(J):
                    len_of_J = len(J)
                    p = [0] + self._indices(I).refinement_splitting_lengths(J).partial_sums()
                    sum_of_elements += prod( (len_of_J - k)**(p[k+1]-p[k]) for k in range(len_of_J) ) * M(J)
            return sum_of_elements

    nM = monomial = Monomial

    class Immaculate(CombinatorialFreeModule, BindableClass):
        r"""
        The immaculate basis of the non-commutative symmetric
        functions.

        The immaculate basis first appears in Berg, Bergeron,
        Saliola, Serrano and Zabrocki's [BBSSZ2012]_. It can be
        described as the family `(\mathfrak{S}_{\alpha})`, where
        `\alpha` runs over all compositions, and
        `\mathfrak{S}_{\alpha}` denotes the immaculate function
        corresponding to `\alpha` (see
        :meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ParentMethods.immaculate_function`).

        If `\alpha` is a composition `(\alpha_1, \alpha_2, \ldots,
        \alpha_m)`, then

        .. MATH::

            \mathfrak{S}_{\alpha}
            = \sum_{\sigma \in S_m} (-1)^{\sigma}
            S_{\alpha_1 + \sigma(1) - 1} S_{\alpha_2 + \sigma(2) - 2}
            \cdots S_{\alpha_m + \sigma(m) - m}.

        .. WARNING::

            This *basis* contains only the immaculate functions
            indexed by compositions (i.e., finite sequences of
            positive integers). To obtain the remaining immaculate
            functions (sensu lato), use the
            :meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ParentMethods.immaculate_function`
            method. Calling the immaculate *basis* with a list
            which is not a composition will currently return
            garbage!

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: I = NCSF.I()
            sage: I([1,3,2])*I([1])
            I[1, 3, 2, 1] + I[1, 3, 3] + I[1, 4, 2] + I[2, 3, 2]
            sage: I([1])*I([1,3,2])
            I[1, 1, 3, 2] - I[2, 2, 1, 2] - I[2, 2, 2, 1] - I[2, 2, 3] - I[3, 2, 2]
            sage: I([1,3])*I([1,1])
            I[1, 3, 1, 1] + I[1, 4, 1] + I[2, 3, 1] + I[2, 4]
            sage: I([3,1])*I([2,1])
            I[3, 1, 2, 1] + I[3, 2, 1, 1] + I[3, 2, 2] + I[3, 3, 1] + I[4, 1, 1, 1] + I[4, 1, 2] + 2*I[4, 2, 1] + I[4, 3] + I[5, 1, 1] + I[5, 2]
            sage: R = NCSF.ribbon()
            sage: I(R[1,3,1])
            I[1, 3, 1] + I[2, 2, 1] + I[2, 3] + I[3, 1, 1] + I[3, 2]
            sage: R(I(R([2,1,3])))
            R[2, 1, 3]
        """

        def __init__(self, NCSF):
            r"""
            TESTS:

            We include a sanity test to verify the conversion to
            and from the complete basis works the way it should::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: I = NonCommutativeSymmetricFunctions(QQ).Immaculate(); I
                Non-Commutative Symmetric Functions over the Rational Field in the Immaculate basis
                sage: all(S(I(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(I(S(I[comp])) == I[comp] for comp in Compositions(5))
                True

            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='I', bracket=False,
                                             category=NCSF.Bases())
            category = self.category()
            S = self.realization_of().complete()
            to_S = self.module_morphism(
                    on_basis = self._to_complete_on_basis,
                    codomain = S,
                    category = category)
            to_S.register_as_coercion()

            from_S = S.module_morphism(
                        on_basis = self._from_complete_on_basis,
                        codomain = self,
                        category = category)
            from_S.register_as_coercion()

        def _realization_name(self):
            r"""
            TESTS::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: I = N.I()
                sage: I._realization_name()
                'Immaculate'
            """
            return "Immaculate"

        def _H(self, alpha):
            r"""
            Return the complete basis element indexed by a list ``alpha`` if
            the list happens to be a composition (possibly with `0`s
            interspersed). Otherwise, return `0`.

            INPUT:

            - ``alpha`` -- a list

            OUTPUT:

            - the complete basis element indexed by ``alpha`` (with any
              zeroes removed) if all entries of ``alpha`` are nonnegative;
              otherwise, `0`

            EXAMPLES::

                sage: I=NonCommutativeSymmetricFunctions(QQ).I()
                sage: I._H([2,0,1])
                S[2, 1]
                sage: I._H([2,0,1,-1])
                0
                sage: I._H([1,0,2])
                S[1, 2]
            """
            S = NonCommutativeSymmetricFunctions(self.base_ring()).complete()
            if any( d < 0 for d in alpha ):
                return S.zero()
            return S( [ d for d in alpha if d > 0 ] )

        @cached_method
        def _to_complete_on_basis(self, alpha):
            r"""
            Return the expansion of an Immaculate basis element in the
            complete basis.

            INPUT:

            - ``alpha`` -- a composition

            OUTPUT:

            - The expansion in the complete basis of the basis element
              of the Immaculate basis ``self`` indexed by the
              composition ``alpha``.

            EXAMPLES::

                sage: I=NonCommutativeSymmetricFunctions(QQ).I()
                sage: I._to_complete_on_basis(Composition([]))
                S[]
                sage: I._to_complete_on_basis(Composition([2,1,3]))
                S[2, 1, 3] - S[2, 2, 2] + S[3, 2, 1] - S[3, 3] - S[4, 1, 1] + S[4, 2]
            """
            alpha_list = alpha._list
            if not alpha_list:
                return self._H([])
            if alpha_list == [1]:
                return self._H([1])
            la = len(alpha_list)
            S = NonCommutativeSymmetricFunctions(self.base_ring()).complete()
            return S.sum( sigma.signature()*self._H( [alpha_list[i]+sigma[i]-(i+1) for i in range(la)] )
                          for sigma in Permutations(la) )

        @cached_method
        def _from_complete_on_basis(self, comp_content):
            r"""
            Return the expansion of a complete basis element in the
            Immaculate basis.

            INPUT:

            - ``comp_content`` -- a composition

            OUTPUT:

            - The expansion in the Immaculate basis of the basis element
              of the complete basis indexed by the composition
              ``comp_content``.

            EXAMPLES::

                sage: I=NonCommutativeSymmetricFunctions(QQ).I()
                sage: I._from_complete_on_basis(Composition([]))
                I[]
                sage: I._from_complete_on_basis(Composition([2,1,3]))
                I[2, 1, 3] + I[2, 2, 2] + I[2, 3, 1] + I[2, 4] + I[3, 1, 2] + I[3, 2, 1] + 2*I[3, 3] + I[4, 1, 1] + 2*I[4, 2] + 2*I[5, 1] + I[6]
            """
            I = NonCommutativeSymmetricFunctions(self.base_ring()).I()
            if not comp_content._list:
                return I([])
            return I.sum_of_terms( ( (comp_shape, number_of_fCT(comp_content, comp_shape))
                                     for comp_shape in Compositions(sum(comp_content)) ),
                                   distinct=True )

        def dual(self):
            r"""
            Return the dual basis to the Immaculate basis of NCSF.

            The basis returned is the dualImmaculate basis of QSym.

            OUTPUT:

            - The dualImmaculate basis of the quasi-symmetric functions.

            EXAMPLES::

                sage: I=NonCommutativeSymmetricFunctions(QQ).Immaculate()
                sage: I.dual()
                Quasisymmetric functions over the Rational Field in the dualImmaculate
                basis
            """
            return self.realization_of().dual().dualImmaculate()

        class Element(CombinatorialFreeModule.Element):
            """
            An element in the Immaculate basis.
            """
            def bernstein_creation_operator(self, n):
                r"""
                Return the image of ``self`` under the `n`-th Bernstein
                creation operator.

                Let `n` be an integer. The `n`-th Bernstein creation
                operator `\mathbb{B}_n` is defined as the endomorphism of
                the space `NSym` of noncommutative symmetric functions
                given by

                .. MATH::

                    \mathbb{B}_n I_{(\alpha_1, \alpha_2, \ldots, \alpha_m)}
                    = I_{(n, \alpha_1, \alpha_2, \ldots, \alpha_m)},

                where `I_{(\alpha_1, \alpha_2, \ldots, \alpha_m)}` is the
                immaculate function associated to the `m`-tuple
                `(\alpha_1, \alpha_2, \ldots, \alpha_m) \in \ZZ^m`.

                This has been introduced in [BBSSZ2012]_, section 3.1, in
                analogy to the Bernstein creation operators on the
                symmetric functions.

                For more information on the `n`-th Bernstein creation
                operator, see
                :meth:`~NonCommutativeSymmetricFunctions.Bases.ElementMethods.bernstein_creation_operator`.

                EXAMPLES::

                    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
                    sage: I = NSym.I()
                    sage: b = I[1,3,2,1]
                    sage: b.bernstein_creation_operator(3)
                    I[3, 1, 3, 2, 1]
                    sage: b.bernstein_creation_operator(5)
                    I[5, 1, 3, 2, 1]
                    sage: elt = b + 3*I[4,1,2]
                    sage: elt.bernstein_creation_operator(1)
                    I[1, 1, 3, 2, 1] + 3*I[1, 4, 1, 2]

                We check that this agrees with the definition on the
                Complete basis::

                    sage: S = NSym.S()
                    sage: S(elt).bernstein_creation_operator(1) == S(elt.bernstein_creation_operator(1))
                    True

                Check on non-positive values of `n`::

                    sage: I[2,2,2].bernstein_creation_operator(-1)
                    I[1, 1, 1, 2] + I[1, 1, 2, 1] + I[1, 2, 1, 1] - I[1, 2, 2]
                    sage: I[2,3,2].bernstein_creation_operator(0)
                    -I[1, 1, 3, 2] - I[1, 2, 2, 2] - I[1, 2, 3, 1] + I[2, 3, 2]
                """
                if n <= 0:
                    return super(NonCommutativeSymmetricFunctions.Immaculate.Element, self).bernstein_creation_operator(n)

                C = Compositions()
                P = self.parent()
                return P.sum_of_terms( (C([n] + list(m)), c) for m,c in self )

    I = Immaculate

    class dualQuasisymmetric_Schur(CombinatorialFreeModule, BindableClass):
        r"""
        The basis of NCSF dual to the Quasisymmetric-Schur basis of QSym.

        The
        :class:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Quasisymmetric_Schur`
        functions are defined in [QSCHUR]_ (see also
        Definition 5.1.1 of [LMvW13]_).  The dual basis in the algebra
        of non-commutative symmetric functions is defined by the following
        formula:

        .. MATH::

            R_\alpha = \sum_{T} dQS_{shape(T)},

        where the sum is over all standard composition tableaux with
        descent composition equal to `\alpha`.
        The
        :class:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Quasisymmetric_Schur`
        basis `QS_\alpha` has the property that

        .. MATH::

            s_\lambda = \sum_{sort(\alpha) = \lambda} QS_\alpha.

        As a consequence the commutative image of a dual
        Quasisymmetric-Schur element in the algebra of symmetric functions
        (the map defined in the method
        :meth:`~sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Bases.ElementMethods.to_symmetric_function`)
        is equal to the Schur function indexed by the decreasing sort of the
        indexing composition.

        .. SEEALSO::

            :class:`~sage.combinat.composition_tableau.CompositionTableaux`,
            :class:`~sage.combinat.composition_tableau.CompositionTableau`.

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: dQS = NCSF.dQS()
            sage: dQS([1,3,2])*dQS([1])
            dQS[1, 2, 4] + dQS[1, 3, 2, 1] + dQS[1, 3, 3] + dQS[3, 2, 2]
            sage: dQS([1])*dQS([1,3,2])
            dQS[1, 1, 3, 2] + dQS[1, 3, 3] + dQS[1, 4, 2] + dQS[2, 3, 2]
            sage: dQS([1,3])*dQS([1,1])
            dQS[1, 3, 1, 1] + dQS[1, 4, 1] + dQS[3, 2, 1] + dQS[4, 2]
            sage: dQS([3,1])*dQS([2,1])
            dQS[1, 1, 4, 1] + dQS[1, 4, 2] + dQS[1, 5, 1] + dQS[2, 4, 1] + dQS[3, 1,
            2, 1] + dQS[3, 2, 2] + dQS[3, 3, 1] + dQS[4, 3] + dQS[5, 2]
            sage: dQS([1,1]).coproduct()
            dQS[] # dQS[1, 1] + dQS[1] # dQS[1] + dQS[1, 1] # dQS[]
            sage: dQS([3,3]).coproduct().monomial_coefficients()[(Composition([1,2]),Composition([1,2]))]
            -1
            sage: S = NCSF.complete()
            sage: dQS(S[1,3,1])
            dQS[1, 3, 1] + dQS[1, 4] + dQS[3, 2] + dQS[4, 1] + dQS[5]
            sage: S(dQS[1,3,1])
            S[1, 3, 1] - S[3, 2] - S[4, 1] + S[5]
            sage: s = SymmetricFunctions(QQ).s()
            sage: s(S(dQS([2,1,3])).to_symmetric_function())
            s[3, 2, 1]
        """

        def __init__(self, NCSF):
            r"""
            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: S = NCSF.complete()
                sage: dQS = NCSF.dualQuasisymmetric_Schur()
                sage: dQS(S(dQS.an_element())) == dQS.an_element()
                True
                sage: S(dQS(S.an_element())) == S.an_element()
                True
                sage: TestSuite(dQS).run() # long time
            """
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='dQS', bracket=False,
                                             category=NCSF.Bases())
            category = self.category()
            self._S = self.realization_of().complete()
            to_S = self.module_morphism(
                    on_basis = self._to_complete_on_basis,
                    codomain = self._S,
                    category = category)
            to_S.register_as_coercion()

            from_S = self._S.module_morphism(
                        on_basis = self._from_complete_on_basis,
                        codomain = self,
                        category = category)
            from_S.register_as_coercion()

        def _realization_name(self):
            r"""
            TESTS::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: dQS = NCSF.dQS()
                sage: dQS._realization_name()
                'dual Quasisymmetric-Schur'
            """
            return "dual Quasisymmetric-Schur"

        @cached_method
        def _to_complete_transition_matrix(self, n):
            r"""
            A matrix representing the transition coefficients to
            the complete basis along with the ordering.

            INPUT:

            - ``n`` -- an integer

            OUTPUT:

            - a pair of a square matrix and the ordered list of compositions

            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: dQS = NCSF.dQS()
                sage: dQS._to_complete_transition_matrix(4)[0]
                [ 1  0  0  0  0  0  0  0]
                [-1  1  0  0  0  0  0  0]
                [-1  0  1  0  0  0  0  0]
                [ 0  0 -1  1  0  0  0  0]
                [ 1 -1  0 -1  1  0  0  0]
                [ 1 -1  0 -1  0  1  0  0]
                [ 1  0 -1 -1  0  0  1  0]
                [-1  1  1  1 -1 -1 -1  1]
            """
            if n == 0:
                return (matrix([[]]), [])
            CO = compositions_order(n)
            # ZZ is faster than over QQ for inverting a matrix
            from sage.rings.integer_ring import ZZ
            MS = MatrixSpace(ZZ, len(CO))
            return (MS([[number_of_SSRCT(al,be) for be in CO] for al in CO]).inverse(),
                    CO)

        @cached_method
        def _to_complete_on_basis(self, comp):
            r"""
            The expansion of the dual Quasisymmetric-Schur basis element
            indexed by ``comp`` in the complete basis.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT:

            - a quasi-symmetric function in the complete basis

            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: dQS = NCSF.dQS()
                sage: dQS._to_complete_on_basis(Composition([1,3,1]))
                S[1, 3, 1] - S[3, 2] - S[4, 1] + S[5]
            """
            if not comp._list:
                return self.one()
            T, comps = self._to_complete_transition_matrix(comp.size())
            i = comps.index(comp)
            return self._S._from_dict({c: T[i,j] for j,c in enumerate(comps)
                                       if T[i,j] != 0},
                                      remove_zeros=False)

        @cached_method
        def _from_complete_on_basis(self, comp_content):
            r"""
            Return the expansion of a complete basis element in the
            dual Quasisymmetric-Schur basis.

            INPUT:

            - ``comp_content`` -- a composition

            OUTPUT:

            - the expansion in the dual Quasisymmetric-Schur basis of
              the basis element of the complete basis indexed by the
              composition ``comp_content``

            EXAMPLES::

                sage: dQS=NonCommutativeSymmetricFunctions(QQ).dQS()
                sage: dQS._from_complete_on_basis(Composition([]))
                dQS[]
                sage: dQS._from_complete_on_basis(Composition([2,1,1]))
                dQS[1, 3] + dQS[2, 1, 1] + dQS[2, 2] + dQS[3, 1] + dQS[4]
            """
            if not comp_content._list:
                return self([])
            return self.sum_of_terms( ( (comp_shape, number_of_SSRCT(comp_content, comp_shape))
                                     for comp_shape in Compositions(sum(comp_content)) ),
                                   distinct=True )

        def dual(self):
            r"""
            The dual basis to the dual Quasisymmetric-Schur basis of NCSF.

            The basis returned is the
            :class:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Quasisymmetric_Schur`
            basis of QSym.

            OUTPUT:

            - the Quasisymmetric-Schur basis of the quasi-symmetric functions

            EXAMPLES::

                sage: dQS=NonCommutativeSymmetricFunctions(QQ).dualQuasisymmetric_Schur()
                sage: dQS.dual()
                Quasisymmetric functions over the Rational Field in the Quasisymmetric
                Schur basis
                sage: dQS.duality_pairing_matrix(dQS.dual(),3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]
            """
            return self.realization_of().dual().Quasisymmetric_Schur()

        def to_symmetric_function_on_basis(self, I):
            r"""
            The commutative image of a dual quasi-symmetric Schur element

            The commutative image of a basis element is obtained by sorting
            the indexing composition of the basis element.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The commutative image of the dual quasi-Schur basis element
              indexed by ``I``. The result is the Schur symmetric function
              indexed by the partition obtained by sorting ``I``.

            EXAMPLES::

                sage: dQS=NonCommutativeSymmetricFunctions(QQ).dQS()
                sage: dQS.to_symmetric_function_on_basis([2,1,3])
                s[3, 2, 1]
                sage: dQS.to_symmetric_function_on_basis([])
                s[]
            """
            s = SymmetricFunctions(self.base_ring()).s()
            return s[Partition(sorted(I,reverse=True))]

    dQS = dualQuasisymmetric_Schur

    class dualYoungQuasisymmetric_Schur(CombinatorialFreeModule, BindableClass):
        r"""
        The basis of NCSF dual to the Young Quasisymmetric-Schur basis of QSym.

        The
        :class:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.YoungQuasisymmetric_Schur`
        functions are given in Definition 5.2.1 of [LMvW13]_.  The dual basis
        in the algebra of non-commutative symmetric functions are related by
        an involution reversing the indexing composition of the complete
        expansion of a quasi-Schur basis element.  This basis has many of the
        same properties as the Quasisymmetric Schur basis and is related to
        that basis by an algebraic transformation.

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: dYQS = NCSF.dYQS()
            sage: dYQS([1,3,2])*dYQS([1])
            dYQS[1, 3, 2, 1] + dYQS[1, 3, 3] + dYQS[1, 4, 2] + dYQS[2, 3, 2]
            sage: dYQS([1])*dYQS([1,3,2])
            dYQS[1, 1, 3, 2] + dYQS[2, 3, 2] + dYQS[3, 3, 1] + dYQS[4, 1, 2]
            sage: dYQS([1,3])*dYQS([1,1])
            dYQS[1, 3, 1, 1] + dYQS[1, 4, 1] + dYQS[2, 3, 1] + dYQS[2, 4]
            sage: dYQS([3,1])*dYQS([2,1])
            dYQS[3, 1, 2, 1] + dYQS[3, 2, 2] + dYQS[3, 3, 1] + dYQS[4, 1, 1, 1]
             + dYQS[4, 1, 2] + dYQS[4, 2, 1] + dYQS[4, 3] + dYQS[5, 1, 1]
             + dYQS[5, 2]
            sage: dYQS([1,1]).coproduct()
            dYQS[] # dYQS[1, 1] + dYQS[1] # dYQS[1] + dYQS[1, 1] # dYQS[]
            sage: dYQS([3,3]).coproduct().monomial_coefficients()[(Composition([1,2]),Composition([2,1]))]
            1
            sage: S = NCSF.complete()
            sage: dYQS(S[1,3,1])
            dYQS[1, 3, 1] + dYQS[1, 4] + dYQS[2, 3] + dYQS[4, 1] + dYQS[5]
            sage: S(dYQS[1,3,1])
            S[1, 3, 1] - S[1, 4] - S[2, 3] + S[5]
            sage: s = SymmetricFunctions(QQ).s()
            sage: s(S(dYQS([2,1,3])).to_symmetric_function())
            s[3, 2, 1]
        """

        def __init__(self, NCSF):
            r"""
            Initialize ``self``.

            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: S = NCSF.complete()
                sage: dYQS = NCSF.dualYoungQuasisymmetric_Schur()
                sage: dYQS(S(dYQS.an_element())) == dYQS.an_element()
                True
                sage: S(dYQS(S.an_element())) == S.an_element()
                True
                sage: TestSuite(dYQS).run() # long time
            """
            category = NCSF.Bases()
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                             prefix='dYQS', bracket=False,
                                             category=category)

            self._S = NCSF.complete()
            self._dQS = NCSF.dualQuasisymmetric_Schur()
            self.module_morphism(on_basis=self._to_complete_on_basis,
                                 codomain=self._S,
                                 category=category).register_as_coercion()

            self._S.module_morphism(on_basis=self._from_complete_on_basis,
                                    codomain=self,
                                    category=category).register_as_coercion()

        def _realization_name(self):
            r"""
            TESTS::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: dYQS = NCSF.dYQS()
                sage: dYQS._realization_name()
                'dual Young Quasisymmetric-Schur'
            """
            return "dual Young Quasisymmetric-Schur"

        def _to_complete_on_basis(self, comp):
            r"""
            The expansion of the dual Quasisymmetric-Schur basis element
            indexed by ``comp`` in the complete basis.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT:

            - a quasi-symmetric function in the complete basis

            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: dYQS = NCSF.dYQS()
                sage: dYQS._to_complete_on_basis(Composition([1,3,1]))
                S[1, 3, 1] - S[1, 4] - S[2, 3] + S[5]
            """
            elt = self._dQS._to_complete_on_basis(comp.reversed())
            return self._S._from_dict({al.reversed(): c for al, c in elt},
                                      coerce=False, remove_zeros=False)

        def _from_complete_on_basis(self, comp):
            r"""
            Return the expansion of a complete basis element in the
            dual Young Quasisymmetric-Schur basis.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT:

            - the expansion in the dual Young Quasisymmetric-Schur basis of
              the basis element in the complete basis indexed by the
              composition ``comp``

            EXAMPLES::

                sage: dYQS=NonCommutativeSymmetricFunctions(QQ).dYQS()
                sage: dYQS._from_complete_on_basis(Composition([]))
                dYQS[]
                sage: dYQS._from_complete_on_basis(Composition([2,1,1]))
                dYQS[2, 1, 1] + dYQS[2, 2] + 2*dYQS[3, 1] + dYQS[4]
            """
            elt = self._dQS._from_complete_on_basis(comp.reversed())
            return self._from_dict({al.reversed(): c for al, c in elt},
                                   coerce=False, remove_zeros=False)

        def dual(self):
            r"""
            The dual basis to the dual Quasisymmetric-Schur basis of NCSF.

            The basis returned is the
            :class:`~sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions.Quasisymmetric_Schur`
            basis of QSym.

            OUTPUT:

            - the Young Quasisymmetric-Schur basis of quasi-symmetric functions

            EXAMPLES::

                sage: dYQS=NonCommutativeSymmetricFunctions(QQ).dualYoungQuasisymmetric_Schur()
                sage: dYQS.dual()
                Quasisymmetric functions over the Rational Field in the Young
                Quasisymmetric Schur basis
                sage: dYQS.duality_pairing_matrix(dYQS.dual(),3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]
            """
            return self.realization_of().dual().Young_Quasisymmetric_Schur()

        def to_symmetric_function_on_basis(self, I):
            r"""
            The commutative image of a dual Young quasi-symmetric
            Schur element.

            The commutative image of a basis element is obtained by sorting
            the indexing composition of the basis element.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The commutative image of the dual Young quasi-Schur basis element
              indexed by ``I``. The result is the Schur symmetric function
              indexed by the partition obtained by sorting ``I``.

            EXAMPLES::

                sage: dYQS=NonCommutativeSymmetricFunctions(QQ).dYQS()
                sage: dYQS.to_symmetric_function_on_basis([2,1,3])
                s[3, 2, 1]
                sage: dYQS.to_symmetric_function_on_basis([])
                s[]
            """
            s = SymmetricFunctions(self.base_ring()).s()
            return s[Partition(sorted(I,reverse=True))]

    dYQS = dualYoungQuasisymmetric_Schur

    class Zassenhaus_left(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of non-commutative symmetric functions in the
        left Zassenhaus basis.

        This basis is the left-version of the basis defined in Section 2.5.1
        of [HLNT09]_.
        It is multiplicative, with `Z_n` defined as the element of `NCSF_n`
        satisfying the equation

        .. MATH::

            \sigma_1 = \cdots exp(Z_n) \cdots exp(Z_2) exp(Z_1),

        where

        .. MATH::

            \sigma_1 = \sum_{n \geq 0} S_n .

        It can be recursively computed by the formula

        .. MATH::

            S_n = \sum_{\lambda \vdash n}
            \frac{1}{m_1(\lambda)! m_2(\lambda)! m_3(\lambda)! \cdots}
            Z_{\lambda_1} Z_{\lambda_2} Z_{\lambda_3} \cdots

        for all `n \geq 0`.
        """
        def __init__(self, NCSF):
            r"""
            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: ZL = NCSF.Zassenhaus_left(); ZL
                Non-Commutative Symmetric Functions over the Rational Field
                 in the Zassenhaus_left basis
                sage: TestSuite(ZL).run()

            TESTS:

            Test coproduct and antipode on the multiplicative identity::

                sage: ZL.one()
                ZL[]
                sage: ZL.one().coproduct()
                ZL[] # ZL[]
                sage: ZL.one().antipode()
                ZL[]

            We include some sanity tests to verify that conversions between
            this basis and other bases work the way they should::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: ZL = NCSF.ZL()
                sage: R = NCSF.ribbon()
                sage: S = NCSF.complete()
                sage: ZL(S[3])
                1/6*ZL[1, 1, 1] + ZL[2, 1] + ZL[3]
                sage: all(S(ZL(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(ZL(S(ZL[comp])) == ZL[comp] for comp in Compositions(5))
                True
                sage: all(R(ZL(R[comp])) == R[comp] for comp in Compositions(5))
                True
                sage: all(ZL(R(ZL[comp])) == ZL[comp] for comp in Compositions(5))
                True
            """
            cat = NCSF.MultiplicativeBasesOnPrimitiveElements()
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                                prefix='ZL', bracket=False,
                                                category=cat)

            # Register coercions
            S = self.realization_of().S()

            to_complete = self.algebra_morphism(self._to_complete_on_generator,
                                                codomain=S)
            to_complete.register_as_coercion()

            from_complete = S.module_morphism(on_basis=self._from_complete_on_basis,
                                              codomain=self)
            from_complete.register_as_coercion()

        def _to_complete_on_generator(self, n):
            r"""
            Expand a (left) Zassenhaus generator of non-commutative symmetric
            functions in the complete basis.

            INPUT:

            - ``n`` -- a positive integer

            OUTPUT:

            The expansion of the (left) Zassenhaus generator indexed by ``n``
            into the complete basis.

            TESTS::

                sage: ZL = NonCommutativeSymmetricFunctions(QQ).Zassenhaus_left()
                sage: ZL._to_complete_on_generator(1)
                S[1]
                sage: ZL._to_complete_on_generator(2)
                -1/2*S[1, 1] + S[2]
            """
            S = self.realization_of().S()
            if n <= 1:
                return S[n]

            from sage.combinat.partitions import ZS1_iterator
            from sage.rings.integer_ring import ZZ
            it = ZS1_iterator(n)
            next(it) # Skip the unique length 1 partition
            res = S[n]
            for p in it:
                d = {}
                for part in p:
                    d[part] = d.get(part, 0) + 1
                coeff = ZZ(prod(factorial(d[l]) for l in d))
                res = res - prod(self._to_complete_on_generator(i) for i in p) / coeff
            return res

        @cached_method
        def _complete_to_zassenhaus_transition_matrix_inverse(self, n):
            r"""
            The change of basis matrix from the S basis to the ZL basis.

            EXAMPLES::

                sage: ZL = NonCommutativeSymmetricFunctions(QQ).Zassenhaus_left()
                sage: ZL._complete_to_zassenhaus_transition_matrix_inverse(3)
                [  1   0   0   0]
                [1/2   1   0   0]
                [1/2   0   1   0]
                [1/6   0   1   1]
            """
            from sage.matrix.constructor import matrix
            S = self.realization_of().S()
            m = []
            for I in Compositions(n):
                x = S(self.basis()[I])
                m.append([x.coefficient(J) for J in Compositions(n)])
            M = matrix(m).inverse()
            M.set_immutable()
            return M

        def _from_complete_on_basis(self, I):
            """
            Convert the Complete basis element indexed by ``I`` to ``self``.

            EXAMPLES::

                sage: ZL = NonCommutativeSymmetricFunctions(QQ).Zassenhaus_left()
                sage: ZL._from_complete_on_basis(Composition([1,3,2]))
                1/12*ZL[1, 1, 1, 1, 1, 1] + 1/6*ZL[1, 1, 1, 1, 2]
                 + 1/2*ZL[1, 2, 1, 1, 1] + ZL[1, 2, 1, 2]
                 + 1/2*ZL[1, 3, 1, 1] + ZL[1, 3, 2]
            """
            n = I.size()
            m = self._complete_to_zassenhaus_transition_matrix_inverse(n)
            C = Compositions(n)
            coeffs = m[C.rank(I)]
            return self._from_dict({J: coeffs[i] for i,J in enumerate(C)})

    ZL = Zassenhaus_left

    class Zassenhaus_right(CombinatorialFreeModule, BindableClass):
        r"""
        The Hopf algebra of non-commutative symmetric functions in
        the right Zassenhaus basis.

        This basis is defined in Section 2.5.1 of [HLNT09]_.
        It is multiplicative, with `Z_n` defined as the element of `NCSF_n`
        satisfying the equation

        .. MATH::

            \sigma_1 = exp(Z_1) exp(Z_2) exp(Z_3) \cdots exp(Z_n) \cdots

        where

        .. MATH::

            \sigma_1 = \sum_{n \geq 0} S_n .

        It can be recursively computed by the formula

        .. MATH::

            S_n = \sum_{\lambda \vdash n}
            \frac{1}{m_1(\lambda)! m_2(\lambda)! m_3(\lambda)! \cdots}
            \cdots Z_{\lambda_3} Z_{\lambda_2} Z_{\lambda_1}

        for all `n \geq 0`.

        Note that there is a variant (called the "noncommutative
        power sum symmetric functions of the third kind")
        in Definition 5.26 of [NCSF2]_ that satisfies:

        .. MATH::

            \sigma_1 = exp(Z_1) exp(Z_2/2) exp(Z_3/3) \cdots exp(Z_n/n) \cdots.
        """
        def __init__(self, NCSF):
            r"""
            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: ZR = NCSF.Zassenhaus_right(); ZR
                Non-Commutative Symmetric Functions over the Rational Field
                 in the Zassenhaus_right basis
                sage: TestSuite(ZR).run()

            TESTS:

            Test coproduct and antipode on the multiplicative identity::

                sage: ZR = NCSF.ZR()
                sage: ZR.one()
                ZR[]
                sage: ZR.one().coproduct()
                ZR[] # ZR[]
                sage: ZR.one().antipode()
                ZR[]

            We include some sanity tests to verify that conversions between
            this basis and other bases work the way they should::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: ZR = NCSF.Zassenhaus_right()
                sage: R = NCSF.ribbon()
                sage: S = NCSF.complete()
                sage: ZR(S[3])
                1/6*ZR[1, 1, 1] + ZR[1, 2] + ZR[3]
                sage: all(S(ZR(S[comp])) == S[comp] for comp in Compositions(5))
                True
                sage: all(ZR(S(ZR[comp])) == ZR[comp] for comp in Compositions(5))
                True
                sage: all(R(ZR(R[comp])) == R[comp] for comp in Compositions(5))
                True
                sage: all(ZR(R(ZR[comp])) == ZR[comp] for comp in Compositions(5))
                True
            """
            cat = NCSF.MultiplicativeBasesOnPrimitiveElements()
            CombinatorialFreeModule.__init__(self, NCSF.base_ring(), Compositions(),
                                                prefix='ZR', bracket=False,
                                                category=cat)

            # Register coercions
            S = self.realization_of().S()
            to_complete = self.algebra_morphism(self._to_complete_on_generator,
                                                codomain=S)
            to_complete.register_as_coercion()

            from_complete = S.module_morphism(on_basis=self._from_complete_on_basis,
                                              codomain=self)
            from_complete.register_as_coercion()

        def _to_complete_on_generator(self, n):
            r"""
            Expand a (right) Zassenhaus generator of non-commutative symmetric
            functions in the complete basis.

            INPUT:

            - ``n`` -- a positive integer

            OUTPUT:

            The expansion of the (right) Zassenhaus generator indexed by ``n``
            into the complete basis.

            TESTS::

                sage: ZR = NonCommutativeSymmetricFunctions(QQ).Zassenhaus_right()
                sage: ZR._to_complete_on_generator(1)
                S[1]
                sage: ZR._to_complete_on_generator(2)
                -1/2*S[1, 1] + S[2]
            """
            S = self.realization_of().S()

            if n <= 1:
                return S[n]

            from sage.combinat.partitions import ZS1_iterator
            from sage.rings.integer_ring import ZZ
            it = ZS1_iterator(n)
            next(it) # Skip the unique length 1 partition
            res = S[n]
            for p in it:
                d = {}
                for part in p:
                    d[part] = d.get(part, 0) + 1
                coeff = ZZ(prod(factorial(d[e]) for e in d))
                res = res - prod(self._to_complete_on_generator(i) for i in reversed(p)) / coeff
            return res

        @cached_method
        def _complete_to_zassenhaus_transition_matrix_inverse(self, n):
            r"""
            The change of basis matrix from the S basis to the ZR basis.

            EXAMPLES::

                sage: ZR = NonCommutativeSymmetricFunctions(QQ).Zassenhaus_right()
                sage: ZR._complete_to_zassenhaus_transition_matrix_inverse(3)
                [  1   0   0   0]
                [1/2   1   0   0]
                [1/2   0   1   0]
                [1/6   1   0   1]
            """
            from sage.matrix.constructor import matrix
            S = self.realization_of().S()
            m = []
            for I in Compositions(n):
                x = S(self.basis()[I])
                m.append([x.coefficient(J) for J in Compositions(n)])
            M = matrix(m).inverse()
            M.set_immutable()
            return M

        def _from_complete_on_basis(self, I):
            """
            Convert the Complete basis element indexed by ``I`` to ``self``.

            EXAMPLES::

                sage: ZR = NonCommutativeSymmetricFunctions(QQ).Zassenhaus_right()
                sage: ZR._from_complete_on_basis(Composition([1,3,2]))
                1/12*ZR[1, 1, 1, 1, 1, 1] + 1/6*ZR[1, 1, 1, 1, 2]
                 + 1/2*ZR[1, 1, 2, 1, 1] + ZR[1, 1, 2, 2]
                 + 1/2*ZR[1, 3, 1, 1] + ZR[1, 3, 2]
            """
            n = I.size()
            m = self._complete_to_zassenhaus_transition_matrix_inverse(n)
            C = Compositions(n)
            coeffs = m[C.rank(I)]
            return self._from_dict({J: coeffs[i] for i,J in enumerate(C)})

    ZR = Zassenhaus_right
