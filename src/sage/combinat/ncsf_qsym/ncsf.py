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
# 1. Make Coersion run faster between multiple bases.
########################################

from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.rings import Rings
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.combinat.composition import Composition, Compositions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.ncsf_qsym.generic_basis_code import BasesOfQSymOrNCSF
from sage.combinat.ncsf_qsym.combinatorics import *
from sage.combinat.partition import Partition
from sage.combinat.permutation import Permutations

class NonCommutativeSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The abstract algebra of non-commutative symmetric functions

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
        Join of Category of graded hopf algebras over Rational Field and Category of monoids with realizations and Category of coalgebras over Rational Field with realizations
        sage: [S[i].degree() for i in range(10)]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    We use the Sage standard renaming idiom to get shorter outputs::

        sage: NCSF.rename("NCSF")
        sage: NCSF
        NCSF

    NCSF has many representations as a concrete algebra. Each of them
    has a distinguished basis, and its elements are expanded in this
    basis. Here is the Psi representation::

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

    .. rubric:: Other concrete representations

    .. TODO:: demonstrate how to customize the basis names

    NCSF admits many other concrete realizations::

        sage: Phi        = NCSF.Phi()
        sage: ribbon     = NCSF.ribbon()
        sage: complete   = NCSF.complete()
        sage: elementary = NCSF.elementary()
        sage: monomial   = NCSF.monomial()

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
    Moebius function for the boolean lattice, the inverse change of
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

    By Moebius inversion on the composition poset, the ribbon
    basis element corresponding to a composition `I` is then the
    alternating sum over all compositions fatter than the
    complement composition of `I` in the elementary basis::

        sage: elementary(ribbon[2,1,2,1])
        L[1, 3, 2] - L[1, 5] - L[4, 2] + L[6]

    .. TODO:: explain the other changes of bases!

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
        Category of hom sets in Category of modules with basis over Rational Field
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

    .. TODO::

        - Bases: monomial, fundamental, forgotten, quasi_schur_dual
          simple() ? (<=> simple modules of HS_n; to be discussed with Florent)
        - Multiplication in:

          - fundamental and monomial (cf. Lenny's code)
          - ribbon (from Mike's code)

        - Coproducts (most done by coercions)
        - some_elements in all bases
        - Systematic coercion checks (in AlgebrasWithBasis().Abstract())
    """
    def __init__(self, R):
        r"""
        TESTS::

            sage: TestSuite(NonCommutativeSymmetricFunctions(QQ)).run()
        """
        assert(R in Rings())
        self._base = R # Won't be needed once CategoryObject won't override base_ring
        Parent.__init__(self, category = GradedHopfAlgebras(R).WithRealizations())

        # COERCION METHODS
        Psi = self.Psi()
        Phi = self.Phi()
        complete = self.complete()
        elementary = self.elementary()
        ribbon = self.ribbon()

        # complete to ribbon, and back
        complete  .module_morphism(ribbon.sum_of_fatter_compositions,               codomain=ribbon    ).register_as_coercion()
        ribbon    .module_morphism(complete.alternating_sum_of_fatter_compositions, codomain=complete  ).register_as_coercion()


        # complete to elementary, and back (should be constructed from _on_generators?)
        complete  .module_morphism(complete._from_elementary_on_basis,              codomain=elementary).register_as_coercion()
        elementary.module_morphism(elementary._from_complete_on_basis,              codomain=complete  ).register_as_coercion()

        # complete to Psi, and back (should be constructed from _on_generators?)
        complete  .module_morphism(Psi._from_complete_on_basis,                     codomain=Psi       ).register_as_coercion()
        Psi       .module_morphism(Psi._to_complete_on_basis,                       codomain=complete  ).register_as_coercion()

        # complete to Phi, and back (should be constructed from _on_generators?)
        complete  .module_morphism(Phi._from_complete_on_basis,                     codomain=Phi       ).register_as_coercion()
        Phi       .module_morphism(Phi._to_complete_on_basis,                       codomain=complete  ).register_as_coercion()

    def _repr_(self): # could be taken care of by the category
        r"""
        EXAMPLES::

            sage: N = NonCommutativeSymmetricFunctions(ZZ)
            sage: N._repr_()
            'Non-Commutative Symmetric Functions over the Integer Ring'
        """
        return "Non-Commutative Symmetric Functions over the %s"%self.base_ring()

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

    _shorthands = tuple(['S', 'R', 'L', 'Phi', 'Psi', 'nM', 'I'])

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
            from generic_basis_code import GradedModulesWithInternalProduct
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

                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S.to_symmetric_function(S[2,1])
                    h[2, 1]
                    sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                    sage: R.to_symmetric_function_on_basis(Composition([2,1]))
                    s[2, 1]
                """
                S = self.realization_of().complete()
                return S.to_symmetric_function(S(self[I]))

            @lazy_attribute
            def to_symmetric_function(self):
                r"""
                Morphism of ``self`` to the algebra of symmetric functions.

                This is constructed by extending the method
                :meth:`to_symmetric_function_on_basis` linearly.

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
                    sage: S = N.complete()
                    sage: S.to_symmetric_function(S[3,1,2])
                    h[3, 2, 1]
                    sage: Phi = N.Phi()
                    sage: Phi.to_symmetric_function(Phi[1,3])
                    h[1, 1, 1, 1] - 3*h[2, 1, 1] + 3*h[3, 1]
                    sage: R.to_symmetric_function
                    Generic morphism:
                      From: Non-Commutative Symmetric Functions over the Rational Field in the Ribbon basis
                      To:   Symmetric Functions over Rational Field in the Schur basis
                """
                on_basis = self.to_symmetric_function_on_basis
                codomain = on_basis(self.one_basis()).parent()
                return self.module_morphism(on_basis=on_basis, codomain=codomain)

            def to_ncsym_on_basis(self, I):
                r"""
                The image of the basis element indexed by ``I`` under the
                map `\kappa` to the symmetric functions in non-commuting
                variables such that for the natural maps `\chi : NCSym \to Sym`
                and `\rho : NSym \to Sym`, we have `\chi \circ \kappa = \rho`.

                This default implementation does a change of basis and
                computes the image in the complete basis.

                INPUT:

                - ``I`` -- a composition

                EXAMPLES::

                    sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S.to_ncsym(S[2,1])
                    1/2*m{{1}, {2}, {3}} + 1/2*m{{1}, {2, 3}} + m{{1, 2}, {3}}
                     + m{{1, 2, 3}} + 1/2*m{{1, 3}, {2}}
                    sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                    sage: R.to_ncsym_on_basis(Composition([2,1]))
                    1/3*m{{1}, {2}, {3}} + 1/6*m{{1}, {2, 3}} + 2/3*m{{1, 2}, {3}} + 1/6*m{{1, 3}, {2}}
                """
                S = self.realization_of().complete()
                return S.to_ncsym(S(self[I]))

            @lazy_attribute
            def to_ncsym(self):
                r"""
                Morphism `\kappa` of ``self`` to the algebra of symmetric
                functions in non-commuting variables that for the natural
                maps `\chi : NCSym \to Sym` and `\rho : NSym \to Sym`, we
                have `\chi \circ \kappa = \rho`.

                This is constructed by extending the method
                :meth:`to_ncsym_on_basis` linearly.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: x = R.an_element(); x
                    2*R[] + 2*R[1] + 3*R[1, 1]
                    sage: R.to_ncsym(x)
                    2*m{} + 2*m{{1}} + 3/2*m{{1}, {2}}
                    sage: S = N.complete()
                    sage: S.to_ncsym(S[1,2])
                    1/2*m{{1}, {2}, {3}} + m{{1}, {2, 3}} + 1/2*m{{1, 2}, {3}}
                     + m{{1, 2, 3}} + 1/2*m{{1, 3}, {2}}
                    sage: Phi = N.Phi()
                    sage: Phi.to_ncsym(Phi[1,3])
                    -1/4*m{{1}, {2}, {3, 4}} - 1/4*m{{1}, {2, 3}, {4}} + m{{1}, {2, 3, 4}}
                     + 1/2*m{{1}, {2, 4}, {3}} - 1/4*m{{1, 2}, {3, 4}} - 1/4*m{{1, 2, 3}, {4}}
                     + m{{1, 2, 3, 4}} + 1/2*m{{1, 2, 4}, {3}} + 1/2*m{{1, 3}, {2, 4}}
                     - 1/4*m{{1, 3, 4}, {2}} - 1/4*m{{1, 4}, {2, 3}}
                    sage: R.to_ncsym
                    Generic morphism:
                      From: Non-Commutative Symmetric Functions over the Rational Field in the Ribbon basis
                      To:   Symmetric functions in non-commuting variables over the Rational Field in the monomial basis
                """
                from sage.combinat.ncsym.ncsym import SymmetricFunctionsNonCommutingVariables
                codomain = SymmetricFunctionsNonCommutingVariables(self.base_ring()).monomial()
                return self.module_morphism(self.to_ncsym_on_basis, codomain=codomain)

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
                    \quad \mathbf{V}_n(\Psi_r) = n \Psi_r,
                    \quad \mathbf{V}_n(\Phi_r) = n \Phi_r

                (where `S_r` denotes the `r`-th complete non-commutative
                symmetric function, `\Lambda_r` denotes the `r`-th elementary
                non-commutative symmetric function, `\Psi_r` denotes the `r`-th
                power-sum non-commutative symmetric function of the first kind,
                and `\Phi_r` denotes the `r`-th power-sum non-commutative
                symmetric function of the second kind). For every positive
                integer `r` with `n \not\mid r`, it satisfes

                .. MATH::

                    \mathbf{V}_n(S_r) = \mathbf{V}_n(\Lambda_r)
                    = \mathbf{V}_n(\Psi_r) = \mathbf{V}_n(\Phi_r) = 0.

                The `n`-th Verschiebung operator is also called the `n`-th
                Verschiebung endomorphism.

                It is a lift of the `n`-th Verschiebung operator on the ring
                of symmetric functions (
                :meth:`sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung`
                ) to the ring of noncommutative symmetric functions.

                The action of the `n`-th Verschiebung operator can also be
                described on the ribbon Schur functions. Namely, every
                composition `I` of size `n \ell` satisfies

                .. MATH::

                    \mathbf{V}_n ( R_I )
                    = (-1)^{\ell(I) - \ell(J)}
                    \cdot R_{J/n},

                where `J` denotes the meet of the compositions `I` and
                `(\underbrace{n, n, \ldots, n}_{\ell\mbox{ times}})`,
                where `\ell(I)` is the length of `I`, and
                where `J / n` denotes the composition obtained by dividing
                every entry of `J` by `n`.
                For a composition `I` of size not divisible by `n`, we have
                `\mathbf{V}_n( R_I ) = 0`.

                .. SEEALSO::

                    :meth:`sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.verschiebung`

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
                dct = {Composition(map(lambda i: i // n, I)): coeff
                       for (I, coeff) in S(self).monomial_coefficients().items()
                       if all(i % n == 0 for i in I)}
                return parent(S._from_dict(dct))

            def to_descent_algebra(self, n):
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
                    sage: (S[1,2,1] - 3 * S[1,1,2]).to_descent_algebra(1)
                    0
                """
                from sage.combinat.descent_algebra import DescentAlgebra
                S = NonCommutativeSymmetricFunctions(self.base_ring()).S()
                S_expansion = S(self)
                B = DescentAlgebra(self.base_ring(), n).B()
                return B.sum(coeff * B[I] for I, coeff in S_expansion.monomial_coefficients().items() if sum(I) == n)

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
                Return the image of ``self`` in the symmetric functions in
                non-commuting variables under the map that fixes the usual
                symmetric functions.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: x = R.an_element(); x
                    2*R[] + 2*R[1] + 3*R[1, 1]
                    sage: x.to_ncsym()
                    2*m{} + 2*m{{1}} + 3/2*m{{1}, {2}}
                    sage: y = N.Phi()[1,2]
                    sage: y.to_ncsym()
                    m{{1}, {2, 3}} + m{{1, 2, 3}}
                """
                return self.parent().to_ncsym(self)

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
                return Family(PositiveIntegers(), lambda i: self.monomial(Composition([i])))

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
                    sage: def double(i) : return Psi[i,i]
                    ...
                    sage: f = Psi.algebra_morphism(double, codomain = Psi)
                    sage: f
                    Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
                    sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
                    2*Psi[] + 3*Psi[1, 1, 3, 3, 2, 2] + Psi[2, 2, 4, 4]
                    sage: f.category()
                    Join of Category of hom sets in Category of modules with basis over Rational Field and Category of hom sets in Category of rings

                When extra properties about the morphism are known, one
                can specify the category of which it is a morphism::

                    sage: def negate(i): return -Psi[i]
                    sage: f = Psi.algebra_morphism(negate, codomain = Psi, category = GradedHopfAlgebrasWithBasis(QQ))
                    sage: f
                    Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
                    sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
                    2*Psi[] - 3*Psi[1, 3, 2] + Psi[2, 4]
                    sage: f.category()
                    Join of Category of hom sets in Category of modules with basis over Rational Field and Category of hom sets in Category of rings

                If ``anti`` is true, this returns an anti-algebra morphism::

                    sage: f = Psi.algebra_morphism(double, codomain = Psi, anti=True)
                    sage: f
                    Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
                    sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
                    2*Psi[] + 3*Psi[2, 2, 3, 3, 1, 1] + Psi[4, 4, 2, 2]
                    sage: f.category()
                    Category of hom sets in Category of modules with basis over Rational Field
                """
                from sage.combinat.ncsf_qsym.generic_basis_code import AlgebraMorphism
                return AlgebraMorphism(self, on_generators, **keywords)

            @lazy_attribute
            def antipode(self):
                r"""
                Return the antipode morphism on the basis ``self``.

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
                    sage: S[3].antipode() # indirect docttest
                    -S[1, 1, 1] + S[1, 2] + S[2, 1] - S[3]
                    sage: S[2,3].coproduct().apply_multilinear_morphism(lambda be,ga: S(be)*S(ga).antipode())
                    0
                    sage: S[2,3].coproduct().apply_multilinear_morphism(lambda be,ga: S(be).antipode()*S(ga))
                    0
                """
                # TODO: avoid this -1^... by using properly

                return (-1)**len(composition) * self.alternating_sum_of_finer_compositions(composition.reversed())

            # @cached_method?
            def coproduct_on_generators(self, i):
                """
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
                    sage: S.coproduct_on_generators(0)
                    'Not a positive integer: 0'
                """
                if i<1:
                    return "Not a positive integer: %s" % `i`
                def C(i): return Composition([i]) if i else Composition([])
                T = self.tensor_square()
                return T.sum_of_monomials( (C(j), C(i-j)) for j in range(0,i+1) )

    class MultiplicativeBasesOnPrimitiveElements(Category_realization_of_parent):
        """
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
                    sage: Psi.antipode_on_generators(0)
                    'Not a positive integer: 0'
                """
                if i<1:
                    return "Not a positive integer: %s" % `i`
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
                    sage: Psi.coproduct_on_generators(0)
                    'Not a positive integer: 0'
                """
                if i<1:
                    return "Not a positive integer: %s" % `i`
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
            if I == []:
                return self.monomial(J)
            elif J == []:
                return self.monomial(I)
            else:
                return self.monomial(Composition(I[:] + J[:])) + \
                       self.monomial(Composition(I[:-1] + [I[-1]+J[0]] + J[1:]))

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
            from sage.combinat.sf.sf import SymmetricFunctions
            s = SymmetricFunctions(self.base_ring()).schur()
            if I == []:
                return s([])
            return s(I.to_skew_partition())

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

        def _from_elementary_on_basis(self, I):
            r"""
            Expand an elementary basis element of non-commutative symmetric
            functions in the complete basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The expansion of the elementary function indexed by ``I`` into
              the complete basis.

            EXAMPLES::

                sage: S=NonCommutativeSymmetricFunctions(QQ).S()
                sage: S._from_elementary_on_basis(Composition([2,1]))
                S[1, 1, 1] - S[2, 1]
                sage: S._from_elementary_on_basis(Composition([]))
                S[]
            """
            n = I.size()
            minus_one = -self.base_ring().one()
            return self.sum_of_terms( (compo, minus_one**(len(compo)-n)) for compo in I.finer() )

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
            The internal product of two non-commutative symmetric complete functions.

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
            The commutative image of a complete non-commutative symmetric function basis
            element. This is obtained by sorting the composition.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The commutative image of the complete basis element indexed by
              ``I``. The result is the complete symmetric function indexed by
              the partition obtained by sorting ``I``.

            EXAMPLES::

                sage: S=NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.to_symmetric_function_on_basis([2,1,3])
                h[3, 2, 1]
                sage: S.to_symmetric_function_on_basis([])
                h[]
            """
            from sage.combinat.sf.sf import SymmetricFunctions
            h = SymmetricFunctions(self.base_ring()).homogeneous()
            return h[Partition(sorted(I,reverse=True))]

        def _to_symmetric_group_algebra_on_basis(self, I):
            r"""
            Return the image of the complete non-commutative symmetric function indexed
            by the composition ``I`` in the symmetric group algebra under the canonical
            embedding of the non-commutative symmetric functions into the symmetric
            group algebra.

            This embedding sends the complete basis element indexed by the composition
            ``I`` to the sum of all permutations whose descent composition is fatter
            than ``I`` (that is, all permutations whose right descent set is contained
            in the subset corresponding to ``I``).

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The sum of all permutations with right descent set contained in ``I``.

            EXAMPLES::

                sage: S=NonCommutativeSymmetricFunctions(QQ).S()
                sage: S._to_symmetric_group_algebra_on_basis(Composition([1,2]))
                [1, 2, 3] + [2, 1, 3] + [3, 1, 2]
                sage: S._to_symmetric_group_algebra_on_basis(Composition([]))
                []
            """
            n = sum(I)
            from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
            from sage.sets.set import Set
            if n == 0:
                return SymmetricGroupAlgebra(self.base_ring(),n).one()
            sga = SymmetricGroupAlgebra(self.base_ring(),n)
            x = sga.zero()
            J = [j-1 for j in I.to_subset()]
            for K in Set(J).subsets():
                for p in Permutations(descents=(K,n)):
                    x += sga(p)
            return x

        def to_ncsym_on_basis(self, I):
            r"""
            Return the image of the complete non-commutative symmetric function
            in the symmetric functions in non-commuting variables under the
            embedding `\mathcal{I}` which fixes the symmetric functions.

            This map is defined by

            .. MATH::

                S_n \mapsto \sum_{A \vdash [n]}
                \frac{\lambda(A)! \lambda(A)^!}{n!} \mathbf{m}_A

            and extended as an algebra homomorphism.

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.to_ncsym_on_basis(Composition([2]))
                1/2*m{{1}, {2}} + m{{1, 2}}
                sage: S.to_ncsym_on_basis(Composition([1,2,1]))
                1/2*m{{1}, {2}, {3}, {4}} + 1/2*m{{1}, {2}, {3, 4}} + m{{1}, {2, 3}, {4}}
                 + m{{1}, {2, 3, 4}} + 1/2*m{{1}, {2, 4}, {3}} + 1/2*m{{1, 2}, {3}, {4}}
                 + 1/2*m{{1, 2}, {3, 4}} + m{{1, 2, 3}, {4}} + m{{1, 2, 3, 4}}
                 + 1/2*m{{1, 2, 4}, {3}} + 1/2*m{{1, 3}, {2}, {4}} + 1/2*m{{1, 3}, {2, 4}}
                 + 1/2*m{{1, 3, 4}, {2}} + 1/2*m{{1, 4}, {2}, {3}} + m{{1, 4}, {2, 3}}
                sage: S.to_ncsym_on_basis(Composition([]))
                m{}

            TESTS:

            Check that the image under `\mathcal{I}` fixes the
            symmetric functions::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: m = SymmetricFunctionsNonCommutingVariables(QQ).monomial()
                sage: mon = SymmetricFunctions(QQ).monomial()
                sage: all(S[c].to_ncsym().to_symmetric_function() == S[c].to_symmetric_function()
                ....:     for i in range(5) for c in Compositions(i))
                True

            We also check that the `NCSym` monomials agree on the homogeneous
            and complete basis::

                sage: h = SymmetricFunctions(QQ).h()
                sage: all(m.from_symmetric_function(h[i]) == S[i].to_ncsym() for i in range(6))
                True
            """
            from sage.combinat.ncsym.ncsym import SymmetricFunctionsNonCommutingVariables
            m = SymmetricFunctionsNonCommutingVariables(self.base_ring()).m()
            if I == []:
                return m.one()

            from sage.combinat.set_partition import SetPartitions
            R = self.base_ring()
            P = SetPartitions()
            c_num = lambda A: prod([factorial(i) for i in A.shape()], R.one())
            return prod(m.sum_of_terms([(P(A), R(c_num(A) / factorial(n))) for A in SetPartitions(n)], distinct=True)
                        for n in I)

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

        # TODO: use alternating_sum_of_finer_compositions, if possible getting rid of this method
        def _from_complete_on_basis(self, I):
            r"""
            Expand a complete basis element of non-commutative symmetric functions
            in the elementary basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The expansion of the complete function indexed by ``I`` into the
              elementary basis.

            EXAMPLES::

                sage: L=NonCommutativeSymmetricFunctions(QQ).L()
                sage: L._from_complete_on_basis(Composition([2,1]))
                L[1, 1, 1] - L[2, 1]
                sage: L._from_complete_on_basis(Composition([]))
                L[]
            """
            n = I.size()
            minus_one = -self.base_ring().one()
            return self.sum_of_terms( (compo, minus_one**(len(compo)-n)) for compo in I.finer() )

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

        EXAMPLES::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: Psi = NCSF.Psi(); Psi
            Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
            sage: Psi.an_element()
            2*Psi[] + 2*Psi[1] + 3*Psi[1, 1]
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

        # TODO: should those be defined using algebra morphism?
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
            I = Composition([n])
            # TODO: I being trivial, there is no refinement going on here, so
            # one can probably be a bit more explicit / fast
            return self.sum_of_terms((J, one/coeff_pi(J,I)) for J in Compositions(n))

        def _from_complete_on_basis(self, I):
            r"""
            Expand a complete basis element of non-commutative symmetric functions
            in the Psi basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The expansion of the complete function indexed by ``I`` in the
              Psi basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi._from_complete_on_basis(Composition([1]))
                Psi[1]
                sage: Psi._from_complete_on_basis(Composition([2]))
                1/2*Psi[1, 1] + 1/2*Psi[2]
                sage: Psi._from_complete_on_basis(Composition([3]))
                1/6*Psi[1, 1, 1] + 1/3*Psi[1, 2] + 1/6*Psi[2, 1] + 1/3*Psi[3]
                sage: Psi._from_complete_on_basis(Composition([2,1]))
                1/2*Psi[1, 1, 1] + 1/2*Psi[2, 1]
                sage: Psi._from_complete_on_basis(Composition([1,2]))
                1/2*Psi[1, 1, 1] + 1/2*Psi[1, 2]
                sage: Psi._from_complete_on_basis(Composition([1,1,1]))
                Psi[1, 1, 1]
            """
            # TODO: make this comment into a reference in the doctest (same thing elsewhere)
            # Proposition 4.5 of NCSF I article
            one = self.base_ring().one()
            return self.sum_of_terms((J, one/coeff_pi(J,I)) for J in I.finer())

        def _to_complete_on_basis(self, I):
            r"""
            Expand a Psi basis element of non-commutative symmetric functions
            in the complete basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The expansion of the Psi function indexed by ``I`` in the
              complete basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi._to_complete_on_basis(Composition([1]))
                S[1]
                sage: Psi._to_complete_on_basis(Composition([2]))
                -S[1, 1] + 2*S[2]
                sage: Psi._to_complete_on_basis(Composition([1,1]))
                S[1, 1]
                sage: Psi._to_complete_on_basis(Composition([3]))
                S[1, 1, 1] - 2*S[1, 2] - S[2, 1] + 3*S[3]
                sage: Psi._to_complete_on_basis(Composition([2,1]))
                -S[1, 1, 1] + 2*S[2, 1]
                sage: Psi._to_complete_on_basis(Composition([1,2]))
                -S[1, 1, 1] + 2*S[1, 2]
                sage: Psi._to_complete_on_basis(Composition([1,1,1]))
                S[1, 1, 1]
            """
            # Proposition 4.5 of NCSF I article
            minus_one = -self.base_ring().one()
            complete = self.realization_of().complete()
            return complete.sum_of_terms((J, minus_one**(len(J)-len(I))*coeff_lp(J,I))
                                            for J in I.finer())

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

        def _from_complete_on_basis(self, I):
            r"""
            Expand a complete basis element of non-commutative symmetric
            functions in the Phi basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The expansion of the complete function indexed by ``I`` in the
              Phi basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Phi = NonCommutativeSymmetricFunctions(QQ).Phi()
                sage: Phi._from_complete_on_basis(Composition([1]))
                Phi[1]
                sage: Phi._from_complete_on_basis(Composition([2]))
                1/2*Phi[1, 1] + 1/2*Phi[2]
                sage: Phi._from_complete_on_basis(Composition([1,1]))
                Phi[1, 1]
                sage: Phi._from_complete_on_basis(Composition([3]))
                1/6*Phi[1, 1, 1] + 1/4*Phi[1, 2] + 1/4*Phi[2, 1] + 1/3*Phi[3]
                sage: Phi._from_complete_on_basis(Composition([2,1]))
                1/2*Phi[1, 1, 1] + 1/2*Phi[2, 1]
                sage: Phi._from_complete_on_basis(Composition([1,2]))
                1/2*Phi[1, 1, 1] + 1/2*Phi[1, 2]
                sage: Phi._from_complete_on_basis(Composition([1,1,1]))
                Phi[1, 1, 1]
            """
            # Proposition 4.9 of NCSF I article
            one = self.base_ring().one()
            return self.sum_of_terms((J, one / coeff_sp(J,I)) for J in I.finer())

        def _to_complete_on_basis(self, I):
            r"""
            Expand a Phi basis element of non-commutative symmetric functions
            in the complete basis.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The expansion of the Phi function indexed by ``I`` in the
              complete basis.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: Phi = NonCommutativeSymmetricFunctions(QQ).Phi()
                sage: Phi._to_complete_on_basis(Composition([1]))
                S[1]
                sage: Phi._to_complete_on_basis(Composition([2]))
                -S[1, 1] + 2*S[2]
                sage: Phi._to_complete_on_basis(Composition([1,1]))
                S[1, 1]
                sage: Phi._to_complete_on_basis(Composition([3]))
                S[1, 1, 1] - 3/2*S[1, 2] - 3/2*S[2, 1] + 3*S[3]
                sage: Phi._to_complete_on_basis(Composition([2,1]))
                -S[1, 1, 1] + 2*S[2, 1]
                sage: Phi._to_complete_on_basis(Composition([1,2]))
                -S[1, 1, 1] + 2*S[1, 2]
                sage: Phi._to_complete_on_basis(Composition([1,1,1]))
                S[1, 1, 1]
            """
            # Proposition 4.9 of NCSF I article
            minus_one = -self.base_ring().one()
            complete = self.realization_of().complete()
            return complete.sum_of_terms((J, minus_one**(len(J)-len(I)) * prod(I) / coeff_ell(J,I))
                                         for J in I.finer())

    class Monomial(CombinatorialFreeModule, BindableClass):
        def __init__(self, NCSF):
            r"""
            The monomial basis defined in Tevlin's paper [Tev2007]_. It
            is the basis denoted by `(M^I)_I` in that paper.

            The Monomial basis is well-defined only when the base ring is a
            `\QQ`-algebra.

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

            - ``self`` - the Monomial basis of non-commutative symmetric functions
            - ``I`` - a composition

            OUTPUT:

            - The expansion of the Monomial function indexed by ``I`` in the complete
              basis.

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
            return sum(m_to_s_stat(self.base_ring(),I,K) * S(K) for K in Compositions(Composition(I).size()))

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
                    p = [0] + Composition(I).refinement_splitting_lengths(J).partial_sums()
                    sum_of_elements += prod( (len_of_J - k)**(p[k+1]-p[k]) for k in range(len_of_J) ) * M(J)
            return sum_of_elements

    nM = monomial = Monomial

    class Immaculate(CombinatorialFreeModule, BindableClass):
        def __init__(self, NCSF):
            r"""
            The immaculate basis of the non-commutative symmetric
            functions. This basis first appears in Berg, Bergeron,
            Saliola, Serrano and Zabrocki's [BBSSZ2012]_.

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

            - ``self`` - The Immaculate basis
            - ``alpha`` - a list

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
                return 0
            return S( [ d for d in alpha if d > 0 ] )

        @cached_method
        def _to_complete_on_basis(self, alpha):
            r"""
            Return the expansion of an Immaculate basis element in the
            complete basis.

            INPUT:

            - ``self`` - the Immaculate basis of the non-commutative
              symmetric functions
            - ``alpha`` - a composition

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
            if alpha == []:
                return self._H([])
            if alpha == [1]:
                return self._H([1])
            return sum( sigma.signature()*self._H( [alpha[i]+sigma[i]-(i+1) for i in range(len(alpha))] ) for sigma in Permutations(len(alpha)))

        @cached_method
        def _from_complete_on_basis(self, comp_content):
            r"""
            Return the expansion of a complete basis element in the
            Immaculate basis.

            INPUT:

            - ``comp_content`` - a composition

            OUTPUT:

            - The expansion in the Immaculate basis of the basis element
              of the complete basis indexed by the composition
              ``comp_content`` .

            EXAMPLES::

                sage: I=NonCommutativeSymmetricFunctions(QQ).I()
                sage: I._from_complete_on_basis(Composition([]))
                I[]
                sage: I._from_complete_on_basis(Composition([2,1,3]))
                I[2, 1, 3] + I[2, 2, 2] + I[2, 3, 1] + I[2, 4] + I[3, 1, 2] + I[3, 2, 1] + 2*I[3, 3] + I[4, 1, 1] + 2*I[4, 2] + 2*I[5, 1] + I[6]
            """
            I = NonCommutativeSymmetricFunctions(self.base_ring()).I()
            if comp_content == []:
                return I([])
            else:
                return sum( number_of_fCT(comp_content,comp_shape) * I(comp_shape) \
                            for comp_shape in Compositions(Composition(comp_content).size()) )

    I = Immaculate
