# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.graphs sage.combinat
r"""
Classical Crystals
"""
# ****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.tensor import TensorProductsCategory


class ClassicalCrystals(Category_singleton):
    """
    The category of classical crystals, that is crystals of finite Cartan type.

    EXAMPLES::

        sage: C = ClassicalCrystals()
        sage: C
        Category of classical crystals
        sage: C.super_categories()
        [Category of regular crystals,
         Category of finite crystals,
         Category of highest weight crystals]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = ClassicalCrystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: ClassicalCrystals().super_categories()
            [Category of regular crystals,
             Category of finite crystals,
             Category of highest weight crystals]
        """
        return [RegularCrystals(), FiniteCrystals(), HighestWeightCrystals()]

    def example(self, n=3):
        """
        Return an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = ClassicalCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        return Crystals().example(n)

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of classical crystals defines no
        additional structure: it only states that its objects are
        `U_q(\mathfrak{g})`-crystals, where `\mathfrak{g}` is of
        finite type.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: ClassicalCrystals().additional_structure()
        """
        return None

    class ParentMethods:

        def demazure_character(self, w, f=None):
            r"""
            Return the Demazure character associated to ``w``.

            INPUT:

            - ``w`` -- an element of the ambient weight lattice
              realization of the crystal, or a reduced word, or an element
              in the associated Weyl group

            OPTIONAL:

            - ``f`` -- a function from the crystal to a module

            This is currently only supported for crystals whose underlying
            weight space is the ambient space.

            The Demazure character is obtained by applying the Demazure operator
            `D_w` (see :meth:`sage.categories.regular_crystals.RegularCrystals.ParentMethods.demazure_operator`)
            to the highest weight element of the classical crystal. The simple
            Demazure operators `D_i` (see
            :meth:`sage.categories.regular_crystals.RegularCrystals.ElementMethods.demazure_operator_simple`)
            do not braid on the level of crystals, but on the level of characters they do.
            That is why it makes sense to input ``w`` either as a weight, a reduced word,
            or as an element of the underlying Weyl group.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape = [2,1])
                sage: e = T.weight_lattice_realization().basis()
                sage: weight = e[0] + 2*e[2]
                sage: weight.reduced_word()
                [2, 1]
                sage: T.demazure_character(weight)
                x1^2*x2 + x1*x2^2 + x1^2*x3 + x1*x2*x3 + x1*x3^2

                sage: T = crystals.Tableaux(['A',3],shape=[2,1])
                sage: T.demazure_character([1,2,3])
                x1^2*x2 + x1*x2^2 + x1^2*x3 + x1*x2*x3 + x2^2*x3
                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([1,2,3])
                sage: T.demazure_character(w)
                x1^2*x2 + x1*x2^2 + x1^2*x3 + x1*x2*x3 + x2^2*x3

                sage: T = crystals.Tableaux(['B',2], shape = [2])
                sage: e = T.weight_lattice_realization().basis()
                sage: weight = -2*e[1]
                sage: T.demazure_character(weight)
                x1^2 + x1*x2 + x2^2 + x1 + x2 + x1/x2 + 1/x2 + 1/x2^2 + 1

                sage: T = crystals.Tableaux("B2",shape=[1/2,1/2])
                sage: b2=WeylCharacterRing("B2",base_ring=QQ).ambient()
                sage: T.demazure_character([1,2],f=lambda x:b2(x.weight()))
                b2(-1/2,1/2) + b2(1/2,-1/2) + b2(1/2,1/2)

            REFERENCES:

            - [De1974]_

            - [Ma2009]_
            """
            from sage.misc.misc_c import prod
            from sage.rings.integer_ring import ZZ
            if hasattr(w, 'reduced_word'):
                word = w.reduced_word()
            else:
                word = w
            n = self.weight_lattice_realization().n
            u = self.algebra(ZZ).sum_of_monomials(self.module_generators)
            u = self.demazure_operator(u, word)
            if f is None:
                from sage.symbolic.ring import SR as P
                x = [P.var('x%s' % (i+1)) for i in range(n)]
                # TODO: use P.linear_combination when PolynomialRing will be a ModulesWithBasis
                return sum((coeff*prod((x[i]**(c.weight()[i]) for i in range(n)), P.one()) for c, coeff in u), P.zero())
            else:
                return sum(coeff * f(c) for c, coeff in u)

        def character(self, R=None):
            """
            Return the character of this crystal.

            INPUT:

            - ``R`` -- a :class:`WeylCharacterRing`
              (default: the default :class:`WeylCharacterRing` for this Cartan type)

            Returns the character of ``self`` as an element of ``R``.

            EXAMPLES::

                sage: C = crystals.Tableaux("A2", shape=[2,1])
                sage: chi = C.character(); chi
                A2(2,1,0)

                sage: T = crystals.TensorProduct(C,C)
                sage: chiT = T.character(); chiT
                A2(2,2,2) + 2*A2(3,2,1) + A2(3,3,0) + A2(4,1,1) + A2(4,2,0)
                sage: chiT == chi^2
                True

            One may specify an alternate :class:`WeylCharacterRing`::

                sage: R = WeylCharacterRing("A2", style='coroots')
                sage: chiT = T.character(R); chiT
                A2(0,0) + 2*A2(1,1) + A2(0,3) + A2(3,0) + A2(2,2)
                sage: chiT in R
                True

            It should have the same Cartan type and use the same
            realization of the weight lattice as ``self``::

                sage: R = WeylCharacterRing("A3", style='coroots')
                sage: T.character(R)
                Traceback (most recent call last):
                ...
                ValueError: Weyl character ring does not have the right Cartan type
            """
            from sage.combinat.root_system.weyl_characters import WeylCharacterRing
            if R is None:
                R = WeylCharacterRing(self.cartan_type())
            if not R.cartan_type() == self.cartan_type():
                raise ValueError("Weyl character ring does not have the right Cartan type")
            assert R.basis().keys() == self.weight_lattice_realization()

            return R.sum_of_monomials( x.weight() for x in self.highest_weight_vectors() )

        def __iter__(self):
            r"""
            Return an iterator over the elements of this crystal.

            This iterator uses little memory, storing only one element
            of the crystal at a time. For details on the complexity, see
            :class:`sage.combinat.crystals.crystals.CrystalBacktracker`.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: [x for x in C]
                [1, 2, 3, 4, 5, 6]

            TESTS::

                sage: C = crystals.Letters(['D',4])
                sage: D = crystals.SpinsPlus(['D',4])
                sage: E = crystals.SpinsMinus(['D',4])
                sage: T = crystals.TensorProduct(D,E,generators=[[D.list()[0],E.list()[0]]])
                sage: U = crystals.TensorProduct(C,E,generators=[[C(1),E.list()[0]]])
                sage: T.cardinality()
                56

                sage: TestSuite(T).run(verbose = True)
                running ._test_an_element() . . . pass
                running ._test_cardinality() . . . pass
                running ._test_category() . . . pass
                running ._test_construction() . . . pass
                running ._test_elements() . . .
                  Running the test suite of self.an_element()
                  running ._test_category() . . . pass
                  running ._test_eq() . . . pass
                  running ._test_new() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                  running ._test_stembridge_local_axioms() . . . pass
                  pass
                running ._test_elements_eq_reflexive() . . . pass
                running ._test_elements_eq_symmetric() . . . pass
                running ._test_elements_eq_transitive() . . . pass
                running ._test_elements_neq() . . . pass
                running ._test_enumerated_set_contains() . . . pass
                running ._test_enumerated_set_iter_cardinality() . . . pass
                running ._test_enumerated_set_iter_list() . . . pass
                running ._test_eq() . . . pass
                running ._test_fast_iter() . . . pass
                running ._test_new() . . . pass
                running ._test_not_implemented_methods() . . . pass
                running ._test_pickling() . . . pass
                running ._test_some_elements() . . . pass
                running ._test_stembridge_local_axioms() . . . pass

                sage: TestSuite(U).run(verbose = True)
                running ._test_an_element() . . . pass
                running ._test_cardinality() . . . pass
                running ._test_category() . . . pass
                running ._test_construction() . . . pass
                running ._test_elements() . . .
                  Running the test suite of self.an_element()
                  running ._test_category() . . . pass
                  running ._test_eq() . . . pass
                  running ._test_new() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                  running ._test_stembridge_local_axioms() . . . pass
                  pass
                running ._test_elements_eq_reflexive() . . . pass
                running ._test_elements_eq_symmetric() . . . pass
                running ._test_elements_eq_transitive() . . . pass
                running ._test_elements_neq() . . . pass
                running ._test_enumerated_set_contains() . . . pass
                running ._test_enumerated_set_iter_cardinality() . . . pass
                running ._test_enumerated_set_iter_list() . . . pass
                running ._test_eq() . . . pass
                running ._test_fast_iter() . . . pass
                running ._test_new() . . . pass
                running ._test_not_implemented_methods() . . . pass
                running ._test_pickling() . . . pass
                running ._test_some_elements() . . . pass
                running ._test_stembridge_local_axioms() . . . pass

            Bump's systematic tests::

                sage: fa3 = lambda a,b,c: crystals.Tableaux(['A',3],shape=[a+b+c,b+c,c])
                sage: fb3 = lambda a,b,c: crystals.Tableaux(['B',3],shape=[a+b+c,b+c,c])
                sage: fc3 = lambda a,b,c: crystals.Tableaux(['C',3],shape=[a+b+c,b+c,c])
                sage: fb4 = lambda a,b,c,d: crystals.Tableaux(['B',4],shape=[a+b+c+d,b+c+d,c+d,d])
                sage: fd4 = lambda a,b,c,d: crystals.Tableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d])
                sage: fd5 = lambda a,b,c,d,e: crystals.Tableaux(['D',5],shape=[a+b+c+d+e,b+c+d+e,c+d+e,d+e,e])
                sage: def fd4spinplus(a,b,c,d):
                ....:     C = crystals.Tableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d])
                ....:     D = crystals.SpinsPlus(['D',4])
                ....:     return crystals.TensorProduct(C,D,generators=[[C[0],D[0]]])
                sage: def fb3spin(a,b,c):
                ....:     C = crystals.Tableaux(['B',3],shape=[a+b+c,b+c,c])
                ....:     D = crystals.Spins(['B',3])
                ....:     return crystals.TensorProduct(C,D,generators=[[C[0],D[0]]])

            .. TODO::

                Choose a good panel of values for `a,b,c ...` both for
                basic systematic tests and for conditionally run,
                computationally involved tests.

            ::

                sage: TestSuite(fb4(1,0,1,0)).run(verbose = True)  # long time (8s on sage.math, 2011)
                running ._test_an_element() . . . pass
                running ._test_cardinality() . . . pass
                running ._test_category() . . . pass
                running ._test_construction() . . . pass
                running ._test_elements() . . .
                  Running the test suite of self.an_element()
                  running ._test_category() . . . pass
                  running ._test_eq() . . . pass
                  running ._test_new() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                  running ._test_stembridge_local_axioms() . . . pass
                  pass
                running ._test_elements_eq_reflexive() . . . pass
                running ._test_elements_eq_symmetric() . . . pass
                running ._test_elements_eq_transitive() . . . pass
                running ._test_elements_neq() . . . pass
                running ._test_enumerated_set_contains() . . . pass
                running ._test_enumerated_set_iter_cardinality() . . . pass
                running ._test_enumerated_set_iter_list() . . . pass
                running ._test_eq() . . . pass
                running ._test_fast_iter() . . . pass
                running ._test_new() . . . pass
                running ._test_not_implemented_methods() . . . pass
                running ._test_pickling() . . . pass
                running ._test_some_elements() . . . pass
                running ._test_stembridge_local_axioms() . . . pass

            ::

                #sage: fb4(1,1,1,1).check() # expensive: the crystal is of size 297297
                #True
            """
            from sage.combinat.crystals.crystals import CrystalBacktracker
            return iter(CrystalBacktracker(self))

        def _test_fast_iter(self, **options):
            r"""
            Test whether the elements returned by :meth:`.__iter__`
            and ``Crystal.list(self)`` are the same (the two
            algorithms are different).

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C._test_fast_iter()
            """
            tester = self._tester(**options)
            S = list(self)
            SS = list(Crystals().parent_class.__iter__(self))
            tester.assertEqual(len(S), len(SS))
            tester.assertEqual(len(S), len(set(S)))
            tester.assertEqual(set(S), set(SS))

        def cardinality(self):
            r"""
            Return the number of elements of the crystal, using Weyl's
            dimension formula on each connected component.

            EXAMPLES::

                sage: C = ClassicalCrystals().example(5)
                sage: C.cardinality()
                6
            """
            return sum(self.weight_lattice_realization().weyl_dimension(x.weight())
                       for x in self.highest_weight_vectors())

    class ElementMethods:

        def lusztig_involution(self):
            r"""
            Return the Lusztig involution on the classical highest weight
            crystal ``self``.

            The Lusztig involution on a finite-dimensional highest weight
            crystal `B(\lambda)` of highest weight `\lambda` maps the
            highest weight vector to the lowest weight vector and the
            Kashiwara operator `f_i` to `e_{i^*}`, where `i^*` is defined as
            `\alpha_{i^*} = -w_0(\alpha_i)`. Here `w_0` is the longest element
            of the Weyl group acting on the `i`-th simple root `\alpha_i`.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',3],shape=[2,1])
                sage: b = B(rows=[[1,2],[4]])
                sage: b.lusztig_involution()
                [[1, 4], [3]]
                sage: b.to_tableau().schuetzenberger_involution(n=4)
                [[1, 4], [3]]

                sage: all(b.lusztig_involution().to_tableau() == b.to_tableau().schuetzenberger_involution(n=4) for b in B)
                True

                sage: B = crystals.Tableaux(['D',4],shape=[1])
                sage: [[b,b.lusztig_involution()] for b in B]
                [[[[1]], [[-1]]], [[[2]], [[-2]]], [[[3]], [[-3]]], [[[4]], [[-4]]], [[[-4]],
                [[4]]], [[[-3]], [[3]]], [[[-2]], [[2]]], [[[-1]], [[1]]]]

                sage: B = crystals.Tableaux(['D',3],shape=[1])
                sage: [[b,b.lusztig_involution()] for b in B]
                [[[[1]], [[-1]]], [[[2]], [[-2]]], [[[3]], [[3]]], [[[-3]], [[-3]]],
                [[[-2]], [[2]]], [[[-1]], [[1]]]]

                sage: C = CartanType(['E',6])
                sage: La = C.root_system().weight_lattice().fundamental_weights()
                sage: T = crystals.HighestWeight(La[1])
                sage: t = T[3]; t
                [(-4, 2, 5)]
                sage: t.lusztig_involution()
                [(-2, -3, 4)]
            """
            hw = self.to_highest_weight()[1]
            hw.reverse()
            A = self.parent().cartan_type().opposition_automorphism()
            hw = [A[i] for i in hw]
            return self.to_lowest_weight()[0].e_string(hw)

    class TensorProducts(TensorProductsCategory):
        """
        The category of classical crystals constructed by tensor
        product of classical crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: ClassicalCrystals().TensorProducts().extra_super_categories()
                [Category of classical crystals]
            """
            return [self.base_category()]
