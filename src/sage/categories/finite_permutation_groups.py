r"""
Finite Permutation Groups
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#                     Nicolas Borie <Nicolas.Borie at u-pusd.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.magmas import Magmas
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.permutation_groups import PermutationGroups

class FinitePermutationGroups(CategoryWithAxiom):
    r"""
    The category of finite permutation groups, i.e. groups concretely
    represented as groups of permutations acting on a finite set.

    It is currently assumed that any finite permutation group comes
    endowed with a distinguished finite set of generators (method
    ``group_generators``); this is the case for all the existing
    implementations in Sage.

    EXAMPLES::

        sage: C = PermutationGroups().Finite(); C
        Category of finite enumerated permutation groups
        sage: C.super_categories()
        [Category of permutation groups,
         Category of finite groups,
         Category of finite finitely generated semigroups]

        sage: C.example()
        Dihedral group of order 6 as a permutation group

    TESTS::

        sage: C is FinitePermutationGroups()
        True
        sage: TestSuite(C).run()

        sage: G = FinitePermutationGroups().example()
        sage: TestSuite(G).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_inverse() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
    """

    def example(self):
        """
        Returns an example of finite permutation group, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: G = FinitePermutationGroups().example(); G
            Dihedral group of order 6 as a permutation group
        """
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        return DihedralGroup(3)

    def extra_super_categories(self):
        """
        Any permutation group is assumed to be endowed with a finite set of generators.

        TESTS:

            sage: PermutationGroups().Finite().extra_super_categories()
            [Category of finitely generated magmas]
        """
        return [Magmas().FinitelyGenerated()]

    class ParentMethods:
        # TODO
        #  - Port features from MuPAD-Combinat, lib/DOMAINS/CATEGORIES/PermutationGroup.mu
        #  - Move here generic code from sage/groups/perm_gps/permgroup.py

        def cycle_index(self, parent = None):
            r"""
            INPUT:

             - ``self`` - a permutation group `G`
             - ``parent`` -- a free module with basis indexed by partitions,
               or behave as such, with a ``term`` and ``sum`` method
               (default: the symmetric functions over the rational field in the p basis)

            Returns the *cycle index* of `G`, which is a gadget counting
            the elements of `G` by cycle type, averaged over the group:

            .. MATH::

                P = \frac{1}{|G|} \sum_{g\in G} p_{ \operatorname{cycle\ type}(g) }

            EXAMPLES:

            Among the permutations of the symmetric group `S_4`, there is
            the identity, 6 cycles of length 2, 3 products of two cycles
            of length 2, 8 cycles of length 3, and 6 cycles of length 4::

                sage: S4 = SymmetricGroup(4)
                sage: P = S4.cycle_index()
                sage: 24 * P
                p[1, 1, 1, 1] + 6*p[2, 1, 1] + 3*p[2, 2] + 8*p[3, 1] + 6*p[4]

            If `l = (l_1,\dots,l_k)` is a partition, ``|G| P[l]`` is the number
            of elements of `G` with cycles of length `(p_1,\dots,p_k)`::

                sage: 24 * P[ Partition([3,1]) ]
                8

            The cycle index plays an important role in the enumeration of
            objects modulo the action of a group (Polya enumeration), via
            the use of symmetric functions and plethysms. It is therefore
            encoded as a symmetric function, expressed in the powersum
            basis::

                sage: P.parent()
                Symmetric Functions over Rational Field in the powersum basis

            This symmetric function can have some nice properties; for
            example, for the symmetric group `S_n`, we get the complete
            symmetric function `h_n`::

                sage: S = SymmetricFunctions(QQ); h = S.h()
                sage: h( P )
                h[4]

            TODO: add some simple examples of Polya enumeration, once it
            will be easy to expand symmetric functions on any alphabet.

            Here are the cycle indices of some permutation groups::

                sage: 6 * CyclicPermutationGroup(6).cycle_index()
                p[1, 1, 1, 1, 1, 1] + p[2, 2, 2] + 2*p[3, 3] + 2*p[6]

                sage: 60 * AlternatingGroup(5).cycle_index()
                p[1, 1, 1, 1, 1] + 15*p[2, 2, 1] + 20*p[3, 1, 1] + 24*p[5]

                sage: for G in TransitiveGroups(5):               # optional - database_gap # long time
                ....:     G.cardinality() * G.cycle_index()
                p[1, 1, 1, 1, 1] + 4*p[5]
                p[1, 1, 1, 1, 1] + 5*p[2, 2, 1] + 4*p[5]
                p[1, 1, 1, 1, 1] + 5*p[2, 2, 1] + 10*p[4, 1] + 4*p[5]
                p[1, 1, 1, 1, 1] + 15*p[2, 2, 1] + 20*p[3, 1, 1] + 24*p[5]
                p[1, 1, 1, 1, 1] + 10*p[2, 1, 1, 1] + 15*p[2, 2, 1] + 20*p[3, 1, 1] + 20*p[3, 2] + 30*p[4, 1] + 24*p[5]

            One may specify another parent for the result::

                sage: F = CombinatorialFreeModule(QQ, Partitions())
                sage: P = CyclicPermutationGroup(6).cycle_index(parent = F)
                sage: 6 * P
                B[[1, 1, 1, 1, 1, 1]] + B[[2, 2, 2]] + 2*B[[3, 3]] + 2*B[[6]]
                sage: P.parent() is F
                True

            This parent should have a ``term`` and ``sum`` method::

                sage: CyclicPermutationGroup(6).cycle_index(parent = QQ)
                Traceback (most recent call last):
                  ...
                AssertionError: `parent` should be (or behave as) a free module with basis indexed by partitions

            REFERENCES:

            - [Ke1991]_

            AUTHORS:

             - Nicolas Borie and Nicolas M. Thiery

            TESTS::

                sage: P = PermutationGroup([]); P
                Permutation Group with generators [()]
                sage: P.cycle_index()
                p[1]
                sage: P = PermutationGroup([[(1)]]); P
                Permutation Group with generators [()]
                sage: P.cycle_index()
                p[1]
            """
            from sage.combinat.permutation import Permutation
            if parent is None:
                 from sage.rings.rational_field import QQ
                 from sage.combinat.sf.sf import SymmetricFunctions
                 parent = SymmetricFunctions(QQ).powersum()
            else:
                assert hasattr(parent, "term") and hasattr(parent, "sum"), \
                    "`parent` should be (or behave as) a free module with basis indexed by partitions"
            base_ring = parent.base_ring()
            # TODO: use self.conjugacy_classes() once available
            from sage.interfaces.gap import gap
            CC = ([Permutation(self(C.Representative())).cycle_type(), base_ring(C.Size())] for C in gap(self).ConjugacyClasses())
            return parent.sum( parent.term( partition, coeff ) for (partition, coeff) in CC)/self.cardinality()


    class ElementMethods:
        # TODO: put abstract_methods for
        #  - cycle_type
        #  - orbit
        #  - ...
        pass
