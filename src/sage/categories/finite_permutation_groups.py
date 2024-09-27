# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.groups
r"""
Finite Permutation Groups
"""
# ****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#                     Nicolas Borie <Nicolas.Borie at u-pusd.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.magmas import Magmas
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.misc.cachefunc import cached_method


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
        running ._test_construction() . . . pass
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
        Return an example of finite permutation group, as per
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

        def cycle_index(self, parent=None):
            r"""
            Return the *cycle index* of ``self``.

            INPUT:

            - ``self`` -- a permutation group `G`
            - ``parent`` -- a free module with basis indexed by partitions,
              or behave as such, with a ``term`` and ``sum`` method
              (default: the symmetric functions over the rational field in the `p` basis)

            The *cycle index* of a permutation group `G`
            (:wikipedia:`Cycle_index`) is a gadget counting the
            elements of `G` by cycle type, averaged over the group:

            .. MATH::

                P = \frac{1}{|G|} \sum_{g\in G} p_{ \operatorname{cycle\ type}(g) }

            EXAMPLES:

            Among the permutations of the symmetric group `S_4`, there is
            the identity, 6 cycles of length 2, 3 products of two cycles
            of length 2, 8 cycles of length 3, and 6 cycles of length 4::

                sage: S4 = SymmetricGroup(4)
                sage: P = S4.cycle_index()                                              # needs sage.combinat
                sage: 24 * P                                                            # needs sage.combinat
                p[1, 1, 1, 1] + 6*p[2, 1, 1] + 3*p[2, 2] + 8*p[3, 1] + 6*p[4]

            If `l = (l_1,\dots,l_k)` is a partition, ``|G| P[l]`` is the number
            of elements of `G` with cycles of length `(p_1,\dots,p_k)`::

                sage: 24 * P[ Partition([3,1]) ]                                        # needs sage.combinat
                8

            The cycle index plays an important role in the enumeration of
            objects modulo the action of a group (Pólya enumeration), via
            the use of symmetric functions and plethysms. It is therefore
            encoded as a symmetric function, expressed in the powersum
            basis::

                sage: P.parent()                                                        # needs sage.combinat
                Symmetric Functions over Rational Field in the powersum basis

            This symmetric function can have some nice properties; for
            example, for the symmetric group `S_n`, we get the complete
            symmetric function `h_n`::

                sage: S = SymmetricFunctions(QQ); h = S.h()
                sage: h( P )                                                            # needs sage.combinat
                h[4]

            .. TODO::

                Add some simple examples of Pólya enumeration, once
                it will be easy to expand symmetric functions on any
                alphabet.

            Here are the cycle indices of some permutation groups::

                sage: 6 * CyclicPermutationGroup(6).cycle_index()                       # needs sage.combinat
                p[1, 1, 1, 1, 1, 1] + p[2, 2, 2] + 2*p[3, 3] + 2*p[6]

                sage: 60 * AlternatingGroup(5).cycle_index()                            # needs sage.combinat
                p[1, 1, 1, 1, 1] + 15*p[2, 2, 1] + 20*p[3, 1, 1] + 24*p[5]

                sage: for G in TransitiveGroups(5):               # long time
                ....:     G.cardinality() * G.cycle_index()
                p[1, 1, 1, 1, 1] + 4*p[5]
                p[1, 1, 1, 1, 1] + 5*p[2, 2, 1] + 4*p[5]
                p[1, 1, 1, 1, 1] + 5*p[2, 2, 1] + 10*p[4, 1] + 4*p[5]
                p[1, 1, 1, 1, 1] + 15*p[2, 2, 1] + 20*p[3, 1, 1] + 24*p[5]
                p[1, 1, 1, 1, 1] + 10*p[2, 1, 1, 1] + 15*p[2, 2, 1] + 20*p[3, 1, 1] + 20*p[3, 2] + 30*p[4, 1] + 24*p[5]

            Permutation groups with arbitrary domains are supported
            (see :issue:`22765`)::

                sage: G = PermutationGroup([['b','c','a']], domain=['a','b','c'])
                sage: G.cycle_index()                                                   # needs sage.combinat
                1/3*p[1, 1, 1] + 2/3*p[3]

            One may specify another parent for the result::

                sage: # needs sage.combinat
                sage: F = CombinatorialFreeModule(QQ, Partitions())
                sage: P = CyclicPermutationGroup(6).cycle_index(parent=F)
                sage: 6 * P
                B[[1, 1, 1, 1, 1, 1]] + B[[2, 2, 2]] + 2*B[[3, 3]] + 2*B[[6]]
                sage: P.parent() is F
                True

            This parent should be a module with basis indexed by partitions::

                sage: CyclicPermutationGroup(6).cycle_index(parent = QQ)
                Traceback (most recent call last):
                  ...
                ValueError: `parent` should be a module with basis indexed by partitions

            REFERENCES:

            - [Ke1991]_

            AUTHORS:

            - Nicolas Borie and Nicolas M. Thiéry

            TESTS::

                sage: P = PermutationGroup([]); P
                Permutation Group with generators [()]
                sage: P.cycle_index()                                                   # needs sage.combinat
                p[1]
                sage: P = PermutationGroup([[(1)]]); P
                Permutation Group with generators [()]
                sage: P.cycle_index()                                                   # needs sage.combinat
                p[1]
            """
            from sage.categories.modules import Modules
            if parent is None:
                from sage.rings.rational_field import QQ
                from sage.combinat.sf.sf import SymmetricFunctions
                parent = SymmetricFunctions(QQ).powersum()
            elif parent not in Modules.WithBasis:
                raise ValueError("`parent` should be a module with basis indexed by partitions")
            base_ring = parent.base_ring()
            return parent.sum_of_terms([C.an_element().cycle_type(), base_ring(C.cardinality())]
                                       for C in self.conjugacy_classes()
                                      ) / self.cardinality()

        @cached_method
        def profile_series(self, variable='z'):
            r"""
            Return the (finite) generating series of the (finite) profile of the group.

            The profile of a permutation group G is the counting function that
            maps each nonnegative integer n onto the number of orbits of the
            action induced by G on the n-subsets of its domain.
            If f is the profile of G, f(n) is thus the number of orbits of
            n-subsets of G.

            INPUT:

            - ``variable`` -- a variable, or variable name as a string (default: `'z'`)

            OUTPUT:

            - A polynomial in ``variable`` with nonnegative integer coefficients.
              By default, a polynomial in z over ZZ.

            .. SEEALSO::

                - :meth:`profile`

            EXAMPLES::

                sage: # needs sage.combinat
                sage: C8 = CyclicPermutationGroup(8)
                sage: C8.profile_series()
                z^8 + z^7 + 4*z^6 + 7*z^5 + 10*z^4 + 7*z^3 + 4*z^2 + z + 1
                sage: D8 = DihedralGroup(8)
                sage: poly_D8 = D8.profile_series()
                sage: poly_D8
                z^8 + z^7 + 4*z^6 + 5*z^5 + 8*z^4 + 5*z^3 + 4*z^2 + z + 1
                sage: poly_D8.parent()
                Univariate Polynomial Ring in z over Rational Field
                sage: D8.profile_series(variable='y')
                y^8 + y^7 + 4*y^6 + 5*y^5 + 8*y^4 + 5*y^3 + 4*y^2 + y + 1
                sage: u = var('u')
                sage: D8.profile_series(u).parent()
                Symbolic Ring
            """
            from sage.rings.integer_ring import ZZ

            if isinstance(variable, str):
                variable = ZZ[variable].gen()
            cycle_poly = self.cycle_index()
            return cycle_poly.expand(2).subs(x0=1, x1=variable)

        profile_polynomial = profile_series

        def profile(self, n, using_polya=True):
            r"""
            Return the value in `n` of the profile of the group ``self``.

            Optional argument ``using_polya`` allows to change the default method.

            INPUT:

            - ``n`` -- nonnegative integer

            - ``using_polya`` -- boolean (default: ``True``); if ``True``, the
              computation uses Pólya enumeration (and all values of the profile
              are cached, so this should be the method used in case several of
              them are needed); if ``False``, uses the GAP interface to compute
              the orbit.

            OUTPUT:

            - A nonnegative integer that is the number of orbits of `n`-subsets
              under the action induced by ``self`` on the subsets of its domain
              (i.e. the value of the profile of ``self`` in `n`)

            .. SEEALSO::

                - :meth:`profile_series`

            EXAMPLES::

                sage: C6 = CyclicPermutationGroup(6)
                sage: C6.profile(2)                                                     # needs sage.combinat
                3
                sage: C6.profile(3)                                                     # needs sage.combinat
                4
                sage: D8 = DihedralGroup(8)
                sage: D8.profile(4, using_polya=False)
                8
            """

            if using_polya:
                return self.profile_polynomial()[n]
            else:
                from sage.libs.gap.libgap import libgap
                subs_n = libgap.Combinations(list(self.domain()), n)
                return len(libgap.Orbits(self, subs_n, libgap.OnSets))

    class ElementMethods:
        # TODO: put abstract_methods for
        #  - cycle_type
        #  - orbit
        #  - ...
        pass
