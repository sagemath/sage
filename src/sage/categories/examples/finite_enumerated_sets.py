r"""
Examples of finite enumerated sets
"""
# ****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing


class Example(UniqueRepresentation, Parent):
    r"""
    An example of a finite enumerated set: `\{1,2,3\}`

    This class provides a minimal implementation of a finite enumerated set.

    See :class:`~sage.sets.finite_enumerated_set.FiniteEnumeratedSet` for a
    full featured implementation.

    EXAMPLES::

        sage: C = FiniteEnumeratedSets().example()
        sage: C.cardinality()
        3
        sage: C.list()
        [1, 2, 3]
        sage: C.an_element()
        1

    This checks that the different methods of the enumerated set `C`
    return consistent results::

        sage: TestSuite(C).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_nonzero_equal() . . . pass
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
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """

    def __init__(self):
        """
        TESTS::

            sage: C = FiniteEnumeratedSets().example()
            sage: C
            An example of a finite enumerated set: {1,2,3}
            sage: C.category()
            Category of facade finite enumerated sets
            sage: TestSuite(C).run()
        """
        self._set = [Integer(_) for _ in [1, 2, 3]]
        Parent.__init__(self, facade=IntegerRing(),
                        category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: FiniteEnumeratedSets().example() # indirect doctest
            An example of a finite enumerated set: {1,2,3}
        """
        return "An example of a finite enumerated set: {1,2,3}"

    def __contains__(self, o):
        """
        EXAMPLES::

            sage: C = FiniteEnumeratedSets().example()
            sage: 1 in C
            True
            sage: 0 in C
            False
        """
        return o in self._set

    def __iter__(self):
        """
        EXAMPLES::

            sage: list(FiniteEnumeratedSets().example()) # indirect doctest
            [1, 2, 3]

        """
        return iter(self._set)


class IsomorphicObjectOfFiniteEnumeratedSet(UniqueRepresentation, Parent):

    def __init__(self, ambient=Example()):
        """
        TESTS::

            sage: C = FiniteEnumeratedSets().IsomorphicObjects().example()
            sage: C
            The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
            sage: C.category()
            Category of facade isomorphic objects of finite enumerated sets
            sage: TestSuite(C).run()
        """
        self._ambient = ambient
        Parent.__init__(self, facade=IntegerRing(),
                        category=FiniteEnumeratedSets().IsomorphicObjects())

    def ambient(self):
        """
        Returns the ambient space for ``self``, as per
        :meth:`Sets.Subquotients.ParentMethods.ambient()
        <sage.categories.sets_cat.Sets.Subquotients.ParentMethods.ambient>`.

        EXAMPLES::

            sage: C = FiniteEnumeratedSets().IsomorphicObjects().example(); C
            The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
            sage: C.ambient()
            An example of a finite enumerated set: {1,2,3}
        """
        return self._ambient

    def lift(self, x):
        """
        INPUT:
         - ``x`` -- an element of ``self``

        Lifts ``x`` to the ambient space for ``self``, as per
        :meth:`Sets.Subquotients.ParentMethods.lift()
        <sage.categories.sets_cat.Sets.Subquotients.ParentMethods.lift>`.

        EXAMPLES::

            sage: C = FiniteEnumeratedSets().IsomorphicObjects().example(); C
            The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
            sage: C.lift(9)
            3
        """
        return x.sqrt()

    def retract(self, x):
        """
        INPUT:
         - ``x`` -- an element of the ambient space for ``self``

        Retracts ``x`` from the ambient space to ``self``, as per
        :meth:`Sets.Subquotients.ParentMethods.retract()
        <sage.categories.sets_cat.Sets.Subquotients.ParentMethods.retract>`.

        EXAMPLES::

            sage: C = FiniteEnumeratedSets().IsomorphicObjects().example(); C
            The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
            sage: C.retract(3)
            9
        """
        return x ** 2

    def __contains__(self, x):
        """
        Membership testing by checking whether the preimage by the
        bijection is in the ambient space.

        EXAMPLES::

            sage: A = FiniteEnumeratedSets().IsomorphicObjects().example(); A
            The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
            sage: list(A)
            [1, 4, 9]
            sage: 4 in A
            True
            sage: 3 in A
            False
            sage: None in A
            False

        TODO: devise a robust implementation, and move it up to
        ``FiniteEnumeratedSets.IsomorphicObjects``.
        """
        try:
            return self.lift(x) in self.ambient()
        except Exception:
            return False
