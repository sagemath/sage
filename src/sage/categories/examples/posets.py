# sage_setup: distribution = sagemath-categories
"""
Examples of posets
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.posets import Posets
from sage.structure.element_wrapper import ElementWrapper
from sage.sets.set import Set, Set_object_enumerated
from sage.sets.positive_integers import PositiveIntegers


class FiniteSetsOrderedByInclusion(UniqueRepresentation, Parent):
    r"""
    An example of a poset: finite sets ordered by inclusion.

    This class provides a minimal implementation of a poset

    EXAMPLES::

        sage: P = Posets().example(); P
        An example of a poset: sets ordered by inclusion

    We conclude by running systematic tests on this poset::

        sage: TestSuite(P).run(verbose = True)
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
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """

    def __init__(self):
        r"""
        EXAMPLES::

            sage: P = Posets().example(); P
            An example of a poset: sets ordered by inclusion
            sage: P.category()
            Category of posets
            sage: type(P)
            <class 'sage.categories.examples.posets.FiniteSetsOrderedByInclusion_with_category'>
            sage: TestSuite(P).run()
        """
        Parent.__init__(self, category=Posets())

    def _repr_(self):
        r"""
        TESTS::

            sage: S = Posets().example()
            sage: S._repr_()
            'An example of a poset: sets ordered by inclusion'
        """
        return "An example of a poset: sets ordered by inclusion"

    def le(self, x, y):
        r"""
        Return whether `x` is a subset of `y`.

        EXAMPLES::

            sage: P = Posets().example()
            sage: P.le( P(Set([1,3])), P(Set([1,2,3])) )
            True
            sage: P.le( P(Set([1,3])), P(Set([1,3])) )
            True
            sage: P.le( P(Set([1,2])), P(Set([1,3])) )
            False
        """
        return x.value.issubset(y.value)

    def an_element(self):
        r"""
        Return an element of this poset.

        EXAMPLES::

            sage: B = Posets().example()
            sage: B.an_element()
            {1, 4, 6}
        """
        return self(Set([1,4,6]))

    class Element(ElementWrapper):

        wrapped_class = Set_object_enumerated


class PositiveIntegersOrderedByDivisibilityFacade(UniqueRepresentation, Parent):
    r"""
    An example of a facade poset: the positive integers ordered by divisibility.

    This class provides a minimal implementation of a facade poset

    EXAMPLES::

        sage: P = Posets().example("facade"); P
        An example of a facade poset: the positive integers ordered by divisibility

        sage: P(5)
        5
        sage: P(0)
        Traceback (most recent call last):
        ...
        ValueError: Can't coerce `0` in any parent `An example of a facade poset: the positive integers ordered by divisibility` is a facade for

        sage: 3 in P
        True
        sage: 0 in P
        False
    """

    element_class = type(Set([]))

    def __init__(self):
        r"""
        EXAMPLES::

            sage: P = Posets().example("facade"); P
            An example of a facade poset: the positive integers ordered by divisibility
            sage: P.category()
            Category of facade posets
            sage: type(P)
            <class 'sage.categories.examples.posets.PositiveIntegersOrderedByDivisibilityFacade_with_category'>
            sage: TestSuite(P).run()
        """
        Parent.__init__(self, facade=(PositiveIntegers(),), category=Posets())

    def _repr_(self):
        r"""
        TESTS::

            sage: S = Posets().example("facade")
            sage: S._repr_()
            'An example of a facade poset: the positive integers ordered by divisibility'
        """
        return "An example of a facade poset: the positive integers ordered by divisibility"

    def le(self, x, y):
        r"""
        Return whether `x` is divisible by `y`.

        EXAMPLES::

            sage: P = Posets().example("facade")
            sage: P.le(3, 6)
            True
            sage: P.le(3, 3)
            True
            sage: P.le(3, 7)
            False
        """
        return x.divides(y)
