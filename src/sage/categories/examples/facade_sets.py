# sage_setup: distribution = sagemath-categories
r"""
Example of facade set
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.sets_cat import Sets
from sage.categories.monoids import Monoids
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet


class PositiveIntegerMonoid(UniqueRepresentation, Parent):
    r"""

    An example of a facade parent: the positive integers viewed as a
    multiplicative monoid

    This class illustrates a minimal implementation of a facade parent
    which models a subset of a set.

    EXAMPLES::

        sage: S = Sets().Facade().example(); S
        An example of facade set: the monoid of positive integers

    TESTS::

        sage: TestSuite(S).run(verbose = True)
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
          running ._test_nonzero_equal() . . . pass
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
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.examples.facade_sets import PositiveIntegerMonoid
            sage: S = PositiveIntegerMonoid(); S
            An example of facade set: the monoid of positive integers

        TESTS::

            sage: TestSuite(S).run()
        """
        Parent.__init__(self, facade=ZZ, category=Monoids())

    def _repr_(self):
        r"""

        EXAMPLES::

            sage: S = Sets().Facade().example()   # indirect doctest
        """
        return "An example of facade set: the monoid of positive integers"

    def _element_constructor_(self, object):
        r"""
        Construction of elements.

        Since ``self`` is a strict subset of the parent it is a facade
        for, it is mandatory to override this method. This method
        indirectly influences membership testing (see
        :meth:`Parent.__contains__`).

        EXAMPLES::

            sage: S = Sets().Facade().example(); S
            An example of facade set: the monoid of positive integers
            sage: S(1)
            1
            sage: S(int(1))
            1
            sage: S(2/1)
            2
            sage: (parent(S(1)) == ZZ, parent(S(int(1))) == ZZ, parent(S(2/1)))
            (True, True, Integer Ring)
            sage: S(1), S(int(1)), S(2/1)
            (1, 1, 2)
            sage: 1 in S
            True
            sage: 2/1 in S, int(1) in S
            (True, True)
            sage: -1 in S, 1/2 in S
            (False, False)
        """
        object = ZZ(object)
        if object > ZZ(0):
            return object
        else:
            raise ValueError("%s should be positive")


class IntegersCompletion(UniqueRepresentation, Parent):
    r"""
    An example of a facade parent: the set of integers completed with
    `+-\infty`

    This class illustrates a minimal implementation of a facade parent
    that models the union of several other parents.

    EXAMPLES::

        sage: S = Sets().Facade().example("union"); S
        An example of a facade set: the integers completed by +-infinity

    TESTS::

        sage: TestSuite(S).run(verbose = True)
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
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.examples.facade_sets import IntegersCompletion
            sage: S = IntegersCompletion(); S
            An example of a facade set: the integers completed by +-infinity

        TESTS::

            sage: TestSuite(S).run()
        """
        # We can't use InfinityRing, because this ring contains 3
        # elements besides +-infinity. We can not use Set either for the
        # moment, because Set([1,2])(1) raises an error
        Parent.__init__(self, facade=(ZZ, FiniteEnumeratedSet([-infinity, +infinity])), category=Sets())

    def _repr_(self):
        r"""

        EXAMPLES::

            sage: S = Sets().Facade().example()   # indirect doctest
        """
        return "An example of a facade set: the integers completed by +-infinity"
