"""
Examples of finite monoids
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.all import Monoids
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

class IntegerModMonoid(UniqueRepresentation, Parent):
    r"""
    An example of a finite monoid: the integers mod `n`

    This class illustrates a minimal implementation of a finite monoid.

    EXAMPLES::

        sage: S = FiniteMonoids().example(); S
        An example of a finite multiplicative monoid: the integers modulo 12

        sage: S.category()
        Category of finitely generated finite enumerated monoids

    We conclude by running systematic tests on this monoid::

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
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
    """

    def __init__(self, n = 12):
        r"""
        EXAMPLES::

            sage: M = FiniteMonoids().example(6); M
            An example of a finite multiplicative monoid: the integers modulo 6

        TESTS::

            sage: TestSuite(M).run()

        """
        self.n = n
        Parent.__init__(self, category=Monoids().Finite().FinitelyGenerated())

    def _repr_(self):
        r"""
        TESTS::

            sage: M = FiniteMonoids().example()
            sage: M._repr_()
            'An example of a finite multiplicative monoid: the integers modulo 12'
        """
        return "An example of a finite multiplicative monoid: the integers modulo %s"%self.n

    def semigroup_generators(self):
        r"""

        Returns a set of generators for ``self``, as per
        :meth:`Semigroups.ParentMethods.semigroup_generators`.
        Currently this returns all integers mod `n`, which is of
        course far from optimal!

        EXAMPLES::

            sage: M = FiniteMonoids().example()
            sage: M.semigroup_generators()
            Family (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
        """
        return Family(tuple(self(ZZ(i)) for i in range(self.n)))

    @cached_method
    def one(self):
        r"""
        Return the one of the monoid, as per :meth:`Monoids.ParentMethods.one`.

        EXAMPLES::

            sage: M = FiniteMonoids().example()
            sage: M.one()
            1

        """
        return self(ZZ.one())

    def product(self, x, y):
        r"""
        Return the product of two elements `x` and `y` of the monoid, as
        per :meth:`Semigroups.ParentMethods.product`.

        EXAMPLES::

            sage: M = FiniteMonoids().example()
            sage: M.product(M(3), M(5))
            3
        """
        return self((x.value * y.value) % self.n)

    def an_element(self):
        r"""
        Returns an element of the monoid, as per :meth:`Sets.ParentMethods.an_element`.

        EXAMPLES::

            sage: M = FiniteMonoids().example()
            sage: M.an_element()
            6
        """
        return self(ZZ(42) % self.n)

    class Element (ElementWrapper):
        wrapped_class = Integer

Example = IntegerModMonoid
