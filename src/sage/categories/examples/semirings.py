# sage_setup: distribution = sagemath-categories
r"""
Examples of semirings
"""
# ****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.categories.semirings import Semirings


# semantic : inside NN ;
# 0 => 0, empty;
# 1 => 1, unique;
# 2 => at least 2
_ADD = [[0, 1, 2], [1, 2, 2], [2, 2, 2]]
_PROD = [[0, 0, 0], [0, 1, 2], [0, 2, 2]]


class Ternary(Element):
    """
    Elements of the Ternary-logic ring.

    The semantic is as follows:

    - 0 -- the integer 0
    - 1 -- the integer 1
    - 2 -- some integer greater than 1

    An alternative semantic is:

    - 0 -- an empty set
    - 1 -- a connected set
    - 2 -- a disconnected set

    The same semantic works for graphs instead of sets.
    """
    def __init__(self, parent, n):
        if n not in [0, 1, 2]:
            raise ValueError
        self._n = n
        Element.__init__(self, parent)

    def _repr_(self):
        return ["0", "1", "many"][self._n]

    def __eq__(self, other):
        if not isinstance(other, Ternary):
            return False
        return self._n == other._n

    def __ne__(self, other):
        return not (self == other)


class TernaryLogic(UniqueRepresentation, Parent):
    r"""
    An example of a semiring.

    This class illustrates a minimal implementation of a semiring.

    EXAMPLES::

        sage: S = Semirings().example(); S
        An example of a semiring: the ternary-logic semiring

    This is the semiring that contains 3 objects::

        sage: S.some_elements()
        [0, 1, many]

    The product rule is as expected::

        sage: S(1) * S(1)
        1
        sage: S(1) + S(1)
        many

    TESTS::

        sage: TestSuite(S).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_distributivity() . . . pass
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
        running ._test_zero() . . . pass
    """
    def __init__(self):
        r"""
        The ternary-logic semiring.

        EXAMPLES::

            sage: S = Semirings().example(); S
            An example of a semiring: the ternary-logic semiring
        """
        Parent.__init__(self, category=Semirings())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: Semirings().example()._repr_()
            'An example of a semiring: the ternary-logic semiring'
        """
        return "An example of a semiring: the ternary-logic semiring"

    def __contains__(self, n) -> bool:
        if isinstance(n, Ternary):
            return True
        return n in [0, 1, 2]

    def summation(self, x, y):
        r"""
        Return the sum of ``x`` and ``y`` in the semiring as per
        :meth:`Semirings.ParentMethods.summation`.

        EXAMPLES::

            sage: S = Semirings().example()
            sage: S(1) + S(1)
            many
        """
        assert x in self
        assert y in self
        return self(_ADD[x._n][y._n])

    def one(self):
        """
        Return the unit of ``self``.
        """
        return self(1)

    def product(self, x, y):
        r"""
        Return the product of ``x`` and ``y`` in the semiring as per
        :meth:`Semirings.ParentMethods.product`.

        EXAMPLES::

            sage: S = Semirings().example()
            sage: S(1) * S(2)
            many
        """
        assert x in self
        assert y in self
        return self(_PROD[x._n][y._n])

    def an_element(self):
        r"""
        Return an element of the semiring.

        EXAMPLES::

            sage: Semirings().example().an_element()
            many
        """
        return self(2)

    def some_elements(self):
        r"""
        Return a list of some elements of the semiring.

        EXAMPLES::

            sage: Semirings().example().some_elements()
            [0, 1, many]
        """
        return [self(i) for i in [0, 1, 2]]

    Element = Ternary
