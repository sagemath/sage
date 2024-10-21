# sage_setup: distribution = sagemath-categories
"""
Positive Integers
"""
# ****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.sets.integer_range import IntegerRangeInfinite
from sage.rings.integer import Integer


class PositiveIntegers(IntegerRangeInfinite):
    r"""
    The enumerated set of positive integers. To fix the ideas,
    we mean `\{1, 2, 3, 4, 5, \dots \}`.

    This class implements the set of positive integers, as an
    enumerated set (see :class:`InfiniteEnumeratedSets
    <sage.categories.infinite_enumerated_sets.InfiniteEnumeratedSets>`).

    This set is an integer range set. The construction is
    therefore done by IntegerRange (see :class:`IntegerRange
    <sage.sets.integer_range.IntegerRange>`).

    EXAMPLES::

        sage: PP = PositiveIntegers()
        sage: PP
        Positive integers
        sage: PP.cardinality()
        +Infinity
        sage: TestSuite(PP).run()
        sage: PP.list()
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot list an infinite set
        sage: it = iter(PP)
        sage: (next(it), next(it), next(it), next(it), next(it))
        (1, 2, 3, 4, 5)
        sage: PP.first()
        1

    TESTS::

        sage: TestSuite(PositiveIntegers()).run()
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: PP = PositiveIntegers()
            sage: PP.category()
            Category of facade infinite enumerated sets
        """
        IntegerRangeInfinite.__init__(self, Integer(1), Integer(1))

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: PositiveIntegers()
            Positive integers
        """
        return "Positive integers"

    def an_element(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: PositiveIntegers().an_element()
            42
        """
        return Integer(42)

    def _sympy_(self):
        r"""
        Return the SymPy set ``Naturals``.

        EXAMPLES::

            sage: PositiveIntegers()._sympy_()                                          # needs sympy
            Naturals
        """
        from sympy import Naturals
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        return Naturals
