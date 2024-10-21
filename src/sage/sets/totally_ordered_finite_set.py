# sage_setup: distribution = sagemath-categories
"""
Totally Ordered Finite Sets

AUTHORS:

- Stepan Starosta (2012): Initial version
"""
# ****************************************************************************
#  Copyright (C) 2012 Stepan Starosta <stepan.starosta@gmail.com>
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
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, rich_to_bool
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.categories.posets import Posets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets


class TotallyOrderedFiniteSetElement(Element):
    """
    Element of a finite totally ordered set.

    EXAMPLES::

        sage: S = TotallyOrderedFiniteSet([2,7], facade=False)
        sage: x = S(2)
        sage: print(x)
        2
        sage: x.parent()
        {2, 7}
    """
    def __init__(self, parent, data):
        r"""
        TESTS::

            sage: T = TotallyOrderedFiniteSet([3,2,1], facade=False)
            sage: TestSuite(T.an_element()).run()
        """
        Element.__init__(self, parent)
        self.value = data

    def __eq__(self, other):
        r"""
        Equality.

        EXAMPLES::

            sage: A = TotallyOrderedFiniteSet(['gaga',1], facade=False)
            sage: A('gaga') == 'gaga' #indirect doctest
            False
            sage: 'gaga' == A('gaga')
            False
            sage: A('gaga') == A('gaga')
            True
        """
        try:
            same_parent = self.parent() is other.parent()
        except AttributeError:
            return False

        if not same_parent:
            return False
        return other.value == self.value

    def __ne__(self, other):
        r"""
        Non-equality.

        EXAMPLES::

            sage: A = TotallyOrderedFiniteSet(['gaga',1], facade=False)
            sage: A('gaga') != 'gaga' #indirect doctest
            True
        """
        return not (self == other)

    def _richcmp_(self, other, op):
        r"""
        Comparison.

        For ``self`` and ``other`` that have the same parent the method compares
        their rank.

        TESTS::

            sage: A = TotallyOrderedFiniteSet([3,2,7], facade=False)
            sage: A(3) < A(2) and A(3) <= A(2) and A(2) <= A(2)
            True
            sage: A(2) > A(3) and A(2) >= A(3) and A(7) >= A(7)
            True
            sage: A(3) >= A(7) or A(2) > A(2)
            False
            sage: A(7) < A(2) or A(2) <= A(3) or A(2) < A(2)
            False
        """
        if self.value == other.value:
            return rich_to_bool(op, 0)
        return richcmp(self.rank(), other.rank(), op)

    def _repr_(self):
        r"""
        String representation.

        TESTS::

            sage: A = TotallyOrderedFiniteSet(['gaga',1], facade=False)
            sage: repr(A('gaga')) #indirect doctest
            "'gaga'"
        """
        return repr(self.value)

    def __str__(self):
        r"""
        String that represents ``self``.

        EXAMPLES::

            sage: A = TotallyOrderedFiniteSet(['gaga',1], facade=False)
            sage: str(A('gaga')) #indirect doctest
            'gaga'
        """
        return str(self.value)


class TotallyOrderedFiniteSet(FiniteEnumeratedSet):
    """
    Totally ordered finite set.

    This is a finite enumerated set assuming that the elements are
    ordered based upon their rank (i.e. their position in the set).

    INPUT:

    - ``elements`` -- list of elements in the set

    - ``facade`` -- boolean (default: ``True``); if ``True``, a facade is used.
      It should be set to ``False`` if the elements do not inherit from
      :class:`~sage.structure.element.Element` or if you want a funny order. See
      examples for more details.

    .. SEEALSO::

        :class:`FiniteEnumeratedSet`

    EXAMPLES::

        sage: S = TotallyOrderedFiniteSet([1,2,3])
        sage: S
        {1, 2, 3}
        sage: S.cardinality()
        3

    By default, totally ordered finite set behaves as a facade::

        sage: S(1).parent()
        Integer Ring

    It makes comparison fails when it is not the standard order::

        sage: T1 = TotallyOrderedFiniteSet([3,2,5,1])
        sage: T1(3) < T1(1)
        False
        sage: T2 = TotallyOrderedFiniteSet([3, x])                                      # needs sage.symbolic
        sage: T2(3) < T2(x)                                                             # needs sage.symbolic
        3 < x

    To make the above example work, you should set the argument facade to
    ``False`` in the constructor. In that case, the elements of the set have a
    dedicated class::

        sage: A = TotallyOrderedFiniteSet([3,2,0,'a',7,(0,0),1], facade=False)
        sage: A
        {3, 2, 0, 'a', 7, (0, 0), 1}
        sage: x = A.an_element()
        sage: x
        3
        sage: x.parent()
        {3, 2, 0, 'a', 7, (0, 0), 1}
        sage: A(3) < A(2)
        True
        sage: A('a') < A(7)
        True
        sage: A(3) > A(2)
        False
        sage: A(1) < A(3)
        False
        sage: A(3) == A(3)
        True

    But then, the equality comparison is always False with elements outside of
    the set::

        sage: A(1) == 1
        False
        sage: 1 == A(1)
        False
        sage: 'a' == A('a')
        False
        sage: A('a') == 'a'
        False

    Since :issue:`16280`, totally ordered sets support elements that do
    not inherit from :class:`sage.structure.element.Element`, whether
    they are facade or not::

        sage: S = TotallyOrderedFiniteSet(['a','b'])
        sage: S('a')
        'a'
        sage: S = TotallyOrderedFiniteSet(['a','b'], facade = False)
        sage: S('a')
        'a'

    Multiple elements are automatically deleted::

        sage: TotallyOrderedFiniteSet([1,1,2,1,2,2,5,4])
        {1, 2, 5, 4}
    """
    Element = TotallyOrderedFiniteSetElement

    @staticmethod
    def __classcall__(cls, iterable, facade=True):
        """
        Standard trick to expand the iterable upon input, and
        guarantees unique representation, independently of the type of
        the iterable. See ``UniqueRepresentation``.

        TESTS::

            sage: S1 = TotallyOrderedFiniteSet([1, 2, 3])
            sage: S2 = TotallyOrderedFiniteSet((1, 2, 3))
            sage: S3 = TotallyOrderedFiniteSet(range(1,4))
            sage: S1 is S2
            True
            sage: S2 is S3
            True
        """
        elements = []
        seen = set()
        for x in iterable:
            if x not in seen:
                elements.append(x)
                seen.add(x)
        return super(FiniteEnumeratedSet, cls).__classcall__(
            cls,
            tuple(elements),
            facade)

    def __init__(self, elements, facade=True):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(TotallyOrderedFiniteSet([1,2,3])).run()
            sage: TestSuite(TotallyOrderedFiniteSet([1,2,3],facade=False)).run()
            sage: TestSuite(TotallyOrderedFiniteSet([1,3,2],facade=False)).run()
            sage: TestSuite(TotallyOrderedFiniteSet([])).run()
        """
        Parent.__init__(self, facade=facade, category=(Posets(), FiniteEnumeratedSets()))
        self._elements = elements
        if facade:
            self._facade_elements = None
        else:
            self._facade_elements = self._elements
            self._elements = [self.element_class(self, x) for x in elements]

    def _element_constructor_(self, data):
        r"""
        Build an element of that set from ``data``.

        EXAMPLES::

            sage: S1 = TotallyOrderedFiniteSet([1,2,3])
            sage: x = S1(1); x  # indirect doctest
            1
            sage: x.parent()
            Integer Ring


            sage: S2 = TotallyOrderedFiniteSet([3,2,1], facade=False)
            sage: y = S2(1); y  # indirect doctest
            1
            sage: y.parent()
            {3, 2, 1}
            sage: y in S2
            True
            sage: S2(y) is y
            True
        """
        if self._facade_elements is None:
            return FiniteEnumeratedSet._element_constructor_(self, data)

        try:
            i = self._facade_elements.index(data)
        except ValueError:
            raise ValueError("%s not in %s" % (data, self))

        return self._elements[i]

    def le(self, x, y):
        r"""
        Return ``True`` if `x \le y` for the order of ``self``.

        EXAMPLES::

            sage: T = TotallyOrderedFiniteSet([1,3,2], facade=False)
            sage: T1, T3, T2 = T.list()
            sage: T.le(T1,T3)
            True
            sage: T.le(T3,T2)
            True
        """
        try:
            return self._elements.index(x) <= self._elements.index(y)
        except Exception:
            raise ValueError("arguments must be elements of the set")
