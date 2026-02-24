r"""
Non-decreasing parking functions

A *non-decreasing parking function* of size `n` is a non-decreasing
function `f` from `\{1,\dots,n\}` to itself such that for all `i`, one
has `f(i) \leq i`.

The number of non-decreasing parking functions of size `n` is the `n`-th
:func:`Catalan number<sage.combinat.combinat.catalan_number>`.

The set of non-decreasing parking functions of size `n` is in bijection with
the set of :mod:`Dyck words<sage.combinat.dyck_word>` of size `n`.

AUTHORS:

- Florent Hivert (2009-04)
- Christian Stump (2012-11) added pretty printing
"""
# ****************************************************************************
#       Copyright (C) 2007 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
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
# ****************************************************************************
from __future__ import annotations

from copy import copy

from .combinat import catalan_number
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_with_grading import SetsWithGrading
from sage.categories.monoids import Monoids
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.rings.integer import Integer
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation


def NonDecreasingParkingFunctions(n=None):
    r"""
    Return the set of Non-Decreasing Parking Functions.

    A *non-decreasing parking function* of size `n` is a non-decreasing
    function `f` from `\{1,\dots,n\}` to itself such that for all `i`,
    one has `f(i) \leq i`.

    EXAMPLES:

    Here are all the-non decreasing parking functions of size 5::

        sage: NonDecreasingParkingFunctions(3).list()
        [[1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 2, 2], [1, 2, 3]]

    If no size is specified, then NonDecreasingParkingFunctions
    returns the set of all non-decreasing parking functions.

    ::

        sage: PF = NonDecreasingParkingFunctions(); PF
        Non-decreasing parking functions
        sage: [] in PF
        True
        sage: [1] in PF
        True
        sage: [2] in PF
        False
        sage: [1,1,3] in PF
        True
        sage: [1,1,4] in PF
        False

    If the size `n` is specified, then NonDecreasingParkingFunctions returns
    the set of all non-decreasing parking functions of size `n`.

    ::

        sage: PF = NonDecreasingParkingFunctions(0)
        sage: PF.list()
        [[]]
        sage: PF = NonDecreasingParkingFunctions(1)
        sage: PF.list()
        [[1]]
        sage: PF = NonDecreasingParkingFunctions(3)
        sage: PF.list()
        [[1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 2, 2], [1, 2, 3]]

        sage: PF3 = NonDecreasingParkingFunctions(3); PF3
        Non-decreasing parking functions of size 3
        sage: [] in PF3
        False
        sage: [1] in PF3
        False
        sage: [1,1,3] in PF3
        True
        sage: [1,1,4] in PF3
        False

    TESTS::

        sage: PF = NonDecreasingParkingFunctions(5)
        sage: len(PF.list()) == PF.cardinality()
        True
        sage: NonDecreasingParkingFunctions("foo")
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'foo' to an integer
    """
    if n is None:
        return NonDecreasingParkingFunctions_all()
    else:
        return NonDecreasingParkingFunctions_n(n)


def is_a(x, n=None) -> bool:
    """
    Check whether a list is a non-decreasing parking function.

    If a size `n` is specified, checks if a list is a non-decreasing
    parking function of size `n`.

    TESTS::

        sage: from sage.combinat.non_decreasing_parking_function import is_a
        sage: is_a([1,1,2])
        True
        sage: is_a([1,1,4])
        False
        sage: is_a([1,1,3], 3)
        True
    """
    if not isinstance(x, (list, tuple)):
        return False
    prev = 1
    for i, elt in enumerate(x):
        if prev > elt or elt > i + 1:
            return False
        prev = elt
    return n is None or n == len(x)


class NonDecreasingParkingFunction(Element):
    r"""
    A *non decreasing parking function* of size `n` is a non-decreasing
    function `f` from `\{1,\dots,n\}` to itself such that for all `i`,
    one has `f(i) \leq i`.

    EXAMPLES::

        sage: NonDecreasingParkingFunction([])
        []
        sage: NonDecreasingParkingFunction([1])
        [1]
        sage: NonDecreasingParkingFunction([2])
        Traceback (most recent call last):
        ...
        ValueError: [2] is not a non-decreasing parking function
        sage: NonDecreasingParkingFunction([1,2])
        [1, 2]
        sage: NonDecreasingParkingFunction([1,1,2])
        [1, 1, 2]
        sage: NonDecreasingParkingFunction([1,1,4])
        Traceback (most recent call last):
        ...
        ValueError: [1, 1, 4] is not a non-decreasing parking function
    """

    def __init__(self, lst):
        """
        TESTS::

            sage: NonDecreasingParkingFunction([1, 1, 2, 2, 5, 6])
            [1, 1, 2, 2, 5, 6]
        """
        if not is_a(lst):
            raise ValueError('%s is not a non-decreasing parking function' % lst)
        parent = NonDecreasingParkingFunctions_n(len(lst))
        Element.__init__(self, parent)
        self._list = lst

    def __getitem__(self, n):
        """
        Return the `n`-th item in the underlying list.

        .. NOTE::

           Note that this is different than the image of `n` under
           function.  It is "off by one".

        EXAMPLES::

            sage: p = NonDecreasingParkingFunction([1, 1, 2, 2, 5, 6])
            sage: p[0]
            1
            sage: p[2]
            2
        """
        return self._list[n]

    def __call__(self, n):
        """
        Return the image of ``n`` under the parking function.

        EXAMPLES::

            sage: p = NonDecreasingParkingFunction([1, 1, 2, 2, 5, 6])
            sage: p(3)
            2
            sage: p(6)
            6
        """
        return self._list[n - 1]

    def _mul_(self, lp) -> NonDecreasingParkingFunction:
        """
        The composition of non-decreasing parking functions.

        EXAMPLES::

            sage: p = NonDecreasingParkingFunction([1,1,3])
            sage: q = NonDecreasingParkingFunction([1,2,2])
            sage: p * q
            [1, 1, 1]
            sage: q * p
            [1, 1, 2]
        """
        lp = lp._list
        sp = self._list
        lp = lp[:] + [i + 1 for i in range(len(lp), len(lp))]  # ?
        sp = sp[:] + [i + 1 for i in range(len(sp), len(lp))]
        return NonDecreasingParkingFunction([sp[i - 1] for i in lp])

    def to_dyck_word(self):
        """
        Implement the bijection to :class:`Dyck
        words<sage.combinat.dyck_word.DyckWords>`, which is defined as follows.
        Take a non decreasing parking function, say [1,1,2,4,5,5], and draw
        its graph::

                     ___
                    |  . 5
                   _|  . 5
               ___|  . . 4
             _|  . . . . 2
            |  . . . . . 1
            |  . . . . . 1

        The corresponding Dyck word [1,1,0,1,0,0,1,0,1,1,0,0] is then read off
        from the sequence of horizontal and vertical steps. The converse
        bijection is :meth:`.from_dyck_word`.

        EXAMPLES::

            sage: NonDecreasingParkingFunction([1,1,2,4,5,5]).to_dyck_word()
            [1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0]
            sage: NonDecreasingParkingFunction([]).to_dyck_word()
            []
            sage: NonDecreasingParkingFunction([1,1,1]).to_dyck_word()
            [1, 1, 1, 0, 0, 0]
            sage: NonDecreasingParkingFunction([1,2,3]).to_dyck_word()
            [1, 0, 1, 0, 1, 0]
            sage: NonDecreasingParkingFunction([1,1,3,3,4,6,6]).to_dyck_word()
            [1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]

        TESTS::

            sage: ndpf = NonDecreasingParkingFunctions(5)
            sage: list(ndpf) == [pf.to_dyck_word().to_non_decreasing_parking_function() for pf in ndpf]
            True
        """
        from sage.combinat.dyck_word import CompleteDyckWords_all
        return CompleteDyckWords_all().from_non_decreasing_parking_function(self)

    def __len__(self) -> int:
        """
        Return the length of ``self``.

        EXAMPLES::

            sage: ndpf = NonDecreasingParkingFunctions(5)
            sage: len(ndpf.random_element())
            5
        """
        return len(self._list)

    grade = __len__  # for the category SetsWithGrading

    def _repr_(self) -> str:
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: NonDecreasingParkingFunction([1,1,1])
            [1, 1, 1]
        """
        return str(self._list)

    def _richcmp_(self, other, op) -> bool:
        """
        Compare ``self`` with ``other``.

        EXAMPLES::

            sage: a = NonDecreasingParkingFunction([1,1,1])
            sage: b = NonDecreasingParkingFunction([1,1,2])
            sage: a == b, a < b, b <= a
            (False, True, False)
        """
        return richcmp(self._list, other._list, op)

    def __hash__(self) -> int:
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: a = NonDecreasingParkingFunction([1,1,1])
            sage: b = NonDecreasingParkingFunction([1,1,2])
            sage: hash(a) == hash(b)
            False
        """
        return hash(tuple(self._list))

    @classmethod
    def from_dyck_word(cls, dw) -> NonDecreasingParkingFunction:
        """
        Bijection from :class:`Dyck
        words<sage.combinat.dyck_word.DyckWords>`. It is the inverse of the
        bijection :meth:`.to_dyck_word`. You can find there the mathematical
        definition.

        EXAMPLES::

            sage: NonDecreasingParkingFunction.from_dyck_word([])
            []
            sage: NonDecreasingParkingFunction.from_dyck_word([1,0])
            [1]
            sage: NonDecreasingParkingFunction.from_dyck_word([1,1,0,0])
            [1, 1]
            sage: NonDecreasingParkingFunction.from_dyck_word([1,0,1,0])
            [1, 2]
            sage: NonDecreasingParkingFunction.from_dyck_word([1,0,1,1,0,1,0,0,1,0])
            [1, 2, 2, 3, 5]

        TESTS::

          sage: ndpf = NonDecreasingParkingFunctions(5)
          sage: list(ndpf) == [NonDecreasingParkingFunction.from_dyck_word(pf.to_dyck_word()) for pf in ndpf]
          True
        """
        res = []
        val = 1
        for i in dw:
            if i == 0:
                val += 1
            else:
                res.append(val)
        return cls(res)


class NonDecreasingParkingFunctions_all(UniqueRepresentation, Parent):
    def __init__(self):
        """
        TESTS::

            sage: PF = NonDecreasingParkingFunctions()
            sage: PF == loads(dumps(PF))
            True
        """
        cat = InfiniteEnumeratedSets() & SetsWithGrading()
        Parent.__init__(self, category=cat)

    def __repr__(self) -> str:
        """
        TESTS::

            sage: repr(NonDecreasingParkingFunctions())
            'Non-decreasing parking functions'
        """
        return "Non-decreasing parking functions"

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: [] in NonDecreasingParkingFunctions()
            True
            sage: [1] in NonDecreasingParkingFunctions()
            True
            sage: [2] in NonDecreasingParkingFunctions()
            False
            sage: [1,1,3] in NonDecreasingParkingFunctions()
            True
            sage: [1,1,4] in NonDecreasingParkingFunctions()
            False
        """
        if isinstance(x, NonDecreasingParkingFunction):
            return True
        return is_a(x)

    def __iter__(self):
        """
        An iterator.

        TESTS::

            sage: it = iter(NonDecreasingParkingFunctions()) # indirect doctest
            sage: [next(it) for i in range(8)]
            [[], [1], [1, 1], [1, 2], [1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 2, 2]]
        """
        for n in NN:
            yield from NonDecreasingParkingFunctions_n(n)

    def graded_component(self, n):
        """
        Return the graded component.

        EXAMPLES::

            sage: P = NonDecreasingParkingFunctions()
            sage: P.graded_component(4) == NonDecreasingParkingFunctions(4)
            True
        """
        return NonDecreasingParkingFunctions_n(n)


class NonDecreasingParkingFunctions_n(UniqueRepresentation, Parent):
    r"""
    The combinatorial class of non-decreasing parking functions of
    size `n`.

    A *non-decreasing parking function* of size `n` is a non-decreasing
    function `f` from `\{1,\dots,n\}` to itself such that for all `i`,
    one has `f(i) \leq i`.

    The number of non-decreasing parking functions of size `n` is the
    `n`-th Catalan number.

    EXAMPLES::

        sage: PF = NonDecreasingParkingFunctions(3)
        sage: PF.list()
        [[1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 2, 2], [1, 2, 3]]
        sage: PF = NonDecreasingParkingFunctions(4)
        sage: PF.list()
        [[1, 1, 1, 1], [1, 1, 1, 2], [1, 1, 1, 3], [1, 1, 1, 4], [1, 1, 2, 2], [1, 1, 2, 3], [1, 1, 2, 4], [1, 1, 3, 3], [1, 1, 3, 4], [1, 2, 2, 2], [1, 2, 2, 3], [1, 2, 2, 4], [1, 2, 3, 3], [1, 2, 3, 4]]
        sage: [ NonDecreasingParkingFunctions(i).cardinality() for i in range(10)]
        [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862]

    .. warning::

       The precise order in which the parking function are generated or
       listed is not fixed, and may change in the future.

    AUTHORS:

    - Florent Hivert
    """

    def __init__(self, n):
        """
        TESTS::

            sage: PF = NonDecreasingParkingFunctions(3)
            sage: PF == loads(dumps(PF))
            True
            sage: TestSuite(PF).run(skip='_test_elements')
        """
        n = Integer(n)
        if n < 0:
            raise ValueError('%s is not a nonnegative integer' % n)
        self.n = n
        Parent.__init__(self, category=Monoids().Enumerated().Finite())

    def __repr__(self) -> str:
        """
        TESTS::

            sage: repr(NonDecreasingParkingFunctions(3))
            'Non-decreasing parking functions of size 3'
        """
        return "Non-decreasing parking functions of size %s" % self.n

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: PF3 = NonDecreasingParkingFunctions(3); PF3
            Non-decreasing parking functions of size 3
            sage: [] in PF3
            False
            sage: [1] in PF3
            False
            sage: [1,1,3] in PF3
            True
            sage: [1,1,1] in PF3
            True
            sage: [1,1,4] in PF3
            False
            sage: all(p in PF3 for p in PF3)
            True
        """
        if isinstance(x, NonDecreasingParkingFunction):
            return True
        return is_a(x, self.n)

    def cardinality(self) -> Integer:
        """
        Return the number of non-decreasing parking functions of size
        `n`.

        This number is the `n`-th :func:`Catalan
        number<sage.combinat.combinat.catalan_number>`.

        EXAMPLES::

            sage: PF = NonDecreasingParkingFunctions(0)
            sage: PF.cardinality()
            1
            sage: PF = NonDecreasingParkingFunctions(1)
            sage: PF.cardinality()
            1
            sage: PF = NonDecreasingParkingFunctions(3)
            sage: PF.cardinality()
            5
            sage: PF = NonDecreasingParkingFunctions(5)
            sage: PF.cardinality()
            42
        """
        return catalan_number(self.n)

    def random_element(self) -> NonDecreasingParkingFunction:
        """
        Return a random parking function of the given size.

        EXAMPLES::

            sage: ndpf = NonDecreasingParkingFunctions(5)
            sage: x = ndpf.random_element(); x  # random
            [1, 2, 2, 4, 5]
            sage: x in ndpf
            True
        """
        from sage.combinat.dyck_word import DyckWords
        n = self.n
        dw = DyckWords(n).random_element()
        return NonDecreasingParkingFunction.from_dyck_word(dw)

    def one(self) -> NonDecreasingParkingFunction:
        """
        Return the unit of this monoid.

        This is the non-decreasing parking function [1, 2, ..., n].

        EXAMPLES::

            sage: ndpf = NonDecreasingParkingFunctions(5)
            sage: x = ndpf.random_element(); x  # random
            sage: e = ndpf.one()
            sage: x == e*x == x*e
            True
        """
        return NonDecreasingParkingFunction(list(range(1, self.n + 1)))

    def __iter__(self):
        """
        Return an iterator for non-decreasing parking functions of size `n`.

        .. warning::

           The precise order in which the parking function are
           generated is not fixed, and may change in the future.

        EXAMPLES::

            sage: PF = NonDecreasingParkingFunctions(0)
            sage: [e for e in PF]      # indirect doctest
            [[]]
            sage: PF = NonDecreasingParkingFunctions(1)
            sage: [e for e in PF]      # indirect doctest
            [[1]]
            sage: PF = NonDecreasingParkingFunctions(3)
            sage: [e for e in PF]      # indirect doctest
            [[1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 2, 2], [1, 2, 3]]
            sage: PF = NonDecreasingParkingFunctions(4)
            sage: [e for e in PF]      # indirect doctest
            [[1, 1, 1, 1], [1, 1, 1, 2], [1, 1, 1, 3], [1, 1, 1, 4], [1, 1, 2, 2], [1, 1, 2, 3], [1, 1, 2, 4], [1, 1, 3, 3], [1, 1, 3, 4], [1, 2, 2, 2], [1, 2, 2, 3], [1, 2, 2, 4], [1, 2, 3, 3], [1, 2, 3, 4]]

        TESTS::

            sage: PF = NonDecreasingParkingFunctions(5)
            sage: [e for e in PF] == PF.list()
            True
            sage: PF = NonDecreasingParkingFunctions(6)
            sage: [e for e in PF] == PF.list()
            True

        Complexity: constant amortized time.
        """
        def iterator_rec(n):
            """
            TESTS::

                sage: PF = NonDecreasingParkingFunctions(2)
                sage: [e for e in PF]      # indirect doctest
                [[1, 1], [1, 2]]
            """
            if n == 0:
                yield []
                return
            if n == 1:
                yield [1]
                return
            for res1 in iterator_rec(n - 1):
                for i in range(res1[-1], n + 1):
                    res = copy(res1)
                    res.append(i)
                    yield res
            return
        for res in iterator_rec(self.n):
            yield NonDecreasingParkingFunction(res)

    Element = NonDecreasingParkingFunction
