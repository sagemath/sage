"""
Integer Range

AUTHORS:

 - Nicolas Borie  (2010-03): First release.
 - Florent Hivert (2010-03): Added a class factory + cardinality method.
 - Vincent Delecroix (2012-02): add methods rank/unrank, make it compliant with
   Python int.
"""
# ****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.infinity import Infinity, MinusInfinity, PlusInfinity


class IntegerRange(UniqueRepresentation, Parent):
    r"""
    The class of :class:`Integer <sage.rings.integer.Integer>` ranges.

    Returns an enumerated set containing an arithmetic progression of integers.

    INPUT:

    - ``begin`` -- integer, Infinity or -Infinity
    - ``end`` -- integer, Infinity or -Infinity
    - ``step`` -- a nonzero integer (default: 1)
    - ``middle_point`` -- integer inside the set (default: ``None``)

    OUTPUT:

    A parent in the category :class:`FiniteEnumeratedSets()
    <sage.categories.finite_enumerated_sets.FiniteEnumeratedSets>` or
    :class:`InfiniteEnumeratedSets()
    <sage.categories.infinite_enumerated_sets.InfiniteEnumeratedSets>`
    depending on the arguments defining ``self``.

    ``IntegerRange(i, j)`` returns the set of `\{i, i+1, i+2, \dots , j-1\}`.
    ``start`` (!) defaults to 0. When ``step`` is given, it specifies the
    increment. The default increment is `1`. IntegerRange allows ``begin`` and
    ``end`` to be infinite.

    ``IntegerRange`` is designed to have similar interface Python
    range. However, whereas ``range`` accept and returns Python ``int``,
    ``IntegerRange`` deals with :class:`Integer <sage.rings.integer.Integer>`.

    If ``middle_point`` is given, then the elements are generated starting
    from it, in a alternating way: `\{m, m+1, m-2, m+2, m-2 \dots \}`.

    EXAMPLES::

        sage: list(IntegerRange(5))
        [0, 1, 2, 3, 4]
        sage: list(IntegerRange(2,5))
        [2, 3, 4]
        sage: I = IntegerRange(2,100,5); I
        {2, 7, ..., 97}
        sage: list(I)
        [2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97]
        sage: I.category()
        Category of facade finite enumerated sets
        sage: I[1].parent()
        Integer Ring

    When ``begin`` and ``end`` are both finite, ``IntegerRange(begin, end,
    step)`` is the set whose list of elements is equivalent to the python
    construction ``range(begin, end, step)``::

        sage: list(IntegerRange(4,105,3)) == list(range(4,105,3))
        True
        sage: list(IntegerRange(-54,13,12)) == list(range(-54,13,12))
        True

    Except for the type of the numbers::

        sage: type(IntegerRange(-54,13,12)[0]), type(list(range(-54,13,12))[0])
        (<... 'sage.rings.integer.Integer'>, <... 'int'>)

    When ``begin`` is finite and ``end`` is +Infinity, ``self`` is the infinite
    arithmetic progression starting from the ``begin`` by step ``step``::

        sage: I = IntegerRange(54,Infinity,3); I
        {54, 57, ...}
        sage: I.category()
        Category of facade infinite enumerated sets
        sage: p = iter(I)
        sage: (next(p), next(p), next(p), next(p), next(p), next(p))
        (54, 57, 60, 63, 66, 69)

        sage: I = IntegerRange(54,-Infinity,-3); I
        {54, 51, ...}
        sage: I.category()
        Category of facade infinite enumerated sets
        sage: p = iter(I)
        sage: (next(p), next(p), next(p), next(p), next(p), next(p))
        (54, 51, 48, 45, 42, 39)

    When ``begin`` and ``end`` are both infinite, you will have to specify the
    extra argument ``middle_point``. ``self`` is then defined by a point
    and a progression/regression setting by ``step``. The enumeration
    is done this way: (let us call `m` the ``middle_point``)
    `\{m, m+step, m-step, m+2step, m-2step, m+3step, \dots \}`::

        sage: I = IntegerRange(-Infinity,Infinity,37,-12); I
        Integer progression containing -12 with increment 37 and bounded with -Infinity and +Infinity
        sage: I.category()
        Category of facade infinite enumerated sets
        sage: -12 in I
        True
        sage: -15 in I
        False
        sage: p = iter(I)
        sage: (next(p), next(p), next(p), next(p), next(p), next(p), next(p), next(p))
        (-12, 25, -49, 62, -86, 99, -123, 136)

    It is also possible to use the argument ``middle_point`` for other cases, finite
    or infinite. The set will be the same as if you didn't give this extra argument
    but the enumeration will begin with this ``middle_point``::

        sage: I = IntegerRange(123,-12,-14); I
        {123, 109, ..., -3}
        sage: list(I)
        [123, 109, 95, 81, 67, 53, 39, 25, 11, -3]
        sage: J = IntegerRange(123,-12,-14,25); J
        Integer progression containing 25 with increment -14 and bounded with 123 and -12
        sage: list(J)
        [25, 11, 39, -3, 53, 67, 81, 95, 109, 123]

    Remember that, like for range, if you define a non empty set, ``begin`` is
    supposed to be included and ``end`` is supposed to be excluded. In the same
    way, when you define a set with a ``middle_point``, the ``begin`` bound will
    be supposed to be included and the ``end`` bound supposed to be excluded::

        sage: I = IntegerRange(-100,100,10,0)
        sage: J = list(range(-100,100,10))
        sage: 100 in I
        False
        sage: 100 in J
        False
        sage: -100 in I
        True
        sage: -100 in J
        True
        sage: list(I)
        [0, 10, -10, 20, -20, 30, -30, 40, -40, 50, -50, 60, -60, 70, -70, 80, -80, 90, -90, -100]


    .. NOTE::

       The input is normalized so that::

        sage: IntegerRange(1, 6, 2) is IntegerRange(1, 7, 2)
        True
        sage: IntegerRange(1, 8, 3) is IntegerRange(1, 10, 3)
        True

    TESTS::

        sage: # Some category automatic tests
        sage: TestSuite(IntegerRange(2,100,3)).run()
        sage: TestSuite(IntegerRange(564,-12,-46)).run()
        sage: TestSuite(IntegerRange(2,Infinity,3)).run()
        sage: TestSuite(IntegerRange(732,-Infinity,-13)).run()
        sage: TestSuite(IntegerRange(-Infinity,Infinity,3,2)).run()
        sage: TestSuite(IntegerRange(56,Infinity,12,80)).run()
        sage: TestSuite(IntegerRange(732,-12,-2743,732)).run()
        sage: # 20 random tests: range and IntegerRange give the same set for finite cases
        sage: for i in range(20):
        ....:     begin = Integer(randint(-300,300))
        ....:     end = Integer(randint(-300,300))
        ....:     step = Integer(randint(-20,20))
        ....:     if step == 0:
        ....:         step = Integer(1)
        ....:     assert list(IntegerRange(begin, end, step)) == list(range(begin, end, step))
        sage: # 20 random tests: range and IntegerRange with middle point for finite cases
        sage: for i in range(20):
        ....:     begin = Integer(randint(-300,300))
        ....:     end = Integer(randint(-300,300))
        ....:     step = Integer(randint(-15,15))
        ....:     if step == 0:
        ....:         step = Integer(-3)
        ....:     I = IntegerRange(begin, end, step)
        ....:     if I.cardinality() == 0:
        ....:         assert len(range(begin, end, step)) == 0
        ....:     else:
        ....:         TestSuite(I).run()
        ....:         L1 = list(IntegerRange(begin, end, step, I.an_element()))
        ....:         L2 = list(range(begin, end, step))
        ....:         L1.sort()
        ....:         L2.sort()
        ....:         assert L1 == L2

    Thanks to :issue:`8543` empty integer range are allowed::

        sage: TestSuite(IntegerRange(0, 5, -1)).run()
    """

    @staticmethod
    def __classcall_private__(cls, begin, end=None, step=Integer(1), middle_point=None):
        """
        TESTS::

            sage: IntegerRange(2,5,0)
            Traceback (most recent call last):
            ...
            ValueError: IntegerRange() step argument must not be zero
            sage: IntegerRange(2) is IntegerRange(0, 2)
            True
            sage: IntegerRange(1.0)                                                     # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: end must be Integer or Infinity, not <... 'sage.rings.real_mpfr.RealLiteral'>
        """
        if isinstance(begin, int):
            begin = Integer(begin)
        if isinstance(end, int):
            end = Integer(end)
        if isinstance(step, int):
            step = Integer(step)

        if end is None:
            end = begin
            begin = Integer(0)
        # check of the arguments
        if not isinstance(begin, (Integer, MinusInfinity, PlusInfinity)):
            raise TypeError("begin must be Integer or Infinity, not %r" % type(begin))
        if not isinstance(end, (Integer, MinusInfinity, PlusInfinity)):
            raise TypeError("end must be Integer or Infinity, not %r" % type(end))
        if not isinstance(step, Integer):
            raise TypeError("step must be Integer, not %r" % type(step))
        if step.is_zero():
            raise ValueError("IntegerRange() step argument must not be zero")

        # If begin and end are infinite, middle_point and step will defined the set.
        if begin == -Infinity and end == Infinity:
            if middle_point is None:
                raise ValueError("Can't iterate over this set, please provide middle_point")

        # If we have a middle point, we go on the special enumeration way...
        if middle_point is not None:
            return IntegerRangeFromMiddle(begin, end, step, middle_point)

        if begin == -Infinity or begin == Infinity:
            raise ValueError("Can't iterate over this set: It is impossible to begin an enumeration with plus/minus Infinity")

        # Check for empty sets
        if step > 0 and begin >= end or step < 0 and begin <= end:
            return IntegerRangeEmpty()

        if end != Infinity and end != -Infinity:
            # Normalize the input
            sgn = 1 if step > 0 else -1
            end = begin + ((end - begin - sgn) // (step) + 1) * step
            return IntegerRangeFinite(begin, end, step)
        else:
            return IntegerRangeInfinite(begin, step)

    def _element_constructor_(self, el):
        """
        TESTS::

            sage: S = IntegerRange(1, 10, 2)
            sage: S(1)      #indirect doctest
            1
            sage: S(0)      #indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: 0 not in {1, 3, 5, 7, 9}
        """
        if el in self:
            if not isinstance(el, Integer):
                return Integer(el)
            return el
        else:
            raise ValueError("%s not in %s" % (el, self))

    element_class = Integer


class IntegerRangeEmpty(IntegerRange, FiniteEnumeratedSet):
    r"""
    A singleton class for empty integer ranges.

    See :class:`IntegerRange` for more details.
    """

    # Needed because FiniteEnumeratedSet.__classcall__ takes an argument.
    @staticmethod
    def __classcall__(cls, *args):
        """
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeEmpty
            sage: I = IntegerRangeEmpty(); I
            {}
            sage: I.category()
            Category of facade finite enumerated sets
            sage: TestSuite(I).run()
            sage: I(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 not in {}
        """
        return FiniteEnumeratedSet.__classcall__(cls, ())


class IntegerRangeFinite(IntegerRange):
    r"""
    The class of finite enumerated sets of integers defined by finite
    arithmetic progressions

    See :class:`IntegerRange` for more details.
    """
    def __init__(self, begin, end, step=Integer(1)):
        r"""
        TESTS::

            sage: I = IntegerRange(123,12,-4)
            sage: I.category()
            Category of facade finite enumerated sets
            sage: TestSuite(I).run()
        """
        self._begin = begin
        self._end = end
        self._step = step
        Parent.__init__(self, facade=IntegerRing(), category=FiniteEnumeratedSets())

    def __contains__(self, elt):
        r"""
        Return ``True`` if ``elt`` is in ``self``.

        EXAMPLES::

            sage: I = IntegerRange(123,12,-4)
            sage: 123 in I
            True
            sage: 127 in I
            False
            sage: 12 in I
            False
            sage: 13 in I
            False
            sage: 14 in I
            False
            sage: 15 in I
            True
            sage: 11 in I
            False
        """
        if not isinstance(elt, Integer):
            try:
                x = Integer(elt)
                if x != elt:
                    return False
                elt = x
            except (ValueError, TypeError):
                return False
        if abs(self._step).divides(Integer(elt)-self._begin):
            return (self._begin <= elt < self._end and self._step > 0) or \
                   (self._begin >= elt > self._end and self._step < 0)
        return False

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: IntegerRange(123,12,-4).cardinality()
            28
            sage: IntegerRange(-57,12,8).cardinality()
            9
            sage: IntegerRange(123,12,4).cardinality()
            0
        """
        return (abs((self._end+self._step-self._begin))-1) // abs(self._step)

    def _repr_(self):
        """
        EXAMPLES::

            sage: IntegerRange(1,2)        #indirect doctest
            {1}
            sage: IntegerRange(1,3)        #indirect doctest
            {1, 2}
            sage: IntegerRange(1,5)        #indirect doctest
            {1, 2, 3, 4}
            sage: IntegerRange(1,6)        #indirect doctest
            {1, ..., 5}
            sage: IntegerRange(123,12,-4)  #indirect doctest
            {123, 119, ..., 15}
            sage: IntegerRange(-57,1,3)    #indirect doctest
            {-57, -54, ..., 0}
        """
        if self.cardinality() < 6:
            return "{" + ", ".join(str(x) for x in self) + "}"
        if self._step == 1:
            return "{%s, ..., %s}" % (self._begin, self._end - self._step)
        return "{%s, %s, ..., %s}" % (self._begin, self._begin + self._step,
                                      self._end - self._step)

    def rank(self, x):
        r"""
        EXAMPLES::

            sage: I = IntegerRange(-57,36,8)
            sage: I.rank(23)
            10
            sage: I.unrank(10)
            23
            sage: I.rank(22)
            Traceback (most recent call last):
            ...
            IndexError: 22 not in self
            sage: I.rank(87)
            Traceback (most recent call last):
            ...
            IndexError: 87 not in self
        """
        if x not in self:
            raise IndexError("%s not in self" % x)
        return Integer((x - self._begin)/self._step)

    def __getitem__(self, i):
        r"""
        Return the i-th element of this integer range.

        EXAMPLES::

            sage: I = IntegerRange(1,13,5)
            sage: I[0], I[1], I[2]
            (1, 6, 11)
            sage: I[3]
            Traceback (most recent call last):
            ...
            IndexError: out of range
            sage: I[-1]
            11
            sage: I[-4]
            Traceback (most recent call last):
            ...
            IndexError: out of range

            sage: I = IntegerRange(13,1,-1)
            sage: l = I.list()
            sage: [I[i] for i in range(I.cardinality())] == l
            True
            sage: l.reverse()
            sage: [I[i] for i in range(-1,-I.cardinality()-1,-1)] == l
            True
        """
        if isinstance(i, slice):
            raise NotImplementedError("not yet")

        if isinstance(i, int):
            i = Integer(i)
        elif not isinstance(i, Integer):
            raise ValueError("argument should be an integer")

        if i < 0:
            if i < -self.cardinality():
                raise IndexError("out of range")
            n = (self._end - self._begin)//(self._step)
            return self._begin + (n+i)*self._step
        else:
            if i >= self.cardinality():
                raise IndexError("out of range")
            return self._begin + i * self._step

    unrank = __getitem__

    def __iter__(self):
        r"""
        Return an iterator over the elements of ``self``.

        EXAMPLES::

            sage: I = IntegerRange(123,12,-4)
            sage: p = iter(I)
            sage: [next(p) for i in range(8)]
            [123, 119, 115, 111, 107, 103, 99, 95]
            sage: I = IntegerRange(-57,12,8)
            sage: p = iter(I)
            sage: [next(p) for i in range(8)]
            [-57, -49, -41, -33, -25, -17, -9, -1]
        """
        n = self._begin
        if self._step > 0:
            while n < self._end:
                yield n
                n += self._step
        else:
            while n > self._end:
                yield n
                n += self._step

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: I = IntegerRange(123,12,-4)
            sage: I.an_element()   #indirect doctest
            115
            sage: I = IntegerRange(-57,12,8)
            sage: I.an_element()   #indirect doctest
            -41
        """
        p = (self._begin + 2*self._step)
        if p in self:
            return p
        else:
            return self._begin


class IntegerRangeInfinite(IntegerRange):
    r""" The class of infinite enumerated sets of integers defined by infinite
    arithmetic progressions.

    See :class:`IntegerRange` for more details.
    """
    def __init__(self, begin, step=Integer(1)):
        r"""
        TESTS::

            sage: I = IntegerRange(-57,Infinity,8)
            sage: I.category()
            Category of facade infinite enumerated sets
            sage: TestSuite(I).run()
        """
        if not isinstance(begin, Integer):
            raise TypeError("begin should be Integer, not %r" % type(begin))
        self._begin = begin
        self._step = step
        Parent.__init__(self, facade=IntegerRing(), category=InfiniteEnumeratedSets())

    def _repr_(self):
        r"""
        TESTS::

            sage: IntegerRange(123,12,-4)            #indirect doctest
            {123, 119, ..., 15}
            sage: IntegerRange(-57,1,3)              #indirect doctest
            {-57, -54, ..., 0}

            sage: IntegerRange(-57,Infinity,8)       #indirect doctest
            {-57, -49, ...}
            sage: IntegerRange(-112,-Infinity,-13)   #indirect doctest
            {-112, -125, ...}
        """
        return "{%s, %s, ...}" % (self._begin, self._begin+self._step)

    def __contains__(self, elt):
        r"""
        Return ``True`` if ``elt`` is in ``self``.

        EXAMPLES::

            sage: I = IntegerRange(-57,Infinity,8)
            sage: -57 in I
            True
            sage: -65 in I
            False
            sage: -49 in I
            True
            sage: 743 in I
            True
        """
        if not isinstance(elt, Integer):
            try:
                elt = Integer(elt)
            except (TypeError, ValueError):
                return False
        if abs(self._step).divides(Integer(elt)-self._begin):
            return (self._step > 0 and elt >= self._begin) or \
                   (self._step < 0 and elt <= self._begin)
        return False

    def rank(self, x):
        r"""
        EXAMPLES::

            sage: I = IntegerRange(-57,Infinity,8)
            sage: I.rank(23)
            10
            sage: I.unrank(10)
            23
            sage: I.rank(22)
            Traceback (most recent call last):
            ...
            IndexError: 22 not in self
        """
        if x not in self:
            raise IndexError("%s not in self" % x)
        return Integer((x - self._begin)/self._step)

    def __getitem__(self, i):
        r"""
        Return the ``i``-th element of ``self``.

        EXAMPLES::

            sage: I = IntegerRange(-8,Infinity,3)
            sage: I.unrank(1)
            -5
        """
        if isinstance(i, slice):
            raise NotImplementedError("not yet")

        if isinstance(i, int):
            i = Integer(i)
        elif not isinstance(i, Integer):
            raise ValueError

        if i < 0:
            raise IndexError("out of range")
        else:
            return self._begin + i * self._step

    unrank = __getitem__

    def __iter__(self):
        r"""
        Return an iterator over the elements of ``self``.

        EXAMPLES::

            sage: I = IntegerRange(-57,Infinity,8)
            sage: p = iter(I)
            sage: [next(p) for i in range(8)]
            [-57, -49, -41, -33, -25, -17, -9, -1]

            sage: I = IntegerRange(-112,-Infinity,-13)
            sage: p = iter(I)
            sage: [next(p) for i in range(8)]
            [-112, -125, -138, -151, -164, -177, -190, -203]
        """
        n = self._begin
        while True:
            yield n
            n += self._step

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: I = IntegerRange(-57,Infinity,8)
            sage: I.an_element()     #indirect doctest
            191

            sage: I = IntegerRange(-112,-Infinity,-13)
            sage: I.an_element()     #indirect doctest
            -515
        """
        return self._begin + 31*self._step


class IntegerRangeFromMiddle(IntegerRange):
    r"""
    The class of finite or infinite enumerated sets defined with
    an inside point, a progression and two limits.

    See :class:`IntegerRange` for more details.
    """
    def __init__(self, begin, end, step=Integer(1), middle_point=Integer(1)):
        r"""
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(-100,100,10,0)
            sage: I.category()
            Category of facade finite enumerated sets
            sage: TestSuite(I).run()
            sage: I = IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)
            sage: I.category()
            Category of facade infinite enumerated sets
            sage: TestSuite(I).run()

            sage: IntegerRange(0, 5, 1, -3)
            Traceback (most recent call last):
            ...
            ValueError: middle_point is not in the interval
        """
        self._begin = begin
        self._end = end
        self._step = step
        self._middle_point = middle_point
        if middle_point not in self:
            raise ValueError("middle_point is not in the interval")

        if (begin != Infinity and begin != -Infinity and
                end != Infinity and end != -Infinity):
            cat = FiniteEnumeratedSets()
        else:
            cat = InfiniteEnumeratedSets()

        Parent.__init__(self, facade=IntegerRing(), category=cat)

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)   #indirect doctest
            Integer progression containing 0 with increment -37 and bounded with +Infinity and -Infinity
            sage: IntegerRangeFromMiddle(-100,100,10,0)              #indirect doctest
            Integer progression containing 0 with increment 10 and bounded with -100 and 100
        """
        return "Integer progression containing %s with increment %s and bounded with %s and %s" % (
            self._middle_point, self._step, self._begin, self._end)

    def __contains__(self, elt):
        r"""
        Return ``True`` if ``elt`` is in ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(-100,100,10,0)
            sage: -110 in I
            False
            sage: -100 in I
            True
            sage: 30 in I
            True
            sage: 90 in I
            True
            sage: 100 in I
            False
        """
        if not isinstance(elt, Integer):
            try:
                elt = Integer(elt)
            except (TypeError, ValueError):
                return False
        if abs(self._step).divides(Integer(elt)-self._middle_point):
            return (self._begin <= elt and elt < self._end) or \
                   (self._begin >= elt and elt > self._end)
        return False

    def next(self, elt):
        r"""
        Return the next element of ``elt`` in ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(-100,100,10,0)
            sage: (I.next(0), I.next(10), I.next(-10), I.next(20), I.next(-100))
            (10, -10, 20, -20, None)
            sage: I = IntegerRangeFromMiddle(-Infinity,Infinity,10,0)
            sage: (I.next(0), I.next(10), I.next(-10), I.next(20), I.next(-100))
            (10, -10, 20, -20, 110)
            sage: I.next(1)
            Traceback (most recent call last):
            ...
            LookupError: 1 not in Integer progression containing 0 with increment 10 and bounded with -Infinity and +Infinity
        """
        if elt not in self:
            raise LookupError('%r not in %r' % (elt, self))
        n = self._middle_point
        if (elt <= n and self._step > 0) or (elt >= n and self._step < 0):
            right = 2*n-elt+self._step
            if right in self:
                return right
            else:
                left = elt-self._step
                if left in self:
                    return left
        else:
            left = 2*n-elt
            if left in self:
                return left
            else:
                right = elt+self._step
                if right in self:
                    return right

    def __iter__(self):
        r"""
        Return an iterator over the elements of ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)
            sage: p = iter(I)
            sage: (next(p), next(p), next(p), next(p), next(p), next(p), next(p), next(p))
            (0, -37, 37, -74, 74, -111, 111, -148)
            sage: I = IntegerRangeFromMiddle(-12,214,10,0)
            sage: p = iter(I)
            sage: (next(p), next(p), next(p), next(p), next(p), next(p), next(p), next(p))
            (0, 10, -10, 20, 30, 40, 50, 60)
        """
        n = self._middle_point
        while n is not None:
            yield n
            n = self.next(n)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

           sage: from sage.sets.integer_range import IntegerRangeFromMiddle
           sage: I = IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)
           sage: I.an_element()    #indirect doctest
           0
           sage: I = IntegerRangeFromMiddle(-12,214,10,0)
           sage: I.an_element()    #indirect doctest
           0
        """
        return self._middle_point
