# -*- coding: utf-8 -*-
"""
Subsets of the Real Line

This module contains subsets of the real line that can be constructed
as the union of a finite set of open and closed intervals.

EXAMPLES::

    sage: RealSet(0,1)
    (0, 1)
    sage: RealSet((0,1), [2,3])
    (0, 1) + [2, 3]
    sage: RealSet(-oo, oo)
    (-oo, +oo)

Brackets must be balanced in Python, so the naive notation for
half-open intervals does not work::

    sage: RealSet([0,1))
    Traceback (most recent call last):
    ...
    SyntaxError: ...

Instead, you can use the following construction functions::

    sage: RealSet.open_closed(0,1)
    (0, 1]
    sage: RealSet.closed_open(0,1)
    [0, 1)
    sage: RealSet.point(1/2)
    {1/2}
    sage: RealSet.unbounded_below_open(0)
    (-oo, 0)
    sage: RealSet.unbounded_below_closed(0)
    (-oo, 0]
    sage: RealSet.unbounded_above_open(1)
    (1, +oo)
    sage: RealSet.unbounded_above_closed(1)
    [1, +oo)

Relations containing symbols and numeric values or constants::

    sage: RealSet(x != 0)
    (-oo, 0) + (0, +oo)
    sage: RealSet(x == pi)
    {pi}
    sage: RealSet(x < 1/2)
    (-oo, 1/2)
    sage: RealSet(1/2 < x)
    (1/2, +oo)
    sage: RealSet(1.5 <= x)
    [1.50000000000000, +oo)

Note that multiple arguments are combined as union::

    sage: RealSet(x >= 0, x < 1)
    (-oo, +oo)
    sage: RealSet(x >= 0, x > 1)
    [0, +oo)
    sage: RealSet(x >= 0, x > -1)
    (-1, +oo)

AUTHORS:

- Laurent Claessens (2010-12-10): Interval and ContinuousSet, posted
  to sage-devel at
  http://www.mail-archive.com/sage-support@googlegroups.com/msg21326.html.

- Ares Ribo (2011-10-24): Extended the previous work defining the
  class RealSet.

- Jordi Saludes (2011-12-10): Documentation and file reorganization.

- Volker Braun (2013-06-22): Rewrite
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.richcmp import richcmp, richcmp_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.rings.all import ZZ
from sage.rings.real_lazy import LazyFieldElement, RLF
from sage.rings.infinity import infinity, minus_infinity


@richcmp_method
class InternalRealInterval(UniqueRepresentation, Parent):
    """
    A real interval.

    You are not supposed to create :class:`RealInterval` objects
    yourself. Always use :class:`RealSet` instead.

    INPUT:

    - ``lower`` -- real or minus infinity; the lower bound of the
      interval.

    - ``lower_closed`` -- boolean; whether the interval is closed
      at the lower bound

    - ``upper`` -- real or (plus) infinity; the upper bound of the
      interval

    - ``upper_closed`` -- boolean; whether the interval is closed
      at the upper bound

    - ``check`` -- boolean; whether to check the other arguments
      for validity
    """
    def __init__(self, lower, lower_closed, upper, upper_closed, check=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RealSet([0, oo])
            Traceback (most recent call last):
            ...
            ValueError: interval cannot be closed at +oo
        """
        self._lower = lower
        self._upper = upper
        self._lower_closed = lower_closed
        self._upper_closed = upper_closed
        if check:
            if not (isinstance(lower, LazyFieldElement) or lower is minus_infinity):
                raise ValueError('lower bound must be in RLF or minus infinity')
            if not (isinstance(upper, LazyFieldElement) or upper is infinity):
                raise ValueError('upper bound must be in RLF or plus infinity')
            if not isinstance(lower_closed, bool):
                raise ValueError('lower_closed must be boolean')
            if not isinstance(upper_closed, bool):
                raise ValueError('upper_closed must be boolean')
            # comparison of infinity with RLF is broken
            if not(lower is minus_infinity or upper is infinity) and lower > upper:
                raise ValueError('lower/upper bounds are not sorted')
            if (lower_closed and lower == minus_infinity):
                raise ValueError('interval cannot be closed at -oo')
            if (upper_closed and upper == infinity):
                raise ValueError('interval cannot be closed at +oo')
            
    def is_empty(self):
        """
        Return whether the interval is empty
        
        The normalized form of :class:`RealSet` has all intervals
        non-empty, so this method usually returns ``False``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.is_empty()
            False
        """
        return (self._lower == self._upper) and not (self._lower_closed and self._upper_closed)

    def is_point(self):
        """
        Return whether the interval consists of a single point

        OUTPUT:
        
        Boolean.

        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.is_point()
            False
        """
        return (self._lower == self._upper) and self._lower_closed and self._upper_closed

    def lower(self):
        """
        Return the lower bound
        
        OUTPUT:

        The lower bound as it was originally specified.
        
        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.lower()
            0
            sage: I.upper()
            1
        """
        if self._lower is minus_infinity:
            return minus_infinity
        else:
            return self._lower._value
    
    def upper(self):
        """
        Return the upper bound

        OUTPUT:

        The upper bound as it was originally specified.
        
        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.lower()
            0
            sage: I.upper()
            1
        """
        if self._upper is infinity:
            return infinity
        else:
            return self._upper._value
                
    def lower_closed(self):
        """
        Return whether the interval is open at the lower bound
        
        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return self._lower_closed

    def upper_closed(self):
        """
        Return whether the interval is closed at the lower bound

        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return self._upper_closed
                
    def lower_open(self):
        """
        Return whether the interval is closed at the upper bound

        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return not self._lower_closed

    def upper_open(self):
        """
        Return whether the interval is closed at the upper bound

        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return not self._upper_closed
                
    def __richcmp__(self, other, op):
        """
        Intervals are sorted by lower bound, then upper bound

        OUTPUT:

        `-1`, `0`, or `+1` depending on how the intervals compare.
        
        EXAMPLES::

            sage: I1 = RealSet.open_closed(1, 3)[0];  I1
            (1, 3]
            sage: I2 = RealSet.open_closed(0, 5)[0];  I2
            (0, 5]
            sage: I1 > I2
            True
            sage: sorted([I1, I2])
            [(0, 5], (1, 3]]

        TESTS:

        Check if a bug in sorting is fixed (:trac:`17714`)::

            sage: RealSet((0, 1),[1, 1],(1, 2))
            (0, 2)
        """
        x = (self._lower, not self._lower_closed, self._upper, self._upper_closed)
        y = (other._lower, not other._lower_closed, other._upper, other._upper_closed)
        return richcmp(x, y, op)
        
    element_class = LazyFieldElement

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: RealSet.open_closed(0, 1)
            (0, 1]
            sage: RealSet.point(0)
            {0}
        """
        if self.is_point():
            return '{' + str(self.lower()) + '}'
        s =  '[' if self._lower_closed else '('
        if self.lower() is minus_infinity:
            s += '-oo'
        else:
            s += str(self.lower())
        s += ', '
        if self.upper() is infinity:
            s += '+oo'
        else:
            s += str(self.upper())
        s +=  ']' if self._upper_closed else ')'
        return s

    def _sympy_condition_(self, variable):
        """
        Convert to a sympy conditional expression.

        INPUT:

        - ``variable`` -- a symbolic variable

        EXAMPLES::

            sage: RealSet(0, 4)._sympy_condition_(x)
            (0 < x) & (x < 4)
        """
        x = variable
        if self.is_point():
            return (x == self.lower())._sympy_()
        true = (x == 0)._sympy_() | True  # trick to get sympy's True
        if self.lower() is not minus_infinity:
            if self._lower_closed:
                lower_condition = (self.lower() <= x)._sympy_()
            else:
                lower_condition = (self.lower() < x)._sympy_()
        else:
            lower_condition = true
        if self.upper() is not infinity:
            if self._upper_closed:
                upper_condition = (x <= self.upper())._sympy_()
            else:
                upper_condition = (x < self.upper())._sympy_()
        else:
            upper_condition = true
        return lower_condition & upper_condition

    def closure(self):
        """
        Return the closure

        OUTPUT:

        The closure as a new :class:`RealInterval`

        EXAMPLES::

            sage: RealSet.open(0,1)[0].closure()
            [0, 1]
            sage: RealSet.open(-oo,1)[0].closure()
            (-oo, 1]
            sage: RealSet.open(0, oo)[0].closure()
            [0, +oo)
        """
        lower_closed = (self._lower != minus_infinity)
        upper_closed = (self._upper != infinity)
        return InternalRealInterval(self._lower, lower_closed, self._upper, upper_closed)
        
    def interior(self):
        """
        Return the interior

        OUTPUT:

        The interior as a new :class:`RealInterval`

        EXAMPLES::

            sage: RealSet.closed(0, 1)[0].interior()
            (0, 1)
            sage: RealSet.open_closed(-oo, 1)[0].interior()
            (-oo, 1)
            sage: RealSet.closed_open(0, oo)[0].interior()
            (0, +oo)
        """
        return InternalRealInterval(self._lower, False, self._upper, False)
        
    def is_connected(self, other):
        """
        Test whether two intervals are connected
        
        OUTPUT:

        Boolean. Whether the set-theoretic union of the two intervals
        has a single connected component.

        EXAMPLES::

            sage: I1 = RealSet.open(0, 1)[0];  I1
            (0, 1)
            sage: I2 = RealSet.closed(1, 2)[0];  I2
            [1, 2]
            sage: I1.is_connected(I2)
            True
            sage: I1.is_connected(I2.interior())
            False
            sage: I1.closure().is_connected(I2.interior())
            True
            sage: I2.is_connected(I1)
            True
            sage: I2.interior().is_connected(I1)
            False
            sage: I2.closure().is_connected(I1.interior())
            True
            sage: I3 = RealSet.closed(1/2, 3/2)[0]; I3
            [1/2, 3/2]
            sage: I1.is_connected(I3)
            True
            sage: I3.is_connected(I1)
            True
        """
        # self is separated and below other
        if self._upper < other._lower:
            return False
        # self is adjacent and below other 
        if self._upper == other._lower:
            return self._upper_closed or other._lower_closed
        # self is separated and above other
        if other._upper < self._lower:
            return False
        # self is adjacent and above other 
        if other._upper == self._lower:
            return self._lower_closed or other._upper_closed
        # They are not separated
        return True
        
    def convex_hull(self, other):
        """
        Return the convex hull of the two intervals

        OUTPUT:

        The convex hull as a new :class:`RealInterval`.

        EXAMPLES::

            sage: I1 = RealSet.open(0, 1)[0];  I1
            (0, 1)
            sage: I2 = RealSet.closed(1, 2)[0];  I2
            [1, 2]
            sage: I1.convex_hull(I2)
            (0, 2]
            sage: I2.convex_hull(I1)
            (0, 2]
            sage: I1.convex_hull(I2.interior())
            (0, 2)
            sage: I1.closure().convex_hull(I2.interior())
            [0, 2)
            sage: I1.closure().convex_hull(I2)
            [0, 2]
            sage: I3 = RealSet.closed(1/2, 3/2)[0]; I3
            [1/2, 3/2]
            sage: I1.convex_hull(I3)
            (0, 3/2]
        """
        if self._lower < other._lower:
            lower = self._lower
            lower_closed = self._lower_closed
        elif self._lower > other._lower:
            lower = other._lower
            lower_closed = other._lower_closed
        else:
            lower = self._lower
            lower_closed = self._lower_closed or other._lower_closed
        if self._upper > other._upper:
            upper = self._upper
            upper_closed = self._upper_closed
        elif self._upper < other._upper:
            upper = other._upper
            upper_closed = other._upper_closed
        else:
            upper = self._upper
            upper_closed = self._upper_closed or other._upper_closed
        return InternalRealInterval(lower, lower_closed, upper, upper_closed)

    def intersection(self, other):
        """
        Return the intersection of the two intervals
        
        INPUT:
        
        - ``other`` -- a :class:`RealInterval`

        OUTPUT:

        The intersection as a new :class:`RealInterval`

        EXAMPLES::

            sage: I1 = RealSet.open(0, 2)[0];  I1
            (0, 2)
            sage: I2 = RealSet.closed(1, 3)[0];  I2
            [1, 3]
            sage: I1.intersection(I2)
            [1, 2)
            sage: I2.intersection(I1)
            [1, 2)
            sage: I1.closure().intersection(I2.interior())
            (1, 2]
            sage: I2.interior().intersection(I1.closure())
            (1, 2]

            sage: I3 = RealSet.closed(10, 11)[0];  I3
            [10, 11]
            sage: I1.intersection(I3)
            (0, 0)
            sage: I3.intersection(I1)
            (0, 0)
        """
        lower = upper = None
        lower_closed = upper_closed = None
        if self._lower < other._lower:
            lower = other._lower
            lower_closed = other._lower_closed
        elif self._lower > other._lower:
            lower = self._lower
            lower_closed = self._lower_closed
        else:
            lower = self._lower
            lower_closed = self._lower_closed and other._lower_closed
        if self._upper > other._upper:
            upper = other._upper
            upper_closed = other._upper_closed
        elif self._upper < other._upper:
            upper = self._upper
            upper_closed = self._upper_closed
        else:
            upper = self._upper
            upper_closed = self._upper_closed and other._upper_closed
        if lower > upper:
            lower = upper = RLF(0)
            lower_closed = upper_closed = False
        return InternalRealInterval(lower, lower_closed, upper, upper_closed)

    def contains(self, x):
        """
        Return whether `x` is contained in the interval

        INPUT:

        - ``x`` -- a real number.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: i = RealSet.open_closed(0,2)[0]; i
            (0, 2]
            sage: i.contains(0)
            False
            sage: i.contains(1)
            True
            sage: i.contains(2)
            True
        """
        if self._lower < x < self._upper:
            return True
        if self._lower == x:
            return self._lower_closed
        if self._upper == x:
            return self._upper_closed
        return False

    def __mul__(self, right):
        r"""
        Scale an interval by a scalar on the left or right.

        If scaled with a negative number, the interval is flipped.

        EXAMPLES::

            sage: i = RealSet.open_closed(0,2)[0]; i
            (0, 2]
            sage: 2 * i
            (0, 4]
            sage: 0 * i
            {0}
            sage: (-2) * i
            [-4, 0)
            sage: i * (-3)
            [-6, 0)
            sage: i * 0
            {0}
            sage: i * 1
            (0, 2]

        TESTS::

            sage: from sage.sets.real_set import InternalRealInterval
            sage: i = InternalRealInterval(RLF(0), False, RLF(0), False)
            sage: (0 * i).is_empty()
            True
        """
        if not isinstance(right, InternalRealInterval):
            right = RLF(right)
            if self.is_empty():
                return self
            lower = self._lower * right
            lower_closed = self._lower_closed
            upper = self._upper * right
            upper_closed = self._upper_closed
            scalar = right
        elif not isinstance(self, InternalRealInterval):
            self = RLF(self)
            if right.is_empty():
                return right
            lower = self * right._lower
            lower_closed = right._lower_closed
            upper = self * right._upper
            upper_closed = right._upper_closed
            scalar = self
        else:
            return NotImplemented
        if scalar == RLF(0):
            return InternalRealInterval(RLF(0), True, RLF(0), True)
        elif scalar < RLF(0):
            lower, lower_closed, upper, upper_closed = upper, upper_closed, lower, lower_closed
        if lower == -infinity:
            lower = -infinity
        if upper == infinity:
            upper = infinity
        return InternalRealInterval(lower, lower_closed,
                                    upper, upper_closed)

    def __rmul__(self, other):
        r"""
        Scale an interval by a scalar on the left.

        If scaled with a negative number, the interval is flipped.

        EXAMPLES::

            sage: i = RealSet.open_closed(0,2)[0]; i
            (0, 2]
            sage: 2 * i
            (0, 4]
            sage: 0 * i
            {0}
            sage: (-2) * i
            [-4, 0)
        """
        return self * other

@richcmp_method
class RealSet(UniqueRepresentation, Parent):

    @staticmethod
    def __classcall__(cls, *args):
        """
        Normalize the input.

        INPUT:
        
        See :class:`RealSet`.
          
        OUTPUT:

        A :class:`RealSet`.

        EXAMPLES::

            sage: R = RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3)); R
            (0, 1] + [2, 3)

        ::

            sage: RealSet(x != 0)
            (-oo, 0) + (0, +oo)
            sage: RealSet(x == pi)
            {pi}
            sage: RealSet(x < 1/2)
            (-oo, 1/2)
            sage: RealSet(1/2 < x)
            (1/2, +oo)
            sage: RealSet(1.5 <= x)
            [1.50000000000000, +oo)
            sage: RealSet(x >= -1)
            [-1, +oo)
            sage: RealSet(x > oo)
            {}
            sage: RealSet(x >= oo)
            {}
            sage: RealSet(x <= -oo)
            {}
            sage: RealSet(x < oo)
            (-oo, +oo)
            sage: RealSet(x > -oo)
            (-oo, +oo)
            sage: RealSet(x != oo)
            (-oo, +oo)
            sage: RealSet(x <= oo)
            Traceback (most recent call last):
            ...
            ValueError: interval cannot be closed at +oo
            sage: RealSet(x == oo)
            Traceback (most recent call last):
            ...
            ValueError: interval cannot be closed at +oo
            sage: RealSet(x >= -oo)
            Traceback (most recent call last):
            ...
            ValueError: interval cannot be closed at -oo

        TESTS::

            sage: TestSuite(R).run()
        """
        from sage.symbolic.expression import Expression
        if len(args) == 1 and isinstance(args[0], RealSet):
            return args[0]   # common optimization
        intervals = []
        if len(args) == 2:
            # allow RealSet(0,1) interval constructor
            try:
                lower, upper = args
                lower.n()
                upper.n()
                args = (RealSet._prep(lower, upper), )
            except (AttributeError, ValueError, TypeError):
                pass
        for arg in args:
            if isinstance(arg, tuple):
                lower, upper = RealSet._prep(*arg)
                intervals.append(InternalRealInterval(lower, False, upper, False))
            elif isinstance(arg, list):
                lower, upper = RealSet._prep(*arg)
                intervals.append(InternalRealInterval(lower, True, upper, True))
            elif isinstance(arg, InternalRealInterval):
                intervals.append(arg)
            elif isinstance(arg, RealSet):
                intervals.extend(arg._intervals)
            elif isinstance(arg, Expression) and arg.is_relational():
                from operator import eq, ne, lt, gt, le, ge
                def rel_to_interval(op, val):
                    """
                    Internal helper function.
                    """
                    oo = infinity
                    try:
                        val = val.pyobject()
                    except AttributeError:
                        pass
                    val = RLF(val)
                    if op == eq:
                        return [InternalRealInterval(val, True, val, True)]
                    elif op == ne:
                        return [InternalRealInterval(-oo, False, val, False),
                                InternalRealInterval(val, False, oo, False)]
                    elif op == gt:
                        return [InternalRealInterval(val, False, oo, False)]
                    elif op == ge:
                        return [InternalRealInterval(val, True, oo, False)]
                    elif op == lt:
                        return [InternalRealInterval(-oo, False, val, False)]
                    else:
                        return [InternalRealInterval(-oo, False, val, True)]

                if (arg.lhs().is_symbol()
                    and (arg.rhs().is_numeric() or arg.rhs().is_constant())
                    and arg.rhs().is_real()):
                    intervals.extend(rel_to_interval(arg.operator(), arg.rhs()))
                elif (arg.rhs().is_symbol()
                    and (arg.lhs().is_numeric() or arg.lhs().is_constant())
                    and arg.lhs().is_real()):
                    op = arg.operator()
                    if op == lt:
                        op = gt
                    elif op == gt:
                        op = lt
                    elif op == le:
                        op = ge
                    elif op == ge:
                        op = le
                    intervals.extend(rel_to_interval(op, arg.lhs()))
                else:
                    raise ValueError(str(arg) + ' does not determine real interval')
            else:
                raise ValueError(str(arg) + ' does not determine real interval')
        intervals = RealSet.normalize(intervals)
        return UniqueRepresentation.__classcall__(cls, *intervals)
                
    def __init__(self, *intervals):
        """
        A subset of the real line

        INPUT:

        Arguments defining a real set. Possibilities are either two
        real numbers to construct an open set or a list/tuple/iterable
        of intervals. The individual intervals can be specified by
        either a :class:`RealInterval`, a tuple of two real numbers
        (constructing an open interval), or a list of two number
        (constructing a closed interval).

        EXAMPLES::

            sage: RealSet(0,1)    # open set from two numbers
            (0, 1)
            sage: i = RealSet(0,1)[0]
            sage: RealSet(i)      # interval
            (0, 1)
            sage: RealSet(i, (3,4))    # tuple of two numbers = open set
            (0, 1) + (3, 4)
            sage: RealSet(i, [3,4])    # list of two numbers = closed set
            (0, 1) + [3, 4]
        """
        Parent.__init__(self, category = Sets())
        self._intervals = intervals
    
    def __richcmp__(self, other, op):
        """
        Intervals are sorted by lower bound, then upper bound

        OUTPUT:

        `-1`, `0`, or `+1` depending on how the intervals compare.
        
        EXAMPLES::

             sage: I1 = RealSet.open_closed(1, 3);  I1
             (1, 3]
             sage: I2 = RealSet.open_closed(0, 5);  I2
             (0, 5]
             sage: I1 > I2
             True
             sage: sorted([I1, I2])
             [(0, 5], (1, 3]]
             sage: I1 == I1
             True
        """
        if not isinstance(other, RealSet):
            return NotImplemented
        # note that the interval representation is normalized into a
        # unique form
        return richcmp(self._intervals, other._intervals, op)

    def __iter__(self):
        """
        Iterate over the component intervals is ascending order
        
        OUTPUT:

        An iterator over the intervals.

        EXAMPLES::

            sage: s = RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3))
            sage: i = iter(s)
            sage: next(i)
            (0, 1]
            sage: next(i)
            [2, 3)
        """
        return iter(self._intervals)

    def n_components(self):
        """
        Return the number of connected components

        See also :meth:`get_interval`

        EXAMPLES::

            sage: s = RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3))
            sage: s.n_components()
            2
        """
        return len(self._intervals)

    def cardinality(self):
        """
        Return the cardinality of the subset of the real line.

        OUTPUT:
        
        Integer or infinity. The size of a discrete set is the number
        of points; the size of a real interval is Infinity.

        EXAMPLES::

           sage: RealSet([0, 0], [1, 1], [3, 3]).cardinality()
           3
           sage: RealSet(0,3).cardinality()
           +Infinity
        """
        n = ZZ(0)
        for interval in self._intervals:
            if interval.is_point():
                n += 1
            else:
                return infinity
        return n

    def is_empty(self):
        """
        Return whether the set is empty
        
        EXAMPLES::

            sage: RealSet(0, 1).is_empty()
            False
            sage: RealSet(0, 0).is_empty()
            True
        """
        return len(self._intervals) == 0

    def get_interval(self, i):
        """
        Return the ``i``-th connected component.

        Note that the intervals representing the real set are always
        normalized, see :meth:`normalize`.

        INPUT:
        
        - ``i`` -- integer.

        OUTPUT:

        The $i$-th connected component as a :class:`RealInterval`.

        EXAMPLES::

            sage: s = RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3))
            sage: s.get_interval(0)
            (0, 1]
            sage: s[0]    # shorthand
            (0, 1]
            sage: s.get_interval(1)
            [2, 3)
            sage: s[0] == s.get_interval(0)
            True
        """
        return self._intervals[i]

    __getitem__ = get_interval

    @staticmethod
    def normalize(intervals):
        """
        Bring a collection of intervals into canonical form

        INPUT:

        - ``intervals`` -- a list/tuple/iterable of intervals.

        OUTPUT:
        
        A tuple of intervals such that

        * they are sorted in ascending order (by lower bound)
        
        * there is a gap between each interval
        
        * all intervals are non-empty

        EXAMPLES::

            sage: i1 = RealSet((0, 1))[0]
            sage: i2 = RealSet([1, 2])[0]
            sage: i3 = RealSet((2, 3))[0]
            sage: RealSet.normalize([i1, i2, i3])
            ((0, 3),)

            sage: RealSet((0, 1), [1, 2], (2, 3))
            (0, 3)
            sage: RealSet((0, 1), (1, 2), (2, 3))
            (0, 1) + (1, 2) + (2, 3)
            sage: RealSet([0, 1], [2, 3])
            [0, 1] + [2, 3]
            sage: RealSet((0, 2), (1, 3))
            (0, 3)
            sage: RealSet(0,0)
            {}
        """
        # sort by lower bound
        intervals = sorted(intervals)
        if not intervals:
            return tuple()
        merged = []
        curr = intervals.pop(0)
        while intervals:
            next = intervals.pop(0)
            if curr._upper > next._lower or (
                    curr._upper == next._lower and
                    (curr._upper_closed or next._lower_closed)):
                curr = curr.convex_hull(next)
            else:
                if not curr.is_empty():
                    merged.append(curr)
                curr = next
        if not curr.is_empty():
            merged.append(curr)
        return tuple(merged)

    def _repr_(self):
        """
        Return a string representation
        
        OUTPUT:

        A string representation.

        EXAMPLES::

            sage: RealSet(0, 1)._repr_()
            '(0, 1)'
        """
        if self.n_components() == 0:
            return '{}'
        else:
            # Switch to u'\u222A' (cup sign) with Python 3
            return ' + '.join(map(repr, self._intervals))
            # return u' ∪ '.join(map(repr, self._intervals)) # py3 only


    def _sympy_condition_(self, variable):
        """
        Convert to a sympy conditional expression.

        INPUT:

        - ``variable`` -- a symbolic variable

        EXAMPLES::

            sage: RealSet(0, 1)._sympy_condition_(x)
            (0 < x) & (x < 1)
            sage: RealSet((0,1), [2,3])._sympy_condition_(x)
            ((2 <= x) & (x <= 3)) | ((0 < x) & (x < 1))
            sage: RealSet.unbounded_below_open(0)._sympy_condition_(x)
            x < 0
            sage: RealSet.unbounded_above_closed(2)._sympy_condition_(x)
            2 <= x

        TESTS::

            sage: RealSet(6,6)._sympy_condition_(x)
            False
            sage: RealSet([6,6])._sympy_condition_(x)
            Eq(x, 6)
        """
        x = variable
        false = (x == 0)._sympy_() & False  # trick to get sympy's False
        if self.n_components() == 0:
            return false
        else:
            cond = false
            for it in self._intervals:
                cond = cond | it._sympy_condition_(x)
            return cond

    @staticmethod
    def _prep(lower, upper=None):
        """
        Helper to prepare the lower and upper bound

        EXAMPLES::

            sage: RealSet._prep(1, 0)
            (0, 1)
            sage: RealSet._prep(oo)
            +Infinity
        """
        if lower == minus_infinity:
            lower = minus_infinity
        if lower == infinity:
            lower = infinity
        else:
            lower = RLF(lower)
        if upper is None:
            return lower
        if upper == minus_infinity:
            upper = minus_infinity
        if upper == infinity:
            upper = infinity
        else:
            upper = RLF(upper)
        if upper is infinity or lower is minus_infinity:
            return lower, upper
        elif lower is infinity or upper is minus_infinity:
            return upper, lower
        elif upper < lower:
            return upper, lower
        else:
            return lower, upper

    @staticmethod
    def open(lower, upper):
        """
        Construct an open interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet`.

        EXAMPLES::

            sage: RealSet.open(1, 0)
            (0, 1)
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(InternalRealInterval(lower, False, upper, False))

    @staticmethod
    def closed(lower, upper):
        """
        Construct a closed interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet`.

        EXAMPLES::

            sage: RealSet.closed(1, 0)
            [0, 1]
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(InternalRealInterval(lower, True, upper, True))

    @staticmethod
    def point(p):
        """
        Construct an interval containing a single point

        INPUT:

        - ``p`` -- a real number.

        OUTPUT:

        A new :class:`RealSet`.

        EXAMPLES::

            sage: RealSet.open(1, 0)
            (0, 1)
        """
        p = RealSet._prep(p)
        return RealSet(InternalRealInterval(p, True, p, True))
    
    @staticmethod
    def open_closed(lower, upper):
        """
        Construct a half-open interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet` that is open at the lower bound and
        closed at the upper bound.

        EXAMPLES::

            sage: RealSet.open_closed(1, 0)
            (0, 1]
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(InternalRealInterval(lower, False, upper, True))

    @staticmethod
    def closed_open(lower, upper):
        """
        Construct an half-open interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet` that is closed at the lower bound and
        open an the upper bound.

        EXAMPLES::

            sage: RealSet.closed_open(1, 0)
            [0, 1)
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(InternalRealInterval(lower, True, upper, False))

    @staticmethod
    def unbounded_below_closed(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from minus infinity to the bound (including).

        EXAMPLES::

            sage: RealSet.unbounded_below_closed(1)
            (-oo, 1]
        """
        bound = RealSet._prep(bound)
        return RealSet(InternalRealInterval(minus_infinity, False, bound, True))

    @staticmethod
    def unbounded_below_open(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from minus infinity to the bound (excluding).

        EXAMPLES::

            sage: RealSet.unbounded_below_open(1)
            (-oo, 1)
        """
        bound = RealSet._prep(bound)
        return RealSet(InternalRealInterval(RLF(minus_infinity), False, RLF(bound), False))

    @staticmethod
    def unbounded_above_closed(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from the bound (including) to plus
        infinity.

        EXAMPLES::

            sage: RealSet.unbounded_above_closed(1)
            [1, +oo)
        """
        bound = RealSet._prep(bound)
        return RealSet(InternalRealInterval(RLF(bound), True, RLF(infinity), False))

    @staticmethod
    def unbounded_above_open(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from the bound (excluding) to plus
        infinity.

        EXAMPLES::

            sage: RealSet.unbounded_above_open(1)
            (1, +oo)
        """
        bound = RealSet._prep(bound)
        return RealSet(InternalRealInterval(RLF(bound), False, RLF(infinity), False))

    def union(self, *other):
        """
        Return the union of the two sets

        INPUT:
        
        - ``other`` -- a :class:`RealSet` or data that defines one.

        OUTPUT:
        
        The set-theoretic union as a new :class:`RealSet`.

        EXAMPLES::

            sage: s1 = RealSet(0,2)
            sage: s2 = RealSet(1,3)
            sage: s1.union(s2)
            (0, 3)
            sage: s1.union(1,3)
            (0, 3)
            sage: s1 | s2    # syntactic sugar
            (0, 3)
            sage: s1 + s2    # syntactic sugar
            (0, 3)
        """
        other = RealSet(*other)
        intervals = self._intervals + other._intervals
        return RealSet(*intervals)
    
    __or__ = union
    __add__ = union

    def intersection(self, *other):
        """
        Return the intersection of the two sets

        INPUT:
        
        - ``other`` -- a :class:`RealSet` or data that defines one.

        OUTPUT:
        
        The set-theoretic intersection as a new :class:`RealSet`.

        EXAMPLES::

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +oo)
            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-oo, -10] + (1, 3)
            sage: s1.intersection(s2)
            (1, 2)
            sage: s1 & s2    # syntactic sugar
            (1, 2)

            sage: s1 = RealSet((0, 1), (2, 3));  s1
            (0, 1) + (2, 3)
            sage: s2 = RealSet([0, 1], [2, 3]);  s2
            [0, 1] + [2, 3]
            sage: s3 = RealSet([1, 2]);  s3
            [1, 2]
            sage: s1.intersection(s2)
            (0, 1) + (2, 3)
            sage: s1.intersection(s3)
            {}
            sage: s2.intersection(s3)
            {1} + {2}
        """
        other = RealSet(*other)
        # TODO: this can be done in linear time since the intervals are already sorted
        intervals = []
        for i1 in self._intervals:
            for i2 in other._intervals:
                intervals.append(i1.intersection(i2))
        return RealSet(*intervals)

    __and__ = intersection

    def inf(self):
        """
        Return the infimum

        OUTPUT:

        A real number or infinity.

        EXAMPLES::

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +oo)
            sage: s1.inf()
            0

            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-oo, -10] + (1, 3)
            sage: s2.inf()
            -Infinity
        """
        if self.n_components() == 0:
            return infinity
        return self._intervals[0].lower()

    def sup(self):
        """
        Return the supremum

        OUTPUT:

        A real number or infinity.

        EXAMPLES::

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +oo)
            sage: s1.sup()
            +Infinity

            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-oo, -10] + (1, 3)
            sage: s2.sup()
            3
        """
        if self.n_components() == 0:
            return minus_infinity
        return self._intervals[-1].upper()

    def complement(self):
        """
        Return the complement

        OUTPUT:

        The set-theoretic complement as a new :class:`RealSet`.
        
        EXAMPLES::

            sage: RealSet(0,1).complement()
            (-oo, 0] + [1, +oo)
       
            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +oo)
            sage: s1.complement()
            (-oo, 0] + [2, 10)

            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-oo, -10] + (1, 3)
            sage: s2.complement()
            (-10, 1] + [3, +oo)
        """
        n = self.n_components()
        if n == 0:
            return RealSet(minus_infinity, infinity)
        intervals = []
        if self.inf() != minus_infinity:
            first = self._intervals[0]
            intervals.append(InternalRealInterval(RLF(minus_infinity), False,
                                          first._lower, first.lower_open()))
        if self.sup() != infinity:
            last = self._intervals[-1]
            intervals.append(InternalRealInterval(last._upper, last.upper_open(),
                                          RLF(infinity), False))
        for i in range(1,n):
            prev = self._intervals[i-1]
            next = self._intervals[i]
            i = InternalRealInterval(prev._upper, prev.upper_open(),
                             next._lower, next.lower_open())
            intervals.append(i)
        return RealSet(*intervals)
                             
    def difference(self, *other):
        """
        Return ``self`` with ``other`` subtracted

        INPUT:
        
        - ``other`` -- a :class:`RealSet` or data that defines one.

        OUTPUT:
        
        The set-theoretic difference of ``self`` with ``other``
        removed as a new :class:`RealSet`.

        EXAMPLES::

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +oo)
            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-oo, -10] + (1, 3)
            sage: s1.difference(s2)
            (0, 1] + [10, +oo)
            sage: s1 - s2    # syntactic sugar
            (0, 1] + [10, +oo)
            sage: s2.difference(s1)
            (-oo, -10] + [2, 3)
            sage: s2 - s1    # syntactic sugar
            (-oo, -10] + [2, 3)
            sage: s1.difference(1,11)
            (0, 1] + [11, +oo)
        """
        other = RealSet(*other)
        return self.intersection(other.complement())

    __sub__ = difference

    def contains(self, x):
        """
        Return whether `x` is contained in the set

        INPUT:

        - ``x`` -- a real number.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: s = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s
            (0, 2) + [10, +oo)
            sage: s.contains(1)
            True
            sage: s.contains(0)
            False
            sage: 10 in s    # syntactic sugar
            True
        """
        x = RLF(x)
        for interval in self._intervals:
            if interval.contains(x):
                return True
        return False
    
    __contains__ = contains
    
    def is_included_in(self, *other):
        r"""
        Tests interval inclusion
            
        INPUT:

        - ``*args`` -- a :class:`RealSet` or something that defines
          one.

        OUTPUT:
        
        Boolean.

        EXAMPLES::

            sage: I = RealSet((1,2))
            sage: J = RealSet((1,3))
            sage: K = RealSet((2,3))
            sage: I.is_included_in(J)
            True
            sage: J.is_included_in(K)
            False
        """
        return RealSet(*other).intersection(self) == self

    def an_element(self):
        """
        Return a point of the set

        OUTPUT:

        A real number. ``ValueError`` if the set is empty.

        EXAMPLES::

            sage: RealSet.open_closed(0, 1).an_element()
            1
            sage: RealSet(0, 1).an_element()
            1/2
            sage: RealSet(-oo,+oo).an_element()
            0
            sage: RealSet(-oo,7).an_element()
            6
            sage: RealSet(7,+oo).an_element()
            8
        """
        from sage.rings.infinity import AnInfinity
        if len(self._intervals) == 0:
            raise ValueError('set is empty')
        i = self._intervals[0]
        if isinstance(i.lower(), AnInfinity):
            if isinstance(i.upper(), AnInfinity):
                return ZZ.zero()
            else:
                return i.upper() - 1
        if isinstance(i.upper(), AnInfinity):
            return i.lower() + 1
        if i.lower_closed():
            return i.lower()
        if i.upper_closed():
            return i.upper()
        return (i.lower() + i.upper())/ZZ(2)

    def is_disjoint_from(self, *other):
        """
        Test whether the two sets are disjoint

        INPUT:

        - ``other`` -- a :class:`RealSet` or data defining one.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: s1 = RealSet((0, 1), (2, 3));  s1
            (0, 1) + (2, 3)
            sage: s2 = RealSet([1, 2]);  s2
            [1, 2]
            sage: s1.is_disjoint_from(s2)
            True
            sage: s1.is_disjoint_from([1, 2])
            True
        """
        other = RealSet(*other)
        return self.intersection(other).is_empty()

    @staticmethod
    def are_pairwise_disjoint(*real_set_collection):
        """
        Test whether sets are pairwise disjoint
       
        INPUT:

        - ``*real_set_collection`` -- a list/tuple/iterable of
          :class:`RealSet`.
 
        OUTPUT:
        
        Boolean.
        
        EXAMPLES::

            sage: s1 = RealSet((0, 1), (2, 3))
            sage: s2 = RealSet((1, 2))
            sage: s3 = RealSet.point(3)
            sage: RealSet.are_pairwise_disjoint(s1, s2, s3)
            True
            sage: RealSet.are_pairwise_disjoint(s1, s2, s3, [10,10])
            True
            sage: RealSet.are_pairwise_disjoint(s1, s2, s3, [-1, 1/2])
            False
        """
        sets = [RealSet(_) for _ in real_set_collection]
        for i in range(len(sets)):
            for j in range(i):
                si = sets[i]
                sj = sets[j]
                if not si.is_disjoint_from(sj):
                    return False
        return True

    def _sage_input_(self, sib, coerced):
        """
        Produce an expression which will reproduce this value when evaluated.

        TESTS::

            sage: sage_input(RealSet())
            RealSet()
            sage: sage_input(RealSet.open(-oo, +oo))
            RealSet(-oo, oo)
            sage: sage_input(RealSet.point(77))
            RealSet.point(77)
            sage: sage_input(RealSet.closed_open(0, +oo))
            RealSet.closed_open(0, oo)
            sage: sage_input(RealSet.open_closed(-oo, 0))
            RealSet.open_closed(-oo, 0)
            sage: sage_input(RealSet.open_closed(-1, 0))
            RealSet.open_closed(-1, 0)
            sage: sage_input(RealSet.closed_open(-1, 0))
            RealSet.closed_open(-1, 0)
            sage: sage_input(RealSet.closed(0, 1))
            RealSet.closed(0, 1)
            sage: sage_input(RealSet.open(0, 1))
            RealSet.open(0, 1)
            sage: sage_input(RealSet.open(0, 1) + RealSet.open(1, 2))
            RealSet.open(0, 1) + RealSet.open(1, 2)
        """

        def interval_input(i):
            lower, upper = i.lower(), i.upper()
            if i.is_point():
                return sib.name('RealSet.point')(lower)
            elif lower == minus_infinity and upper == infinity:
                return sib.name('RealSet')(sib(minus_infinity), sib(infinity))
            else:
                if i.lower_closed():
                    if i.upper_closed():
                        t = 'RealSet.closed'
                    else:
                        t = 'RealSet.closed_open'
                else:
                    if i.upper_closed():
                        t = 'RealSet.open_closed'
                    else:
                        t = 'RealSet.open'
                return sib.name(t)(sib(lower), sib(upper))

        if self.is_empty():
            return sib.name('RealSet')()
        else:
            return sib.sum(interval_input(i) for i in self)

    def __mul__(self, right):
        r"""
        Scale a real set by a scalar on the left or right.

        EXAMPLES::

            sage: A = RealSet([0, 1/2], (2, infinity)); A
            [0, 1/2] + (2, +oo)
            sage: 2 * A
            [0, 1] + (4, +oo)
            sage: A * 100
            [0, 50] + (200, +oo)
            sage: 1.5 * A
            [0.000000000000000, 0.750000000000000] + (3.00000000000000, +oo)
            sage: (-2) * A
            (-oo, -4) + [-1, 0]
        """
        if not isinstance(right, RealSet):
            return RealSet(*[e * right for e in self])
        elif not isinstance(self, RealSet):
            return RealSet(*[self * e for e in right])
        else:
            return NotImplemented

    def __rmul__(self, other):
        r"""
        Scale a real set by a scalar on the left.

        TESTS::

            sage: A = RealSet([0, 1/2], RealSet.unbounded_above_closed(2)); A
            [0, 1/2] + [2, +oo)
            sage: pi * A
            [0, 1/2*pi] + [2*pi, +oo)
        """
        return self * other
