r"""
Tools for generating lists of integers in lexicographic order

IMPORTANT NOTE (2009/02):
The internal functions in this file will be deprecated soon.
Please only use them through :class:`IntegerListsLex`.

AUTHORS:

- Mike Hansen

- Travis Scrimshaw (2012-05-12): Fixed errors when returning ``None`` from
  first. Added checks to make sure ``max_slope`` is satisfied.

- Travis Scrimshaw (2012-10-29): Made ``IntegerListsLex`` into a parent with
  the element class ``IntegerListsLexElement``.
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from sage.arith.all import binomial
from sage.rings.integer_ring import ZZ
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.list_clone import ClonableArray
from sage.misc.stopgap import stopgap


def first(n, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Return the lexicographically smallest valid composition of ``n``
    satisfying the conditions.

    .. warning::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    .. TODO::

       Move this into Cython.

    Preconditions:

    - ``floor`` and ``ceiling`` need to satisfy the slope constraints,
      e.g. be obtained ``fromcomp2floor`` or ``comp2ceil``

    - ``floor`` must be below ``ceiling`` to ensure
      the existence a valid composition

    TESTS::

        sage: import sage.combinat.integer_list_old as integer_list
        sage: f = lambda l: lambda i: l[i-1]
        sage: f([0,1,2,3,4,5])(1)
        0
        sage: integer_list.first(12, 4, 4, f([0,0,0,0]), f([4,4,4,4]), -1, 1)
        [4, 3, 3, 2]
        sage: integer_list.first(36, 9, 9, f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -1, 1)
        [7, 6, 5, 5, 4, 3, 3, 2, 1]
        sage: integer_list.first(25, 9, 9, f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        [7, 6, 5, 4, 2, 1, 0, 0, 0]
        sage: integer_list.first(36, 9, 9, f([3,3,3,2,1,4,2,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        [7, 6, 5, 5, 5, 4, 3, 1, 0]

    ::

        sage: I = integer_list.IntegerListsLex(6, max_slope=2, min_slope=2)
        sage: list(I)
        [[6], [2, 4], [0, 2, 4]]
    """
    stopgap("First uses the old implementation of IntegerListsLex, which does not allow for arbitrary input;"
            " non-allowed input can return wrong results,"
            " please see the documentation for IntegerListsLex for details.",
            17548)
    # Check trivial cases, and standardize min_length to be at least 1
    if n < 0:
        return None
    if max_length <= 0:
        if n == 0:
            return []
        return None
    if min_length <= 0:
        if n == 0:
            return []
        min_length = 1

    # Increase min_length until n <= sum([ceiling(i) for i in range(min_length)])
    # This may run forever!
    # Find the actual length the list needs to be
    N = 0
    for i in range(1, min_length + 1):
        ceil = ceiling(i)
        if ceil < floor(i):
            return None
        N += ceil
    while N < n:
        min_length += 1
        if min_length > max_length:
            return None

        ceil = ceiling(min_length)
        if ceil == 0 and max_slope <= 0 or ceil < floor(min_length):
            return None

        N += ceil

    # Trivial case
    if min_length == 1:
        if n < floor(1):
            return None
        return [n]

    if max_slope < min_slope:
        return None

    # Compute the minimum values
    # We are constrained below by the max slope
    result = [floor(min_length)]
    n -= floor(min_length)
    for i in reversed(range(1, min_length)):
        result.insert(0, max(floor(i), result[0] - max_slope))
        n -= result[0]
        if n < 0:
            return None

    if n == 0:  # There is nothing more to do
        return result

    if min_slope == float('-inf'):
        for i in range(1, min_length+1):
            if n <= ceiling(i) - result[i-1]:  # -1 for indexing
                result[i-1] += n
                break
            else:
                n -= ceiling(i) - result[i-1]
                result[i-1] = ceiling(i)
    else:
        low_x = 1
        low_y = result[0]
        high_x = 1
        high_y = result[0]

        while n > 0:
            # invariant after each iteration of the loop:
            # [low_x, low_y] is the coordinate of the rightmost point of the
            # current diagonal s.t. result[low_x] < low_y
            low_y += 1
            while low_x < min_length and low_y + min_slope > result[low_x]:
                low_x += 1
                low_y += min_slope

            high_y += 1
            while high_y > ceiling(high_x):
                high_x += 1
                high_y += min_slope

            n -= low_x - high_x + 1

        for j in range(1, high_x):
            result[j-1] = ceiling(j)
        for i in range(0, -n):
            result[high_x+i-1] = high_y + min_slope * i - 1
        for i in range(-n, low_x-high_x+1):
            result[high_x+i-1] = high_y + min_slope * i

    # Special check for equal slopes
    if min_slope == max_slope and any(val + min_slope != result[i + 1]
                                      for i, val in enumerate(result[:-1])):
            return None

    return result


def lower_regular(comp, min_slope, max_slope):
    """
    Return the uppest regular composition below ``comp``

    TESTS::

        sage: import sage.combinat.integer_list_old as integer_list
        sage: integer_list.lower_regular([4,2,6], -1, 1)
        [3, 2, 3]
        sage: integer_list.lower_regular([4,2,6], -1, infinity)
        [3, 2, 6]
        sage: integer_list.lower_regular([1,4,2], -1, 1)
        [1, 2, 2]
        sage: integer_list.lower_regular([4,2,6,3,7], -2, 1)
        [4, 2, 3, 3, 4]
        sage: integer_list.lower_regular([4,2,infinity,3,7], -2, 1)
        [4, 2, 3, 3, 4]
        sage: integer_list.lower_regular([1, infinity, 2], -1, 1)
        [1, 2, 2]
        sage: integer_list.lower_regular([infinity, 4, 2], -1, 1)
        [4, 3, 2]
    """
    new_comp = comp[:]
    for i in range(1, len(new_comp)):
        new_comp[i] = min(new_comp[i], new_comp[i-1] + max_slope)

    for i in reversed(range(len(new_comp)-1)):
        new_comp[i] = min(new_comp[i], new_comp[i+1] - min_slope)

    return new_comp


def rightmost_pivot(comp, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    TESTS::

        sage: import sage.combinat.integer_list_old as integer_list
        sage: f = lambda l: lambda i: l[i-1]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9, f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -1, 0)
        [7, 2]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 0)
        [7, 1]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 4)
        [8, 1]
        sage: integer_list.rightmost_pivot([7,6,5,5,4,3,3,2,1], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        [8, 1]
        sage: integer_list.rightmost_pivot([7,6,5,5,5,5,5,4,4], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        sage: integer_list.rightmost_pivot([3,3,3,2,1,1,0,0,0], 9, 9,f([3,3,3,2,1,1,0,0,0]), f([7,6,5,5,5,5,5,4,4]), -2, 1)
        sage: g = lambda x: lambda i: x
        sage: integer_list.rightmost_pivot([1],1,1,g(0),g(2),-10, 10)
        sage: integer_list.rightmost_pivot([1,2],2,2,g(0),g(2),-10, 10)
        sage: integer_list.rightmost_pivot([1,2],2,2,g(1),g(2), -10, 10)
        sage: integer_list.rightmost_pivot([1,2],2,3,g(1),g(2), -10, 10)
        [2, 1]
        sage: integer_list.rightmost_pivot([2,2],2,3,g(2),g(2),-10, 10)
        sage: integer_list.rightmost_pivot([2,3],2,3,g(2),g(2),-10,+10)
        sage: integer_list.rightmost_pivot([3,2],2,3,g(2),g(2),-10,+10)
        sage: integer_list.rightmost_pivot([3,3],2,3,g(2),g(2),-10,+10)
        [1, 2]
        sage: integer_list.rightmost_pivot([6],1,3,g(0),g(6),-1,0)
        [1, 0]
        sage: integer_list.rightmost_pivot([6],1,3,g(0),g(6),-2,0)
        [1, 0]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(0),g(10),-1,10)
        [2, 6]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(5),g(10),-10,10)
        [3, 5]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(5),g(10),-1,10)
        [2, 6]
        sage: integer_list.rightmost_pivot([7,9,8,7],1,5,g(4),g(10),-2,10)
        [3, 7]
        sage: integer_list.rightmost_pivot([9,8,7],1,4,g(4),g(10),-2,0)
        [1, 4]
        sage: integer_list.rightmost_pivot([1,3],1,5,lambda i: i,g(10),-10,10)
        sage: integer_list.rightmost_pivot([1,4],1,5,lambda i: i,g(10),-10,10)
        sage: integer_list.rightmost_pivot([2,4],1,5,lambda i: i,g(10),-10,10)
        [1, 1]
    """
    if max_slope < min_slope:
        return None

    x = len(comp)
    if x == 0:
        return None

    y = len(comp) + 1
    while y <= max_length:
        if ceiling(y) > 0:
            break
        if max_slope <= 0:
            y = max_length + 1
            break
        y += 1

    ceilingx_x = comp[x - 1]-1
    floorx_x = floor(x)
    if x > 1:
        floorx_x = max(floorx_x, comp[x - 2] + min_slope)

    F = comp[x - 1] - floorx_x
    G = ceilingx_x - comp[x - 1]  # this is -1

    highX = x
    lowX = x

    while not (ceilingx_x >= floorx_x and
               (G >= 0 or
               (y < max_length + 1 and
                F - max(floor(y), floorx_x + (y - x) * min_slope) >= 0 and
                G + min(ceiling(y), ceilingx_x + (y - x) * max_slope) >= 0))):

        if x == 1:
            return None

        x -= 1

        oldfloorx_x = floorx_x
        ceilingx_x = comp[x-1] - 1
        floorx_x = floor(x)
        if x > 1:
            floorx_x = max(floorx_x, comp[x-2] + min_slope)

        min_slope_lowX = min_slope * (lowX - x)
        max_slope_highX = max_slope * (highX - x)

        # Update G
        if max_slope == float('+inf'):
            # In this case, we have
            #  -- ceiling_x(i) = ceiling(i) for i > x
            #  --G >= 0 or G = -1
            G += ceiling(x + 1) - comp[x]
        else:
            G += (highX - x) * ((comp[x-1] + max_slope) - comp[x]) - 1
            temp = (ceilingx_x + max_slope_highX) - ceiling(highX)
            while highX > x and temp >= 0:
                G -= temp
                highX -= 1
                max_slope_highX = max_slope * (highX-x)
                temp = (ceilingx_x + max_slope_highX) - ceiling(highX)

        if G >= 0 and comp[x-1] > floorx_x:
            # By case 1, x is at the rightmost pivot position
            break

        # Update F
        if y < max_length+1:
            F += comp[x-1] - floorx_x
            if min_slope != float('-inf'):
                F += (lowX - x) * (oldfloorx_x - (floorx_x + min_slope))
                temp = floor(lowX) - (floorx_x + min_slope_lowX)
                while lowX > x and temp >= 0:
                    F -= temp
                    lowX -= 1
                    min_slope_lowX = min_slope * (lowX-x)
                    temp = floor(lowX) - (floorx_x + min_slope_lowX)

    return [x, floorx_x]


def next(comp, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Return the next integer list after ``comp`` that satisfies the
    constraints.

    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    EXAMPLES::

        sage: from sage.combinat.integer_list_old import next
        sage: IV = sage.combinat.integer_list_old.IntegerListsLex(n=2,length=3,min_slope=0)
        sage: next([0,1,1], 3, 3, lambda i: 0, lambda i: 5, 0, 10)
        [0, 0, 2]
    """
    stopgap("Next uses the old implementation of IntegerListsLex, which does not allow for arbitrary input;"
            " non-allowed input can return wrong results,"
            " please see the documentation for IntegerListsLex for details.",
            17548)
    x = rightmost_pivot(comp, min_length, max_length, floor, ceiling, min_slope, max_slope)
    if x is None:
        return None
    [x, low] = x
    high = comp[x-1] - 1

#     // Build wrappers around floor and ceiling to take into
#     // account the new constraints on the value of compo[x].
#     //
#     // Efficiency note: they are not wrapped more than once, since
#     // the method Next calls first, but not the converse.

    if min_slope == float('-inf'):
        new_floor = lambda i: floor(x+(i-1))
    else:
        new_floor = lambda i: max(floor(x+(i-1)), low+(i-1)*min_slope)

    if max_slope == float('+inf'):
        new_ceiling = lambda i: comp[x-1] - 1 if i == 1 else ceiling(x+(i-1))
    else:
        new_ceiling = lambda i: min(ceiling(x+(i-1)), high+(i-1)*max_slope)

    res = []
    res += comp[:x-1]
    f = first(sum(comp[x-1:]), max(min_length-x+1, 0), max_length-x+1,
                 new_floor, new_ceiling, min_slope, max_slope)
    if f is None:  # Check to make sure it is valid
        return None
    res += f
    return res


def iterator(n, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    EXAMPLES::

        sage: from sage.combinat.integer_list_old import iterator
        sage: IV = sage.combinat.integer_list_old.IntegerListsLex(n=2,length=3,min_slope=0)
        sage: list(iterator(2, 3, 3, lambda i: 0, lambda i: 5, 0, 10))
        [[0, 1, 1], [0, 0, 2]]
    """
    stopgap("Iterator uses the old implementation of IntegerListsLex, which does not allow for arbitrary input;"
            " non-allowed input can return wrong results,"
            " please see the documentation for IntegerListsLex for details.",
            17548)
    succ = lambda x: next(x, min_length, max_length, floor, ceiling, min_slope, max_slope)

    # Handle the case where n is a list of integers
    if isinstance(n, list):
        for i in range(n[0], min(n[1]+1, upper_bound(min_length, max_length, floor, ceiling, min_slope, max_slope))):
            for el in iterator(i, min_length, max_length, floor, ceiling, min_slope, max_slope):
                yield el
    else:
        f = first(n, min_length, max_length, floor, ceiling, min_slope, max_slope)
        while f is not None:
            yield f
            f = succ(f)


def upper_regular(comp, min_slope, max_slope):
    """
    Return the uppest regular composition above ``comp``.

    TESTS::

        sage: import sage.combinat.integer_list_old as integer_list
        sage: integer_list.upper_regular([4,2,6],-1,1)
        [4, 5, 6]
        sage: integer_list.upper_regular([4,2,6],-2, 1)
        [4, 5, 6]
        sage: integer_list.upper_regular([4,2,6,3,7],-2, 1)
        [4, 5, 6, 6, 7]
        sage: integer_list.upper_regular([4,2,6,1], -2, 1)
        [4, 5, 6, 4]
    """

    new_comp = comp[:]
    for i in range(1, len(new_comp)):
        new_comp[i] = max(new_comp[i], new_comp[i-1] + min_slope)

    for i in reversed(range(len(new_comp)-1)):
        new_comp[i] = max(new_comp[i], new_comp[i+1] - max_slope)

    return new_comp


def comp2floor(f, min_slope, max_slope):
    """
    Given a composition, returns the lowest regular function N->N above
    it.

    EXAMPLES::

        sage: from sage.combinat.integer_list_old import comp2floor
        sage: f = comp2floor([2, 1, 1],-1,0)
        sage: [f(i) for i in range(10)]
        [2, 1, 1, 1, 2, 3, 4, 5, 6, 7]
    """
    if not f:
        return lambda i: 0
    floor = upper_regular(f, min_slope, max_slope)
    return lambda i: floor[i] if i < len(floor) else max(0, floor[-1]-(i-len(floor))*min_slope)


def comp2ceil(c, min_slope, max_slope):
    """
    Given a composition, return the lowest regular function N->N below it.

    EXAMPLES::

        sage: from sage.combinat.integer_list_old import comp2ceil
        sage: f = comp2ceil([2, 1, 1],-1,0)
        sage: [f(i) for i in range(10)]
        [2, 1, 1, 1, 2, 3, 4, 5, 6, 7]
    """
    if not c:
        return lambda i: 0
    ceil = lower_regular(c, min_slope, max_slope)
    return lambda i: ceil[i] if i < len(ceil) else max(0, ceil[-1]-(i-len(ceil))*min_slope)


def upper_bound(min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Compute a coarse upper bound on the size of a vector satisfying the
    constraints.

    TESTS::

        sage: import sage.combinat.integer_list_old as integer_list
        sage: f = lambda x: lambda i: x
        sage: integer_list.upper_bound(0,4,f(0), f(1),-infinity,infinity)
        4
        sage: integer_list.upper_bound(0, infinity, f(0), f(1), -infinity, infinity)
        inf
        sage: integer_list.upper_bound(0, infinity, f(0), f(1), -infinity, -1)
        1
        sage: integer_list.upper_bound(0, infinity, f(0), f(5), -infinity, -1)
        15
        sage: integer_list.upper_bound(0, infinity, f(0), f(5), -infinity, -2)
        9
    """
    from sage.functions.all import floor as flr
    if max_length < float('inf'):
        return sum([ceiling(j) for j in range(max_length)])
    elif max_slope < 0 and ceiling(1) < float('inf'):
        maxl = flr(-ceiling(1) / max_slope)
        return ceiling(1) * (maxl + 1) + binomial(maxl + 1, 2) * max_slope
    # FIXME: only checking the first 10000 values, but that should generally
    # be enough
    elif [ceiling(j) for j in range(10000)] == [0] * 10000:
        return 0
    else:
        return float('inf')


def is_a(comp, min_length, max_length, floor, ceiling, min_slope, max_slope):
    """
    Return ``True`` if ``comp`` meets the constraints imposed by the
    arguments.

    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    EXAMPLES::

        sage: from sage.combinat.integer_list_old import is_a
        sage: IV = sage.combinat.integer_list_old.IntegerListsLex(n=2,length=3,min_slope=0)
        sage: all(is_a(iv, 3, 3, lambda i: 0, lambda i: 5, 0, 10) for iv in IV)
        True
    """
    if not(min_length <= len(comp) <= max_length):
        return False
    for i in range(len(comp)):
        if not(floor(i + 1) <= comp[i] <= ceiling(i + 1)):
            return False
    for i in range(len(comp) - 1):
        slope = comp[i + 1] - comp[i]
        if not(min_slope <= slope <= max_slope):
            return False
    return True


class IntegerListsLexElement(ClonableArray):
    """
    Element class for :class:`IntegerListsLex`.
    """
    def check(self):
        """
        Check to make sure this is a valid element in its
        :class:`IntegerListsLex` parent.

        .. TODO:: Placeholder. Implement a proper check.

        EXAMPLES::

            sage: C = IntegerListsLex(4)
            sage: C([4]).check()
            True
        """
        return True


class IntegerListsLex(Parent):
    r"""
    A combinatorial class `C` for integer lists satisfying certain
    sum, length, upper/lower bound and regularity constraints. The
    purpose of this tool is mostly to provide a Constant Amortized
    Time iterator through those lists, in lexicographic order.

    INPUT:

    - ``n`` -- a non negative integer
    - ``min_length`` -- a non negative integer
    - ``max_length`` -- an integer or `\infty`
    - ``length`` -- an integer; overrides min_length and max_length if
      specified
    - ``min_part`` -- the minimum value of each part; defaults to ``0``
    - ``max_part`` -- the maximum value of each part; defaults to `+\infty`
    - ``floor`` -- a function `f` (or list);    defaults to
      ``lambda i: min_part``
    - ``ceiling`` -- a function `f` (or list);  defaults to
      ``lambda i: max_part``
    - ``min_slope`` -- an integer or `-\infty`; defaults to `-\infty`
    - ``max_slope`` -- an integer or `+\infty`; defaults to `+\infty`

    An *integer list* is a list `l` of nonnegative integers, its
    *parts*. The *length* of `l` is the number of its parts;
    the *sum* of `l` is the sum of its parts.

    .. NOTE::

       Two valid integer lists are considered equivalent if they only
       differ by trailing zeroes. In this case, only the list with the
       least number of trailing zeroes will be produced.

    The constraints on the lists are as follow:

    - Sum: `sum(l) == n`

    - Length: ``min_length <= len(l) <= max_length``

    - Lower and upper bounds: ``floor(i) <= l[i] <= ceiling(i)``, for
      ``i`` from 0 to ``len(l)``

    - Regularity condition: ``minSlope <= l[i+1]-l[i] <= maxSlope``,
      for ``i`` from 0 to ``len(l)-1``

    This is a generic low level tool. The interface has been designed
    with efficiency in mind. It is subject to incompatible changes in
    the future. More user friendly interfaces are provided by high
    level tools like :class:`Partitions` or :class:`Compositions`.

    EXAMPLES:

    We create the combinatorial class of lists of length 3 and sum 2::

        sage: import sage.combinat.integer_list_old as integer_list
        sage: C = integer_list.IntegerListsLex(2, length=3)
        sage: C
        Integer lists of sum 2 satisfying certain constraints
        sage: C.cardinality()
        6
        sage: [p for p in C]
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

        sage: [2, 0, 0] in C
        True
        sage: [2, 0, 1] in C
        False
        sage: "a" in C
        False
        sage: ["a"] in C
        False

        sage: C.first()
        [2, 0, 0]

    One can specify lower and upper bound on each part::

        sage: list(integer_list.IntegerListsLex(5, length = 3, floor = [1,2,0], ceiling = [3,2,3]))
        [[3, 2, 0], [2, 2, 1], [1, 2, 2]]

    Using the slope condition, one can generate integer partitions
    (but see :mod:`sage.combinat.partition.Partitions`)::

        sage: list(integer_list.IntegerListsLex(4, max_slope=0))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

    This is the list of all partitions of `7` with parts at least `2`::

        sage: list(integer_list.IntegerListsLex(7, max_slope = 0, min_part = 2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]

    This is the list of all partitions of `5` and length at most 3
    which are bounded below by [2,1,1]::

        sage: list(integer_list.IntegerListsLex(5, max_slope = 0, max_length = 3, floor = [2,1,1]))
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]

    Note that ``[5]`` is considered valid, because the lower bound
    constraint only apply to existing positions in the list. To
    obtain instead the partitions containing ``[2,1,1]``, one need to
    use ``min_length``::

        sage: list(integer_list.IntegerListsLex(5, max_slope = 0, min_length = 3, max_length = 3, floor = [2,1,1]))
        [[3, 1, 1], [2, 2, 1]]

    This is the list of all partitions of `5` which are contained in
    ``[3,2,2]``::

        sage: list(integer_list.IntegerListsLex(5, max_slope = 0, max_length = 3, ceiling = [3,2,2]))
        [[3, 2], [3, 1, 1], [2, 2, 1]]

    This is the list of all compositions of `4` (but see Compositions)::

        sage: list(integer_list.IntegerListsLex(4, min_part = 1))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    This is the list of all integer vectors of sum `4` and length `3`::

        sage: list(integer_list.IntegerListsLex(4, length = 3))
        [[4, 0, 0], [3, 1, 0], [3, 0, 1], [2, 2, 0], [2, 1, 1],
         [2, 0, 2], [1, 3, 0], [1, 2, 1], [1, 1, 2], [1, 0, 3],
         [0, 4, 0], [0, 3, 1], [0, 2, 2], [0, 1, 3], [0, 0, 4]]

    There are all the lists of sum 4 and length 4 such that l[i] <= i::

        sage: list(integer_list.IntegerListsLex(4, length=4, ceiling=lambda i: i))
        [[0, 1, 2, 1], [0, 1, 1, 2], [0, 1, 0, 3], [0, 0, 2, 2], [0, 0, 1, 3]]

    This is the list of all monomials of degree `4` which divide the
    monomial `x^3y^1z^2` (a monomial being identified with its
    exponent vector)::

        sage: R.<x,y,z> = QQ[]
        sage: m = [3,1,2]
        sage: def term(exponents):
        ....:     return x^exponents[0] * y^exponents[1] * z^exponents[2]
        sage: list( integer_list.IntegerListsLex(4, length = len(m), ceiling = m, element_constructor = term) )
        [x^3*y, x^3*z, x^2*y*z, x^2*z^2, x*y*z^2]

    Note the use of the element_constructor feature.

    In general, the complexity of the iteration algorithm is constant
    time amortized for each integer list produced.  There is one
    degenerate case though where the algorithm may run forever without
    producing anything. If max_length is `+\infty` and max_slope `>0`,
    testing whether there exists a valid integer list of sum `n` may
    be only semi-decidable. In the following example, the algorithm
    will enter an infinite loop, because it needs to decide whether
    `ceiling(i)` is nonzero for some `i`::

        sage: list( integer_list.IntegerListsLex(1, ceiling = lambda i: 0) ) # todo: not implemented

    .. NOTE::

       Caveat: counting is done by brute force generation. In some
       special cases, it would be possible to do better by counting
       techniques for integral point in polytopes.

    .. NOTE::

       Caveat: with the current implementation, the constraints should
       satisfy the following conditions:

       - The upper and lower bounds themselves should satisfy the
         slope constraints.

       - The maximal and minimal part values should not be equal.

    Those conditions are not always checked by the algorithm, and the
    result may be completely incorrect if they are not satisfied:

    In the following example, the floor conditions do not satisfy the
    slope conditions since the floor for the third part is also 3::

        sage: I = integer_list.IntegerListsLex(16, min_length=2, min_part=3, max_slope=-1, floor=[5,3])
        Traceback (most recent call last):
        ...
        ValueError: floor does not satisfy the max slope condition

    Compare this with the following input, which is equivalent
    but it bypasses the checks because the floor is a function::

        sage: f = lambda x: 5 if x == 0 else 3
        sage: I = integer_list.IntegerListsLex(16, min_length=2, max_slope=-1, floor=f)
        sage: list(I)
        [[13, 3], [12, 4], [11, 5], [10, 6]]

    With some work, this could be fixed without affecting the overall
    complexity and efficiency. Also, the generation algorithm could be
    extended to deal with non-constant slope constraints and with
    negative parts, as well as to accept a range parameter instead of
    a single integer for the sum `n` of the lists (the later was
    readily implemented in MuPAD-Combinat). Encouragements,
    suggestions, and help are welcome.

    .. TODO::

        Integrate all remaining tests from
        http://mupad-combinat.svn.sourceforge.net/viewvc/mupad-combinat/trunk/MuPAD-Combinat/lib/COMBINAT/TEST/MachineIntegerListsLex.tst

    TESTS::

        sage: g = lambda x: lambda i: x
        sage: list(integer_list.IntegerListsLex(0, floor = g(1), min_slope = 0))
        [[]]
        sage: list(integer_list.IntegerListsLex(0, floor = g(1), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(integer_list.IntegerListsLex(0, max_length=0, floor = g(1), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(integer_list.IntegerListsLex(0, max_length=0, floor = g(0), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(integer_list.IntegerListsLex(0, min_part = 1, min_slope = 0))
        [[]]
        sage: list(integer_list.IntegerListsLex(1, min_part = 1, min_slope = 0))
        [[1]]
        sage: list(integer_list.IntegerListsLex(0, min_length = 1, min_part = 1, min_slope = 0))
        []
        sage: list(integer_list.IntegerListsLex(0, min_length = 1, min_slope = 0))
        [[0]]
        sage: list(integer_list.IntegerListsLex(3, max_length=2, ))
        [[3], [2, 1], [1, 2], [0, 3]]
        sage: partitions = {"min_part": 1, "max_slope": 0}
        sage: partitions_min_2 = {"floor": g(2), "max_slope": 0}
        sage: compositions = {"min_part": 1}
        sage: integer_vectors = lambda l: {"length": l}
        sage: lower_monomials = lambda c: {"length": c, "floor": lambda i: c[i]}
        sage: upper_monomials = lambda c: {"length": c, "ceiling": lambda i: c[i]}
        sage: constraints = { "min_part":1, "min_slope": -1, "max_slope": 0}
        sage: list(integer_list.IntegerListsLex(6, **partitions))
        [[6],
         [5, 1],
         [4, 2],
         [4, 1, 1],
         [3, 3],
         [3, 2, 1],
         [3, 1, 1, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(integer_list.IntegerListsLex(6, **constraints))
        [[6],
         [3, 3],
         [3, 2, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(integer_list.IntegerListsLex(1, **partitions_min_2))
        []
        sage: list(integer_list.IntegerListsLex(2, **partitions_min_2))
        [[2]]
        sage: list(integer_list.IntegerListsLex(3, **partitions_min_2))
        [[3]]
        sage: list(integer_list.IntegerListsLex(4, **partitions_min_2))
        [[4], [2, 2]]
        sage: list(integer_list.IntegerListsLex(5, **partitions_min_2))
        [[5], [3, 2]]
        sage: list(integer_list.IntegerListsLex(6, **partitions_min_2))
        [[6], [4, 2], [3, 3], [2, 2, 2]]
        sage: list(integer_list.IntegerListsLex(7, **partitions_min_2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]
        sage: list(integer_list.IntegerListsLex(9, **partitions_min_2))
        [[9], [7, 2], [6, 3], [5, 4], [5, 2, 2], [4, 3, 2], [3, 3, 3], [3, 2, 2, 2]]
        sage: list(integer_list.IntegerListsLex(10, **partitions_min_2))
        [[10],
         [8, 2],
         [7, 3],
         [6, 4],
         [6, 2, 2],
         [5, 5],
         [5, 3, 2],
         [4, 4, 2],
         [4, 3, 3],
         [4, 2, 2, 2],
         [3, 3, 2, 2],
         [2, 2, 2, 2, 2]]
        sage: list(integer_list.IntegerListsLex(4, **compositions))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: list(integer_list.IntegerListsLex(6, min_length=1, floor=[7]))
        []

    Noted on :trac:`17898`::

        sage: list(integer_list.IntegerListsLex(4, min_part=1, length=3, min_slope=1))
        []
        sage: integer_list.IntegerListsLex(6, ceiling=[4,2], floor=[3,3]).list()
        []
        sage: integer_list.IntegerListsLex(6, min_part=1, max_part=3, max_slope=-4).list()
        []
    """
    def __init__(self,
                 n,
                 length = None, min_length=0, max_length=float('+inf'),
                 floor=None, ceiling = None,
                 min_part = 0, max_part = float('+inf'),
                 min_slope=float('-inf'), max_slope=float('+inf'),
                 name = None,
                 element_constructor = None,
                 element_class = None,
                 global_options = None):
        """
        Initialize ``self``.

        TESTS::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: C == loads(dumps(C)) # this did fail at some point, really!
            True
            sage: C is loads(dumps(C)) # todo: not implemented
            True
            sage: C.cardinality().parent() is ZZ
            True
            sage: TestSuite(C).run()
        """
        stopgap("The old implementation of IntegerListsLex does not allow for arbitrary input;"
                " non-allowed input can return wrong results,"
                " please see the documentation for IntegerListsLex for details.",
                17548)
        # Convert to float infinity
        from sage.rings.infinity import infinity
        if max_slope == infinity:
            max_slope = float('+inf')
        if min_slope == -infinity:
            min_slope = float('-inf')
        if max_length == infinity:
            max_length = float('inf')
        if max_part == infinity:
            max_part = float('+inf')

        if floor is None:
            self.floor_list = []
        else:
            try:
                # Is ``floor`` an iterable?
                # Not ``floor[:]`` because we want ``self.floor_list``
                #    mutable, and applying [:] to a tuple gives a tuple.
                self.floor_list = list(floor)
                # Make sure the floor list will make the list satisfy the constraints
                if min_slope != float('-inf'):
                    for i in range(1, len(self.floor_list)):
                        self.floor_list[i] = max(self.floor_list[i], self.floor_list[i-1] + min_slope)

                # Some input checking
                for i in range(1, len(self.floor_list)):
                    if self.floor_list[i] - self.floor_list[i-1] > max_slope:
                        raise ValueError("floor does not satisfy the max slope condition")
                if self.floor_list and min_part - self.floor_list[-1] > max_slope:
                    raise ValueError("floor does not satisfy the max slope condition")
            except TypeError:
                self.floor = floor
        if ceiling is None:
            self.ceiling_list = []
        else:
            try:
                # Is ``ceiling`` an iterable?
                self.ceiling_list = list(ceiling)
                # Make sure the ceiling list will make the list satisfy the constraints
                if max_slope != float('+inf'):
                    for i in range(1, len(self.ceiling_list)):
                        self.ceiling_list[i] = min(self.ceiling_list[i], self.ceiling_list[i-1] + max_slope)

                # Some input checking
                for i in range(1, len(self.ceiling_list)):
                    if self.ceiling_list[i] - self.ceiling_list[i-1] < min_slope:
                        raise ValueError("ceiling does not satisfy the min slope condition")
                if self.ceiling_list and max_part - self.ceiling_list[-1] < min_slope:
                    raise ValueError("ceiling does not satisfy the min slope condition")
            except TypeError:
                # ``ceiling`` is not an iterable.
                self.ceiling = ceiling
        if name is not None:
            self.rename(name)
        if n in ZZ:
            self.n = n
            self.n_range = [n]
        else:
            self.n_range = n
        if length is not None:
            min_length = length
            max_length = length
        self.min_length = min_length
        self.max_length = max_length
        self.min_part = min_part
        self.max_part = max_part
        # FIXME: the internal functions currently assume that floor and ceiling start at 1
        # this is a workaround
        self.max_slope = max_slope
        self.min_slope = min_slope
        if element_constructor is not None:
            self._element_constructor_ = element_constructor
        if element_class is not None:
            self.Element = element_class
        if global_options is not None:
            self.global_options = global_options
        Parent.__init__(self, category=FiniteEnumeratedSets())

    Element = IntegerListsLexElement

    def _element_constructor_(self, lst):
        """
        Construct an element with ``self`` as parent.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(4)
            sage: C([4])
            [4]
        """
        return self.element_class(self, lst)

    def __eq__(self, x):
        """
        Compare two different :class:`IntegerListsLex`.

        For now, the comparison is done just on their repr's which is
        not robust!

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: D = integer_list.IntegerListsLex(4, length=3)
            sage: repr(C) == repr(D)
            False
            sage: C == D
            False
        """
        return repr(self) == repr(x)

    def __ne__(self, other):
        """
        Compare two different :class:`IntegerListsLex`.

        For now, the comparison is done just on their repr's which is
        not robust!

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: D = integer_list.IntegerListsLex(4, length=3)
            sage: C != D
            True
        """
        return not self.__eq__(other)

    def __hash__(self):
        """
        Compute a hash for ``self``.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: h = hash(C)
        """
        return hash(repr(self)) ^ 53397379531

    def _repr_(self):
        """
        Return the name of this combinatorial class.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: C # indirect doctest
            Integer lists of sum 2 satisfying certain constraints

            sage: C = integer_list.IntegerListsLex([1,2,4], length=3)
            sage: C # indirect doctest
            Integer lists of sum in [1, 2, 4] satisfying certain constraints

            sage: C = integer_list.IntegerListsLex([1,2,4], length=3, name="A given name")
            sage: C
            A given name
        """
        if hasattr(self, "n"):
            return "Integer lists of sum %s satisfying certain constraints" % self.n

        return "Integer lists of sum in %s satisfying certain constraints" % self.n_range

    def floor(self, i):
        """
        Return the minimum part that can appear at the `i^{th}` position of
        any list produced.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(4, length=2, min_part=1)
            sage: C.floor(0)
            1
            sage: C = integer_list.IntegerListsLex(4, length=2, floor=[1,2])
            sage: C.floor(0)
            1
            sage: C.floor(1)
            2
        """
        if i < len(self.floor_list):
            return max(self.min_part, self.floor_list[i])
        if self.min_slope != float('-inf') and self.min_slope > 0:
            return self.min_part + (i - len(self.floor_list)) * self.min_slope
        return self.min_part

    def ceiling(self, i):
        """
        Return the maximum part that can appear in the `i^{th}`
        position in any list produced.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(4, length=2, max_part=3)
            sage: C.ceiling(0)
            3
            sage: C = integer_list.IntegerListsLex(4, length=2, ceiling=[3,2])
            sage: C.ceiling(0)
            3
            sage: C.ceiling(1)
            2
        """
        if i < len(self.ceiling_list):
            return min(self.max_part, self.ceiling_list[i])
        if self.max_slope != float('inf') and self.max_slope < 0:
            return self.max_part + (i - len(self.ceiling_list)) * self.max_slope
        return self.max_part

    # Temporary adapter to use the preexisting list/iterator/is_a function above.
    # FIXME: fix their specs so that floor and ceiling start from 0 instead of 1...
    # FIXME: integrate them as methods of this class
    def build_args(self):
        """
        Return a list of arguments that can be passed into the pre-existing
        ``first``, ``next``, ``is_a``, ... functions in this module.

        ``n`` is currently not included in this list.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: C.build_args()
            [3,
             3,
             <function ...<lambda> at 0x...>,
             <function ...<lambda> at 0x...>,
             -inf,
             inf]
        """
        return [self.min_length, self.max_length,
                lambda i: self.floor(i - 1), lambda i: self.ceiling(i - 1),
                self.min_slope, self.max_slope]

    def first(self):
        """
        Return the lexicographically maximal element in ``self``.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: C.first()
            [2, 0, 0]
        """
        # Make sure we have a valid return
        f = first(self.n_range[0], *(self.build_args()))
        if f is None:
            return None
        return self._element_constructor_(f)

    def __iter__(self):
        """
        Return an iterator for the elements of ``self``.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: list(C) #indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
        """
        args = self.build_args()
        for n in self.n_range:
            l = first(n, *args)
            while l is not None:
                yield self._element_constructor_(l)
                l = next(l, *args)

    def count(self):
        """
        Default brute force implementation of count by iteration
        through all the objects.

        Note that this skips the call to ``_element_constructor_``,
        unlike the default implementation.

        .. TODO::

            Do the iteration in place to save on copying time

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: C.cardinality() == C.count()
            True
        """
        args = self.build_args()
        c = ZZ(0)
        for n in self.n_range:
            l = first(n, *args)
            while l is not None:
                c += 1
                l = next(l, *args)
        return c

    def __contains__(self, v):
        """
        Return ``True`` if and only if ``v`` is in ``self``.

        EXAMPLES::

            sage: import sage.combinat.integer_list_old as integer_list
            sage: C = integer_list.IntegerListsLex(2, length=3)
            sage: [2, 0, 0] in C
            True
            sage: [2, 0] in C
            False
            sage: [3, 0, 0] in C
            False
            sage: all(v in C for v in C)
            True
        """
        if isinstance(v, self.element_class) or isinstance(v, list):
            return is_a(v, *(self.build_args())) and sum(v) in self.n_range
        return False
