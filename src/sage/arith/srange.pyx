# sage_setup: distribution = sagemath-categories
"""
Ranges and the ``[1,2,..,n]`` notation

AUTHORS:

- Jeroen Demeyer (2016-02-22): moved here from ``misc.py`` and cleaned
  up.
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from libc.math cimport ceil
from cysignals.signals cimport sig_check

from sage.rings.integer cimport Integer
from sage.structure.element cimport parent
from sage.structure.sequence import Sequence


def xsrange(start, end=None, step=1, universe=None, *, coerce=True, bint include_endpoint=False, endpoint_tolerance=1e-5):
    """
    Return an iterator over numbers
    ``start, start+step, ..., start+k*step``,
    where ``start+k*step < end`` and ``start+(k+1)*step >= end``.

    This provides one way to iterate over Sage integers as opposed to
    Python int's.  It also allows you to specify step sizes for such
    an iteration.

    INPUT:

    - ``start`` -- number (default: 0)

    - ``end`` -- number

    - ``step`` -- number (default: 1)

    - ``universe`` -- parent or type where all the elements should live
      (default: deduce from inputs)

    - ``coerce`` -- convert ``start``, ``end`` and ``step`` to the same
      universe (either the universe given in ``universe`` or the
      automatically detected universe)

    - ``include_endpoint`` -- whether or not to include the endpoint
      (default: ``False``). This is only relevant if ``end`` is actually of
      the form ``start + k*step`` for some integer `k`.

    ` ``endpoint_tolerance`` -- used to determine whether or not the
      endpoint is hit for inexact rings (default: 1e-5)

    OUTPUT: iterator

    Unlike :func:`range`, ``start`` and ``end`` can be any type of
    numbers, and the resulting iterator involves numbers of that type.

    .. warning::

       You need to be careful when using this function over inexact
       rings: the elements are computed via repeated addition rather
       than multiplication, which may produce slightly different
       results. For example::

           sage: sum([1.1] * 10) == 1.1 * 10
           False

       Also, the question of whether the endpoint is hit exactly for
       a given ``start + k*step`` is fuzzy for an inexact ring. If
       ``start + k*step = end`` for some `k` within
       ``endpoint_tolerance`` of being integral, it is considered an
       exact hit, thus avoiding spurious values falling just below the
       endpoint.

    EXAMPLES::

        sage: xsrange(10)
        <...generator object at 0x...>
        sage: for i in xsrange(1,5):
        ....:     print(i)
        1
        2
        3
        4

    See :func:`srange` for more examples.

    TESTS:

    Ranges can be very large, see :issue:`20094`::

        sage: it = xsrange(10^30, 10^100)
        sage: for i in range(5):
        ....:     print(next(it))
        1000000000000000000000000000000
        1000000000000000000000000000001
        1000000000000000000000000000002
        1000000000000000000000000000003
        1000000000000000000000000000004
    """
    if end is None:
        end = start
        start = 0

    if coerce:
        if universe is None:
            universe = Sequence([start, end, step]).universe()
        start, end, step = universe(start), universe(end), universe(step)

    if not step:
        raise ValueError("step argument must not be zero")

    count = (end - start) / step

    # If count is exact, set endpoint_tolerance to zero
    cdef bint count_is_exact
    try:
        count_is_exact = parent(count).is_exact()
    except Exception:
        count_is_exact = False

    # Round count up (with tolerance if applicable) and check if
    # endpoint is of the form start + k*step
    if count_is_exact:
        try:
            icount = count.ceil()
        except AttributeError:
            try:
                # Double minus to divide with ceil
                icount = -((start - end) // step)
            except TypeError:
                icount = ceil(count)
        if icount != count:
            include_endpoint = False
    else:
        fcount = count - endpoint_tolerance
        try:
            icount = fcount.ceil()
        except AttributeError:
            icount = ceil(fcount)
        if abs(count - icount) > endpoint_tolerance:
            include_endpoint = False

    icount = Integer(icount)

    if icount < 0:
        return

    cur = start
    # yield in chunks of 1024
    cdef long k
    while icount > 1024:
        for k in range(1024):
            yield cur
            cur += step
        icount -= 1024
        sig_check()

    for k in range(icount):
        yield cur
        cur += step
    if include_endpoint:
        yield end


def srange(*args, **kwds):
    r"""
    Return a list of numbers
    ``start, start+step, ..., start+k*step``,
    where ``start+k*step < end`` and ``start+(k+1)*step >= end``.

    This provides one way to iterate over Sage integers as opposed to
    Python int's.  It also allows you to specify step sizes for such
    an iteration.

    INPUT:

    - ``start`` -- number (default: 0)

    - ``end`` -- number

    - ``step`` -- number (default: 1)

    - ``universe -- parent or type where all the elements should live
      (default: deduce from inputs). This is only used if ``coerce`` is
      true.

    - ``coerce`` -- convert ``start``, ``end`` and ``step`` to the same
      universe (either the universe given in ``universe`` or the
      automatically detected universe)

    - ``include_endpoint`` -- whether or not to include the endpoint
      (default: ``False``). This is only relevant if ``end`` is actually of
      the form ``start + k*step`` for some integer `k`.

    ` ``endpoint_tolerance`` -- used to determine whether or not the
      endpoint is hit for inexact rings (default: ``1e-5``)

    OUTPUT: list

    .. NOTE::

       This function is called ``srange`` to distinguish
       it from the built-in Python ``range`` command.  The s
       at the beginning of the name stands for "Sage".

    .. SEEALSO::

        :func:`xsrange` -- iterator which is used to implement :func:`srange`.

    EXAMPLES::

        sage: v = srange(5); v
        [0, 1, 2, 3, 4]
        sage: type(v[2])
        <class 'sage.rings.integer.Integer'>
        sage: srange(1, 10)
        [1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: srange(10, 1, -1)
        [10, 9, 8, 7, 6, 5, 4, 3, 2]
        sage: srange(10,1,-1, include_endpoint=True)
        [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        sage: srange(1, 10, universe=RDF)
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

        sage: srange(1, 10, 1/2)
        [1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 5, 11/2, 6, 13/2, 7, 15/2, 8, 17/2, 9, 19/2]

        sage: # needs sage.rings.real_mpfr
        sage: srange(1, 5, 0.5)
        [1.00000000000000, 1.50000000000000, 2.00000000000000, 2.50000000000000,
         3.00000000000000, 3.50000000000000, 4.00000000000000, 4.50000000000000]
        sage: srange(0, 1, 0.4)
        [0.000000000000000, 0.400000000000000, 0.800000000000000]
        sage: srange(1.0, 5.0, include_endpoint=True)
        [1.00000000000000, 2.00000000000000, 3.00000000000000, 4.00000000000000,
         5.00000000000000]
        sage: srange(1.0, 1.1)
        [1.00000000000000]
        sage: srange(1.0, 1.0)
        []

        sage: V = VectorSpace(QQ, 2)                                                    # needs sage.modules
        sage: srange(V([0,0]), V([5,5]), step=V([2,2]))                                 # needs sage.modules
        [(0, 0), (2, 2), (4, 4)]

    Including the endpoint::

        sage: srange(0, 10, step=2, include_endpoint=True)
        [0, 2, 4, 6, 8, 10]
        sage: srange(0, 10, step=3, include_endpoint=True)
        [0, 3, 6, 9]

    Try some inexact rings::

        sage: srange(0.5, 1.1, 0.1, universe=RDF, include_endpoint=False)
        [0.5, 0.6, 0.7, 0.7999999999999999, 0.8999999999999999, 0.9999999999999999]
        sage: srange(0.5, 1, 0.1, universe=RDF, include_endpoint=False)
        [0.5, 0.6, 0.7, 0.7999999999999999, 0.8999999999999999]
        sage: srange(0.5, 0.9, 0.1, universe=RDF, include_endpoint=False)
        [0.5, 0.6, 0.7, 0.7999999999999999]
        sage: srange(0, 1.1, 0.1, universe=RDF, include_endpoint=True)
        [0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6, 0.7,
         0.7999999999999999, 0.8999999999999999, 0.9999999999999999, 1.1]
        sage: srange(0, 0.2, 0.1, universe=RDF, include_endpoint=True)
        [0.0, 0.1, 0.2]
        sage: srange(0, 0.3, 0.1, universe=RDF, include_endpoint=True)
        [0.0, 0.1, 0.2, 0.3]

    More examples::

        sage: Q = RationalField()
        sage: srange(1, 10, Q('1/2'))
        [1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 5, 11/2, 6, 13/2, 7, 15/2, 8, 17/2, 9, 19/2]
        sage: srange(1, 5, 0.5)                                                         # needs sage.rings.real_mpfr
        [1.00000000000000, 1.50000000000000, 2.00000000000000, 2.50000000000000,
         3.00000000000000, 3.50000000000000, 4.00000000000000, 4.50000000000000]
        sage: srange(0, 1, 0.4)                                                         # needs sage.rings.real_mpfr
        [0.000000000000000, 0.400000000000000, 0.800000000000000]

    Negative steps are also allowed::

        sage: srange(4, 1, -1)
        [4, 3, 2]
        sage: srange(4, 1, -1/2)
        [4, 7/2, 3, 5/2, 2, 3/2]

    TESTS:

    These are doctests from :issue:`6409`::

        sage: srange(1,QQ(0),include_endpoint=True)
        []
        sage: srange(1,QQ(0),-1,include_endpoint=True)
        [1, 0]

    Test :issue:`11753`::

        sage: srange(1,1,0)
        Traceback (most recent call last):
        ...
        ValueError: step argument must not be zero

    No problems with large lists::

        sage: srange(10^5) == list(range(10^5))
        True
    """
    return [x for x in xsrange(*args, **kwds)]


def ellipsis_iter(*args, step=None):
    """
    Same as ellipsis_range, but as an iterator (and may end with an
    Ellipsis).

    See also ellipsis_range.

    Use (1,2,...) notation.

    EXAMPLES::

        sage: A = ellipsis_iter(1,2,Ellipsis)
        sage: [next(A) for _ in range(10)]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: next(A)
        11
        sage: A = ellipsis_iter(1,3,5,Ellipsis)
        sage: [next(A) for _ in range(10)]
        [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
        sage: A = ellipsis_iter(1,2,Ellipsis,5,10,Ellipsis)
        sage: [next(A) for _ in range(10)]
        [1, 2, 3, 4, 5, 10, 11, 12, 13, 14]

    TESTS:

    These were carefully chosen tests, only to be changed if the
    semantics of ellipsis ranges change. In other words, if they don't
    pass, it's probably a bug in the implementation, not in the
    doctest.

    ::

        sage: list(1,..,10)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: list(1,3,..,10)
        [1, 3, 5, 7, 9]
        sage: list(1,..,10,..,20)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        sage: list(1,3,..,10,..,20)
        [1, 3, 5, 7, 9, 10, 12, 14, 16, 18, 20]
        sage: list(1,3,..,10,10,..,20)
        [1, 3, 5, 7, 9, 10, 12, 14, 16, 18, 20]
        sage: list(0,2,..,10,10,..,20,20,..,25)
        [0, 2, 4, 6, 8, 10, 10, 12, 14, 16, 18, 20, 20, 22, 24]
        sage: list(10,..,1)
        []
        sage: list(10,11,..,1)
        []
        sage: list(10,9,..,1)
        [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        sage: list(100,..,10,..,20)
        [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        sage: list(0,..,10,..,-20)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: list(100,..,10,..,-20)
        []
        sage: list(100,102,..,10,..,20)
        [10, 12, 14, 16, 18, 20]
    """
    step_magic = 0
    if step is None:
        step = 1
        if Ellipsis in args:
            i = list(args).index(Ellipsis)
            if i > 1:
                step = args[i-1]-args[i-2]
                step_magic = i

    S = Sequence([a for a in args if a is not Ellipsis] + [step])
    universe = S.universe()
    args = [Ellipsis if a is Ellipsis else universe(a) for a in args]
    step = universe(step)

    # this is a bit more complicated because we can't pop what's already been yielded
    next_ = None
    skip = False
    last_end = None
    # first we handle step_magic (which may require two pops if the range is empty)
    if step_magic:
        for i in range(step_magic-2):
            yield args[i]
        if len(args) > step_magic+1:
            i = step_magic
            more = xsrange(args[i-2], args[i+1], step, coerce=False, include_endpoint=True)
            for a in more:
                last_end = a
                yield a
            skip = True
            next_ = None
            step_magic += 1
        else:
            yield args[step_magic-2]

    # now onto the rest
    for i in range(step_magic, len(args)):
        if skip:
            skip = False
        elif args[i] is Ellipsis:
            if i == len(args)-1:
                # continue forever
                cur = args[i-1]
                if last_end != cur:
                    yield cur
                while True:
                    cur += step
                    yield cur
            start, end = args[i-1], args[i+1]
            if i < 2 or args[i-2] is not Ellipsis:
                next_ = None
            more = xsrange(start, end, step, coerce=False, include_endpoint=True)
            try:
                first = next(more)
                if last_end != first:
                    yield first
                for a in more:
                    last_end = a
                    yield a
            except StopIteration:
                last_end = None
            skip = True
            next_ = None
        else:
            if next_ is not None:
                yield next_
            next_ = args[i]
            last_end = None


def ellipsis_range(*args, step=None):
    """
    Return arithmetic sequence determined by the numeric arguments and
    ellipsis. Best illustrated by examples.

    Use [1,2,..,n] notation.

    EXAMPLES::

        sage: ellipsis_range(1,Ellipsis,11,100)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 100]
        sage: ellipsis_range(0,2,Ellipsis,10,Ellipsis,20)
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        sage: ellipsis_range(0,2,Ellipsis,11,Ellipsis,20)
        [0, 2, 4, 6, 8, 10, 11, 13, 15, 17, 19]
        sage: ellipsis_range(0,2,Ellipsis,11,Ellipsis,20, step=3)
        [0, 2, 5, 8, 11, 14, 17, 20]
        sage: ellipsis_range(10,Ellipsis,0)
        []

    TESTS:

    These were carefully chosen tests, only to be changed if the
    semantics of ellipsis ranges change. In other words, if they don't
    pass it's probably a bug in the implementation, not in the
    doctest.

    Note 10 only appears once (though it is in both ranges).

    ::

        sage: ellipsis_range(0,Ellipsis,10,Ellipsis,20,step=2)
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

    Sometimes one or more ranges is empty.

    ::

        sage: ellipsis_range(100,Ellipsis,10,Ellipsis,20,step=2)
        [10, 12, 14, 16, 18, 20]
        sage: ellipsis_range(0,Ellipsis,10,Ellipsis,-20,step=2)
        [0, 2, 4, 6, 8, 10]
        sage: ellipsis_range(100,Ellipsis,10,Ellipsis,-20,step=2)
        []

    We always start on the leftmost point of the range.

    ::

        sage: ellipsis_range(0,Ellipsis,10,Ellipsis,20,step=3)
        [0, 3, 6, 9, 10, 13, 16, 19]
        sage: ellipsis_range(100,Ellipsis,10,Ellipsis,20,step=3)
        [10, 13, 16, 19]
        sage: ellipsis_range(0,Ellipsis,10,Ellipsis,-20,step=3)
        [0, 3, 6, 9]
        sage: ellipsis_range(100,Ellipsis,10,Ellipsis,-20,step=3)
        []
        sage: ellipsis_range(0,1,Ellipsis,-10)
        []
        sage: ellipsis_range(0,1,Ellipsis,-10,step=1)
        [0]
        sage: ellipsis_range(100,0,1,Ellipsis,-10)
        [100]

    Note the duplicate 5 in the output.

    ::

        sage: ellipsis_range(0,Ellipsis,5,5,Ellipsis,10)
        [0, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10]

    Examples in which the step determines the parent of the elements::

        sage: [1..3, step=0.5]                                                          # needs sage.rings.real_mpfr
        [1.00000000000000, 1.50000000000000, 2.00000000000000, 2.50000000000000, 3.00000000000000]
        sage: v = [1..5, step=1/1]; v
        [1, 2, 3, 4, 5]
        sage: parent(v[2])
        Rational Field
    """
    step_magic = 0
    if step is None:
        step = 1
        if Ellipsis in args:
            i = args.index(Ellipsis)
            if i > 1:
                step = args[i-1]-args[i-2]
                step_magic = i

    S = Sequence([a for a in args if a is not Ellipsis] + [step])
    universe = S.universe()
    args = [Ellipsis if a is Ellipsis else universe(a) for a in args]
    step = universe(step)

    skip = False
    last_end = None
    L = []
    for i in range(len(args)):
        if skip:
            skip = False
        elif args[i] is Ellipsis:
            if len(args) == i+1:
                raise IndexError("ellipsis range must have an endpoint, use (n..) for infinite sequence")
            start, end = args[i-1], args[i+1]
            if i < 2 or args[i-2] is not Ellipsis:
                L.pop()
                if i == step_magic:
                    L.pop()
                    start = args[i-2]
            more = srange(start, end, step, coerce=False, include_endpoint=True)
            if more:
                if last_end == more[0]:
                    L.pop()
                last_end = more[-1]
                L += more
            else:
                last_end = None
            skip = True
        else:
            L.append(args[i])
            last_end = None
    return L
