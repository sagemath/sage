r"""
Enumerated set of lists of integers with constraints: base classes

- :class:`IntegerListsBackend`: base class for the Cython back-end of
  an enumerated set of lists of integers with specified constraints.

- :class:`Envelope`: a utility class for upper (lower) envelope of a
  function under constraints.
"""

# ****************************************************************************
#       Copyright (C) 2015 Bryan Gillespie <Brg008@gmail.com>
#                          Nicolas M. Thiery <nthiery at users.sf.net>
#                          Anne Schilling <anne@math.ucdavis.edu>
#                          Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from cpython.object cimport Py_LE, Py_EQ, Py_NE, Py_GE
from sage.misc.constant_function import ConstantFunction
from sage.structure.element cimport RingElement
from sage.rings.integer cimport Integer

Infinity = float('+inf')
MInfinity = float('-inf')


cdef class IntegerListsBackend():
    """
    Base class for the Cython back-end of an enumerated set of lists of
    integers with specified constraints.

    This base implements the basic operations, including checking for
    containment using :meth:`_contains`, but not iteration. For
    iteration, subclass this class and implement an ``_iter()`` method.

    EXAMPLES::

        sage: from sage.combinat.integer_lists.base import IntegerListsBackend
        sage: L = IntegerListsBackend(6, max_slope=-1)
        sage: L._contains([3,2,1])
        True
    """
    def __init__(self,
                 n=None, length=None, *,
                 min_length=0, max_length=Infinity,
                 floor=None, ceiling=None,
                 min_part=0, max_part=Infinity,
                 min_slope=MInfinity, max_slope=Infinity,
                 min_sum=0, max_sum=Infinity):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.integer_lists.base import IntegerListsBackend
            sage: C = IntegerListsBackend(2, length=3)
            sage: C = IntegerListsBackend(min_sum=1.4)                                  # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
            sage: C = IntegerListsBackend(min_sum=Infinity)
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce <class 'sage.rings.infinity.PlusInfinity'> to an integer
        """
        if n is not None:
            min_sum = n
            max_sum = n
        self.min_sum = Integer(min_sum) if min_sum != -Infinity else -Infinity
        self.max_sum = Integer(max_sum) if max_sum != Infinity else Infinity

        if length is not None:
            min_length = length
            max_length = length
        self.min_length = Integer(max(min_length, 0))
        self.max_length = Integer(max_length) if max_length != Infinity else Infinity

        self.min_slope = Integer(min_slope) if min_slope != -Infinity else -Infinity
        self.max_slope = Integer(max_slope) if max_slope != Infinity else Infinity

        self.min_part = Integer(min_part) if min_part != -Infinity else -Infinity
        self.max_part = Integer(max_part) if max_part != Infinity else Infinity

        if isinstance(floor, Envelope):
            self.floor = floor
        else:
            if floor is None:
                floor = -Infinity
            elif isinstance(floor, (list, tuple)):
                floor = tuple(Integer(i) for i in floor)
            elif callable(floor):
                pass
            else:
                raise TypeError("floor should be a list, tuple, or function")
            self.floor = Envelope(floor, sign=-1,
                    min_part=self.min_part, max_part=self.max_part,
                    min_slope=self.min_slope, max_slope=self.max_slope,
                    min_length=self.min_length)

        if isinstance(ceiling, Envelope):
            self.ceiling = ceiling
        else:
            if ceiling is None:
                ceiling = Infinity
            elif isinstance(ceiling, (list, tuple)):
                ceiling = tuple(Integer(i) if i != Infinity else Infinity
                                for i in ceiling)
            elif callable(ceiling):
                pass
            else:
                raise ValueError("Unable to parse value of parameter ceiling")
            self.ceiling = Envelope(ceiling, sign=1,
                    min_part=self.min_part, max_part=self.max_part,
                    min_slope=self.min_slope, max_slope=self.max_slope,
                    min_length=self.min_length)

    def __richcmp__(self, other, int op):
        r"""
        Basic comparison function, supporting only checking for
        equality.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3).backend
            sage: D = IntegerListsLex(2, length=3).backend; L = list(D._iter())
            sage: E = IntegerListsLex(2, min_length=3).backend
            sage: G = IntegerListsLex(4, length=3).backend
            sage: C >= C
            True
            sage: C == D
            True
            sage: C != D
            False
            sage: C == E
            False
            sage: C != E
            True
            sage: C == None
            False
            sage: C == G
            False
            sage: C <= G
            Traceback (most recent call last):
            ...
            TypeError: IntegerListsBackend can only be compared for equality
        """
        cdef IntegerListsBackend left = <IntegerListsBackend>self
        cdef IntegerListsBackend right = <IntegerListsBackend>other
        equal = (type(left) is type(other) and
            left.min_length == right.min_length and
            left.max_length == right.max_length and
            left.min_sum == right.min_sum and
            left.max_sum == right.max_sum and
            left.min_slope == right.min_slope and
            left.max_slope == right.max_slope and
            left.floor == right.floor and
            left.ceiling == right.ceiling)
        if equal:
            return op == Py_EQ or op == Py_LE or op == Py_GE
        if op == Py_EQ:
            return False
        if op == Py_NE:
            return True
        else:
            raise TypeError("IntegerListsBackend can only be compared for equality")

    def _repr_(self):
        """
        Return the name of this enumerated set.

        EXAMPLES::

            sage: from sage.combinat.integer_lists.base import IntegerListsBackend
            sage: C = IntegerListsBackend(2, length=3)
            sage: C._repr_()
            'Integer lists of sum 2 satisfying certain constraints'
        """
        if self.min_sum == self.max_sum:
            return "Integer lists of sum {} satisfying certain constraints".format(self.min_sum)
        elif self.max_sum == Infinity:
            if self.min_sum == 0:
                return "Integer lists with arbitrary sum satisfying certain constraints"
            else:
                return "Integer lists of sum at least {} satisfying certain constraints".format(self.min_sum)
        else:
            return "Integer lists of sum between {} and {} satisfying certain constraints".format(self.min_sum, self.max_sum)

    def _contains(self, comp):
        """
        Return ``True`` if ``comp`` meets the constraints imposed
        by the arguments.

        EXAMPLES::

            sage: C = IntegerListsLex(n=2, max_length=3, min_slope=0)
            sage: all(l in C for l in C)  # indirect doctest
            True
        """
        if len(comp) < self.min_length or len(comp) > self.max_length:
            return False
        n = sum(comp)
        if n < self.min_sum or n > self.max_sum:
            return False
        for i in range(len(comp)):
            if comp[i] < self.floor(i):
                return False
            if comp[i] > self.ceiling(i):
                return False
        for i in range(len(comp)-1):
            slope = comp[i+1] - comp[i]
            if slope < self.min_slope or slope > self.max_slope:
                return False
        return True

    def __getstate__(self):
        """
        Pickle ``self``.

        EXAMPLES::

            sage: from sage.combinat.integer_lists.base import IntegerListsBackend
            sage: C = IntegerListsBackend(2, length=3)
            sage: C.__getstate__()
            {'ceiling': <sage.combinat.integer_lists.base.Envelope object at ...>,
             'floor': <sage.combinat.integer_lists.base.Envelope object at ...>,
             'max_length': 3,
             'max_part': inf,
             'max_slope': inf,
             'max_sum': 2,
             'min_length': 3,
             'min_part': 0,
             'min_slope': -inf,
             'min_sum': 2}
        """
        return {"min_sum": self.min_sum,
                "max_sum": self.max_sum,
                "min_length": self.min_length,
                "max_length": self.max_length,
                "min_part": self.min_part,
                "max_part": self.max_part,
                "min_slope": self.min_slope,
                "max_slope": self.max_slope,
                "floor": self.floor,
                "ceiling": self.ceiling}

    def __setstate__(self, state):
        """
        Unpickle ``self`` from the state ``state``.

        EXAMPLES::

            sage: from sage.combinat.integer_lists.base import IntegerListsBackend
            sage: C = IntegerListsBackend(2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: C == loads(dumps(C)) # this did fail at some point, really!
            True
            sage: C is loads(dumps(C)) # todo: not implemented
            True
        """
        self.__init__(**state)


cdef class Envelope():
    r"""
    The (currently approximated) upper (lower) envelope of a function
    under the specified constraints.

    INPUT:

    - ``f`` -- a function, list, or tuple; if ``f`` is a list, it is
      considered as the function ``f(i)=f[i]``, completed for larger
      `i` with ``f(i)=max_part``.

    - ``min_part``, ``max_part``, ``min_slope``, ``max_slope``, ...
      as for :class:`IntegerListsLex` (please consult for details).

    - ``sign`` -- (+1 or -1) multiply the input values with ``sign``
      and multiply the output with ``sign``. Setting this to `-1` can
      be used to implement a lower envelope.

    The *upper envelope* `U(f)` of `f` is the (pointwise) largest
    function which is bounded above by `f` and satisfies the
    ``max_part`` and ``max_slope`` conditions. Furthermore, for
    ``i,i+1<min_length``, the upper envelope also satisfies the
    ``min_slope`` condition.

    Upon computing `U(f)(i)`, all the previous values
    for `j\leq i` are computed and cached; in particular `f(i)` will
    be computed at most once for each `i`.

    .. TODO::

        - This class is a good candidate for Cythonization, especially
          to get the critical path in ``__call__`` super fast.

        - To get full envelopes, we would want both the ``min_slope``
          and ``max_slope`` conditions to always be satisfied. This is
          only properly defined for the restriction of `f` to a finite
          interval `0,..,k`, and depends on `k`.

        - This is the core "data structure" of
          ``IntegerListsLex``. Improving the lookahead there
          essentially depends on having functions with a good
          complexity to compute the area below an envelope; and in
          particular how it evolves when increasing the length.

    EXAMPLES::

        sage: from sage.combinat.integer_lists import Envelope

    Trivial upper and lower envelopes::

        sage: f = Envelope([3,2,2])
        sage: [f(i) for i in range(10)]
        [3, 2, 2, inf, inf, inf, inf, inf, inf, inf]
        sage: f = Envelope([3,2,2], sign=-1)
        sage: [f(i) for i in range(10)]
        [3, 2, 2, 0, 0, 0, 0, 0, 0, 0]

    A more interesting lower envelope::

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1)
        sage: [f(i) for i in range(10)]
        [4, 3, 5, 4, 5, 4, 3, 2, 2, 2]

    Currently, adding ``max_slope`` has no effect::

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0)
        sage: [f(i) for i in range(10)]
        [4, 3, 5, 4, 5, 4, 3, 2, 2, 2]

    unless ``min_length`` is large enough::

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0, min_length=2)
        sage: [f(i) for i in range(10)]
        [4, 3, 5, 4, 5, 4, 3, 2, 2, 2]

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0, min_length=4)
        sage: [f(i) for i in range(10)]
        [5, 5, 5, 4, 5, 4, 3, 2, 2, 2]

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0, min_length=5)
        sage: [f(i) for i in range(10)]
        [5, 5, 5, 5, 5, 4, 3, 2, 2, 2]

    A non trivial upper envelope::

        sage: f = Envelope([9,1,5,4], max_part=7, max_slope=2)
        sage: [f(i) for i in range(10)]
        [7, 1, 3, 4, 6, 7, 7, 7, 7, 7]

    TESTS::

        sage: f = Envelope(3, min_slope=1)
        sage: [f(i) for i in range(10)]
        [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

        sage: f = Envelope(3, min_slope=1, min_length=5)
        sage: [f(i) for i in range(10)]
        [-1, 0, 1, 2, 3, 3, 3, 3, 3, 3]

        sage: f = Envelope(3, sign=-1, min_slope=1)
        sage: [f(i) for i in range(10)]
        [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

        sage: f = Envelope(3, sign=-1, max_slope=-1, min_length=4)
        sage: [f(i) for i in range(10)]
        [6, 5, 4, 3, 3, 3, 3, 3, 3, 3]
    """
    def __init__(self, f, *,
                 min_part=0, max_part=Infinity,
                 min_slope=MInfinity, max_slope=Infinity,
                 min_length=0, max_length=Infinity, sign=1):
        r"""
        Initialize this envelope.

        TESTS::

            sage: from sage.combinat.integer_lists import Envelope
            sage: f = Envelope(3, sign=-1, max_slope=-1, min_length=4)
            sage: f.sign
            -1
            sage: f.max_part
            -3
            sage: f.max_slope
            inf
            sage: f.min_slope
            1
            sage: TestSuite(f).run(skip='_test_pickling')
            sage: Envelope(3, sign=1/3, max_slope=-1, min_length=4)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: Envelope(3, sign=-2, max_slope=-1, min_length=4)
            Traceback (most recent call last):
            ...
            ValueError: sign should be +1 or -1
        """
        # self.sign = sign for the output values (the sign change for
        # f is handled here in __init__)
        self.sign = Integer(sign)
        if self.sign == 1:
            self.max_part = max_part
            self.min_slope = min_slope
            self.max_slope = max_slope
            if max_part == 0:
                # This uses that all entries are nonnegative.
                # This is not for speed optimization but for
                # setting the limit start and avoid hangs.
                # See #17979: comment 389
                f = Integer(0)
        elif self.sign == -1:
            self.max_part = -min_part
            self.min_slope = -max_slope
            self.max_slope = -min_slope
        else:
            raise ValueError("sign should be +1 or -1")

        # Handle different types of f and multiply f with sign
        if isinstance(f, RingElement) or f == Infinity or f == -Infinity:
            limit_start = 0
            self.max_part = min(self.sign * f, self.max_part)
            f = ConstantFunction(Infinity)
        elif isinstance(f, (list, tuple)):
            limit_start = len(f)
            f_tab = [self.sign * i for i in f]
            f = lambda k: f_tab[k] if k < len(f_tab) else Infinity
        else:
            g = f
            f = lambda k: self.sign * g(k)
            # At this point, this is not really used
            limit_start = Infinity

        self.f = f
        # For i >= limit_start, f is constant
        # This does not necessarily means that self is constant!
        self.f_limit_start = limit_start
        self.precomputed = []

        if min_length > 0:
            self(min_length-1)
            for i in range(min_length-1,0,-1):
                self.precomputed[i-1] = min(self.precomputed[i-1], self.precomputed[i] - self.min_slope)

    def __richcmp__(self, other, int op):
        r"""
        Basic comparison function, supporting only checking for
        equality.

        EXAMPLES::

            sage: from sage.combinat.integer_lists import Envelope
            sage: f = Envelope([3,2,2])
            sage: g = Envelope([3,2,2])
            sage: h = Envelope([3,2,2], min_part=2)
            sage: f == f, f == h, f == None
            (True, False, False)
            sage: f < f, f != h, f != None
            (False, True, True)

        This would be desirable::

            sage: f == g          # todo: not implemented
            True
        """
        cdef Envelope left = <Envelope>self
        cdef Envelope right = <Envelope>other
        equal = (type(left) is type(other) and
            left.sign == right.sign and
            left.f == right.f and
            left.f_limit_start == right.f_limit_start and
            left.max_part == right.max_part and
            left.min_slope == right.min_slope and
            left.max_slope == right.max_slope)
        if equal:
            return op == Py_EQ or op == Py_LE or op == Py_GE
        if op == Py_EQ:
            return False
        if op == Py_NE:
            return True
        else:
            raise TypeError("Envelopes can only be compared for equality")

    def limit_start(self):
        """
        Return from which `i` on the bound returned by ``limit`` holds.

        .. SEEALSO:: :meth:`limit` for the precise specifications.

        EXAMPLES::

            sage: from sage.combinat.integer_lists import Envelope
            sage: Envelope([4,1,5]).limit_start()
            3
            sage: Envelope([4,1,5], sign=-1).limit_start()
            3

            sage: Envelope([4,1,5], max_part=2).limit_start()
            3

            sage: Envelope(4).limit_start()
            0
            sage: Envelope(4, sign=-1).limit_start()
            0

            sage: Envelope(lambda x: 3).limit_start() == Infinity
            True
            sage: Envelope(lambda x: 3, max_part=2).limit_start() == Infinity
            True

            sage: Envelope(lambda x: 3, sign=-1, min_part=2).limit_start() == Infinity
            True
        """
        return self.f_limit_start

    def limit(self):
        r"""
        Return a bound on the limit of ``self``.

        OUTPUT: nonnegative integer or `\infty`

        This returns some upper bound for the accumulation points of
        this upper envelope. For a lower envelope, a lower bound is
        returned instead.

        In particular this gives a bound for the value of ``self`` at
        `i` for `i` large enough. Special case: for a lower envelop,
        and when the limit is `\infty`, the envelope is guaranteed to
        tend to `\infty` instead.

        When ``s=self.limit_start()`` is finite, this bound is
        guaranteed to be valid for `i>=s`.

        Sometimes it's better to have a loose bound that starts early;
        sometimes the converse holds. At this point which specific
        bound and starting point is returned is not set in stone, in
        order to leave room for later optimizations.

        EXAMPLES::

            sage: from sage.combinat.integer_lists import Envelope
            sage: Envelope([4,1,5]).limit()
            inf
            sage: Envelope([4,1,5], max_part=2).limit()
            2
            sage: Envelope([4,1,5], max_slope=0).limit()
            1
            sage: Envelope(lambda x: 3, max_part=2).limit()
            2

        Lower envelopes::

            sage: Envelope(lambda x: 3, min_part=2, sign=-1).limit()
            2
            sage: Envelope([4,1,5], min_slope=0, sign=-1).limit()
            5
            sage: Envelope([4,1,5], sign=-1).limit()
            0

        .. SEEALSO:: :meth:`limit_start`
        """
        if self.limit_start() < Infinity and self.max_slope <= 0:
            return self(self.limit_start())
        else:
            return self.max_part * self.sign

    def __call__(self, Py_ssize_t k):
        """
        Return the value of this envelope at `k`.

        EXAMPLES::

            sage: from sage.combinat.integer_lists import Envelope
            sage: f = Envelope([4,1,5,3,5])
            sage: f.__call__(2)
            5
            sage: [f(i) for i in range(10)]
            [4, 1, 5, 3, 5, inf, inf, inf, inf, inf]

        .. NOTE::

            See the documentation of :class:`Envelope` for tests and
            examples.
        """
        if k >= len(self.precomputed):
            for i in range(len(self.precomputed), k+1):
                value = min(self.f(i), self.max_part)
                if i > 0:
                    value = min(value, self.precomputed[i-1] + self.max_slope)
                self.precomputed.append(value)
        return self.precomputed[k] * self.sign

    def adapt(self, m, j):
        """
        Return this envelope adapted to an additional local constraint.

        INPUT:

        - ``m`` -- nonnegative integer (starting value)

        - ``j`` -- nonnegative integer (position)

        This method adapts this envelope to the additional local
        constraint imposed by having a part `m` at position `j`.
        Namely, this returns a function which computes, for any `i>j`,
        the minimum of the ceiling function and the value restriction
        given by the slope conditions.

        EXAMPLES::

            sage: from sage.combinat.integer_lists import Envelope
            sage: f = Envelope(3)
            sage: g = f.adapt(1,1)
            sage: g is f
            True
            sage: [g(i) for i in range(10)]
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

            sage: f = Envelope(3, max_slope=1)
            sage: g = f.adapt(1,1)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

        Note that, in both cases above, the adapted envelope is only
        guaranteed to be valid for `i>j`! This is to leave potential
        room in the future for sharing similar adapted envelopes::

            sage: g = f.adapt(0,0)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

            sage: g = f.adapt(2,2)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

            sage: g = f.adapt(3,3)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

        Now with a lower envelope::

            sage: f = Envelope(1, sign=-1, min_slope=-1)
            sage: g = f.adapt(2,2)
            sage: [g(i) for i in range(10)]
            [4, 3, 2, 1, 1, 1, 1, 1, 1, 1]
            sage: g = f.adapt(1,3)
            sage: [g(i) for i in range(10)]
            [4, 3, 2, 1, 1, 1, 1, 1, 1, 1]
        """
        if self.max_slope == Infinity:
            return self
        m *= self.sign
        m = m - j * self.max_slope
        return lambda i: self.sign * min(m + i*self.max_slope, self.sign*self(i) )

    def __reduce__(self):
        """
        Pickle ``self``.

        EXAMPLES::

            sage: from sage.combinat.integer_lists import Envelope
            sage: h = Envelope(3, min_part=2)
            sage: loads(dumps(h)) == h
            True
        """
        args = (type(self),
                self.sign, self.f, self.f_limit_start, self.precomputed,
                self.max_part, self.min_slope, self.max_slope)
        return _unpickle_Envelope, args


def _unpickle_Envelope(type t, _sign, _f, _f_limit_start, _precomputed,
        _max_part, _min_slope, _max_slope):
    """
    Internal function to support pickling for :class:`Envelope`.

    EXAMPLES::

        sage: from sage.combinat.integer_lists.base import Envelope, _unpickle_Envelope
        sage: _unpickle_Envelope(Envelope,
        ....:     1, lambda i:i, Infinity, [], 4, -1, 3)
        <sage.combinat.integer_lists.base.Envelope object at ...>
    """
    cdef Envelope self = t.__new__(t)
    self.sign = _sign
    self.f = _f
    self.f_limit_start = _f_limit_start
    self.precomputed = _precomputed
    self.max_part = _max_part
    self.min_slope = _min_slope
    self.max_slope = _max_slope
    return self
