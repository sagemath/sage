"""
Lazy Power Series

This file provides an implementation of lazy univariate power
series, which uses the stream class for its internal data
structure. The lazy power series keep track of their approximate
order as much as possible without forcing the computation of any
additional coefficients. This is required for recursively defined
power series.

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevant section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatse9.html.
"""
# ****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
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

import builtins

from .stream import Stream, Stream_class
from .series_order import  bounded_decrement, increment, inf, unk
from sage.rings.integer import Integer
from sage.misc.misc_c import prod
from functools import partial
from sage.misc.misc import is_iterator
from sage.misc.repr import repr_lincomb
from sage.misc.cachefunc import cached_method

from sage.rings.ring import Algebra
from sage.structure.parent import Parent
from sage.categories.all import Rings
from sage.structure.element import Element, parent, AlgebraElement


class LazyPowerSeriesRing(Algebra):
    def __init__(self, R, names=None, element_class=None):
        """
        TESTS::

            sage: from sage.combinat.species.series import LazyPowerSeriesRing
            sage: L = LazyPowerSeriesRing(QQ)

        Equality testing is undecidable in general, and not much
        efforts are done at this stage to implement equality when
        possible. Hence the failing tests below::

            sage: TestSuite(L).run()
            Failure in ...
            The following tests failed: _test_additive_associativity, _test_associativity, _test_distributivity, _test_elements, _test_one, _test_prod, _test_zero

        ::

            sage: LazyPowerSeriesRing(QQ, 'z').gen()
            z
            sage: LazyPowerSeriesRing(QQ, ['z']).gen()
            z
            sage: LazyPowerSeriesRing(QQ, ['x', 'z'])
            Traceback (most recent call last):
            ...
            NotImplementedError: only univariate lazy power series rings are supported
        """
        #Make sure R is a ring with unit element
        if R not in Rings():
            raise TypeError("Argument R must be a ring.")

        #Take care of the names
        if names is None:
            names = 'x'
        elif isinstance(names, (list, tuple)):
            if len(names) != 1:
                raise NotImplementedError(
                    'only univariate lazy power series rings are supported')
            names = names[0]
        else:
            names = str(names)

        self._element_class = element_class if element_class is not None else LazyPowerSeries
        self._order = None
        self._name = names
        self._zero_base_ring = R.zero()
        Parent.__init__(self, R, category=Rings())

    def ngens(self):
        """
        EXAMPLES::

            sage: LazyPowerSeriesRing(QQ).ngens()
            1
        """
        return 1

    def __repr__(self):
        """
        EXAMPLES::

            sage: LazyPowerSeriesRing(QQ)
            Lazy Power Series Ring over Rational Field
        """
        return "Lazy Power Series Ring over %s" % self.base_ring()

    def __eq__(self, x):
        """
        Check whether ``self`` is equal to ``x``.

        EXAMPLES::

            sage: LQ = LazyPowerSeriesRing(QQ)
            sage: LZ = LazyPowerSeriesRing(ZZ)
            sage: LQ == LQ
            True
            sage: LZ == LQ
            False
        """
        if not isinstance(x, LazyPowerSeriesRing):
            return False
        return self.base_ring() == x.base_ring()

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: LQ = LazyPowerSeriesRing(QQ)
            sage: LZ = LazyPowerSeriesRing(ZZ)
            sage: LQ != LQ
            False
            sage: LZ != LQ
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: LQ = LazyPowerSeriesRing(QQ)
            sage: LZ = LazyPowerSeriesRing(ZZ)
            sage: hash(LQ) == hash(LQ)
            True
            sage: hash(LZ) == hash(LQ)
            False
        """
        # with a random number, so that the hash is not that of the base ring
        return hash((16079305, self.base_ring()))

    def _coerce_impl(self, x):
        """
        EXAMPLES::

            sage: L1 = LazyPowerSeriesRing(QQ)
            sage: L2 = LazyPowerSeriesRing(RR)
            sage: L2.has_coerce_map_from(L1)
            True
            sage: L1.has_coerce_map_from(L2)
            False

        ::

            sage: a = L1([1]) + L2([1])
            sage: a.coefficients(3)
            [2.00000000000000, 2.00000000000000, 2.00000000000000]
        """
        return self(x)

    def __call__(self, x=None, order=unk):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: L()
            Uninitialized lazy power series
            sage: L(1)
            1
            sage: L(ZZ).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]
            sage: L(iter(ZZ)).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]
            sage: L(Stream(ZZ)).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]

        ::

            sage: a = L([1,2,3])
            sage: a.coefficients(3)
            [1, 2, 3]
            sage: L(a) is a
            True
            sage: L_RR = LazyPowerSeriesRing(RR)
            sage: b = L_RR(a)
            sage: b.coefficients(3)
            [1.00000000000000, 2.00000000000000, 3.00000000000000]
            sage: L(b)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to coerce ... into self

        TESTS::

            sage: L(pi)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to coerce pi into self
        """
        cls = self._element_class
        BR = self.base_ring()

        if x is None:
            res = cls(self, stream=None, order=unk, aorder=unk,
                      aorder_changed=True, is_initialized=False)
            res.compute_aorder = uninitialized
            return res

        if isinstance(x, LazyPowerSeries):
            x_parent = x.parent()
            if x_parent.__class__ != self.__class__:
                raise ValueError

            if x_parent.base_ring() == self.base_ring():
                return x
            else:
                if self.base_ring().has_coerce_map_from(x_parent.base_ring()):
                    return x._new(partial(x._change_ring_gen, self.base_ring()), lambda ao: ao, x, parent=self)


        if BR.has_coerce_map_from(parent(x)):
            x = BR(x)
            return self.term(x, 0)

        if hasattr(x, "__iter__") and not isinstance(x, Stream_class):
            x = iter(x)

        if is_iterator(x):
            x = Stream(x)

        if isinstance(x, Stream_class):
            aorder = order if order != unk else 0
            return cls(self, stream=x, order=order, aorder=aorder,
                       aorder_changed=False, is_initialized=True)
        elif not isinstance(x, Element):
            x = BR(x)
            return self.term(x, 0)

        raise TypeError("do not know how to coerce %s into self" % x)

    @cached_method
    def zero(self):
        """
        Return the zero power series.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.zero()
            0
        """
        return self.term(self._zero_base_ring, 0)

    def identity_element(self):
        """
        Return the one power series.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.identity_element()
            1
        """
        return self(self.base_ring()(1))

    def gen(self, i=0):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.gen().coefficients(5)
            [0, 1, 0, 0, 0]
        """
        res = self._new_initial(1, Stream([0,1,0]))
        res._name = self._name
        return res

    def term(self, r, n):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.term(0,0)
            0
            sage: L.term(3,2).coefficients(5)
            [0, 0, 3, 0, 0]
        """
        if n < 0:
            raise ValueError("n must be non-negative")
        BR = self.base_ring()
        if r == 0:
            res = self._new_initial(inf, Stream([0]))
            res._name = "0"
        else:
            zero = BR(0)
            s = [zero]*n+[BR(r),zero]
            res = self._new_initial(n, Stream(s))

            if n == 0:
                res._name = repr(r)
            elif n == 1:
                res._name = repr(r) + "*" + self._name
            else:
                res._name = "%s*%s^%s" % (repr(r), self._name, n)

        return res

    def _new_initial(self, order, stream):
        """
        Return a new power series with specified order.

        INPUT:


        -  ``order`` - a non-negative integer

        -  ``stream`` - a Stream object


        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: L._new_initial(0, Stream([1,2,3,0])).coefficients(5)
            [1, 2, 3, 0, 0]
        """
        return self._element_class(self, stream=stream, order=order, aorder=order,
                                   aorder_changed=False, is_initialized=True)


    def _sum_gen(self, series_list):
        """
        Return a generator for the coefficients of the sum of the lazy
        power series in series_list.

        INPUT:


        -  ``series_list`` - a list of lazy power series


        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: series_list = [ L([1]), L([0,1]), L([0,0,1]) ]
            sage: g = L._sum_gen(series_list)
            sage: [next(g) for i in range(5)]
            [1, 2, 3, 3, 3]
        """
        last_index = len(series_list) - 1
        assert last_index >= 0
        n = 0
        while True:
            r = sum( [f.coefficient(n) for f in series_list] )
            yield r
            n += 1

    def sum(self, a):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: l = [L(ZZ)]*3
            sage: L.sum(l).coefficients(10)
            [0, 3, -3, 6, -6, 9, -9, 12, -12, 15]
        """
        return self( self._sum_gen(a) )

    #Potentially infinite sum
    def _sum_generator_gen(self, g):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([1])
            sage: def f():
            ....:     while True:
            ....:         yield s
            sage: g = L._sum_generator_gen(f())
            sage: [next(g) for i in range(10)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        s = Stream(g)
        n = 0
        while True:
            r = s[n].coefficient(n)
            for i in range(len(s)-1):
                r += s[i].coefficient(n)
            yield r
            n += 1

    def sum_generator(self, g):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: g = [L([1])]*6 + [L(0)]
            sage: t = L.sum_generator(g)
            sage: t.coefficients(10)
            [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]

        ::

            sage: s = L([1])
            sage: def g():
            ....:     while True:
            ....:         yield s
            sage: t = L.sum_generator(g())
            sage: t.coefficients(9)
            [1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return self(self._sum_generator_gen(g))

    #Potentially infinite product
    def _product_generator_gen(self, g):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import _integers_from
            sage: L = LazyPowerSeriesRing(QQ)
            sage: g = (L([1]+[0]*i+[1]) for i in _integers_from(0))
            sage: g2 = L._product_generator_gen(g)
            sage: [next(g2) for i in range(10)]
            [1, 1, 2, 4, 7, 12, 20, 33, 53, 84]
        """
        z = next(g)
        yield z.coefficient(0)
        yield z.coefficient(1)

        n = 2

        for x in g:
            z = z * x
            yield z.coefficient(n)
            n += 1

        while True:
            yield z.coefficient(n)
            n += 1

    def product_generator(self, g):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s1 = L([1,1,0])
            sage: s2 = L([1,0,1,0])
            sage: s3 = L([1,0,0,1,0])
            sage: s4 = L([1,0,0,0,1,0])
            sage: s5 = L([1,0,0,0,0,1,0])
            sage: s6 = L([1,0,0,0,0,0,1,0])
            sage: s = [s1, s2, s3, s4, s5, s6]
            sage: def g():
            ....:     for a in s:
            ....:         yield a
            sage: p = L.product_generator(g())
            sage: p.coefficients(26)
            [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 4, 4, 4, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0]

        ::

            sage: def m(n):
            ....:     yield 1
            ....:     while True:
            ....:         for i in range(n-1):
            ....:             yield 0
            ....:         yield 1
            sage: def s(n):
            ....:     q = 1/n
            ....:     yield 0
            ....:     while True:
            ....:         for i in range(n-1):
            ....:             yield 0
            ....:         yield q

        ::

            sage: def lhs_gen():
            ....:     n = 1
            ....:     while True:
            ....:         yield L(m(n))
            ....:         n += 1

        ::

            sage: def rhs_gen():
            ....:     n = 1
            ....:     while True:
            ....:         yield L(s(n))
            ....:         n += 1
            sage: lhs = L.product_generator(lhs_gen())
            sage: rhs = L.sum_generator(rhs_gen()).exponential()
            sage: lhs.coefficients(10)
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
            sage: rhs.coefficients(10)
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        return self(self._product_generator_gen(g))



class LazyPowerSeries(AlgebraElement):
    def __init__(self, A, stream=None, order=None, aorder=None, aorder_changed=True, is_initialized=False, name=None):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L()
            sage: loads(dumps(f))
            Uninitialized lazy power series
        """
        AlgebraElement.__init__(self, A)
        self._stream = stream
        self.order = unk if order is None else order
        self.aorder = unk if aorder is None else aorder
        if self.aorder == inf:
            self.order = inf
        self.aorder_changed = aorder_changed
        self.is_initialized = is_initialized
        self._name = name

    def compute_aorder(*args, **kwargs):
        """
        The default compute_aorder does nothing.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L(1)
            sage: a.compute_aorder() is None
            True
        """
        return None

    def _get_repr_info(self, x):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1,2,3])
            sage: a.compute_coefficients(5)
            sage: a._get_repr_info('x')
            [('1', 1), ('x', 2), ('x^2', 3)]
        """
        n = len(self._stream)
        m = ['1', x]
        m += [x + "^" + str(i) for i in range(2, n)]
        c = [self._stream[i] for i in range(n)]
        return [(mo, co) for mo, co in zip(m, c) if co != 0]

    def __repr__(self):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L(); s._name = 's'; s
            s

        ::

            sage: L()
            Uninitialized lazy power series

        ::

            sage: a = L([1,2,3])
            sage: a
            O(1)
            sage: a.compute_coefficients(2)
            sage: a
            1 + 2*x + 3*x^2 + O(x^3)
            sage: a.compute_coefficients(4)
            sage: a
            1 + 2*x + 3*x^2 + 3*x^3 + 3*x^4 + 3*x^5 + ...

        ::

            sage: a = L([1,2,3,0])
            sage: a.compute_coefficients(5)
            sage: a
            1 + 2*x + 3*x^2
        """
        if self._name is not None:
            return self._name

        if self.is_initialized:
            n = len(self._stream)
            x = self.parent()._name
            baserepr = repr_lincomb(self._get_repr_info(x))
            if self._stream.is_constant():
                if self._stream[n-1] == 0:
                    l = baserepr
                else:
                    l = baserepr + " + " + repr_lincomb([(x+"^"+str(i), self._stream[n-1]) for i in range(n, n+3)]) + " + ..."
            else:
                l = baserepr + " + O(x^%s)" % n if n > 0 else "O(1)"
        else:
            l = 'Uninitialized lazy power series'
        return l


    def refine_aorder(self):
        """
        Refines the approximate order of self as much as possible without
        computing any coefficients.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([0,0,0,0,1])
            sage: a.aorder
            0
            sage: a.coefficient(2)
            0
            sage: a.aorder
            0
            sage: a.refine_aorder()
            sage: a.aorder
            3

        ::

            sage: a = L([0,0])
            sage: a.aorder
            0
            sage: a.coefficient(5)
            0
            sage: a.refine_aorder()
            sage: a.aorder
            Infinite series order

        ::

            sage: a = L([0,0,1,0,0,0])
            sage: a[4]
            0
            sage: a.refine_aorder()
            sage: a.aorder
            2
        """
        #If we already know the order, then we don't have
        #to worry about the approximate order
        if self.order != unk:
            return

        #aorder can never be infinity since order would have to
        #be infinity as well
        assert self.aorder != inf

        if self.aorder == unk or not self.is_initialized:
            self.compute_aorder()
        else:
            #Try to improve the approximate order
            ao = self.aorder
            c = self._stream
            n = c.number_computed()


            if ao == 0 and n > 0:
                while ao < n:
                    if self._stream[ao] == 0:
                        self.aorder += 1
                        ao += 1
                    else:
                        break

            #Try to recognize the zero series
            if ao == n:
                #For non-constant series, we cannot do anything
                if not c.is_constant():
                    return
                if c[n-1] == 0:
                    self.aorder = inf
                    self.order  = inf
                    return

            if ao < n:
                self.order = ao


        if hasattr(self, '_reference') and self._reference is not None:
            self._reference._copy(self)

    def initialize_coefficient_stream(self, compute_coefficients):
        """
        Initializes the coefficient stream.

        INPUT: compute_coefficients

        TESTS::

            sage: from sage.combinat.species.series_order import inf, unk
            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L()
            sage: compute_coefficients = lambda ao: iter(ZZ)
            sage: f.order = inf
            sage: f.aorder = inf
            sage: f.initialize_coefficient_stream(compute_coefficients)
            sage: f.coefficients(5)
            [0, 0, 0, 0, 0]

        ::

            sage: f = L()
            sage: compute_coefficients = lambda ao: iter(ZZ)
            sage: f.order = 1
            sage: f.aorder = 1
            sage: f.initialize_coefficient_stream(compute_coefficients)
            sage: f.coefficients(5)
            [0, 1, -1, 2, -2]
        """
        ao = self.aorder
        assert ao != unk

        if ao == inf:
            self.order = inf
            self._stream = Stream(0)
        else:
            self._stream = Stream(compute_coefficients(ao))

        self.is_initialized = True

    def compute_coefficients(self, i):
        """
        Computes all the coefficients of self up to i.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1,2,3])
            sage: a.compute_coefficients(5)
            sage: a
            1 + 2*x + 3*x^2 + 3*x^3 + 3*x^4 + 3*x^5 + ...
        """
        self.coefficient(i)

    def coefficients(self, n):
        """
        Return the first n coefficients of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,2,3,0])
            sage: f.coefficients(5)
            [1, 2, 3, 0, 0]
        """
        return [self.coefficient(i) for i in range(n)]

    def is_zero(self):
        """
        Return True if and only if self is zero.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([0,2,3,0])
            sage: s.is_zero()
            False

        ::

            sage: s = L(0)
            sage: s.is_zero()
            True

        ::

            sage: s = L([0])
            sage: s.is_zero()
            False
            sage: s.coefficient(0)
            0
            sage: s.coefficient(1)
            0
            sage: s.is_zero()
            True
        """
        self.refine_aorder()
        return self.order == inf

    def set_approximate_order(self, new_order):
        """
        Sets the approximate order of self and returns True if the
        approximate order has changed otherwise it will return False.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([0,0,0,3,2,1,0])
            sage: f.get_aorder()
            0
            sage: f.set_approximate_order(3)
            True
            sage: f.set_approximate_order(3)
            False
        """
        self.aorder_changed = ( self.aorder != new_order )
        self.aorder = new_order
        return self.aorder_changed

    def _copy(self, x):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L.term(2, 2)
            sage: g = L()
            sage: g._copy(f)
            sage: g.order
            2
            sage: g.aorder
            2
            sage: g.is_initialized
            True
            sage: g.coefficients(4)
            [0, 0, 2, 0]
        """
        self.order  = x.order
        self.aorder = x.aorder
        self.aorder_changed = x.aorder_changed
        self.compute_aorder = x.compute_aorder
        self.is_initialized = x.is_initialized
        self._stream = x._stream

    def define(self, x):
        """
        EXAMPLES: Test Recursive 0

        ::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: one = L(1)
            sage: monom = L.gen()
            sage: s = L()
            sage: s._name = 's'
            sage: s.define(one+monom*s)
            sage: s.aorder
            0
            sage: s.order
            Unknown series order
            sage: [s.coefficient(i) for i in range(6)]
            [1, 1, 1, 1, 1, 1]

        Test Recursive 1

        ::

            sage: s = L()
            sage: s._name = 's'
            sage: s.define(one+monom*s*s)
            sage: s.aorder
            0
            sage: s.order
            Unknown series order
            sage: [s.coefficient(i) for i in range(6)]
            [1, 1, 2, 5, 14, 42]

        Test Recursive 1b

        ::

            sage: s = L()
            sage: s._name = 's'
            sage: s.define(monom + s*s)
            sage: s.aorder
            1
            sage: s.order
            Unknown series order
            sage: [s.coefficient(i) for i in range(7)]
            [0, 1, 1, 2, 5, 14, 42]

        Test Recursive 2

        ::

            sage: s = L()
            sage: s._name = 's'
            sage: t = L()
            sage: t._name = 't'
            sage: s.define(one+monom*t*t*t)
            sage: t.define(one+monom*s*s)
            sage: [s.coefficient(i) for i in range(9)]
            [1, 1, 3, 9, 34, 132, 546, 2327, 10191]
            sage: [t.coefficient(i) for i in range(9)]
            [1, 1, 2, 7, 24, 95, 386, 1641, 7150]

        Test Recursive 2b

        ::

            sage: s = L()
            sage: s._name = 's'
            sage: t = L()
            sage: t._name = 't'
            sage: s.define(monom + t*t*t)
            sage: t.define(monom + s*s)
            sage: [s.coefficient(i) for i in range(9)]
            [0, 1, 0, 1, 3, 3, 7, 30, 63]
            sage: [t.coefficient(i) for i in range(9)]
            [0, 1, 1, 0, 2, 6, 7, 20, 75]

        Test Recursive 3

        ::

            sage: s = L()
            sage: s._name = 's'
            sage: s.define(one+monom*s*s*s)
            sage: [s.coefficient(i) for i in range(10)]
            [1, 1, 3, 12, 55, 273, 1428, 7752, 43263, 246675]
        """
        self._copy(x)
        x._reference = self

    def coefficient(self, n):
        """
        Return the coefficient of xn in self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(ZZ)
            sage: [f.coefficient(i) for i in range(5)]
            [0, 1, -1, 2, -2]
        """
        # The following line must not be written n < self.get_aorder()
        # because comparison of Integer and OnfinityOrder is not implemented.
        if self.get_aorder() > n:
            return self.parent()._zero_base_ring

        assert self.is_initialized

        return self._stream[n]

    def get_aorder(self):
        """
        Return the approximate order of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L.gen()
            sage: a.get_aorder()
            1
        """
        self.refine_aorder()
        return self.aorder

    def get_order(self):
        """
        Return the order of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L.gen()
            sage: a.get_order()
            1
        """
        self.refine_aorder()
        return self.order

    def get_stream(self):
        """
        Return self's underlying Stream object.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L.gen()
            sage: s = a.get_stream()
            sage: [s[i] for i in range(5)]
            [0, 1, 0, 0, 0]
        """
        self.refine_aorder()
        return self._stream

    def _approximate_order(self, compute_coefficients, new_order, *series):
        if self.is_initialized:
            return

        ochanged = self.aorder_changed

        ao = new_order(*[s.aorder for s in series])
        ao = inf if ao == unk else ao

        tchanged = self.set_approximate_order(ao)

        if len(series) == 0:
            must_initialize_coefficient_stream = True
            tchanged = ochanged = False
        elif len(series) == 1 or len(series) == 2:
            must_initialize_coefficient_stream = ( self.aorder == unk or self.is_initialized is False)
        else:
            raise ValueError

        if ochanged or tchanged:
            for s in series:
                s.compute_aorder()
            ao = new_order(*[s.aorder for s in series])
            tchanged = self.set_approximate_order(ao)

        if must_initialize_coefficient_stream:
            self.initialize_coefficient_stream(compute_coefficients)

        if hasattr(self, '_reference') and self._reference is not None:
            self._reference._copy(self)

    def _new(self, compute_coefficients, order_op, *series, **kwds):
        parent = kwds['parent'] if 'parent' in kwds else self.parent()
        new_fps = self.__class__(parent, stream=None, order=unk, aorder=self.aorder,
                                 aorder_changed=True, is_initialized=False)

        new_fps.compute_aorder = lambda: new_fps._approximate_order(compute_coefficients, order_op, *series)
        return new_fps

    def _add_(self, y):
        """
        EXAMPLES: Test Plus 1

        ::

            sage: from sage.combinat.species.series import *
            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: gs0 = L([0])
            sage: gs1 = L([1])
            sage: sum1 = gs0 + gs1
            sage: sum2 = gs1 + gs1
            sage: sum3 = gs1 + gs0
            sage: [gs0.coefficient(i) for i in range(11)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: [gs1.coefficient(i) for i in range(11)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [sum1.coefficient(i) for i in range(11)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [sum2.coefficient(i) for i in range(11)]
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            sage: [sum3.coefficient(i) for i in range(11)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        Test Plus 2

        ::

            sage: gs1 = L([1,2,4,8,0])
            sage: gs2 = L([-1, 0,-1,-9,22,0])
            sage: sum = gs1 + gs2
            sage: sum2 = gs2 + gs1
            sage: [ sum.coefficient(i) for i in range(5) ]
            [0,  2, 3, -1, 22]
            sage: [ sum.coefficient(i) for i in range(5, 11) ]
            [0, 0, 0, 0, 0, 0]
            sage: [ sum2.coefficient(i) for i in range(5) ]
            [0,  2, 3, -1, 22]
            sage: [ sum2.coefficient(i) for i in range(5, 11) ]
            [0, 0, 0, 0, 0, 0]
        """
        return self._new(partial(self._plus_gen, y), min, self, y)

    add = _add_



    def _plus_gen(self, y, ao):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: gs1 = L([1])
            sage: g = gs1._plus_gen(gs1, 0)
            sage: [next(g) for i in range(5)]
            [2, 2, 2, 2, 2]

        ::

            sage: g = gs1._plus_gen(gs1, 2)
            sage: [next(g) for i in range(5)]
            [0, 0, 2, 2, 2]
        """
        base_ring = self.parent().base_ring()
        zero = base_ring(0)
        for n in range(ao):
            yield zero
        n = ao
        while True:
            yield self._stream[n] + y._stream[n]
            n += 1

    def _mul_(self, y):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: gs0 = L(0)
            sage: gs1 = L([1])

        ::

            sage: prod0 = gs0 * gs1
            sage: [prod0.coefficient(i) for i in range(11)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        ::

            sage: prod1 = gs1 * gs0
            sage: [prod1.coefficient(i) for i in range(11)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        ::

            sage: prod2 = gs1 * gs1
            sage: [prod2.coefficient(i) for i in range(11)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

        ::

            sage: gs1 = L([1,2,4,8,0])
            sage: gs2 = L([-1, 0,-1,-9,22,0])

        ::

            sage: prod1 = gs1 * gs2
            sage: [prod1.coefficient(i) for i in range(11)]
            [-1, -2, -5, -19, 0, 0, 16, 176, 0, 0, 0]

        ::

            sage: prod2 = gs2 * gs1
            sage: [prod2.coefficient(i) for i in range(11)]
            [-1, -2, -5, -19, 0, 0, 16, 176, 0, 0, 0]
        """

        return self._new(partial(self._times_gen, y), lambda a,b:a+b, self, y)

    times = _mul_

    def _times_gen(self, y, ao):
        r"""
        Return an iterator for the coefficients of self \* y.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,1,0])
            sage: g = f._times_gen(f,0)
            sage: [next(g) for i in range(5)]
            [1, 2, 1, 0, 0]
        """
        base_ring = self.parent().base_ring()
        zero = base_ring(0)

        for n in range(ao):
            yield zero

        n = ao
        while True:
            low = self.aorder
            high = n - y.aorder
            nth_coefficient = zero

            #Handle the zero series
            if low == inf or high == inf:
                yield zero
                n += 1
                continue

            for k in range(low, high+1):
                cx = self._stream[k]
                if cx == 0:
                    continue
                nth_coefficient += cx * y._stream[n-k]
            yield nth_coefficient
            n += 1

    def __pow__(self, n):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,1,0])  # 1+x
            sage: g = f^3
            sage: g.coefficients(4)
            [1, 3, 3, 1]

        ::

            sage: f^0
            1
        """
        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError("n must be a nonnegative integer")
        return prod([self]*n, self.parent().identity_element())

    def __invert__(self):
        """
        Return 1 over this power series, i.e. invert this power series.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: x = L.gen()

        Geometric series::

            sage: a = ~(1-x); a.compute_coefficients(10); a
            1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + O(x^11)

        (Shifted) Fibonacci numbers::

            sage: b = ~(1-x-x^2); b.compute_coefficients(10); b
            1 + x + 2*x^2 + 3*x^3 + 5*x^4 + 8*x^5
            + 13*x^6 + 21*x^7 + 34*x^8 + 55*x^9 + 89*x^10 + O(x^11)

        Series whose constant coefficient is `0` cannot be inverted::

            sage: ~x
            Traceback (most recent call last):
            ....
            ZeroDivisionError: cannot invert x because constant coefficient is 0
        """
        if self.get_aorder() > 0:
            raise ZeroDivisionError(
                'cannot invert {} because '
                'constant coefficient is 0'.format(self))
        return self._new(self._invert_gen, lambda a: 0, self)

    invert = __invert__

    def _invert_gen(self, ao):
        r"""
        Return an iterator for the coefficients of 1 over this power series.

        TESTS::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1, -1, 0])
            sage: g = f._invert_gen(0)
            sage: [next(g) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        from itertools import count

        assert ao == 0

        ic0 = ~self.coefficient(0)
        yield ic0
        if self.order == 0:
            return

        one = self.parent()(1)
        base = one - ic0 * self
        base.coefficient(0)
        ao_base = base.get_aorder()
        assert ao_base >= 1

        current = one + base
        k = 1
        for n in count(1):
            while ao_base*k < n:
                current = one + base * current
                k += 1
                current.coefficient(n)  # make sure new current is initialized
            ao_base = base.get_aorder()  # update this so that while above is faster
            yield current.coefficient(n) * ic0

    def _div_(self, other):
        """
        Divide this power series by ``other``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: x = L.gen()

        Fibonacci numbers::

            sage: b = x / (1-x-x^2); b.compute_coefficients(10); b
            x + x^2 + 2*x^3 + 3*x^4 + 5*x^5 + 8*x^6
            + 13*x^7 + 21*x^8 + 34*x^9 + 55*x^10 + O(x^11)
        """
        return self * ~other

    div = _div_

    def __call__(self, y):
        """
        Return the composition of this power series and the power series y.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([1])
            sage: t = L([0,0,1])
            sage: u = s(t)
            sage: u.coefficients(11)
            [1, 0, 1, 1, 2, 3, 5, 8, 13, 21, 34]

        Test Compose 2

        ::

            sage: s = L([1])
            sage: t = L([0,0,1,0])
            sage: u = s(t)
            sage: u.aorder
            0
            sage: u.order
            Unknown series order
            sage: u.coefficients(10)
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
            sage: u.aorder
            0
            sage: u.order
            0

        Test Compose 3 s = 1/(1-x), t = x/(1-x) s(t) = (1-x)/(1-2x)

        ::

            sage: s = L([1])
            sage: t = L([0,1])
            sage: u = s(t)
            sage: u.coefficients(14)
            [1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
        """
        return self._new(partial(self._compose_gen, y), lambda a,b:a*b, self, y)

    composition = __call__

    def _compose_gen(self, y, ao):
        """
        Return a iterator for the coefficients of the composition of this
        power series with the power series y.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([1])
            sage: t = L([0,1])
            sage: g = s._compose_gen(t, 0)
            sage: [next(g) for i in range(10)]
            [1, 1, 2, 4, 8, 16, 32, 64, 128, 256]
        """
        assert y.coefficient(0) == 0
        yield self._stream[0]
        z = self.tail().compose(y) * y
        z.coefficient(1)
        n = 1
        while True:
            yield z._stream[n]
            n += 1

    def tail(self):
        """
        Return the power series whose coefficients obtained by subtracting
        the constant term from this series and then dividing by x.

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(range(20))
            sage: g = f.tail()
            sage: g.coefficients(10)
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        return self._new(lambda a0: self.iterator(1), bounded_decrement, self)

    def iterator(self, n=0, initial=None):
        """
        Return an iterator for the coefficients of self starting at n.

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(range(10))
            sage: g = f.iterator(2)
            sage: [next(g) for i in range(5)]
            [2, 3, 4, 5, 6]
            sage: g = f.iterator(2, initial=[0,0])
            sage: [next(g) for i in range(5)]
            [0, 0, 2, 3, 4]
        """
        if initial is not None:
            for x in initial:
                yield x
        while True:
            yield self._stream[n]
            n += 1

    compose = __call__

    def _power_gen(self):
        """
        Return a generator for all the powers self^k starting with k = 1.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,1,0])
            sage: g = f._power_gen()
            sage: next(g).coefficients(5)
            [1, 1, 0, 0, 0]
            sage: next(g).coefficients(5)
            [1, 2, 1, 0, 0]
            sage: next(g).coefficients(5)
            [1, 3, 3, 1, 0]
        """
        z = self
        while True:
            yield z
            z = z*self

    def derivative(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: one = L(1)
            sage: monom = L.gen()
            sage: s = L([1])
            sage: u = s.derivative()
            sage: u.coefficients(10)
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        ::

            sage: s = L()
            sage: s._name = 's'
            sage: s.define(one+monom*s*s)
            sage: u = s.derivative()
            sage: u.coefficients(5) #[1*1, 2*2, 3*5, 4*14, 5*42]
            [1, 4, 15, 56, 210]

        ::

            sage: s = L([1])
            sage: t = L([0,1])
            sage: u = s(t).derivative()
            sage: v = (s.derivative().compose(t))*t.derivative()
            sage: u.coefficients(11)
            [1, 4, 12, 32, 80, 192, 448, 1024, 2304, 5120, 11264]
            sage: v.coefficients(11)
            [1, 4, 12, 32, 80, 192, 448, 1024, 2304, 5120, 11264]

        ::

            sage: s = L(); s._name='s'
            sage: t = L(); t._name='t'
            sage: s.define(monom+t*t*t)
            sage: t.define(monom+s*s)
            sage: u = (s*t).derivative()
            sage: v = s.derivative()*t + s*t.derivative()
            sage: u.coefficients(10)
            [0, 2, 3, 4, 30, 72, 133, 552, 1791, 4260]
            sage: v.coefficients(10)
            [0, 2, 3, 4, 30, 72, 133, 552, 1791, 4260]
            sage: u.coefficients(10) == v.coefficients(10)
            True

        ::

            sage: f = L._new_initial(2, Stream([0,0,4,5,6,0]))
            sage: d = f.derivative()
            sage: d.get_aorder()
            1
            sage: d.coefficients(5)
            [0, 8, 15, 24, 0]
        """
        return self._new(self._diff_gen, bounded_decrement, self)

    def _diff_gen(self, ao):
        """
        Return an iterator for the coefficients of the derivative of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1])
            sage: g = f._diff_gen(0)
            sage: [next(g) for i in range(10)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        n = 1
        while True:
            yield n*self._stream[n]
            n += 1

    ###########
    #Integrals#
    ###########
    def integral(self, integration_constant = 0):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: zero = L(0)
            sage: s = zero
            sage: t = s.integral()
            sage: t.is_zero()
            True

        ::

            sage: s = zero
            sage: t = s.integral(1)
            sage: t.coefficients(6)
            [1, 0, 0, 0, 0, 0]
            sage: t._stream.is_constant()
            True

        ::

            sage: s = L.term(1, 0)
            sage: t = s.integral()
            sage: t.coefficients(6)
            [0, 1, 0, 0, 0, 0]
            sage: t._stream.is_constant()
            True

        ::

            sage: s = L.term(1,0)
            sage: t = s.integral(1)
            sage: t.coefficients(6)
            [1, 1, 0, 0, 0, 0]
            sage: t._stream.is_constant()
            True

        ::

            sage: s = L.term(1, 4)
            sage: t = s.integral()
            sage: t.coefficients(10)
            [0, 0, 0, 0, 0, 1/5, 0, 0, 0, 0]

        ::

            sage: s = L.term(1,4)
            sage: t = s.integral(1)
            sage: t.coefficients(10)
            [1, 0, 0, 0, 0, 1/5, 0, 0, 0, 0]

        TESTS::

            sage: from sage.combinat.species.stream import Stream
            sage: f = L._new_initial(2, Stream([0,0,4,5,6,0]))
            sage: i = f.derivative().integral()
            sage: i.get_aorder()
            2
            sage: i.coefficients(5)
            [0, 0, 4, 5, 6]
            sage: i = f.derivative().integral(1)
            sage: i.get_aorder()
            0
            sage: i.coefficients(5)
            [1, 0, 4, 5, 6]
        """
        if integration_constant == 0:
            return self._new(self._integral_zero_gen, increment, self)
        else:
            L = self.parent()
            return L._new_initial(0, Stream(self._integral_nonzero_gen(integration_constant)))

    def _integral_zero_gen(self, ao):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L.gen()
            sage: g = s._integral_zero_gen(1)
            sage: [next(g) for i in range(5)]
            [0, 0, 1/2, 0, 0]
        """
        for n in range(ao):
            yield self.parent().zero()
        n = ao
        while True:
            #Check to see if the stream is finite
            if self.is_finite(n-1):
                yield self._stream[n-1]
                break
            else:
                yield (Integer(1)/Integer(n))*self._stream[n-1]
                n += 1


    def _integral_nonzero_gen(self, integration_constant):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L._new_initial(2, Stream([0,0,4,5,6,0])).derivative()
            sage: g = f._integral_nonzero_gen(1)
            sage: [next(g) for i in range(5)]
            [1, 0, 4, 5, 6]
        """
        yield integration_constant
        ao = self.aorder
        assert ao != unk

        if ao == inf:
            yield self.parent()._zero_base_ring
        else:
            for _ in range(ao-1):
                yield self.parent()._zero_base_ring

            n = max(1, ao)
            while True:
                self.coefficient(n - 1)

                #Check to see if the stream is finite
                if self.is_finite(n-1):
                    yield self.coefficient(n-1)
                    break
                else:
                    yield (Integer(1)/Integer(n))*self.coefficient(n-1)
                    n += 1

    def is_finite(self, n=None):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([0,0,1,0,0]); a
            O(1)
            sage: a.is_finite()
            False
            sage: c = a[4]
            sage: a.is_finite()
            False
            sage: a.is_finite(4)
            False
            sage: c = a[5]
            sage: a.is_finite()
            True
            sage: a.is_finite(4)
            True
        """
        if self.order is inf:
            return True

        s = self._stream

        if n is None:
            n = len(s)

        if s.is_constant() and all(s[i] == 0 for i in range(n-1, max(n,len(s)))):
            return True

        return False

    def exponential(self):
        """
        TESTS::

            sage: def inv_factorial():
            ....:     q = 1
            ....:     yield 0
            ....:     yield q
            ....:     n = 2
            ....:     while True:
            ....:         q = q / n
            ....:         yield q
            ....:         n += 1
            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(inv_factorial()) #e^(x)-1
            sage: u = f.exponential()
            sage: g = inv_factorial()
            sage: z1 = [1,1,2,5,15,52,203,877,4140,21147,115975]
            sage: l1 = [z*next(g) for z in z1]
            sage: l1 = [1] + l1[1:]
            sage: u.coefficients(11)
            [1, 1, 1, 5/6, 5/8, 13/30, 203/720, 877/5040, 23/224, 1007/17280, 4639/145152]
            sage: l1 == u.coefficients(11)
            True
        """
        base_ring = self.parent().base_ring()
        s = self.parent()()
        s.define( (self.derivative()*s).integral(base_ring(1)) )
        return s

    def __getitem__(self, i):
        """
        Return the ith coefficient of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,2,3,0])
            sage: [f[i] for i in range(5)]
            [1, 2, 3, 0, 0]
        """
        return self.coefficient(i)


    #########################
    #Min and max restriction#
    #########################
    def restricted(self, min=None, max=None):
        """
        Return the power series restricted to the coefficients starting at
        ``min`` and going up to, but not including ``max``.

        If ``min`` is not specified, then it is assumed to be zero. If
        ``max`` is not specified, then it is assumed to be infinity.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1])
            sage: a.restricted().coefficients(10)
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: a.restricted(min=2).coefficients(10)
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: a.restricted(max=5).coefficients(10)
            [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
            sage: a.restricted(min=2, max=6).coefficients(10)
            [0, 0, 1, 1, 1, 1, 0, 0, 0, 0]
        """

        if ((min is None and max is None) or
            (max is None and self.get_aorder() >= min)):
            return self

        if min is None:
            min = 0
        return self._new(partial(self._restricted_gen, min, max),
                         lambda ao: builtins.max(ao, min), self)

    def _restricted_gen(self, mn, mx, ao):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1])
            sage: g = a._restricted_gen(None, None, 2)
            sage: [next(g) for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: g = a._restricted_gen(1, None, 2)
            sage: [next(g) for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: g = a._restricted_gen(3, None, 2)
            sage: [next(g) for i in range(10)]
            [0, 0, 0, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: g = a._restricted_gen(1, 5, 2)
            sage: [next(g) for i in range(6)]
            [0, 0, 1, 1, 1, 0]
        """
        BR = self.parent().base_ring()
        for n in range(max(mn,ao)):
            yield BR(0)

        n = max(mn, ao)
        while True:
            if mx is not None and n >= mx:
                yield BR(0)
                break
            else:
                yield self._stream[n]
            n += 1


    #############
    #Change Ring#
    #############
    def _change_ring_gen(self, R, ao):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L2 = LazyPowerSeriesRing(RR)
            sage: a = L([1])
            sage: b = L2(a)
            sage: b.parent()
            Lazy Power Series Ring over Real Field with 53 bits of precision
            sage: b.coefficients(3)
            [1.00000000000000, 1.00000000000000, 1.00000000000000]
        """
        for n in range(ao):
            yield R(0)

        n = ao
        while True:
            yield R(self._stream[n])
            n += 1

#################################



def uninitialized():
    """
    EXAMPLES::

        sage: from sage.combinat.species.series import uninitialized
        sage: uninitialized()
        Traceback (most recent call last):
        ...
        RuntimeError: we should never be here
    """
    raise RuntimeError("we should never be here")
