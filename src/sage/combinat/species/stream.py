"""
Streams or Infinite Arrays

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevant section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatse12.html.
"""
import types
from sage.structure.sage_object import SageObject


def _integers_from(n):
    """
    Returns a generator for the integers starting at n.

    EXAMPLES::

        sage: from sage.combinat.species.stream import _integers_from
        sage: g = _integers_from(5)
        sage: [next(g) for i in range(5)]
        [5, 6, 7, 8, 9]
    """
    while True:
        yield n
        n += 1

def _apply_function(func, list):
    """
    Returns an iterator for func(i) for i in list.

    EXAMPLES::

        sage: from sage.combinat.species.stream import _apply_function
        sage: def square(l):
        ....:     l.append(l[-1]^2)
        ....:     return l[-1]
        ...
        sage: l = [2]
        sage: g = _apply_function(square, l)
        sage: [next(g) for i in range(5)]
        [4, 16, 256, 65536, 4294967296]
    """
    while True:
        try:
            yield func(list)
        except Exception:
            break

def Stream(x=None, const=None):
    """
    Returns a stream.

    EXAMPLES: We can create a constant stream by just passing a

    ::

        sage: from sage.combinat.species.stream import Stream
        sage: s = Stream(const=0)
        sage: [s[i] for i in range(10)]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    """
    if const is not None:
        return Stream_class(const=const)
    elif hasattr(x, '__iter__'):
        return Stream_class(iter(x))
    elif isinstance(x, (types.FunctionType, types.LambdaType)):
        return Stream_class(func=x)

    return Stream_class(iter([x,0]))

class Stream_class(SageObject):
    """
    EXAMPLES::

        sage: from sage.combinat.species.stream import Stream
        sage: from builtins import zip
        sage: s = Stream(const=0)
        sage: len(s)
        1
        sage: [x for (x,i) in zip(s, range(4))]
        [0, 0, 0, 0]
        sage: len(s)
        1

    ::

        sage: s = Stream(const=4)
        sage: g = iter(s)
        sage: l1 = [x for (x,i) in zip(g, range(10))]
        sage: l = [4 for k in range(10)]
        sage: l == l1
        True

    ::

        sage: h = lambda l: 1 if len(l) < 2 else l[-1] + l[-2]
        sage: fib = Stream(h)
        sage: [x for (x,i) in zip(fib, range(11))]
        [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]

    ::

        sage: r = [4, 3, 5, 2, 6, 1, 1, 1, 1, 1]
        sage: l = [4, 3, 5, 2, 6, 1]
        sage: s = Stream(l)
        sage: s[3] = -1
        sage: [x for (x,i) in zip(s, r)]
        [4, 3, 5, -1, 6, 1, 1, 1, 1, 1]
        sage: s[5] = -2
        sage: [x for (x,i) in zip(s, r)]
        [4, 3, 5, -1, 6, -2, 1, 1, 1, 1]
        sage: s[6] = -3
        sage: [x for (x,i) in zip(s, r)]
        [4, 3, 5, -1, 6, -2, -3, 1, 1, 1]
        sage: s[8] = -4
        sage: [x for (x,i) in zip(s, r)]
        [4, 3, 5, -1, 6, -2, -3, 1, -4, 1]
        sage: a = Stream(const=0)
        sage: a[2] = 3
        sage: [x for (x,i) in zip(a, range(4))]
        [0, 0, 3, 0]
    """

    def __init__(self, gen=None, const=None, func=None):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream_class, Stream
            sage: s = Stream_class(const=4)
            sage: loads(dumps(s))
            <sage.combinat.species.stream.Stream_class object at ...>

        ::

            sage: sorted(s.__dict__.items())
            [('_constant', 4),
             ('_gen', None),
             ('_last_index', 0),
             ('_list', [4]),
             ('end_reached', True)]

        ::

            sage: s = Stream(ZZ)
            sage: sorted(s.__dict__.items())
            [('_constant', None),
             ('_gen', <generator object at 0x...>),
             ('_last_index', -1),
             ('_list', []),
             ('end_reached', False)]
        """
        #We define self._list up here so that
        #_apply_function can make use of it if
        #it needs to.
        self._list = []


        if func is not None:
            if gen is not None:
                raise ValueError("you cannot specify both a function and a generator")
            gen = _apply_function(func, self._list)

        #Constant stream
        if const is not None:
            self._list = [const]
            self._last_index = 0      # last_index == len(self._list) - 1
            self._gen = None
            self._constant = const
            self.end_reached = True
        else:
            self._last_index = -1     # last_index == len(self._list) - 1
            self._gen = gen
            self._constant = const
            self.end_reached = False

    def __setitem__(self, i, t):
        """
        Set the i-th entry of self to t.

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream

        ::

            sage: s = Stream(const=0)
            sage: s[5]
            0
            sage: s.data()
            [0]
            sage: s[5] = 5
            sage: s[5]
            5
            sage: s.data()
            [0, 0, 0, 0, 0, 5]

        ::

            sage: s = Stream(ZZ)
            sage: s[10]
            -5
            sage: s.data()
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5]
            sage: s[10] = 10
            sage: s.data()
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 10]
        """
        # Compute all of the coefficients up to (and including) the ith one
        self[i]

        if i < len(self._list):
            #If we are here, we can just change the entry in self._list
            self._list[i] = t
        else:
            #If we are here, then the stream has become constant.  We just
            #extend self._list with self._constant and then change the
            #last entry.
            self._list += [ self._constant ] * (i+1 - len(self._list))
            self._last_index = i
            self._list[i] = t

    def set_gen(self, gen):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: from builtins import zip
            sage: fib = Stream()
            sage: def g():
            ....:        yield 1
            ....:        yield 1
            ....:        n = 0
            ....:        while True:
            ....:            yield fib[n] + fib[n+1]
            ....:            n += 1

        ::

            sage: fib.set_gen(g())
            sage: [x for (x,i) in zip(fib, range(11))]
            [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]

        ::

            sage: l = [4,3,5,2,6,1]
            sage: s = Stream(l)
            sage: s[3]
            2
            sage: len(s)
            4
            sage: g = iter(l)
            sage: s.set_gen(g)
            sage: s[5]
            3
            sage: len(s)
            6
        """
        self._gen = gen
        self.end_reached = False

    def map(self, f):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: s = Stream(ZZ)
            sage: square = lambda x: x^2
            sage: ss = s.map(square)
            sage: [ss[i] for i in range(10)]
            [0, 1, 1, 4, 4, 9, 9, 16, 16, 25]

        TESTS::

            sage: from builtins import zip
            sage: f = lambda l: 0 if len(l) == 0 else l[-1] + 1
            sage: o = Stream(f)
            sage: [x for (x,i) in zip(o, range(10))]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: double = lambda z: 2*z
            sage: t = o.map(double)
            sage: [x for (x,i) in zip(t, range(10))]
            [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

        ::

            sage: double = lambda z: 2*z
            sage: o = Stream([0,1,2,3])
            sage: [x for (x,i) in zip(o, range(6))]
            [0, 1, 2, 3, 3, 3]
            sage: t = o.map(double)
            sage: [x for (x,i) in zip(t, range(6))]
            [0, 2, 4, 6, 6, 6]
        """
        return Stream((f(x) for x in self))

    def __getitem__(self, i):
        """
        Returns the ith entry of self.

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: s = Stream(ZZ)
            sage: [s[i] for i in range(10)]
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]
            sage: s[1]
            1

        ::

            sage: s = Stream([1,2,3])
            sage: [s[i] for i in range(10)]
            [1, 2, 3, 3, 3, 3, 3, 3, 3, 3]

        ::

            sage: s = Stream(QQ)
            sage: s[10]
            -3
        """
        if i <= self._last_index:
            return self._list[i]
        elif self.end_reached:
            if self._constant is not False:
                return self._constant
            else:
                raise IndexError("out of position")
        else:
            while self._last_index < i:
                try:
                    self._list.append(next(self._gen))
                    self._last_index += 1
                except StopIteration:
                    self.end_reached = True
                    self._constant = self._list[-1]
                    return self[i]

            return self._list[i]


    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: s = Stream([1,2,3])
            sage: g = iter(s)
            sage: [next(g) for i in range(5)]
            [1, 2, 3, 3, 3]
        """
        i = 0
        while True:
            try:
                yield self[i]
            except IndexError:
                break
            i += 1

    def __len__(self):
        """
        Returns the number of coefficients computed so far.

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: l = [4,3,5,7,4,1,9,7]
            sage: s = Stream(l)
            sage: s[3]
            7
            sage: len(s)
            4
            sage: s[3]
            7
            sage: len(s)
            4
            sage: s[1]
            3
            sage: len(s)
            4
            sage: s[4]
            4
            sage: len(s)
            5

        TESTS::

            sage: l = ['Hello', ' ', 'World', '!']
            sage: s = Stream(l)
            sage: len(s)
            0
            sage: s[2]
            'World'
            sage: len(s)
            3
            sage: u = ""
            sage: for i in range(len(s)): u += s[i]
            sage: u
            'Hello World'
            sage: v = ""
            sage: for i in range(10): v += s[i]
            sage: v
            'Hello World!!!!!!!'
            sage: len(s)
            4
        """
        return len(self._list)

    number_computed = __len__

    def data(self):
        """
        Returns a list of all the coefficients computed so far.

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream, _integers_from
            sage: s = Stream(_integers_from(3))
            sage: s.data()
            []
            sage: s[5]
            8
            sage: s.data()
            [3, 4, 5, 6, 7, 8]
        """
        return self._list

    def is_constant(self):
        """
        Returns True if and only if

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: s = Stream([1,2,3])
            sage: s.is_constant()
            False
            sage: s[3]
            3
            sage: s.data()
            [1, 2, 3]
            sage: s.is_constant()
            True

        TESTS::

            sage: l = [2,3,5,7,11,0]
            sage: s = Stream(l)
            sage: s.is_constant()
            False
            sage: s[3]
            7
            sage: s.is_constant()
            False
            sage: s[5]
            0
            sage: s.is_constant()
            False
            sage: s[6]
            0
            sage: s.is_constant()
            True

        ::

            sage: s = Stream(const='I am constant.')
            sage: s.is_constant()
            True
        """
        return self.end_reached


    def stretch(self, k):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: s = Stream(range(1, 10))
            sage: s2 = s.stretch(2)
            sage: [s2[i] for i in range(10)]
            [1, 0, 2, 0, 3, 0, 4, 0, 5, 0]
        """
        return Stream(self._stretch_gen(k))

    def _stretch_gen(self, k):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: s = Stream(range(1, 10))
            sage: g = s._stretch_gen(2)
            sage: [next(g) for i in range(10)]
            [1, 0, 2, 0, 3, 0, 4, 0, 5, 0]
        """
        yield self[0]
        for i in _integers_from(1):
            for j in range(k-1):
                yield 0
            yield self[i]
