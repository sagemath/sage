"""
Series Order

This file provides some utility classes which are useful when
working with unknown, known, and infinite series orders for
univariate power series.

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevant section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatsu30.html.
"""
from sage.rings.integer import Integer

class SeriesOrderElement:
    pass

class InfiniteSeriesOrder(SeriesOrderElement):
    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: o = InfiniteSeriesOrder(); o
            Infinite series order
        """
        return "Infinite series order"

    def __add__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: o = InfiniteSeriesOrder()
            sage: o + 2
            Infinite series order
            sage: o + o
            Infinite series order

        ::

            sage: u = UnknownSeriesOrder()
            sage: o + u
            Unknown series order

        TESTS::

            sage: o + -1
            Traceback (most recent call last):
            ...
            ValueError: x must be a positive integer
        """
        if isinstance(x, (int, Integer)):
            if x < 0:
                raise ValueError("x must be a positive integer")
            return self

        if isinstance(x, InfiniteSeriesOrder):
            return self

        if isinstance(x, UnknownSeriesOrder):
            return x

        raise TypeError("x must be a positive integer or a SeriesOrderElement")

    __radd__ = __add__


    def __mul__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: o = InfiniteSeriesOrder()
            sage: o * 2
            Infinite series order
            sage: o * o
            Infinite series order

        ::

            sage: u = UnknownSeriesOrder()
            sage: o * u
            Unknown series order

        TESTS::

            sage: o * -1
            Traceback (most recent call last):
            ...
            ValueError: x must be a positive integer
        """
        if isinstance(x, (int, Integer)):
            if x < 0:
                raise ValueError("x must be a positive integer")
            elif x == 0:
                return x
            return self

        if isinstance(x, InfiniteSeriesOrder):
            return self

        if isinstance(x, UnknownSeriesOrder):
            return x

        raise TypeError("x must be a positive integer or a SeriesOrderElement")

    __rmul__ = __mul__

    def __lt__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: o = InfiniteSeriesOrder()
            sage: o < 2
            False
            sage: o < o
            False

        ::

            sage: u = UnknownSeriesOrder()
            sage: o < u
            False
            sage: 2 < o   # TODO: Not Implemented
            True
        """
        if isinstance(x, (int, Integer)):
            if x < 0:
                raise ValueError("x must be a positive integer")
            return False

        if isinstance(x, InfiniteSeriesOrder):
            return False

        if isinstance(x, UnknownSeriesOrder):
            return False


        raise TypeError("x must be a positive integer or a SeriesOrderElement")


    def __gt__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: o = InfiniteSeriesOrder()
            sage: o > 2
            True
        """
        return True

class UnknownSeriesOrder(SeriesOrderElement):
    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: u = UnknownSeriesOrder(); u
            Unknown series order
        """
        return "Unknown series order"

    def __add__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: u = UnknownSeriesOrder()
            sage: u + 2
            Unknown series order
            sage: u + u
            Unknown series order
        """

        if isinstance(x, (int, Integer)):
            if x < 0:
                raise ValueError("x must be a positive integer")
            return self

        if isinstance(x, SeriesOrderElement):
            return self

        raise TypeError("x must be a positive integer or a SeriesOrderElement")

    __radd__ = __add__


    def __mul__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: u = UnknownSeriesOrder()
            sage: u * 2
            Unknown series order
            sage: u * u
            Unknown series order
        """
        if isinstance(x, (int, Integer)):
            if x < 0:
                raise ValueError("x must be a positive integer")
            return self

        if isinstance(x, SeriesOrderElement):
            return self

        raise TypeError("x must be a positive integer or a SeriesOrderElement")

    __rmul__ = __mul__

    def __lt__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: u = UnknownSeriesOrder()
            sage: u < 2
            True
            sage: o = InfiniteSeriesOrder()
            sage: u < o
            True
        """
        if isinstance(x, (int, Integer)):
            if x < 0:
                raise ValueError("x must be a positive integer")
            return True

        if isinstance(x, SeriesOrderElement):
            return True

        raise TypeError("x must be a positive integer or a SeriesOrderElement")

    def __gt__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.species.series_order import *
            sage: u = UnknownSeriesOrder()
            sage: u > 2
            False
        """
        return False

def bounded_decrement(x):
    """
    EXAMPLES::

        sage: from sage.combinat.species.series_order import *
        sage: u = UnknownSeriesOrder()
        sage: bounded_decrement(u)
        Unknown series order
        sage: bounded_decrement(4)
        3
        sage: bounded_decrement(0)
        0
    """
    if isinstance(x, SeriesOrderElement):
        return x

    if isinstance(x, (int, Integer)):
        if x < 0:
            raise ValueError("x must be a positive integer")
        return max(0, x - 1)

    raise TypeError("x must be a positive integer or a SeriesOrderElement")


def increment(x):
    """
    EXAMPLES::

        sage: from sage.combinat.species.series_order import *
        sage: u = UnknownSeriesOrder()
        sage: increment(u)
        Unknown series order
        sage: increment(2)
        3
    """
    if isinstance(x, SeriesOrderElement):
        return x + 1

    if isinstance(x, (int, Integer)):
        if x < 0:
            raise ValueError("x must be a positive integer")
        return x+1

    raise TypeError("x must be a positive integer or a SeriesOrderElement")

inf = InfiniteSeriesOrder()
unk = UnknownSeriesOrder()
