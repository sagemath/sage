r"""
Interval Exchange Transformations and Linear Involution

An interval exchage transformation is a map defined on an interval (see
help(iet.IntervalExchangeTransformation) for a more complete help.

EXAMPLES:

Initialization of a simple iet with integer lengths::

    sage: T = iet.IntervalExchangeTransformation(Permutation([3,2,1]), [3,1,2])
    sage: T
    Interval exchange transformation of [0, 6[ with permutation
    1 2 3
    3 2 1

Rotation corresponds to iet with two intervals::

    sage: p = iet.Permutation('a b', 'b a')
    sage: T = iet.IntervalExchangeTransformation(p, [1, (sqrt(5)-1)/2])
    sage: print(T.in_which_interval(0))
    a
    sage: print(T.in_which_interval(T(0)))
    a
    sage: print(T.in_which_interval(T(T(0))))
    b
    sage: print(T.in_which_interval(T(T(T(0)))))
    a

There are two plotting methods for iet::

    sage: p = iet.Permutation('a b c','c b a')
    sage: T = iet.IntervalExchangeTransformation(p, [1, 2, 3])

.. plot the domain and the range of T::

    sage: T.plot_two_intervals()
    Graphics object consisting of 12 graphics primitives

.. plot T as a function::

    sage: T.plot_function()
    Graphics object consisting of 3 graphics primitives
"""
from __future__ import print_function
from __future__ import absolute_import

from copy import copy
from sage.structure.sage_object import SageObject

from .template import side_conversion, interval_conversion

class IntervalExchangeTransformation(SageObject):
    r"""
    Interval exchange transformation

    INPUT:

    - ``permutation`` - a permutation (LabelledPermutationIET)

    - ``lengths`` - the list of lengths

    EXAMPLES:

    Direct initialization::

        sage: p = iet.IET(('a b c','c b a'),{'a':1,'b':1,'c':1})
        sage: p.permutation()
        a b c
        c b a
        sage: p.lengths()
        [1, 1, 1]

    Initialization from a iet.Permutation::

        sage: perm = iet.Permutation('a b c','c b a')
        sage: l = [0.5,1,1.2]
        sage: t = iet.IET(perm,l)
        sage: t.permutation() == perm
        True
        sage: t.lengths() == l
        True

    Initialization from a Permutation::

        sage: p = Permutation([3,2,1])
        sage: iet.IET(p, [1,1,1])
        Interval exchange transformation of [0, 3[ with permutation
        1 2 3
        3 2 1

    If it is not possible to convert lengths to real values an error is raised::

        sage: iet.IntervalExchangeTransformation(('a b','b a'),['e','f'])
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'e' to a float

    The value for the lengths must be positive::

        sage: iet.IET(('a b','b a'),[-1,-1])
        Traceback (most recent call last):
        ...
        ValueError: lengths must be positive
    """
    def __init__(self,permutation=None,lengths=None):
        r"""
        INPUT:

        - ``permutation`` - a permutation (LabelledPermutationIET)

        - ``lengths`` - the list of lengths

        TEST::

            sage: p=iet.IntervalExchangeTransformation(('a','a'),[1])
            sage: p == loads(dumps(p))
            True
        """
        from .labelled import LabelledPermutationIET
        if permutation is None or lengths is None:
            self._permutation = LabelledPermutationIET()
            self._lengths = []
        else:
            self._permutation = permutation
            self._lengths = lengths

    def permutation(self):
        r"""
        Returns the permutation associated to this iet.

        OUTPUT:

        permutation -- the permutation associated to this iet

        EXAMPLES::

            sage: perm = iet.Permutation('a b c','c b a')
            sage: p = iet.IntervalExchangeTransformation(perm,(1,2,1))
            sage: p.permutation() == perm
            True
        """
        return copy(self._permutation)

    def length(self):
        r"""
        Returns the total length of the interval.

        OUTPUT:

        real -- the length of the interval

        EXAMPLES::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t.length()
            2
        """
        return sum(self._lengths)

    def lengths(self):
        r"""
        Returns the list of lengths associated to this iet.

        OUTPUT:

        list -- the list of lengths of subinterval

        EXAMPLES::

            sage: p = iet.IntervalExchangeTransformation(('a b','b a'),[1,3])
            sage: p.lengths()
            [1, 3]
        """
        return copy(self._lengths)

    def normalize(self, total=1):
        r"""
        Returns a interval exchange transformation of normalized lengths.

        The normalization consists in multiplying all lengths by a
        constant in such way that their sum is given by ``total``
        (default is 1).

        INPUT:

        - ``total`` - (default: 1) The total length of the interval

        OUTPUT:

        iet -- the normalized iet

        EXAMPLES::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'), [1,3])
            sage: t.length()
            4
            sage: s = t.normalize(2)
            sage: s.length()
            2
            sage: s.lengths()
            [1/2, 3/2]

        TESTS::

           sage: s = t.normalize('bla')
           Traceback (most recent call last):
           ...
           TypeError: unable to convert 'bla' to a float
           sage: s = t.normalize(-691)
           Traceback (most recent call last):
           ...
           ValueError: the total length must be positive
        """
        try:
            float(total)
        except ValueError:
            raise TypeError("unable to convert {!r} to a float".format(total))

        if total <= 0:
            raise ValueError("the total length must be positive")

        res = copy(self)
        coeff = total / res.length()
        res._multiply_lengths(coeff)
        return res

    def _multiply_lengths(self, x):
        r"""
        Multiplies the lengths of self by x (no verification on x).

        INPUT:

        - ``x`` - a positive number

        TESTS::

            sage: t = iet.IET(("a","a"), [1])
            sage: t.lengths()
            [1]
            sage: t._multiply_lengths(2)
            sage: t.lengths()
            [2]
        """
        self._lengths = [t*x for t in self._lengths]

    def _repr_(self):
        r"""
        A representation string.

        EXAMPLES::

            sage: a = iet.IntervalExchangeTransformation(('a','a'),[1])
            sage: a   # indirect doctest
            Interval exchange transformation of [0, 1[ with permutation
            a
            a
        """
        interval = "[0, %s["%self.length()
        s = "Interval exchange transformation of %s "%interval
        s += "with permutation\n%s"%self._permutation
        return s

    def is_identity(self):
        r"""
        Returns True if self is the identity.

        OUTPUT:

        boolean -- the answer

        EXAMPLES::

            sage: p = iet.Permutation("a b","b a")
            sage: q = iet.Permutation("c d","d c")
            sage: s = iet.IET(p, [1,5])
            sage: t = iet.IET(q, [5,1])
            sage: (s*t).is_identity()
            True
            sage: (t*s).is_identity()
            True
        """
        return self._permutation.is_identity()

    def inverse(self):
        r"""
        Returns the inverse iet.

        OUTPUT:

        iet -- the inverse interval exchange transformation

        EXAMPLES::

            sage: p = iet.Permutation("a b","b a")
            sage: s = iet.IET(p, [1,sqrt(2)-1])
            sage: t = s.inverse()
            sage: t.permutation()
            b a
            a b
            sage: t.lengths()
            [1, sqrt(2) - 1]
            sage: t*s
            Interval exchange transformation of [0, sqrt(2)[ with permutation
            aa bb
            aa bb

        We can verify with the method .is_identity()::

            sage: p = iet.Permutation("a b c d","d a c b")
            sage: s = iet.IET(p, [1, sqrt(2), sqrt(3), sqrt(5)])
            sage: (s * s.inverse()).is_identity()
            True
            sage: (s.inverse() * s).is_identity()
            True
        """
        res = copy(self)
        res._permutation._inversed()
        return res

    def __mul__(self, other):
        r"""
        Composition of iet.

        The domain (i.e. the length) of the two iets must be the same). The
        alphabet choosen depends on the permutation.

        TESTS:

        ::

            sage: p = iet.Permutation("a b", "a b")
            sage: t = iet.IET(p, [1,1])
            sage: r = t*t
            sage: r.permutation()
            aa bb
            aa bb
            sage: r.lengths()
            [1, 1]

        ::

            sage: p = iet.Permutation("a b","b a")
            sage: t = iet.IET(p, [1,1])
            sage: r = t*t
            sage: r.permutation()
            ab ba
            ab ba
            sage: r.lengths()
            [1, 1]

        ::

            sage: p = iet.Permutation("1 2 3 4 5","5 4 3 2 1")
            sage: q = iet.Permutation("a b","b a")
            sage: s = iet.IET(p, [1]*5)
            sage: t = iet.IET(q, [1/2, 9/2])
            sage: r = s*t
            sage: r.permutation()
            a5 b1 b2 b3 b4 b5
            b5 a5 b4 b3 b2 b1
            sage: r.lengths()
            [1/2, 1, 1, 1, 1, 1/2]
            sage: r = t*s
            sage: r.permutation()
            1b 2b 3b 4b 5a 5b
            5b 4b 3b 2b 1b 5a
            sage: r.lengths()
            [1, 1, 1, 1, 1/2, 1/2]
            sage: t = iet.IET(q, [3/2, 7/2])
            sage: r = s*t
            sage: r.permutation()
            a4 a5 b1 b2 b3 b4
            a5 b4 a4 b3 b2 b1
            sage: r.lengths()
            [1/2, 1, 1, 1, 1, 1/2]
            sage: t = iet.IET(q, [5/2,5/2])
            sage: r = s*t
            sage: r.permutation()
            a3 a4 a5 b1 b2 b3
            a5 a4 b3 a3 b2 b1
            sage: r = t*s
            sage: r.permutation()
            1b 2b 3a 3b 4a 5a
            3b 2b 1b 5a 4a 3a

        ::

            sage: p = iet.Permutation("a b","b a")
            sage: s = iet.IET(p, [4,2])
            sage: q = iet.Permutation("c d","d c")
            sage: t = iet.IET(q, [3, 3])
            sage: r1 = t * s
            sage: r1.permutation()
            ac ad bc
            ad bc ac
            sage: r1.lengths()
            [1, 3, 2]
            sage: r2 = s * t
            sage: r2.permutation()
            ca cb da
            cb da ca
            sage: r2.lengths()
            [1, 2, 3]

        ::

            sage: r * s
            Traceback (most recent call last):
            ...
            ValueError: self and other are not IET of the same length
        """
        if not(isinstance(other, IntervalExchangeTransformation) and
               self.length() == other.length()):
            raise ValueError("self and other are not IET of the same length")

        from .labelled import LabelledPermutationIET

        other_sg = other.range_singularities()[1:]
        self_sg = self.domain_singularities()[1:]

        n_other = len(other._permutation)
        n_self = len(self._permutation)

        interval_other = other._permutation._intervals[1]
        interval_self = self._permutation._intervals[0]

        d_other = dict([(i,[]) for i in interval_other])
        d_self = dict([(i,[]) for i in interval_self])

        i_other = 0
        i_self = 0

        x = 0
        l_lengths = []
        while i_other < n_other and i_self < n_self:
            j_other = interval_other[i_other]
            j_self = interval_self[i_self]

            d_other[j_other].append(j_self)
            d_self[j_self].append(j_other)

            if other_sg[i_other] < self_sg[i_self]:
                l = other_sg[i_other] - x
                x = other_sg[i_other]
                i_other += 1
            elif other_sg[i_other] > self_sg[i_self]:
                l = self_sg[i_self] - x
                x = self_sg[i_self]
                i_self += 1
            else:
                l = self_sg[i_self] - x
                x = self_sg[i_self]
                i_other += 1
                i_self += 1

            l_lengths.append(((j_other,j_self),l))

        alphabet_other = other._permutation.alphabet()
        alphabet_self = self._permutation.alphabet()

        d_lengths = dict(l_lengths)

        l_lengths = []
        top_interval = []
        for i in other._permutation._intervals[0]:
            for j in d_other[i]:
                a = alphabet_other.unrank(i)
                b = alphabet_self.unrank(j)
                top_interval.append(str(a)+str(b))
                l_lengths.append(d_lengths[(i,j)])

        bottom_interval = []
        for i in self._permutation._intervals[1]:
            for j in d_self[i]:
                a = alphabet_other.unrank(j)
                b = alphabet_self.unrank(i)
                bottom_interval.append(str(a)+str(b))

        p = LabelledPermutationIET((top_interval,bottom_interval))
        return IntervalExchangeTransformation(p,l_lengths)

    def __eq__(self, other):
        r"""
        Tests equality

        TESTS::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t == t
            True
        """
        return (
            type(self) is type(other) and
            self._permutation == other._permutation and
            self._lengths == other._lengths)

    def __ne__(self, other):
        r"""
        Tests difference

        TESTS::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t != t
            False
        """
        return (
            type(self) is not type(other) or
            self._permutation != other._permutation or
            self._lengths != other._lengths)

    def in_which_interval(self, x, interval=0):
        r"""
        Returns the letter for which x is in this interval.

        INPUT:

        - ``x`` - a positive number

        - ``interval`` - (default: 'top') 'top' or 'bottom'


        OUTPUT:

        label -- a label corresponding to an interval

        TEST:

        ::

            sage: t = iet.IntervalExchangeTransformation(('a b c','c b a'),[1,1,1])
            sage: t.in_which_interval(0)
            'a'
            sage: t.in_which_interval(0.3)
            'a'
            sage: t.in_which_interval(1)
            'b'
            sage: t.in_which_interval(1.9)
            'b'
            sage: t.in_which_interval(2)
            'c'
            sage: t.in_which_interval(2.1)
            'c'
            sage: t.in_which_interval(3)
            Traceback (most recent call last):
            ...
            ValueError: your value does not lie in [0;l[

        .. and for the bottom interval::

            sage: t.in_which_interval(0,'bottom')
            'c'
            sage: t.in_which_interval(1.2,'bottom')
            'b'
            sage: t.in_which_interval(2.9,'bottom')
            'a'

        TESTS::

            sage: t.in_which_interval(-2.9,'bottom')
            Traceback (most recent call last):
            ...
            ValueError: your value does not lie in [0;l[
        """
        interval = interval_conversion(interval)

        if x < 0 or x >= self.length():
            raise ValueError("your value does not lie in [0;l[")

        i = 0

        while x >= 0:
            x -= self._lengths[self._permutation._intervals[interval][i]]
            i += 1

        i -= 1
        x += self._lengths[self._permutation._intervals[interval][i]]

        j = self._permutation._intervals[interval][i]
        return self._permutation._alphabet.unrank(j)

    def singularities(self):
        r"""
        The list of singularities of `T` and `T^{-1}`.

        OUTPUT:

        list -- two lists of positive numbers which corresponds to extremities
            of subintervals

        EXAMPLE::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1/2,3/2])
            sage: t.singularities()
            [[0, 1/2, 2], [0, 3/2, 2]]
        """
        return [self.domain_singularities(), self.range_singularities()]

    def domain_singularities(self):
        r"""
        Returns the list of singularities of T

        OUTPUT:

        list -- positive reals that corresponds to singularities in the top
            interval

        EXAMPLES::

            sage: t = iet.IET(("a b","b a"), [1, sqrt(2)])
            sage: t.domain_singularities()
            [0, 1, sqrt(2) + 1]
        """
        l = [0]
        for j in self._permutation._intervals[0]:
            l.append(l[-1] + self._lengths[j])
        return l

    def range_singularities(self):
        r"""
        Returns the list of singularities of `T^{-1}`

        OUTPUT:

        list -- real numbers that are singular for `T^{-1}`


        EXAMPLES::

            sage: t = iet.IET(("a b","b a"), [1, sqrt(2)])
            sage: t.range_singularities()
            [0, sqrt(2), sqrt(2) + 1]
        """
        l = [0]
        for j in self._permutation._intervals[1]:
            l.append(l[-1] + self._lengths[j])
        return l

    def __call__(self, value):
        r"""
        Return the image of value by this transformation

        EXAMPLES::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1/2,3/2])
            sage: t(0)
            3/2
            sage: t(1/2)
            0
            sage: t(1)
            1/2
            sage: t(3/2)
            1

        TESTS::

            sage: t(-3/2)
            Traceback (most recent call last):
            ...
            ValueError: value must positive and smaller than length
        """
        if not(value >= 0 and value < self.length()):
            raise ValueError("value must positive and smaller than length")

        dom_sg = self.domain_singularities()
        im_sg = self.range_singularities()

        a = self.in_which_interval(value)

        i0 = self._permutation[0].index(a)
        i1 = self._permutation[1].index(a)

        return value - dom_sg[i0] + im_sg[i1]

    def rauzy_move(self, side='right', iterations=1):
        r"""
        Performs a Rauzy move.

        INPUT:

        - ``side`` - 'left' (or 'l' or 0) or 'right' (or 'r' or 1)

        - ``iterations`` - integer (default :1) the number of iteration of Rauzy
           moves to perform

        OUTPUT:

        iet -- the Rauzy move of self

        EXAMPLES::

            sage: phi = QQbar((sqrt(5)-1)/2)
            sage: t1 = iet.IntervalExchangeTransformation(('a b','b a'),[1,phi])
            sage: t2 = t1.rauzy_move().normalize(t1.length())
            sage: l2 = t2.lengths()
            sage: l1 = t1.lengths()
            sage: l2[0] == l1[1] and l2[1] == l1[0]
            True
        """
        side = side_conversion(side)

        res = copy(self)
        for i in range(iterations):
            res = res._rauzy_move(side)
        return res

    def _rauzy_move(self,side=-1):
        r"""
        Performs a Rauzy move

        INPUT:

        - ``side`` - must be 0 or -1 (no verification)

        TEST::

            sage: t = iet.IntervalExchangeTransformation(('a b c','c b a'),[1,1,3])
            sage: t
            Interval exchange transformation of [0, 5[ with permutation
            a b c
            c b a
            sage: t1 = t.rauzy_move()   #indirect doctest
            sage: t1
            Interval exchange transformation of [0, 4[ with permutation
            a b c
            c a b
            sage: t2 = t1.rauzy_move()   #indirect doctest
            sage: t2
            Interval exchange transformation of [0, 3[ with permutation
            a b c
            c b a
            sage: t2.rauzy_move()   #indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: top and bottom extrem intervals have equal lengths
        """
        top = self._permutation._intervals[0][side]
        bottom = self._permutation._intervals[1][side]

        length_top = self._lengths[top]
        length_bottom = self._lengths[bottom]

        if length_top > length_bottom:
            winner = 0
            winner_interval = top
            loser_interval = bottom
        elif length_top < length_bottom:
            winner = 1
            winner_interval = bottom
            loser_interval = top
        else:
            raise ValueError("top and bottom extrem intervals have equal lengths")

        res = IntervalExchangeTransformation(([],[]),{})
        res._permutation = self._permutation.rauzy_move(winner=winner,side=side)
        res._lengths = self._lengths[:]
        res._lengths[winner_interval] -= res._lengths[loser_interval]

        return res

    def __copy__(self):
        r"""
        Returns a copy of this interval exchange transformation.

        EXAMPLES::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: s = copy(t)
            sage: s == t
            True
            sage: s is t
            False
        """
        res = self.__class__()
        res._permutation = copy(self._permutation)
        res._lengths = copy(self._lengths)
        return res

    def plot_function(self,**d):
        r"""
        Return a plot of the interval exchange transformation as a
        function.

        INPUT:

        - Any option that is accepted by line2d

        OUTPUT:

        2d plot -- a plot of the iet as a function

        EXAMPLES::

            sage: t = iet.IntervalExchangeTransformation(('a b c d','d a c b'),[1,1,1,1])
            sage: t.plot_function(rgbcolor=(0,1,0))
            Graphics object consisting of 4 graphics primitives
        """
        from sage.plot.all import Graphics
        from sage.plot.plot import line2d

        G = Graphics()
        l = self.singularities()
        t = self._permutation._twin

        for i in range(len(self._permutation)):
            j = t[0][i]
            G += line2d([(l[0][i],l[1][j]),(l[0][i+1],l[1][j+1])],**d)

        return G

    def plot_two_intervals(self,
                           position=(0,0),
                           vertical_alignment='center',
                           horizontal_alignment='left',
                           interval_height=0.1,
                           labels_height=0.05,
                           fontsize=14,
                           labels=True,
                           colors=None):
        r"""
        Returns a picture of the interval exchange transformation.

        INPUT:

        - ``position`` - a 2-uple of the position

        - ``horizontal_alignment`` - left (default), center or right

        - ``labels`` - boolean (default: True)

        - ``fontsize`` - the size of the label


        OUTPUT:

        2d plot -- a plot of the two intervals (domain and range)

        EXAMPLES::

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t.plot_two_intervals()
            Graphics object consisting of 8 graphics primitives
        """
        from sage.plot.all import Graphics
        from sage.plot.plot import line2d
        from sage.plot.plot import text
        from sage.plot.colors import rainbow

        G = Graphics()

        lengths = [float(_) for _ in self._lengths]
        total_length = sum(lengths)

        if colors is None:
            colors = rainbow(len(self._permutation), 'rgbtuple')

        if horizontal_alignment == 'left':
            s = position[0]
        elif horizontal_alignment == 'center':
            s = position[0] - total_length / 2
        elif horizontal_alignment == 'right':
            s = position[0] - total_length
        else:
            raise ValueError("horizontal_alignement must be left, center or right")

        top_height = position[1] + interval_height
        for i in self._permutation._intervals[0]:
            G += line2d([(s,top_height), (s+lengths[i],top_height)],
                        rgbcolor=colors[i])
            if labels:
                G += text(str(self._permutation._alphabet.unrank(i)),
                          (s+float(lengths[i])/2, top_height+labels_height),
                          horizontal_alignment='center',
                          rgbcolor=colors[i],
                          fontsize=fontsize)

            s += lengths[i]

        if horizontal_alignment == 'left':
            s = position[0]
        elif horizontal_alignment == 'center':
            s = position[0] - total_length / 2
        elif horizontal_alignment == 'right':
            s = position[0] - total_length
        else:
            raise ValueError("horizontal_alignement must be left, center or right")

        bottom_height = position[1] - interval_height
        for i in self._permutation._intervals[1]:
            G += line2d([(s,bottom_height), (s+lengths[i],bottom_height)],
                        rgbcolor=colors[i])
            if labels:
                G += text(str(self._permutation._alphabet.unrank(i)),
                          (s+float(lengths[i])/2, bottom_height-labels_height),
                          horizontal_alignment='center',
                          rgbcolor=colors[i],
                          fontsize=fontsize)
            s += lengths[i]

        return G

    plot = plot_two_intervals

    def show(self):
        r"""
        Shows a picture of the interval exchange transformation

        EXAMPLES::

            sage: phi = QQbar((sqrt(5)-1)/2)
            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,phi])
            sage: t.show()
        """
        self.plot_two_intervals().show(axes=False)

#TODO
# class LinearInvolution(SageObject):
#     r"""_
#     Linear involutions
#     """
#     pass
