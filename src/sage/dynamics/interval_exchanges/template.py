r"""
Permutations template

This file define high level operations on permutations (alphabet,
the different rauzy induction, ...) shared by reduced and labeled
permutations.

AUTHORS:

- Vincent Delecroix (2008-12-20): initial version

.. TODO::

    - construct as options different string representations for a permutation

      - the two intervals: str
      - the two intervals on one line: str_one_line
      - the separatrix diagram: str_separatrix_diagram
      - twin[0] and twin[1] for reduced permutation
      - nothing (useful for Rauzy diagram)
"""
#*****************************************************************************
#       Copyright (C) 2008 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.structure.sage_object import SageObject

from copy import copy

from sage.rings.integer import Integer
from sage.combinat.words.alphabet import Alphabet
from sage.graphs.graph import DiGraph
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.nested_class import NestedClassMetaclass


def interval_conversion(interval=None):
    r"""
    Converts the argument in 0 or 1.

    INPUT:

    - ``winner`` - 'top' (or 't' or 0) or bottom (or 'b' or 1)

    OUTPUT:

    integer -- 0 or 1

    TESTS:

    ::

        sage: from sage.dynamics.interval_exchanges.template import interval_conversion
        sage: interval_conversion('top')
        0
        sage: interval_conversion('t')
        0
        sage: interval_conversion(0)
        0
        sage: interval_conversion('bottom')
        1
        sage: interval_conversion('b')
        1
        sage: interval_conversion(1)
        1

    .. Non admissible strings raise a ValueError::

        sage: interval_conversion('')
        Traceback (most recent call last):
        ...
        ValueError: the interval can not be the empty string
        sage: interval_conversion('right')
        Traceback (most recent call last):
        ...
        ValueError: 'right' can not be converted to interval
        sage: interval_conversion('top_right')
        Traceback (most recent call last):
        ...
        ValueError: 'top_right' can not be converted to interval
    """
    if isinstance(interval, (int, Integer)):
        if interval != 0 and interval != 1:
            raise ValueError("interval must be 0 or 1")
        return interval

    if isinstance(interval,str):
        if interval == '':
            raise ValueError("the interval can not be the empty string")
        if 'top'.startswith(interval): return 0
        if 'bottom'.startswith(interval): return 1
        raise ValueError("'%s' can not be converted to interval" % (interval))

    raise TypeError("'%s' is not an admissible type" % (str(interval)))

def side_conversion(side=None):
    r"""
    Converts the argument in 0 or -1.

    INPUT:

    - ``side`` - either 'left' (or 'l' or 0) or 'right' (or 'r' or -1)

    OUTPUT:

    integer -- 0 or -1

    TESTS:

    ::

        sage: from sage.dynamics.interval_exchanges.template import side_conversion
        sage: side_conversion('left')
        0
        sage: side_conversion('l')
        0
        sage: side_conversion(0)
        0
        sage: side_conversion('right')
        -1
        sage: side_conversion('r')
        -1
        sage: side_conversion(1)
        -1
        sage: side_conversion(-1)
        -1

    .. Non admissible strings raise a ValueError::

        sage: side_conversion('')
        Traceback (most recent call last):
        ...
        ValueError: no empty string for side
        sage: side_conversion('top')
        Traceback (most recent call last):
        ...
        ValueError: 'top' can not be converted to a side
    """
    if side is None: return -1

    if isinstance(side,str):
        if side == '':
            raise ValueError("no empty string for side")
        if 'left'.startswith(side): return 0
        if 'right'.startswith(side): return -1
        raise ValueError("'%s' can not be converted to a side" % (side))

    if isinstance(side, (int,Integer)):
        if side != 0 and side != 1 and side != -1:
            raise ValueError("side must be 0 or 1")
        if side == 0: return 0
        return -1

    raise TypeError("'%s' is not an admissible type" % (str(side)))

def twin_list_iet(a=None):
    r"""
    Returns the twin list of intervals.

    The twin intervals is the correspondance between positions of labels in such
    way that a[interval][position] is a[1-interval][twin[interval][position]]

    INPUT:

    - ``a`` - two lists of labels

    OUTPUT:

    list -- a list of two lists of integers

    TESTS::

        sage: from sage.dynamics.interval_exchanges.template import twin_list_iet
        sage: twin_list_iet([['a','b','c'],['a','b','c']])
        [[0, 1, 2], [0, 1, 2]]
        sage: twin_list_iet([['a','b','c'],['a','c','b']])
        [[0, 2, 1], [0, 2, 1]]
        sage: twin_list_iet([['a','b','c'],['b','a','c']])
        [[1, 0, 2], [1, 0, 2]]
        sage: twin_list_iet([['a','b','c'],['b','c','a']])
        [[2, 0, 1], [1, 2, 0]]
        sage: twin_list_iet([['a','b','c'],['c','a','b']])
        [[1, 2, 0], [2, 0, 1]]
        sage: twin_list_iet([['a','b','c'],['c','b','a']])
        [[2, 1, 0], [2, 1, 0]]
    """
    if a is None : return [[],[]]

    twin = [[0]*len(a[0]), [0]*len(a[1])]

    for i in range(len(twin[0])) :
        c = a[0][i]
        j = a[1].index(c)
        twin[0][i] = j
        twin[1][j] = i

    return twin

def twin_list_li(a=None):
    r"""
    Returns the twin list of intervals

    INPUT:

    - ``a`` - two lists of labels

    OUTPUT:

    list -- a list of two lists of couples of integers

    TESTS::

        sage: from sage.dynamics.interval_exchanges.template import twin_list_li
        sage: twin_list_li([['a','a','b','b'],[]])
        [[(0, 1), (0, 0), (0, 3), (0, 2)], []]
        sage: twin_list_li([['a','a','b'],['b']])
        [[(0, 1), (0, 0), (1, 0)], [(0, 2)]]
        sage: twin_list_li([['a','a'],['b','b']])
        [[(0, 1), (0, 0)], [(1, 1), (1, 0)]]
        sage: twin_list_li([['a'], ['a','b','b']])
        [[(1, 0)], [(0, 0), (1, 2), (1, 1)]]
        sage: twin_list_li([[], ['a','a','b','b']])
        [[], [(1, 1), (1, 0), (1, 3), (1, 2)]]
    """
    if a is None: return [[],[]]

    twin = [
        [(0,j) for j in range(len(a[0]))],
        [(1,j) for j in range(len(a[1]))]]

    for i in (0,1):
        for j in range(len(twin[i])) :
            if twin[i][j] == (i,j) :
                if a[i][j] in a[i][j+1:] :
                # two up or two down
                    j2 = (a[i][j+1:]).index(a[i][j]) + j + 1
                    twin[i][j] = (i,j2)
                    twin[i][j2] = (i,j)
                else :
                    # one up, one down (here i=0)
                    j2 = a[1].index(a[i][j])
                    twin[0][j] = (1,j2)
                    twin[1][j2] = (0,j)

    return twin

class Permutation(SageObject):
    r"""
    Template for all permutations.

    .. WARNING::

        Internal class! Do not use directly!

    This class implement generic algorithm (stratum, connected component, ...)
    and unfies all its children.
    """
    def _repr_(self):
        r"""
        Representation method of self.

        Apply the function str to _repr_type(_repr_options) if _repr_type is
        callable and _repr_type else.

        TESTS:

        ::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p._repr_type = 'str'
            sage: p._repr_options = ('\n',)
            sage: p   #indirect doctest
            a b c
            c b a
            sage: p._repr_options = (' / ',)
            sage: p   #indirect doctest
            a b c / c b a

        ::

            sage: p._repr_type = 'separatrix_diagram'
            sage: p._repr_options = (False,)
            sage: p   #indirect doctest
            [[('c', 'a'), 'b'], ['b', ('c', 'a')]]
            sage: p._repr_options = (True,)
            sage: p
            [[(('c', 'a'), 'L'), ('b', 'R')], [('b', 'L'), (('c', 'a'), 'R')]]

        ::

            sage: p._repr_type = '_twin'
            sage: p   #indirect doctest
            [[2, 1, 0], [2, 1, 0]]
        """
        if self._repr_type is None:
            return ''

        elif self._repr_type == 'reduced':
            return ''.join(map(str,self[1]))

        else:
            f = getattr(self, self._repr_type)
            if callable(f):
                return str(f(*self._repr_options))
            else:
                return str(f)

    def str(self, sep= "\n"):
        r"""
        A string representation of the generalized permutation.

        INPUT:

        - ``sep`` - (default: '\n') a separator for the two intervals

        OUTPUT:

        string -- the string that represents the permutation


        EXAMPLES:

        For permutations of iet::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.str()
            'a b c\nc b a'
            sage: p.str(sep=' | ')
            'a b c | c b a'

        ..the permutation can be rebuilt from the standard string::

            sage: p == iet.Permutation(p.str())
            True

        For permutations of li::

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.str()
            'a b b\nc c a'
            sage: p.str(sep=' | ')
            'a b b | c c a'

        ..the generalized permutation can be rebuilt from the standard string::

            sage: p == iet.GeneralizedPermutation(p.str())
            True

        """
        l = self.list()
        s0 = ' '.join(map(str,l[0]))
        s1 = ' '.join(map(str,l[1]))
        return s0 + sep + s1

    _repr_type = 'str'
    _repr_options = ("\n",)

    def _set_alphabet(self, alphabet):
        r"""
        Sets the alphabet of self.

        TESTS:

            sage: p = iet.GeneralizedPermutation('a a','b b')
            sage: p.alphabet([0,1])   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 0
            1 1
            sage: p.alphabet("cd")   #indirect doctest
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c c
            d d

        Tests with reduced permutations::

            sage: p = iet.Permutation('a b','b a',reduced=True)
            sage: p.alphabet([0,1])   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 1
            1 0
            sage: p.alphabet("cd")   #indirect doctest
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c d
            d c

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)
            sage: p.alphabet([0,1])   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 0
            1 1
            sage: p.alphabet("cd")   #indirect doctest
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c c
            d d
        """
        alphabet = Alphabet(alphabet)
        if alphabet.cardinality() < len(self):
            raise ValueError("Your alphabet has not enough letters")
        self._alphabet = alphabet

    def alphabet(self, data=None):
        r"""
        Manages the alphabet of self.

        If there is no argument, the method returns the alphabet used. If the
        argument could be converted to an alphabet, this alphabet will be used.

        INPUT:

        - ``data`` - None or something that could be converted to an alphabet


        OUTPUT:

        -- either None or the current alphabet


        EXAMPLES::

            sage: p = iet.Permutation('a b','a b')
            sage: p.alphabet([0,1])
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 1
            0 1
            sage: p.alphabet("cd")
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c d
            c d
        """
        if data is None:
            return self._alphabet
        else:
            self._set_alphabet(data)

    def letters(self):
        r"""
        Returns the list of letters of the alphabet used for representation.

        The letters used are not necessarily the whole alphabet (for example if
        the alphabet is infinite).

        OUTPUT:

        -- a list of labels


        EXAMPLES::

            sage: p = iet.Permutation([1,2],[2,1])
            sage: p.alphabet(Alphabet(name="NN"))
            sage: p
            0 1
            1 0
            sage: p.letters()
            [0, 1]
        """
        return [self._alphabet.unrank(_) for _ in range(len(self))]

    def left_right_inverse(self):
        r"""
        Returns the left-right inverse.

        You can also use the shorter .lr_inverse()

        OUTPUT:

        -- a permutation


        EXAMPLES::

            sage: p = iet.Permutation('a b c','c a b')
            sage: p.left_right_inverse()
            c b a
            b a c
            sage: p = iet.Permutation('a b c d','c d a b')
            sage: p.left_right_inverse()
            d c b a
            b a d c

        ::

             sage: p = iet.GeneralizedPermutation('a a','b b c c')
             sage: p.left_right_inverse()
             a a
             c c b b

        ::

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: p.left_right_inverse() == p
            True
            sage: p = iet.Permutation('a b c','c a b',reduced=True)
            sage: q = p.left_right_inverse()
            sage: q == p
            False
            sage: q
            a b c
            b c a

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: p.left_right_inverse() == p
            True
            sage: p = iet.GeneralizedPermutation('a b b','c c a',reduced=True)
            sage: q = p.left_right_inverse()
            sage: q == p
            False
            sage: q
            a a b
            b c c

        TESTS::

            sage: p = iet.GeneralizedPermutation('a a','b b')
            sage: p.left_right_inverse()
            a a
            b b
            sage: p is p.left_right_inverse()
            False
            sage: p == p.left_right_inverse()
            True
        """
        result = copy(self)
        result._reversed()
        return result

    lr_inverse = left_right_inverse
    vertical_inverse = left_right_inverse

    def top_bottom_inverse(self):
        r"""
        Returns the top-bottom inverse.

        You can use also use the shorter .tb_inverse().

        OUTPUT:

        -- a permutation


        EXAMPLES::

            sage: p = iet.Permutation('a b','b a')
            sage: p.top_bottom_inverse()
            b a
            a b
            sage: p = iet.Permutation('a b','b a',reduced=True)
            sage: p.top_bottom_inverse() == p
            True

        ::

            sage: p = iet.Permutation('a b c d','c d a b')
            sage: p.top_bottom_inverse()
            c d a b
            a b c d

        TESTS::

            sage: p = iet.Permutation('a b','a b')
            sage: p == p.top_bottom_inverse()
            True
            sage: p is p.top_bottom_inverse()
            False
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)
            sage: p == p.top_bottom_inverse()
            True
            sage: p is p.top_bottom_inverse()
            False
        """
        result = copy(self)
        result._inversed()
        return result

    tb_inverse = top_bottom_inverse
    horizontal_inverse = top_bottom_inverse

    def symmetric(self):
        r"""
        Returns the symmetric permutation.

        The symmetric permutation is the composition of the top-bottom
        inversion and the left-right inversion (which are geometrically
        orientation reversing).

        OUTPUT:

        -- a permutation


        EXAMPLES::

            sage: p = iet.Permutation("a b c","c b a")
            sage: p.symmetric()
            a b c
            c b a
            sage: q = iet.Permutation("a b c d","b d a c")
            sage: q.symmetric()
            c a d b
            d c b a

        ::

            sage: p = iet.Permutation('a b c d','c a d b')
            sage: q = p.symmetric()
            sage: q1 = p.tb_inverse().lr_inverse()
            sage: q2 = p.lr_inverse().tb_inverse()
            sage: q == q1
            True
            sage: q == q2
            True


        TESTS:

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: q = p.symmetric()
            sage: q1 = p.tb_inverse().lr_inverse()
            sage: q2 = p.lr_inverse().tb_inverse()
            sage: q == q1
            True
            sage: q == q2
            True

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True,flips='a')
            sage: q = p.symmetric()
            sage: q1 = p.tb_inverse().lr_inverse()
            sage: q2 = p.lr_inverse().tb_inverse()
            sage: q == q1
            True
            sage: q == q2
            True

        """
        res = copy(self)
        res._inversed()
        res._reversed()
        return res

    def has_rauzy_move(self, winner='top', side=None):
        r"""
        Tests the legality of a Rauzy move.

        INPUT:

        - ``winner`` - 'top' or 'bottom' corresponding to the interval

        - ``side`` - 'left' or 'right' (default)


        OUTPUT:

        -- a boolean


        EXAMPLES::

            sage: p = iet.Permutation('a b','a b')
            sage: p.has_rauzy_move('top','right')
            False
            sage: p.has_rauzy_move('bottom','right')
            False
            sage: p.has_rauzy_move('top','left')
            False
            sage: p.has_rauzy_move('bottom','left')
            False

        ::

            sage: p = iet.Permutation('a b c','b a c')
            sage: p.has_rauzy_move('top','right')
            False
            sage: p.has_rauzy_move('bottom', 'right')
            False
            sage: p.has_rauzy_move('top','left')
            True
            sage: p.has_rauzy_move('bottom','left')
            True

        ::

            sage: p = iet.Permutation('a b','b a')
            sage: p.has_rauzy_move('top','right')
            True
            sage: p.has_rauzy_move('bottom','right')
            True
            sage: p.has_rauzy_move('top','left')
            True
            sage: p.has_rauzy_move('bottom','left')
            True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)

        if side == -1:
            return self.has_right_rauzy_move(winner)
        elif side == 0:
            return self.lr_inverse().has_right_rauzy_move(winner)

    def rauzy_move(self, winner, side='right', iteration=1):
        r"""
        Returns the permutation after a Rauzy move.

        INPUT:

        - ``winner`` - 'top' or 'bottom' interval

        - ``side`` - 'right' or 'left' (default: 'right') corresponding
          to the side on which the Rauzy move must be performed.

        - ``iteration`` - a non negative integer


        OUTPUT:

        - a permutation


        TESTS::

            sage: p = iet.Permutation('a b','b a')
            sage: p.rauzy_move(winner=0, side='right') == p
            True
            sage: p.rauzy_move(winner=1, side='right') == p
            True
            sage: p.rauzy_move(winner=0, side='left') == p
            True
            sage: p.rauzy_move(winner=1, side='left') == p
            True

        ::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.rauzy_move(winner=0, side='right')
            a b c
            c a b
            sage: p.rauzy_move(winner=1, side='right')
            a c b
            c b a
            sage: p.rauzy_move(winner=0, side='left')
            a b c
            b c a
            sage: p.rauzy_move(winner=1, side='left')
            b a c
            c b a
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)

        if side == -1:
            tmp = self
            for k in range(iteration):
                tmp = tmp.right_rauzy_move(winner)
            return tmp

        elif side == 0:
            tmp = self
            for k in range(iteration):
                tmp = tmp.left_rauzy_move(winner)
            return tmp

class PermutationIET(Permutation):
    """
    Template for permutation from Interval Exchange Transformation.

    .. WARNING::

        Internal class! Do not use directly!

    AUTHOR:

    - Vincent Delecroix (2008-12-20): initial version
    """
    def _init_alphabet(self,a) :
        r"""
        Initializes the alphabet from intervals.

        TESTS::

            sage: p = iet.Permutation('a b c d','d c a b')   #indirect doctest
            sage: p.alphabet() == Alphabet(['a', 'b', 'c', 'd'])
            True
            sage: p = iet.Permutation([0,1,2],[1,0,2],reduced=True)   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1,2])
            True
        """
        self._alphabet = Alphabet(a[0])

    def is_irreducible(self, return_decomposition=False) :
        r"""
        Tests the irreducibility.

        An abelian permutation p = (p0,p1) is reducible if:
        set(p0[:i]) = set(p1[:i]) for an i < len(p0)

        OUTPUT:

        - a boolean

        EXAMPLES::

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.is_irreducible()
            True

            sage: p = iet.Permutation('a b c', 'b a c')
            sage: p.is_irreducible()
            False
        """
        s0, s1 = 0, 0
        for i in range(len(self)-1) :
            s0 += i
            s1 += self._twin[0][i]
            if s0 == s1 :
                if return_decomposition :
                    return False, (self[0][:i+1], self[0][i+1:], self[1][:i+1], self[1][i+1:])
                return False
        if return_decomposition:
            return True, (self[0],[],self[1],[])
        return True

    def decompose(self):
        r"""
        Returns the decomposition of self.

        OUTPUT:

        -- a list of permutations


        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a').decompose()[0]
            sage: p
            a b c
            c b a

        ::

            sage: p1,p2,p3 = iet.Permutation('a b c d e','b a c e d').decompose()
            sage: p1
            a b
            b a
            sage: p2
            c
            c
            sage: p3
            d e
            e d
        """
        l = []
        test, t = self.is_irreducible(return_decomposition=True)
        l.append(self.__class__((t[0],t[2])))

        while not test:
            q = self.__class__((t[1],t[3]))
            test, t = q.is_irreducible(return_decomposition=True)
            l.append(self.__class__((t[0],t[2])))

        return l

    def intersection_matrix(self):
        r"""
        Returns the intersection matrix.

        This `d*d` antisymmetric matrix is given by the rule :

        .. math::

            m_{ij} = \begin{cases}
                1 & \text{$i < j$ and $\pi(i) > \pi(j)$} \\
                -1 & \text{$i > j$ and $\pi(i) < \pi(j)$} \\
                0 & \text{else}
                \end{cases}

        OUTPUT:

        - a matrix

        EXAMPLES::

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p.intersection_matrix()
            [ 0  1  1  1]
            [-1  0  1  1]
            [-1 -1  0  1]
            [-1 -1 -1  0]

        ::

            sage: p = iet.Permutation('1 2 3 4 5','5 3 2 4 1')
            sage: p.intersection_matrix()
            [ 0  1  1  1  1]
            [-1  0  1  0  1]
            [-1 -1  0  0  1]
            [-1  0  0  0  1]
            [-1 -1 -1 -1  0]
        """
        n = self.length_top()
        m = matrix(n)
        for i in range(n):
            for j in range(i,n):
                if self._twin[0][i] > self._twin[0][j]:
                    m[i,j] = 1
                    m[j,i] = -1
        return m

    def attached_out_degree(self):
        r"""
        Returns the degree of the singularity at the left of the interval.

        OUTPUT:

        -- a positive integer


        EXAMPLES::

            sage: p1 = iet.Permutation('a b c d e f g','d c g f e b a')
            sage: p2 = iet.Permutation('a b c d e f g','e d c g f b a')
            sage: p1.attached_out_degree()
            3
            sage: p2.attached_out_degree()
            1
        """
        left_corner = ((self[1][0], self[0][0]), 'L')
        for s in self.separatrix_diagram(side=True):
            if left_corner in s:
                return len(s)//2 - 1

    def attached_in_degree(self):
        r"""
        Returns the degree of the singularity at the right of the interval.

        OUTPUT:

        -- a positive integer


        EXAMPLES::

            sage: p1 = iet.Permutation('a b c d e f g','d c g f e b a')
            sage: p2 = iet.Permutation('a b c d e f g','e d c g f b a')
            sage: p1.attached_in_degree()
            1
            sage: p2.attached_in_degree()
            3
        """
        right_corner = ((self[0][-1], self[1][-1]), 'R')

        for s in self.separatrix_diagram(side=True):
            if right_corner in s:
                return len(s)//2 - 1

    def attached_type(self):
        r"""
        Return the singularity degree attached on the left and the right.

        OUTPUT:

        ``([degre], angle_parity)`` -- if the same singularity is attached on the left and right

        ``([left_degree, right_degree], 0)`` -- the degrees at the left and the right which are different singularitites

        EXAMPLES:

        With two intervals::

            sage: p = iet.Permutation('a b','b a')
            sage: p.attached_type()
            ([0], 1)

        With three intervals::

            sage: p = iet.Permutation('a b c','b c a')
            sage: p.attached_type()
            ([0], 1)

            sage: p = iet.Permutation('a b c','c a b')
            sage: p.attached_type()
            ([0], 1)

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.attached_type()
            ([0, 0], 0)

        With four intervals::

            sage: p = iet.Permutation('1 2 3 4','4 3 2 1')
            sage: p.attached_type()
            ([2], 0)
        """
        left_corner = ((self[1][0], self[0][0]), 'L')
        right_corner = ((self[0][-1], self[1][-1]), 'R')

        l = self.separatrix_diagram(side=True)

        for s in l:
            if left_corner in s and right_corner in s:
                i1 = s.index(left_corner)
                i2 = s.index(right_corner)
                return ([len(s)//2 - 1], ((i2-i1+1)//2) % 2)
            elif left_corner in s:
                left_degree = len(s)//2 - 1
            elif right_corner in s:
                right_degree = len(s)//2 - 1

        return ([left_degree,right_degree], 0)

    def separatrix_diagram(self,side=False):
        r"""
        Returns the separatrix diagram of the permutation.

        INPUT:

        - ``side`` - boolean


        OUTPUT:

        -- a list of lists


        EXAMPLES::

            sage: iet.Permutation([0, 1], [1, 0]).separatrix_diagram()
            [[(1, 0), (1, 0)]]

        ::

            sage: iet.Permutation('a b c d','d c b a').separatrix_diagram()
            [[('d', 'a'), 'b', 'c', ('d', 'a'), 'b', 'c']]
        """
        separatrices = range(len(self)) # bottom intervals
        labels = self[1] # their labels

        singularities = []


        twin = self._twin
        n = len(self)-1

        while separatrices != []:
            start = separatrices.pop(0)
            separatrix = start
            if side:
                singularity = [(labels[start],'L')]
            else:
                singularity = [labels[start]]

            while True:
                if separatrix == 0:
                    separatrix = twin[0][0]
                    if side:
                        a = singularity.pop()[0]
                    else:
                        a = singularity.pop()
                    if side:
                        singularity.append(((a,labels[separatrix]), 'L'))
                    else:
                        singularity.append((a,labels[separatrix]))

                    if separatrix == start:
                        singularities.append(singularity)
                        break

                    del separatrices[separatrices.index(separatrix)]

                else:
                    separatrix -= 1
                    if side:
                        singularity.append((labels[separatrix],'R'))
                    else:
                        singularity.append(labels[separatrix])

                    if separatrix == twin[0][n] :
                        separatrix = n
                        if side:
                            a = singularity.pop()[0]
                        else:
                            a = singularity.pop()
                        if side:
                            singularity.append(((a,labels[separatrix]),'R'))
                        else:
                            singularity.append((a,labels[separatrix]))

                    separatrix = twin[0][twin[1][separatrix]+1]

                if separatrix == start:
                    singularities.append(singularity)
                    break

                elif separatrix != twin[0][0]:
                    del separatrices[separatrices.index(separatrix)]
                    if side:
                        singularity.append((labels[separatrix],'L'))
                    else:
                        singularity.append(labels[separatrix])

        return singularities

    def stratum(self, marked_separatrix='no'):
        r"""
        Returns the strata in which any suspension of this permutation lives.

        OUTPUT:

        - a stratum of Abelian differentials

        EXAMPLES::

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.stratum()
            H(0, 0)

            sage: p = iet.Permutation('a b c d', 'd a b c')
            sage: p.stratum()
            H(0, 0, 0)

            sage: p = iet.Permutation(range(9), [8,5,2,7,4,1,6,3,0])
            sage: p.stratum()
            H(1, 1, 1, 1)

        You can specify that you want to attach the singularity on the left (or
        on the right) with the option marked_separatrix::

            sage: a = 'a b c d e f g h i j'
            sage: b3 = 'd c g f e j i h b a'
            sage: b2 = 'd c e g f j i h b a'
            sage: b1 = 'e d c g f h j i b a'
            sage: p3 = iet.Permutation(a, b3)
            sage: p3.stratum()
            H(3, 2, 1)
            sage: p3.stratum(marked_separatrix='out')
            H^out(3, 2, 1)
            sage: p2 = iet.Permutation(a, b2)
            sage: p2.stratum()
            H(3, 2, 1)
            sage: p2.stratum(marked_separatrix='out')
            H^out(2, 3, 1)
            sage: p1 = iet.Permutation(a, b1)
            sage: p1.stratum()
            H(3, 2, 1)
            sage: p1.stratum(marked_separatrix='out')
            H^out(1, 3, 2)

        AUTHORS:
            - Vincent Delecroix (2008-12-20)
        """
        from sage.dynamics.flat_surfaces.strata import AbelianStratum

        if not self.is_irreducible():
            return [x.stratum(marked_separatrix) for x in self.decompose()]

        if len(self) == 1:
            return AbelianStratum([])

        singularities = [len(x)//2 - 1 for x in self.separatrix_diagram()]

        return AbelianStratum(singularities,marked_separatrix=marked_separatrix)

    def genus(self) :
        r"""
        Returns the genus corresponding to any suspension of the permutation.

        OUTPUT:

        -- a positive integer

        EXAMPLES::

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.genus()
            1

        ::

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p.genus()
            2

        REFERENCES:
            Veech
        """
        return self.stratum().genus()

    def arf_invariant(self):
        r"""
        Returns the Arf invariant of the suspension of self.

        OUTPUT:

        integer -- 0 or 1

        EXAMPLES:

        Permutations from the odd and even component of H(2,2,2)::

            sage: a = range(10)
            sage: b1 = [3,2,4,6,5,7,9,8,1,0]
            sage: b0 = [6,5,4,3,2,7,9,8,1,0]
            sage: p1 = iet.Permutation(a,b1)
            sage: p1.arf_invariant()
            1
            sage: p0 = iet.Permutation(a,b0)
            sage: p0.arf_invariant()
            0

        Permutations from the odd and even component of H(4,4)::

            sage: a = range(11)
            sage: b1 = [3,2,5,4,6,8,7,10,9,1,0]
            sage: b0 = [5,4,3,2,6,8,7,10,9,1,0]
            sage: p1 = iet.Permutation(a,b1)
            sage: p1.arf_invariant()
            1
            sage: p0 = iet.Permutation(a,b0)
            sage: p0.arf_invariant()
            0

        REFERENCES:

        [Jo80] D. Johnson, "Spin structures and quadratic forms on surfaces", J.
        London Math. Soc (2), 22, 1980, 365-373

        [KoZo03] M. Kontsevich, A. Zorich "Connected components of the moduli
        spaces of Abelian differentials with prescribed singularities",
        Inventiones Mathematicae, 153, 2003, 631-678
        """
        M = self.intersection_matrix()
        F, C = M.symplectic_form()

        g = F.rank()/2
        n = F.ncols()

        s = 0
        for i in range(g):
            a = C.row(i)

            a_indices = []
            for k in xrange(n):
                if a[k] != 0: a_indices.append(k)

            t_a = len(a_indices) % 2
            for j1 in xrange(len(a_indices)):
                for j2 in xrange(j1+1,len(a_indices)):
                    t_a = (t_a + M[a_indices[j1], a_indices[j2]]) % 2

            b = C.row(g+i)

            b_indices = []
            for k in xrange(n):
                if b[k] != 0: b_indices.append(k)

            t_b = len(b_indices) % 2
            for j1 in xrange(len(b_indices)):
                for j2 in xrange(j1+1,len(b_indices)):
                    t_b = (t_b + M[b_indices[j1],b_indices[j2]]) % 2

            s = (s + t_a * t_b) % 2

        return s

    def connected_component(self,marked_separatrix='no'):
        r"""
        Returns a connected components of a stratum.

        EXAMPLES:

        Permutations from the stratum H(6)::

            sage: a = range(8)
            sage: b_hyp = [7,6,5,4,3,2,1,0]
            sage: b_odd = [3,2,5,4,7,6,1,0]
            sage: b_even = [5,4,3,2,7,6,1,0]
            sage: p_hyp = iet.Permutation(a, b_hyp)
            sage: p_odd = iet.Permutation(a, b_odd)
            sage: p_even = iet.Permutation(a, b_even)
            sage: p_hyp.connected_component()
            H_hyp(6)
            sage: p_odd.connected_component()
            H_odd(6)
            sage: p_even.connected_component()
            H_even(6)

        Permutations from the stratum H(4,4)::

            sage: a = range(11)
            sage: b_hyp = [10,9,8,7,6,5,4,3,2,1,0]
            sage: b_odd = [3,2,5,4,6,8,7,10,9,1,0]
            sage: b_even = [5,4,3,2,6,8,7,10,9,1,0]
            sage: p_hyp = iet.Permutation(a,b_hyp)
            sage: p_odd = iet.Permutation(a,b_odd)
            sage: p_even = iet.Permutation(a,b_even)
            sage: p_hyp.stratum() == AbelianStratum(4,4)
            True
            sage: p_hyp.connected_component()
            H_hyp(4, 4)
            sage: p_odd.stratum() == AbelianStratum(4,4)
            True
            sage: p_odd.connected_component()
            H_odd(4, 4)
            sage: p_even.stratum() == AbelianStratum(4,4)
            True
            sage: p_even.connected_component()
            H_even(4, 4)

        As for stratum you can specify that you want to attach the singularity
        on the left of the interval using the option marked_separatrix::

            sage: a = [1,2,3,4,5,6,7,8,9]
            sage: b4_odd = [4,3,6,5,7,9,8,2,1]
            sage: b4_even = [6,5,4,3,7,9,8,2,1]
            sage: b2_odd = [4,3,5,7,6,9,8,2,1]
            sage: b2_even = [7,6,5,4,3,9,8,2,1]
            sage: p4_odd = iet.Permutation(a,b4_odd)
            sage: p4_even = iet.Permutation(a,b4_even)
            sage: p2_odd = iet.Permutation(a,b2_odd)
            sage: p2_even = iet.Permutation(a,b2_even)
            sage: p4_odd.connected_component(marked_separatrix='out')
            H_odd^out(4, 2)
            sage: p4_even.connected_component(marked_separatrix='out')
            H_even^out(4, 2)
            sage: p2_odd.connected_component(marked_separatrix='out')
            H_odd^out(2, 4)
            sage: p2_even.connected_component(marked_separatrix='out')
            H_even^out(2, 4)
            sage: p2_odd.connected_component() == p4_odd.connected_component()
            True
            sage: p2_odd.connected_component('out') == p4_odd.connected_component('out')
            False
        """
        from sage.dynamics.flat_surfaces.strata import (HypCCA,
                                                        OddCCA, EvenCCA)

        if not self.is_irreducible():
            return [x.connected_component(marked_separatrix) for x in self.decompose()]

        stratum = self.stratum(marked_separatrix=marked_separatrix)
        cc = stratum._cc

        if len(cc) == 1:
            return stratum.connected_components()[0]

        if HypCCA in cc:
            if self.is_hyperelliptic():
                return HypCCA(stratum)
            else:
                cc = cc[1:]

        if len(cc) == 1:
            return cc[0](stratum)

        else:
            spin = self.arf_invariant()
            if spin == 0:
                return EvenCCA(stratum)
            else:
                return OddCCA(stratum)

    def order_of_rauzy_action(self, winner, side=None):
        r"""
        Returns the order of the action of a Rauzy move.

        INPUT:

        - ``winner`` - string ``'top'`` or ``'bottom'``

        - ``side`` - string ``'left'`` or ``'right'``

        OUTPUT:

        An integer corresponding to the order of the Rauzy action.

        EXAMPLES::

            sage: p = iet.Permutation('a b c d','d a c b')
            sage: p.order_of_rauzy_action('top', 'right')
            3
            sage: p.order_of_rauzy_action('bottom', 'right')
            2
            sage: p.order_of_rauzy_action('top', 'left')
            1
            sage: p.order_of_rauzy_action('bottom', 'left')
            3
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)

        if side == -1:
            return self.length(winner) - self._twin[winner][-1] - 1
        elif side == 0:
            return self._twin[winner][0]

    def erase_marked_points(self):
        r"""
        Returns a permutation equivalent to self but without marked points.

        EXAMPLES::

            sage: a = iet.Permutation('a b1 b2 c d', 'd c b1 b2 a')
            sage: a.erase_marked_points()
            a b1 c d
            d c b1 a
        """
        res = copy(self)

        l = res.list()
        left_corner = ((l[1][0], l[0][0]), 'L')
        right_corner = ((l[0][-1], l[1][-1]), 'R')

        s = res.separatrix_diagram(side=True)
        lengths = [len(_) for _ in s]

        while 2 in lengths:
            if lengths == [2]:
                return res

            i = lengths.index(2)
            t = s[i]
            if t[0] == left_corner or t[0] == right_corner:
                letter = t[1][0]
            else:
                letter = t[0][0]

            res = res.erase_letter(letter)

            l = res.list()

            s = res.separatrix_diagram(side=True)
            lengths = [len(_) for _ in s]

        return res

    def is_hyperelliptic(self):
        r"""
        Returns True if the permutation is in the class of the symmetric
        permutations (with eventual marked points).

        This is equivalent to say that the suspension lives in an hyperelliptic
        stratum of Abelian differentials H_hyp(2g-2) or H_hyp(g-1, g-1) with
        some marked points.

        EXAMPLES::

            sage: iet.Permutation('a b c d','d c b a').is_hyperelliptic()
            True
            sage: iet.Permutation('0 1 2 3 4 5','5 2 1 4 3 0').is_hyperelliptic()
            False

        REFERENCES:

        Gerard Rauzy, "Echanges d'intervalles et transformations induites",
        Acta Arith. 34, no. 3, 203-212, 1980

        M. Kontsevich, A. Zorich "Connected components of the moduli space
        of Abelian differentials with prescripebd singularities" Invent. math.
        153, 631-678 (2003)
        """
        test = self.erase_marked_points()

        n = test.length_top()
        cylindric = test.cylindric()
        return cylindric._twin[0] == range(n-1,-1,-1)

    def cylindric(self):
        r"""
        Returns a permutation in the Rauzy class such that

            twin[0][-1] == 0
            twin[1][-1] == 0

        TESTS::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.cylindric() == p
            True
            sage: p = iet.Permutation('a b c d','b d a c')
            sage: q = p.cylindric()
            sage: q[0][0] == q[1][-1]
            True
            sage: q[1][0] == q[1][0]
            True
        """
        tmp = copy(self)
        n = self.length(0)

        a0 = tmp._twin[0][-1]
        a1 = tmp._twin[1][-1]
        p_min = min(a0,a1)

        while p_min > 0:
            if p_min == a0:
                k_min = min(tmp._twin[1][a0+1:])
                k = n - tmp._twin[1].index(k_min) - 1

                tmp = tmp.rauzy_move(0, iteration=k)

            else:
                k_min = min(tmp._twin[0][a1+1:])
                k = n - tmp._twin[0].index(k_min) - 1

                tmp = tmp.rauzy_move(1, iteration=k)

            a0 = tmp._twin[0][-1]
            a1 = tmp._twin[1][-1]
            p_min = min(a0,a1)

        if a0 == 0:
            k = n - tmp._twin[1].index(0) - 1
            tmp = tmp.rauzy_move(0, iteration = k)

        else:
            k = n - tmp._twin[0].index(0) - 1
            tmp = tmp.rauzy_move(1, iteration=k)

        return tmp

    def is_cylindric(self):
        r"""
        Returns True if the permutation is Rauzy_1n.

        A permutation is cylindric if 1 and n are exchanged.

        EXAMPLES::

            sage: iet.Permutation('1 2 3','3 2 1').is_cylindric()
            True
            sage: iet.Permutation('1 2 3','2 1 3').is_cylindric()
            False
        """
        return self._twin[0][-1] == 0 and self._twin[1][-1] == 0

    def to_permutation(self):
        r"""
        Returns the permutation as an element of the symmetric group.

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.to_permutation()
            [3, 2, 1]

        ::

            sage: p = Permutation([2,4,1,3])
            sage: q = iet.Permutation(p)
            sage: q.to_permutation() == p
            True
        """
        from sage.combinat.permutation import Permutation
        return Permutation([x+1 for x in self._twin[1]])

class PermutationLI(Permutation):
    r"""
    Template for quadratic permutation.

    .. WARNING::

        Internal class! Do not use directly!

    AUTHOR:

    - Vincent Delecroix (2008-12-20): initial version
    """
    def _init_twin(self,a):
        r"""
        Initialization of the twin data.

        TEST::

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)   #indirect doctest
            sage: p._twin
            [[(0, 1), (0, 0)], [(1, 1), (1, 0)]]
        """
        # creation of the twin
        self._twin = [[],[]]
        l = [[(0,j) for j in range(len(a[0]))],[(1,j) for j in range(len(a[1]))]]
        for i in range(2) :
            for j in range(len(l[i])) :
                if l[i][j] == (i,j) :
                    if a[i][j] in a[i][j+1:] :
                        # two up or two down
                        j2 = (a[i][j+1:]).index(a[i][j]) + j + 1
                        l[i][j] = (i,j2)
                        l[i][j2] = (i,j)
                    else :
                        # one up, one down (here i=0)
                        j2 = a[1].index(a[i][j])
                        l[0][j] = (1,j2)
                        l[1][j2] = (0,j)

        self._twin[0] = l[0]
        self._twin[1] = l[1]

    def _init_alphabet(self, intervals) :
        r"""
        Intialization procedure of the alphabet of self from intervals list

        TEST::

            sage: p = iet.GeneralizedPermutation('a a','b b')   #indirect doctest
            sage: p.alphabet()
            {'a', 'b'}
        """
        tmp_alphabet = []
        for letter in intervals[0] + intervals[1] :
            if letter not in tmp_alphabet :
                tmp_alphabet.append(letter)

        self._alphabet = Alphabet(tmp_alphabet)

    def is_irreducible(self, return_decomposition=False):
        r"""
        Test of reducibility

        A quadratic (or generalized) permutation is *reducible* if there exists
        a decomposition

        .. math::

           A1 u B1 | ... | B1 u A2

           A1 u B2 | ... | B2 u A2

        where no corners is empty, or exactly one corner is empty
        and it is on the left, or two and they are both on the
        right or on the left. The definition is due to [BL08]_ where they prove
        that the property of being irreducible is stable under Rauzy induction.

        INPUT:

        -  ``return_decomposition`` - boolean (default: False) - if True, and
           the permutation is reducible, returns also the blocs A1 u B1, B1 u
           A2, A1 u B2 and B2 u A2 of a decomposition as above.

        OUTPUT:

        If return_decomposition is True, returns a 2-uple
        (test,decomposition) where test is the preceding test and
        decomposition is a 4-uple (A11,A12,A21,A22) where:

        A11 = A1 u B1
        A12 = B1 u A2
        A21 = A1 u B2
        A22 = B2 u A2

        EXAMPLES::

            sage: GP = iet.GeneralizedPermutation

            sage: GP('a a','b b').is_irreducible()
            False
            sage: GP('a a b','b c c').is_irreducible()
            True
            sage: GP('1 2 3 4 5 1','5 6 6 4 3 2').is_irreducible()
            True

        TESTS:

        Test reducible permutations with no empty corner::

            sage: GP('1 4 1 3','4 2 3 2').is_irreducible(True)
            (False, (['1', '4'], ['1', '3'], ['4', '2'], ['3', '2']))

        Test reducible permutations with one left corner empty::

            sage: GP('1 2 2 3 1','4 4 3').is_irreducible(True)
            (False, (['1'], ['3', '1'], [], ['3']))
            sage: GP('4 4 3','1 2 2 3 1').is_irreducible(True)
            (False, ([], ['3'], ['1'], ['3', '1']))

        Test reducible permutations with two left corners empty::

            sage: GP('1 1 2 3','4 2 4 3').is_irreducible(True)
            (False, ([], ['3'], [], ['3']))

        Test reducible permutations with two right corners empty::

            sage: GP('1 2 2 3 3','1 4 4').is_irreducible(True)
            (False, (['1'], [], ['1'], []))
            sage: GP('1 2 2','1 3 3').is_irreducible(True)
            (False, (['1'], [], ['1'], []))
            sage: GP('1 2 3 3','2 1 4 4 5 5').is_irreducible(True)
            (False, (['1', '2'], [], ['2', '1'], []))

        AUTHORS:

            - Vincent Delecroix (2008-12-20)
        """
        l0 = self.length_top()
        l1 = self.length_bottom()
        s0,s1 = self.list()

        # testing two corners empty on the right (i12 = i22 = 0)
        A11, A21, A12, A22 = [],[],[],[]

        for i11 in range(1, l0):
            if s0[i11-1] in A11:
                break
            A11 = s0[:i11]

            for i21 in range(1, l1):
                if s1[i21-1] in A21:
                    break
                A21 = s1[:i21]

                if sorted(A11) == sorted(A21):
                    if return_decomposition:
                        return False,(A11,A12,A21,A22)
                    return False
            A21 = []

        # testing no corner empty but one or two on the left
        A11, A12, A21, A22 = [], [], [], []
        for i11 in range(0, l0):
            if i11 > 0 and s0[i11-1] in A11:
                break
            A11 = s0[:i11]

            for i21 in xrange(0, l1) :
                if i21 > 0 and s1[i21-1] in A21:
                    break
                A21 = s1[:i21]

                for i12 in xrange(l0 - 1, i11 - 1, -1) :
                    if s0[i12] in A12 or s0[i12] in A21:
                        break
                    A12 = s0[i12:]

                    for i22 in xrange(l1 - 1, i21 - 1, -1) :
                        if s1[i22] in A22 or s1[i22] in A11:
                            break
                        A22 = s1[i22:]

                        if sorted(A11 + A22) == sorted(A12 + A21) :
                            if return_decomposition :
                                return False, (A11,A12,A21,A22)
                            return False
                    A22 = []
                A12 = []
            A21 = []


        if return_decomposition:
            return True, ()
        return True

    def has_right_rauzy_move(self, winner):
        r"""
        Test of Rauzy movability (with an eventual specified choice of winner)

        A quadratic (or generalized) permutation is rauzy_movable type
        depending on the possible length of the last interval. It's
        dependent of the length equation.

        INPUT:

        - ``winner`` - the integer 'top' or 'bottom'

        EXAMPLES::

            sage: p = iet.GeneralizedPermutation('a a','b b')
            sage: p.has_right_rauzy_move('top')
            False
            sage: p.has_right_rauzy_move('bottom')
            False

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.has_right_rauzy_move('top')
            True
            sage: p.has_right_rauzy_move('bottom')
            True

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: p.has_right_rauzy_move('top')
            True
            sage: p.has_right_rauzy_move('bottom')
            False

        ::

            sage: p = iet.GeneralizedPermutation('a a b b','c c')
            sage: p.has_right_rauzy_move('top')
            False
            sage: p.has_right_rauzy_move('bottom')
            True
        """
        winner = interval_conversion(winner)
        loser = 1 - winner

        # the same letter at the right-end (False)
        if (self._twin[0][-1][0] == 1 and
            self._twin[0][-1][1] == self.length_bottom() - 1):
            return False

        # the winner (or loser) letter is repeated on the other interval (True)
        if self._twin[winner][-1][0] == loser:
            return True
        if self._twin[loser][-1][0] == winner:
            return True

        # the loser letters is the only letter repeated in
        # the loser interval (False)
        if [i for i,_ in self._twin[loser]].count(loser) == 2:
            return False

        return True

def labelize_flip(couple):
    r"""
    Returns a string from a 2-uple couple of the form (name, flip).

    TESTS::

        sage: from sage.dynamics.interval_exchanges.template import labelize_flip
        sage: labelize_flip((0,1))
        ' 0'
        sage: labelize_flip((0,-1))
        '-0'
    """
    if couple[1] == -1: return '-' + str(couple[0])
    return ' ' + str(couple[0])

class FlippedPermutation(Permutation):
    r"""
    Template for flipped generalized permutations.

    .. WARNING::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version
    """
    def _init_flips(self,intervals,flips):
        r"""
        Initialize the flip list

        TESTS:

            sage: iet.Permutation('a b','b a',flips='a').flips() #indirect doctest
            ['a']
            sage: iet.Permutation('a b','b a',flips='b').flips() #indirect doctest
            ['b']
            sage: iet.Permutation('a b','b a',flips='ab').flips() #indirect doctest
            ['a', 'b']

        ::

            sage: iet.GeneralizedPermutation('a a','b b',flips='a').flips()
            ['a']
            sage: iet.GeneralizedPermutation('a a','b b',flips='b').flips()
            ['b']
            sage: iet.GeneralizedPermutation('a a','b b',flips='ab').flips()
            ['a', 'b']
        """
        self._flips = [[1]*self.length_top(), [1]*self.length_bottom()]
        for interval in (0,1):
            for i,letter in enumerate(intervals[interval]):
                if letter in flips:
                    self._flips[interval][i] = -1

    def str(self, sep="\n"):
        r"""
        String representation.

        TESTS::

            sage: p = iet.GeneralizedPermutation('a a','b b',flips='a')
            sage: print(p.str())
            -a -a
             b  b
             sage: print(p.str('/'))
             -a -a/ b  b
        """
        l = self.list(flips=True)
        return (' '.join(map(labelize_flip, l[0]))
                + sep
                + ' '.join(map(labelize_flip, l[1])))


class FlippedPermutationIET(FlippedPermutation, PermutationIET):
    r"""
    Template for flipped Abelian permutations.

    .. WARNING::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version
    """
    def flips(self):
        r"""
        Returns the list of flips.

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a',flips='ac')
            sage: p.flips()
            ['a', 'c']
        """
        result = []
        l = self.list(flips=False)
        for i,f in enumerate(self._flips[0]):
            if f == -1:
                result.append(l[0][i])
        return result

class FlippedPermutationLI(FlippedPermutation, PermutationLI):
    r"""
    Template for flipped quadratic permutations.

    .. WARNING::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version
    """
    def flips(self):
        r"""
        Returns the list of flipped intervals.

        EXAMPLES::

            sage: p = iet.GeneralizedPermutation('a a','b b',flips='a')
            sage: p.flips()
            ['a']
            sage: p = iet.GeneralizedPermutation('a a','b b',flips='b',reduced=True)
            sage: p.flips()
            ['b']
        """
        res = []
        l = self.list(flips=False)
        for i,f in enumerate(self._flips[0]):
            if f == -1:
                res.append(l[0][i])
        for i,f in enumerate(self._flips[1]):
            if f == -1:
                res.append(l[1][i])
        return list(set(res))


class RauzyDiagram(SageObject):
    r"""
    Template for Rauzy diagrams.

    .. WARNING::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version
    """
    # TODO: pickle problem of Path (it does not understand what is its parent)
    __metaclass__ = NestedClassMetaclass

    class Path(SageObject):
        r"""
        Path in Rauzy diagram.

            A path in a Rauzy diagram corresponds to a subsimplex of the simplex of
            lengths. This correspondance is obtained via the Rauzy induction. To a
            idoc IET we can associate a unique path in a Rauzy diagram. This
            establishes a correspondance between infinite full path in Rauzy diagram
            and equivalence topologic class of IET.
        """
        def __init__(self, parent, *data):
            r"""
            Constructor of the path.

            TEST::

                sage: p = iet.Permutation('a b c', 'c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p, 0, 1, 0); g
                Path of length 3 in a Rauzy diagram

            Check for :trac:`8388`::

                sage: loads(dumps(g)) == g
                True
            """
            self._parent = parent

            if len(data) == 0:
                raise ValueError("No empty data")

            start = data[0]
            if start not in self._parent:
                raise ValueError("Starting point not in this Rauzy diagram")

            self._start = self._parent._permutation_to_vertex(start)

            cur_vertex = self._start
            self._edge_types = []

            n = len(self._parent._edge_types)
            for i in data[1:]:
                if not isinstance(i, (int,Integer)): # try parent method
                    i = self._parent.edge_types_index(i)

                if i < 0 or i > n:
                    raise ValueError("indices must be integer between 0 and %d" % (n))
                neighbours = self._parent._succ[cur_vertex]
                if neighbours[i] is None:
                    raise ValueError("Invalid path")

                cur_vertex = neighbours[i]
                self._edge_types.append(i)

            self._end = cur_vertex

        def _repr_(self):
            r"""
            Returns a representation of the path.

            TEST::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p)   #indirect doctest
                Path of length 0 in a Rauzy diagram
                sage: r.path(p,'top')   #indirect doctest
                Path of length 1 in a Rauzy diagram
                sage: r.path(p,'bottom')   #indirect doctest
                Path of length 1 in a Rauzy diagram
            """
            return "Path of length %d in a Rauzy diagram" % (len(self))

        def start(self):
            r"""
            Returns the first vertex of the path.

            EXAMPLES::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p, 't', 'b')
                sage: g.start() == p
                True
            """
            return self._parent._vertex_to_permutation(self._start)

        def end(self):
            r"""
            Returns the last vertex of the path.

            EXAMPLES::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g1 = r.path(p, 't', 'b', 't')
                sage: g1.end() == p
                True
                sage: g2 = r.path(p, 'b', 't', 'b')
                sage: g2.end() == p
                True
            """
            return self._parent._vertex_to_permutation(self._end)

        def edge_types(self):
            r"""
            Returns the edge types of the path.

            EXAMPLES::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p, 0, 1)
                sage: g.edge_types()
                [0, 1]
            """
            return copy(self._edge_types)

        def __eq__(self, other):
            r"""
            Tests equality

            TEST::

                sage: p1 = iet.Permutation('a b','b a')
                sage: r1 = p1.rauzy_diagram()
                sage: p2 = p1.reduced()
                sage: r2 = p2.rauzy_diagram()
                sage: r1.path(p1,0,1) == r2.path(p2,0,1)
                False
                sage: r1.path(p1,0) == r1.path(p1,0)
                True
                sage: r1.path(p1,1) == r1.path(p1,0)
                False
            """
            return (
                type(self) is type(other) and
                self._parent == other._parent and
                self._start == other._start and
                self._edge_types == other._edge_types)

        def __ne__(self,other):
            r"""
            Tests inequality

            TEST::

                sage: p1 = iet.Permutation('a b','b a')
                sage: r1 = p1.rauzy_diagram()
                sage: p2 = p1.reduced()
                sage: r2 = p2.rauzy_diagram()
                sage: r1.path(p1,0,1) != r2.path(p2,0,1)
                True
                sage: r1.path(p1,0) != r1.path(p1,0)
                False
                sage: r1.path(p1,1) != r1.path(p1,0)
                True
            """
            return (
                type(self) is not type(other) or
                self._parent != other._parent or
                self._start != other._start or
                self._edge_types != other._edge_types)

        def __copy__(self):
            r"""
            Returns a copy of the path.

            TESTS::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g1 = r.path(p,0,1,0,0)
                sage: g2 = copy(g1)
                sage: g1 is g2
                False
            """
            res = self.__class__(self._parent, self.start())
            res._edge_types = copy(self._edge_types)
            res._end = copy(self._end)
            return res

        def pop(self):
            r"""
            Pops the queue of the path

            OUTPUT:

            a path corresponding to the last edge

            EXAMPLES::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p,0,1,0)
                sage: g0,g1,g2,g3 = g[0], g[1], g[2], g[3]
                sage: g.pop() == r.path(g2,0)
                True
                sage: g == r.path(g0,0,1)
                True
                sage: g.pop() == r.path(g1,1)
                True
                sage: g == r.path(g0,0)
                True
                sage: g.pop() == r.path(g0,0)
                True
                sage: g == r.path(g0)
                True
                sage: g.pop() == r.path(g0)
                True
            """
            if len(self) == 0:
                return self._parent.path(self.start())

            else:
                e = self._edge_types.pop()
                self._end = self._parent._pred[self._end][e]
                return self._parent.path(self.end(),e)

        def append(self, edge_type):
            r"""
            Append an edge to the path.

            EXAMPLES::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p)
                sage: g.append('top')
                sage: g
                Path of length 1 in a Rauzy diagram
                sage: g.append('bottom')
                sage: g
                Path of length 2 in a Rauzy diagram
            """
            if not isinstance(edge_type, (int,Integer)):
                edge_type = self._parent.edge_types_index(edge_type)

            elif edge_type < 0 or edge_type >= len(self._parent._edge_types):
                raise ValueError("Edge type not valid")

            if self._parent._succ[self._end][edge_type] is None:
                raise ValueError("%d is not a valid edge" % (edge_type))

            self._edge_types.append(edge_type)
            self._end = self._parent._succ[self._end][edge_type]

        def _fast_append(self, edge_type):
            r"""
            Append an edge to the path without verification.

            EXAMPLES::

                sage: p = iet.GeneralizedPermutation('a a','b b c c')
                sage: r = p.rauzy_diagram()

            .. try to add 1 with append::

                sage: g = r.path(p)
                sage: r[p][1] is None
                True
                sage: g.append(1)
                Traceback (most recent call last):
                ...
                ValueError: 1 is not a valid edge

            .. the same with fast append::

                sage: g = r.path(p)
                sage: r[p][1] is None
                True
                sage: g._fast_append(1)
            """
            self._edge_types.append(edge_type)
            self._end = self._parent._succ[self._end][edge_type]

        def extend(self, path):
            r"""
            Extends self with another path.

            EXAMPLES::

                sage: p = iet.Permutation('a b c d','d c b a')
                sage: r = p.rauzy_diagram()
                sage: g1 = r.path(p,'t','t')
                sage: g2 = r.path(p.rauzy_move('t',iteration=2),'b','b')
                sage: g = r.path(p,'t','t','b','b')
                sage: g == g1 + g2
                True
                sage: g = copy(g1)
                sage: g.extend(g2)
                sage: g == g1 + g2
                True
            """
            if self._parent != path._parent:
                raise ValueError("Not on the same Rauzy diagram")

            if self._end != path._start:
                raise ValueError("The end of the first path must the start of the second")

            self._edge_types.extend(path._edge_types)
            self._end = path._end

        def _fast_extend(self, path):
            r"""
            Extension with no verification.

            EXAMPLES::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: p0, p1 = r[p]
                sage: g = r.path(p)
                sage: g._fast_extend(r.path(p0))
                sage: g
                Path of length 0 in a Rauzy diagram
                sage: g._fast_extend(r.path(p1))
                sage: g
                Path of length 0 in a Rauzy diagram
            """
            self._edge_types.extend(path._edge_types)
            self._end = path._end

        def __len__(self):
            r"""
            Returns the length of the path.

            TEST::


                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: len(r.path(p))
                0
                sage: len(r.path(p,0))
                1
                sage: len(r.path(p,1))
                1
            """
            return len(self._edge_types)

        def __getitem__(self, i):
            r"""
            TESTS::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p,'t','b')
                sage: g[0] == p
                True
                sage: g[1] == p.rauzy_move('t')
                True
                sage: g[2] == p.rauzy_move('t').rauzy_move('b')
                True
                sage: g[-1] == g[2]
                True
                sage: g[-2] == g[1]
                True
                sage: g[-3] == g[0]
                True
            """
            if i > len(self) or i < -len(self)-1:
                raise IndexError("path index out of range")

            if i == 0: return self.start()
            if i < 0: i = i + len(self) + 1

            v = self._start
            for k in range(i):
                v = self._parent._succ[v][self._edge_types[k]]
            return self._parent._vertex_to_permutation(v)

        def __add__(self, other):
            r"""
            Concatenation of paths.

            EXAMPLES::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p) + r.path(p,'b') == r.path(p,'b')
                True
                sage: r.path(p,'b') + r.path(p) == r.path(p,'b')
                True
                sage: r.path(p,'t') + r.path(p,'b') == r.path(p,'t','b')
                True
            """
            if self._end != other._start:
                raise ValueError("The end of the first path is not the start of the second")

            res = copy(self)
            res._fast_extend(other)
            return res

        def __mul__(self, n):
            r"""
            Multiple of a loop.

            EXAMPLES::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: l = r.path(p,'b')
                sage: l * 2 == r.path(p,'b','b')
                True
                sage: l * 3 == r.path(p,'b','b','b')
                True
            """
            if not self.is_loop():
                raise ValueError("Must be a loop to have multiple")

            if not isinstance(n, (Integer,int)):
                raise TypeError("The multiplier must be an integer")

            if n < 0:
                raise ValueError("The multiplier must be non negative")

            res = copy(self)
            for i in range(n-1):
                res += self
            return res

        def is_loop(self):
            r"""
            Tests whether the path is a loop (start point = end point).

            EXAMPLES::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p).is_loop()
                True
                sage: r.path(p,0,1,0,0).is_loop()
                True
            """
            return self._start == self._end

        def winners(self):
            r"""
            Returns the winner list associated to the edge of the path.

            EXAMPLES::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p).winners()
                []
                sage: r.path(p,0).winners()
                ['b']
                sage: r.path(p,1).winners()
                ['a']
            """
            return self.composition(
                self._parent.edge_to_winner,
                list.__add__)

        def losers(self):
            r"""
            Returns a list of the loosers on the path.

            EXAMPLES::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g0 = r.path(p,'t','b','t')
                sage: g0.losers()
                ['a', 'c', 'b']
                sage: g1 = r.path(p,'b','t','b')
                sage: g1.losers()
                ['c', 'a', 'b']
            """
            return self.composition(
                self._parent.edge_to_loser,
                list.__add__)

        def __iter__(self):
            r"""
            Iterator over the permutations of the path.

            EXAMPLES::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p)
                sage: for q in g:
                ....:     print(p)
                a b c
                c b a
                sage: g = r.path(p, 't', 't')
                sage: for q in g:
                ....:     print(q)
                ....:     print("*****")
                a b c
                c b a
                *****
                a b c
                c a b
                *****
                a b c
                c b a
                *****
                sage: g = r.path(p,'b','t')
                sage: for q in g:
                ....:     print(q)
                ....:     print("*****")
                a b c
                c b a
                *****
                a c b
                c b a
                *****
                a c b
                c b a
                *****
            """
            i = self._start

            for edge_type in self._edge_types:
                yield self._parent._vertex_to_permutation(i)
                i = self._parent._succ[i][edge_type]

            yield self.end()

        def composition(self, function, composition = None):
            r"""
            Compose an edges function on a path

            INPUT:

            - ``path`` - either a Path or a tuple describing a path

            - ``function`` - function must be of the form

            - ``composition`` - the composition function

            AUTHOR:

            - Vincent Delecroix (2009-09-29)

            EXAMPLES::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: def f(i,t):
                ....:     if t is None: return []
                ....:     return [t]
                sage: g = r.path(p)
                sage: g.composition(f,list.__add__)
                []
                sage: g = r.path(p,0,1)
                sage: g.composition(f, list.__add__)
                [0, 1]
            """
            result = function(None,None)
            cur_vertex = self._start
            p = self._parent._element

            if composition is None: composition = result.__class__.__mul__

            for i in self._edge_types:
                self._parent._set_element(cur_vertex)
                result = composition(result, function(p,i))
                cur_vertex = self._parent._succ[cur_vertex][i]

            return result

        def right_composition(self, function, composition = None) :
            r"""
            Compose an edges function on a path

            INPUT:

            - ``function`` - function must be of the form (indice,type) -> element. Moreover function(None,None) must be an identity element for initialization.

            - ``composition`` - the composition function for the function. * if None (default None)

            TEST::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: def f(i,t):
                ....:     if t is None: return []
                ....:     return [t]
                sage: g = r.path(p)
                sage: g.right_composition(f,list.__add__)
                []
                sage: g = r.path(p,0,1)
                sage: g.right_composition(f, list.__add__)
                [1, 0]
            """
            result = function(None,None)
            p = self._parent._element
            cur_vertex = self._start

            if composition is None: composition = result.__class__.__mul__

            for i in self._edge_types:
                self._parent._set_element(cur_vertex)
                result = composition(function(p,i),result)
                cur_vertex = self._parent._succ[cur_vertex][i]

            return result

    def __init__(self, p,
                 right_induction=True,
                 left_induction=False,
                 left_right_inversion=False,
                 top_bottom_inversion=False,
                 symmetric=False):
        r"""
        - ``self._succ`` contains successors

        - ``self._pred`` contains predecessors

        - ``self._element_class`` is the class of elements of ``self``

        - ``self._element`` is an instance of this class (hence
          contains the alphabet, the representation mode, ...). It is
          used to store data about property of permutations and also as
          a fast iterator.

        INPUT:

        - ``right_induction`` - boolean or 'top' or 'bottom': consider the
          right induction

        - ``left_induction`` - boolean or 'top' or 'bottom': consider the
          left induction

        - ``left_right_inversion`` - consider the left right inversion

        - ``top_bottom_inversion`` - consider the top bottom inversion

        - ``symmetric`` - consider the symmetric

        TESTS::

            sage: r1 = iet.RauzyDiagram('a b','b a')
            sage: r2 = loads(dumps(r1))
        """
        self._edge_types = []
        self._index = {}

        if right_induction is True:
            self._index['rt_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(0,-1)))
            self._index['rb_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(1,-1)))

        elif isinstance(right_induction, str):
            if right_induction == '':
                raise ValueError("right_induction can not be empty string")

            elif 'top'.startswith(right_induction):
                self._index['rt_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(0,-1)))

            elif 'bottom'.startswith(right_induction):
                self._index['rb_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(1,-1)))

            else:
                raise ValueError("%s is not valid for right_induction" % (right_induction))

        if left_induction is True:
            self._index['lt_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(0,0)))
            self._index['lb_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(1,0)))

        elif isinstance(left_induction,str):
            if left_induction == '':
                raise ValueError("left_induction can not be empty string")

            elif 'top'.startswith(left_induction):
                self._index['lt_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move', (0,0)))

            elif 'bottom'.startswith(left_induction):
                self._index['lb_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move', (1,0)))

            else:
                raise ValueError("%s is not valid for left_induction" % (right_induction))

        if left_right_inversion is True:
            self._index['lr_inverse'] = len(self._edge_types)
            self._edge_types.append(('left_right_inverse', ()))

        if top_bottom_inversion is True:
            self._index['tb_inverse'] = len(self._edge_types)
            self._edge_types.append(('top_bottom_inverse', ()))

        if symmetric is True:
            self._index['symmetric'] = len(self._edge_types)
            self._edge_types.append(('symmetric', ()))

        self._n = len(p)
        self._element_class = p.__class__
        self._element = copy(p)
        self._alphabet = self._element._alphabet

        self._pred = {}
        self._succ = {}

        self.complete(p)

    def __eq__(self, other):
        r"""
        Tests equality.

        TESTS:

        ::

            sage: iet.RauzyDiagram('a b','b a') == iet.RauzyDiagram('a b c','c b a')
            False
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r1 = iet.RauzyDiagram('a c b','c b a', alphabet='abc')
            sage: r2 = iet.RauzyDiagram('a b c','c a b', alphabet='abc')
            sage: r == r1
            True
            sage: r == r2
            True
            sage: r1 == r2
            True

        ::

            sage: r = iet.RauzyDiagram('a b c d','d c b a')
            sage: for p in r:
            ....:     p.rauzy_diagram() == r
            True
            True
            True
            True
            True
            True
            True
        """
        return (
            type(self) is type(other) and
            self._edge_types == other._edge_types and
            self._succ.keys()[0] in other._succ)

    def __ne__(self, other):
        r"""
        Tests difference.


        TEST::

            sage: iet.RauzyDiagram('a b','b a') != iet.RauzyDiagram('a b c','c b a')
            True
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r1 = iet.RauzyDiagram('a c b','c b a', alphabet='abc')
            sage: r2 = iet.RauzyDiagram('a b c','c a b', alphabet='abc')
            sage: r != r1
            False
            sage: r != r2
            False
            sage: r1 != r2
            False
        """
        return (
            type(self) is not type(other) or
            self._edge_types != other._edge_types or
            self._succ.keys()[0] not in other._succ)

    def vertices(self):
        r"""
        Returns a list of the vertices.

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r.vertices(): print(p)
            a b
            b a
        """
        return [self._vertex_to_permutation(x) for x in self._succ.keys()]

    def vertex_iterator(self):
        r"""
        Returns an iterator over the vertices

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r.vertex_iterator(): print(p)
            a b
            b a

        ::

            sage: r = iet.RauzyDiagram('a b c d','d c b a')
            sage: from itertools import ifilter
            sage: r_1n = ifilter(lambda x: x.is_cylindric(), r)
            sage: for p in r_1n: print(p)
            a b c d
            d c b a
        """
        from itertools import imap
        return imap(
            lambda x: self._vertex_to_permutation(x),
            self._succ.keys())

    def edges(self,labels=True):
        r"""
        Returns a list of the edges.

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: len(r.edges())
            2
        """
        return list(self.edge_iterator())

    def edge_iterator(self):
        r"""
        Returns an iterator over the edges of the graph.

        EXAMPLES::

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: for e in r.edge_iterator():
            ....:     print(e[0].str(sep='/') + ' --> ' + e[1].str(sep='/'))
            a b/b a --> a b/b a
            a b/b a --> a b/b a
        """
        for x in self._succ.keys():
            for i,y in enumerate(self._succ[x]):
                if y is not None:
                    yield(
                        (self._vertex_to_permutation(x),
                         self._vertex_to_permutation(y),
                         i))

    def edge_types_index(self, data):
        r"""
        Try to convert the data as an edge type.

        INPUT:

        - ``data`` - a string

        OUTPUT:

        integer

        EXAMPLES:

        For a standard Rauzy diagram (only right induction) the 0 index
        corresponds to the 'top' induction and the index 1 corresponds to the
        'bottom' one::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: r.edge_types_index('top')
            0
            sage: r[p][0] == p.rauzy_move('top')
            True
            sage: r.edge_types_index('bottom')
            1
            sage: r[p][1] == p.rauzy_move('bottom')
            True

        The special operations (inversion and symmetry) always appears after the
        different Rauzy inductions::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(symmetric=True)
            sage: r.edge_types_index('symmetric')
            2
            sage: r[p][2] == p.symmetric()
            True

        This function always try to resolve conflictuous name. If it's
        impossible a ValueError is raised::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(left_induction=True)
            sage: r.edge_types_index('top')
            Traceback (most recent call last):
            ...
            ValueError: left and right inductions must be differentiated
            sage: r.edge_types_index('top_right')
            0
            sage: r[p][0] == p.rauzy_move(0)
            True
            sage: r.edge_types_index('bottom_left')
            3
            sage: r[p][3] == p.rauzy_move('bottom', 'left')
            True

        ::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(left_right_inversion=True,top_bottom_inversion=True)
            sage: r.edge_types_index('inversion')
            Traceback (most recent call last):
            ...
            ValueError: left-right and top-bottom inversions must be differentiated
            sage: r.edge_types_index('lr_inverse')
            2
            sage: p.lr_inverse() == r[p][2]
            True
            sage: r.edge_types_index('tb_inverse')
            3
            sage: p.tb_inverse() == r[p][3]
            True

        Short names are accepted::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(right_induction='top',top_bottom_inversion=True)
            sage: r.edge_types_index('top_rauzy_move')
            0
            sage: r.edge_types_index('t')
            0
            sage: r.edge_types_index('tb')
            1
            sage: r.edge_types_index('inversion')
            1
            sage: r.edge_types_index('inverse')
            1
            sage: r.edge_types_index('i')
            1
         """
        if not isinstance(data,str):
            raise ValueError("the edge type must be a string")

        if 'top_rauzy_move'.startswith(data) or 't_rauzy_move'.startswith(data):
            if 'lt_rauzy' in self._index:
                if 'rt_rauzy' in self._index:
                    raise ValueError("left and right inductions must "
                                     "be differentiated")
                return self._index['lt_rauzy']

            if 'rt_rauzy' in self._index:
                return self._index['rt_rauzy']

            raise ValueError("no top induction in this Rauzy diagram")

        if ('bottom_rauzy_move'.startswith(data) or
            'b_rauzy_move'.startswith(data)):
            if 'lb_rauzy' in self._index:
                if 'rb_rauzy' in self._index:
                    raise ValueError("left and right inductions must "
                                     "be differentiated")
                return self._index['lb_rauzy']

            if 'rb_rauzy' in self._index:
                return self._index['rb_rauzy']

            raise ValueError("no bottom Rauzy induction in this diagram")

        if ('left_rauzy_move'.startswith(data) or
            'l_rauzy_move'.startswith(data)):
            if 'lt_rauzy' in self._index:
                if 'lb_rauzy' in self._index:
                    raise ValueError("top and bottom inductions must be differentiated")
                return self._index['lt_rauzy']

            if 'lb_rauzy' in self._index:
                return self._index('lb_rauzy')

            raise ValueError("no left Rauzy induction in this diagram")

        if ('lt_rauzy_move'.startswith(data) or
            'tl_rauzy_move'.startswith(data) or
            'left_top_rauzy_move'.startswith(data) or
            'top_left_rauzy_move'.startswith(data)):
            if not 'lt_rauzy' in self._index:
                raise ValueError("no top-left Rauzy induction in this diagram")
            else:
                return self._index['lt_rauzy']

        if ('lb_rauzy_move'.startswith(data) or
            'bl_rauzy_move'.startswith(data) or
            'left_bottom_rauzy_move'.startswith(data) or
            'bottom_left_rauzy_move'.startswith(data)):
            if not 'lb_rauzy' in self._index:
                raise ValueError("no bottom-left Rauzy induction in this diagram")
            else:
                return self._index['lb_rauzy']

        if 'right'.startswith(data):
            raise ValueError("ambiguity with your edge name: %s" % (data))

        if ('rt_rauzy_move'.startswith(data) or
            'tr_rauzy_move'.startswith(data) or
            'right_top_rauzy_move'.startswith(data) or
            'top_right_rauzy_move'.startswith(data)):
            if not 'rt_rauzy' in self._index:
                raise ValueError("no top-right Rauzy induction in this diagram")
            else:
                return self._index['rt_rauzy']

        if ('rb_rauzy_move'.startswith(data) or
            'br_rauzy_move'.startswith(data) or
            'right_bottom_rauzy_move'.startswith(data) or
            'bottom_right_rauzy_move'.startswith(data)):
            if not 'rb_rauzy' in self._index:
                raise ValueError("no bottom-right Rauzy induction in this diagram")
            else:
                return self._index['rb_rauzy']

        if 'symmetric'.startswith(data):
            if not 'symmetric' in self._index:
                raise ValueError("no symmetric in this diagram")
            else:
                return self._index['symmetric']

        if 'inversion'.startswith(data) or data == 'inverse':
            if 'lr_inverse' in self._index:
                if 'tb_inverse' in self._index:
                    raise ValueError("left-right and top-bottom inversions must be differentiated")
                return self._index['lr_inverse']

            if 'tb_inverse' in self._index:
                return self._index['tb_inverse']

            raise ValueError("no inversion in this diagram")

        if ('lr_inversion'.startswith(data) or
            data == 'lr_inverse' or
            'left_right_inversion'.startswith(data) or
            data == 'left_right_inverse'):
            if not 'lr_inverse' in self._index:
                raise ValueError("no left-right inversion in this diagram")
            else:
                return self._index['lr_inverse']

        if ('tb_inversion'.startswith(data) or
            data == 'tb_inverse' or
            'top_bottom_inversion'.startswith(data)
            or data == 'top_bottom_inverse'):
            if not 'tb_inverse' in self._index:
                raise ValueError("no top-bottom inversion in this diagram")
            else:
                return self._index['tb_inverse']

        raise ValueError("this edge type does not exist: %s" % (data))

    def edge_types(self):
        r"""
        Print information about edges.

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b', 'b a')
            sage: r.edge_types()
            0: rauzy_move(0, -1)
            1: rauzy_move(1, -1)

        ::

            sage: r = iet.RauzyDiagram('a b', 'b a', left_induction=True)
            sage: r.edge_types()
            0: rauzy_move(0, -1)
            1: rauzy_move(1, -1)
            2: rauzy_move(0, 0)
            3: rauzy_move(1, 0)

        ::

            sage: r = iet.RauzyDiagram('a b',' b a',symmetric=True)
            sage: r.edge_types()
            0: rauzy_move(0, -1)
            1: rauzy_move(1, -1)
            2: symmetric()
        """
        for i,(edge_type,t) in enumerate(self._edge_types):
            print(str(i) + ": " + edge_type + str(t))

    def alphabet(self, data=None):
        r"""
        TESTS::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.alphabet() == Alphabet(['a','b'])
            True
            sage: r = iet.RauzyDiagram([0,1],[1,0])
            sage: r.alphabet() == Alphabet([0,1])
            True
        """
        if data is None:
            return self._element._alphabet
        else:
            self._element._set_alphabet(data)

    def letters(self):
        r"""
        Returns the letters used by the RauzyDiagram.

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.alphabet()
            {'a', 'b'}
            sage: r.letters()
            ['a', 'b']
            sage: r.alphabet('ABCDEF')
            sage: r.alphabet()
            {'A', 'B', 'C', 'D', 'E', 'F'}
            sage: r.letters()
            ['A', 'B']
        """
        return self._element.letters()

    def _vertex_to_permutation(self,data=None):
        r"""
        Converts the (vertex) data to a permutation.

        TESTS:

            sage: r = iet.RauzyDiagram('a b','b a')   #indirect doctest
        """
        if data is not None:
            self._set_element(data)
            return copy(self._element)

    def edge_to_matrix(self, p=None, edge_type=None):
        r"""
        Return the corresponding matrix

        INPUT:

        - ``p`` - a permutation

        - ``edge_type`` - 0 or 1 corresponding to the type of the edge

        OUTPUT:

        A matrix

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a')
            sage: d = p.rauzy_diagram()
            sage: d.edge_to_matrix(p,1)
            [1 0 1]
            [0 1 0]
            [0 0 1]
        """
        if p is None and edge_type is None:
            return identity_matrix(self._n)

        function_name = self._edge_types[edge_type][0] + '_matrix'

        if not hasattr(self._element_class,function_name):
            return identity_matrix(self._n)

        arguments = self._edge_types[edge_type][1]

        return getattr(p,function_name)(*arguments)

    def edge_to_winner(self, p=None, edge_type=None):
        r"""
        Return the corresponding winner

        TEST::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.edge_to_winner(None,None)
            []
        """
        if p is None and edge_type is None:
            return []

        function_name = self._edge_types[edge_type][0] + '_winner'

        if not hasattr(self._element_class, function_name):
            return [None]

        arguments = self._edge_types[edge_type][1]

        return [getattr(p,function_name)(*arguments)]

    def edge_to_loser(self, p=None, edge_type=None):
        r"""
        Return the corresponding loser

        TEST::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.edge_to_loser(None,None)
            []
        """
        if p is None and edge_type is None:
            return []

        function_name = self._edge_types[edge_type][0] + '_loser'

        if not hasattr(self._element_class, function_name):
            return [None]

        arguments = self._edge_types[edge_type][1]

        return [getattr(p,function_name)(*arguments)]

    def _all_npath_extension(self, g, length=0):
        r"""
        Returns an iterator over all extension of fixed length of p.

        INPUT:

        - ``p`` - a path

        - ``length`` - a non-negative integer

        TESTS:

        ::

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: for g in r._all_npath_extension(g0,0):
            ....:     print(g)
            Path of length 0 in a Rauzy diagram
            sage: for g in r._all_npath_extension(g0,1):
            ....:     print(g)
            Path of length 1 in a Rauzy diagram
            Path of length 1 in a Rauzy diagram
            sage: for g in r._all_npath_extension(g0,2):
            ....:     print(g)
            Path of length 2 in a Rauzy diagram
            Path of length 2 in a Rauzy diagram
            Path of length 2 in a Rauzy diagram
            Path of length 2 in a Rauzy diagram

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: len(list(r._all_npath_extension(g0,0)))
            1
            sage: len(list(r._all_npath_extension(g0,1)))
            1
            sage: len(list(r._all_npath_extension(g0,2)))
            2
            sage: len(list(r._all_npath_extension(g0,3)))
            3
            sage: len(list(r._all_npath_extension(g0,4)))
            5
        """
        if length == 0:
            yield g

        else:
            for i in range(len(self._edge_types)):
                if self._succ[g._end][i] is not None:
                    g._fast_append(i)
                    for h in self._all_npath_extension(g,length-1): yield h
                    g.pop()

    def _all_path_extension(self, g, length=0):
        r"""
        Returns an iterator over all path extension of p.

        TESTS:

        ::

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: for g in r._all_path_extension(g0,0):
            ....:     print(g)
            Path of length 0 in a Rauzy diagram
            sage: for g in r._all_path_extension(g0, 1):
            ....:     print(g)
            Path of length 0 in a Rauzy diagram
            Path of length 1 in a Rauzy diagram
            Path of length 1 in a Rauzy diagram

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: len(list(r._all_path_extension(g0,0)))
            1
            sage: len(list(r._all_path_extension(g0,1)))
            2
            sage: len(list(r._all_path_extension(g0,2)))
            4
            sage: len(list(r._all_path_extension(g0,3)))
            7
        """
        yield g

        if length > 0:
            for i in range(len(self._edge_types)):
                if self._succ[g._end][i] is not None:
                    g._fast_append(i)
                    for h in self._all_path_extension(g,length-1): yield h
                    g.pop()

    def __iter__(self):
        r"""
        Iterator over the permutations of the Rauzy diagram.

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r: print(p)
            a b
            b a
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: for p in r: print(p.stratum())
            H(0, 0)
            H(0, 0)
            H(0, 0)
        """
        for data in self._succ.iterkeys():
            yield self._vertex_to_permutation(data)

    def __contains__(self, element):
        r"""
        Containance test.

        TESTS::

            sage: p = iet.Permutation('a b c d','d c b a',reduced=True)
            sage: q = iet.Permutation('a b c d','d b c a',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: s = q.rauzy_diagram()
            sage: p in r
            True
            sage: p in s
            False
            sage: q in r
            False
            sage: q in s
            True
        """
        for p in self._succ.iterkeys():
            if self._vertex_to_permutation(p) == element:
                return True

        return False

    def _repr_(self):
        r"""
        Returns a representation of self

        TEST::

            sage: iet.RauzyDiagram('a b','b a')   #indirect doctest
            Rauzy diagram with 1 permutation
            sage: iet.RauzyDiagram('a b c','c b a')   #indirect doctest
            Rauzy diagram with 3 permutations
        """
        if len(self._succ) == 0:
            return "Empty Rauzy diagram"
        elif len(self._succ) == 1:
            return "Rauzy diagram with 1 permutation"
        else:
            return "Rauzy diagram with %d permutations" % (len(self._succ))

    def __getitem__(self,p):
        r"""
        Returns the neighbors of p.

        Just use the function vertex_to_permutation that must be defined
        in each child.

        INPUT:

        - ``p`` - a permutation in the Rauzy diagram

        TESTS::


            sage: p = iet.Permutation('a b c','c b a')
            sage: p0 = iet.Permutation('a b c','c a b',alphabet="abc")
            sage: p1 = iet.Permutation('a c b','c b a',alphabet="abc")
            sage: r = p.rauzy_diagram()
            sage: r[p] == [p0, p1]
            True
            sage: r[p1] == [p1, p]
            True
            sage: r[p0] == [p, p0]
            True
        """
        if not isinstance(p, self._element_class):
            raise ValueError("Your element does not have the good type")

        perm = self._permutation_to_vertex(p)
        return [self._vertex_to_permutation(x) for x in self._succ[perm]]

    def __len__(self):
        r"""
        Deprecated use cardinality.

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.cardinality()
            1
        """
        return self.cardinality()

    def cardinality(self):
        r"""
        Returns the number of permutations in this Rauzy diagram.

        OUTPUT:

        - `integer` - the number of vertices in the diagram

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.cardinality()
            1
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r.cardinality()
            3
            sage: r = iet.RauzyDiagram('a b c d','d c b a')
            sage: r.cardinality()
            7
        """
        return len(self._succ)

    def complete(self, p):
        r"""
        Completion of the Rauzy diagram.

        Add to the Rauzy diagram all permutations that are obtained by
        successive operations defined by edge_types(). The permutation must be
        of the same type and the same length as the one used for the creation.

        INPUT:

        - ``p`` - a permutation of Interval exchange transformation

        Rauzy diagram is the reunion of all permutations that could be
        obtained with successive rauzy moves. This function just use the
        functions __getitem__ and has_rauzy_move and rauzy_move which must
        be defined for child and their corresponding permutation types.

        TEST::

            sage: r = iet.RauzyDiagram('a b c','c b a')   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',left_induction=True) #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',symmetric=True)   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',lr_inversion=True)   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',tb_inversion=True)   #indirect doctest
        """
        if p.__class__ is not self._element_class:
            raise ValueError("your permutation is not of good type")

        if len(p) != self._n:
            raise ValueError("your permutation has not the good length")

        pred = self._pred
        succ = self._succ
        p = self._permutation_to_vertex(p)
        perm = self._element
        l = []

        if not p in succ:
            succ[p] = [None] * len(self._edge_types)
            pred[p] = [None] * len(self._edge_types)
            l.append(p)

        while(l != []):
            p = l.pop()
            self._set_element(p)

            for t, edge in enumerate(self._edge_types):
                if (not hasattr(perm, 'has_'+edge[0]) or
                    getattr(perm, 'has_'+edge[0])(*(edge[1]))):
                    q = getattr(perm,edge[0])(*(edge[1]))
                    q = self._permutation_to_vertex(q)
                    if not q in succ:
                        succ[q] = [None] * len(self._edge_types)
                        pred[q] = [None] * len(self._edge_types)
                        l.append(q)
                    succ[p][t] = q
                    pred[q][t] = p

    def path(self, *data):
        r"""
        Returns a path over this Rauzy diagram.

        INPUT:

        - ``initial_vertex`` - the initial vertex (starting point of the path)

        - ``data`` - a sequence of edges

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: g = r.path(p, 'top', 'bottom')
        """
        if len(data) == 0:
            raise TypeError("Must be non empty")
        elif len(data) == 1 and isinstance(data[0], self.Path):
                return copy(data[0])
        return self.Path(self,*data)

    def graph(self):
        r"""
        Returns the Rauzy diagram as a Graph object

        The graph returned is more precisely a DiGraph (directed graph) with
        loops and multiedges allowed.

        EXAMPLES::

            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r
            Rauzy diagram with 3 permutations
            sage: r.graph()
            Looped multi-digraph on 3 vertices

        """
        G = DiGraph(loops=True,multiedges=True)

        for p,neighbours in self._succ.iteritems():
            p = self._vertex_to_permutation(p)
            for i,n in enumerate(neighbours):
                if n is not None:
                    q = self._vertex_to_permutation(n)
                    G.add_edge(p,q,i)

        return G

class FlippedRauzyDiagram(RauzyDiagram):
    r"""
    Template for flipped Rauzy diagrams.

    .. WARNING::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2009-09-29): initial version
    """
    def complete(self, p, reducible=False):
        r"""
        Completion of the Rauzy diagram

        Add all successors of p for defined operations in edge_types. Could be
        used for generating non (strongly) connected Rauzy diagrams. Sometimes,
        for flipped permutations, the maximal connected graph in all
        permutations is not strongly connected. Finding such components needs to
        call most than once the .complete() method.

        INPUT:

        - ``p`` - a permutation

        - ``reducible`` - put or not reducible permutations

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a',flips='a')
            sage: d = p.rauzy_diagram()
            sage: d
            Rauzy diagram with 3 permutations
            sage: p = iet.Permutation('a b c','c b a',flips='b')
            sage: d.complete(p)
            sage: d
            Rauzy diagram with 8 permutations
            sage: p = iet.Permutation('a b c','c b a',flips='a')
            sage: d.complete(p)
            sage: d
            Rauzy diagram with 8 permutations
        """
        if p.__class__ is not self._element_class:
            raise ValueError("your permutation is not of good type")

        if len(p) != self._n:
            raise ValueError("your permutation has not the good length")

        pred = self._pred
        succ = self._succ
        p = self._permutation_to_vertex(p)
        l = []

        if not p in succ:
            succ[p] = [None] * len(self._edge_types)
            pred[p] = [None] * len(self._edge_types)
            l.append(p)

        while(l != []):
            p = l.pop()

            for t, edge_type in enumerate(self._edge_types):
                perm = self._vertex_to_permutation(p)

                if (not hasattr(perm,'has_' + edge_type[0]) or
                    getattr(perm, 'has_' + edge_type[0])(*(edge_type[1]))):
                    q = perm.rauzy_move(t)
                    q = self._permutation_to_vertex(q)
                    if reducible or perm.is_irreducible():
                        if not q in succ:
                            succ[q] = [None] * len(self._edge_types)
                            pred[q] = [None] * len(self._edge_types)
                            l.append(q)

                        succ[p][t] = q
                        pred[q][t] = p
