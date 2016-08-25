r"""
Reduced permutations

A reduced (generalized) permutation is better suited to study strata of Abelian
(or quadratic) holomorphic forms on Riemann surfaces. The Rauzy diagram is an
invariant of such a component. Corentin Boissy proved the identification of
Rauzy diagrams with connected components of stratas. But the geometry of the
diagram and the relation with the strata is not yet totally understood.

AUTHORS:

- Vincent Delecroix (2000-09-29): initial version

TESTS::

    sage: from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
    sage: ReducedPermutationIET([['a','b'],['b','a']])
    a b
    b a
    sage: ReducedPermutationIET([[1,2,3],[3,1,2]])
    1 2 3
    3 1 2
    sage: from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
    sage: ReducedPermutationLI([[1,1],[2,2,3,3,4,4]])
    1 1
    2 2 3 3 4 4
    sage: ReducedPermutationLI([['a','a','b','b','c','c'],['d','d']])
    a a b b c c
    d d
    sage: from sage.dynamics.interval_exchanges.reduced import FlippedReducedPermutationIET
    sage: FlippedReducedPermutationIET([[1,2,3],[3,2,1]],flips=[1,2])
    -1 -2  3
     3 -2 -1
    sage: FlippedReducedPermutationIET([['a','b','c'],['b','c','a']],flips='b')
     a -b  c
    -b  c  a
    sage: from sage.dynamics.interval_exchanges.reduced import FlippedReducedPermutationLI
    sage: FlippedReducedPermutationLI([[1,1],[2,2,3,3,4,4]], flips=[1,4])
    -1 -1
     2  2  3  3 -4 -4
    sage: FlippedReducedPermutationLI([['a','a','b','b'],['c','c']],flips='ac')
    -a -a  b  b
    -c -c
    sage: from sage.dynamics.interval_exchanges.reduced import ReducedRauzyDiagram
    sage: p = ReducedPermutationIET([[1,2,3],[3,2,1]])
    sage: d = ReducedRauzyDiagram(p)
"""
#*****************************************************************************
#       Copyright (C) 2008 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from sage.structure.sage_object import SageObject

from copy import copy

from sage.combinat.words.alphabet import Alphabet
from sage.rings.integer import Integer

from .template import PermutationIET, PermutationLI   # permutations
from .template import FlippedPermutationIET, FlippedPermutationLI   # flipped permutations
from .template import twin_list_iet, twin_list_li
from .template import RauzyDiagram, FlippedRauzyDiagram

from .template import interval_conversion, side_conversion

class ReducedPermutation(SageObject) :
    r"""
    Template for reduced objects.

    .. warning::

        Internal class! Do not use directly!

    INPUT:

    - ``intervals`` - a list of two list of labels

    - ``alphabet`` - (default: None) any object that can be used to initialize an Alphabet or None. In this latter case, the letter of the intervals are used to generate one.
   """
    def __init__(self,intervals=None,alphabet=None):
        r"""
        TESTS::

            sage: from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            sage: p = ReducedPermutationIET()
            sage: loads(dumps(p)) == p
            True
            sage: p = ReducedPermutationIET([['a','b'],['b','a']])
            sage: loads(dumps(p)) == p
            True
            sage: from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            sage: p = ReducedPermutationLI()
            sage: loads(dumps(p)) == p
            True
            sage: p = ReducedPermutationLI([['a','a'],['b','b']])
            sage: loads(dumps(p)) == p
            True
        """
        self._hash = None

        if intervals is None:
            self._twin = [[],[]]
            self._alphabet = alphabet

        else:
            self._init_twin(intervals)

            if alphabet is None:
                self._init_alphabet(intervals)
            else:
                alphabet = Alphabet(alphabet)
                if alphabet.cardinality() < len(self):
                    raise TypeError("The alphabet is too short")
                self._alphabet = alphabet

    def __len__(self):
        r"""
        TESTS::

            sage: p = iet.Permutation('a b','b a',reduced=True)
            sage: len(p)
            2
            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: len(p)
            3
            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: len(p)
            3
        """
        return (len(self._twin[0]) + len(self._twin[1])) / 2

    def length_top(self):
        r"""
        Returns the number of intervals in the top segment.

        OUTPUT:

        integer -- the length of the top segment

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.length_top()
            3
            sage: p = iet.GeneralizedPermutation('a a b','c d c b d')
            sage: p.length_top()
            3
            sage: p = iet.GeneralizedPermutation('a b c b d c d', 'e a e')
            sage: p.length_top()
            7
        """
        return len(self._twin[0])

    def length_bottom(self):
        r"""
        Returns the number of intervals in the bottom segment.

        OUTPUT:

        integer -- the length of the bottom segment

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.length_bottom()
            3
            sage: p = iet.GeneralizedPermutation('a a b','c d c b d')
            sage: p.length_bottom()
            5
        """
        return len(self._twin[1])

    def length(self, interval=None):
        r"""
        Returns the 2-uple of lengths.

        p.length() is identical to (p.length_top(), p.length_bottom())
        If an interval is specified, it returns the length of the specified
        interval.

        INPUT:

        - ``interval`` - None, 'top' (or 't' or 0) or 'bottom' (or 'b' or 1)

        OUTPUT:

        integer or 2-uple of integers -- the corresponding lengths

        EXAMPLES::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.length()
            (3, 3)
            sage: p = iet.GeneralizedPermutation('a a b','c d c b d')
            sage: p.length()
            (3, 5)
        """
        if interval is None :
            return len(self._twin[0]),len(self._twin[1])
        else :
            interval = interval_conversion(interval)
            return len(self._twin[interval])

    def __getitem__(self,i):
        r"""
        TESTS::

            sage: p = iet.Permutation('a b', 'b a', reduced=True)
            sage: p[0]
            ['a', 'b']
            sage: p[1]
            ['b', 'a']
            sage: p.alphabet([0,1])
            sage: p[0]
            [0, 1]
            sage: p[1]
            [1, 0]
        """
        return self.list()[i]

    def __copy__(self) :
        r"""
        Returns a copy of self.

        EXAMPLES::

            sage: p = iet.Permutation('a b c', 'c b a', reduced=True)
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
            sage: p = iet.GeneralizedPermutation('a b b','c c a', reduced=True)
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
        """
        q = self.__class__()

        q._twin = [self._twin[0][:], self._twin[1][:]]
        q._alphabet = self._alphabet
        q._repr_type = self._repr_type
        q._repr_options = self._repr_options

        return q

    def erase_letter(self, letter):
        r"""
        Erases a letter.

        INPUT:

        - ``letter`` - a letter which is a label of an interval of self

        EXAMPLES:

        ::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.erase_letter('a')
            b c
            c b

        ::

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.erase_letter('a')
            b b
            c c
        """
        l = self.list()

        del l[0][l[0].index(letter)]
        del l[1][l[1].index(letter)]

        return self.__class__(l)

    def right_rauzy_move(self, winner):
        r"""
        Performs a Rauzy move on the right.

        EXAMPLES:

        ::

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: p.right_rauzy_move(0)
            a b c
            c a b
            sage: p.right_rauzy_move(1)
            a b c
            b c a

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: p.right_rauzy_move(0)
            a b b
            c c a
        """
        winner = interval_conversion(winner)

        result = copy(self)

        loser_to = result._get_loser_to(winner)
        # beware here, loser_to can contain 2 or 3 items
        # (depending on the presence of flip)

        result._twin_rauzy_move(winner, loser_to)

        return result

    def left_rauzy_move(self, winner):
        r"""
        Performs a Rauzy move on the left.

        EXAMPLES:

        ::

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: p.left_rauzy_move(0)
            a b c
            b c a
            sage: p.right_rauzy_move(1)
            a b c
            b c a

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: p.left_rauzy_move(0)
            a a b
            b c c
        """
        winner = interval_conversion(winner)

        result = copy(self)
        result._reversed()
        result = result.right_rauzy_move(winner)
        result._reversed()
        return result

_RP = ReducedPermutation

def alphabetized_atwin(twin, alphabet):
    """
    Alphabetization of a twin of iet.

    TESTS::

        sage: from sage.dynamics.interval_exchanges.reduced import alphabetized_atwin

    ::

        sage: twin = [[0,1],[0,1]]
        sage: alphabet = Alphabet("ab")
        sage: alphabetized_atwin(twin, alphabet)
        [['a', 'b'], ['a', 'b']]

    ::

        sage: twin = [[1,0],[1,0]]
        sage: alphabet = Alphabet([0,1])
        sage: alphabetized_atwin(twin, alphabet)
        [[0, 1], [1, 0]]

    ::

        sage: twin = [[1,2,3,0],[3,0,1,2]]
        sage: alphabet = Alphabet("abcd")
        sage: alphabetized_atwin(twin,alphabet)
        [['a', 'b', 'c', 'd'], ['d', 'a', 'b', 'c']]
    """
    l = [[],[]]

    l[0] = [alphabet.unrank(x) for x in range(len(twin[0]))]
    l[1] = [alphabet.unrank(x) for x in twin[1]]

    return l

def ReducedPermutationsIET_iterator(
    nintervals=None,
    irreducible=True,
    alphabet=None):
    r"""
    Returns an iterator over reduced permutations

    INPUT:

    - ``nintervals`` - integer or None

    - ``irreducible`` - boolean

    - ``alphabet`` - something that should be converted to an alphabet
      of at least nintervals letters

    TESTS::

        sage: for p in iet.Permutations_iterator(3,reduced=True,alphabet="abc"):
        ....:     print(p)  #indirect doctest
        a b c
        b c a
        a b c
        c a b
        a b c
        c b a
    """
    from itertools import imap,ifilter
    from sage.combinat.permutation import Permutations

    if irreducible is False:
        if nintervals is None:
            raise ValueError("please choose a number of intervals")

        nintervals = Integer(nintervals)

        if not(nintervals > 0):
            raise ValueError('number of intervals must be positive')

        a0 = range(1,nintervals+1)
        f = lambda x: ReducedPermutationIET([a0,list(x)],
            alphabet=alphabet)
        return imap(f, Permutations(nintervals))
    else:
        return ifilter(lambda x: x.is_irreducible(),
        ReducedPermutationsIET_iterator(nintervals,False,alphabet))

class ReducedPermutationIET(ReducedPermutation, PermutationIET):
    """
    Reduced permutation from iet

    Permutation from iet without numerotation of intervals. For initialization,
    you should use GeneralizedPermutation which is the class factory for all
    permutation types.

    EXAMPLES:

    Equality testing (no equality of letters but just of ordering)::

        sage: p = iet.Permutation('a b c', 'c b a', reduced = True)
        sage: q = iet.Permutation('p q r', 'r q p', reduced = True)
        sage: p == q
        True

    Reducibility testing::

        sage: p = iet.Permutation('a b c', 'c b a', reduced = True)
        sage: p.is_irreducible()
        True

    ::

        sage: q = iet.Permutation('a b c d', 'b a d c', reduced = True)
        sage: q.is_irreducible()
        False


    Rauzy movability and Rauzy move::

        sage: p = iet.Permutation('a b c', 'c b a', reduced = True)
        sage: p.has_rauzy_move(1)
        True
        sage: p.rauzy_move(1)
        a b c
        b c a

    Rauzy diagrams::

        sage: p = iet.Permutation('a b c d', 'd a b c')
        sage: p_red = iet.Permutation('a b c d', 'd a b c', reduced = True)
        sage: d = p.rauzy_diagram()
        sage: d_red = p_red.rauzy_diagram()
        sage: p.rauzy_move(0) in d
        True
        sage: d.cardinality(), d_red.cardinality()
        (12, 6)
    """
    def _init_twin(self, a):
        r"""
        Initializes the _twin attribute

        TESTS::

            sage: p = iet.Permutation('a b','b a',reduced=True)   #indirect doctest
            sage: p._twin
            [[1, 0], [1, 0]]
        """
        self._twin = twin_list_iet(a)

    def list(self):
        r"""
        Returns a list of two list that represents the permutation.

        EXAMPLES::

            sage: p = iet.GeneralizedPermutation('a b','b a',reduced=True)
            sage: p.list() == [p[0], p[1]]
            True
            sage: p.list() == [['a', 'b'], ['b', 'a']]
            True

        ::

            sage: p = iet.GeneralizedPermutation('a b c', 'b c a',reduced=True)
            sage: iet.GeneralizedPermutation(p.list(),reduced=True) == p
            True
        """
        a0 = [self._alphabet.unrank(_) for _ in range(0,len(self))]
        a1 = [self._alphabet.unrank(_) for _ in self._twin[1]]
        return [a0,a1]

    def is_identity(self):
        r"""
        Returns True if self is the identity.

        EXAMPLES::

            sage: iet.Permutation("a b","a b",reduced=True).is_identity()
            True
            sage: iet.Permutation("a b","b a",reduced=True).is_identity()
            False
        """
        for i in range(len(self)):
            if self._twin[0][i] != i:
                return False
        return True

    def __hash__(self):
        r"""
        Returns a hash value (does not depends of the alphabet).

        TESTS::

            sage: p = iet.Permutation([1,2],[1,2], reduced=True)
            sage: q = iet.Permutation([1,2],[2,1], reduced=True)
            sage: r = iet.Permutation([2,1],[1,2], reduced=True)
            sage: hash(p) == hash(q)
            False
            sage: hash(q) == hash(r)
            True
        """
        if self._hash is None:
            self._hash = hash(tuple(self._twin[0]))
        return self._hash

    def __cmp__(self, other):
        r"""
        Defines a natural lexicographic order.

        TESTS::

            sage: p = iet.GeneralizedPermutation('a b','a b',reduced=True)
            sage: q = copy(p)
            sage: q.alphabet([0,1])
            sage: p == q
            True
            sage: p0 = iet.GeneralizedPermutation('a b', 'a b', reduced=True)
            sage: p1 = iet.GeneralizedPermutation('a b', 'b a', reduced=True)
            sage: p0 < p1 and p1 > p0
            True
            sage: q0 = iet.GeneralizedPermutation('a b c','a b c',reduced=True)
            sage: q1 = iet.GeneralizedPermutation('a b c','a c b',reduced=True)
            sage: q2 = iet.GeneralizedPermutation('a b c','b a c',reduced=True)
            sage: q3 = iet.GeneralizedPermutation('a b c','b c a',reduced=True)
            sage: q4 = iet.GeneralizedPermutation('a b c','c a b',reduced=True)
            sage: q5 = iet.GeneralizedPermutation('a b c','c b a',reduced=True)
            sage: p0 < q0 and q0 > p0 and p1 < q0 and q0 > p1
            True
            sage: q0 < q1 and q1 > q0
            True
            sage: q1 < q2 and q2 > q1
            True
            sage: q2 < q3 and q3 > q2
            True
            sage: q3 < q4 and q4 > q3
            True
            sage: q4 < q5 and q5 > q4
            True
        """
        if type(self) is not type(other):
            raise ValueError("Permutations must be of the same type")

        if len(self) > len(other):
            return 1
        elif len(self) < len(other):
            return -1

        n = len(self)
        j = 0
        while (j < n and self._twin[1][j] == other._twin[1][j]):
            j += 1

        if j != n:
            if self._twin[1][j] > other._twin[1][j]: return 1
            else: return -1

        return 0

    def _reversed(self):
        r"""
        Reverses the permutation.

        EXAMPLES:

        ::

            sage: p = iet.Permutation('a b c','c a b',reduced=True)
            sage: p
            a b c
            c a b
            sage: p._reversed()
            sage: p
            a b c
            b c a

        ::

            sage: p = iet.Permutation('a b c d','d a b c',reduced=True)
            sage: p
            a b c d
            d a b c
            sage: p._reversed()
            sage: p
            a b c d
            b c d a
        """
        tmp = [self._twin[0][:], self._twin[1][:]]

        n = self.length_top()
        for i in (0, 1):
            for j in range(n):
                tmp[i][n - 1 - j] = n - 1 - self._twin[i][j]

        self._twin = tmp

    def _inversed(self):
        r"""
        Inverses the permutation.

        EXAMPLES:

        ::

            sage: p = iet.Permutation('a b c','c a b',reduced=True)
            sage: p
            a b c
            c a b
            sage: p._inversed()
            sage: p
            a b c
            b c a

        ::

            sage: p = iet.Permutation('a b c d','d a b c',reduced=True)
            sage: p
            a b c d
            d a b c
            sage: p._inversed()
            sage: p
            a b c d
            b c d a
        """
        self._twin = [self._twin[1], self._twin[0]]

    def has_rauzy_move(self, winner, side='right'):
        r"""
        Tests if the permutation is rauzy_movable on the left.

        EXAMPLES:

        ::

            sage: p = iet.Permutation('a b c','a c b',reduced=True)
            sage: p.has_rauzy_move(0,'right')
            True
            sage: p.has_rauzy_move(0,'left')
            False
            sage: p.has_rauzy_move(1,'right')
            True
            sage: p.has_rauzy_move(1,'left')
            False

        ::

            sage: p = iet.Permutation('a b c d','c a b d',reduced=True)
            sage: p.has_rauzy_move(0,'right')
            False
            sage: p.has_rauzy_move(0,'left')
            True
            sage: p.has_rauzy_move(1,'right')
            False
            sage: p.has_rauzy_move(1,'left')
            True
        """
        side = side_conversion(side)
        winner = interval_conversion(winner)

        return self._twin[winner][side] % len(self) != side % len(self)

    def rauzy_move_relabel(self, winner, side='right'):
        r"""
        Returns the relabelization obtained from this move.

        EXAMPLE::

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: q = p.reduced()
            sage: p_t = p.rauzy_move('t')
            sage: q_t = q.rauzy_move('t')
            sage: s_t = q.rauzy_move_relabel('t')
            sage: s_t
            WordMorphism: a->a, b->b, c->c, d->d
            sage: map(s_t, p_t[0]) == map(Word, q_t[0])
            True
            sage: map(s_t, p_t[1]) == map(Word, q_t[1])
            True
            sage: p_b = p.rauzy_move('b')
            sage: q_b = q.rauzy_move('b')
            sage: s_b = q.rauzy_move_relabel('b')
            sage: s_b
            WordMorphism: a->a, b->d, c->b, d->c
            sage: map(s_b, q_b[0]) == map(Word, p_b[0])
            True
            sage: map(s_b, q_b[1]) == map(Word, p_b[1])
            True
        """
        from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
        from sage.combinat.words.morphism import WordMorphism

        winner = interval_conversion(winner)
        side = side_conversion(side)

        p = LabelledPermutationIET(self.list())

        l0_q = p.rauzy_move(winner, side).list()[0]

        d = dict([(self._alphabet[i],l0_q[i]) for i in range(len(self))])

        return WordMorphism(d)

    def _get_loser_to(self, winner) :
        r"""
        This function return the position of the future loser position.

        TESTS:

        ::

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: p._get_loser_to(0)
            (1, 1)
            sage: p._get_loser_to(1)
            (0, 1)

        ::

            sage: p = iet.Permutation('a b c','c a b',reduced=True)
            sage: p._get_loser_to(0)
            (1, 1)
            sage: p._get_loser_to(1)
            (0, 2)

        ::

            sage: p = iet.Permutation('a b c','b c a',reduced=True)
            sage: p._get_loser_to(0)
            (1, 2)
            sage: p._get_loser_to(1)
            (0, 1)
        """
        return (1-winner, self._twin[winner][-1]+1)

    def _twin_rauzy_move(self, winner_interval, loser_to) :
        r"""
        Do a Rauzy move for this choice of winner.

        TESTS::

            sage: p = iet.Permutation('a b','b a',reduced=True)
            sage: p.rauzy_move(0) == p   #indirect doctest
            True
            sage: p.rauzy_move(1) == p   #indirect doctest
            True
        """
        loser_interval = 1 - winner_interval

        loser_twin_interval = winner_interval
        loser_twin_position = self._twin[loser_interval][-1]

        loser_interval_to, loser_position_to = loser_to

        # move the loser
        del self._twin[loser_interval][-1]
        self._twin[loser_interval_to].insert(loser_position_to, loser_twin_position)
        self._twin[loser_twin_interval][loser_twin_position] = loser_position_to

        # increment the twins in the winner interval
        for j in range(loser_position_to + 1, self.length(loser_interval_to)) :
            self._twin[winner_interval][self._twin[loser_interval_to][j]] += 1

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        OUTPUT:

        A Rauzy diagram

        EXAMPLES:

        ::

            sage: p = iet.Permutation('a b c d', 'd a b c',reduced=True)
            sage: d = p.rauzy_diagram()
            sage: p.rauzy_move(0) in d
            True
            sage: p.rauzy_move(1) in d
            True

        For more information, try help RauzyDiagram
        """
        return ReducedRauzyDiagram(self, **kargs)

def alphabetized_qtwin(twin, alphabet):
    """
    Alphabetization of a qtwin.

    TESTS::

        sage: from sage.dynamics.interval_exchanges.reduced import alphabetized_qtwin

    ::

        sage: twin = [[(1,0),(1,1)],[(0,0),(0,1)]]
        sage: alphabet = Alphabet("ab")
        sage: alphabetized_qtwin(twin,alphabet)
        [['a', 'b'], ['a', 'b']]

    ::

        sage: twin = [[(1,1), (1,0)],[(0,1), (0,0)]]
        sage: alphabet=Alphabet("AB")
        sage: alphabetized_qtwin(twin,alphabet)
        [['A', 'B'], ['B', 'A']]
        sage: alphabet=Alphabet("BA")
        sage: alphabetized_qtwin(twin,alphabet)
        [['B', 'A'], ['A', 'B']]

    ::

        sage: twin = [[(0,1),(0,0)],[(1,1),(1,0)]]
        sage: alphabet=Alphabet("ab")
        sage: alphabetized_qtwin(twin,alphabet)
        [['a', 'a'], ['b', 'b']]

    ::

        sage: twin = [[(0,2),(1,1),(0,0)],[(1,2),(0,1),(1,0)]]
        sage: alphabet=Alphabet("abc")
        sage: alphabetized_qtwin(twin,alphabet)
        [['a', 'b', 'a'], ['c', 'b', 'c']]
    """
    i_a = 0
    l = [[False]*len(twin[0]), [False]*len(twin[1])]
    # False means empty here
    for i in range(2) :
        for j in range(len(l[i])) :
            if l[i][j] is False :
                l[i][j] = alphabet[i_a]
                l[twin[i][j][0]][twin[i][j][1]] = alphabet[i_a]
                i_a += 1
    return l

class ReducedPermutationLI(ReducedPermutation, PermutationLI):
    r"""
    Reduced quadratic (or generalized) permutation.

    EXAMPLES:

    Reducibility testing::

        sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
        sage: p.is_irreducible()
        True

    ::

        sage: p = iet.GeneralizedPermutation('a b c a', 'b d d c', reduced = True)
        sage: p.is_irreducible()
        False
        sage: test, decomposition = p.is_irreducible(return_decomposition = True)
        sage: test
        False
        sage: decomposition
        (['a'], ['c', 'a'], [], ['c'])

    Rauzy movavability and Rauzy move::

        sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
        sage: p.has_rauzy_move(0)
        True
        sage: p.rauzy_move(0)
        a a b b
        c c
        sage: p.rauzy_move(0).has_rauzy_move(0)
        False
        sage: p.rauzy_move(1)
        a b b
        c c a

    Rauzy diagrams::

        sage: p_red = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
        sage: d_red = p_red.rauzy_diagram()
        sage: d_red.cardinality()
        4
    """
    def _init_twin(self, a):
        r"""
        Initializes the _twin attribute

        TESTS::

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)   #indirect doctest
            sage: p._twin
            [[(0, 1), (0, 0)], [(1, 1), (1, 0)]]
        """
        self._twin = twin_list_li(a)

    def list(self) :
        r"""
        The permutations as a list of two lists.

        EXAMPLES::

            sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
            sage: list(p)
            [['a', 'b', 'b'], ['c', 'c', 'a']]
        """
        return alphabetized_qtwin(self._twin, self.alphabet())

    def __eq__(self, other) :
        r"""
        Tests equality.

        Two reduced permutations are equal if they have the same order of
        apparition of intervals. Non necessarily the same alphabet.

        TESTS::

            sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
            sage: q = iet.GeneralizedPermutation('b a a', 'c c b', reduced = True)
            sage: r = iet.GeneralizedPermutation('t s s', 'w w t', reduced = True)
            sage: p == q
            True
            sage: p == r
            True
        """
        return type(self) is type(other) and self._twin == other._twin

    def __ne__(self, other) :
        """
        Tests difference.

        TESTS::

            sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
            sage: q = iet.GeneralizedPermutation('b b a', 'c c a', reduced = True)
            sage: r = iet.GeneralizedPermutation('i j j', 'k k i', reduced = True)
            sage: p != q
            True
            sage: p != r
            False
        """
        return type(self) is not type(other) or (self._twin != other._twin)

    def _get_loser_to(self, winner) :
        r"""
        This function return the position of the future loser position.

        TESTS::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: p._get_loser_to(0)
            (1, 1)
            sage: p._get_loser_to(1)
            (1, 1)
        """
        loser = 1 - winner

        if self._twin[winner][-1][0] == loser:
            return (loser, self._twin[winner][-1][1] + 1)
        else:
            return (winner, self._twin[winner][-1][1])

    def _twin_rauzy_move(self, winner_interval, loser_to) :
        r"""
        Rauzy move on the twin data

        TESTS:

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: p.rauzy_move(0)   #indirect doctest
            a a b
            b c c
            sage: p.rauzy_move(1)   #indirect doctest
            a a
            b b c c

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: p.rauzy_move(0)   #indirect doctest
            a a b
            b c c
            sage: p.rauzy_move(1)   #indirect doctest
            a a
            b b c c
        """
        loser_interval = 1 - winner_interval

        loser_interval_to, loser_position_to = loser_to
        loser_twin_interval, loser_twin_position = self._twin[loser_interval][-1]

        # increment the twins in the winner interval
        interval = [(self._twin[loser_interval_to][j], j) for j in range(loser_position_to, len(self._twin[loser_interval_to]))]
        for (i,j),k in interval : self._twin[i][j] = (loser_interval_to, k+1)

        # prepare the loser new position in its twin
        self._twin[loser_twin_interval][loser_twin_position] = loser_to

        # move the loser
        loser_twin = self._twin[loser_interval][-1]
        self._twin[loser_interval_to].insert(loser_position_to, loser_twin)
        del self._twin[loser_interval][-1]

    def _reversed(self):
        r"""
        Reverses the permutation.

        EXAMPLES:

        ::

            sage: p = iet.GeneralizedPermutation('a b b','c c a',reduced=True)
            sage: p._reversed()
            sage: p
            a a b
            b c c

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: p._reversed()
            sage: p
            a a
            b b c c
        """
        tmp = [self._twin[0][:], self._twin[1][:]]

        n = self.length()

        for i in (0,1):
            for j in range(n[i]):
                interval, position = self._twin[i][j]
                tmp[i][n[i] - 1 - j] = (
                    interval,
                    n[interval] - 1 - position)

        self._twin = tmp

    def _inversed(self):
        r"""
        Inverses the permutation.

        EXAMPLES:

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)
            sage: p._inversed()
            sage: p
            a a
            b b

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: p._inversed()
            sage: p
            a b b
            c c a

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: p._inversed()
            sage: p
            a a b b
            c c
        """
        self._twin = [self._twin[1], self._twin[0]]

        for interval in (0,1):
            for j in xrange(self.length(interval)):
                self._twin[interval][j] = (1-self._twin[interval][j][0],
                    self._twin[interval][j][1])

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        The Rauzy diagram of a permutation corresponds to all permutations
        that we could obtain from this one by Rauzy move. The set obtained
        is a labelled Graph. The label of vertices being 0 or 1 depending
        on the type.

        OUTPUT:

        Rauzy diagram -- the graph of permutations obtained by rauzy induction

        EXAMPLES::

            sage: p = iet.Permutation('a b c d', 'd a b c')
            sage: d = p.rauzy_diagram()
        """
        return ReducedRauzyDiagram(self, **kargs)

def labelize_flip(couple):
    r"""
    Returns a string from a 2-uple couple of the form (name, flip).

    TESTS::

        sage: from sage.dynamics.interval_exchanges.reduced import labelize_flip
        sage: labelize_flip((4,1))
        ' 4'
        sage: labelize_flip(('a',-1))
        '-a'
    """
    if couple[1] == -1: return '-' + str(couple[0])
    return ' ' + str(couple[0])

class FlippedReducedPermutation(ReducedPermutation):
    r"""
    Flipped Reduced Permutation.

    .. warning::

        Internal class! Do not use directly!

    INPUT:

    - ``intervals`` - a list of two lists

    - ``flips`` - the flipped letters

    - ``alphabet`` - an alphabet
    """
    def __init__(self, intervals=None, flips=None, alphabet=None):
        r"""
        TESTS::

            sage: p = iet.Permutation('a b','b a',reduced=True,flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.Permutation('a b','b a',reduced=True,flips='b')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.Permutation('a b','b a',reduced=True,flips='ab')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='b')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='ab')
            sage: p == loads(dumps(p))
            True
        """
        self._hash = None

        if intervals is None:
            self._twin = [[],[]]
            self._flips = [[],[]]
            self._alphabet = None

        else:
            if flips is None: flips = []

            if alphabet is None : self._init_alphabet(intervals)
            else : self._alphabet = Alphabet(alphabet)

            self._init_twin(intervals)
            self._init_flips(intervals, flips)

            self._hash = None

    def __copy__(self):
        r"""
        Returns a copy of self.

        TESTS::

            sage: p = iet.GeneralizedPermutation('a a b', 'c c b',reduced=True, flips='a')
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
        """
        p = self.__class__()

        p._twin = [self._twin[0][:], self._twin[1][:]]
        p._flips = [self._flips[0][:], self._flips[1][:]]
        p._alphabet = self._alphabet
        p._repr_type = self._repr_type
        p._repr_options = self._repr_options

        return p

    def right_rauzy_move(self, winner):
        r"""
        Performs a Rauzy move on the right.

        EXAMPLE::

            sage: p = iet.Permutation('a b c','c b a',reduced=True,flips='c')
            sage: p.right_rauzy_move('top')
            -a  b -c
            -a -c  b
        """
        winner = interval_conversion(winner)

        result = copy(self)

        loser_to = result._get_loser_to(winner)

        result._flip_rauzy_move(winner, loser_to)
        result._twin_rauzy_move(winner, loser_to)

        return result

    def _reversed(self):
        r"""
        Reverses the permutation

        TESTS:

        ::

            sage: p = iet.Permutation('a b c d','c d b a',reduced=True,flips='a')
            sage: p
            -a  b  c  d
             c  d  b -a
            sage: p._reversed()
            sage: p
             a  b  c -d
            -d  c  a  b
            sage: p._reversed()
            sage: p
            -a  b  c  d
             c  d  b -a

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True,flips='a')
            sage: p
            -a -a  b
             b  c  c
            sage: p._reversed()
            sage: p
             a -b -b
             c  c  a
            sage: p._reversed()
            sage: p
            -a -a  b
             b  c  c
        """
        super(FlippedReducedPermutation, self)._reversed()
        self._flips[0].reverse()
        self._flips[1].reverse()

    def _inversed(self):
        r"""
        Inverses the permutation.

        TESTS:

        ::

            sage: p = iet.Permutation('a b c d','c d b a',flips='a')
            sage: p
            -a  b  c  d
             c  d  b -a
            sage: p._inversed()
            sage: p
             c  d  b -a
            -a  b  c  d
            sage: p._inversed()
            sage: p
            -a  b  c  d
             c  d  b -a

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True,flips='a')
            sage: p
            -a -a  b
             b  c  c
            sage: p._inversed()
            sage: p
             a  b  b
            -c -c  a
            sage: p._inversed()
            sage: p
            -a -a  b
             b  c  c
        """
        super(FlippedReducedPermutation, self)._inversed()
        self._flips.reverse()

class FlippedReducedPermutationIET(
    FlippedReducedPermutation,
    FlippedPermutationIET,
    ReducedPermutationIET):
    r"""
    Flipped Reduced Permutation from iet

    EXAMPLES

    ::

        sage: p = iet.Permutation('a b c', 'c b a', flips=['a'], reduced=True)
        sage: p.rauzy_move(1)
        -a -b  c
        -a  c -b

    TESTS::

        sage: p = iet.Permutation('a b','b a',flips=['a'])
        sage: p == loads(dumps(p))
        True
    """
    def __cmp__(self, other):
        r"""
        Defines a natural lexicographic order.

        TESTS::

            sage: p = iet.Permutation('a b','a b',reduced=True,flips='a')
            sage: q = copy(p)
            sage: q.alphabet([0,1])
            sage: p == q
            True
            sage: l0 = ['a b','a b']
            sage: l1 = ['a b','b a']
            sage: l2 = ['b a', 'a b']
            sage: p0 = iet.Permutation(l0, reduced=True, flips='ab')
            sage: p1 = iet.Permutation(l1, reduced=True, flips='a')
            sage: p2 = iet.Permutation(l2, reduced=True, flips='b')
            sage: p3 = iet.Permutation(l1, reduced=True, flips='ab')
            sage: p4 = iet.Permutation(l2 ,reduced=True,flips='ab')
            sage: p0 == p1 or p0 == p2 or p0 == p3 or p0 == p4
            False
            sage: p0 != p1 and p0 != p2 and p0 != p3 and p0 != p4
            True
            sage: p1 == p2 and p3 == p4
            True
            sage: p1 != p2 or p3 != p4
            False
            sage: p1 == p3 or p1 == p4 or p2 == p3 or p2 == p4
            False
            sage: p1 != p3 and p1 != p4 and p2 != p3 and p2 != p4
            True
            sage: p1 = iet.Permutation(l1,reduced=True, flips='a')
            sage: p2 = iet.Permutation(l1,reduced=True, flips='b')
            sage: p3 = iet.Permutation(l1,reduced=True, flips='ab')
            sage: p2 > p3 and p3 < p2
            True
            sage: p1 > p2 and p2 < p1
            True
            sage: p1 > p3 and p3 < p1
            True
            sage: q1 = iet.Permutation(l0, reduced=True, flips='a')
            sage: q2 = iet.Permutation(l0, reduced=True, flips='b')
            sage: q3 = iet.Permutation(l0, reduced=True, flips='ab')
            sage: q2 > q1 and q2 > q3 and q1 < q2 and q3 < q2
            True
            sage: q1 > q3
            True
            sage: q3 < q1
            True
            sage: r = iet.Permutation('a b c','a b c', reduced=True, flips='a')
            sage: r > p1 and r > p2 and r > p3
            True
            sage: p1 < r and p2 < r and p3 < r
            True
        """
        if type(self) is not type(other):
            return -1

        if len(self) > len(other):
            return 1
        elif len(self) < len(other):
            return -1

        n = len(self)
        j = 0
        while (j < n and
            self._twin[1][j] == other._twin[1][j] and
            self._flips[1][j] == other._flips[1][j]):
            j += 1

        if j != n:
            if self._twin[1][j] > other._twin[1][j]: return 1
            elif self._twin[1][j] < other._twin[1][j]: return -1
            else: return self._flips[1][j]

        return 0

    def list(self, flips=False):
        r"""
        Returns a list representation of self.

        INPUT:

        - ``flips`` - boolean (default: False) if True the output contains
           2-uple of (label, flip)

        EXAMPLES::

            sage: p = iet.Permutation('a b','b a',reduced=True,flips='b')
            sage: p.list(flips=True)
            [[('a', 1), ('b', -1)], [('b', -1), ('a', 1)]]
            sage: p.list(flips=False)
            [['a', 'b'], ['b', 'a']]
            sage: p.alphabet([0,1])
            sage: p.list(flips=True)
            [[(0, 1), (1, -1)], [(1, -1), (0, 1)]]
            sage: p.list(flips=False)
            [[0, 1], [1, 0]]

        One can recover the initial permutation from this list::

            sage: p = iet.Permutation('a b','b a',reduced=True,flips='a')
            sage: iet.Permutation(p.list(), flips=p.flips(), reduced=True) == p
            True
        """
        if flips:
            a0 = zip([self.alphabet().unrank(_) for _ in range(0,len(self))], self._flips[0])
            a1 = zip([self.alphabet().unrank(_) for _ in self._twin[1]], self._flips[1])

        else:
            a0 = [self.alphabet().unrank(_) for _ in range(0,len(self))]
            a1 = [self.alphabet().unrank(_) for _ in self._twin[1]]

        return [a0,a1]

    def _get_loser_to(self, winner) :
        r"""
        This function return the position of the future loser position.

        TESTS:

        ::

            sage: p = iet.GeneralizedPermutation('a b b','c c a', reduced=True, flips='ab')
            sage: p._get_loser_to(0)
            (0, 2)
            sage: p._get_loser_to(1)
            (0, 0)

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True,flips='a')
            sage: p._get_loser_to(0)
            (0, 1)
        """
        if self._flips[winner][-1] == 1 :
            return (1-winner, self._twin[winner][-1]+1)
        else :
            return (1-winner, self._twin[winner][-1])

    def _flip_rauzy_move(self, winner, loser_to):
        r"""
        Performs a Rauzy move on flips.

        TESTS::

            sage: p = iet.GeneralizedPermutation('a b b','c c a', reduced=True, flips='ab')
            sage: p
            -a -b -b
             c  c -a
            sage: p.rauzy_move(1)   #indirect doctest
             a -b  a
             c  c -b
            sage: p.rauzy_move(0)   #indirect doctest
             a  -b  a -b
             c  c
        """
        loser = 1 - winner

        loser_twin_interval, loser_twin_position = winner, self._twin[loser][-1]
        loser_interval_to, loser_position_to = loser_to

        flip = self._flips[winner][-1] * self._flips[loser][-1]

        self._flips[loser_twin_interval][loser_twin_position] = flip

        del self._flips[loser][-1]
        self._flips[loser_interval_to].insert(loser_position_to, flip)

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        EXAMPLES::

            sage: p = iet.Permutation('a b','b a',reduced=True,flips='a')
            sage: r = p.rauzy_diagram()
            sage: p in r
            True
        """
        return FlippedReducedRauzyDiagram(self, **kargs)

class FlippedReducedPermutationLI(
    FlippedReducedPermutation,
    FlippedPermutationLI,
    ReducedPermutationLI):
    r"""
    Flipped Reduced Permutation from li

    EXAMPLES:

    Creation using the GeneralizedPermutation function::

        sage: p = iet.GeneralizedPermutation('a a b', 'b c c', reduced=True, flips='a')
    """
    def _get_loser_to(self, winner) :
        r"""
        This function return the position of the future loser position.

        TESTS:

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True,flips='b')
            sage: p._get_loser_to(0)
            (1, 0)
            sage: p._get_loser_to(1)
            (1, 1)

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True,flips='c')
            sage: p._get_loser_to(0)
            (1, 1)
            sage: p._get_loser_to(1)
            (1, 2)
        """
        loser = 1 - winner

        if self._twin[winner][-1][0] == loser :
            if self._flips[winner][-1] == 1 :
                return (loser, self._twin[winner][-1][1] + 1)
            else :
                return (loser, self._twin[winner][-1][1])
        else :
            if self._flips[winner][-1] == 1 :
                return (winner, self._twin[winner][-1][1])
            else :
                return (winner, self._twin[winner][-1][1] + 1)

    def _flip_rauzy_move(self, winner, loser_to):
        r"""
        Performs a Rauzy move on flips

        TESTS::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',flips='b',reduced=True)
            sage: p
             a  a -b
            -b  c  c
            sage: p.rauzy_move('top')   #indirect doctest
             a  a -b
            -c -b -c
        """
        loser = 1 - winner

        loser_twin_interval, loser_twin_position = self._twin[loser][-1]
        loser_interval_to, loser_position_to = loser_to

        flip = self._flips[winner][-1] * self._flips[loser][-1]

        self._flips[loser_twin_interval][loser_twin_position] = flip

        del self._flips[loser][-1]
        self._flips[loser_interval_to].insert(loser_position_to, flip)

    def list(self, flips=False):
        r"""
        Returns a list representation of self.

        INPUT:

        - ``flips`` - boolean (default: False) return the list with flips

        EXAMPLES:

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='a')
            sage: p.list(flips=True)
            [[('a', -1), ('a', -1)], [('b', 1), ('b', 1)]]
            sage: p.list(flips=False)
            [['a', 'a'], ['b', 'b']]

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True,flips='abc')
            sage: p.list(flips=True)
            [[('a', -1), ('a', -1), ('b', -1)], [('b', -1), ('c', -1), ('c', -1)]]
            sage: p.list(flips=False)
            [['a', 'a', 'b'], ['b', 'c', 'c']]

        one can rebuild the permutation from the list::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',flips='a',reduced=True)
            sage: iet.GeneralizedPermutation(p.list(),flips=p.flips(),reduced=True) == p
            True
        """
        i_a = 0
        l = [[False]*len(self._twin[0]), [False]*len(self._twin[1])]

        if flips:
            for i in range(2):  # False means empty here
                for j in range(len(l[i])):
                    if l[i][j] is False:
                        l[i][j] = (self._alphabet.unrank(i_a), self._flips[i][j])
                        l[self._twin[i][j][0]][self._twin[i][j][1]] = l[i][j]
                        i_a += 1

        else:
            for i in range(2):  # False means empty here
                for j in range(len(l[i])):
                    if l[i][j] is False:
                        l[i][j] = self._alphabet.unrank(i_a)
                        l[self._twin[i][j][0]][self._twin[i][j][1]] = l[i][j]
                        i_a += 1
        return l

    def __eq__(self, other) :
        r"""
        TESTS::

            sage: a0 = [0,0,1]
            sage: a1 = [1,2,2]
            sage: p = iet.GeneralizedPermutation(a0,a1,reduced=True,flips=[0])
            sage: q = copy(p)
            sage: q.alphabet("abc")
            sage: p == q
            True
            sage: b0 = [1,0,0]
            sage: b1 = [2,2,1]
            sage: r = iet.GeneralizedPermutation(b0,b1,reduced=True,flips=[0])
            sage: p == r or q == r
            False
        """
        return (type(self) is type(other) and
            self._twin == other._twin and
            self._flips == other._flips)

    def __ne__(self, other) :
        r"""
        TESTS::

            sage: a0 = [0,0,1]
            sage: a1 = [1,2,2]
            sage: p = iet.GeneralizedPermutation(a0,a1,reduced=True,flips=[0])
            sage: q = copy(p)
            sage: q.alphabet("abc")
            sage: p != q
            False
            sage: b0 = [1,0,0]
            sage: b1 = [2,2,1]
            sage: r = iet.GeneralizedPermutation(b0,b1,reduced=True,flips=[0])
            sage: p != r and q != r
            True
        """
        return (type(self) is not type(other) or
            self._twin != other._twin or
            self._flips != other._flips)

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        For more explanation and a list of arguments try help(iet.RauzyDiagram)

        EXAMPLES::

            sage: p = iet.GeneralizedPermutation('a a b','c c b',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: p in r
            True
        """
        return FlippedReducedRauzyDiagram(self, **kargs)

class ReducedRauzyDiagram(RauzyDiagram):
    r"""
    Rauzy diagram of reduced permutations
    """
    def _permutation_to_vertex(self, p):
        r"""
        The necessary data to store the permutation.

        TESTS::


            sage: p = iet.Permutation('a b c','c b a',reduced=True)   #indirect doctest
            sage: r = p.rauzy_diagram()
            sage: p in r
            True
        """
        return (tuple(p._twin[0]), tuple(p._twin[1]))

    def _set_element(self, data=None):
        r"""
        Sets self._element with data.

        TESTS::

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: p in r   #indirect doctest
            True
        """
        self._element._twin = [list(data[0]), list(data[1])]

class FlippedReducedRauzyDiagram(FlippedRauzyDiagram, ReducedRauzyDiagram):
    r"""
    Rauzy diagram of flipped reduced permutations.
    """
    def _permutation_to_vertex(self, p):
        r"""
        TESTS::

            sage: p = iet.GeneralizedPermutation('a b b','c c a',flips='a',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: p in r   #indirect doctest
            True
        """
        return ((tuple(p._twin[0]), tuple(p._twin[1])),
                (tuple(p._flips[0]), tuple(p._flips[1])))

    def _set_element(self, data=None):
        r"""
        Sets self._element with data.
        TESTS::

            sage: r = iet.RauzyDiagram('a b c','c b a',flips='b',reduced=True)   #indirect doctest
        """
        self._element._twin = [list(data[0][0]), list(data[0][1])]
        self._element._flips = [list(data[1][0]), list(data[1][1])]
