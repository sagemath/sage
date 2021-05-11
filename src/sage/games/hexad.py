r"""
Hexads in S(5,6,12)

This module completes a 5-element subset of a 12-set `X`
into a hexad in a Steiner system `S(5, 6, 12)` using Curtis
and Conway's "kitten method".  The labeling is either the
"modulo 11" labeling or the "shuffle" labeling.

The  main functions implemented in this file are
:meth:`Minimog.blackjack_move` and :meth:`Minimog.find_hexad`.

Enter ``blackjack_move?``
for help to play blackjack (i.e., the rules of the game), or
``find_hexad?`` for help finding hexads of `S(5, 6, 12)` in
the shuffle labeling.

This picture is the kitten in the "shuffle" labeling::


                                     6

                                     9

                                10     8

                             7     2      5

                         9     4     11     9

                      10     8     3      10     8

                 1                                      0


The corresponding MINIMOG is::

             +-----+-----+-----+-----+
             |  6  |  3  |  0  |  9  |
             +-----+-----+-----+-----+
             |  5  |  2  |  7  | 10  |
             +-----+-----+-----+-----+
             |  4  |  1  |  8  | 11  |
             +-----+-----+-----+-----+

which is specified by the global variable ``minimog_shuffle``.

See the docstrings for :meth:`Minimog.find_hexad` and
:meth:`Minimog.blackjack_move` for further details and examples.

AUTHOR:

David Joyner (2006-05)

REFERENCES:

- [Cu1984]_

- [Co1984]_

- [CS1986]_

- [KR2001]_

Some details are also online at:  http://www.permutationpuzzles.org/hexad/
"""
# ****************************************************************************
#       Copyright (C) 2005 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your choice)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.infinity import infinity
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.calculus.calculus import SR


def view_list(L):
    """
    This provides a printout of the lines, crosses and squares
    of the MINIMOG, as in Curtis' paper [Cu1984]_.

    EXAMPLES::

        sage: from sage.games.hexad import *
        sage: M = Minimog(type="shuffle")
        sage: view_list(M.line[1])
        <BLANKLINE>
        [0 0 0]
        [1 1 1]
        [0 0 0]
        sage: view_list(M.cross[1])
        <BLANKLINE>
        [1 1 1]
        [0 0 1]
        [0 0 1]
        sage: view_list(M.square[1])
        <BLANKLINE>
        [0 0 0]
        [1 1 0]
        [1 1 0]
    """
    return matrix(GF(2), 3, 3, lambda x, y: 1 if (x,y) in L else 0)



def picture_set(A, L):
    """
    This is needed in the :meth:`Minimog.find_hexad` function below.

    EXAMPLES::

        sage: from sage.games.hexad import *
        sage: M = Minimog(type="shuffle")
        sage: picture_set(M.picture00, M.cross[2])
        {5, 7, 8, 9, 10}
        sage: picture_set(M.picture02, M.square[7])
        {2, 3, 5, 8}
    """
    return set([A[x] for x in L])


class Minimog(object):
    r"""
    This implements the Conway/Curtis minimog idea for describing
    the Steiner triple system `S(5, 6, 12)`.

    EXAMPLES::

        sage: from sage.games.hexad import *
        sage: Minimog(type="shuffle")
        Minimog of type shuffle
        sage: M = Minimog(type = "modulo11")
        sage: M.minimog
        [        0         3 +Infinity         2]
        [        5         9         8        10]
        [        4         1         6         7]
    """
    def __init__(self, type="shuffle"):
        self.type = type
        MS34 = MatrixSpace(SR,3,4)
        minimog_modulo11 = MS34([[0,3,infinity,2],[5,9,8,10],[4,1,6,7]])
        minimog_shuffle = MS34([[6,3,0,9],[5,2,7,10],[4,1,8,11]])
        if type == "shuffle":
            self.minimog = minimog_shuffle
        elif type == "modulo11":
            self.minimog = minimog_modulo11
        else:
            raise ValueError("that Minimog type is not implemented")
        # This initializes the variables in the game.
        MS34 = MatrixSpace(SR,3,4)
        A = self.minimog
        MS33 = MatrixSpace(SR,3,3)
        self.picture00 = MS33([[A[(1,0)],A[(2,3)],A[(0,1)]],[A[(2,2)],A[(1,1)],A[(2,0)]],[A[(0,3)],A[(1,3)],A[(1,2)]]])
        #######     self.picture00 is the "picture at 6"
        self.picture02 = MS33([[A[(1,0)],A[(2,3)],A[(0,1)]],[A[(1,1)],A[(2,0)],A[(2,2)]],[A[(1,2)],A[(0,3)],A[(1,3)]]])
        #######     self.picture02 is the "picture at 1"
        self.picture21 = MS33([[A[(2,2)],A[(1,3)],A[(0,1)]],[A[(0,3)],A[(2,3)],A[(2,0)]],[A[(1,0)],A[(1,1)],A[(1,2)]]])
        #######     self.picture21 is the "picture at 0"

        self.line = list(range(12))
        self.line[0] = set([(0,0),(0,1),(0,2)])
        self.line[1] = set([(1,0),(1,1),(1,2)])
        self.line[2] = set([(2,0),(2,1),(2,2)])
        self.line[3] = set([(0,2),(1,2),(2,2)])
        self.line[4] = set([(0,1),(1,1),(2,1)])
        self.line[5] = set([(0,0),(1,0),(2,0)])
        self.line[6] = set([(0,0),(1,1),(2,2)])
        self.line[7] = set([(2,0),(0,1),(1,2)])
        self.line[8] = set([(0,2),(1,0),(2,1)])
        self.line[9] = set([(2,0),(1,1),(0,2)])
        self.line[10] = set([(0,0),(1,2),(2,1)])
        self.line[11] = set([(1,0),(0,1),(2,2)])

        self.cross = list(range(18))
        self.cross[0] = set([(0,0),(0,1),(0,2),(1,0),(2,0)])
        self.cross[1] = set([(0,0),(0,1),(0,2),(1,2),(2,2)])
        self.cross[2] = set([(0,0),(1,0),(2,0),(2,1),(2,2)])
        self.cross[3] = set([(2,0),(2,1),(2,2),(0,2),(1,2)])
        self.cross[4] = set([(0,0),(0,1),(0,2),(1,1),(2,1)])
        self.cross[5] = set([(0,0),(1,0),(2,0),(1,1),(1,2)])
        self.cross[6] = set([(1,0),(1,1),(1,2),(0,2),(2,2)])
        self.cross[7] = set([(0,1),(1,1),(2,1),(2,0),(2,2)])
        self.cross[8] = set([(0,0),(0,1),(1,0),(1,1),(2,2)])
        self.cross[9] = set([(0,0),(1,1),(1,2),(2,1),(2,2)])
        self.cross[10] = set([(2,0),(2,1),(1,0),(1,1),(0,2)])
        self.cross[11] = set([(0,1),(0,2),(1,1),(1,2),(2,0)])
        self.cross[12] = set([(0,0),(1,0),(0,2),(1,2),(2,1)])
        self.cross[13] = set([(1,0),(0,1),(0,2),(2,1),(2,2)])
        self.cross[14] = set([(0,1),(1,0),(1,2),(2,0),(2,2)])
        self.cross[15] = set([(0,0),(0,1),(1,2),(2,0),(2,1)])
        self.cross[16] = set([(1,0),(1,1),(1,2),(0,1),(2,1)])
        self.cross[17] = set([(0,0),(0,2),(1,1),(2,0),(2,2)])
        self.box = set([(i,j) for i in range(3) for j in range(3)])
        self.square = [set([]) for i in range(18)]
        for i in range(18):
            self.square[i] = self.box - self.cross[i]

        MS34_GF3 = MatrixSpace(GF(2), 3, 4)
        cols = {}
        cols[1] = MS34_GF3([[1,0,0,0],[1,0,0,0],[1,0,0,0]])
        cols[2] = MS34_GF3([[0,1,0,0],[0,1,0,0],[0,1,0,0]])
        cols[3] = MS34_GF3([[0,0,1,0],[0,0,1,0],[0,0,1,0]])
        cols[4] = MS34_GF3([[0,0,0,1],[0,0,0,1],[0,0,0,1]])
        self.col = cols

        tets = {}
        tets[1] = MS34_GF3([[1,1,1,1],[0,0,0,0],[0,0,0,0]])
        tets[2] = MS34_GF3([[1,0,0,0],[0,1,1,1],[0,0,0,0]])
        tets[3] = MS34_GF3([[1,0,0,0],[0,0,0,0],[0,1,1,1]])
        tets[4] = MS34_GF3([[0,1,0,0],[1,0,1,0],[0,0,0,1]])
        tets[5] = MS34_GF3([[0,0,0,1],[1,1,0,0],[0,0,1,0]])
        tets[6] = MS34_GF3([[0,0,1,0],[1,0,0,1],[0,1,0,0]])
        tets[7] = MS34_GF3([[0,1,0,0],[0,0,0,1],[1,0,1,0]])
        tets[8] = MS34_GF3([[0,0,1,0],[0,1,0,0],[1,0,0,1]])
        tets[9] = MS34_GF3([[0,0,0,1],[0,0,1,0],[1,1,0,0]])
        self.tet = tets

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = Minimog(type="modulo11")
            sage: M
            Minimog of type modulo11
        """
        return "Minimog of type %s" % self.type

    def __str__(self):
        """
        EXAMPLES::

            sage: M = Minimog(type="modulo11")
            sage: print(M)
            Minimog of type modulo11 associated to
            [        0         3 +Infinity         2]
            [        5         9         8        10]
            [        4         1         6         7]

        """
        return "Minimog of type %s associated to\n %s" % (self.type,
                                                          self.minimog)

    def _latex_(self):
        r"""
        Prints latex code.

        EXAMPLES::

            sage: M = Minimog(type="modulo11")
            sage: latex(M)
            Minimog of type modulo11 associated to
             $\left(\begin{array}{rrrr}
            0 & 3 & +\infty & 2 \\
            5 & 9 & 8 & 10 \\
            4 & 1 & 6 & 7
            \end{array}\right)$
        """
        from sage.misc.latex import latex
        return "Minimog of type %s associated to\n $%s$" % (self.type,
                                                            latex(self.minimog))

    def print_kitten(self):
        """
        This simply prints the "kitten" (expressed as a triangular diagram
        of symbols).

        EXAMPLES::

            sage: from sage.games.hexad import *
            sage: M = Minimog("shuffle")
            sage: M.print_kitten()
                     0
            <BLANKLINE>
                     8
                   9  10
                  5  11  3
                 8  2  4  8
               9  10  7  9  10
            <BLANKLINE>
            6                   1
            sage: M = Minimog("modulo11")
            sage: M.print_kitten()
                     +Infinity
            <BLANKLINE>
                     6
                   2  10
                  5  7  3
                 6  9  4  6
               2  10  8  2  10
            <BLANKLINE>
            0                   1

        """
        MINIMOG = self.minimog
        kitten = '         {}'.format(MINIMOG[0][2])
        kitten += '\n           '
        kitten += '\n         {}'.format(MINIMOG[2][2])
        kitten += '\n       {}  {}'.format(MINIMOG[0][3], MINIMOG[1][3])
        kitten += '\n      {}  {}  {}'.format(MINIMOG[1][0], MINIMOG[2][3], MINIMOG[0][1])
        kitten += '\n     {0}  {1}  {2}  {0}'.format(MINIMOG[2][2], MINIMOG[1][1], MINIMOG[2][0])
        kitten += '\n   {0}  {1}  {2}  {0}  {1}'.format(MINIMOG[0][3], MINIMOG[1][3], MINIMOG[1][2])
        kitten += '\n           \n'
        kitten += '{}                   {}'.format(MINIMOG[0][0], MINIMOG[2][1])
        print(kitten)

    def find_hexad0(self, pts):
        """
        Find a hexad of type 0.

        INPUT:

        - ``pts`` -- a set of 2 distinct elements of MINIMOG, but not
          including the "points at infinity"

        OUTPUT:

        hexad containing ``pts`` and of type 0 (the 3 points "at
        infinity" union a line)

        .. NOTE::

            The 3 points "at infinity" are
            ``{MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}``.

        EXAMPLES::

            sage: from sage.games.hexad import *
            sage: M = Minimog(type="shuffle")
            sage: M.find_hexad0(set([2,4]))
            ([0, 1, 2, 4, 6, 8], ['line 1', 'picture 1'])

        """
        MINIMOG = self.minimog
        L = set(pts)
        H = set([MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]])
        for i in range(12):
            if L <= picture_set(self.picture02, self.line[i]):
                WHAT = ["line " + str(i), "picture " + str(1)]
                H = H | picture_set(self.picture02, self.line[i])
                return list(H), WHAT
            if L <= picture_set(self.picture21, self.line[i]):
                WHAT = ["line " + str(i), "picture " + str(0)]
                H = H | picture_set(self.picture21, self.line[i])
                return list(H), WHAT
            if L <= picture_set(self.picture00, self.line[i]):
                WHAT = ["line " + str(i), "picture " + str(6)]
                H = H | picture_set(self.picture00, self.line[i])
                return list(H), WHAT
        return [], []

    def find_hexad1(self, pts):
        """
        Find a hexad of type 1.

        INPUT:

        - ``pts`` -- a set of 5 distinct elements of MINIMOG

        OUTPUT:

        hexad containing ``pts`` and of type 1 (union of 2 parallel
        lines -- *no* points "at infinity")

        .. NOTE::

            The 3 points "at infinity" are
            ``{MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}``.

        EXAMPLES::

            sage: from sage.games.hexad import *
            sage: M = Minimog(type="shuffle")
            sage: M.find_hexad1(set([2,3,4,5,8]))
            ([2, 3, 4, 5, 8, 11], ['lines (1, 2)', 'picture 1'])
        """
        H = set(pts)
        L = set(pts)
        linez = [(1, 2), (1, 3), (2, 3), (4, 5), (4, 6), (5, 6),
                 (7, 8), (7, 9), (8, 9), (10, 11), (10, 12), (11, 12)]
        for x in linez:
            x1 = int(x[0] - 1)
            x2 = int(x[1] - 1)      # (recall | is union)
            if L <= (picture_set(self.picture02, self.line[x1]) | picture_set(self.picture02,self.line[x2])):
                WHAT = ["lines " + str(x), "picture " + str(1)]
                H = picture_set(self.picture02, self.line[x1]) | picture_set(self.picture02,self.line[x2])
                return list(H), WHAT
            if L <= (picture_set(self.picture21, self.line[x1]) | picture_set(self.picture21,self.line[x2])):
                WHAT = ["lines " + str(x), "picture " + str(0)]
                H = picture_set(self.picture21, self.line[x1]) | picture_set(self.picture21,self.line[x2])
                return list(H), WHAT
            if L <= (picture_set(self.picture00, self.line[x1]) | picture_set(self.picture00,self.line[x2])):
                WHAT = ["lines " + str(x), "picture " + str(6)]
                H = picture_set(self.picture00, self.line[x1]) | picture_set(self.picture00,self.line[x2])
                return list(H), WHAT
        return [],[]

    def find_hexad2(self, pts, x0):
        r"""
        Find a hexad of type 2.

        INPUT:

        - ``pts`` -- a list S of 4 elements of MINIMOG, not including
          any "points at infinity"
        - ``x0`` -- in ``{MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}``

        OUTPUT:

        hexad containing `S \cup \{x0\}` of type 2

        EXAMPLES::

            sage: from sage.games.hexad import *
            sage: M = Minimog(type="shuffle")
            sage: M.find_hexad2([2,3,4,5],1)
            ([], [])

        The above output indicates that there is no hexad of type 2
        containing `\{2,3,4,5\}`. However, there is one containing
        `\{2,3,4,8\}`::

            sage: M.find_hexad2([2,3,4,8],0)
            ([0, 2, 3, 4, 8, 9], ['cross 12', 'picture 0'])
        """
        MINIMOG = self.minimog
        L = set(pts)
        H = set([x0])
        for i in range(18):
            if (x0 == MINIMOG[2][1] and L <= picture_set(self.picture02,self.cross[i])):
                WHAT = ["cross " + str(i), "picture " + str(1)]
                H = H | picture_set(self.picture02, self.cross[i])
                return list(H), WHAT
            if (x0 == MINIMOG[0][2] and L <= picture_set(self.picture21,self.cross[i])):
                WHAT = ["cross " + str(i), "picture " + str(MINIMOG[0][2])]
                H = H | picture_set(self.picture21, self.cross[i])
                return list(H), WHAT
            if (x0 == MINIMOG[0][0] and L <= picture_set(self.picture00,self.cross[i])):
                WHAT = ["cross " + str(i), "picture " + str(6)]
                H = H | picture_set(self.picture00, self.cross[i])
                return list(H), WHAT
        return [],[]

    def find_hexad3(self, pts, x0, x1):
        r"""
        Find a hexad of type 3.

        INPUT:

        - ``pts`` -- a list of 3 elements of MINIMOG, not including any
          "points at infinity"
        - ``x0``, ``x1``  -- in ``{MINIMOG[0][2], MINIMOG[2][1],
          MINIMOG[0][0]}``

        OUTPUT:

        hexad containing pts union \{x0,x1\} of type 3 (square at
        picture of "omitted point at infinity")

        EXAMPLES::

            sage: from sage.games.hexad import *
            sage: M = Minimog(type="shuffle")
            sage: M.find_hexad3([2,3,4],0,1)
            ([0, 1, 2, 3, 4, 11], ['square 2', 'picture 6'])

        """
        MINIMOG = self.minimog
        L = set(pts)
        H = set([x0, x1])
        for i in range(18):
            if (not (MINIMOG[0][2] in H) and L <= picture_set(self.picture21, self.square[i])):
                WHAT = ["square " + str(i), "picture " + str(MINIMOG[0][2])]
                H = H | picture_set(self.picture21, self.square[i])
                return list(H), WHAT
            if (not (MINIMOG[2][1] in H) and L <= picture_set(self.picture02, self.square[i])):
                WHAT = ["square " + str(i), "picture " + str(MINIMOG[2][1])]
                H = H | picture_set(self.picture02, self.square[i])
                return list(H), WHAT
            if (not (MINIMOG[0][0] in H) and L <= picture_set(self.picture00, self.square[i])):
                WHAT = ["square " + str(i), "picture " + str(MINIMOG[0][0])]
                H = H | picture_set(self.picture00, self.square[i])
                return list(H), WHAT
        return [], []

    def find_hexad(self, pts):
        r"""
        Find a hexad of some type.

        INPUT:

        - ``pts`` -- a list S of 5 elements of MINIMOG

        OUTPUT:

        hexad containing `S \cup \{x0\}` of some type

        .. NOTE::

            The 3 "points at infinity" are
            ``{MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}``.

        Theorem ([Cu1984]_,  [Co1984]_): Each hexads is of exactly one
        of the following types:

        0. {3 "points at infinity"} union {any line},
        1. the union of any two (distinct) parallel lines in the same
           picture,
        2. one "point at infinity" union a cross in the corresponding
           picture, or
        3. two "points at infinity" union a square in the picture
           corresponding to the omitted point at infinity.

        More precisely, there are 132 such hexads (12 of type 0,
        12 of type 1, 54 of type 2, and 54 of type 3).
        They form a Steiner system of type `(5,6,12)`.

        EXAMPLES::

            sage: from sage.games.hexad import *
            sage: M = Minimog(type="shuffle")
            sage: M.find_hexad([0,1,2,3,4])
            ([0, 1, 2, 3, 4, 11], ['square 2', 'picture 6'])
            sage: M.find_hexad([1,2,3,4,5])
            ([1, 2, 3, 4, 5, 6], ['square 8', 'picture 0'])
            sage: M.find_hexad([2,3,4,5,8])
            ([2, 3, 4, 5, 8, 11], ['lines (1, 2)', 'picture 1'])
            sage: M.find_hexad([0,1,2,4,6])
            ([0, 1, 2, 4, 6, 8], ['line 1', 'picture 1'])
            sage: M = Minimog(type="modulo11")
            sage: M.find_hexad([1,2,3,4,SR(infinity)]) # random (machine dependent?) order
            ([+Infinity, 2, 3, 4, 1, 10], ['square 8', 'picture 0'])

        AUTHOR:

        David Joyner (2006-05)

        REFERENCES: [Cu1984]_,  [Co1984]_
        """
        MINIMOG = self.minimog
        L = set(pts)
        LL = L.copy()
        pts_at_infty = set([MINIMOG[0][2],MINIMOG[2][1],MINIMOG[0][0]])
        # recall & means intersection
        L2 = LL & pts_at_infty
        if len(L2) == 3:  # must be type 0 (line + pts at infty)
            H, WHAT = self.find_hexad0(LL - pts_at_infty)
            return H, WHAT
        if len(L2) == 2:  # type 0 or 3
            if (MINIMOG[0][2] in LL and MINIMOG[2][1] in LL):
                H, WHAT = self.find_hexad3(LL - set([MINIMOG[0][2],MINIMOG[2][1]]),
                                           MINIMOG[0][2], MINIMOG[2][1])
                if H != []:   # must be type 3
                    return list(H), WHAT
                if H == []:   # could be type 0
                    H, WHAT = self.find_hexad0(LL - L2)
                    if H != []:   # must be type 0
                        return list(H), WHAT
            if (MINIMOG[2][1] in LL and MINIMOG[0][0] in LL):
                H, WHAT = self.find_hexad3(LL - set([MINIMOG[2][1],MINIMOG[0][0]]),
                                           MINIMOG[2][1], MINIMOG[0][0])
                if H != []:   # must be type 3
                    return list(H), WHAT
                if H == []:   # could be type 0
                    H, WHAT = self.find_hexad0(LL - L2)
                    if H != []:   # must be type 0
                        return list(H), WHAT
            if (MINIMOG[0][2] in LL and MINIMOG[0][0] in LL):
                H, WHAT = self.find_hexad3(LL - set([MINIMOG[0][2],MINIMOG[0][0]]),
                                           MINIMOG[0][2], MINIMOG[0][0])
                if H != []:   # must be type 3
                    return list(H), WHAT
                if H == []:   # could be type 0
                    H, WHAT = self.find_hexad0(LL - L2)
                    if H != []:   # must be type 0
                        return list(H), WHAT
        if len(L2) == 1:
            H, WHAT = self.find_hexad2(LL - L2, list(L2)[0])
            if H == []:   # not a cross in picture at infinity
                if list(L2)[0] == MINIMOG[2][1]:
                    L1 = LL - L2
                    H, WHAT = self.find_hexad3(L1, MINIMOG[0][0], MINIMOG[2][1])
                    if H != []:
                        return list(H), WHAT
                    L1 = LL - L2
                    H, WHAT = self.find_hexad3(L1, MINIMOG[0][2], MINIMOG[2][1])
                    if H != []:
                        return list(H), WHAT
                if list(L2)[0] == MINIMOG[0][0]:
                    L1 = (LL - L2)
                    H, WHAT = self.find_hexad3(L1, MINIMOG[0][0], MINIMOG[2][1])
                    if H != []:
                        return list(H), WHAT
                    L1 = (LL - L2)
                    H, WHAT = self.find_hexad3(L1, MINIMOG[0][0], MINIMOG[0][2])
                    if H != []:
                        return list(H), WHAT
                if list(L2)[0] == MINIMOG[0][2]:
                    L1 = (LL - L2)
                    H, WHAT = self.find_hexad3(L1, MINIMOG[0][0], MINIMOG[0][2])
                    if H != []:
                        return list(H), WHAT
                    L1 = (LL - L2)
                    H, WHAT = self.find_hexad3(L1, MINIMOG[2][1], MINIMOG[0][2])
                    if H != []:
                        return list(H), WHAT
            return list(H), WHAT
            # a cross in a pic at infty
        if not L2:  # L is either a union of 2 lines or a cross
            for i in LL:
                for j in pts_at_infty:
                    H, WHAT = self.find_hexad2(LL - set([i]),j)
                    if (H != [] and i in H):
                        return list(H), WHAT  # L is in a cross
            H, WHAT = self.find_hexad1(LL)  # L is a union of lines
            return H, WHAT

    def blackjack_move(self, L0):
        r"""
        Perform a blackjack move.

        INPUT:

        - ``L0`` -- a list of cards of length 6, taken
          from `\{0, 1, ..., 11\}`

        .. RUBRIC:: MATHEMATICAL BLACKJACK

        Mathematical blackjack is played with 12 cards, labeled `0, ..., 11`
        (for example:  king,  ace, `2`, `3`, ..., `10`, jack, where the
        king is `0` and the  jack is `11`). Divide the 12 cards into two
        piles of `6` (to be fair, this should be done randomly). Each of
        the `6` cards of one of these piles are to be placed face up on
        the table. The remaining cards are in a stack which is shared
        and visible to both players. If the sum of the cards face up on
        the table is less than 21 then no legal move is possible so you
        must shuffle the cards and deal a new game. (Conway calls such
        a game `*={0|0}`, where `0={|}`; in this game the first player
        automatically wins.)

        * Players alternate moves.
        * A move consists of exchanging a card on the table with a
          lower card from the other pile.
        * The player whose move makes the sum of the cards on the table
          under 21 loses.

        The winning strategy (given below) for this game is due to
        Conway and Ryba. There is a Steiner system `S(5,6,12)` of hexads
        in the set `\{0, 1, ..., 11\}`. This Steiner system is associated
        to the MINIMOG of in the "shuffle numbering" rather than the
        "modulo `11` labeling".

        **Proposition** ([KR2001]_) For this Steiner system, the
        winning strategy is to choose a move which is a hexad from
        this system.

        EXAMPLES::

            sage: M = Minimog(type="modulo11")
            sage: M.blackjack_move([0,2,3,6,1,10])
            '6 --> 5. The total went from 22 to 21.'
            sage: M = Minimog(type="shuffle")
            sage: M.blackjack_move([0,2,4,6,7,11])
            '4 --> 3. The total went from 30 to 29.'

        Is this really a hexad? ::

            sage: M.find_hexad([11,2,3,6,7])
            ([0, 2, 3, 6, 7, 11], ['square 9', 'picture 1'])

        So, yes it is, but here is further confirmation::

            sage: M.blackjack_move([0,2,3,6,7,11])
            This is a hexad.
            There is no winning move, so make a random legal move.
            [0, 2, 3, 6, 7, 11]

        Now, suppose player 2 replaced the 11 by a 9. Your next move::

            sage: M.blackjack_move([0,2,3,6,7,9])
            '7 --> 1. The total went from 27 to 21.'

        You have now won. Sage will even tell you so::

            sage: M.blackjack_move([0,2,3,6,1,9])
            'No move possible. Shuffle the deck and redeal.'

        AUTHOR:

        David Joyner (2006-05)

        REFERENCES: [CS1986]_, [KR2001]_
        """
        total = sum(L0)
        if total < 22:
            return "No move possible. Shuffle the deck and redeal."
        L = set(L0)
        for x in L:
            h, WHAT = self.find_hexad(L - set([x]))
            if list(L0) == list(h):
                print("      This is a hexad. \n      There is no winning move, so make a random legal move.")
                return L0
            y = list(set(h) - (L - set([x])))[0]
            if y < x:
                return str(x) + ' --> ' + str(y) + ". The total went from " + str(total) + " to " + str(total - x + y) + "."
        print("This is a hexad. \n There is no winning move, so make a random legal move.")
        return L0

