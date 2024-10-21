# sage_setup: distribution = sagemath-categories
"""
Dancing links C++ wrapper
"""
# ****************************************************************************
#       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>,
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

# OneExactCover and AllExactCovers are almost exact copies of the
# functions with the same name in sage/combinat/dlx.py by Tom Boothby.

from .dancing_links import dlx_solver


def DLXCPP(rows):
    """
    Solve the Exact Cover problem by using the Dancing Links algorithm
    described by Knuth.

    Consider a matrix M with entries of 0 and 1, and compute a subset
    of the rows of this matrix which sum to the vector of all 1s.

    The dancing links algorithm works particularly well for sparse
    matrices, so the input is a list of lists of the form::

       [
        [i_11,i_12,...,i_1r]
        ...
        [i_m1,i_m2,...,i_ms]
       ]

    where M[j][i_jk] = 1.

    The first example below corresponds to the matrix::

       1110
       1010
       0100
       0001

    which is exactly covered by::

       1110
       0001

    and

    ::

       1010
       0100
       0001

    If soln is a solution given by DLXCPP(rows) then

    [ rows[soln[0]], rows[soln[1]], ... rows[soln[len(soln)-1]] ]

    is an exact cover.

    Solutions are given as a list.

    EXAMPLES::

        sage: rows = [[0,1,2]]
        sage: rows+= [[0,2]]
        sage: rows+= [[1]]
        sage: rows+= [[3]]
        sage: [x for x in DLXCPP(rows)]
        [[3, 0], [3, 1, 2]]
    """
    if not rows:
        return

    x = dlx_solver(rows)

    while x.search():
        yield x.get_solution()


def AllExactCovers(M):
    """
    Solve the exact cover problem on the matrix M (treated as a dense
    binary matrix).

    EXAMPLES: No exact covers::

        sage: M = Matrix([[1,1,0],[1,0,1],[0,1,1]])                                     # needs sage.modules
        sage: [cover for cover in AllExactCovers(M)]                                    # needs sage.modules
        []

    Two exact covers::

        sage: M = Matrix([[1,1,0],[1,0,1],[0,0,1],[0,1,0]])                             # needs sage.modules
        sage: [cover for cover in AllExactCovers(M)]                                    # needs sage.modules
        [[(1, 1, 0), (0, 0, 1)], [(1, 0, 1), (0, 1, 0)]]
    """
    rows = []
    for R in M.rows():
        row = []
        for i in range(len(R)):
            if R[i]:
                row.append(i)
        rows.append(row)
    for s in DLXCPP(rows):
        yield [M.row(i) for i in s]


def OneExactCover(M):
    """
    Solve the exact cover problem on the matrix M (treated as a dense
    binary matrix).

    EXAMPLES::

        sage: # needs sage.modules
        sage: M = Matrix([[1,1,0],[1,0,1],[0,1,1]])  # no exact covers
        sage: print(OneExactCover(M))
        None
        sage: M = Matrix([[1,1,0],[1,0,1],[0,0,1],[0,1,0]]) # two exact covers
        sage: OneExactCover(M)
        [(1, 1, 0), (0, 0, 1)]
    """

    for s in AllExactCovers(M):
        return s
