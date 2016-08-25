r"""
Alternating Sign Matrices

AUTHORS:

- Mike Hansen (2007): Initial version
- Pierre Cange, Luis Serrano (2012): Added monotone triangles
- Travis Scrimshaw (2013-28-03): Added element class for ASM's and made
  :class:`MonotoneTriangles` inherit from :class:`GelfandTsetlinPatterns`
- Jessica Striker (2013): Added additional methods
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2012 Pierre Cagne <pierre.cagne@ens.fr>,
#                          Luis Serrano <luisgui.serrano@gmail.com>
#                     2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#                     2013 Jessica Striker <jessicapalencia@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
# python3
from __future__ import division

import itertools
import copy
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.flatten import flatten
from sage.misc.all import prod
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.misc.all import cached_method
from sage.rings.all import ZZ
from sage.arith.all import factorial
from sage.rings.integer import Integer
from sage.combinat.posets.lattices import LatticePoset
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPatternsTopRow
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.non_decreasing_parking_function import NonDecreasingParkingFunction
from sage.combinat.permutation import Permutation
from sage.combinat.six_vertex_model import SquareIceModel
from sage.misc.decorators import rename_keyword


class AlternatingSignMatrix(Element):
    r"""
    An alternating sign matrix.

    An alternating sign matrix is a square matrix of `0`'s, `1`'s and `-1`'s
    such that the sum of each row and column is `1` and the non-zero
    entries in each row and column alternate in sign.

    These were introduced in [MiRoRu]_.

    REFERENCES:

    .. [MiRoRu] \W. H. Mills, David P Robbins, Howard Rumsey Jr.,
       *Alternating sign matrices and descending plane partitions*,
       Journal of Combinatorial Theory, Series A,
       Volume 34, Issue 3, May 1983, Pages 340--359.
       http://www.sciencedirect.com/science/article/pii/0097316583900687
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, asm):
        """
        Create an ASM.

        EXAMPLES::

            sage: AlternatingSignMatrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        asm = matrix(asm)
        if not asm.is_square():
            raise ValueError("The alternating sign matrices must be square")
        P = AlternatingSignMatrices(asm.nrows())
        if asm not in P:
            raise ValueError("Invalid alternating sign matrix")
        return P(asm)

    def __init__(self, parent, asm):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: elt = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: TestSuite(elt).run()
        """
        self._matrix = asm
        Element.__init__(self, parent)

    def __hash__(self):
        r"""
        TESTS::

            sage: A = AlternatingSignMatrices(3)
            sage: elt = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: hash(elt)
            12
        """
        return hash(self._matrix)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return repr(self._matrix)

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: M == A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            True
            sage: M == A([[1, 0, 0],[0, 0, 1],[0, 1, 0]])
            False
        """
        if isinstance(other, AlternatingSignMatrix):
            return self._matrix == other._matrix
        return self._matrix == other

    def __ne__(self, other):
        """
        Check not equals. This is needed, see :trac:`14762`.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: M != A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            False
            sage: M != A([[1, 0, 0],[0, 0, 1],[0, 1, 0]])
            True
        """
        return not self == other

    def __le__(self, other):
        """
        Check less than or equal to. This is needed, see :trac:`15372`.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: M <= A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            True
            sage: M <= A([[1, 0, 0],[0, 0, 1],[0, 1, 0]])
            False
        """
        if isinstance(other, AlternatingSignMatrix):
            return self._matrix <= other._matrix
        return False #return False if other is not an ASM

    def __lt__(self, other):
        """
        Check less than. This is needed, see :trac:`15372`.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: M < A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            False
        """
        if isinstance(other, AlternatingSignMatrix):
            return self._matrix < other._matrix
        return False #return False if other is not an ASM

    def __ge__(self, other):
        """
        Check greater than or equal to. This is needed, see :trac:`15372`.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: M >= A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            True
            sage: M >= A([[1, 0, 0],[0, 0, 1],[0, 1, 0]])
            True
        """
        if isinstance(other, AlternatingSignMatrix):
            return self._matrix >= other._matrix
        return False #return False if other is not an ASM

    def __gt__(self, other):
        """
        Check greater than. This is needed, see :trac:`15372`.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: M > A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            False
        """
        if isinstance(other, AlternatingSignMatrix):
            return self._matrix > other._matrix
        return False #return False if other is not an ASM

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: latex(A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]))
            \left(\begin{array}{rrr}
            1 & 0 & 0 \\
            0 & 1 & 0 \\
            0 & 0 & 1
            \end{array}\right)
        """
        return self._matrix._latex_()

    def to_matrix(self):
        """
        Return ``self`` as a regular matrix.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: asm = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: m = asm.to_matrix(); m
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: m.parent()
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
        """
        return copy.copy(self._matrix)

    @combinatorial_map(name='to monotone triangle')
    def to_monotone_triangle(self):
        r"""
        Return a monotone triangle from ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).to_monotone_triangle()
            [[3, 2, 1], [2, 1], [1]]
            sage: asm = A([[0, 1, 0],[1, -1, 1],[0, 1, 0]])
            sage: asm.to_monotone_triangle()
            [[3, 2, 1], [3, 1], [2]]
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.to_monotone_triangle()
            [[3, 2, 1], [3, 1], [3]]
            sage: A.from_monotone_triangle(asm.to_monotone_triangle()) == asm
            True
        """
        n = self._matrix.nrows()
        triangle = [None]*n
        prev = [0]*n
        for j, row in enumerate(self._matrix):
            add_row = [a+b for (a,b) in itertools.izip(row, prev)]
            line = [i+1 for (i,val) in enumerate(add_row) if val==1]
            triangle[n-1-j] = list(reversed(line))
            prev = add_row
        return MonotoneTriangles(n)(triangle)

    @combinatorial_map(name='rotate counterclockwise')
    def rotate_ccw(self):
        r"""
        Return the counterclockwise quarter turn rotation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).rotate_ccw()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.rotate_ccw()
            [1 0 0]
            [0 0 1]
            [0 1 0]
        """
        l = list(self._matrix.transpose())
        l.reverse()
        return AlternatingSignMatrix(matrix(l))

    def inversion_number(self):
        r"""
        Return the inversion number of ``self``.

        If we denote the entries of the alternating sign matrix as `a_{i,j}`,
        the inversion number is defined as `\sum_{i>k}\sum_{j<l}a_{i,j}a_{k,l}`.
        When restricted to permutation matrices, this gives the usual inversion
        number of the permutation.

        This definition is equivalent to the one given in [MiRoRu]_.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).inversion_number()
            0
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.inversion_number()
            2
            sage: asm = A([[0, 1, 0],[1, -1, 1],[0, 1, 0]])
            sage: asm.inversion_number()
            2
            sage: P=Permutations(5)
            sage: all(p.number_of_inversions()==AlternatingSignMatrix(p.to_matrix()).inversion_number() for p in P)
            True
        """
        inversion_num = 0
        asm_matrix = self.to_matrix()
        nonzero_cells = asm_matrix.nonzero_positions()
        for (i,j) in nonzero_cells:
            for (k,l) in nonzero_cells:
                if i > k and j < l:
                    inversion_num += asm_matrix[i][j]*asm_matrix[k][l]
        return inversion_num

    @combinatorial_map(name='rotate clockwise')
    def rotate_cw(self):
        r"""
        Return the clockwise quarter turn rotation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).rotate_cw()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.rotate_cw()
            [0 1 0]
            [1 0 0]
            [0 0 1]
        """
        l = list(self._matrix.transpose())
        l.reverse()
        return AlternatingSignMatrix(matrix(l).transpose().antitranspose())

    @combinatorial_map(name='transpose')
    def transpose(self):
        r"""
        Return the counterclockwise quarter turn rotation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).transpose()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.transpose()
            [0 1 0]
            [0 0 1]
            [1 0 0]
        """
        return AlternatingSignMatrix(self._matrix.transpose())

    def corner_sum_matrix(self):
        r"""
        Return the corner sum matrix from ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).corner_sum_matrix()
            [0 0 0 0]
            [0 1 1 1]
            [0 1 2 2]
            [0 1 2 3]
            sage: asm = A([[0, 1, 0],[1, -1, 1],[0, 1, 0]])
            sage: asm.corner_sum_matrix()
            [0 0 0 0]
            [0 0 1 1]
            [0 1 1 2]
            [0 1 2 3]
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.corner_sum_matrix()
            [0 0 0 0]
            [0 0 0 1]
            [0 1 1 2]
            [0 1 2 3]

        TESTS:

        Some non-symmetric tests::

            sage: A = AlternatingSignMatrices(3)
            sage: asm = A([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
            sage: asm.corner_sum_matrix()
            [0 0 0 0]
            [0 0 1 1]
            [0 0 1 2]
            [0 1 2 3]
            sage: B = AlternatingSignMatrices(4)
            sage: asm = B([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, -1, 1], [0, 0, 1, 0]])
            sage: asm.corner_sum_matrix()
            [0 0 0 0 0]
            [0 0 0 1 1]
            [0 1 1 2 2]
            [0 1 2 2 3]
            [0 1 2 3 4]

        """
        asm = self.to_matrix()
        n = asm.nrows() + 1
        return matrix([[nw_corner_sum(asm,i,j) for j in range(n)] for i in range(n)])

    def height_function(self):
        r"""
        Return the height function from ``self``. A height function
        corresponding to an `n \times n` ASM is an `(n+1) \times (n+1)` matrix
        such that the first row is `0, 1, \ldots, n`, the last row is
        `n, n-1, \ldots, 1, 0`, and the difference between adjacent entries
        is 1.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).height_function()
            [0 1 2 3]
            [1 0 1 2]
            [2 1 0 1]
            [3 2 1 0]
            sage: asm = A([[0, 1, 0],[1, -1, 1],[0, 1, 0]])
            sage: asm.height_function()
            [0 1 2 3]
            [1 2 1 2]
            [2 1 2 1]
            [3 2 1 0]
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.height_function()
            [0 1 2 3]
            [1 2 1 2]
            [2 3 2 1]
            [3 2 1 0]
        """
        asm = self.to_matrix()
        n = asm.nrows() + 1
        return matrix([[i+j-2*nw_corner_sum(asm,i,j) for i in range(n)] for j in range(n)])

    def to_six_vertex_model(self):
        r"""
        Return the six vertex model configuration from ``self``.
        This method calls :meth:`sage.combinat.six_vertex_model.from_alternating_sign_matrix`.

        EXAMPLES::

            sage: asm = AlternatingSignMatrix([[0,1,0],[1,-1,1],[0,1,0]])
            sage: asm.to_six_vertex_model()
                ^    ^    ^
                |    |    |
            --> # -> # <- # <--
                ^    |    ^
                |    V    |
            --> # <- # -> # <--
                |    ^    |
                V    |    V
            --> # -> # <- # <--
                |    |    |
                V    V    V

        TESTS::

            sage: ASM = AlternatingSignMatrices(5)
            sage: all((x.to_six_vertex_model()).to_alternating_sign_matrix() == x
            ....:     for x in ASM)
            True
        """

        asm = self.to_matrix()
        n = asm.nrows()
        M = SquareIceModel(n)
        return M.from_alternating_sign_matrix(self)

    def to_fully_packed_loop(self):
        r"""
        Return the fully packed loop configuration from ``self``.

        .. SEEALSO::

            :class:FullyPackedLoop

        EXAMPLES::

            sage: asm = AlternatingSignMatrix([[1,0,0],[0,1,0],[0,0,1]])
            sage: fpl = asm.to_fully_packed_loop()
            sage: fpl
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |

        """
        from sage.combinat.fully_packed_loop import FullyPackedLoop
        return FullyPackedLoop(self)

    def link_pattern(self):
        """
        Return the link pattern corresponding to the fully packed loop 
        corresponding to self.

        EXAMPLES:

        We can extract the underlying link pattern (a non-crossing
        partition) from a fully packed loop::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: A.link_pattern()
            [(1, 2), (3, 6), (4, 5)]

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: B.link_pattern()
            [(1, 6), (2, 5), (3, 4)]
        """
        return self.to_fully_packed_loop().link_pattern()

    @combinatorial_map(name='gyration')
    def gyration(self):
        r"""
        Return the alternating sign matrix obtained by applying gyration
        to the height function in bijection with ``self``.

        Gyration acts on height functions as follows. Go through the entries of
        the matrix, first those for which the sum of the row and column indices
        is even, then for those for which it is odd, and increment or decrement
        the squares by 2 wherever possible such that the resulting matrix is
        still a height function. Gyration was first defined in [Wieland00]_ as
        an action on fully-packed loops.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).gyration()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: asm = A([[0, 1, 0],[1, -1, 1],[0, 1, 0]])
            sage: asm.gyration()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: asm = A([[0, 0, 1],[1, 0, 0],[0, 1, 0]])
            sage: asm.gyration()
            [0 1 0]
            [0 0 1]
            [1 0 0]
            sage: A = AlternatingSignMatrices(3)
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).gyration().gyration()
            [ 0  1  0]
            [ 1 -1  1]
            [ 0  1  0]
            sage: A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]).gyration().gyration().gyration()
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: A = AlternatingSignMatrices(4)
            sage: M = A([[0,0,1,0],[1,0,0,0],[0,1,-1,1],[0,0,1,0]])
            sage: for i in range(5):
            ....:     M = M.gyration()
            sage: M
            [1 0 0 0]
            [0 0 0 1]
            [0 1 0 0]
            [0 0 1 0]
        """
        A = self.parent()
        hf = list(self.height_function())
        k = len(hf) - 1
        for i in range(1,k):
            for j in range(1,k):
                if (i+j) % 2 == 0 \
                        and hf[i-1][j] == hf[i+1][j] == hf[i][j+1] == hf[i][j-1]:
                    if hf[i][j] < hf[i+1][j]:
                        hf[i][j] += 2
                    else:
                        hf[i][j] -= 2
        for i in range(1,k):
            for j in range(1,k):
                if (i+j) % 2 == 1 \
                        and hf[i-1][j] == hf[i+1][j] == hf[i][j+1] == hf[i][j-1]:
                    if hf[i][j] < hf[i+1][j]:
                        hf[i][j] += 2
                    else:
                        hf[i][j] -= 2
        return A.from_height_function(matrix(hf))

    def gyration_orbit(self):
        r"""
        Return the gyration orbit of ``self`` (including ``self``)

        EXAMPLES::

            sage: AlternatingSignMatrix([[0,1,0],[1,-1,1],[0,1,0]]).gyration_orbit()
            [
            [ 0  1  0]  [1 0 0]  [0 0 1]
            [ 1 -1  1]  [0 1 0]  [0 1 0]
            [ 0  1  0], [0 0 1], [1 0 0]
            ]

            sage: AlternatingSignMatrix([[0,1,0,0],[1,-1,1,0],[0,1,-1,1],[0,0,1,0]]).gyration_orbit()
            [
            [ 0  1  0  0]  [1 0 0 0]  [ 0  0  1  0]  [0 0 0 1]
            [ 1 -1  1  0]  [0 1 0 0]  [ 0  1 -1  1]  [0 0 1 0]
            [ 0  1 -1  1]  [0 0 1 0]  [ 1 -1  1  0]  [0 1 0 0]
            [ 0  0  1  0], [0 0 0 1], [ 0  1  0  0], [1 0 0 0]
            ]

            sage: len(AlternatingSignMatrix([[0,1,0,0,0,0],[0,0,1,0,0,0],[1,-1,0,0,0,1],\
            [0,1,0,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0]]).gyration_orbit())
            12
        """
        cyc = [self]
        B = self.gyration()
        while self != B:
            cyc.append(B)
            B = B.gyration()
        return cyc

    def ASM_compatible(self, B):
        r"""
        Return ``True`` if ``self`` and ``B`` are compatible alternating sign
        matrices in the sense of [EKLP92]_. (If ``self`` is of size `n`, ``B``
        must  be of size `n+1`.)

        In [EKLP92]_, there is a notion of a pair of ASM's with sizes differing
        by 1 being compatible, in the sense that they can be combined to encode
        a tiling of the Aztec Diamond.

        REFERENCES:

        .. [EKLP92] \N. Elkies, G. Kuperberg, M. Larsen, J. Propp,
           *Alternating-Sign Matrices and Domino Tilings*, Journal of Algebraic
           Combinatorics, volume 1 (1992), p. 111-132.

        EXAMPLES::

            sage: A = AlternatingSignMatrix(matrix([[0,0,1,0],[0,1,-1,1],[1,0,0,0],[0,0,1,0]]))
            sage: B = AlternatingSignMatrix(matrix([[0,0,1,0,0],[0,0,0,1,0],[1,0,0,-1,1],[0,1,0,0,0],[0,0,0,1,0]]))
            sage: A.ASM_compatible(B)
            True
            sage: A = AlternatingSignMatrix(matrix([[0,1,0],[1,-1,1],[0,1,0]]))
            sage: B = AlternatingSignMatrix(matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]))
            sage: A.ASM_compatible(B)
            False
        """
        if B.parent()._n - self.parent()._n != 1:
            raise ValueError("mismatched sizes")

        AA = self.corner_sum_matrix()
        BB = B.corner_sum_matrix()
        for i in range(0, len(AA[0])):
            for j in range(0, len(AA[0])):
                if not (AA[i,j]>=BB[i,j] and AA[i,j]>=BB[i+1,j+1]-1 \
                        and AA[i,j]<=BB[i+1,j] and AA[i,j]<=BB[i,j+1]):
                    return False
        return True

    def ASM_compatible_bigger(self):
        r"""
        Return all ASM's compatible with ``self`` that are of size one greater
        than ``self``.

        Given an `n \times n` alternating sign matrix `A`, there are as many
        ASM's of size `n+1` compatible with `A` as 2 raised to the power of
        the number of 1's in `A` [EKLP92]_.

        EXAMPLES::

            sage: A = AlternatingSignMatrix(matrix([[1,0],[0,1]]))
            sage: A.ASM_compatible_bigger()
            [
            [ 0  1  0]  [1 0 0]  [0 1 0]  [1 0 0]
            [ 1 -1  1]  [0 0 1]  [1 0 0]  [0 1 0]
            [ 0  1  0], [0 1 0], [0 0 1], [0 0 1]
            ]
            sage: B = AlternatingSignMatrix(matrix([[0,1],[1,0]]))
            sage: B.ASM_compatible_bigger()
            [
            [0 0 1]  [0 0 1]  [0 1 0]  [ 0  1  0]
            [0 1 0]  [1 0 0]  [0 0 1]  [ 1 -1  1]
            [1 0 0], [0 1 0], [1 0 0], [ 0  1  0]
            ]
        """
        n = self.parent()._n + 1
        M = AlternatingSignMatrices(n)
        sign = []
        asm = self.to_matrix()
        B = matrix(n+1)
        A = matrix([[2*(i+j-2*nw_corner_sum(asm,i,j))+1 for i in range(n)]
                    for j in range(n)])
        for a in range(n+1):
            B[a,0] = 2*a
            B[0,a] = 2*a
            B[a,n] = 2*(n-a)
            B[n,a] = 2*(n-a)

        for i in range(1,n):
            for j in range(1,n):
                if A[i-1,j-1] == A[i,j] == A[i-1,j]-2 == A[i,j-1]-2:
                    B[i,j] = -A[i,j]
                    sign.append([i,j])
                else:
                    B[i,j] = list({A[i-1,j-1]-1,A[i-1,j-1]+3} & {A[i-1,j]-3,A[i-1,j]+1} & {A[i,j-1]-3,A[i,j-1]+1} & {A[i,j]-1,A[i,j]+3})[0]

        output = [B]
        for b in range(len(sign)):
            N = len(output)
            for c in range(N):
                d = copy.copy(output[c])
                output[c][sign[b][0],sign[b][1]] = -output[c][sign[b][0], sign[b][1]] + 3
                d[sign[b][0],sign[b][1]] = -d[sign[b][0], sign[b][1]]-1
                output.append(d)

        for k in range(len(output)):
            output[k] = M.from_height_function(output[k]/2)
        return(output)

    def ASM_compatible_smaller(self):
        r"""
        Return the list of all ASMs compatible with ``self`` that are of size
        one smaller than ``self``.

        Given an alternating sign matrix `A` of size `n`, there are as many
        ASM's of size `n-1` compatible with it as 2 raised to the power of
        the number of `-1`'s in `A` [EKLP92]_.

        EXAMPLES::

            sage: A = AlternatingSignMatrix(matrix([[0,0,1,0],[0,1,-1,1],[1,0,0,0],[0,0,1,0]]))
            sage: A.ASM_compatible_smaller()
            [
            [0 0 1]  [ 0  1  0]
            [1 0 0]  [ 1 -1  1]
            [0 1 0], [ 0  1  0]
            ]
            sage: B = AlternatingSignMatrix(matrix([[1,0,0],[0,0,1],[0,1,0]]))
            sage: B.ASM_compatible_smaller()
            [
            [1 0]
            [0 1]
            ]

        """
        n = self.parent()._n
        M = AlternatingSignMatrices(n)
        A = matrix(n)
        asm = self.to_matrix()
        B = matrix([[2*(i+j-2*nw_corner_sum(asm,i,j)) for i in range(n)] for j in range(n)])
        sign = []
        for a in range(n):
            A[a,0] = 2*a + 1
            A[0,a] = 2*a + 1
            A[n-1,a] = 2*(n-a) - 1
            A[a,n-1] = 2*(n-a) - 1

        for i in range(n-1):
            for j in range(n-1):
                if B[i+1,j+1] == B[i,j] == B[i,j+1]+2 == B[i+1,j]+2:
                    A[i,j] = -B[i,j]
                    sign.append([i,j])
                else:
                    A[i,j] = list({B[i,j]+1,B[i,j]-3} & {B[i,j+1]+3,B[i,j+1]-1} & {B[i+1,j]+3,B[i+1,j]-1} & {B[i+1,j+1]+1,B[i+1,j+1]-3})[0]

        output = [A]
        for b in range(len(sign)):
            N = len(output)
            for c in range(N):
                d = copy.copy(output[c])
                output[c][sign[b][0],sign[b][1]] = -output[c][sign[b][0], sign[b][1]]+1
                d[sign[b][0],sign[b][1]] = -d[sign[b][0], sign[b][1]]-3
                output.append(d)
        for k in range(0,len(output)):
            output[k] = M.from_height_function((output[k]-matrix.ones(n,n))/2)
        return(output)

    @combinatorial_map(name='to Dyck word')
    @rename_keyword(deprecation=19572, method='algorithm')
    def to_dyck_word(self, algorithm):
        r"""
        Return a Dyck word determined by the specified algorithm.
 
        The algorithm 'last_diagonal' uses the last diagonal of the monotone 
        triangle corresponding to ``self``. The algorithm 'link_pattern' returns 
        the Dyck word in bijection with the link pattern of the fully packed 
        loop.

        Note that these two algorithms in general yield different Dyck words for a
        given alternating sign matrix.

        INPUT:

        - ``algorithm``  - 

          - ``'last_diagonal'`` 
          - ``'link_pattern'`` 

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[0,1,0],[1,0,0],[0,0,1]]).to_dyck_word(algorithm = 'last_diagonal')
            [1, 1, 0, 0, 1, 0]
            sage: d = A([[0,1,0],[1,-1,1],[0,1,0]]).to_dyck_word(algorithm = 'last_diagonal'); d
            [1, 1, 0, 1, 0, 0]
            sage: parent(d)
            Complete Dyck words
            sage: A = AlternatingSignMatrices(3)
            sage: asm = A([[0,1,0],[1,0,0],[0,0,1]])
            sage: asm.to_dyck_word(algorithm = 'link_pattern')
            [1, 0, 1, 0, 1, 0]
            sage: asm = A([[0,1,0],[1,-1,1],[0,1,0]])
            sage: asm.to_dyck_word(algorithm = 'link_pattern')
            [1, 0, 1, 1, 0, 0]
            sage: A = AlternatingSignMatrices(4)
            sage: asm = A([[0,0,1,0],[1,0,0,0],[0,1,-1,1],[0,0,1,0]])
            sage: asm.to_dyck_word(algorithm = 'link_pattern')
            [1, 1, 1, 0, 1, 0, 0, 0]
            sage: asm.to_dyck_word()
            Traceback (most recent call last):
            ...
            TypeError: to_dyck_word() takes exactly 2 arguments (1 given)
            sage: asm.to_dyck_word(algorithm = 'notamethod')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'notamethod'
        """
        if algorithm == 'last_diagonal':
            MT = self.to_monotone_triangle()
            nplus = self._matrix.nrows() + 1
            parkfn = [nplus - row[0] for row in list(MT) if len(row) > 0]
            return NonDecreasingParkingFunction(parkfn).to_dyck_word().reverse()
        
        elif algorithm == 'link_pattern':
            from sage.combinat.perfect_matching import PerfectMatching        
            from sage.combinat.dyck_word import DyckWords        
            p = PerfectMatching(self.link_pattern()).to_non_crossing_set_partition()
            asm = self.to_matrix()
            n = asm.nrows()
            d = DyckWords(n)
            return d.from_noncrossing_partition(p)

        raise ValueError("unknown algorithm '%s'" % algorithm)

    def number_negative_ones(self):
        """
        Return the number of entries in ``self`` equal to -1.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: asm = A([[0,1,0],[1,0,0],[0,0,1]])
            sage: asm.number_negative_ones()
            0
            sage: asm = A([[0,1,0],[1,-1,1],[0,1,0]])
            sage: asm.number_negative_ones()
            1
        """
        a = self._matrix
        return sum(1 for (i,j) in a.nonzero_positions() if a[i,j] == -1)

    def is_permutation(self):
        """
        Return ``True`` if ``self`` is a permutation matrix
        and ``False`` otherwise.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: asm = A([[0,1,0],[1,0,0],[0,0,1]])
            sage: asm.is_permutation()
            True
            sage: asm = A([[0,1,0],[1,-1,1],[0,1,0]])
            sage: asm.is_permutation()
            False
        """
        return self.number_negative_ones() == 0

    def to_permutation(self):
        """
        Return the corresponding permutation if ``self`` is a permutation
        matrix.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: asm = A([[0,1,0],[1,0,0],[0,0,1]])
            sage: p = asm.to_permutation(); p
            [2, 1, 3]
            sage: parent(p)
            Standard permutations
            sage: asm = A([[0,1,0],[1,-1,1],[0,1,0]])
            sage: asm.to_permutation()
            Traceback (most recent call last):
            ...
            ValueError: Not a permutation matrix
        """
        if not self.is_permutation():
            raise ValueError('Not a permutation matrix')
        asm_matrix = self.to_matrix()
        return Permutation([ j+1 for (i,j) in asm_matrix.nonzero_positions() ])

    @combinatorial_map(name='to semistandard tableau')
    def to_semistandard_tableau(self):
        """
        Return the semistandard tableau corresponding the monotone triangle
        corresponding to ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[0,0,1],[1,0,0],[0,1,0]]).to_semistandard_tableau()
            [[1, 1, 3], [2, 3], [3]]
            sage: t = A([[0,1,0],[1,-1,1],[0,1,0]]).to_semistandard_tableau(); t
            [[1, 1, 2], [2, 3], [3]]
            sage: parent(t)
            Semistandard tableaux
            """
        from sage.combinat.tableau import SemistandardTableau
        mt = self.to_monotone_triangle()
        ssyt = [[0]*(len(mt) - j) for j in range(len(mt))]
        for i in range(len(mt)):
            for j in range(len(mt[i])):
                ssyt[i][j] = mt[j][-(i+1)]
        return SemistandardTableau(ssyt)



    def left_key(self):
        r"""
        Return the left key of the alternating sign matrix ``self``.

        The left key of an alternating sign matrix was defined by Lascoux
        in [LascouxPreprint]_ and is obtained by successively removing all the
        `-1`'suntil what remains is a permutation matrix. This notion
        corresponds to the notion of left key for semistandard tableaux. So
        our algorithm proceeds as follows: we map ``self`` to its
        corresponding monotone triangle, view that monotone triangle as a
        semistandard tableaux, take its left key, and then map back through
        monotone triangles to the permutation matrix which is the left key.

        REFERENCES:

        .. [Aval07] \J.-C. Aval. *Keys and alternating sign matrices*.
           Sem. Lothar. Combin. 59 (2007/10), Art. B59f, 13 pp.

        .. [LascouxPreprint] \A. Lascoux. *Chern and Yang through ice*.
           Preprint.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[0,0,1],[1,0,0],[0,1,0]]).left_key()
            [0 0 1]
            [1 0 0]
            [0 1 0]
            sage: t = A([[0,1,0],[1,-1,1],[0,1,0]]).left_key(); t
            [1 0 0]
            [0 0 1]
            [0 1 0]
            sage: parent(t)
            Alternating sign matrices of size 3
        """
        lkey = self.to_semistandard_tableau().left_key_tableau()
        mt = [[0]*(len(lkey) - j) for j in range(len(lkey))]
        for i in range(len(lkey)):
            for j in range(len(lkey[i])):
                mt[i][j] = lkey[len(lkey[i])-j-1][i]
        A = AlternatingSignMatrices(len(lkey))
        return A.from_monotone_triangle(mt)

    @combinatorial_map(name='to left key permutation')
    def left_key_as_permutation(self):
        """
        Return the permutation of the left key of ``self``.

        .. SEEALSO::

            - :meth:`left_key()`

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A([[0,0,1],[1,0,0],[0,1,0]]).left_key_as_permutation()
            [3, 1, 2]
            sage: t = A([[0,1,0],[1,-1,1],[0,1,0]]).left_key_as_permutation(); t
            [1, 3, 2]
            sage: parent(t)
            Standard permutations
        """
        return self.left_key().to_permutation()

class AlternatingSignMatrices(UniqueRepresentation, Parent):
    r"""
    Class of all `n \times n` alternating sign matrices.

    An alternating sign matrix of size `n` is an `n \times n` matrix of `0`'s,
    `1`'s and `-1`'s such that the sum of each row and column is `1` and the
    non-zero entries in each row and column alternate in sign.

    Alternating sign matrices of size `n` are in bijection with
    :class:`monotone triangles <MonotoneTriangles>` with `n` rows.

    INPUT:

    - `n` -- an integer, the size of the matrices.

    - ``use_monotone_triangle`` -- deprecated

    EXAMPLES:

    This will create an instance to manipulate the alternating sign
    matrices of size 3::

        sage: A = AlternatingSignMatrices(3)
        sage: A
        Alternating sign matrices of size 3
        sage: A.cardinality()
        7

    Notably, this implementation allows to make a lattice of it::

        sage: L = A.lattice()
        sage: L
        Finite lattice containing 7 elements
        sage: L.category()
        Join of Category of finite lattice posets
         and Category of finite enumerated sets
         and Category of facade sets
    """
    def __init__(self, n, use_monotone_triangles=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: A = AlternatingSignMatrices(4)
            sage: TestSuite(A).run()
        """
        self._n = n
        self._matrix_space = MatrixSpace(ZZ, n)
        if use_monotone_triangles is not None:
            from sage.misc.superseded import deprecation
            deprecation(18208, 'use_monotone_triangles is deprecated')
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: A = AlternatingSignMatrices(4); A
            Alternating sign matrices of size 4
        """
        return "Alternating sign matrices of size %s" % self._n

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A._repr_option('element_ascii_art')
            True
        """
        return self._matrix_space._repr_option(key)

    def __contains__(self, asm):
        """
        Check if ``asm`` is in ``self``.

        TESTS::

            sage: A = AlternatingSignMatrices(3)
            sage: [[0,1,0],[1,0,0],[0,0,1]] in A
            True
            sage: [[0,1,0],[1,-1,1],[0,1,0]] in A
            True
            sage: [[0, 1],[1,0]] in A
            False
            sage: [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]] in A
            False
            sage: [[-1, 1, 1],[1,-1,1],[1,1,-1]] in A
            False
        """
        if isinstance(asm, AlternatingSignMatrix):
            return asm._matrix.nrows() == self._n
        try:
            asm = self._matrix_space(asm)
        except (TypeError, ValueError):
            return False
        for row in asm:
            pos = False
            for val in row:
                if val > 0:
                    if pos:
                        return False
                    else:
                        pos = True
                elif val < 0:
                    if pos:
                        pos = False
                    else:
                        return False
            if not pos:
                return False
        if any(sum(row) != 1 for row in asm.columns()):
            return False
        return True

    def _element_constructor_(self, asm):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: elt = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]]); elt
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: elt.parent() is A
            True
            sage: A([[3, 2, 1], [2, 1], [1]])
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if isinstance(asm, AlternatingSignMatrix):
            if asm.parent() is self:
                return asm
            raise ValueError("Cannot convert between alternating sign matrices of different sizes")
        if asm in MonotoneTriangles(self._n):
            return self.from_monotone_triangle(asm)
        m = self._matrix_space(asm)
        m.set_immutable()
        return self.element_class(self, m)

    Element = AlternatingSignMatrix

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A._an_element_()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return self.element_class(self, self._matrix_space.identity_matrix())

    def random_element(self):
        r"""
        Return a uniformly random alternating sign matrix.

        EXAMPLES::

            sage: AlternatingSignMatrices(7).random_element()  # random
            [ 0  0  0  0  1  0  0]
            [ 0  0  1  0 -1  0  1]
            [ 0  0  0  0  1  0  0]
            [ 0  1 -1  0  0  1  0]
            [ 1 -1  1  0  0  0  0]
            [ 0  0  0  1  0  0  0]
            [ 0  1  0  0  0  0  0]
            sage: a = AlternatingSignMatrices(5).random_element()
            sage: bool(a.number_negative_ones()) or a.is_permutation()
            True

        This is done using a modified version of Propp and Wilson's "coupling
        from the past" algorithm. It creates a uniformly random Gelfand-Tsetlin
        triangle with top row `[n, n-1, \ldots 2, 1]`, and then converts it to
        an alternating sign matrix.
        """
        from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPatterns
        n = self._n
        toprow = [n-i for i in range(n)]
        gt = GelfandTsetlinPatterns(top_row = toprow, strict = True)
        randomgt = gt.random_element()
        A = AlternatingSignMatrices(n)
        return A.from_monotone_triangle(randomgt)

    def from_monotone_triangle(self, triangle):
        r"""
        Return an alternating sign matrix from a monotone triangle.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A.from_monotone_triangle([[3, 2, 1], [2, 1], [1]])
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: A.from_monotone_triangle([[3, 2, 1], [3, 2], [3]])
            [0 0 1]
            [0 1 0]
            [1 0 0]
        """
        n = len(triangle)
        if n != self._n:
            raise ValueError("Incorrect size")
        asm = []

        prev = [0]*n
        for line in reversed(triangle):
            v = [1 if j+1 in reversed(line) else 0 for j in range(n)]
            row = [a-b for (a, b) in zip(v, prev)]
            asm.append(row)
            prev = v

        m = self._matrix_space(asm)
        m.set_immutable()
        return self.element_class(self, m)

    def from_corner_sum(self, corner):
        r"""
        Return an alternating sign matrix from a corner sum matrix.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A.from_corner_sum(matrix([[0,0,0,0],[0,1,1,1],[0,1,2,2],[0,1,2,3]]))
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: A.from_corner_sum(matrix([[0,0,0,0],[0,0,1,1],[0,1,1,2],[0,1,2,3]]))
            [ 0  1  0]
            [ 1 -1  1]
            [ 0  1  0]
        """
        asm_list=[]
        n = len(list(corner)) - 1
        for k in range(n):
            asm_list.append([])
        for i in range(n):
            for j in range(n):
                y = corner[i+1][j+1] \
                     - sum([sum([asm_list[i2][j2] for i2 in range(i)])
                            for j2 in range(j)]) \
                     - sum([asm_list[i2][j] for i2 in range(i)]) \
                     - sum([asm_list[i][j2] for j2 in range(j)])
                asm_list[i].append(y)
        return AlternatingSignMatrix(asm_list)

    def from_height_function(self,height):
        r"""
        Return an alternating sign matrix from a height function.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A.from_height_function(matrix([[0,1,2,3],[1,2,1,2],[2,3,2,1],[3,2,1,0]]))
            [0 0 1]
            [1 0 0]
            [0 1 0]
            sage: A.from_height_function(matrix([[0,1,2,3],[1,2,1,2],[2,1,2,1],[3,2,1,0]]))
            [ 0  1  0]
            [ 1 -1  1]
            [ 0  1  0]
        """
        return self.from_corner_sum(matrix( [[((i+j-height[i][j])//2)
                                              for i in range(len(list(height)))]
                                             for j in range(len(list(height)))] ))

    def from_contre_tableau(self, comps):
        r"""
        Return an alternating sign matrix from a contre-tableau.

        EXAMPLES::

            sage: ASM = AlternatingSignMatrices(3)
            sage: ASM.from_contre_tableau([[1, 2, 3], [1, 2], [1]])
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: ASM.from_contre_tableau([[1, 2, 3], [2, 3], [3]])
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        n = len(comps)
        M = [ [0 for _ in range(n)] for _ in range(n) ]

        previous_set = set([])
        for col in range(n-1, -1, -1):
            s = set( comps[col] )
            for x in s.difference(previous_set):
                M[x-1][col] = 1
            for x in previous_set.difference(s):
                M[x-1][col] = -1

            previous_set = s

        return AlternatingSignMatrix(M)

    def size(self):
        r"""
        Return the size of the matrices in ``self``.

        TESTS::

            sage: A = AlternatingSignMatrices(4)
            sage: A.size()
            4
        """
        return self._n

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of `n \times n` alternating sign matrices is equal to

        .. MATH::

            \prod_{k=0}^{n-1} \frac{(3k+1)!}{(n+k)!} = \frac{1! 4! 7! 10!
            \cdots (3n-2)!}{n! (n+1)! (n+2)! (n+3)! \cdots (2n-1)!}

        EXAMPLES::

            sage: [AlternatingSignMatrices(n).cardinality() for n in range(0, 11)]
            [1, 1, 2, 7, 42, 429, 7436, 218348, 10850216, 911835460, 129534272700]
        """
        return Integer(prod( [ factorial(3*k+1)/factorial(self._n+k)
                       for k in range(self._n)] ))

    def matrix_space(self):
        """
        Return the underlying matrix space.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: A.matrix_space()
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
        """
        return self._matrix_space

    def __iter__(self):
        r"""
        Iterator on the alternating sign matrices of size `n`.

        TESTS::

            sage: A = AlternatingSignMatrices(4)
            sage: len(list(A))
            42
        """
        for t in MonotoneTriangles(self._n):
            yield self.from_monotone_triangle(t)

    def _lattice_initializer(self):
        r"""
        Return a 2-tuple to use in argument of ``LatticePoset``.

        For more details about the cover relations, see
        ``MonotoneTriangles``. Notice that the returned matrices are
        made immutable to ensure their hashability required by
        ``LatticePoset``.

        EXAMPLES:

        Proof of the lattice property for alternating sign matrices of
        size 3::

            sage: A = AlternatingSignMatrices(3)
            sage: P = Poset(A._lattice_initializer())
            sage: P.is_lattice()
            True
        """
        (mts, rels) = MonotoneTriangles(self._n)._lattice_initializer()
        bij = dict((t, self.from_monotone_triangle(t)) for t in mts)
        asms, rels = bij.itervalues(), [(bij[a], bij[b]) for (a,b) in rels]
        return (asms, rels)

    def cover_relations(self):
        r"""
        Iterate on the cover relations between the alternating sign
        matrices.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: for (a,b) in A.cover_relations():
            ....:   eval('a, b')
            (
            [1 0 0]  [0 1 0]
            [0 1 0]  [1 0 0]
            [0 0 1], [0 0 1]
            )
            (
            [1 0 0]  [1 0 0]
            [0 1 0]  [0 0 1]
            [0 0 1], [0 1 0]
            )
            (
            [0 1 0]  [ 0  1  0]
            [1 0 0]  [ 1 -1  1]
            [0 0 1], [ 0  1  0]
            )
            (
            [1 0 0]  [ 0  1  0]
            [0 0 1]  [ 1 -1  1]
            [0 1 0], [ 0  1  0]
            )
            (
            [ 0  1  0]  [0 0 1]
            [ 1 -1  1]  [1 0 0]
            [ 0  1  0], [0 1 0]
            )
            (
            [ 0  1  0]  [0 1 0]
            [ 1 -1  1]  [0 0 1]
            [ 0  1  0], [1 0 0]
            )
            (
            [0 0 1]  [0 0 1]
            [1 0 0]  [0 1 0]
            [0 1 0], [1 0 0]
            )
            (
            [0 1 0]  [0 0 1]
            [0 0 1]  [0 1 0]
            [1 0 0], [1 0 0]
            )

        """
        (_, rels) = self._lattice_initializer()
        return (_ for _ in rels)

    def lattice(self):
        r"""
        Return the lattice of the alternating sign matrices of size
        `n`, created by ``LatticePoset``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: L = A.lattice()
            sage: L
            Finite lattice containing 7 elements

        """
        return LatticePoset(self._lattice_initializer(), cover_relations=True)

    @cached_method
    def gyration_orbits(self):
        r"""
        Return the list of gyration orbits of ``self``.

        EXAMPLES::

            sage: AlternatingSignMatrices(3).gyration_orbits()
            ((
              [1 0 0]  [0 0 1]  [ 0  1  0]
              [0 1 0]  [0 1 0]  [ 1 -1  1]
              [0 0 1], [1 0 0], [ 0  1  0]
             ),
             (
              [0 1 0]  [1 0 0]
              [1 0 0]  [0 0 1]
              [0 0 1], [0 1 0]
             ),
             (
              [0 0 1]  [0 1 0]
              [1 0 0]  [0 0 1]
              [0 1 0], [1 0 0]
             ))
        """
        ASMs = list(self)
        perm = Permutation([ASMs.index(asm.gyration())+1 for asm in ASMs])
        return tuple([tuple([ASMs[i-1] for i in cyc])
                      for cyc in perm.cycle_tuples()])

    def gyration_orbit_sizes(self):
        r"""
        Return the sizes of gyration orbits of ``self``.

        EXAMPLES::

            sage: AlternatingSignMatrices(3).gyration_orbit_sizes()
            [3, 2, 2]
            sage: AlternatingSignMatrices(4).gyration_orbit_sizes()
            [4, 8, 2, 8, 8, 8, 2, 2]

            sage: A = AlternatingSignMatrices(5)
            sage: li = [5,10,10,10,10,10,2,5,10,10,10,10,10,10,10,10,10,10,10,10,\
            4,10,10,10,10,10,10,4,5,10,10,10,10,10,10,10,2,4,5,10,10,10,10,10,10,\
            4,5,10,10,2,2]
            sage: A.gyration_orbit_sizes() == li
            True
        """
        return [len(orbit) for orbit in self.gyration_orbits()]

class MonotoneTriangles(GelfandTsetlinPatternsTopRow):
    r"""
    Monotone triangles with `n` rows.

    A monotone triangle is a number triangle `(a_{i,j})_{1 \leq i \leq
    n , 1 \leq j \leq i}` on `\{1, \dots, n\}` such that:

    - `a_{i,j} < a_{i,j+1}`

    - `a_{i+1,j} < a_{i,j} \leq a_{i+1,j+1}`

    This notably requires that the bottom column is ``[1,...,n]``.

    Alternatively a monotone triangle is a strict Gelfand-Tsetlin pattern with
    top row `(n, \ldots, 2, 1)`.

    INPUT:

    - ``n`` -- The number of rows in the monotone triangles

    EXAMPLES:

    This represents the monotone triangles with base ``[3,2,1]``::

        sage: M = MonotoneTriangles(3)
        sage: M
        Monotone triangles with 3 rows
        sage: M.cardinality()
        7

    The monotone triangles are a lattice::

        sage: M.lattice()
        Finite lattice containing 7 elements

    Monotone triangles can be converted to alternating sign matrices
    and back::

        sage: M = MonotoneTriangles(5)
        sage: A = AlternatingSignMatrices(5)
        sage: all(A.from_monotone_triangle(m).to_monotone_triangle() == m for m in M)
        True
    """
    def __init__(self, n):
        r"""
        Initialize ``self``.

        TESTS::

            sage: M = MonotoneTriangles(4)
            sage: TestSuite(M).run()
            sage: M2 = MonotoneTriangles(int(4))
            sage: M is M2
            True
        """
        GelfandTsetlinPatternsTopRow.__init__(self, tuple(reversed(range(1, n+1))), True)

    def _repr_(self):
        r"""
        String representation.

        TESTS::

            sage: M = MonotoneTriangles(4)
            sage: M
            Monotone triangles with 4 rows
        """
        return "Monotone triangles with %s rows" % self._n

    def cardinality(self):
        r"""
        Cardinality of ``self``.

        The number of monotone triangles with `n` rows is equal to

        .. MATH::

            \prod_{k=0}^{n-1} \frac{(3k+1)!}{(n+k)!} = \frac{1! 4! 7! 10!
            \cdots (3n-2)!}{n! (n+1)! (n+2)! (n+3)! \cdots (2n-1)!}

        EXAMPLES::

            sage: M = MonotoneTriangles(4)
            sage: M.cardinality()
            42
        """
        return Integer(prod( [ factorial(3*k+1)/factorial(self._n+k)
                       for k in range(self._n)] ))

    def _lattice_initializer(self):
        r"""
        Return a 2-tuple to use in argument of ``LatticePoset``.

        This couple is composed by the set of the monotone triangles
        with `n` rows and the cover relations. Specializing this
        function allows to generate the monotone triangles just once,
        and so to speed up the computation in comparison of
        ``(list(self), self.cover_relations())``. Notice that the
        function also switch the representation of monotone triangles
        from list of list to tuple of tuple in order to make them
        hashable (required to make a poset with them).

        EXAMPLES::

            sage: M = MonotoneTriangles(3)
            sage: P = Poset(M._lattice_initializer())
            sage: P.is_lattice()
            True
        """
        # get a list of the elements and switch to a tuple
        # representation
        set_ = [tuple(tuple(_) for _ in x) for x in list(self)]
        return (set_, [(a,b) for a in set_ for b in set_ if _is_a_cover(a,b)])

    def cover_relations(self):
        r"""
        Iterate on the cover relations in the set of monotone triangles
        with `n` rows.

        EXAMPLES::

            sage: M = MonotoneTriangles(3)
            sage: for (a,b) in M.cover_relations():
            ....:   eval('a, b')
            ([[3, 2, 1], [2, 1], [1]], [[3, 2, 1], [2, 1], [2]])
            ([[3, 2, 1], [2, 1], [1]], [[3, 2, 1], [3, 1], [1]])
            ([[3, 2, 1], [2, 1], [2]], [[3, 2, 1], [3, 1], [2]])
            ([[3, 2, 1], [3, 1], [1]], [[3, 2, 1], [3, 1], [2]])
            ([[3, 2, 1], [3, 1], [2]], [[3, 2, 1], [3, 1], [3]])
            ([[3, 2, 1], [3, 1], [2]], [[3, 2, 1], [3, 2], [2]])
            ([[3, 2, 1], [3, 1], [3]], [[3, 2, 1], [3, 2], [3]])
            ([[3, 2, 1], [3, 2], [2]], [[3, 2, 1], [3, 2], [3]])
        """
        set_ = list(self)
        return ((a,b) for a in set_ for b in set_ if _is_a_cover(a,b))

    def lattice(self):
        r"""
        Return the lattice of the monotone triangles with `n` rows.

        EXAMPLES::

            sage: M = MonotoneTriangles(3)
            sage: P = M.lattice()
            sage: P
            Finite lattice containing 7 elements

        """
        return LatticePoset(self._lattice_initializer(), cover_relations=True)

def _is_a_cover(mt0, mt1):
    r"""
    Define the cover relations.

    Return ``True`` if and only if the second argument is a cover of
    the first one.

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm._is_a_cover([[1,2,3],[1,2],[1]], [[1,2,3],[1,3],[1]])
        True
        sage: asm._is_a_cover([[1,2,3],[1,3],[2]], [[1,2,3],[1,2],[1]])
        False
    """
    diffs = 0
    for (a,b) in itertools.izip(flatten(mt0), flatten(mt1)):
        if a != b:
            if a+1 == b:
                diffs += 1
            else:
                return False
        if diffs > 1:
            return False
    return diffs == 1

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.alternating_sign_matrix', 'AlternatingSignMatrices_n', AlternatingSignMatrices)
register_unpickle_override('sage.combinat.alternating_sign_matrix', 'MonotoneTriangles_n', MonotoneTriangles)

class ContreTableaux(Parent):
    """
    Factory class for the combinatorial class of contre tableaux of size `n`.

    EXAMPLES::

        sage: ct4 = ContreTableaux(4); ct4
        Contre tableaux of size 4
        sage: ct4.cardinality()
        42
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, n, **kwds):
        r"""
        Factory pattern.

        Check properties on arguments, then call the appropriate class.

        EXAMPLES::

            sage: C = ContreTableaux(4)
            sage: type(C)
            <class 'sage.combinat.alternating_sign_matrix.ContreTableaux_n'>

        """
        assert(isinstance(n, (int, Integer)))
        return ContreTableaux_n(n, **kwds)


class ContreTableaux_n(ContreTableaux):
    def __init__(self, n):
        """
        TESTS::

            sage: ct2 = ContreTableaux(2); ct2
            Contre tableaux of size 2
            sage: ct2 == loads(dumps(ct2))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS::

            sage: repr(ContreTableaux(2))
            'Contre tableaux of size 2'
        """
        return "Contre tableaux of size %s"%self.n

    def __eq__(self, other):
        """
        TESTS::

            sage: C = ContreTableaux(4)
            sage: C == loads(dumps(C))
            True

        """
        return self.n == other.n

    def cardinality(self):
        """
        EXAMPLES::

            sage: [ ContreTableaux(n).cardinality() for n in range(0, 11)]
            [1, 1, 2, 7, 42, 429, 7436, 218348, 10850216, 911835460, 129534272700]
        """
        return prod( [ factorial(3*k+1)/factorial(self.n+k) for k in range(self.n)] )

    def _iterator_rec(self, i):
        """
        EXAMPLES::

            sage: c = ContreTableaux(2)
            sage: list(c._iterator_rec(0))
            [[]]
            sage: list(c._iterator_rec(1))
            [[[1, 2]]]
            sage: list(c._iterator_rec(2))
            [[[1, 2], [1]], [[1, 2], [2]]]
        """
        if i == 0:
            yield []
        elif i == 1:
            yield [range(1, self.n+1)]
        else:
            for columns in self._iterator_rec(i-1):
                previous_column = columns[-1]
                for column in _next_column_iterator(previous_column, len(previous_column)-1):
                    yield columns + [ column ]

    def __iter__(self):
        """
        EXAMPLES::

            sage: list(ContreTableaux(0))
            [[]]
            sage: list(ContreTableaux(1))
            [[[1]]]
            sage: list(ContreTableaux(2))
            [[[1, 2], [1]], [[1, 2], [2]]]
            sage: list(ContreTableaux(3))
            [[[1, 2, 3], [1, 2], [1]],
             [[1, 2, 3], [1, 2], [2]],
             [[1, 2, 3], [1, 3], [1]],
             [[1, 2, 3], [1, 3], [2]],
             [[1, 2, 3], [1, 3], [3]],
             [[1, 2, 3], [2, 3], [2]],
             [[1, 2, 3], [2, 3], [3]]]
        """
        for z in self._iterator_rec(self.n):
            yield z


def _next_column_iterator(previous_column, height, i = None):
    """
    Return a generator for all columns of height ``height``
    properly filled from row 1 to ``i``.
    ("Properly filled" means strictly increasing and having
    the property that the `k`-th entry is `\geq` to the `k`-th
    entry of ``previous_column`` for each `k`.)

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: list(asm._next_column_iterator([1], 0))
        [[]]
        sage: list(asm._next_column_iterator([1,5],1))
        [[1], [2], [3], [4], [5]]
        sage: list(asm._next_column_iterator([1,4,5],2))
        [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5], [4, 5]]
    """
    if i is None:
        i = height
    if i == 0:
        yield [-1]*height
    else:
        for column in _next_column_iterator(previous_column, height, i-1):
            min_value = previous_column[i-1]
            if i > 1:
                min_value = max(min_value, column[i-2]+1)
            for value in range(min_value, previous_column[i]+1):
                c = column[:]
                c[i-1] = value
                yield c


def _previous_column_iterator(column, height, max_value):
    """
    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: list(asm._previous_column_iterator([2,3], 3, 4))
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
    """
    new_column = [1] + column + [ max_value ] * (height - len(column))
    return _next_column_iterator(new_column, height)


class TruncatedStaircases(Parent):
    """
    Factory class for the combinatorial class of truncated staircases
    of size ``n`` with last column ``last_column``.

    EXAMPLES::

        sage: t4 = TruncatedStaircases(4, [2,3]); t4
        Truncated staircases of size 4 with last column [2, 3]
        sage: t4.cardinality()
        4
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, n, last_column, **kwds):
        r"""
        Factory pattern.

        Check properties on arguments, then call the appropriate class.

        TESTS::

            sage: T = TruncatedStaircases(4, [2,3])
            sage: type(T)
            <class 'sage.combinat.alternating_sign_matrix.TruncatedStaircases_nlastcolumn'>

        """
        assert(isinstance(n, (int, Integer)))
        return TruncatedStaircases_nlastcolumn(n, last_column, **kwds)


class TruncatedStaircases_nlastcolumn(TruncatedStaircases):
    def __init__(self, n, last_column):
        """
        TESTS::

            sage: t4 = TruncatedStaircases(4, [2,3]); t4
            Truncated staircases of size 4 with last column [2, 3]
            sage: t4 == loads(dumps(t4))
            True
        """
        self.n = n
        self.last_column = last_column

    def __repr__(self):
        """
        TESTS::

            sage: repr(TruncatedStaircases(4, [2,3]))
            'Truncated staircases of size 4 with last column [2, 3]'
        """
        return "Truncated staircases of size %s with last column %s"%(self.n, self.last_column)

    def _iterator_rec(self, i):
        """
        EXAMPLES::

            sage: t = TruncatedStaircases(3, [2,3])
            sage: list(t._iterator_rec(1))
            []
            sage: list(t._iterator_rec(2))
            [[[2, 3]]]
            sage: list(t._iterator_rec(3))
            [[[1, 2, 3], [2, 3]]]
        """
        if i < len(self.last_column):
            return
        elif i == len(self.last_column):
            yield [self.last_column]
        else:
            for columns in self._iterator_rec(i-1):
                previous_column = columns[0]
                for column in _previous_column_iterator(previous_column, len(previous_column)+1, self.n):
                    yield [column] + columns

    def __iter__(self):
        """
        EXAMPLES::

            sage: list(TruncatedStaircases(4, [2,3]))
            [[[4, 3, 2, 1], [3, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 2], [3, 2]]]
        """
        for z in self._iterator_rec(self.n):
            yield [list(reversed(x)) for x in z]

    def __eq__(self, other):
        r"""
        TESTS::

            sage: T = TruncatedStaircases(4, [2,3])
            sage: T == loads(dumps(T))
            True
        """
        return ((self.n == other.n) and
                (self.last_column == other.last_column))

    def cardinality(self):
        r"""
        EXAMPLES::

            sage: T = TruncatedStaircases(4, [2,3])
            sage: T.cardinality()
            4
        """
        c = 0
        for _ in self:
            c += 1
        return c

def nw_corner_sum(M,i,j):
    r"""
    Return the sum of entries to the northwest of `(i,j)` in matrix.

    EXAMPLES::

        sage: from sage.combinat.alternating_sign_matrix import nw_corner_sum
        sage: A = matrix.ones(3,3)
        sage: nw_corner_sum(A,2,2)
        4
    """
    if i >= 0 and j >= 0:
        return sum([sum(M[i2][:j]) for i2 in range(i)])
    return 0

