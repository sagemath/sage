r"""
Gelfand-Tsetlin Patterns

AUTHORS:

- Travis Scrimshaw (2013-15-03): Initial version

REFERENCES:

.. [BBF] \B. Brubaker, D. Bump, and S. Friedberg.
   Weyl Group Multiple Dirichlet Series: Type A Combinatorial Theory.
   Ann. of Math. Stud., vol. 175, Princeton Univ. Press, New Jersey, 2011.

.. [GC50] \I. M. Gelfand and M. L. Cetlin.
   Finite-Dimensional Representations of the Group of Unimodular Matrices.
   Dokl. Akad. Nauk SSSR **71**, pp. 825--828, 1950.

.. [Tok88] \T. Tokuyama.
   A Generating Function of Strict Gelfand Patterns and Some Formulas on
   Characters of General Linear Groups.
   J. Math. Soc. Japan **40** (4), pp. 671--685, 1988.

"""
# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
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
from __future__ import annotations

from sage.structure.parent import Parent
from sage.structure.list_clone import ClonableArray
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.combinat.partition import Partitions
from sage.combinat.tableau import Tableau, SemistandardTableaux
from sage.combinat.combinatorial_map import combinatorial_map
from sage.misc.misc_c import prod


class GelfandTsetlinPattern(ClonableArray,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A Gelfand-Tsetlin (sometimes written as Gelfand-Zetlin or Gelfand-Cetlin)
    pattern.  They were originally defined in [GC50]_.

    A Gelfand-Tsetlin pattern is a triangular array:

    .. MATH::

        \begin{array}{ccccccccc}
        a_{1,1} & & a_{1,2} & & a_{1,3} & & \cdots & & a_{1,n} \\
        & a_{2,2} & & a_{2,3} & & \cdots & & a_{2,n} \\
        & & a_{3,3} & &  \cdots & & a_{3,n} \\
        & & & \ddots \\
        & & & & a_{n,n}
        \end{array}

    such that `a_{i,j} \geq a_{i+1,j+1} \geq a_{i,j+1}`.

    Gelfand-Tsetlin patterns are in bijection with semistandard Young tableaux
    by the following algorithm. Let `G` be a Gelfand-Tsetlin pattern with
    `\lambda^{(k)}` being the `(n-k+1)`-st row (note that this is a partition).
    The definition of `G` implies

    .. MATH::

        \lambda^{(0)} \subseteq \lambda^{(1)} \subseteq \cdots \subseteq
        \lambda^{(n)},

    where `\lambda^{(0)}` is the empty partition, and each skew shape
    `\lambda^{(k)}/\lambda^{(k-1)}` is a horizontal strip. Thus define `T(G)`
    by inserting `k` into the squares of the skew shape
    `\lambda^{(k)}/ \lambda^{(k-1)}`, for `k=1,\dots,n`.

    To each entry in a Gelfand-Tsetlin pattern, one may attach a decoration of
    a circle or a box (or both or neither).  These decorations appear in the
    study of Weyl group multiple Dirichlet series, and are implemented here
    following the exposition in [BBF]_.

    .. NOTE::

        We use the "right-hand" rule for determining circled and boxed entries.

    .. WARNING::

        The entries in Sage are 0-based and are thought of as flushed to the
        left in a matrix. In other words, the coordinates of entries in the
        Gelfand-Tsetlin patterns are thought of as the matrix:

        .. MATH::

            \begin{bmatrix}
            g_{0,0} & g_{0,1} & g_{0,2} & \cdots & g_{0,n-2} & g_{n-1,n-1} \\
            g_{1,0} & g_{1,1} & g_{1,2} & \cdots & g_{1,n-2} \\
            g_{2,0} & g_{2,1} & g_{2,2} & \cdots \\
            \vdots & \vdots & \vdots \\
            g_{n-2,0} & g_{n-2,1} \\
            g_{n-1,0}
            \end{bmatrix}.

        However, in the discussions, we will be using the **standard**
        numbering system.

    EXAMPLES::

        sage: G = GelfandTsetlinPattern([[3, 2, 1], [2, 1], [1]]); G
        [[3, 2, 1], [2, 1], [1]]
        sage: G.pp()
          3     2     1
             2     1
                1
        sage: G = GelfandTsetlinPattern([[7, 7, 4, 0], [7, 7, 3], [7, 5], [5]]); G.pp()
          7     7     4     0
             7     7     3
                7     5
                   5
        sage: G.to_tableau().pp()
          1  1  1  1  1  2  2
          2  2  2  2  2  3  3
          3  3  3  4
    """
    # Note that the width == height, so len(gt) == len(gt[0]) except
    #   we don't have to check if it is the entry GT pattern
    @staticmethod
    def __classcall_private__(self, gt):
        """
        Return ``gt`` as a proper element of :class:`GelfandTsetlinPatterns`.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[3,2,1],[2,1],[1]])
            sage: G.parent()
            Gelfand-Tsetlin patterns
            sage: TestSuite(G).run()
        """
        return GelfandTsetlinPatterns()(gt)

    def check(self):
        """
        Check that this is a valid Gelfand-Tsetlin pattern.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: G([[3,2,1],[2,1],[1]]).check()
        """
        assert all(self[i - 1][j] >= self[i][j] >= self[i - 1][j + 1]
                   for i in range(1, len(self)) for j in range(len(self[i])))

    def _hash_(self) -> int:
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: gt = G([[3,2,1],[2,1],[1]])
            sage: hash(gt) == hash(gt)
            True

        Check that :trac:`14717` is fixed::

            sage: GT = GelfandTsetlinPattern([[2, 1, 0], [2, 0], [1]])
            sage: GT in {}
            False
        """
        return hash(tuple(map(tuple, self)))

    def _repr_diagram(self) -> str:
        """
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: print(G([[3,2,1],[2,1],[1]])._repr_diagram())
              3     2     1
                 2     1
                    1
        """
        ret = ''
        for i, row in enumerate(self):
            if i != 0:
                ret += '\n'
            ret += '   '*i
            ret += '   '.join('%3s' % val for val in row)
        return ret

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: G([[3,2,1],[2,1],[1]]).pp()
              3     2     1
                 2     1
                    1
        """
        print(self._repr_diagram())

    def _latex_(self) -> str:
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: latex(G([[3,2,1],[2,1],[1]]))
            \begin{array}{ccccc}
            3 & & 2 & & 1 \\
            & 2 & & 1 & \\
            & & 1 & &
            \end{array}
            sage: latex(G([]))
            \emptyset
        """
        n = len(self)
        if n == 0:
            return "\\emptyset"
        ret = "\\begin{array}{" + 'c'*(n*2-1) + "}\n"
        for i, row in enumerate(self):
            if i > 0:
                ret += " \\\\\n"
            ret += "& "*i
            ret += " & & ".join(repr(val) for val in row)
            ret += " &"*i
        return ret + "\n\\end{array}"

    @combinatorial_map(name='to semistandard tableau')
    def to_tableau(self):
        r"""
        Return ``self`` as a semistandard Young tableau.

        The conversion from a Gelfand-Tsetlin pattern to a semistandard Young
        tableaux is as follows. Let `G` be a Gelfand-Tsetlin pattern with
        `\lambda^{(k)}` being the `(n-k+1)`-st row (note that this is a
        partition).  The definition of `G` implies

        .. MATH::

            \lambda^{(0)} \subseteq \lambda^{(1)} \subseteq \cdots \subseteq
            \lambda^{(n)},

        where `\lambda^{(0)}` is the empty partition, and each skew shape
        `\lambda^{(k)} / \lambda^{(k-1)}` is a horizontal strip. Thus define
        `T(G)` by inserting `k` into the squares of the skew shape
        `\lambda^{(k)} / \lambda^{(k-1)}`, for `k=1,\dots,n`.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: elt = G([[3,2,1],[2,1],[1]])
            sage: T = elt.to_tableau(); T
            [[1, 2, 3], [2, 3], [3]]
            sage: T.pp()
              1  2  3
              2  3
              3
            sage: G(T) == elt
            True
        """
        ret = []
        for i, row in enumerate(reversed(self)):
            for j, val in enumerate(row):
                if j >= len(ret):
                    if val == 0:
                        break
                    ret.append([i+1]*val)
                else:
                    ret[j].extend([i+1]*(val-len(ret[j])))
        S = SemistandardTableaux(max_entry=len(self))
        return S(ret)

    @cached_method
    def boxed_entries(self) -> tuple:
        """
        Return the position of the boxed entries of ``self``.

        Using the *right-hand* rule, an entry `a_{i,j}` is boxed if
        `a_{i,j} = a_{i-1,j-1}`; i.e., `a_{i,j}` has the same value as its
        neighbor to the northwest.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[3,2,1],[3,1],[1]])
            sage: G.boxed_entries()
            ((1, 0),)
        """
        ret = []
        for i in range(1, len(self)):
            for j in range(len(self[i])):
                if self[i][j] == self[i-1][j]:
                    ret.append((i, j))
        return tuple(ret)

    @cached_method
    def circled_entries(self) -> tuple:
        """
        Return the circled entries of ``self``.

        Using the *right-hand* rule, an entry `a_{i,j}` is circled if
        `a_{i,j} = a_{i-1,j}`; i.e., `a_{i,j}` has the same value as its
        neighbor to the northeast.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[3,2,1],[3,1],[1]])
            sage: G.circled_entries()
            ((1, 1), (2, 0))
        """
        ret = []
        for i in range(1, len(self)):
            for j in range(len(self[i])):
                if self[i][j] == self[i-1][j+1]:
                    ret.append((i, j))
        return tuple(ret)

    @cached_method
    def special_entries(self) -> tuple:
        """
        Return the special entries.

        An entry `a_{i,j}` is special if `a_{i-1,j-1} > a_{i,j} > a_{i-1,j}`,
        that is to say, the entry is neither boxed nor circled and is **not**
        in the first row. The name was coined by [Tok88]_.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[3,2,1],[3,1],[1]])
            sage: G.special_entries()
            ()
            sage: G = GelfandTsetlinPattern([[4,2,1],[4,1],[2]])
            sage: G.special_entries()
            ((2, 0),)
        """
        ret = []
        for i in range(1, len(self)):
            for j in range(len(self[i])):
                if self[i-1][j] > self[i][j] and self[i][j] > self[i-1][j+1]:
                    ret.append((i, j))
        return tuple(ret)

    def number_of_boxes(self) -> int:
        """
        Return the number of boxed entries. See :meth:`boxed_entries()`.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[3,2,1],[3,1],[1]])
            sage: G.number_of_boxes()
            1
        """
        return len(self.boxed_entries())

    def number_of_circles(self) -> int:
        """
        Return the number of boxed entries. See :meth:`circled_entries()`.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[3,2,1],[3,1],[1]])
            sage: G.number_of_circles()
            2
        """
        return len(self.circled_entries())

    def number_of_special_entries(self) -> int:
        """
        Return the number of special entries. See :meth:`special_entries()`.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[4,2,1],[4,1],[2]])
            sage: G.number_of_special_entries()
            1
        """
        return len(self.special_entries())

    def is_strict(self) -> bool:
        """
        Return ``True`` if ``self`` is a strict Gelfand-Tsetlin pattern.

        A Gelfand-Tsetlin pattern is said to be *strict* if every row is
        strictly decreasing.

        EXAMPLES::

            sage: GelfandTsetlinPattern([[7,3,1],[6,2],[4]]).is_strict()
            True
            sage: GelfandTsetlinPattern([[3,2,1],[3,1],[1]]).is_strict()
            True
            sage: GelfandTsetlinPattern([[6,0,0],[3,0],[2]]).is_strict()
            False
        """
        for row in self:
            if any(row[i] == row[i+1] for i in range(len(row)-1)):
                return False
        return True

    def row_sums(self) -> list:
        r"""
        Return the list of row sums.

        For a Gelfand-Tsetlin pattern `G`, the `i`-th row sum `d_i` is

        .. MATH::

            d_i = d_i(G) = \sum_{j=i}^{n} a_{i,j}.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[5,3,2,1,0],[4,3,2,0],[4,2,1],[3,2],[3]])
            sage: G.row_sums()
            [11, 9, 7, 5, 3]
            sage: G = GelfandTsetlinPattern([[3,2,1],[3,1],[2]])
            sage: G.row_sums()
            [6, 4, 2]
        """
        return [sum(self[i][j] for j in range(len(self[i])))
                for i in range(len(self))]

    def weight(self) -> tuple:
        r"""
        Return the weight of ``self``.

        Define the weight of `G` to be the content of the tableau to which `G`
        corresponds under the bijection between Gelfand-Tsetlin patterns and
        semistandard tableaux.  More precisely,

        .. MATH::

            \mathrm{wt}(G) = (d_n, d_{n-1}-d_n, \dots, d_1-d_2),

        where the `d_i` are the row sums.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[2,1,0],[1,0],[1]])
            sage: G.weight()
            (1, 0, 2)
            sage: G = GelfandTsetlinPattern([[4,2,1],[3,1],[2]])
            sage: G.weight()
            (2, 2, 3)
        """
        wt = [self.row_sums()[-1]] + [self.row_sums()[i-1]-self.row_sums()[i] for i in reversed(range(1, len(self[0])))]
        return tuple(wt)

    def Tokuyama_coefficient(self, name='t'):
        r"""
        Return the Tokuyama coefficient attached to ``self``.

        Following the exposition of [BBF]_, Tokuyama's formula asserts

        .. MATH::

            \sum_{G} (t+1)^{s(G)} t^{l(G)}
            z_1^{d_{n+1}} z_2^{d_{n}-d_{n+1}} \cdots z_{n+1}^{d_1-d_2}
            =
            s_{\lambda}(z_1,\dots,z_{n+1}) \prod_{i<j} (z_j+tz_i),

        where the sum is over all strict Gelfand-Tsetlin patterns with fixed
        top row `\lambda + \rho`, with `\lambda` a partition with at most
        `n+1` parts and `\rho = (n, n-1, \ldots, 1, 0)`, and `s_\lambda` is a
        Schur function.

        INPUT:

        - ``name`` -- (Default: ``'t'``) An alternative name for the
          variable `t`.

        EXAMPLES::

            sage: P = GelfandTsetlinPattern([[3,2,1],[2,2],[2]])
            sage: P.Tokuyama_coefficient()
            0
            sage: G = GelfandTsetlinPattern([[3,2,1],[3,1],[2]])
            sage: G.Tokuyama_coefficient()
            t^2 + t
            sage: G = GelfandTsetlinPattern([[2,1,0],[1,1],[1]])
            sage: G.Tokuyama_coefficient()
            0
            sage: G = GelfandTsetlinPattern([[5,3,2,1,0],[4,3,2,0],[4,2,1],[3,2],[3]])
            sage: G.Tokuyama_coefficient()
            t^8 + 3*t^7 + 3*t^6 + t^5
        """
        R = PolynomialRing(ZZ, name)
        t = R.gen(0)
        if not self.is_strict():
            return R.zero()
        return (t+1)**(self.number_of_special_entries()) * t**(self.number_of_boxes())

    @combinatorial_map(order=2, name='Bender-Knuth involution')
    def bender_knuth_involution(self, i) -> GelfandTsetlinPattern:
        r"""
        Return the image of ``self`` under the `i`-th Bender-Knuth involution.

        If the triangle ``self`` has size `n` then this is defined for `0 < i < n`.

        The entries of ``self`` can take values in any ordered ring. Usually,
        this will be the integers but can also be the rationals or the real numbers.

        This implements the construction of the Bender-Knuth involution using toggling
        due to Berenstein-Kirillov.

        This agrees with the Bender-Knuth involution on semistandard tableaux.

        EXAMPLES::

            sage: G = GelfandTsetlinPattern([[5,3,2,1,0],[4,3,2,0],[4,2,1],[3,2],[3]])
            sage: G.bender_knuth_involution(2)
            [[5, 3, 2, 1, 0], [4, 3, 2, 0], [4, 2, 1], [4, 1], [3]]

            sage: G = GelfandTsetlinPattern([[3,2,0],[2.2,0],[2]])
            sage: G.bender_knuth_involution(2)
            [[3, 2, 0], [2.80000000000000, 2], [2]]

        TESTS::

            sage: all(all( G.bender_knuth_involution(i).to_tableau() == G.to_tableau().bender_knuth_involution(i)
            ....:       for i in range(1,len(G)) ) for G in GelfandTsetlinPatterns(top_row=[3,3,3,0,0]))
            True

            sage: G = GelfandTsetlinPattern([[2,1,0],[1,0],[0]])
            sage: G.bender_knuth_involution(0)
            Traceback (most recent call last):
            ...
            ValueError: must have 0 < 0 < 3
            sage: G.bender_knuth_involution(3)
            Traceback (most recent call last):
            ...
            ValueError: must have 0 < 3 < 3

        """
        n = len(self)

        def toggle(i, j):
            """
            Return the toggle of entry 'G[i][j]' in a Gelfand-Tsetlin pattern, 'G'.
            """
            if i == n-1:
                return self[n-2][0]+self[n-2][1]-self[n-1][0]

            if j == 0:
                left = self[i-1][0]
            else:
                left = min(self[i-1][j], self[i+1][j-1])
            if j == n-i-1:
                right = self[i-1][j+1]
            else:
                right = max(self[i-1][j+1], self[i+1][j])

            return left + right - self[i][j]

        if not 0 < i < n:
            raise ValueError(f"must have 0 < {i} < {n}")
        r = n-i
        P = self.parent()
        data = [list(row) for row in self]
        data[r] = [toggle(r, s) for s in range(i)]
        return P.element_class(P, data)


class GelfandTsetlinPatterns(UniqueRepresentation, Parent):
    """
    Gelfand-Tsetlin patterns.

    INPUT:

    - ``n`` -- The width or depth of the array, also known as the rank

    - ``k`` -- (Default: ``None``) If specified, this is the maximum value that
      can occur in the patterns

    - ``top_row`` -- (Default: ``None``) If specified, this is the fixed top
      row of all patterns

    - ``strict`` -- (Default: ``False``) Set to ``True`` if all patterns are
      strict patterns

    TESTS:

    Check that the number of Gelfand-Tsetlin patterns is equal to the number
    of semistandard Young tableaux::

        sage: G = GelfandTsetlinPatterns(3,3)
        sage: c = 0
        sage: from sage.combinat.crystals.kirillov_reshetikhin import partitions_in_box
        sage: for p in partitions_in_box(3,3):
        ....:    S = SemistandardTableaux(p, max_entry=3)
        ....:    c += S.cardinality()
        sage: c == G.cardinality()
        True

    Note that the top row in reverse of the Gelfand-Tsetlin pattern is the
    shape of the corresponding semistandard Young tableau under the bijection
    described in :meth:`GelfandTsetlinPattern.to_tableau()`::

        sage: G = GelfandTsetlinPatterns(top_row=[2,2,1])
        sage: S = SemistandardTableaux([2,2,1], max_entry=3)
        sage: G.cardinality() == S.cardinality()
        True
    """
    @staticmethod
    def __classcall_private__(cls, n=None, k=None, strict=False, top_row=None):
        """
        Return the correct parent based upon the inputs.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: G2 = GelfandTsetlinPatterns()
            sage: G is G2
            True
            sage: G = GelfandTsetlinPatterns(3,4, strict=True)
            sage: G2 = GelfandTsetlinPatterns(int(3),int(4), strict=True)
            sage: G is G2
            True
            sage: G = GelfandTsetlinPatterns(top_row=[3,1,1])
            sage: G2 = GelfandTsetlinPatterns(top_row=(3,1,1))
            sage: G is G2
            True
        """
        if top_row is not None:
            top_row = tuple(top_row)
            if any(top_row[i] < top_row[i+1] for i in range(len(top_row)-1)):
                raise ValueError("The top row must be weakly decreasing")
            if n is not None and n != len(top_row):
                raise ValueError("n must be the length of the specified top row")
            return GelfandTsetlinPatternsTopRow(top_row, strict)
        return super(GelfandTsetlinPatterns, cls).__classcall__(cls, n, k, strict)

    def __init__(self, n, k, strict):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: TestSuite(G).run()
            sage: G = GelfandTsetlinPatterns(3)
            sage: TestSuite(G).run()
            sage: G = GelfandTsetlinPatterns(3, 3)
            sage: TestSuite(G).run()
            sage: G = GelfandTsetlinPatterns(3, 3, strict=True)
            sage: TestSuite(G).run()
        """
        self._n = n
        self._k = k
        self._strict = strict
        # Note - if a top row is given, then n and k are not None
        if k is not None and (n is not None or strict):
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=InfiniteEnumeratedSets())

    def __contains__(self, gt):
        """
        Check to see if ``gt`` is in ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns()
            sage: [[3, 1],[2]] in G
            True
            sage: [[2, 3],[4]] in G
            False
            sage: [[3, 1],[0]] in G
            False
            sage: [] in G
            True
            sage: G = GelfandTsetlinPatterns(3,2)
            sage: [] in G
            False
            sage: [[2,0,0],[1,0],[1]] in G
            True
            sage: [[0,0],[0]] in G
            False
            sage: [[3,0,0],[2,0],[0]] in G
            False
            sage: G = GelfandTsetlinPatterns(3,strict=True)
            sage: [[2,1,0],[2,1],[1]] in G
            True
            sage: [[3,0,0],[3,0],[0]] in G
            False
        """
        if not isinstance(gt, (list, tuple, GelfandTsetlinPattern)):
            return False
        # Check if it has the correct width/depth (if applicable)
        if self._n is not None and len(gt) != self._n:
            return False
        # Check if it has the correct maximum value
        if self._k is not None and any( val > self._k for row in gt for val in row ):
            return False
        # Check if it is a GT pattern
        if not all( gt[i-1][j] >= gt[i][j] >= gt[i-1][j+1]
                    for i in range(1, len(gt)) for j in range(len(gt[i])) ):
            return False
        # Check if it is strict if applicable
        if self._strict and any( gt[i][j] == gt[i][j-1] for i in range(len(gt))
                for j in range(1, len(gt[i])) ):
            return False
        return True

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: GelfandTsetlinPatterns(4)
            Gelfand-Tsetlin patterns of width 4
            sage: GelfandTsetlinPatterns(4, 3, strict=True)
            Strict Gelfand-Tsetlin patterns of width 4 and max value 3
            sage: G = GelfandTsetlinPatterns(k=3, strict=True); G
            Strict Gelfand-Tsetlin patterns with max value 3
        """
        base = "Gelfand-Tsetlin patterns"
        if self._strict:
            base = "Strict " + base
        if self._n is not None:
            if self._k is not None:
                return base + " of width %s and max value %s" % (self._n, self._k)
            return base + " of width %s" % self._n
        if self._k is not None:
            return base + " with max value %s" % self._k
        return base

    def _element_constructor_(self, gt):
        """
        Construct an element of ``self`` from ``gt``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(3, 3, strict=True); G
            Strict Gelfand-Tsetlin patterns of width 3 and max value 3
            sage: elt = G([[3,2,1],[2,1],[1]]); elt.pp()
              3     2     1
                 2     1
                    1
            sage: elt.parent()
            Strict Gelfand-Tsetlin patterns of width 3 and max value 3
        """
        if isinstance(gt, GelfandTsetlinPattern) and gt.parent() == self:
            return gt
        if isinstance(gt, Tableau):
            gt = [list(x) for x in reversed(gt.to_chain()[1:])]
            n = len(gt)
            for i in range(n):
                while len(gt[i]) < n-i:
                    gt[i].append(0)
            if self._n is not None:
                if len(gt) == 0:
                    gt = [[0]]
                while self._n != len(gt):
                    gt.insert(0, gt[0][:] + [0])
            return self.element_class(self, gt)
        return self.element_class(self, list(gt))

    Element = GelfandTsetlinPattern

    def _coerce_map_from_(self, S):
        """
        TESTS::

            sage: t = GelfandTsetlinPattern([[1]])
            sage: t == 0
            False
            sage: t == GelfandTsetlinPattern([[1]])
            True

        Check that :trac:`25919` is fixed::

            sage: t = GelfandTsetlinPattern([[1]])
            sage: u = GelfandTsetlinPatterns()[1]
            sage: v = GelfandTsetlinPatterns(top_row=(1,))[0]
            sage: t == u
            True
            sage: u == t
            True
            sage: t == v
            True
            sage: v == t
            True
            sage: u == v
            True
            sage: v == u
            True
        """
        if isinstance(S, GelfandTsetlinPatternsTopRow):
            return True

    def __iter__(self):
        """
        Iterate through ``self`` by using a backtracing algorithm.

        EXAMPLES::

            sage: L = list(GelfandTsetlinPatterns(3,3))
            sage: c = 0
            sage: from sage.combinat.crystals.kirillov_reshetikhin import partitions_in_box
            sage: for p in partitions_in_box(3,3):
            ....:    S = SemistandardTableaux(p, max_entry=3)
            ....:    c += S.cardinality()
            sage: c == len(L)
            True
            sage: G = GelfandTsetlinPatterns(3, 3, strict=True)
            sage: all(x.is_strict() for x in G)
            True
            sage: G = GelfandTsetlinPatterns(k=3, strict=True)
            sage: all(x.is_strict() for x in G)
            True

        Checking iterator when the set is infinite::

            sage: T = GelfandTsetlinPatterns()
            sage: it = T.__iter__()
            sage: [next(it) for i in range(10)]
            [[],
             [[1]],
             [[2]],
             [[1, 1], [1]],
             [[3]],
             [[2, 1], [1]],
             [[2, 1], [2]],
             [[1, 1, 1], [1, 1], [1]],
             [[4]],
             [[3, 1], [1]]]
            sage: T = GelfandTsetlinPatterns(k=1)
            sage: it = T.__iter__()
            sage: [next(it) for i in range(10)]
            [[],
             [[0]],
             [[1]],
             [[0, 0], [0]],
             [[1, 0], [0]],
             [[1, 0], [1]],
             [[1, 1], [1]],
             [[0, 0, 0], [0, 0], [0]],
             [[1, 0, 0], [0, 0], [0]],
             [[1, 0, 0], [1, 0], [0]]]

        Check that :trac:`14718` is fixed::

            sage: T = GelfandTsetlinPatterns(1,3)
            sage: list(T)
            [[[0]],
             [[1]],
             [[2]],
             [[3]]]
        """
        # Special cases
        if self._n is None:
            yield self.element_class(self, [])
            if self._k is None:
                # Since both `n` and `k` are none, we need special consideration
                #   while iterating, so we do so by specifying the top row by
                #   using the iterator for partitions
                n = 1
                while True:
                    if self._strict:
                        P = Partitions(n, max_slope=-1)
                    else:
                        P = Partitions(n)
                    for p in P:
                        for x in GelfandTsetlinPatterns(top_row=tuple(p), strict=self._strict):
                            yield self.element_class(self, list(x))
                    n += 1
            for x in range(self._k+1):
                yield self.element_class(self, [[x]])
            n = 2
            while not self._strict or n <= self._k+1:
                for x in self._list_iter(n):
                    yield self.element_class(self, x)
                n += 1
            return
        if self._n < 0:
            return
        if self._n == 0:
            yield self.element_class(self, [])
            return
        if self._n == 1:
            if self._k is not None:
                for x in range(self._k+1):
                    yield self.element_class(self, [[x]])
            else:
                k = 1
                while True:
                    yield self.element_class(self, [[k]])
                    k += 1
            return
        for x in self._list_iter(self._n):
            yield self.element_class(self, x)

    def _list_iter(self, n):
        """
        Fast iterator which returns Gelfand-Tsetlin patterns of width ``n`` as
        lists of lists.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(3, 1)
            sage: L = [x for x in G._list_iter(3)]
            sage: len(L) == G.cardinality()
            True
            sage: type(L[0])
            <class 'list'>
        """
        # Setup the first row
        iters = [None] * n
        ret = [None] * n
        iters[0] = self._top_row_iter(n)
        ret[0] = next(iters[0])
        min_pos = 0
        iters[1] = self._row_iter(ret[0])
        pos = 1
        while pos >= min_pos:
            try:
                ret[pos] = next(iters[pos])
                pos += 1
                # If we've reached 0 width, yield and backstep
                if pos == n:
                    yield ret[:]
                    pos -= 1
                    continue
                iters[pos] = self._row_iter(ret[pos-1])
            except StopIteration:
                pos -= 1

    def _top_row_iter(self, n):
        """
        Helper iterator for the top row.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(3, 1)
            sage: for x in G._top_row_iter(3): x
            [0, 0, 0]
            [1, 0, 0]
            [1, 1, 0]
            [1, 1, 1]
            sage: G = GelfandTsetlinPatterns(3, 2, strict=True)
            sage: for x in G._top_row_iter(3): x
            [2, 1, 0]
        """
        row = [-1] * n
        pos = 0
        while pos >= 0:
            if pos == n:
                yield row[:]
                pos -= 1
                continue
            # If it would create an invalid entry, backstep
            if (pos > 0 and (row[pos] >= row[pos-1]
                    or (self._strict and row[pos] == row[pos-1]-1)) ) \
                    or (self._k is not None and row[pos] >= self._k):
                row[pos] = -1
                pos -= 1
                continue
            row[pos] += 1
            pos += 1

    def _row_iter(self, upper_row):
        """
        Helper iterator for any row with a row above it.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(3, 4)
            sage: for x in G._row_iter([4,2,1]): x
            [2, 1]
            [2, 2]
            [3, 1]
            [3, 2]
            [4, 1]
            [4, 2]
            sage: G = GelfandTsetlinPatterns(3, 2, strict=True)
            sage: for x in G._row_iter([2, 1, 0]): x
            [1, 0]
            [2, 0]
            [2, 1]
        """
        row = [x-1 for x in upper_row[1:]]
        row_len = len(row)
        pos = 0
        while pos >= 0:
            if pos == row_len:
                yield row[:]
                pos -= 1
                continue
            # If it would create an invalid entry, backstep
            if ( pos > 0 and (row[pos] >= row[pos-1] \
                    or (self._strict and row[pos] == row[pos-1]-1)) ) \
                    or row[pos] >= upper_row[pos] \
                    or (self._k is not None and row[pos] >= self._k):
                row[pos] = upper_row[pos+1] - 1
                pos -= 1
                continue
            row[pos] += 1
            pos += 1

    def _toggle_markov_chain(self, chain_state, row, col, direction):
        """
        Helper for coupling from the past. Advance the Markov chain one step.

        INPUT:

        - ``chain_state`` -- A GelfandTsetlin pattern represented as a list of lists
        - ``row`` -- The row of the cell being modified
        - ``col`` -- The column of the cell being modified
        - ``direction`` -- The direction to change the cell 1 = increase, 0 = decrease

        OUTPUT:

        ``chain_state`` is possibly modified.

        TESTS:

            sage: G=GelfandTsetlinPatterns(3,4)
            sage: state = [[3,2,1],[3,1],[2]]
            sage: G._toggle_markov_chain(state, 0, 0, 1)
            sage: state
            [[4, 2, 1], [3, 1], [2]]
            sage: G._toggle_markov_chain(state, 1, 1, 1)
            sage: state
            [[4, 2, 1], [3, 2], [2]]
            sage: G._toggle_markov_chain(state, 0, 2, 1)
            sage: state
            [[4, 2, 2], [3, 2], [2]]
            sage: G._toggle_markov_chain(state, 0, 2, 1)
            sage: state
            [[4, 2, 2], [3, 2], [2]]
            sage: G._toggle_markov_chain(state, 0, 2, 0)
            sage: state
            [[4, 2, 1], [3, 2], [2]]
            sage: G._toggle_markov_chain(state, 0, 2, 0)
            sage: state
            [[4, 2, 0], [3, 2], [2]]
            sage: G._toggle_markov_chain(state, 0, 2, 0)
            sage: state
            [[4, 2, 0], [3, 2], [2]]

            """
        if direction == 1:
            upbound = self._k
            if row != 0:
                upbound = min(upbound, chain_state[row - 1][col])
            if self._strict and col > 0:
                upbound = min(upbound, chain_state[row][col - 1] - 1)
            if row < self._n and col > 0:
                upbound = min(upbound, chain_state[row + 1][col - 1])
            if chain_state[row][col] < upbound:
                chain_state[row][col] += 1
        else:
            lobound = 0
            if row != 0:
                lobound = max(lobound, chain_state[row - 1][col + 1])
            if self._strict and col < self._n - row - 1:
                lobound = max(lobound, chain_state[row][col + 1] + 1)
            if row < self._n and col < self._n - row - 1:
                lobound = max(lobound, chain_state[row + 1][col])
            if chain_state[row][col] > lobound:
                chain_state[row][col] -= 1

    def _cftp_upper(self):
        """
        Return the largest member of the poset of Gelfand-Tsetlin patterns having the given ``n`` and ``k``.

        TESTS::

            sage: GelfandTsetlinPatterns(3, 5)._cftp_upper()
            [[5, 5, 5], [5, 5], [5]]
            sage: GelfandTsetlinPatterns(3, 5, strict=True)._cftp_upper()
            [[5, 4, 3], [5, 4], [5]]
        """
        if self._strict:
            return [[self._k - j for j in range(self._n - i)] for i in range(self._n)]
        else:
            return [[self._k for j in range(self._n - i)] for i in range(self._n)]

    def _cftp_lower(self):
        """
        Return the smallest member of the poset of Gelfand-Tsetlin patterns having the given ``n`` and ``k``.

        TESTS::

            sage: GelfandTsetlinPatterns(3, 5)._cftp_lower()
            [[0, 0, 0], [0, 0], [0]]
            sage: GelfandTsetlinPatterns(3, 5, strict=True)._cftp_lower()
            [[2, 1, 0], [1, 0], [0]]
        """
        if self._strict:
            return [[self._n - j - i - 1 for j in range(self._n - i)] for i in range(self._n)]
        else:
            return [[0 for j in range(self._n - i)] for i in range(self._n)]

    def _cftp(self, start_row):
        """
        Implement coupling from the past.

        ALGORITHM:

        The set of Gelfand-Tsetlin patterns can partially ordered by
        elementwise domination.  The partial order has unique maximum
        and minimum elements that are computed by the methods
        :meth:`_cftp_upper` and :meth:`_cftp_lower`. We then run the Markov
        chain that randomly toggles each element up or down from the
        past until the state reached from the upper and lower start
        points coalesce as described in [Propp1997]_.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(3, 5)
            sage: G._cftp(0)  # random
            [[5, 3, 2], [4, 2], [3]]
            sage: G._cftp(0) in G
            True
        """
        from sage.misc.randstate import current_randstate
        from sage.misc.randstate import seed
        from sage.misc.randstate import random

        count = self._n * self._k
        seedlist = [(current_randstate().long_seed(), count)]
        upper = []
        lower = []
        while True:
            upper = self._cftp_upper()
            lower = self._cftp_lower()
            for currseed, count in seedlist:
                with seed(currseed):
                    for _ in range(count):
                        for row in range(start_row, self._n):
                            for col in range(self._n - row):
                                direction = random() % 2
                                self._toggle_markov_chain(upper, row, col, direction)
                                self._toggle_markov_chain(lower, row, col, direction)
            if all(all(x == y for x, y in zip(l1, l2)) for l1, l2 in zip(upper, lower)):
                break
            count = seedlist[0][1] * 2
            seedlist.insert(0, (current_randstate().long_seed(), count))
        return GelfandTsetlinPattern(upper)

    def random_element(self) -> GelfandTsetlinPattern:
        """
        Return a uniformly random Gelfand-Tsetlin pattern.

        EXAMPLES::

            sage: g = GelfandTsetlinPatterns(4, 5)
            sage: x = g.random_element()
            sage: x in g
            True
            sage: len(x)
            4
            sage: all(y in range(5+1) for z in x for y in z)
            True
            sage: x.check()

        ::

            sage: g = GelfandTsetlinPatterns(4, 5, strict=True)
            sage: x = g.random_element()
            sage: x in g
            True
            sage: len(x)
            4
            sage: all(y in range(5+1) for z in x for y in z)
            True
            sage: x.check()
            sage: x.is_strict()
            True
        """
        if self._n is not None and self._k is not None:
            if self._strict and self._k+1 < self._n:
                raise ValueError('Cannot sample from empty set')
            elif self._k < 0:
                raise ValueError('Cannot sample from empty set')
            else:
                return self._cftp(0)
        else:
            raise ValueError('Cannot sample from infinite set')


class GelfandTsetlinPatternsTopRow(GelfandTsetlinPatterns):
    """
    Gelfand-Tsetlin patterns with a fixed top row.
    """
    def __init__(self, top_row, strict):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(top_row=[4,4,3,1])
            sage: TestSuite(G).run()

        TESTS:

        Check a border case in :trac:`14765`::

            sage: G = GelfandTsetlinPatterns(top_row=[])
            sage: list(G)
            [[]]
        """
        self._row = top_row
        n = len(top_row)
        if n == 0:
            k = 0
        else:
            k = top_row[0]
        GelfandTsetlinPatterns.__init__(self, n, k, strict)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: GelfandTsetlinPatterns(top_row=[4,4,3,1])
            Gelfand-Tsetlin patterns with top row [4, 4, 3, 1]
            sage: GelfandTsetlinPatterns(top_row=[5,4,3,1], strict=True)
            Strict Gelfand-Tsetlin patterns with top row [5, 4, 3, 1]
        """
        base = "Gelfand-Tsetlin patterns with top row %s" % list(self._row)
        if self._strict:
            base = "Strict " + base
        return base

    def __contains__(self, gt) -> bool:
        """
        Check if ``gt`` is in ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(top_row=[4,4,1])
            sage: [[4,4,1], [4,2], [3]] in G
            True
            sage: [[4,3,1], [4,2], [3]] in G
            False
        """
        # Check if the top row matches (if applicable)
        if gt and tuple(gt[0]) != self._row:
            return False
        return GelfandTsetlinPatterns.__contains__(self, gt)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(top_row=[4,2,1])
            sage: list(G)
            [[[4, 2, 1], [2, 1], [1]],
             [[4, 2, 1], [2, 1], [2]],
             [[4, 2, 1], [2, 2], [2]],
             [[4, 2, 1], [3, 1], [1]],
             [[4, 2, 1], [3, 1], [2]],
             [[4, 2, 1], [3, 1], [3]],
             [[4, 2, 1], [3, 2], [2]],
             [[4, 2, 1], [3, 2], [3]],
             [[4, 2, 1], [4, 1], [1]],
             [[4, 2, 1], [4, 1], [2]],
             [[4, 2, 1], [4, 1], [3]],
             [[4, 2, 1], [4, 1], [4]],
             [[4, 2, 1], [4, 2], [2]],
             [[4, 2, 1], [4, 2], [3]],
             [[4, 2, 1], [4, 2], [4]]]
        """
        # If we enforce strictness, check to see if a specified top row is strict
        if self._strict and any(self._row[i] == self._row[i+1] for i in range(self._n-1)):
            return
        if self._n == 0:
            yield self.element_class(self, [])
            return
        if self._n == 1:
            yield self.element_class(self, [list(self._row)])
            return
        # Setup the first row
        iters = [None]*self._n
        ret = [None]*self._n
        ret[0] = list(self._row)
        min_pos = 1
        iters[1] = self._row_iter(ret[0])
        pos = 1
        while pos >= min_pos:
            try:
                ret[pos] = next(iters[pos])
                pos += 1
                # If we've reached 0 width, yield and backstep
                if pos == self._n:
                    yield self.element_class(self, ret[:])
                    pos -= 1
                    continue
                iters[pos] = self._row_iter(ret[pos-1])
            except StopIteration:
                pos -= 1

    def top_row(self):
        """
        Return the top row of ``self``.

        EXAMPLES::

            sage: G = GelfandTsetlinPatterns(top_row=[4,4,3,1])
            sage: G.top_row()
            (4, 4, 3, 1)
        """
        return self._row

    def Tokuyama_formula(self, name='t'):
        r"""
        Return the Tokuyama formula of ``self``.

        Following the exposition of [BBF]_, Tokuyama's formula asserts

        .. MATH::

            \sum_{G} (t+1)^{s(G)} t^{l(G)}
            z_1^{d_{n+1}} z_2^{d_{n}-d_{n+1}} \cdots z_{n+1}^{d_1-d_2}
            = s_{\lambda} (z_1, \ldots, z_{n+1}) \prod_{i<j} (z_j+tz_i),

        where the sum is over all strict Gelfand-Tsetlin patterns with fixed
        top row `\lambda+\rho`, with `\lambda` a partition with at most
        `n+1` parts and `\rho = (n,n-1,\dots,1,0)`, and `s_{\lambda}` is a Schur
        function.

        INPUT:

        - ``name`` -- (Default: ``'t'``) An alternative name for the
          variable `t`.

        EXAMPLES::

            sage: GT = GelfandTsetlinPatterns(top_row=[2,1,0],strict=True)
            sage: GT.Tokuyama_formula()
            t^3*x1^2*x2 + t^2*x1*x2^2 + t^2*x1^2*x3 + t^2*x1*x2*x3 + t*x1*x2*x3 + t*x2^2*x3 + t*x1*x3^2 + x2*x3^2
            sage: GT = GelfandTsetlinPatterns(top_row=[3,2,1],strict=True)
            sage: GT.Tokuyama_formula()
            t^3*x1^3*x2^2*x3 + t^2*x1^2*x2^3*x3 + t^2*x1^3*x2*x3^2 + t^2*x1^2*x2^2*x3^2 + t*x1^2*x2^2*x3^2 + t*x1*x2^3*x3^2 + t*x1^2*x2*x3^3 + x1*x2^2*x3^3
            sage: GT = GelfandTsetlinPatterns(top_row=[1,1,1],strict=True)
            sage: GT.Tokuyama_formula()
            0
        """
        n = self._n
        variables = [name] + ["x%d" % i for i in range(1, n+1)]
        R = PolynomialRing(ZZ, names=variables)
        t = R.gen(0)
        x = R.gens()[1:]
        GT = GelfandTsetlinPatterns(top_row=self._row, strict=True)
        return sum((t+1)**(gt.number_of_special_entries()) * t**(gt.number_of_boxes()) * prod(x[i]**gt.weight()[i] for i in range(n)) for gt in GT)

    def _cftp_upper(self) -> list:
        """
        Return the largest member of the poset of Gelfand-Tsetlin patterns having the given ``top_row``.

        TESTS::

            sage: GelfandTsetlinPatterns(top_row = [5, 4, 3])._cftp_upper()
            [[5, 4, 3], [5, 4], [5]]
        """
        return [[self._row[j] for j in range(self._n - i)] for i in range(self._n)]

    def _cftp_lower(self) -> list:
        """
        Return the smallest member of the poset of Gelfand-Tsetlin patterns having the given ``top_row``.

        TESTS::

            sage: GelfandTsetlinPatterns(top_row = [5, 4, 3])._cftp_lower()
            [[5, 4, 3], [4, 3], [3]]
        """
        return [[self._row[i + j] for j in range(self._n - i)] for i in range(self._n)]

    def random_element(self) -> GelfandTsetlinPattern:
        """
        Return a uniformly random Gelfand-Tsetlin pattern with specified top row.

        EXAMPLES::

            sage: g = GelfandTsetlinPatterns(top_row = [4, 3, 1, 1])
            sage: x = g.random_element()
            sage: x in g
            True
            sage: x[0] == [4, 3, 1, 1]
            True
            sage: x.check()

            sage: g = GelfandTsetlinPatterns(top_row=[4, 3, 2, 1], strict=True)
            sage: x = g.random_element()
            sage: x in g
            True
            sage: x[0] == [4, 3, 2, 1]
            True
            sage: x.is_strict()
            True
            sage: x.check()
        """
        if self._strict:
            return self._cftp(1)
        l = [i for i in self._row if i > 0]
        return SemistandardTableaux(l, max_entry=self._n).random_element().to_Gelfand_Tsetlin_pattern()  # type:ignore
