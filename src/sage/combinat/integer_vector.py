"""
(Non-negative) Integer vectors

AUTHORS:

* Mike Hansen (2007) - original module
* Nathann Cohen, David Joyner (2009-2010) - Gale-Ryser stuff
* Nathann Cohen, David Joyner (2011) - Gale-Ryser bugfix
* Travis Scrimshaw (2012-05-12) - Updated doc-strings to tell the user of
  that the class's name is a misnomer (that they only contains non-negative
  entries).
* Federico Poloni (2013) - specialized ``rank()``
* Travis Scrimshaw (2013-02-04) - Refactored to use ``ClonableIntArray``
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from itertools import product
from collections.abc import Sequence
import numbers

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.misc.classcall_metaclass import ClasscallMetaclass

from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.rings.infinity import PlusInfinity
from sage.arith.misc import binomial
from sage.rings.integer_ring import ZZ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.rings.integer import Integer


def is_gale_ryser(r, s):
    r"""
    Test whether the given sequences satisfy the condition
    of the Gale-Ryser theorem.

    Given a binary matrix `B` of dimension `n\times m`, the
    vector of row sums is defined as the vector whose
    `i^{\mbox{th}}` component is equal to the sum of the `i^{\mbox{th}}`
    row in `A`. The vector of column sums is defined similarly.

    If, given a binary matrix, these two vectors are easy to compute,
    the Gale-Ryser theorem lets us decide whether, given two
    nonnegative vectors `r,s`, there exists a binary matrix
    whose row/column sums vectors are `r` and `s`.

    This functions answers accordingly.

    INPUT:

    - ``r``, ``s`` -- lists of nonnegative integers

    ALGORITHM:

    Without loss of generality, we can assume that:

    - The two given sequences do not contain any `0` ( which would
      correspond to an empty column/row )

    - The two given sequences are ordered in decreasing order
      (reordering the sequence of row (resp. column) sums amounts to
      reordering the rows (resp. columns) themselves in the matrix,
      which does not alter the columns (resp. rows) sums.

    We can then assume that `r` and `s` are partitions
    (see the corresponding class :class:`Partition`)

    If `r^*` denote the conjugate of `r`, the Gale-Ryser theorem
    asserts that a binary Matrix satisfying the constraints exists
    if and only if `s \preceq r^*`, where `\preceq` denotes
    the domination order on partitions.

    EXAMPLES::

        sage: from sage.combinat.integer_vector import is_gale_ryser
        sage: is_gale_ryser([4,2,2], [3,3,1,1])                                         # needs sage.combinat
        True
        sage: is_gale_ryser([4,2,1,1], [3,3,1,1])                                       # needs sage.combinat
        True
        sage: is_gale_ryser([3,2,1,1], [3,3,1,1])                                       # needs sage.combinat
        False

    REMARK: In the literature, what we are calling a
    Gale-Ryser sequence sometimes goes by the (rather
    generic-sounding) term ''realizable sequence''.
    """

    # The sequences only contain nonnegative integers
    if [x for x in r if x < 0] or [x for x in s if x < 0]:
        return False

    # builds the corresponding partitions, i.e.
    # removes the 0 and sorts the sequences
    from sage.combinat.partition import Partition
    r2 = Partition(sorted([x for x in r if x > 0], reverse=True))
    s2 = Partition(sorted([x for x in s if x > 0], reverse=True))

    # If the two sequences only contained zeroes
    if len(r2) == 0 and len(s2) == 0:
        return True

    rstar = Partition(r2).conjugate()

    #                                same number of 1s           domination
    return len(rstar) <= len(s2) and sum(r2) == sum(s2) and rstar.dominates(s)


def gale_ryser_theorem(p1, p2, algorithm='gale',
                       *, solver=None, integrality_tolerance=1e-3):
    r"""
    Return the binary matrix given by the Gale-Ryser theorem.

    The Gale Ryser theorem asserts that if `p_1,p_2` are two
    partitions of `n` of respective lengths `k_1,k_2`, then there is
    a binary `k_1\times k_2` matrix `M` such that `p_1` is the vector
    of row sums and `p_2` is the vector of column sums of `M`, if
    and only if the conjugate of `p_2` dominates `p_1`.

    INPUT:

    - ``p1``, ``p2`` -- list of integers representing the vectors
      of row/column sums

    - ``algorithm`` -- two possible string values:

      - ``'ryser'`` implements the construction due to Ryser [Ryser63]_.
      - ``'gale'`` (default) implements the construction due to Gale [Gale57]_.

    - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear Programming
      (MILP) solver to be used. If set to ``None``, the default one is used. For
      more information on MILP solvers and which default solver is used, see
      the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
      of the class
      :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``integrality_tolerance`` -- parameter for use with MILP solvers over an
      inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`

    OUTPUT: a binary matrix if it exists, ``None`` otherwise

    Gale's Algorithm:

    (Gale [Gale57]_): A matrix satisfying the constraints of its
    sums can be defined as the solution of the following
    Linear Program, which Sage knows how to solve.

    .. MATH::

        \forall i&\sum_{j=1}^{k_2} b_{i,j}=p_{1,j}\\
        \forall i&\sum_{j=1}^{k_1} b_{j,i}=p_{2,j}\\
        &b_{i,j}\mbox{ is a binary variable}

    Ryser's Algorithm:

    (Ryser [Ryser63]_): The construction of an `m \times n` matrix
    `A=A_{r,s}`, due to Ryser, is described as follows. The
    construction works if and only if have `s\preceq r^*`.

    * Construct the `m \times n` matrix `B` from `r` by defining
      the `i`-th row of `B` to be the vector whose first `r_i`
      entries are `1`, and the remainder are 0s, `1 \leq i \leq m`.
      This maximal matrix `B` with row sum `r` and ones left
      justified has column sum `r^{*}`.

    * Shift the last `1` in certain rows of `B` to column `n` in
      order to achieve the sum `s_n`.  Call this `B` again.

      * The `1`'s in column `n` are to appear in those
        rows in which `A` has the largest row sums, giving
        preference to the bottom-most positions in case of ties.
      * Note: When this step automatically "fixes" other columns,
        one must skip ahead to the first column index
        with a wrong sum in the step below.

    * Proceed inductively to construct columns `n-1`, ..., `2`, `1`.
      Note: when performing the induction on step `k`, we consider
      the row sums of the first `k` columns.

    * Set `A = B`. Return `A`.

    EXAMPLES:

    Computing the matrix for `p_1=p_2=2+2+1`::

        sage: # needs sage.combinat sage.modules
        sage: from sage.combinat.integer_vector import gale_ryser_theorem
        sage: p1 = [2,2,1]
        sage: p2 = [2,2,1]
        sage: print(gale_ryser_theorem(p1, p2))         # not tested
        [1 1 0]
        [1 0 1]
        [0 1 0]
        sage: A = gale_ryser_theorem(p1, p2)
        sage: rs = [sum(x) for x in A.rows()]
        sage: cs = [sum(x) for x in A.columns()]
        sage: p1 == rs; p2 == cs
        True
        True

    Or for a non-square matrix with `p_1=3+3+2+1` and `p_2=3+2+2+1+1`,
    using Ryser's algorithm::

        sage: # needs sage.combinat sage.modules
        sage: from sage.combinat.integer_vector import gale_ryser_theorem
        sage: p1 = [3,3,1,1]
        sage: p2 = [3,3,1,1]
        sage: gale_ryser_theorem(p1, p2, algorithm='ryser')
        [1 1 1 0]
        [1 1 0 1]
        [1 0 0 0]
        [0 1 0 0]
        sage: p1 = [4,2,2]
        sage: p2 = [3,3,1,1]
        sage: gale_ryser_theorem(p1, p2, algorithm='ryser')
        [1 1 1 1]
        [1 1 0 0]
        [1 1 0 0]
        sage: p1 = [4,2,2,0]
        sage: p2 = [3,3,1,1,0,0]
        sage: gale_ryser_theorem(p1, p2, algorithm='ryser')
        [1 1 1 1 0 0]
        [1 1 0 0 0 0]
        [1 1 0 0 0 0]
        [0 0 0 0 0 0]
        sage: p1 = [3,3,2,1]
        sage: p2 = [3,2,2,1,1]
        sage: print(gale_ryser_theorem(p1, p2, algorithm='gale'))       # not tested
        [1 1 1 0 0]
        [1 1 0 0 1]
        [1 0 1 0 0]
        [0 0 0 1 0]

    With `0` in the sequences, and with unordered inputs::

        sage: from sage.combinat.integer_vector import gale_ryser_theorem
        sage: gale_ryser_theorem([3,3,0,1,1,0], [3,1,3,1,0], algorithm='ryser')         # needs sage.combinat sage.modules
        [1 1 1 0 0]
        [1 0 1 1 0]
        [0 0 0 0 0]
        [1 0 0 0 0]
        [0 0 1 0 0]
        [0 0 0 0 0]
        sage: p1 = [3,1,1,1,1]; p2 = [3,2,2,0]
        sage: gale_ryser_theorem(p1, p2, algorithm='ryser')                             # needs sage.combinat sage.modules
        [1 1 1 0]
        [1 0 0 0]
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]

    TESTS:

    This test created a random bipartite graph on `n+m` vertices. Its
    adjacency matrix is binary, and it is used to create some
    "random-looking" sequences which correspond to an existing matrix. The
    ``gale_ryser_theorem`` is then called on these sequences, and the output
    checked for correction.::

        sage: def test_algorithm(algorithm, low=10, high=50):
        ....:    n,m = randint(low,high), randint(low,high)
        ....:    g = graphs.RandomBipartite(n, m, .3)
        ....:    s1 = sorted(g.degree([(0,i) for i in range(n)]), reverse = True)
        ....:    s2 = sorted(g.degree([(1,i) for i in range(m)]), reverse = True)
        ....:    m = gale_ryser_theorem(s1, s2, algorithm = algorithm)
        ....:    ss1 = sorted(map(lambda x : sum(x) , m.rows()), reverse = True)
        ....:    ss2 = sorted(map(lambda x : sum(x) , m.columns()), reverse = True)
        ....:    if ((ss1 != s1) or (ss2 != s2)):
        ....:        print("Algorithm %s failed with this input:" % algorithm)
        ....:        print(s1, s2)

        sage: for algorithm in ["gale", "ryser"]:       # long time                     # needs sage.combinat sage.modules
        ....:    for i in range(50):
        ....:       test_algorithm(algorithm, 3, 10)

    Null matrix::

        sage: gale_ryser_theorem([0,0,0],[0,0,0,0], algorithm='gale')                   # needs sage.combinat sage.modules
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: gale_ryser_theorem([0,0,0],[0,0,0,0], algorithm='ryser')                  # needs sage.combinat sage.modules
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]

    REFERENCES:

    ..  [Ryser63] \H. J. Ryser, Combinatorial Mathematics,
        Carus Monographs, MAA, 1963.
    ..  [Gale57] \D. Gale, A theorem on flows in networks, Pacific J. Math.
        7(1957)1073-1082.
    """
    from sage.matrix.constructor import matrix

    if not is_gale_ryser(p1, p2):
        return False

    if algorithm == "ryser":  # ryser's algorithm
        from sage.combinat.permutation import Permutation

        # Sorts the sequences if they are not, and remembers the permutation
        # applied
        tmp = sorted(enumerate(p1), reverse=True, key=lambda x: x[1])
        r = [x[1] for x in tmp]
        r_permutation = [x-1 for x in Permutation([x[0]+1 for x in tmp]).inverse()]
        m = len(r)

        tmp = sorted(enumerate(p2), reverse=True, key=lambda x: x[1])
        s = [x[1] for x in tmp]
        s_permutation = [x-1 for x in Permutation([x[0]+1 for x in tmp]).inverse()]

        # This is the partition equivalent to the sliding algorithm
        cols = []
        for t in reversed(s):
            c = [0] * m
            i = 0
            while t:
                k = i + 1
                while k < m and r[i] == r[k]:
                    k += 1
                if t >= k - i:  # == number rows of the same length
                    for j in range(i, k):
                        r[j] -= 1
                        c[j] = 1
                    t -= k - i
                else:  # Remove the t last rows of that length
                    for j in range(k-t, k):
                        r[j] -= 1
                        c[j] = 1
                    t = 0
                i = k
            cols.append(c)

        # We added columns to the back instead of the front
        A0 = matrix(list(reversed(cols))).transpose()

        # Applying the permutations to get a matrix satisfying the
        # order given by the input
        A0 = A0.matrix_from_rows_and_columns(r_permutation, s_permutation)
        return A0

    elif algorithm == "gale":
        from sage.numerical.mip import MixedIntegerLinearProgram
        k1, k2 = len(p1), len(p2)
        p = MixedIntegerLinearProgram(solver=solver)
        b = p.new_variable(binary=True)
        for (i, c) in enumerate(p1):
            p.add_constraint(p.sum([b[i, j] for j in range(k2)]) == c)
        for (i, c) in enumerate(p2):
            p.add_constraint(p.sum([b[j, i] for j in range(k1)]) == c)
        p.set_objective(None)
        p.solve()
        b = p.get_values(b, convert=ZZ, tolerance=integrality_tolerance)
        M = [[0]*k2 for i in range(k1)]
        for i in range(k1):
            for j in range(k2):
                M[i][j] = b[i, j]
        return matrix(M)

    else:
        raise ValueError('the only two algorithms available are "gale" and "ryser"')


def _default_function(l, default, i):
    """
    EXAMPLES::

        sage: from sage.combinat.integer_vector import _default_function
        sage: import functools
        sage: f = functools.partial(_default_function, [1,2,3], 99)
        sage: f(-1)
        99
        sage: f(0)
        1
        sage: f(1)
        2
        sage: f(2)
        3
        sage: f(3)
        99
    """
    try:
        if i < 0:
            return default
        return l[i]
    except IndexError:
        return default


def list2func(l, default=None):
    """
    Given a list ``l``, return a function that takes in a value ``i`` and
    return ``l[i]``. If default is not ``None``, then the function will
    return the default value for out of range ``i``'s.

    EXAMPLES::

        sage: f = sage.combinat.integer_vector.list2func([1,2,3])
        sage: f(0)
        1
        sage: f(1)
        2
        sage: f(2)
        3
        sage: f(3)
        Traceback (most recent call last):
        ...
        IndexError: list index out of range

    ::

        sage: f = sage.combinat.integer_vector.list2func([1,2,3], 0)
        sage: f(2)
        3
        sage: f(3)
        0
    """
    if default is None:
        return lambda i: l[i]
    else:
        from functools import partial
        return partial(_default_function, l, default)


class IntegerVector(ClonableArray):
    """
    An integer vector.
    """

    def check(self):
        """
        Check to make sure this is a valid integer vector by making sure
        all entries are nonnegative.

        EXAMPLES::

            sage: IV = IntegerVectors()
            sage: elt = IV([1,2,1])
            sage: elt.check()

        Check :issue:`34510`::

            sage: IV3 = IntegerVectors(n=3)
            sage: IV3([2,2])
            Traceback (most recent call last):
            ...
            ValueError: [2, 2] doesn't satisfy correct constraints
            sage: IVk3 = IntegerVectors(k=3)
            sage: IVk3([2,2])
            Traceback (most recent call last):
            ...
            ValueError: [2, 2] doesn't satisfy correct constraints
            sage: IV33 = IntegerVectors(n=3, k=3)
            sage: IV33([2,2])
            Traceback (most recent call last):
            ...
            ValueError: [2, 2] doesn't satisfy correct constraints
        """
        if any(x < 0 for x in self):
            raise ValueError("all entries must be nonnegative")
        if self not in self.parent():
            raise ValueError(f"{self} doesn't satisfy correct constraints")

    def trim(self):
        """
        Remove trailing zeros from the integer vector.

        EXAMPLES::

            sage: IV = IntegerVectors()
            sage: IV([5,3,5,1,0,0]).trim()
            [5, 3, 5, 1]
            sage: IV([5,0,5,1,0]).trim()
            [5, 0, 5, 1]
            sage: IV([4,3,3]).trim()
            [4, 3, 3]
            sage: IV([0,0,0]).trim()
            []

            sage: IV = IntegerVectors(k=4)
            sage: v = IV([4,3,2,0]).trim(); v
            [4, 3, 2]
            sage: v.parent()
            Integer vectors
        """
        P = IntegerVectors()
        v = list(self)
        if all(i == 0 for i in v):
            return P.element_class(P, [], check=False)
        while not v[-1]:
            v = v[:-1]
        return P.element_class(P, v, check=False)

    def specht_module(self, base_ring=None):
        r"""
        Return the Specht module corresponding to ``self``.

        EXAMPLES::

            sage: SM = IntegerVectors()([2,0,1,0,2]).specht_module(QQ); SM              # needs sage.combinat sage.modules
            Specht module of [(0, 0), (0, 1), (2, 0), (4, 0), (4, 1)] over Rational Field
            sage: s = SymmetricFunctions(QQ).s()                                        # needs sage.combinat sage.modules
            sage: s(SM.frobenius_image())                                               # needs sage.combinat sage.modules
            s[2, 2, 1]
        """
        from sage.combinat.specht_module import SpechtModule
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ
        R = SymmetricGroupAlgebra(base_ring, sum(self))
        return SpechtModule(R, self)

    def specht_module_dimension(self, base_ring=None):
        r"""
        Return the dimension of the Specht module corresponding to ``self``.

        INPUT:

        - ``BR`` -- (default: `\QQ`) the base ring

        EXAMPLES::

            sage: IntegerVectors()([2,0,1,0,2]).specht_module_dimension()               # needs sage.combinat sage.modules
            5
            sage: IntegerVectors()([2,0,1,0,2]).specht_module_dimension(GF(2))          # needs sage.combinat sage.modules sage.rings.finite_rings
            5
        """
        from sage.combinat.specht_module import specht_module_rank
        return specht_module_rank(self, base_ring)


class IntegerVectors(Parent, metaclass=ClasscallMetaclass):
    """
    The class of (nonnegative) integer vectors.

    INPUT:

    - ``n`` -- if set to an integer, returns the combinatorial class
      of integer vectors whose sum is ``n``; if set to ``None``
      (default), no such constraint is defined

    - ``k`` -- the length of the vectors; set to ``None`` (default) if
      you do not want such a constraint

    .. NOTE::

        The entries are nonnegative integers.

    EXAMPLES:

    If ``n`` is not specified, it returns the class of all integer vectors::

        sage: IntegerVectors()
        Integer vectors
        sage: [] in IntegerVectors()
        True
        sage: [1,2,1] in IntegerVectors()
        True
        sage: [1, 0, 0] in IntegerVectors()
        True

    Entries are nonnegative::

        sage: [-1, 2] in IntegerVectors()
        False

    If ``n`` is specified, then it returns the class of all integer vectors
    which sum to ``n``::

        sage: IV3 = IntegerVectors(3); IV3
        Integer vectors that sum to 3

    Note that trailing zeros are ignored so that ``[3, 0]`` does not show
    up in the following list (since ``[3]`` does)::

        sage: IntegerVectors(3, max_length=2).list()
        [[3], [2, 1], [1, 2], [0, 3]]

    If ``n`` and ``k`` are both specified, then it returns the class
    of integer vectors that sum to ``n`` and are of length ``k``::

        sage: IV53 = IntegerVectors(5,3); IV53
        Integer vectors of length 3 that sum to 5
        sage: IV53.cardinality()
        21
        sage: IV53.first()
        [5, 0, 0]
        sage: IV53.last()
        [0, 0, 5]
        sage: IV53.random_element().parent() is IV53
        True

    Further examples::

        sage: IntegerVectors(-1, 0, min_part=1).list()
        []
        sage: IntegerVectors(-1, 2, min_part=1).list()
        []
        sage: IntegerVectors(0, 0, min_part=1).list()
        [[]]
        sage: IntegerVectors(3, 0, min_part=1).list()
        []
        sage: IntegerVectors(0, 1, min_part=1).list()
        []
        sage: IntegerVectors(2, 2, min_part=1).list()
        [[1, 1]]
        sage: IntegerVectors(2, 3, min_part=1).list()
        []
        sage: IntegerVectors(4, 2, min_part=1).list()
        [[3, 1], [2, 2], [1, 3]]

    ::

        sage: IntegerVectors(0, 3, outer=[0,0,0]).list()
        [[0, 0, 0]]
        sage: IntegerVectors(1, 3, outer=[0,0,0]).list()
        []
        sage: IntegerVectors(2, 3, outer=[0,2,0]).list()
        [[0, 2, 0]]
        sage: IntegerVectors(2, 3, outer=[1,2,1]).list()
        [[1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1]]
        sage: IntegerVectors(2, 3, outer=[1,1,1]).list()
        [[1, 1, 0], [1, 0, 1], [0, 1, 1]]
        sage: IntegerVectors(2, 5, outer=[1,1,1,1,1]).list()
        [[1, 1, 0, 0, 0],
         [1, 0, 1, 0, 0],
         [1, 0, 0, 1, 0],
         [1, 0, 0, 0, 1],
         [0, 1, 1, 0, 0],
         [0, 1, 0, 1, 0],
         [0, 1, 0, 0, 1],
         [0, 0, 1, 1, 0],
         [0, 0, 1, 0, 1],
         [0, 0, 0, 1, 1]]

    ::

        sage: iv = [ IntegerVectors(n,k) for n in range(-2, 7) for k in range(7) ]
        sage: all(map(lambda x: x.cardinality() == len(x.list()), iv))
        True
        sage: essai = [[1,1,1], [2,5,6], [6,5,2]]
        sage: iv = [ IntegerVectors(x[0], x[1], max_part = x[2]-1) for x in essai ]
        sage: all(map(lambda x: x.cardinality() == len(x.list()), iv))
        True

    An example showing the same output by using IntegerListsLex::

        sage: IntegerVectors(4, max_length=2).list()
        [[4], [3, 1], [2, 2], [1, 3], [0, 4]]
        sage: list(IntegerListsLex(4, max_length=2))
        [[4], [3, 1], [2, 2], [1, 3], [0, 4]]

    .. SEEALSO::

        :class:`sage.combinat.integer_lists.invlex.IntegerListsLex`
    """
    @staticmethod
    def __classcall_private__(cls, n=None, k=None, **kwargs):
        """
        Choose the correct parent based upon input.

        EXAMPLES::

            sage: IV1 = IntegerVectors(3, 2)
            sage: IV2 = IntegerVectors(3, 2)
            sage: IV1 is IV2
            True

        TESTS::

            sage: IV2 = IntegerVectors(3, 2, length=2)
            Traceback (most recent call last):
            ...
            ValueError: k and length both specified

        :issue:`29524`::

            sage: IntegerVectors(3, 3/1)
            Traceback (most recent call last):
            ...
            TypeError: 'k' must be an integer or a tuple, got Rational
        """
        if 'length' in kwargs:
            if k is not None:
                raise ValueError("k and length both specified")
            k = kwargs.pop('length')
        if kwargs:
            return IntegerVectorsConstraints(n, k, **kwargs)

        if k is None:
            if n is None:
                return IntegerVectors_all()
            return IntegerVectors_n(n)
        if n is None:
            return IntegerVectors_k(k)

        if isinstance(k, numbers.Integral):
            return IntegerVectors_nk(n, k)
        elif isinstance(k, (tuple, list)):
            return IntegerVectors_nnondescents(n, tuple(k))
        else:
            raise TypeError("'k' must be an integer or a tuple, got {}".format(type(k).__name__))

    def __init__(self, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: IV = IntegerVectors()
            sage: TestSuite(IV).run()
        """
        if category is None:
            category = EnumeratedSets()
        Parent.__init__(self, category=category)

    def _element_constructor_(self, lst):
        """
        Construct an element of ``self`` from ``lst``.

        EXAMPLES::

            sage: IV = IntegerVectors()
            sage: elt = IV([3, 1, 0, 3, 2]); elt
            [3, 1, 0, 3, 2]
            sage: elt.parent()
            Integer vectors

            sage: IV9 = IntegerVectors(9)
            sage: elt9 = IV9(elt)
            sage: elt9.parent()
            Integer vectors that sum to 9
        """
        return self.element_class(self, lst)

    Element = IntegerVector

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [] in IntegerVectors()
            True
            sage: [3,2,2,1] in IntegerVectors()
            True
        """
        if isinstance(x, IntegerVector):
            return True

        if not isinstance(x, Sequence):
            return False

        for i in x:
            if i not in ZZ:
                return False
            if i < 0:
                return False
        return True

    def _unrank_helper(self, x, rtn):
        """
        Return the element at rank ``x`` by iterating through all integer vectors beginning with ``rtn``.

        INPUT:

        - ``x`` -- nonnegative integer
        - ``rtn`` -- list of nonnegative integers

        EXAMPLES::

            sage: IV = IntegerVectors(k=5)
            sage: IV._unrank_helper(10, [2,0,0,0,0])
            [1, 0, 0, 0, 1]

            sage: IV = IntegerVectors(n=7)
            sage: IV._unrank_helper(100, [7,0,0,0])
            [2, 0, 0, 5]

            sage: IV = IntegerVectors(n=12, k=7)
            sage: IV._unrank_helper(1000, [12,0,0,0,0,0,0])
            [5, 3, 1, 1, 1, 1, 0]
        """
        ptr = 0
        while True:
            current_rank = self.rank(rtn)
            if current_rank < x:
                rtn[ptr+1] = rtn[ptr]
                rtn[ptr] = 0
                ptr += 1
            elif current_rank > x:
                rtn[ptr] -= 1
                rtn[ptr-1] += 1
            else:
                return self._element_constructor_(rtn)


class IntegerVectors_all(UniqueRepresentation, IntegerVectors):
    """
    Class of all integer vectors.
    """

    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: IV = IntegerVectors()
            sage: TestSuite(IV).run()
        """
        IntegerVectors.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: IntegerVectors()
            Integer vectors
        """
        return "Integer vectors"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: IV = IntegerVectors()
            sage: it = IV.__iter__()
            sage: [next(it) for x in range(10)]
            [[], [1], [2], [2, 0], [1, 1], [0, 2], [3], [3, 0], [2, 1], [1, 2]]
        """
        yield self.element_class(self, [])
        n = 1
        while True:
            for k in range(1, n + 1):
                for v in integer_vectors_nk_fast_iter(n, k):
                    yield self.element_class(self, v, check=False)
            n += 1


class IntegerVectors_n(UniqueRepresentation, IntegerVectors):
    """
    Integer vectors that sum to `n`.
    """

    def __init__(self, n):
        """
        TESTS::

            sage: IV = IntegerVectors(3)
            sage: TestSuite(IV).run()
        """
        self.n = n
        if self.n == 0:
            IntegerVectors.__init__(self, category=EnumeratedSets())
        else:
            IntegerVectors.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: IV = IntegerVectors(3)
            sage: IV
            Integer vectors that sum to 3
        """
        return "Integer vectors that sum to {}".format(self.n)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = IntegerVectors(3).__iter__()
            sage: [next(it) for x in range(10)]
            [[3],
             [3, 0],
             [2, 1],
             [1, 2],
             [0, 3],
             [3, 0, 0],
             [2, 1, 0],
             [2, 0, 1],
             [1, 2, 0],
             [1, 1, 1]]
        """
        if not self.n:
            yield self.element_class(self, [])

        k = 1
        while True:
            for iv in integer_vectors_nk_fast_iter(self.n, k):
                yield self.element_class(self, iv, check=False)
            k += 1

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [0] in IntegerVectors(0)
            True
            sage: [3] in IntegerVectors(3)
            True
            sage: [3] in IntegerVectors(2)
            False
            sage: [3,2,2,1] in IntegerVectors(9)
            False
            sage: [3,2,2,1] in IntegerVectors(8)
            True
        """
        if not IntegerVectors.__contains__(self, x):
            return False
        return sum(x) == self.n

    def rank(self, x):
        """
        Return the rank of a given element.

        INPUT:

        - ``x`` -- list with ``sum(x) == n``

        EXAMPLES::

            sage: IntegerVectors(n=5).rank([5,0])
            1
            sage: IntegerVectors(n=5).rank([3,2])
            3
        """
        if sum(x) != self.n:
            raise ValueError("argument is not a member of IntegerVectors({},{})".format(self.n, None))

        n, k, s = self.n, len(x), 0
        r = binomial(k + n - 1, n + 1)
        for i in range(k - 1):
            s += x[k - 1 - i]
            r += binomial(s + i, i + 1)
        return r

    def unrank(self, x):
        """
        Return the element at given rank x.

        INPUT:

        - ``x`` -- integer

        EXAMPLES::

            sage: IntegerVectors(n=5).unrank(2)
            [4, 1]
            sage: IntegerVectors(n=10).unrank(10)
            [1, 9]
        """
        rtn = [self.n]
        while self.rank(rtn) <= x:
            rtn.append(0)
        rtn.pop()

        return IntegerVectors._unrank_helper(self, x, rtn)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: IntegerVectors(n=0).cardinality()
            1
            sage: IntegerVectors(n=10).cardinality()
            +Infinity
        """
        if self.n == 0:
            return Integer(1)
        return PlusInfinity()


class IntegerVectors_k(UniqueRepresentation, IntegerVectors):
    """
    Integer vectors of length `k`.
    """

    def __init__(self, k):
        """
        TESTS::

            sage: IV = IntegerVectors(k=2)
            sage: TestSuite(IV).run()
        """
        self.k = k
        if self.k == 0:
            IntegerVectors.__init__(self, category=EnumeratedSets())
        else:
            IntegerVectors.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: IV = IntegerVectors(k=2)
            sage: IV
            Integer vectors of length 2
        """
        return "Integer vectors of length {}".format(self.k)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = IntegerVectors(k=2).__iter__()
            sage: [next(it) for x in range(10)]
            [[0, 0],
             [1, 0],
             [0, 1],
             [2, 0],
             [1, 1],
             [0, 2],
             [3, 0],
             [2, 1],
             [1, 2],
             [0, 3]]

        TESTS:

        Check corner case::

            sage: IV = IntegerVectors(k=0)
            sage: list(IV)
            [[]]
        """
        n = 0
        if self.k == 0:  # special case
            yield self.element_class(self, [], check=False)
            return
        while True:
            for iv in integer_vectors_nk_fast_iter(n, self.k):
                yield self.element_class(self, iv, check=False)
            n += 1

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [] in IntegerVectors(k=0)
            True
            sage: [3] in IntegerVectors(k=1)
            True
            sage: [3] in IntegerVectors(k=2)
            False
            sage: [3,2,2,1] in IntegerVectors(k=3)
            False
            sage: [3,2,2,1] in IntegerVectors(k=4)
            True
        """
        if not IntegerVectors.__contains__(self, x):
            return False
        return len(x) == self.k

    def rank(self, x):
        """
        Return the rank of a given element.

        INPUT:

        - ``x`` -- list with ``len(x) == k``

        EXAMPLES::

            sage: IntegerVectors(k=5).rank([0,0,0,0,0])
            0
            sage: IntegerVectors(k=5).rank([1,1,0,0,0])
            7
        """
        if len(x) != self.k:
            raise ValueError("argument is not a member of IntegerVectors({},{})".format(None, self.k))

        n, k, s = sum(x), self.k, 0
        r = binomial(n + k - 1, k)
        for i in range(k - 1):
            s += x[k - 1 - i]
            r += binomial(s + i, i + 1)
        return r

    def unrank(self, x):
        """
        Return the element at given rank x.

        INPUT:

        - ``x`` -- integer such that ``x < self.cardinality()``

        EXAMPLES::

            sage: IntegerVectors(k=5).unrank(10)
            [1, 0, 0, 0, 1]
            sage: IntegerVectors(k=5).unrank(15)
            [0, 0, 2, 0, 0]
            sage: IntegerVectors(k=0).unrank(0)
            []
        """
        if self.k == 0 and x != 0:
            raise IndexError(f"Index {x} is out of range for the IntegerVector.")
        rtn = [0]*self.k
        if self.k == 0 and x == 0:
            return rtn

        while self.rank(rtn) <= x:
            rtn[0] += 1
        rtn[0] -= 1

        return IntegerVectors._unrank_helper(self, x, rtn)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: IntegerVectors(k=0).cardinality()
            1
            sage: IntegerVectors(k=10).cardinality()
            +Infinity
        """
        if self.k == 0:
            return Integer(1)
        return PlusInfinity()


class IntegerVectors_nk(UniqueRepresentation, IntegerVectors):
    """
    Integer vectors of length `k` that sum to `n`.

    AUTHORS:

    - Martin Albrecht
    - Mike Hansen
    """

    def __init__(self, n, k):
        """
        TESTS::

            sage: IV = IntegerVectors(2, 3)
            sage: TestSuite(IV).run()
        """
        self.n = n
        self.k = k
        IntegerVectors.__init__(self, category=FiniteEnumeratedSets())

    def _list_rec(self, n, k):
        """
        Return a list of a exponent tuples of length ``size`` such
        that the degree of the associated monomial is `D`.

        INPUT:

        - ``n`` -- degree (must be 0)

        - ``k`` -- length of exponent tuples (must be 0)

        EXAMPLES::

            sage: IV = IntegerVectors(2,3)
            sage: list(IV._list_rec(2,3))
            [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
        """
        if k == 1:
            yield (n,)
            return

        for nbar in range(n + 1):
            n_diff = n - nbar
            for rest in self._list_rec(nbar, k - 1):
                yield (n_diff,) + rest

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: IV = IntegerVectors(2, 3)
            sage: list(IV)
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
            sage: list(IntegerVectors(3, 0))
            []
            sage: list(IntegerVectors(3, 1))
            [[3]]
            sage: list(IntegerVectors(0, 1))
            [[0]]
            sage: list(IntegerVectors(0, 2))
            [[0, 0]]
            sage: list(IntegerVectors(2, 2))
            [[2, 0], [1, 1], [0, 2]]
            sage: IntegerVectors(0, 0).list()
            [[]]
            sage: IntegerVectors(1, 0).list()
            []
            sage: IntegerVectors(0, 1).list()
            [[0]]
            sage: IntegerVectors(2, 2).list()
            [[2, 0], [1, 1], [0, 2]]
            sage: IntegerVectors(-1,0).list()
            []
            sage: IntegerVectors(-1,2).list()
            []
        """
        if self.n < 0:
            return

        if not self.k:
            if not self.n:
                yield self.element_class(self, [], check=False)
            return
        elif self.k == 1:
            yield self.element_class(self, [self.n], check=False)
            return

        for nbar in range(self.n+1):
            n = self.n - nbar
            for rest in integer_vectors_nk_fast_iter(nbar, self.k - 1):
                yield self.element_class(self, [n] + rest, check=False)

    def _repr_(self):
        """
        TESTS::

            sage: IV = IntegerVectors(2,3)
            sage: IV
            Integer vectors of length 3 that sum to 2
        """
        return "Integer vectors of length {} that sum to {}".format(self.k,
                                                                    self.n)

    def __contains__(self, x):
        """
        TESTS::

            sage: IV = IntegerVectors(2, 3)
            sage: all(i in IV for i in IV)
            True
            sage: [0,1,2] in IV
            False
            sage: [2.0, 0, 0] in IV
            True
            sage: [0,1,0,1] in IV
            False
            sage: [0,1,1] in IV
            True
            sage: [-1,2,1] in IV
            False

            sage: [0] in IntegerVectors(0, 1)
            True
            sage: [] in IntegerVectors(0, 0)
            True
            sage: [] in IntegerVectors(0, 1)
            False
            sage: [] in IntegerVectors(1, 0)
            False
            sage: [3] in IntegerVectors(2, 1)
            False
            sage: [3] in IntegerVectors(3, 1)
            True
            sage: [3,2,2,1] in IntegerVectors(9, 5)
            False
            sage: [3,2,2,1] in IntegerVectors(8, 5)
            False
            sage: [3,2,2,1] in IntegerVectors(8, 4)
            True

        Check :issue:`34510`::

            sage: IV33 = IntegerVectors(n=3, k=3)
            sage: IV33([0])
            Traceback (most recent call last):
            ...
            ValueError: [0] doesn't satisfy correct constraints
        """
        if not IntegerVectors.__contains__(self, x):
            return False

        if len(x) != self.k:
            return False

        if sum(x) != self.n:
            return False

        if len(x) > 0 and min(x) < 0:
            return False

        return True

    def rank(self, x):
        """
        Return the rank of a given element.

        INPUT:

        - ``x`` -- list with ``sum(x) == n`` and ``len(x) == k``

        TESTS::

            sage: IV = IntegerVectors(4,5)
            sage: list(range(IV.cardinality())) == [IV.rank(x) for x in IV]
            True
        """
        if x not in self:
            raise ValueError("argument is not a member of IntegerVectors({},{})".format(self.n, self.k))

        k, s, r = self.k, 0, 0
        for i in range(k - 1):
            s += x[k - 1 - i]
            r += binomial(s + i, i + 1)
        return r

    def unrank(self, x):
        """
        Return the element at given rank x.

        INPUT:

        - ``x`` -- integer such that ``x < self.cardinality()``

        EXAMPLES::

            sage: IntegerVectors(4,5).unrank(30)
            [1, 0, 1, 0, 2]
            sage: IntegerVectors(2,3).unrank(5)
            [0, 0, 2]
        """
        if x >= self.cardinality():
            raise IndexError(f"Index {x} is out of range for the IntegerVector.")
        rtn = [0]*self.k
        rtn[0] = self.n
        return IntegerVectors._unrank_helper(self, x, rtn)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: IntegerVectors(3,5).cardinality()
            35
            sage: IntegerVectors(99, 3).cardinality()
            5050
            sage: IntegerVectors(10^9 - 1, 3).cardinality()
            500000000500000000
        """
        n, k = self.n, self.k
        return Integer(binomial(n + k - 1, n))


class IntegerVectors_nnondescents(UniqueRepresentation, IntegerVectors):
    r"""
    Integer vectors graded by two parameters.

    The grading parameters on the integer vector `v` are:

    - ``n`` -- the sum of the parts of `v`

    - ``c`` -- the non descents composition of `v`

    In other words: the length of `v` equals `c_1 + \cdots + c_k`, and `v`
    is decreasing in the consecutive blocs of length `c_1, \ldots, c_k`,

    INPUT:

    - ``n`` -- the positive integer `n`
    - ``comp`` -- the composition `c`

    Those are the integer vectors of sum `n` that are lexicographically
    maximal (for the natural left-to-right reading) in their orbit by the
    Young subgroup `S_{c_1} \times \cdots \times S_{c_k}`. In particular,
    they form a set of orbit representative of integer vectors with
    respect to this Young subgroup.
    """
    @staticmethod
    def __classcall_private__(cls, n, comp):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: IntegerVectors(4, [2,1]) is IntegerVectors(int(4), (2,1))
            True
        """
        return super().__classcall__(cls, n, tuple(comp))

    def __init__(self, n, comp):
        """
        EXAMPLES::

            sage: IV = IntegerVectors(4, [2])
            sage: TestSuite(IV).run()
        """
        self.n = n
        self.comp = comp
        IntegerVectors.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: IntegerVectors(4, [2])
            Integer vectors of 4 with non-descents composition [2]
        """
        return "Integer vectors of {} with non-descents composition {}".format(self.n, list(self.comp))

    def __iter__(self):
        """
        TESTS::

            sage: IntegerVectors(0, []).list()
            [[]]
            sage: IntegerVectors(5, []).list()
            []
            sage: IntegerVectors(0, [1]).list()
            [[0]]
            sage: IntegerVectors(4, [1]).list()
            [[4]]
            sage: IntegerVectors(4, [2]).list()
            [[4, 0], [3, 1], [2, 2]]
            sage: IntegerVectors(4, [2,2]).list()
            [[4, 0, 0, 0],
             [3, 1, 0, 0],
             [2, 2, 0, 0],
             [3, 0, 1, 0],
             [2, 1, 1, 0],
             [2, 0, 2, 0],
             [2, 0, 1, 1],
             [1, 1, 2, 0],
             [1, 1, 1, 1],
             [1, 0, 3, 0],
             [1, 0, 2, 1],
             [0, 0, 4, 0],
             [0, 0, 3, 1],
             [0, 0, 2, 2]]
            sage: IntegerVectors(5, [1,1,1]).list()
            [[5, 0, 0],
             [4, 1, 0],
             [4, 0, 1],
             [3, 2, 0],
             [3, 1, 1],
             [3, 0, 2],
             [2, 3, 0],
             [2, 2, 1],
             [2, 1, 2],
             [2, 0, 3],
             [1, 4, 0],
             [1, 3, 1],
             [1, 2, 2],
             [1, 1, 3],
             [1, 0, 4],
             [0, 5, 0],
             [0, 4, 1],
             [0, 3, 2],
             [0, 2, 3],
             [0, 1, 4],
             [0, 0, 5]]
            sage: IntegerVectors(0, [2,3]).list()
            [[0, 0, 0, 0, 0]]
        """
        for iv in IntegerVectors(self.n, len(self.comp)):
            blocks = [IntegerVectors(iv[i], val, max_slope=0).list()
                      for i, val in enumerate(self.comp)]
            for parts in product(*blocks):
                res = []
                for part in parts:
                    res += part
                yield self.element_class(self, res, check=False)


class IntegerVectorsConstraints(IntegerVectors):
    """
    Class of integer vectors subject to various constraints.
    """

    def __init__(self, n=None, k=None, **constraints):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: TestSuite(IntegerVectors(min_slope=0)).run()
            sage: TestSuite(IntegerVectors(3, max_slope=0)).run()
            sage: TestSuite(IntegerVectors(3, max_length=4)).run()
            sage: TestSuite(IntegerVectors(k=2, max_part=4)).run()
            sage: TestSuite(IntegerVectors(k=2, min_part=2, max_part=4)).run()
            sage: TestSuite(IntegerVectors(3, 2, max_slope=0)).run()
        """
        self.n = n
        self.k = k
        if k is not None and self.k >= 0:
            constraints['length'] = self.k
        if 'outer' in constraints:
            constraints['ceiling'] = constraints['outer']
            del constraints['outer']
        if 'inner' in constraints:
            constraints['floor'] = constraints['inner']
            del constraints['inner']
        self.constraints = constraints

        if n is not None:
            if k is not None or 'max_length' in constraints:
                category = FiniteEnumeratedSets()
            else:
                category = EnumeratedSets()
        elif k is not None and 'max_part' in constraints:  # n is None
            category = FiniteEnumeratedSets()
        else:
            category = EnumeratedSets()
        IntegerVectors.__init__(self, category=category)  # placeholder category

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: IntegerVectors(min_slope=0)
            Integer vectors with constraints: min_slope=0

            sage: IntegerVectors(3, max_length=2)
            Integer vectors that sum to 3 with constraints: max_length=2

            sage: IntegerVectors(2, 3, min_slope=0)
            Integer vectors that sum to 2 with constraints: length=3, min_slope=0
        """
        if self.n is not None:
            base = "Integer vectors that sum to {} with constraints: ".format(self.n)
        else:
            base = "Integer vectors with constraints: "
        return base + ", ".join("{}={}".format(key, self.constraints[key])
                                for key in sorted(self.constraints))

    def __eq__(self, rhs):
        """
        EXAMPLES::

            sage: IntegerVectors(min_slope=0) == IntegerVectors(min_slope=0)
            True
            sage: IntegerVectors(2, min_slope=0) == IntegerVectors(2, min_slope=0)
            True
            sage: IntegerVectors(2, 3, min_slope=0) == IntegerVectors(2, 3, min_slope=0)
            True
        """
        if isinstance(rhs, IntegerVectorsConstraints):
            return self.n == rhs.n and self.k == rhs.k and self.constraints == rhs.constraints
        return False

    def __ne__(self, rhs):
        """
        EXAMPLES::

            sage: IntegerVectors(min_slope=0) != IntegerVectors(min_slope=3)
            True
        """
        return not self.__eq__(rhs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: hash(IntegerVectors(min_slope=0)) == hash(IntegerVectors(min_slope=0))
            True
            sage: hash(IntegerVectors(2, min_slope=0)) == hash(IntegerVectors(2, min_slope=0))
            True
            sage: hash(IntegerVectors(2, 3, min_slope=0)) == hash(IntegerVectors(2, 3, min_slope=0))
            True
            sage: hash(IntegerVectors(min_slope=0)) != hash(IntegerVectors(min_slope=3))
            True
        """
        return hash((self.n, self.k, tuple(self.constraints.items())))

    def __contains__(self, x):
        """
        TESTS::

            sage: [3,2,2,1] in IntegerVectors(8, 4, min_part=1)                         # needs sage.combinat
            True
            sage: [3,2,2,1] in IntegerVectors(8, 4, min_part=2)                         # needs sage.combinat
            False

            sage: [0,3,0,1,2] in IntegerVectors(6, max_length=3)                        # needs sage.combinat
            False
        """
        if isinstance(x, IntegerVector) and x.parent() is self:
            return True

        if not IntegerVectors.__contains__(self, x):
            return False

        if self.k is not None and len(x) != self.k:
            return False

        if self.n is not None and sum(x) != self.n:
            return False

        from sage.combinat.misc import check_integer_list_constraints
        return check_integer_list_constraints(x, singleton=True, **self.constraints)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: IntegerVectors(3, 3, min_part=1).cardinality()
            1
            sage: IntegerVectors(5, 3, min_part=1).cardinality()
            6
            sage: IntegerVectors(13, 4, max_part=4).cardinality()
            20
            sage: IntegerVectors(k=4, max_part=3).cardinality()
            256
            sage: IntegerVectors(k=3, min_part=2, max_part=4).cardinality()
            27
            sage: IntegerVectors(13, 4, min_part=2, max_part=4).cardinality()
            16
        """
        if self.k is None:
            if self.n is None:
                return PlusInfinity()
            if ('max_length' not in self.constraints
                    and self.constraints.get('min_part', 0) <= 0):
                return PlusInfinity()
        elif ('max_part' in self.constraints
                and self.constraints['max_part'] != PlusInfinity()):
            if (self.n is None and len(self.constraints) == 2
                    and 'min_part' in self.constraints
                    and self.constraints['min_part'] >= 0):
                num = self.constraints['max_part'] - self.constraints['min_part'] + 1
                return Integer(num ** self.k)
            if len(self.constraints) == 1:
                m = self.constraints['max_part']
                if self.n is None:
                    return Integer((m + 1) ** self.k)
                if m >= self.n:
                    return Integer(binomial(self.n + self.k - 1, self.n))
                # do by inclusion / exclusion on the number
                # i of parts greater than m
                n, k = self.n, self.k
                return Integer(sum(
                    (-1)**i * binomial(n + k - 1 - i * (m + 1), k - 1)
                    * binomial(k, i) for i in range(self.n // (m + 1) + 1)))
        return ZZ.sum(ZZ.one() for x in self)

    def __iter__(self):
        """
        EXAMPLES::

            sage: IntegerVectors(-1, 0, min_part=1).list()
            []
            sage: IntegerVectors(-1, 2, min_part=1).list()
            []
            sage: IntegerVectors(0, 0, min_part=1).list()
            [[]]
            sage: IntegerVectors(3, 0, min_part=1).list()
            []
            sage: IntegerVectors(0, 1, min_part=1).list()
            []
            sage: IntegerVectors(2, 2, min_part=1).list()
            [[1, 1]]
            sage: IntegerVectors(2, 3, min_part=1).list()
            []
            sage: IntegerVectors(4, 2, min_part=1).list()
            [[3, 1], [2, 2], [1, 3]]

        ::

            sage: IntegerVectors(0, 3, outer=[0,0,0]).list()
            [[0, 0, 0]]
            sage: IntegerVectors(1, 3, outer=[0,0,0]).list()
            []
            sage: IntegerVectors(2, 3, outer=[0,2,0]).list()
            [[0, 2, 0]]
            sage: IntegerVectors(2, 3, outer=[1,2,1]).list()
            [[1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1]]
            sage: IntegerVectors(2, 3, outer=[1,1,1]).list()
            [[1, 1, 0], [1, 0, 1], [0, 1, 1]]
            sage: IntegerVectors(2, 5, outer=[1,1,1,1,1]).list()
            [[1, 1, 0, 0, 0],
             [1, 0, 1, 0, 0],
             [1, 0, 0, 1, 0],
             [1, 0, 0, 0, 1],
             [0, 1, 1, 0, 0],
             [0, 1, 0, 1, 0],
             [0, 1, 0, 0, 1],
             [0, 0, 1, 1, 0],
             [0, 0, 1, 0, 1],
             [0, 0, 0, 1, 1]]

        ::

            sage: iv = [ IntegerVectors(n, k) for n in range(-2, 7) for k in range(7) ]
            sage: all(map(lambda x: x.cardinality() == len(x.list()), iv))
            True
            sage: essai = [[1,1,1], [2,5,6], [6,5,2]]
            sage: iv = [ IntegerVectors(x[0], x[1], max_part=x[2]-1) for x in essai ]
            sage: all(map(lambda x: x.cardinality() == len(x.list()), iv))
            True
        """
        from sage.combinat.integer_lists import IntegerListsLex

        if self.n is None:
            if self.k is not None and 'max_part' in self.constraints:
                n_list = range((self.constraints['max_part'] + 1) * self.k)
            else:
                n_list = NN
        else:
            n_list = [self.n]
        for n in n_list:
            for x in IntegerListsLex(n, check=False, **self.constraints):
                yield self.element_class(self, x, check=False)


def integer_vectors_nk_fast_iter(n, k):
    """
    A fast iterator for integer vectors of ``n`` of length ``k`` which
    yields Python lists filled with Sage Integers.

    EXAMPLES::

        sage: from sage.combinat.integer_vector import integer_vectors_nk_fast_iter
        sage: list(integer_vectors_nk_fast_iter(3, 2))
        [[3, 0], [2, 1], [1, 2], [0, 3]]
        sage: list(integer_vectors_nk_fast_iter(2, 2))
        [[2, 0], [1, 1], [0, 2]]
        sage: list(integer_vectors_nk_fast_iter(1, 2))
        [[1, 0], [0, 1]]

    We check some corner cases::

        sage: list(integer_vectors_nk_fast_iter(5, 1))
        [[5]]
        sage: list(integer_vectors_nk_fast_iter(1, 1))
        [[1]]
        sage: list(integer_vectors_nk_fast_iter(2, 0))
        []
        sage: list(integer_vectors_nk_fast_iter(0, 2))
        [[0, 0]]
        sage: list(integer_vectors_nk_fast_iter(0, 0))
        [[]]
    """
    # "bad" input
    if n < 0 or k < 0:
        return

    # Check some corner cases first
    if not k:
        if not n:
            yield []
        return
    n = Integer(n)
    if k == 1:
        yield [n]
        return

    zero = ZZ.zero()
    one = ZZ.one()
    k = int(k)

    pos = 0  # Current position
    rem = zero  # Amount remaining
    cur = [n] + [zero] * (k - 1)  # Current list
    yield list(cur)
    while pos >= 0:
        if not cur[pos]:
            pos -= 1
            continue
        cur[pos] -= one
        rem += one
        if not rem:
            yield list(cur)
        elif pos == k - 2:
            cur[pos + 1] = rem
            yield list(cur)
            cur[pos + 1] = zero
        else:
            pos += 1
            cur[pos] = rem  # Guaranteed to be at least 1
            rem = zero
            yield list(cur)


# October 2012: fixing outdated pickles which use classes being deprecated
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.integer_vector', 'IntegerVectors_nconstraints', IntegerVectorsConstraints)
register_unpickle_override('sage.combinat.integer_vector', 'IntegerVectors_nkconstraints', IntegerVectorsConstraints)
