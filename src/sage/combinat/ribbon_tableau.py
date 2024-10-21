# sage_setup: distribution = sagemath-combinat
r"""
Ribbon Tableaux
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
import functools

from sage.structure.parent import Parent
from sage.structure.element import parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_cat import Sets
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.combinat.combinat import CombinatorialElement
from sage.combinat.skew_partition import SkewPartition, SkewPartitions
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux, SemistandardSkewTableaux
from sage.combinat.tableau import Tableaux
from sage.combinat.partition import Partition, _Partitions
from sage.combinat.permutation import to_standard
from sage.misc.persist import register_unpickle_override
from . import permutation


class RibbonTableau(SkewTableau):
    r"""
    A ribbon tableau.

    A ribbon is a connected skew shape which does not contain any
    `2 \times 2` boxes.  A ribbon tableau is a skew tableau
    whose shape is partitioned into ribbons, each of which is filled
    with identical entries.

    EXAMPLES::

        sage: rt = RibbonTableau([[None, 1],[2,3]]); rt
        [[None, 1], [2, 3]]
        sage: rt.inner_shape()
        [1]
        sage: rt.outer_shape()
        [2, 2]

        sage: rt = RibbonTableau([[None, None, 0, 0, 0], [None, 0, 0, 2], [1, 0, 1]]); rt.pp()
          .  .  0  0  0
          .  0  0  2
          1  0  1

    In the previous example, each ribbon is uniquely determined by a
    nonzero entry.  The 0 entries are used to fill in the rest of the
    skew shape.

    .. NOTE::

        Sanity checks are not performed; lists can contain any object.

    ::

        sage: RibbonTableau(expr=[[1,1],[[5],[3,4],[1,2]]])
        [[None, 1, 2], [None, 3, 4], [5]]

    TESTS::

        sage: RibbonTableau([[0, 0, 3, 0], [1, 1, 0], [2, 0, 4]]).evaluation()
        [2, 1, 1, 1]
    """
    # The following method is private and will only get called
    # when calling RibbonTableau() directly, and not via element_class
    @staticmethod
    def __classcall_private__(cls, rt=None, expr=None):
        """
        Return a ribbon tableau object.

        EXAMPLES::

            sage: rt = RibbonTableau([[None, 1],[2,3]]); rt
            [[None, 1], [2, 3]]
            sage: TestSuite(rt).run()
        """
        if expr is not None:
            return RibbonTableaux().from_expr(expr)

        try:
            rt = [tuple(row) for row in rt]
        except TypeError:
            raise TypeError("each element of the ribbon tableau must be an iterable")
        if not all(row for row in rt):
            raise TypeError("a ribbon tableau cannot have empty rows")
        # calls the inherited __init__ method (of SkewTableau )
        return RibbonTableaux()(rt)

    def length(self):
        """
        Return the length of the ribbons into a ribbon tableau.

        EXAMPLES::

            sage: RibbonTableau([[None, 1],[2,3]]).length()
            1
            sage: RibbonTableau([[1,0],[2,0]]).length()
            2
        """
        if self.to_expr() == [[], []]:
            return 0

        tableau = self.to_expr()[1]
        l = 0
        t = 0
        for k in range(len(tableau)):
            t += len([x for x in tableau[k] if x is not None and x > -1])
            l += len([x for x in tableau[k] if x is not None and x > 0])

        if l == 0:
            return t
        else:
            return t // l

    def to_word(self):
        """
        Return a word obtained from a row reading of ``self``.

        .. WARNING::

            Unlike the ``to_word`` method on skew tableaux (which are a
            superclass of this), this method does not filter out
            ``None`` entries.

        EXAMPLES::

            sage: R = RibbonTableau([[0, 0, 3, 0], [1, 1, 0], [2, 0, 4]])
            sage: R.to_word()
            word: 2041100030
        """
        from sage.combinat.words.word import Word
        return Word([letter for row in reversed(self) for letter in row])

# ===================
#   Ribbon Tableaux
# ===================


class RibbonTableaux(UniqueRepresentation, Parent):
    r"""
    Ribbon tableaux.

    A ribbon tableau is a skew tableau whose skew shape ``shape`` is
    tiled by ribbons of length ``length``. The weight ``weight`` is
    calculated from the labels on the ribbons.

    .. NOTE::

        Here we impose the condition that the ribbon tableaux are semistandard.

    INPUT(Optional):

    - ``shape`` -- skew shape as a list of lists or an object of type
      SkewPartition

    - ``length`` -- integer; ``shape`` is partitioned into ribbons of
      length ``length``

    - ``weight`` -- list of integers; computed from the values of
      nonzero entries labeling the ribbons

    EXAMPLES::

        sage: RibbonTableaux([[2,1],[]], [1,1,1], 1)
        Ribbon tableaux of shape [2, 1] / [] and weight [1, 1, 1] with 1-ribbons

        sage: R = RibbonTableaux([[5,4,3],[2,1]], [2,1], 3)
        sage: for i in R: i.pp(); print("\n")
          .  .  0  0  0
          .  0  0  2
          1  0  1
        <BLANKLINE>
          .  .  1  0  0
          .  0  0  0
          1  0  2
        <BLANKLINE>
          .  .  0  0  0
          .  1  0  1
          2  0  0
        <BLANKLINE>

    REFERENCES:

    .. [vanLeeuwen91] Marc. A. A. van Leeuwen, *Edge sequences, ribbon tableaux,
       and an action of affine permutations*. Europe J. Combinatorics. **20**
       (1999). http://wwwmathlabo.univ-poitiers.fr/~maavl/pdf/edgeseqs.pdf
    """
    @staticmethod
    def __classcall_private__(cls, shape=None, weight=None, length=None):
        """
        Return the correct parent object.

        EXAMPLES::

            sage: R = RibbonTableaux([[2,1],[]],[1,1,1],1)
            sage: R2 = RibbonTableaux(SkewPartition([[2,1],[]]),(1,1,1),1)
            sage: R is R2
            True
        """
        if shape is None and weight is None and length is None:
            return super().__classcall__(cls)

        return RibbonTableaux_shape_weight_length(shape, weight, length)

    def __init__(self):
        """
        EXAMPLES::

            sage: R = RibbonTableaux()
            sage: TestSuite(R).run()
        """
        Parent.__init__(self, category=Sets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RibbonTableaux()
            Ribbon tableaux
        """
        return "Ribbon tableaux"

    def _element_constructor_(self, rt):
        """
        Construct an element of ``self`` from ``rt``.

        EXAMPLES::

            sage: R = RibbonTableaux()
            sage: elt = R([[0, 0, 3, 0], [1, 1, 0], [2, 0, 4]]); elt
            [[0, 0, 3, 0], [1, 1, 0], [2, 0, 4]]
            sage: elt.parent() is R
            True
        """
        return self.element_class(self, rt)

    def from_expr(self, l):
        """
        Return a :class:`RibbonTableau` from a MuPAD-Combinat expr for a skew
        tableau. The first list in ``expr`` is the inner shape of the skew
        tableau. The second list are the entries in the rows of the skew
        tableau from bottom to top.

        Provided primarily for compatibility with MuPAD-Combinat.

        EXAMPLES::

            sage: RibbonTableaux().from_expr([[1,1],[[5],[3,4],[1,2]]])
            [[None, 1, 2], [None, 3, 4], [5]]
        """
        return self.element_class(self, SkewTableaux().from_expr(l))

    Element = RibbonTableau
    options = Tableaux.options


class RibbonTableaux_shape_weight_length(RibbonTableaux):
    """
    Ribbon tableaux of a given shape, weight, and length.
    """
    @staticmethod
    def __classcall_private__(cls, shape, weight, length):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: R = RibbonTableaux([[2,1],[]],[1,1,1],1)
            sage: R2 = RibbonTableaux(SkewPartition([[2,1],[]]),(1,1,1),1)
            sage: R is R2
            True
        """
        if shape in _Partitions:
            shape = _Partitions(shape)
            shape = SkewPartition([shape, shape.core(length)])
        else:
            shape = SkewPartition(shape)

        if shape.size() != length * sum(weight):
            raise ValueError("incompatible shape and weight")

        return super().__classcall__(cls, shape, tuple(weight), length)

    def __init__(self, shape, weight, length):
        """
        EXAMPLES::

            sage: R = RibbonTableaux([[2,1],[]],[1,1,1],1)
            sage: TestSuite(R).run()
        """
        self._shape = shape
        self._weight = weight
        self._length = length
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __iter__(self):
        """
        EXAMPLES::

            sage: RibbonTableaux([[2,1],[]],[1,1,1],1).list()
            [[[1, 3], [2]], [[1, 2], [3]]]
            sage: RibbonTableaux([[2,2],[]],[1,1],2).list()
            [[[0, 0], [1, 2]], [[1, 0], [2, 0]]]
        """
        for x in graph_implementation_rec(self._shape, self._weight, self._length, list_rec):
            yield self.from_expr(x)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RibbonTableaux([[2,1],[]], [1,1,1], 1)
            Ribbon tableaux of shape [2, 1] / [] and weight [1, 1, 1] with 1-ribbons
        """
        return "Ribbon tableaux of shape %s and weight %s with %s-ribbons" % (repr(self._shape), list(self._weight), self._length)

    def __contains__(self, x):
        """
        Note that this just checks to see if ``x`` appears in ``self``.

        This should be improved to provide actual checking.

        EXAMPLES::

            sage: r = RibbonTableaux([[2,2],[]],[1,1],2)
            sage: [[0, 0], [1, 2]] in r
            True
            sage: [[1, 0], [2, 0]] in r
            True
            sage: [[0, 1], [2, 0]] in r
            False
        """
        try:
            x = RibbonTableau(x)
        except (ValueError, TypeError):
            return False
        return x in self.list()

        # return x.is_ribbon() and x.shape() == self._shape \
        # and tuple(x.weight()) == self._weight and x in list(self)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: RibbonTableaux([[2,1],[]],[1,1,1],1).cardinality()
            2
            sage: RibbonTableaux([[2,2],[]],[1,1],2).cardinality()
            2
            sage: RibbonTableaux([[4,3,3],[]],[2,1,1,1],2).cardinality()
            5

        TESTS::

            sage: RibbonTableaux([6,6,6], [4,2], 3).cardinality()
            6
            sage: RibbonTableaux([3,3,3,2,1], [3,1], 3).cardinality()
            1
            sage: RibbonTableaux([3,3,3,2,1], [2,2], 3).cardinality()
            2
            sage: RibbonTableaux([3,3,3,2,1], [2,1,1], 3).cardinality()
            5
            sage: RibbonTableaux([3,3,3,2,1], [1,1,1,1], 3).cardinality()
            12
            sage: RibbonTableaux([5,4,3,2,1], [2,2,1], 3).cardinality()
            10

        ::

            sage: RibbonTableaux([8,7,6,5,1,1], [3,2,2,1], 3).cardinality()
            85
            sage: RibbonTableaux([5,4,3,2,1,1,1], [2,2,1], 3).cardinality()
            10

        ::

            sage: RibbonTableaux([7,7,7,2,1,1], [3,2,0,1,1], 3).cardinality()
            25

        Weights with some zeros in the middle and end::

            sage: RibbonTableaux([3,3,3], [0,1,0,2,0], 3).cardinality()
            3
            sage: RibbonTableaux([3,3,3], [1,0,1,0,1,0,0,0], 3).cardinality()
            6
        """
        # Strip zeros for graph_implementation_rec
        wt = [i for i in self._weight if i != 0]
        return Integer(graph_implementation_rec(self._shape, wt, self._length, count_rec)[0])


def insertion_tableau(skp, perm, evaluation, tableau, length):
    """
    INPUT:

    - ``skp`` -- skew partitions

    - ``perm, evaluation`` -- nonnegative integers

    - ``tableau`` -- skew tableau

    - ``length`` -- integer

    TESTS::

        sage: from sage.combinat.ribbon_tableau import insertion_tableau
        sage: insertion_tableau([[1], []], [1], 1, [[], []], 1)
        [[], [[1]]]
        sage: insertion_tableau([[2, 1], []], [1, 1], 2, [[], [[1]]], 1)
        [[], [[2], [1, 2]]]
        sage: insertion_tableau([[2, 1], []], [0, 0], 3, [[], [[2], [1, 2]]], 1)
        [[], [[2], [1, 2]]]
        sage: insertion_tableau([[1, 1], []], [1], 2, [[], [[1]]], 1)
        [[], [[2], [1]]]
        sage: insertion_tableau([[2], []], [0, 1], 2, [[], [[1]]], 1)
        [[], [[1, 2]]]
        sage: insertion_tableau([[2, 1], []], [0, 1], 3, [[], [[2], [1]]], 1)
        [[], [[2], [1, 3]]]
        sage: insertion_tableau([[1, 1], []], [2], 1, [[], []], 2)
        [[], [[1], [0]]]
        sage: insertion_tableau([[2], []], [2, 0], 1, [[], []], 2)
        [[], [[1, 0]]]
        sage: insertion_tableau([[2, 2], []], [0, 2], 2, [[], [[1], [0]]], 2)
        [[], [[1, 2], [0, 0]]]
        sage: insertion_tableau([[2, 2], []], [2, 0], 2, [[], [[1, 0]]], 2)
        [[], [[2, 0], [1, 0]]]
        sage: insertion_tableau([[2, 2], [1]], [3, 0], 1, [[], []], 3)
        [[1], [[1, 0], [0]]]
    """
    psave = Partition(skp[1])
    partc = skp[1] + [0] * (len(skp[0]) - len(skp[1]))

    tableau = SkewTableau(expr=tableau).to_expr()[1]

    for k in range(len(tableau)):
        tableau[-(k + 1)] += [0] * (skp[0][k] - partc[k] - len(tableau[-(k + 1)]))

    # We construct a tableau from the southwest corner to the northeast one
    tableau = [[0] * (skp[0][k] - partc[k])
               for k in reversed(range(len(tableau), len(skp[0])))] + tableau

    tableau = SkewTableaux().from_expr([skp[1], tableau]).conjugate()
    tableau = tableau.to_expr()[1]

    skp = SkewPartition(skp).conjugate().to_list()
    skp[1].extend([0] * (len(skp[0]) - len(skp[1])))

    if len(perm) > len(skp[0]):
        return None

    lp = len(perm)
    for k in range(lp):
        if perm[-(k + 1)] != 0:
            tableau[len(tableau) - lp + k][skp[0][lp - (k + 1)] - skp[1][lp - (k + 1)] - 1] = evaluation

    return SkewTableau(expr=[psave.conjugate(), tableau]).conjugate().to_expr()


def count_rec(nexts, current, part, weight, length):
    """
    INPUT:

    - ``nexts, current, part`` -- skew partitions

    - ``weight`` -- nonnegative integer list

    - ``length`` -- integer

    TESTS::

        sage: from sage.combinat.ribbon_tableau import count_rec
        sage: count_rec([], [], [[2, 1, 1], []], [2], 2)
        [0]
        sage: count_rec([[0], [1]], [[[2, 1, 1], [0, 0, 2, 0]], [[4], [2, 0, 0, 0]]], [[4, 1, 1], []], [2, 1], 2)
        [1]
        sage: count_rec([], [[[], [2, 2]]], [[2, 2], []], [2], 2)
        [1]
        sage: count_rec([], [[[], [2, 0, 2, 0]]], [[4], []], [2], 2)
        [1]
        sage: count_rec([[1], [1]], [[[2, 2], [0, 0, 2, 0]], [[4], [2, 0, 0, 0]]], [[4, 2], []], [2, 1], 2)
        [2]
        sage: count_rec([[1], [1], [2]], [[[2, 2, 2], [0, 0, 2, 0]], [[4, 1, 1], [0, 2, 0, 0]], [[4, 2], [2, 0, 0, 0]]], [[4, 2, 2], []], [2, 1, 1], 2)
        [4]
        sage: count_rec([[4], [1]], [[[4, 2, 2], [0, 0, 2, 0]], [[4, 3, 1], [0, 2, 0, 0]]], [[4, 3, 3], []], [2, 1, 1, 1], 2)
        [5]
    """
    if not current:
        return [0]
    if nexts:
        return [sum(j for i in nexts for j in i)]
    else:
        return [len(current)]


def list_rec(nexts, current, part, weight, length):
    """
    INPUT:

    - ``nexts, current, part`` -- skew partitions

    - ``weight`` -- nonnegative integer list

    - ``length`` -- integer

    TESTS::

        sage: from sage.combinat.ribbon_tableau import list_rec
        sage: list_rec([], [[[], [1]]], [[1], []], [1], 1)
        [[[], [[1]]]]
        sage: list_rec([[[[], [[1]]]]], [[[1], [1, 1]]], [[2, 1], []], [1, 2], 1)
        [[[], [[2], [1, 2]]]]
        sage: list_rec([], [[[1], [3, 0]]], [[2, 2], [1]], [1], 3)
        [[[1], [[1, 0], [0]]]]
        sage: list_rec([[[[], [[2]]]]], [[[1], [1, 1]]], [[2, 1], []], [0, 1, 2], 1)
        [[[], [[3], [2, 3]]]]
        sage: list_rec([], [[[], [2]]], [[1, 1], []], [1], 2)
        [[[], [[1], [0]]]]
        sage: list_rec([], [[[], [2, 0]]], [[2], []], [1], 2)
        [[[], [[1, 0]]]]
        sage: list_rec([[[[], [[1], [0]]]], [[[], [[1, 0]]]]], [[[1, 1], [0, 2]], [[2], [2, 0]]], [[2, 2], []], [1, 1], 2)
        [[[], [[1, 2], [0, 0]]], [[], [[2, 0], [1, 0]]]]
        sage: list_rec([], [[[], [2, 2]]], [[2, 2], []], [2], 2)
        [[[], [[1, 1], [0, 0]]]]
        sage: list_rec([], [[[], [1, 1]]], [[2], []], [2], 1)
        [[[], [[1, 1]]]]
        sage: list_rec([[[[], [[1, 1]]]]], [[[2], [1, 1]]], [[2, 2], []], [2, 2], 1)
        [[[], [[2, 2], [1, 1]]]]
    """
    if current == [] and nexts == [] and weight == []:
        return [[part[1], []]]

    # Test if the current nodes is not an empty node
    if not current:
        return []

    # Test if the current nodes drive us to new solutions
    if nexts:
        return [insertion_tableau(part, curr_i[1], len(weight),
                                  nexts_ij, length)
                for nexts_i, curr_i in zip(nexts, current)
                for nexts_ij in nexts_i]

    # The current nodes are at the bottom of the tree
    return [insertion_tableau(part, curr_i[1],
                              len(weight), [[], []], length)
            for curr_i in current]


# ===============================
#   Spin and Cospin Polynomials
# ===============================
def spin_rec(t, nexts, current, part, weight, length):
    """
    Routine used for constructing the spin polynomial.

    INPUT:

    - ``weight`` -- list of nonnegative integers

    - ``length`` -- the length of the ribbons we're tiling with

    - ``t`` -- the variable

    EXAMPLES::

        sage: from sage.combinat.ribbon_tableau import spin_rec
        sage: sp = SkewPartition
        sage: t = ZZ['t'].gen()
        sage: spin_rec(t, [], [[[], [3, 3]]], sp([[2, 2, 2], []]), [2], 3)
        [t^4]
        sage: spin_rec(t, [[0], [t^4]], [[[2, 1, 1, 1, 1], [0, 3]], [[2, 2, 2], [3, 0]]], sp([[2, 2, 2, 2, 1], []]), [2, 1], 3)
        [t^5]
        sage: spin_rec(t, [], [[[], [3, 3, 0]]], sp([[3, 3], []]), [2], 3)
        [t^2]
        sage: spin_rec(t, [[t^4], [t^3], [t^2]], [[[2, 2, 2], [0, 0, 3]], [[3, 2, 1], [0, 3, 0]], [[3, 3], [3, 0, 0]]], sp([[3, 3, 3], []]), [2, 1], 3)
        [t^6 + t^4 + t^2]
        sage: spin_rec(t, [[t^5], [t^4], [t^6 + t^4 + t^2]], [[[2, 2, 2, 2, 1], [0, 0, 3]], [[3, 3, 1, 1, 1], [0, 3, 0]], [[3, 3, 3], [3, 0, 0]]], sp([[3, 3, 3, 2, 1], []]), [2, 1, 1], 3)
        [2*t^7 + 2*t^5 + t^3]
    """
    if not current:
        return [parent(t).zero()]

    tmp = []
    partp = part[0].conjugate()
    ell = len(partp)

    # compute the contribution of the ribbons added at
    # the current node
    for val in current:
        perms = val[1]
        perm = [partp[i] + ell - (i + 1) - perms[i] for i in reversed(range(ell))]
        perm = to_standard(perm)
        tmp.append(weight[-1] * (length - 1) - perm.number_of_inversions())

    if nexts:
        return [sum(sum(t**tval * nval for nval in nexts[i])
                    for i, tval in enumerate(tmp))]
    return [sum(t**val for val in tmp)]


def spin_polynomial_square(part, weight, length):
    r"""
    Return the spin polynomial associated with ``part``, ``weight``, and
    ``length``, with the substitution `t \to t^2` made.

    EXAMPLES::

        sage: from sage.combinat.ribbon_tableau import spin_polynomial_square
        sage: spin_polynomial_square([6,6,6],[4,2],3)
        t^12 + t^10 + 2*t^8 + t^6 + t^4
        sage: spin_polynomial_square([6,6,6],[4,1,1],3)
        t^12 + 2*t^10 + 3*t^8 + 2*t^6 + t^4
        sage: spin_polynomial_square([3,3,3,2,1], [2,2], 3)
        t^7 + t^5
        sage: spin_polynomial_square([3,3,3,2,1], [2,1,1], 3)
        2*t^7 + 2*t^5 + t^3
        sage: spin_polynomial_square([3,3,3,2,1], [1,1,1,1], 3)
        3*t^7 + 5*t^5 + 3*t^3 + t
        sage: spin_polynomial_square([5,4,3,2,1,1,1], [2,2,1], 3)
        2*t^9 + 6*t^7 + 2*t^5
        sage: spin_polynomial_square([[6]*6, [3,3]], [4,4,2], 3)
        3*t^18 + 5*t^16 + 9*t^14 + 6*t^12 + 3*t^10
    """
    R = ZZ['t']

    if part in _Partitions:
        part = SkewPartition([part, _Partitions([])])
    elif part in SkewPartitions():
        part = SkewPartition(part)

    if part == [[], []] and not weight:
        return R.one()

    t = R.gen()
    return R(graph_implementation_rec(part, weight, length,
                                      functools.partial(spin_rec, t))[0])


def spin_polynomial(part, weight, length):
    """
    Return the spin polynomial associated to ``part``, ``weight``, and
    ``length``.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: from sage.combinat.ribbon_tableau import spin_polynomial
        sage: spin_polynomial([6,6,6],[4,2],3)
        t^6 + t^5 + 2*t^4 + t^3 + t^2
        sage: spin_polynomial([6,6,6],[4,1,1],3)
        t^6 + 2*t^5 + 3*t^4 + 2*t^3 + t^2
        sage: spin_polynomial([3,3,3,2,1], [2,2], 3)
        t^(7/2) + t^(5/2)
        sage: spin_polynomial([3,3,3,2,1], [2,1,1], 3)
        2*t^(7/2) + 2*t^(5/2) + t^(3/2)
        sage: spin_polynomial([3,3,3,2,1], [1,1,1,1], 3)
        3*t^(7/2) + 5*t^(5/2) + 3*t^(3/2) + sqrt(t)
        sage: spin_polynomial([5,4,3,2,1,1,1], [2,2,1], 3)
        2*t^(9/2) + 6*t^(7/2) + 2*t^(5/2)
        sage: spin_polynomial([[6]*6, [3,3]], [4,4,2], 3)
        3*t^9 + 5*t^8 + 9*t^7 + 6*t^6 + 3*t^5
    """
    from sage.symbolic.ring import SR
    sp = spin_polynomial_square(part, weight, length)
    t = SR.var('t')
    coeffs = sp.list()
    return sum(c * t**(ZZ(i) / 2) for i, c in enumerate(coeffs))


def cospin_polynomial(part, weight, length):
    """
    Return the cospin polynomial associated to ``part``, ``weight``, and
    ``length``.

    EXAMPLES::

        sage: from sage.combinat.ribbon_tableau import cospin_polynomial
        sage: cospin_polynomial([6,6,6],[4,2],3)
        t^4 + t^3 + 2*t^2 + t + 1
        sage: cospin_polynomial([3,3,3,2,1], [3,1], 3)
        1
        sage: cospin_polynomial([3,3,3,2,1], [2,2], 3)
        t + 1
        sage: cospin_polynomial([3,3,3,2,1], [2,1,1], 3)
        t^2 + 2*t + 2
        sage: cospin_polynomial([3,3,3,2,1], [1,1,1,1], 3)
        t^3 + 3*t^2 + 5*t + 3
        sage: cospin_polynomial([5,4,3,2,1,1,1], [2,2,1], 3)
        2*t^2 + 6*t + 2
        sage: cospin_polynomial([[6]*6, [3,3]], [4,4,2], 3)
        3*t^4 + 6*t^3 + 9*t^2 + 5*t + 3
    """
    R = ZZ['t']

    # The power in the spin polynomial are all half integers
    # or all integers.  Manipulation of expressions need to
    # separate cases
    sp = spin_polynomial_square(part, weight, length)
    if sp == 0:
        return R.zero()

    coeffs = [c for c in sp.list() if c]
    d = len(coeffs) - 1
    t = R.gen()
    return R(sum(c * t**(d - i) for i, c in enumerate(coeffs)))


#   //////////////////////////////////////////////////////////////////////////////////////////
#   // Generic function for driving into the graph of partitions coding all ribbons
#   // tableaux of a given shape and weight
#   //////////////////////////////////////////////////////////////////////////////////////////
#   // This function constructs the graph of the set of k-ribbon tableaux
#   // of a given skew shape and a given weight.
#   // The first argument is always a skew partition.
#   // In the case where the inner partition is empty, there is no branch without solutions.
#   // In the other cases, there is in average a lot of branches without solutions.
#   /////////////////////////////////////////////////////////////////////////////////////////


def graph_implementation_rec(skp, weight, length, function):
    """
    TESTS::

        sage: from sage.combinat.ribbon_tableau import graph_implementation_rec, list_rec
        sage: graph_implementation_rec(SkewPartition([[1], []]), [1], 1, list_rec)
        [[[], [[1]]]]
        sage: graph_implementation_rec(SkewPartition([[2, 1], []]), [1, 2], 1, list_rec)
        [[[], [[2], [1, 2]]]]
        sage: graph_implementation_rec(SkewPartition([[], []]), [0], 1, list_rec)
        [[[], []]]
    """
    if sum(weight) == 0:
        weight = []

    partp = skp[0].conjugate()
    ell = len(partp)
    outer = skp[1]
    outer_len = len(outer)

    # Some tests in order to know if the shape and the weight are compatible.
    if weight and weight[-1] <= len(partp):
        perms = permutation.Permutations([0] * (len(partp) - weight[-1]) + [length] * (weight[-1])).list()
    else:
        return function([], [], skp, weight, length)

    selection = []

    for j in range(len(perms)):
        retire = [(val + ell - (i + 1) - perms[j][i]) for i, val in enumerate(partp)]
        retire.sort(reverse=True)
        retire = [val - ell + (i + 1) for i, val in enumerate(retire)]

        if retire[-1] >= 0 and retire == sorted(retire, reverse=True):
            retire = Partition(retire).conjugate()

            # Cutting branches if the retired partition has a line strictly included into the inner one
            if len(retire) >= outer_len:
                append = True
                for k in range(outer_len):
                    if retire[k] - outer[k] < 0:
                        append = False
                        break
                if append:
                    selection.append([retire, perms[j]])

    # selection contains the list of current nodes

    if len(weight) == 1:
        return function([], selection, skp, weight, length)
    else:
        # The recursive calls permit us to construct the list of the sons
        # of all current nodes in selection
        a = [graph_implementation_rec([p[0], outer], weight[:-1], length, function)
             for p in selection]
        return function(a, selection, skp, weight, length)


class MultiSkewTableau(CombinatorialElement):
    """
    A multi skew tableau which is a tuple of skew tableaux.

    EXAMPLES::

        sage: s = MultiSkewTableau([ [[None,1],[2,3]], [[1,2],[2]] ])
        sage: s.size()
        6
        sage: s.weight()
        [2, 3, 1]
        sage: s.shape()
        [[2, 2] / [1], [2, 1] / []]

    TESTS::

        sage: mst = MultiSkewTableau([ [[None,1],[2,3]], [[1,2],[2]] ])
        sage: TestSuite(mst).run()
    """
    @staticmethod
    def __classcall_private__(cls, x):
        """
        Construct a multi skew tableau.

        EXAMPLES::

            sage: s = MultiSkewTableau([ [[None,1],[2,3]], [[1,2],[2]] ])
        """
        if isinstance(x, MultiSkewTableau):
            return x
        return MultiSkewTableaux()([SkewTableau(i) for i in x])

    def size(self):
        """
        Return the size of ``self``.

        This is the sum of the sizes of the skew
        tableaux in ``self``.

        EXAMPLES::

            sage: s = SemistandardSkewTableaux([[2,2],[1]]).list()
            sage: a = MultiSkewTableau([s[0],s[1],s[2]])
            sage: a.size()
            9
        """
        return sum(x.size() for x in self)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: s = SemistandardSkewTableaux([[2,2],[1]]).list()
            sage: a = MultiSkewTableau([s[0],s[1],s[2]])
            sage: a.weight()
            [5, 3, 1]
        """
        weights = [x.weight() for x in self]
        m = max([len(x) for x in weights])
        weight = [0] * m
        for w in weights:
            for i in range(len(w)):
                weight[i] += w[i]
        return weight

    def shape(self):
        """
        Return the shape of ``self``.

        EXAMPLES::

            sage: s = SemistandardSkewTableaux([[2,2],[1]]).list()
            sage: a = MultiSkewTableau([s[0],s[1],s[2]])
            sage: a.shape()
            [[2, 2] / [1], [2, 2] / [1], [2, 2] / [1]]
        """
        return [x.shape() for x in self]

    def inversion_pairs(self):
        """
        Return a list of the inversion pairs of ``self``.

        EXAMPLES::

            sage: s = MultiSkewTableau([ [[2,3],[5,5]], [[1,1],[3,3]], [[2],[6]] ])
            sage: s.inversion_pairs()
            [((0, (0, 0)), (1, (0, 0))),
             ((0, (1, 0)), (1, (0, 1))),
             ((0, (1, 1)), (1, (0, 0))),
             ((0, (1, 1)), (1, (1, 1))),
             ((0, (1, 1)), (2, (0, 0))),
             ((1, (0, 1)), (2, (0, 0))),
             ((1, (1, 1)), (2, (0, 0)))]
        """
        inv = []
        for k in range(len(self)):
            for b in self[k].cells():
                inv += self._inversion_pairs_from_position(k, b)
        return inv

    def inversions(self):
        """
        Return the number of inversion pairs of ``self``.

        EXAMPLES::

            sage: t1 = SkewTableau([[1]])
            sage: t2 = SkewTableau([[2]])
            sage: MultiSkewTableau([t1,t1]).inversions()
            0
            sage: MultiSkewTableau([t1,t2]).inversions()
            0
            sage: MultiSkewTableau([t2,t2]).inversions()
            0
            sage: MultiSkewTableau([t2,t1]).inversions()
            1
            sage: s = MultiSkewTableau([ [[2,3],[5,5]], [[1,1],[3,3]], [[2],[6]] ])
            sage: s.inversions()
            7
        """
        return len(self.inversion_pairs())

    def _inversion_pairs_from_position(self, k, ij):
        """
        Return the number of inversions at the cell position `(i,j)` in the
        ``k``-th tableaux in ``self``.

        EXAMPLES::

            sage: s = MultiSkewTableau([ [[2,3],[5,5]], [[1,1],[3,3]], [[2],[6]] ])
            sage: s._inversion_pairs_from_position(0, (1,1))
            [((0, (1, 1)), (1, (0, 0))),
             ((0, (1, 1)), (1, (1, 1))),
             ((0, (1, 1)), (2, (0, 0)))]
            sage: s._inversion_pairs_from_position(1, (0,1))
            [((1, (0, 1)), (2, (0, 0)))]
        """
        pk = k
        pi, pj = ij
        c = pi - pj
        value = self[pk][pi][pj]
        pk_cells = self[pk].cells_by_content(c)
        same_diagonal = [t.cells_by_content(c) for t in self[pk + 1:]]
        above_diagonal = [t.cells_by_content(c + 1) for t in self[pk + 1:]]

        res = []
        for i, j in pk_cells:
            if pi < i and value > self[pk][i][j]:
                res.append(((pk, (pi, pj)), (pk, (i, j))))
        for k in range(len(same_diagonal)):
            for i, j in same_diagonal[k]:
                if value > self[pk + k + 1][i][j]:
                    res.append(((pk, (pi, pj)), (pk + k + 1, (i, j))))
        for k in range(len(above_diagonal)):
            for i, j in above_diagonal[k]:
                if value < self[pk + k + 1][i][j]:
                    res.append(((pk, (pi, pj)), (pk + k + 1, (i, j))))
        return res


class MultiSkewTableaux(UniqueRepresentation, Parent):
    r"""
    Multiskew tableaux.
    """

    def __init__(self, category=None):
        """
        EXAMPLES::

            sage: R = MultiSkewTableaux()
            sage: TestSuite(R).run()
        """
        if category is None:
            category = Sets()
        Parent.__init__(self, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: MultiSkewTableaux()
            Multi Skew Tableaux tableaux
        """
        return "Multi Skew Tableaux tableaux"

    def _element_constructor_(self, rt):
        """
        Construct an element of ``self`` from ``rt``.

        EXAMPLES::

            sage: R = MultiSkewTableaux()
            sage: R([[[1, 1], [2]], [[None, 2], [3, 3]]])
            [[[1, 1], [2]], [[None, 2], [3, 3]]]
        """
        return self.element_class(self, rt)

    Element = MultiSkewTableau


class SemistandardMultiSkewTableaux(MultiSkewTableaux):
    """
    Semistandard multi skew tableaux.

    A multi skew tableau is a `k`-tuple of skew tableaux of
    given shape with a specified total weight.

    EXAMPLES::

        sage: S = SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2]); S
        Semistandard multi skew tableaux of shape [[2, 1] / [], [2, 2] / [1]] and weight [2, 2, 2]
        sage: S.list()
        [[[[1, 1], [2]], [[None, 2], [3, 3]]],
         [[[1, 2], [2]], [[None, 1], [3, 3]]],
         [[[1, 3], [2]], [[None, 2], [1, 3]]],
         [[[1, 3], [2]], [[None, 1], [2, 3]]],
         [[[1, 1], [3]], [[None, 2], [2, 3]]],
         [[[1, 2], [3]], [[None, 2], [1, 3]]],
         [[[1, 2], [3]], [[None, 1], [2, 3]]],
         [[[2, 2], [3]], [[None, 1], [1, 3]]],
         [[[1, 3], [3]], [[None, 1], [2, 2]]],
         [[[2, 3], [3]], [[None, 1], [1, 2]]]]
    """
    @staticmethod
    def __classcall_private__(cls, shape, weight):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S1 = SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2])
            sage: shape_alt = ( SkewPartition([[2,1],[]]), SkewPartition([[2,2],[1]]) )
            sage: S2 = SemistandardMultiSkewTableaux(shape_alt, (2,2,2))
            sage: S1 is S2
            True
        """
        shape = tuple(SkewPartition(x) for x in shape)
        weight = Partition(weight)

        if sum(weight) != sum(s.size() for s in shape):
            raise ValueError("the sum of weight must be the sum of the sizes of shape")

        return super().__classcall__(cls, shape, weight)

    def __init__(self, shape, weight):
        """
        TESTS::

            sage: S = SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2])
            sage: TestSuite(S).run()
        """
        self._shape = shape
        self._weight = weight
        MultiSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2])
            Semistandard multi skew tableaux of shape [[2, 1] / [], [2, 2] / [1]] and weight [2, 2, 2]
        """
        return "Semistandard multi skew tableaux of shape %s and weight %s" % (list(self._shape), self._weight)

    def __contains__(self, x):
        """
        TESTS::

            sage: s = SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2])
            sage: all(i in s for i in s)
            True
        """
        try:
            x = MultiSkewTableau(x)
        except TypeError:
            return False
        if x.weight() != list(self._weight):
            return False
        if x.shape() != list(self._shape):
            return False
        return all(xi.is_semistandard() for xi in x)

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: sp = SkewPartitions(3).list()
            sage: SemistandardMultiSkewTableaux([SkewPartition([[1, 1, 1], []]), SkewPartition([[3], []])],[2,2,2]).list()
            [[[[1], [2], [3]], [[1, 2, 3]]]]

            sage: a = SkewPartition([[8,7,6,5,1,1],[2,1,1]])
            sage: weight = [3,3,2]
            sage: k = 3
            sage: s = SemistandardMultiSkewTableaux(a.quotient(k),weight)
            sage: len(s.list())
            34
            sage: RibbonTableaux(a,weight,k).cardinality()
            34

        TESTS:

        Check that :issue:`36196` is fixed::

            sage: shapes = [[[1], [0]], [[1], [0]], [[1], [0]]]
            sage: weight = [1, 1, 1]
            sage: SMST = SemistandardMultiSkewTableaux(shapes, weight)
            sage: list(SMST)
            [[[[1]], [[2]], [[3]]],
             [[[2]], [[1]], [[3]]],
             [[[1]], [[3]], [[2]]],
             [[[2]], [[3]], [[1]]],
             [[[3]], [[1]], [[2]]],
             [[[3]], [[2]], [[1]]]]
        """
        parts = self._shape
        mu = self._weight

        # Splitting the partition
        s = [p.size() for p in parts]
        parts = [p.to_list() for p in parts]

        # Gluing the partitions
        parttmp = parts[0]
        for i in range(1, len(parts)):
            trans = parttmp[0][0]
            current_part = parts[i]
            current_part[1] += [0] * (len(current_part[0]) - len(current_part[1]))
            inner_current = [trans + j for j in current_part[1]]
            outer_current = [trans + j for j in current_part[0]]
            parttmp = [outer_current + parttmp[0], inner_current + parttmp[1]]

        # List the corresponding skew tableaux
        l = (st.to_word() for st in SemistandardSkewTableaux(parttmp, mu))

        S = SkewTableaux()
        for lk in l:
            pos = 0  # Double check this
            lk = list(lk)
            w = lk[:s[0]]
            restmp = [S.from_shape_and_word(parts[0], w)]
            for i in range(1, len(parts)):
                pos += s[i-1]
                w = lk[pos: pos + s[i]]
                restmp.append(S.from_shape_and_word(parts[i], w))
            yield self.element_class(self, restmp)


class RibbonTableau_class(RibbonTableau):
    """
    This exists solely for unpickling ``RibbonTableau_class`` objects.
    """

    def __setstate__(self, state):
        r"""
        Unpickle old ``RibbonTableau_class`` objects.

        TESTS::

            sage: loads(b'x\x9c5\xcc\xbd\x0e\xc2 \x14@\xe1\xb4Z\x7f\xd0\x07\xc1\x85D}\x8f\x0e\x8d\x1d\t\xb9\x90\x1bJ\xa44\x17\xe8h\xa2\x83\xef-\xda\xb8\x9do9\xcf\xda$\xb0(\xcc4j\x17 \x8b\xe8\xb4\x9e\x82\xca\xa0=\xc2\xcc\xba\x1fo\x8b\x94\xf1\x90\x12\xa3\xea\xf4\xa2\xfaA+\xde7j\x804\xd0\xba-\xe5]\xca\xd4H\xdapI[\xde.\xdf\xe8\x82M\xc2\x85\x8c\x16#\x1b\xe1\x8e\xea\x0f\xda\xf5\xd5\xf9\xdd\xd1\x1e%1>\x14]\x8a\x0e\xdf\xb8\x968"\xceZ|\x00x\xef5\x11')
            [[None, 1], [2, 3]]
            sage: loads(dumps( RibbonTableau([[None, 1],[2,3]]) ))
            [[None, 1], [2, 3]]
        """
        self.__class__ = RibbonTableau
        self.__init__(RibbonTableaux(), state['_list'])


register_unpickle_override('sage.combinat.ribbon_tableau', 'RibbonTableau_class', RibbonTableau_class)
register_unpickle_override('sage.combinat.ribbon_tableau', 'RibbonTableaux_shapeweightlength', RibbonTableaux)
register_unpickle_override('sage.combinat.ribbon_tableau', 'SemistandardMultiSkewTtableaux_shapeweight', SemistandardMultiSkewTableaux)
