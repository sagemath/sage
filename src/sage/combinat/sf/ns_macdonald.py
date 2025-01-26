# sage.doctest: needs sage.combinat sage.modules
"""
Non-symmetric Macdonald polynomials
"""
import copy

from sage.combinat.combinat import CombinatorialObject
from sage.combinat.words.word import Word
from sage.combinat.combination import Combinations
from sage.combinat.permutation import Permutation
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.misc.misc_c import prod
from sage.combinat.backtrack import GenericBacktracker
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation


class LatticeDiagram(CombinatorialObject):
    def boxes(self):
        """
        EXAMPLES::

            sage: a = LatticeDiagram([3,0,2])
            sage: a.boxes()
            [(1, 1), (1, 2), (1, 3), (3, 1), (3, 2)]
            sage: a = LatticeDiagram([2, 1, 3, 0, 0, 2])
            sage: a.boxes()
            [(1, 1), (1, 2), (2, 1), (3, 1), (3, 2), (3, 3), (6, 1), (6, 2)]
        """
        res = []
        for i in range(1, len(self) + 1):
            res += [(i, j + 1) for j in range(self[i])]
        return res

    def __getitem__(self, i):
        """
        Return the `i`-th entry of ``self``.

        Note that the indexing for lattice diagrams starts at `1`.

        EXAMPLES::

            sage: a = LatticeDiagram([3,0,2])
            sage: a[1]
            3
            sage: a[0]
            Traceback (most recent call last):
            ...
            ValueError: indexing starts at 1
            sage: a[-1]
            2
        """
        if i == 0:
            raise ValueError("indexing starts at 1")
        elif i < 0:
            i += 1
        return self._list[i - 1]

    def leg(self, i, j):
        """
        Return the leg of the box ``(i,j)`` in ``self``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.leg(5,2)
            [(5, 3)]
        """
        return [(i, k) for k in range(j + 1, self[i] + 1)]

    def arm_left(self, i, j):
        """
        Return the left arm of the box ``(i,j)`` in ``self``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.arm_left(5,2)
            [(1, 2), (3, 2)]
        """
        return [(ip, j) for ip in range(1, i) if j <= self[ip] <= self[i]]

    def arm_right(self, i, j):
        """
        Return the right arm of the box ``(i,j)`` in ``self``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.arm_right(5,2)
            [(8, 1)]
        """
        return [(ip, j - 1) for ip in range(i + 1, len(self) + 1)
                if j - 1 <= self[ip] < self[i]]

    def arm(self, i, j):
        """
        Return the arm of the box ``(i,j)`` in ``self``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.arm(5,2)
            [(1, 2), (3, 2), (8, 1)]
        """
        return self.arm_left(i, j) + self.arm_right(i, j)

    def l(self, i, j):
        """
        Return ``self[i] - j``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.l(5,2)
            1
        """
        return self[i] - j

    def a(self, i, j):
        """
        Return the length of the arm of the box ``(i,j)`` in ``self``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.a(5,2)
            3
        """
        return len(self.arm(i, j))

    def size(self):
        """
        Return the number of boxes in ``self``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.size()
            22
        """
        return sum(self._list)

    def flip(self):
        """
        Return the flip of ``self``, where flip is defined as follows. Let
        ``r = max(self)``. Then ``self.flip()[i] = r - self[i]``.

        EXAMPLES::

            sage: a = LatticeDiagram([3,0,2])
            sage: a.flip()
            [0, 3, 1]
        """
        r = max(self)
        return LatticeDiagram([r - i for i in self])

    def boxes_same_and_lower_right(self, ii, jj):
        """
        Return an iterator of the boxes of ``self`` that are in row ``jj``
        but not identical with ``(ii, jj)``, or lie in the row
        ``jj - 1`` (the row directly below ``jj``; this might be the
        basement) and strictly to the right of ``(ii, jj)``.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a = a.shape()
            sage: list(a.boxes_same_and_lower_right(1,1))
            [(2, 1), (3, 1), (6, 1), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0)]
            sage: list(a.boxes_same_and_lower_right(1,2))
            [(3, 2), (6, 2), (2, 1), (3, 1), (6, 1)]
            sage: list(a.boxes_same_and_lower_right(3,3))
            [(6, 2)]
            sage: list(a.boxes_same_and_lower_right(2,3))
            [(3, 3), (3, 2), (6, 2)]
        """
        # Add all of the boxes in the same row
        for i in range(1, len(self) + 1):
            if self[i] >= jj and i != ii:
                yield (i, jj)

        for i in range(ii + 1, len(self) + 1):
            if self[i] >= jj - 1:
                yield (i, jj - 1)


class AugmentedLatticeDiagramFilling(CombinatorialObject):
    def __init__(self, l, pi=None):
        """
        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a == loads(dumps(a))
            True
            sage: pi = Permutation([2,3,1]).to_permutation_group_element()
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]],pi)
            sage: a == loads(dumps(a))
            True
        """
        if pi is None:
            pi = Permutation([1]).to_permutation_group_element()
        self._list = [[pi(i + 1)] + li for i, li in enumerate(l)]

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a[0]
            Traceback (most recent call last):
            ...
            ValueError: indexing starts at 1
            sage: a[1,0]
            1
            sage: a[2,0]
            2
            sage: a[3,2]
            4
            sage: a[3][2]
            4
        """
        if isinstance(i, tuple):
            i, j = i
            return self._list[i - 1][j]
        if i < 1:
            raise ValueError("indexing starts at 1")
        return self._list[i - 1]

    def shape(self):
        """
        Return the shape of ``self``.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.shape()
            [2, 1, 3, 0, 0, 2]
        """
        return LatticeDiagram([max(0, len(self[i]) - 1)
                               for i in range(1, len(self) + 1)])

    def __contains__(self, ij):
        """
        Return ``True`` if the box ``(i,j) (= ij)`` is in ``self``. Note that this
        does not include the basement row.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: (1,1) in a
            True
            sage: (1,0) in a
            False
        """
        i, j = ij
        return 0 < i <= len(self) and 0 < j <= len(self[i])

    def are_attacking(self, i, j, ii, jj):
        """
        Return ``True`` if the boxes ``(i,j)`` and ``(ii,jj)`` in ``self`` are attacking.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: all( a.are_attacking(i,j,ii,jj) for (i,j),(ii,jj) in a.attacking_boxes())
            True
            sage: a.are_attacking(1,1,3,2)
            False
        """
        # If the two boxes are at the same height,
        # then they are attacking
        if j == jj:
            return True

        # Make it so that the lower box is in position i,j
        if jj < j:
            i, j, ii, jj = ii, jj, i, j

        # If the lower box is one row below and
        # strictly to the right, then they are
        # attacking.
        return j == jj - 1 and i > ii

    def boxes(self):
        """
        Return a list of the coordinates of the boxes of ``self``, including
        the 'basement row'.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.boxes()
            [(1, 1),
             (1, 2),
             (2, 1),
             (3, 1),
             (3, 2),
             (3, 3),
             (6, 1),
             (6, 2),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0)]
        """
        return self.shape().boxes() + [(i, 0)
                                       for i in range(1, len(self.shape()) + 1)]

    def attacking_boxes(self):
        """
        Return a list of pairs of boxes in ``self`` that are attacking.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.attacking_boxes()[:5]
            [((1, 1), (2, 1)),
             ((1, 1), (3, 1)),
             ((1, 1), (6, 1)),
             ((1, 1), (2, 0)),
             ((1, 1), (3, 0))]
        """
        boxes = self.boxes()
        res = []
        for (i, j), (ii, jj) in Combinations(boxes, 2):
            if self.are_attacking(i, j, ii, jj):
                res.append(((i, j), (ii, jj)))
        return res

    def is_non_attacking(self):
        """
        Return ``True`` if ``self`` is non-attacking.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.is_non_attacking()
            True
            sage: a = AugmentedLatticeDiagramFilling([[1, 1, 1], [2, 3], [3]])
            sage: a.is_non_attacking()
            False
            sage: a = AugmentedLatticeDiagramFilling([[2,2],[1]])
            sage: a.is_non_attacking()
            False
            sage: pi = Permutation([2,1]).to_permutation_group_element()
            sage: a = AugmentedLatticeDiagramFilling([[2,2],[1]],pi)
            sage: a.is_non_attacking()
            True
        """
        for a, b in self.attacking_boxes():
            if self[a] == self[b]:
                return False
        return True

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.weight()
            [1, 2, 1, 1, 2, 1]
        """
        ed = self.reading_word().evaluation_dict()
        entries = list(ed)
        m = max(entries) + 1 if entries else -1
        return [ed.get(k, 0) for k in range(1, m)]

    def descents(self):
        """
        Return a list of the descents of ``self``.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.descents()
            [(1, 2), (3, 2)]
        """
        res = []
        for i, j in self.shape().boxes():
            if self[i, j] > self[i, j - 1]:
                res.append((i, j))
        return res

    def maj(self):
        """
        Return the major index of ``self``.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.maj()
            3
        """
        res = 0
        shape = self.shape()
        for i, j in self.descents():
            res += shape.l(i, j) + 1
        return res

    def reading_order(self):
        """
        Return a list of coordinates of the boxes in ``self``, starting from
        the top right, and reading from right to left.

        Note that this includes the 'basement row' of ``self``.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.reading_order()
            [(3, 3),
             (6, 2),
             (3, 2),
             (1, 2),
             (6, 1),
             (3, 1),
             (2, 1),
             (1, 1),
             (6, 0),
             (5, 0),
             (4, 0),
             (3, 0),
             (2, 0),
             (1, 0)]
        """
        boxes = self.boxes()

        def fn(ij):
            return (-ij[1], -ij[0])
        boxes.sort(key=fn)
        return boxes

    def reading_word(self):
        """
        Return the reading word of ``self``, obtained by reading the boxes
        entries of ``self`` from right to left, starting in the upper right.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.reading_word()
            word: 25465321
        """
        w = [self[i, j] for i, j in self.reading_order() if j > 0]
        return Word(w)

    def inversions(self):
        """
        Return a list of the inversions of ``self``.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.inversions()[:5]
            [((6, 2), (3, 2)),
             ((1, 2), (6, 1)),
             ((1, 2), (3, 1)),
             ((1, 2), (2, 1)),
             ((6, 1), (3, 1))]
            sage: len(a.inversions())
            25
        """
        atboxes = [set(x) for x in self.attacking_boxes()]
        res = []
        order = self.reading_order()
        for a in range(len(order)):
            i, j = order[a]
            for b in range(a + 1, len(order)):
                ii, jj = order[b]
                if self[i, j] > self[ii, jj] and set(((i, j), (ii, jj))) in atboxes:
                    res.append(((i, j), (ii, jj)))
        return res

    def _inv_aux(self):
        """
        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a._inv_aux()
            7
        """
        res = 0
        shape = self.shape()
        for i in range(1, len(self) + 1):
            for j in range(i + 1, len(self) + 1):
                if shape[i] <= shape[j]:
                    res += 1
        return res

    def inv(self):
        """
        Return ``self``'s inversion statistic.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.inv()
            15
        """
        res = len(self.inversions())
        res -= sum(self.shape().a(i, j) for i, j in self.descents())
        res -= self._inv_aux()
        return res

    def coinv(self):
        """
        Return ``self``'s co-inversion statistic.

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.coinv()
            2
        """
        shape = self.shape()
        return sum(shape.a(i, j) for i, j in shape.boxes()) - self.inv()

    def coeff(self, q, t):
        """
        Return the coefficient in front of ``self`` in the HHL formula for the
        expansion of the non-symmetric Macdonald polynomial
        E(self.shape()).

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: q,t = var('q,t')                                                      # needs sage.symbolic
            sage: a.coeff(q,t)                                                          # needs sage.symbolic
            (t - 1)^4/((q^2*t^3 - 1)^2*(q*t^2 - 1)^2)
        """
        res = 1
        shape = self.shape()
        for i, j in shape.boxes():
            if self[i, j] != self[i, j - 1]:
                res *= (1 - t) / (1 - q**(shape.l(i, j) + 1)
                                  * t**(shape.a(i, j) + 1))
        return res

    def coeff_integral(self, q, t):
        """
        Return the coefficient in front of ``self`` in the HHL formula for the
        expansion of the integral non-symmetric Macdonald polynomial
        E(self.shape())

        EXAMPLES::

            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: q,t = var('q,t')                                                      # needs sage.symbolic
            sage: a.coeff_integral(q,t)                                                 # needs sage.symbolic
            (q^2*t^3 - 1)^2*(q*t^2 - 1)^2*(t - 1)^4
        """
        res = 1
        shape = self.shape()
        for i, j in shape.boxes():
            if self[i, j] != self[i, j - 1]:
                res *= (1 - q**(shape.l(i, j) + 1) * t**(shape.a(i, j) + 1))
        for i, j in shape.boxes():
            if self[i, j] == self[i, j - 1]:
                res *= (1 - t)
        return res

    def permuted_filling(self, sigma):
        """
        EXAMPLES::

            sage: pi=Permutation([2,1,4,3]).to_permutation_group_element()
            sage: fill=[[2],[1,2,3],[],[3,1]]
            sage: AugmentedLatticeDiagramFilling(fill).permuted_filling(pi)
            [[2, 1], [1, 2, 1, 4], [4], [3, 4, 2]]
        """
        new_filling = []
        for col in self:
            nc = [sigma(x) for x in col]
            nc.pop(0)
            new_filling.append(nc)
        return AugmentedLatticeDiagramFilling(new_filling, sigma)


def NonattackingFillings(shape, pi=None):
    """
    Returning the finite set of nonattacking fillings of a
    given shape.

    EXAMPLES::

        sage: NonattackingFillings([0,1,2])
        Nonattacking fillings of [0, 1, 2]
        sage: NonattackingFillings([0,1,2]).list()
        [[[1], [2, 1], [3, 2, 1]],
         [[1], [2, 1], [3, 2, 2]],
         [[1], [2, 1], [3, 2, 3]],
         [[1], [2, 1], [3, 3, 1]],
         [[1], [2, 1], [3, 3, 2]],
         [[1], [2, 1], [3, 3, 3]],
         [[1], [2, 2], [3, 1, 1]],
         [[1], [2, 2], [3, 1, 2]],
         [[1], [2, 2], [3, 1, 3]],
         [[1], [2, 2], [3, 3, 1]],
         [[1], [2, 2], [3, 3, 2]],
         [[1], [2, 2], [3, 3, 3]]]
    """
    return NonattackingFillings_shape(tuple(shape), pi)


class NonattackingFillings_shape(Parent, UniqueRepresentation):
    def __init__(self, shape, pi=None):
        """
        EXAMPLES::

            sage: n = NonattackingFillings([0,1,2])
            sage: n == loads(dumps(n))
            True
        """
        self.pi = pi
        self._shape = LatticeDiagram(shape)
        self._name = "Nonattacking fillings of %s" % list(shape)
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __repr__(self):
        """
        EXAMPLES::

            sage: NonattackingFillings([0,1,2])
            Nonattacking fillings of [0, 1, 2]
        """
        return self._name

    def flip(self):
        """
        Return the nonattacking fillings of the flipped shape.

        EXAMPLES::

            sage: NonattackingFillings([0,1,2]).flip()
            Nonattacking fillings of [2, 1, 0]
        """
        return NonattackingFillings(list(self._shape.flip()), self.pi)

    def __iter__(self):
        """
        EXAMPLES::

            sage: NonattackingFillings([0,1,2]).list()  #indirect doctest
            [[[1], [2, 1], [3, 2, 1]],
             [[1], [2, 1], [3, 2, 2]],
             [[1], [2, 1], [3, 2, 3]],
             [[1], [2, 1], [3, 3, 1]],
             [[1], [2, 1], [3, 3, 2]],
             [[1], [2, 1], [3, 3, 3]],
             [[1], [2, 2], [3, 1, 1]],
             [[1], [2, 2], [3, 1, 2]],
             [[1], [2, 2], [3, 1, 3]],
             [[1], [2, 2], [3, 3, 1]],
             [[1], [2, 2], [3, 3, 2]],
             [[1], [2, 2], [3, 3, 3]]]
            sage: len(_)
            12

        TESTS::

            sage: NonattackingFillings([3,2,1,1]).cardinality()
            3
            sage: NonattackingFillings([3,2,1,2]).cardinality()
            4
            sage: NonattackingFillings([1,2,3]).cardinality()
            12
            sage: NonattackingFillings([2,2,2]).cardinality()
            1
            sage: NonattackingFillings([1,2,3,2]).cardinality()
            24
        """
        if sum(self._shape) == 0:
            yield AugmentedLatticeDiagramFilling([[] for _ in self._shape],
                                                 self.pi)
            return

        for z in NonattackingBacktracker(self._shape, self.pi):
            yield AugmentedLatticeDiagramFilling(z, self.pi)


class NonattackingBacktracker(GenericBacktracker):
    def __init__(self, shape, pi=None):
        """
        EXAMPLES::

            sage: from sage.combinat.sf.ns_macdonald import NonattackingBacktracker
            sage: n = NonattackingBacktracker(LatticeDiagram([0,1,2]))
            sage: n._ending_position
            (3, 2)
            sage: n._initial_state
            (2, 1)
        """
        self._shape = shape
        self._n = sum(shape)
        self._initial_data = [[None] * s for s in shape]
        if pi is None:
            pi = Permutation([1]).to_permutation_group_element()
        self.pi = pi

        # The ending position will be at the highest box
        # which is farthest right
        ending_row = max(shape)
        ending_col = len(shape) - list(reversed(list(shape))).index(ending_row)
        self._ending_position = (ending_col, ending_row)

        # Get the lowest box that is farthest left
        starting_row = 1
        nonzero = [i for i in shape if i != 0]
        starting_col = list(shape).index(nonzero[0]) + 1

        GenericBacktracker.__init__(self, self._initial_data, (starting_col, starting_row))

    def _rec(self, obj, state):
        """
        EXAMPLES::

            sage: from sage.combinat.sf.ns_macdonald import NonattackingBacktracker
            sage: n = NonattackingBacktracker(LatticeDiagram([0,1,2]))
            sage: len(list(n))
            12
            sage: obj = [ [], [None], [None, None]]
            sage: state = 2, 1
            sage: list(n._rec(obj, state))
            [([[], [1], [None, None]], (3, 1), False),
             ([[], [2], [None, None]], (3, 1), False)]
        """
        # We need to set the i,j^th entry.
        i, j = state

        # Get the next state
        new_state = self.get_next_pos(i, j)
        yld = bool(new_state is None)

        for k in range(1, len(self._shape) + 1):
            # We check to make sure that k does not
            # violate any of the attacking conditions
            if j == 1 and any(self.pi(x + 1) == k
                              for x in range(i, len(self._shape))):
                continue
            if any(obj[ii - 1][jj - 1] == k for ii, jj in
                    self._shape.boxes_same_and_lower_right(i, j) if jj != 0):
                continue

            # Fill in the in the i,j box with k+1
            obj[i - 1][j - 1] = k

            # Yield the object
            yield copy.deepcopy(obj), new_state, yld

    def get_next_pos(self, ii, jj):
        """
        EXAMPLES::

            sage: from sage.combinat.sf.ns_macdonald import NonattackingBacktracker
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: n = NonattackingBacktracker(a.shape())
            sage: n.get_next_pos(1, 1)
            (2, 1)
            sage: n.get_next_pos(6, 1)
            (1, 2)
            sage: n = NonattackingBacktracker(LatticeDiagram([2,2,2]))
            sage: n.get_next_pos(3, 1)
            (1, 2)
        """
        if (ii, jj) == self._ending_position:
            return None

        for i in range(ii + 1, len(self._shape) + 1):
            if self._shape[i] >= jj:
                return i, jj

        for i in range(1, ii + 1):
            if self._shape[i] >= jj + 1:
                return i, jj + 1

        raise ValueError("we should never be here")


def _check_muqt(mu, q, t, pi=None):
    """
    EXAMPLES::

        sage: from sage.combinat.sf.ns_macdonald import _check_muqt
        sage: P, q, t, n, R, x = _check_muqt([0,0,1],None,None)
        sage: P
        Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: q
        q
        sage: t
        t
        sage: n
        Nonattacking fillings of [0, 0, 1]
        sage: R
        Multivariate Polynomial Ring in x0, x1, x2 over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: x
        (x0, x1, x2)

    ::

        sage: q,t = var('q,t')                                                          # needs sage.symbolic
        sage: P, q, t, n, R, x = _check_muqt([0,0,1],q,None)                            # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: you must specify either both q and t or neither of them

    ::

        sage: P, q, t, n, R, x = _check_muqt([0,0,1],q,2)                               # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: the parents of q and t must be the same
    """
    if q is None and t is None:
        P = PolynomialRing(QQ, 'q,t').fraction_field()
        q, t = P.gens()
    elif q is not None and t is not None:
        if q.parent() != t.parent():
            raise ValueError("the parents of q and t must be the same")
        P = q.parent()
    else:
        raise ValueError("you must specify either both q and t or neither of them")
    n = NonattackingFillings(mu, pi)
    R = PolynomialRing(P, len(n._shape), 'x')
    x = R.gens()
    return P, q, t, n, R, x


def E(mu, q=None, t=None, pi=None):
    """
    Return the non-symmetric Macdonald polynomial in type A
    corresponding to a shape ``mu``, with basement permuted according to
    ``pi``.

    Note that if both `q` and `t` are specified, then they must have
    the same parent.

    REFERENCE:

    - J. Haglund, M. Haiman, N. Loehr.
      *A combinatorial formula for non-symmetric Macdonald polynomials*.
      :arxiv:`math/0601693v3`.

    .. SEEALSO::

        :class:`~sage.combinat.root_system.non_symmetric_macdonald_polynomials.NonSymmetricMacdonaldPolynomials`
        for a type free implementation where the polynomials are
        constructed recursively by the application of intertwining
        operators.

    EXAMPLES::

        sage: from sage.combinat.sf.ns_macdonald import E
        sage: E([0,0,0])
        1
        sage: E([1,0,0])
        x0
        sage: E([0,1,0])
        (t - 1)/(q*t^2 - 1)*x0 + x1
        sage: E([0,0,1])
        (t - 1)/(q*t - 1)*x0 + (t - 1)/(q*t - 1)*x1 + x2
        sage: E([1,1,0])
        x0*x1
        sage: E([1,0,1])
        (t - 1)/(q*t^2 - 1)*x0*x1 + x0*x2
        sage: E([0,1,1])
        (t - 1)/(q*t - 1)*x0*x1 + (t - 1)/(q*t - 1)*x0*x2 + x1*x2
        sage: E([2,0,0])
        x0^2 + (q*t - q)/(q*t - 1)*x0*x1 + (q*t - q)/(q*t - 1)*x0*x2
        sage: E([0,2,0])
        (t - 1)/(q^2*t^2 - 1)*x0^2 + (q^2*t^3 - q^2*t^2 + q*t^2 - 2*q*t + q - t + 1)/(q^3*t^3 - q^2*t^2 - q*t + 1)*x0*x1 + x1^2 + (q*t^2 - 2*q*t + q)/(q^3*t^3 - q^2*t^2 - q*t + 1)*x0*x2 + (q*t - q)/(q*t - 1)*x1*x2
    """
    P, q, t, n, R, x = _check_muqt(mu, q, t, pi)
    res = 0
    for a in n:
        weight = a.weight()
        res += q**a.maj() * t**a.coinv() * a.coeff(q, t) * prod(x[i]**weight[i] for i in range(len(weight)))
    return res


def E_integral(mu, q=None, t=None, pi=None):
    """
    Return the integral form for the non-symmetric Macdonald
    polynomial in type A corresponding to a shape mu.

    Note that if both q and t are specified, then they must have the
    same parent.

    REFERENCE:

    - J. Haglund, M. Haiman, N. Loehr.
      *A combinatorial formula for non-symmetric Macdonald polynomials*.
      :arxiv:`math/0601693v3`.

    EXAMPLES::

        sage: from sage.combinat.sf.ns_macdonald import E_integral
        sage: E_integral([0,0,0])
        1
        sage: E_integral([1,0,0])
        (-t + 1)*x0
        sage: E_integral([0,1,0])
        (-q*t^2 + 1)*x0 + (-t + 1)*x1
        sage: E_integral([0,0,1])
        (-q*t + 1)*x0 + (-q*t + 1)*x1 + (-t + 1)*x2
        sage: E_integral([1,1,0])
        (t^2 - 2*t + 1)*x0*x1
        sage: E_integral([1,0,1])
        (q*t^3 - q*t^2 - t + 1)*x0*x1 + (t^2 - 2*t + 1)*x0*x2
        sage: E_integral([0,1,1])
        (q^2*t^3 + q*t^4 - q*t^3 - q*t^2 - q*t - t^2 + t + 1)*x0*x1 + (q*t^2 - q*t - t + 1)*x0*x2 + (t^2 - 2*t + 1)*x1*x2
        sage: E_integral([2,0,0])
        (t^2 - 2*t + 1)*x0^2 + (q^2*t^2 - q^2*t - q*t + q)*x0*x1 + (q^2*t^2 - q^2*t - q*t + q)*x0*x2
        sage: E_integral([0,2,0])
        (q^2*t^3 - q^2*t^2 - t + 1)*x0^2 + (q^4*t^3 - q^3*t^2 - q^2*t + q*t^2 - q*t + q - t + 1)*x0*x1 + (t^2 - 2*t + 1)*x1^2 + (q^4*t^3 - q^3*t^2 - q^2*t + q)*x0*x2 + (q^2*t^2 - q^2*t - q*t + q)*x1*x2
    """
    P, q, t, n, R, x = _check_muqt(mu, q, t, pi)
    res = 0
    for a in n:
        weight = a.weight()
        res += q**a.maj() * t**a.coinv() * a.coeff_integral(q, t) * prod(x[i]**weight[i] for i in range(len(weight)))
    return res


def Ht(mu, q=None, t=None, pi=None):
    """
    Return the symmetric Macdonald polynomial using the Haiman,
    Haglund, and Loehr formula.

    Note that if both `q` and `t` are specified, then they must have the
    same parent.

    REFERENCE:

    - J. Haglund, M. Haiman, N. Loehr.
      *A combinatorial formula for non-symmetric Macdonald polynomials*.
      :arxiv:`math/0601693v3`.

    EXAMPLES::

        sage: from sage.combinat.sf.ns_macdonald import Ht
        sage: HHt = SymmetricFunctions(QQ['q','t'].fraction_field()).macdonald().Ht()
        sage: Ht([0,0,1])
        x0 + x1 + x2
        sage: HHt([1]).expand(3)
        x0 + x1 + x2
        sage: Ht([0,0,2])
        x0^2 + (q + 1)*x0*x1 + x1^2 + (q + 1)*x0*x2 + (q + 1)*x1*x2 + x2^2
        sage: HHt([2]).expand(3)
        x0^2 + (q + 1)*x0*x1 + x1^2 + (q + 1)*x0*x2 + (q + 1)*x1*x2 + x2^2
    """
    P, q, t, n, R, x = _check_muqt(mu, q, t, pi)
    res = 0
    for a in n:
        weight = a.weight()
        res += q**a.maj() * t**a.inv() * prod(x[i]**weight[i] for i in range(len(weight)))
    return res
