r"""
Kleshchev partitions
====================

A partition (tuple) `\mu` is Kleshchev if it can be recursively
obtained by adding a sequence of good nodes to the empty
:class:`PartitionTuple` of the same :meth:`~PartitionTuple.level`
and *multicharge*. In this way, the set of Kleshchev multipartitions becomes
a realization of a Kashiwara crystal :mod:`sage.combinat.crystals.crystals`
for a irreducible integral highest weight representation of
`U_q(\widehat{\mathfrak{sl}}_e)`.

The Kleshchev multipartitions first appeared in the work of Ariki and Mathas
[AM2000]_ where it was shown that they index the irreducible representations
of the cyclotomic Hecke algebras of type `A` [AK1994]_. Soon afterwards Ariki
[Ariki2001]_ showed that the set of Kleshchev multipartitions naturally label
the irreducible representations of these algebras.  As a far reaching
generalization of these ideas the Ariki-Brundan-Kleshchev categorification
theorem [Ariki1996]_ [BK2009]_ says that these algebras categorify the
irreducible integral highest weight representations of the quantum group
`U_q(\widehat{\mathfrak{sl}}_e)` of the affine special linear group. Under
this categorification, `q` corresponds to the grading shift on the
cyclotomic Hecke algebras, where the grading from the Brundan-Kleshchev
graded isomorphism theorem to the *KLR algebras* of type `A` [BK2009]_.

The group algebras of the symmetric group in characteristic `p` are an
important special case of the cyclotomic Hecke algebras of type `A`.
In this case, depending on your prefer convention, the set of Kleshchev
partitions is the set of *`p`-regular* or *`p`-restricted*
:class:`~sage.combinat.partition.Partitions`. In this case, Kleshchev
[Kle1995]_ proved that the *modular branching rules* were given by adding
and removing *good nodes*; see :meth:`~KleshchevPartition.good_cells`.
Lascoux, Leclerc and Thibon [LLT1996]_ noticed that Kleshchev's branching
rules coincided with Kashiwara's crystal operators for the fundamental
representation of `L(\Lambda_0)` of `U_q(\widehat{\mathfrak{sl}}_p)`
and their celebrated *LLT conjecture* said that decomposition matrices of
the :class:`sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra` of
the symmetric group should be computable using the canonical basis of
`L(\Lambda_0)`. This was proved and generalised to all cyclotomic Hecke
algebras of type `A` by Ariki [Ariki1996]_ and then further generalized
to the graded setting by Brundan and Kleshchev [BK2009]_.

The main class for accessing Kleshchev partition (tuples) is
:class:`~KleshchevPartitions`. Unfortunately, just as with the
symmetric group, different authors use different conventions when
defining Kleshchev partitions, which depends on whether you read
components from left to right, or right to left, and whether you
read the nodes in the partition in each component from top to bottom
or bottom to top. The :class:`~KleshchevPartitions` class supports
these four different conventions::

    sage: KleshchevPartitions(2, [0,0], size=2, convention='left regular')[:]
    [([1], [1]), ([2], [])]
    sage: KleshchevPartitions(2, [0,0], size=2, convention='left restricted')[:]
    [([1], [1]), ([], [1, 1])]
    sage: KleshchevPartitions(2, [0,0], size=2, convention='right regular')[:]
    [([1], [1]), ([], [2])]
    sage: KleshchevPartitions(2, [0,0], size=2, convention='right restricted')[:]
    [([1], [1]), ([1, 1], [])]

By default, the ``left restricted`` convention is used. As a shorthand,
``LG``, ``LS``, ``RG`` and ``RS``, respectively, can be used to specify the
``convention`` With the ``left`` convention the partition tuples should be
ordered with the most dominant partitions in the partition tuple on the left
and with the ``right`` convention the most dominant partition is on the right.

The :class:`~KleshchevPartitions` class can automatically convert between
these four different conventions::

    sage: KPlg = KleshchevPartitions(2, [0,0], size=2, convention='left regular')
    sage: KPls = KleshchevPartitions(2, [0,0], size=2, convention='left restricted')
    sage: [KPlg(mu) for mu in KPls] # indirect doc test
    [([1], [1]), ([2], [])]

AUTHORS:

- Andrew Mathas and Travis Scrimshaw (2018-05-1): Initial version
"""

from .partition import Partition, Partitions
from .partition_tuple import PartitionTuple, PartitionTuples

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.all import NN, ZZ, IntegerModRing
from sage.cpython.getattr import getattr_from_other_class

from collections import defaultdict

#--------------------------------------------------
# Kleshchev partition - element classes
#--------------------------------------------------
class KleshchevPartition(Partition):
    r"""
    Abstract base class for Kleshchev partitions. See
    :class:`~KleshchevPartitions`.
    """

    def conormal_cells(self, i=None):
        r"""
        Return a dictionary of the cells of ``self`` which are conormal.

        Following [Kle1995]_, the *conormal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and addable
        cells (in that order) from this list. The addable `i`-cells that
        remain at the end of the this process are the conormal `i`-cells.

        When computing conormal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of conormal cells
        is returned, which gives the conormal cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention="regular")
            sage: KP([5,4,4,3,2]).conormal_cells()
            {0: [(1, 4)], 1: [(5, 0), (4, 2)]}
            sage: KP([5,4,4,3,2]).conormal_cells(0)
            [(1, 4)]
            sage: KP([5,4,4,3,2]).conormal_cells(1)
            [(5, 0), (4, 2)]
            sage: KP = KleshchevPartitions(3, convention="restricted")
            sage: KP([5,4,4,3,2]).conormal_cells()
            {0: [(1, 4), (3, 3)], 2: [(0, 5)]}
        """
        # We use a dictionary for the conormal nodes as the indexing set is Z when e=0
        conormals = defaultdict(list)   # the conormal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells) - #(addable cells)

        # determine if we read up or down the partition
        KP = self.parent()
        rows = list(range(len(self)+1))
        if KP._convention[1] == 'G':
            rows.reverse()

        # work through the rows
        for row in rows:
            if row == len(self): # addable cell at bottom of partition
                res = KP._multicharge[0] - row
                if carry[res] == 0:
                    conormals[res].append((row, 0))
                else:
                    carry[res] += 1
            else:
                res = KP._multicharge[0] + self[row] - row - 1
                if row == len(self)-1 or self[row] > self[row+1]: # removable cell
                    carry[res] -= 1
                if row == 0 or self[row-1] > self[row]:               #addable cell
                    if carry[res+1] >= 0:
                        conormals[res+1].append((row, self[row]))
                    else:
                        carry[res+1] += 1

        # finally return the result
        return dict(conormals) if i is None else conormals[i]

    def cogood_cells(self, i=None):
        r"""
        Return a list of the cells of ``self`` that are cogood.

        The cogood `i`-cell is the 'last' conormal `i`-cell. As with the
        conormal cells we can choose to read either up or down the partition as
        specified by :meth:`~KleshchevPartitions.convention`.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of cogood cells
        is returned, which gives the cogood cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention="regular")
            sage: KP([5,4,4,3,2]).cogood_cells()
            {0: (1, 4), 1: (4, 2)}
            sage: KP([5,4,4,3,2]).cogood_cells(0)
            (1, 4)
            sage: KP([5,4,4,3,2]).cogood_cells(1)
            (4, 2)
            sage: KP = KleshchevPartitions(4, convention='restricted')
            sage: KP([5,4,4,3,2]).cogood_cells()
            {1: (0, 5), 2: (4, 2), 3: (1, 4)}
            sage: KP([5,4,4,3,2]).cogood_cells(0)
            sage: KP([5,4,4,3,2]).cogood_cells(2)
            (4, 2)
        """
        conormal_cells = self.conormal_cells(i)
        if i is None:
            return {i: conormal_cells[i][-1] for i in conormal_cells}
        elif not conormal_cells:
            return None

        return conormal_cells[-1]

    def normal_cells(self, i=None):
        r"""
        Return a dictionary of the cells of the partition that are normal.

        Following [Kle1995]_, the *normal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and
        addable cells (in that order) from this list. The removable
        `i`-cells that remain at the end of the this process are the
        normal `i`-cells.

        When computing normal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of normal cells
        is returned, which gives the normal cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([5,4,4,3,2]).normal_cells()
            {1: [(2, 3), (0, 4)]}
            sage: KP([5,4,4,3,2]).normal_cells(1)
            [(2, 3), (0, 4)]
            sage: KP = KleshchevPartitions(3, convention='restricted')
            sage: KP([5,4,4,3,2]).normal_cells()
            {0: [(4, 1)], 2: [(3, 2)]}
            sage: KP([5,4,4,3,2]).normal_cells(2)
            [(3, 2)]
        """
        # We use a dictionary for the normal nodes as the indexing set is Z when e=0
        normals = defaultdict(list)     # the normal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)

        # determine if we read up or down the partition
        KP = self.parent()
        rows = list(range(len(self)+1))
        if KP._convention[1] == 'S':
            rows.reverse()

        # work through the rows
        for row in rows:
            if row == len(self): # addable cell at bottom of partition
                carry[KP._multicharge[0]-row] += 1
            else:
                res = KP._multicharge[0] + self[row] - row - 1
                if row == len(self) - 1 or self[row] > self[row+1]: # removable cell
                    if carry[res] == 0:
                        normals[res].insert(0, (row, self[row]-1))
                    else:
                        carry[res] -= 1
                if row == 0 or self[row-1] > self[row]:              # addable cell
                    carry[res+1] += 1

        # finally return the result
        return dict(normals) if i is None else normals[i]

    def good_cells(self, i=None):
        """
        Return a list of the cells of ``self`` that are good.

        The good `i`-cell is the 'first' normal `i`-cell. As with the normal
        cells we can choose to read either up or down the partition as
        specified by :meth:`~KleshchevPartitions.convention`.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of good cells
        is returned, which gives the good cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP3 = KleshchevPartitions(3, convention='regular')
            sage: KP3([5,4,4,3,2]).good_cells()
            {1: (2, 3)}
            sage: KP3([5,4,4,3,2]).good_cells(1)
            (2, 3)
            sage: KP4 = KleshchevPartitions(4, convention='restricted')
            sage: KP4([5,4,4,3,2]).good_cells()
            {1: (2, 3)}
            sage: KP4([5,4,4,3,2]).good_cells(0)
            sage: KP4([5,4,4,3,2]).good_cells(1)
            (2, 3)
        """
        normal_cells = self.normal_cells(i)
        if i is None:
            return {j: normal_cells[j][0] for j in normal_cells}
        elif not normal_cells:
            return None

        return normal_cells[0]

    def good_residue_sequence(self):
        """
        Return a sequence of good nodes from the empty partition
        to ``self``, or ``None`` if no such sequence exists.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([5,4,4,3,2]).good_residue_sequence()
            [0, 2, 1, 1, 0, 2, 0, 2, 1, 1, 0, 2, 0, 2, 2, 0, 1, 1]
            sage: KP = KleshchevPartitions(3, convention='restricted')
            sage: KP([5,4,4,3,2]).good_residue_sequence()
            [0, 1, 2, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 0]
        """
        if not self:
            return []

        good_cells = self.good_cells()
        assert good_cells

        res = sorted(good_cells)[0]
        r, c = good_cells[res]
        good_seq = type(self)(self.parent(), self.remove_cell(r,c)).good_residue_sequence()
        good_seq.append(self.parent()._index_set(res))
        return good_seq

    def good_cell_sequence(self):
        """
        Return a sequence of good nodes from the empty partition
        to ``self``, or ``None`` if no such sequence exists.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([5,4,4,3,2]).good_cell_sequence()
            [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
             (3, 0), (2, 1), (1, 2), (3, 1), (0, 3), (1, 3),
             (2, 2), (3, 2), (4, 0), (4, 1), (0, 4), (2, 3)]
            sage: KP = KleshchevPartitions(3, convention='restricted')
            sage: KP([5,4,4,3,2]).good_cell_sequence()
            [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0),
             (0, 3), (2, 1), (1, 2), (1, 3), (3, 0), (3, 1),
             (2, 2), (4, 0), (2, 3), (3, 2), (0, 4), (4, 1)]
        """
        if not self:
            return []

        good_cells = self.good_cells()
        assert good_cells

        cell = good_cells[sorted(good_cells)[0]]
        good_seq = type(self)(self.parent(), self.remove_cell(*cell)).good_cell_sequence()
        good_seq.append(cell)
        return good_seq

    def mullineux_conjugate(self):
        r"""
        Return the partition tuple that is the Mullineux conjugate of ``self``.

        It follows from results in [BK2009]_, [Mat2015]_ that if `\nu` is the
        Mullineux conjugate of the Kleshchev partition tuple `\mu` then the
        simple module `D^\nu =(D^\mu)^{\text{sgn}}` is obtained from `D^\mu`
        by twisting by the `\text{sgn}`-automorphism with is the
        Iwahori-Hecke algebra analogue of tensoring with the one
        dimensional sign representation.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([5,4,4,3,2]).mullineux_conjugate()
            [9, 7, 1, 1]
            sage: KP = KleshchevPartitions(3, convention='restricted')
            sage: KP([5,4,4,3,2]).mullineux_conjugate()
            [3, 2, 2, 2, 2, 2, 2, 1, 1, 1]
            sage: KP = KleshchevPartitions(3, [2], convention='regular')
            sage: mc = KP([5,4,4,3,2]).mullineux_conjugate(); mc
            [9, 7, 1, 1]
            sage: mc.parent().multicharge()
            (1,)
            sage: KP = KleshchevPartitions(3, [2], convention='restricted')
            sage: mc = KP([5,4,4,3,2]).mullineux_conjugate(); mc
            [3, 2, 2, 2, 2, 2, 2, 1, 1, 1]
            sage: mc.parent().multicharge()
            (1,)
        """
        P = self.parent()
        if not self:
            size = None
            if isinstance(P, KleshchevPartitions_size):
                size = P._size
            KP = KleshchevPartitions(P._e, [-c for c in P._multicharge],
                                     size=size, convention=P._convention)
            return KP.element_class(KP, [])

        good_cells = self.good_cells()
        assert good_cells

        r, c = sorted(good_cells.values())[0]
        # This is technically wrong when the parent has a fixed size because
        #   the resulting Kleshchev partition after removing a cell has abs
        #   smaller size. However, this is useful to avoid constructing
        #   transient parents.
        mu = P.element_class(P, self.remove_cell(r, c)).mullineux_conjugate()
        # add back on a cogood cell of residue -residue(k,r,c)
        KP = mu.parent()
        return KP.element_class(KP, mu.add_cell(*mu.cogood_cells( r-c-self.parent()._multicharge[0]) ))

    def is_regular(self):
        r"""
        Return ``True`` if ``self`` is a `e`-regular partition tuple.

        A partition tuple is `e`-regular if we can get to the empty partition
        tuple by successively removing a sequence of good cells in the down
        direction. Equivalently, all partitions are `0`-regular and if `e > 0`
        then a partition is `e`-regular if no `e` non-zero parts of ``self``
        are equal.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2)
            sage: KP([2,1,1]).is_regular()
            False
            sage: KP = KleshchevPartitions(3)
            sage: KP([2,1,1]).is_regular()
            True
            sage: KP([]).is_regular()
            True
        """
        if self.size() == 0 or self.parent()._e == 0:
            return True
        KP = self.parent()
        return super(KleshchevPartition, self).is_regular(KP._e, KP._multicharge)

    def is_restricted(self):
        r"""
        Return ``True`` if ``self`` is an `e`-restricted partition tuple.

        A partition tuple is `e`-restricted if we can get to the empty
        partition tuple by successively removing a sequence of good cells in
        the up direction. Equivalently, all partitions are `0`-restricted and
        if `e > 0` then a partition is `e`-restricted if the difference of
        successive parts of ``self`` are always strictly less than `e`.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2, convention='regular')
            sage: KP([3,1]).is_restricted()
            False
            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([3,1]).is_restricted()
            True
            sage: KP([]).is_restricted()
            True
        """
        if self.size() == 0 or self.parent()._e == 0:
            return True
        KP = self.parent()
        return super(KleshchevPartition, self).is_restricted(KP._e, KP._multicharge)

class KleshchevPartitionTuple(PartitionTuple):
    r"""
    Abstract base class for Kleshchev partition tuples. See
    :class:`~KleshchevPartitions`.
    """

    def conormal_cells(self, i=None):
        r"""
        Return a dictionary of the cells of the partition that are conormal.

        Following [Kle1995]_, the *conormal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and addable
        cells (in that order) from this list. The addable `i`-cells that
        remain at the end of the this process are the conormal `i`-cells.

        When computing conormal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of conormal cells
        is returned, which gives the conormal cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, [0,1], convention="left regular")
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells()
            {0: [(1, 2, 1), (1, 1, 3), (1, 0, 5)],
            1: [(1, 3, 0), (0, 2, 0), (0, 1, 2), (0, 0, 4)]}
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells(1)
            [(1, 3, 0), (0, 2, 0), (0, 1, 2), (0, 0, 4)]
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells(2)
            []
            sage: KP = KleshchevPartitions(3, [0,1], convention="right restricted")
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells(0)
            [(1, 0, 5), (1, 1, 3), (1, 2, 1)]
        """
        # We use a dictionary for the conormal nodes as the indexing set is Z when e=0
        conormals = defaultdict(list)   # the conormal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)

        part_lens = [len(part) for part in self]  # so we don't repeatedly call these
        # the indices for the rows ending in addable nodes
        KP = self.parent()
        if KP._convention[0] == 'L':
            rows = [(k,r) for k,ell in enumerate(part_lens) for r in range(ell+1)]
        else:
            rows = [(k,r) for k,ell in reversed(list(enumerate(part_lens))) for r in range(ell+1)]
        if KP._convention[1] == 'G':
            rows.reverse()

        for row in rows:
            k,r = row
            if r == part_lens[k]: # addable cell at bottom of a component
                res = KP._multicharge[k] - r
                if carry[res] == 0:
                    conormals[res].append((k, r, 0))
                else:
                    carry[res] += 1
            else:
                part = self[k]
                res = KP._multicharge[k] + (part[r] - r - 1)
                if r == part_lens[k] - 1 or part[r] > part[r+1]: # removable cell
                    carry[res] -= 1
                if r == 0 or part[r-1] > part[r]:                # addable cell
                    if carry[res+1] == 0:
                        conormals[res+1].append((k, r, part[r]))
                    else:
                        carry[res+1] += 1

        # finally return the result
        if i is None:
            return dict(conormals)
        return conormals[i]

    def cogood_cells(self, i=None):
        r"""
        Return a list of the cells of the partition that are cogood.

        The cogood `i`-cell is the 'last' conormal `i`-cell. As with the
        conormal cells we can choose to read either up or down the partition
        as specified by :meth:`~KleshchevPartitions.convention`.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of cogood cells
        is returned, which gives the cogood cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, [0,1])
            sage: pt = KP([[4, 2], [5, 3, 1]])
            sage: pt.cogood_cells()
            {0: (1, 2, 1), 1: (1, 3, 0)}
            sage: pt.cogood_cells(0)
            (1, 2, 1)
            sage: KP = KleshchevPartitions(4, [0,1], convention="left regular")
            sage: pt = KP([[5, 2, 2], [6, 1, 1]])
            sage: pt.cogood_cells()
            {1: (0, 0, 5), 2: (1, 3, 0)}
            sage: pt.cogood_cells(0) is None
            True
            sage: pt.cogood_cells(1) is None
            False
        """
        conormal_cells = self.conormal_cells(i)
        if i is None:
            return {j: conormal_cells[j][-1] for j in conormal_cells}
        elif not conormal_cells:
            return None

        return conormal_cells[-1]

    def normal_cells(self, i=None):
        r"""
        Return a dictionary of the removable cells of the partition that
        are normal.

        Following [Kle1995]_, the *normal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and
        addable cells (in that order) from this list. The removable
        `i`-cells that remain at the end of the this process are the
        normal `i`-cells.

        When computing normal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of normal cells
        is returned, which gives the normal cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, [0,1], convention="left restricted")
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells()
            {2: [(1, 0, 4), (1, 1, 2), (1, 2, 0)]}
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells(1)
            []
            sage: KP = KleshchevPartitions(3, [0,1], convention="left regular")
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells()
            {0: [(0, 1, 1), (0, 0, 3)], 2: [(1, 2, 0), (1, 1, 2), (1, 0, 4)]}
            sage: KP = KleshchevPartitions(3, [0,1], convention="right regular")
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells()
            {2: [(1, 2, 0), (1, 1, 2), (1, 0, 4)]}
            sage: KP = KleshchevPartitions(3, [0,1], convention="right restricted")
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells()
            {0: [(0, 0, 3), (0, 1, 1)], 2: [(1, 0, 4), (1, 1, 2), (1, 2, 0)]}
        """
        # We use a dictionary for the normal nodes as the indexing set is Z when e=0
        normals = defaultdict(list)     # the normal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)

        part_lens = [len(part) for part in self]  # so we don't repeatedly call these
        KP = self.parent()
        if KP._convention[0] == 'L':
            rows = [(k,r) for k,ell in enumerate(part_lens) for r in range(ell+1)]
        else:
            rows = [(k,r) for k,ell in reversed(list(enumerate(part_lens))) for r in range(ell+1)]
        if KP._convention[1] == 'S':
            rows.reverse()

        for row in rows:
            k,r = row
            if r == part_lens[k]: # addable cell at bottom of a component
                carry[KP._multicharge[k]-r] += 1
            else:
                part = self[k]
                res = KP._multicharge[k] + (part[r] - r - 1)
                if r == part_lens[k]-1 or part[r] > part[r+1]: # removable cell
                    if carry[res] == 0:
                        normals[res].insert(0, (k, r, part[r]-1))
                    else:
                        carry[res] -= 1
                if r == 0 or part[r-1] > part[r]:               #addable cell
                    carry[res+1] += 1

        # finally return the result
        if i is None:
            return dict(normals)    # change the defaultdict into a dict

        return normals[i]

    def good_cells(self, i=None):
        r"""
        Return a list of the cells of the partition tuple which are good.

        The good `i`-cell is the 'first' normal `i`-cell. As with the normal
        cells we can choose to read either up or down the partition as specified
        by :meth:`~KleshchevPartitions.convention`.

        INPUT:

        - ``i`` -- (optional) a residue

        OUTPUT:

        If no residue ``i`` is specified then a dictionary of good cells
        is returned, which gives the good cells for ``0 <= i < e``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, [0,1])
            sage: pt = KP([[4, 2], [5, 3, 1]])
            sage: pt.good_cells()
            {2: (1, 0, 4)}
            sage: pt.good_cells(2)
            (1, 0, 4)
            sage: KP = KleshchevPartitions(4, [0,1], convention="left regular")
            sage: pt = KP([[5, 2, 2], [6, 2, 1]])
            sage: pt.good_cells()
            {0: (0, 0, 4), 2: (1, 0, 5), 3: (0, 2, 1)}
            sage: pt.good_cells(1) is None
            True
        """
        normal_cells = self.normal_cells(i)
        if i is None:
            return {j: normal_cells[j][0] for j in normal_cells}
        elif not normal_cells:
            return None

        return normal_cells[0]

    def good_residue_sequence(self):
        """
        Return a sequence of good nodes from the empty partition to ``self``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, [0,1])
            sage: KP([[4, 2], [5, 3, 1]]).good_residue_sequence()
            [0, 1, 2, 1, 2, 0, 1, 0, 2, 2, 0, 1, 0, 2, 2]
        """
        if self.size() == 0:
            return []
        good_cells = self.good_cells()
        assert good_cells

        res = sorted(good_cells.keys())[0]
        k, r, c = good_cells[res]
        good_seq = type(self)(self.parent(), self.remove_cell(k,r,c)).good_residue_sequence()
        good_seq.append( self.parent()._index_set(res) )
        return good_seq

    def good_cell_sequence(self):
        """
        Return a sequence of good nodes from the empty partition to ``self``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3,[0,1])
            sage: KP([[4, 2], [5, 3, 1]]).good_cell_sequence()
            [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1), (0, 1, 0),
             (1, 1, 0), (1, 1, 1), (1, 0, 2), (1, 2, 0), (0, 0, 2),
             (0, 1, 1), (1, 0, 3), (0, 0, 3), (1, 1, 2), (1, 0, 4)]
        """
        if self.size() == 0:
            return []
        good_cells = self.good_cells()
        assert good_cells

        cell = good_cells[sorted(good_cells)[0]]
        good_seq = type(self)(self.parent(), self.remove_cell(*cell)).good_cell_sequence()
        good_seq.append(cell)
        return good_seq

    def mullineux_conjugate(self):
        r"""
        Return the partition that is the Mullineux conjugate of ``self``.

        It follows from results in [Kle1996]_ [Bru1998]_ that if `\nu` is the
        Mullineux conjugate of the Kleshchev partition tuple `\mu` then the
        simple module `D^\nu =(D^\mu)^{\text{sgn}}` is obtained from `D^\mu`
        by twisting by the `\text{sgn}`-automorphism with is the Hecke algebra
        analogue of tensoring with the one dimensional sign representation.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, [0,1])
            sage: mc = KP([[4, 2], [5, 3, 1]]).mullineux_conjugate(); mc
            ([2, 2, 1, 1], [3, 2, 2, 1, 1])
            sage: mc.parent()
            Kleshchev partitions with e=3 and multicharge=(0,2)

        """
        P = self.parent()
        if self.size() == 0:
            size = None
            if isinstance(P, KleshchevPartitions_size):
                size = P._size
            KP = KleshchevPartitions(P._e, [-c for c in P._multicharge],
                                     size=size, convention=P._convention)
            return KP.element_class(KP, [[]]*P._level)

        good_cells = self.good_cells()
        assert good_cells

        k,r,c = sorted(good_cells.values())[0]
        # This is technically wrong when the parent has a fixed size because
        #   the resulting Kleshchev partition after removing a cell has abs
        #   smaller size. However, this is useful to avoid constructing
        #   transient parents.
        mu = P.element_class(P, self.remove_cell(k,r,c)).mullineux_conjugate()
        # add back on a cogood cell of residue -residue(k,r,c)
        KP = mu.parent()
        return KP.element_class(KP, mu.add_cell(*mu.cogood_cells( r-c-self.parent()._multicharge[k])))

    def is_regular(self):
        r"""
        Return ``True`` if ``self`` is a `e`-regular partition tuple.

        A partition tuple is `e`-regular if we can get to the
        empty partition tuple by successively removing a sequence
        of good cells in the down direction.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2, [0,2], convention="right restricted")
            sage: KP([[3,2,1], [2,1,1]]).is_regular()
            False
            sage: KP = KleshchevPartitions(4, [0,2], convention="right restricted")
            sage: KP([[3,2,1], [2,1,1]]).is_regular()
            True
            sage: KP([[], []]).is_regular()
            True
        """
        if self.size() == 0:
            return True
        KP = self.parent()
        return _is_regular(self.to_list(), KP._multicharge, KP._convention)

    def is_restricted(self):
        r"""
        Return ``True`` if ``self`` is an `e`-restricted partition tuple.

        A partition tuple is `e`-restricted if we can get to the
        empty partition tuple by successively removing a sequence
        of good cells in the up direction.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2, [0,2], convention="left regular")
            sage: KP([[3,2,1], [3,1]]).is_restricted()
            False
            sage: KP = KleshchevPartitions(3, [0,2], convention="left regular")
            sage: KP([[3,2,1], [3,1]]).is_restricted()
            True
            sage: KP([[], []]).is_restricted()
            True
        """
        if self.size() == 0:
            return True
        KP = self.parent()
        return _is_restricted(self.to_list(), KP._multicharge, KP._convention)

class KleshchevCrystalMixin(object):
    """
    Mixin class for the crystal structure of a Kleshchev partition.
    """
    def epsilon(self, i):
        r"""
        Return the Kashiwara crystal operator `\varepsilon_i` applied to ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: [x.epsilon(i) for i in C.index_set()]
            [0, 3, 0]
        """
        return len(self.normal_cells(i))

    def phi(self, i):
        r"""
        Return the Kashiwara crystal operator `\varphi_i` applied to ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: [x.phi(i) for i in C.index_set()]
            [3, 2, 0]
        """
        return len(self.conormal_cells(i))

    def Epsilon(self):
        r"""
        Return `\varepsilon` of ``self``.

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.Epsilon()
            3*Lambda[1]
        """
        P = self.parent()
        WLR = P.weight_lattice_realization()
        La = WLR.fundamental_weights()
        n = self.normal_cells()
        return WLR.sum(len(n[i])*La[i] for i in P.index_set() if i in n)

    def Phi(self):
        r"""
        Return `\phi` of ``self``.

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.Phi()
            3*Lambda[0] + 2*Lambda[1]
        """
        P = self.parent()
        WLR = P.weight_lattice_realization()
        La = WLR.fundamental_weights()
        c = self.conormal_cells()
        return WLR.sum(len(c[i])*La[i] for i in P.index_set() if i in c)

    def weight(self):
        r"""
        Return the weight of ``self``.

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: x = C([[5,4,1], [3,2,1,1]])
            sage: x.weight()
            3*Lambda[0] - Lambda[1] - 5*delta
            sage: x.Phi() - x.Epsilon()
            3*Lambda[0] - Lambda[1]

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="right regular")
            sage: y = C([[5,1,1], [4,2,2,1,1]])
            sage: y.weight()
            6*Lambda[0] - 4*Lambda[1] - 4*delta
            sage: y.Phi() - y.Epsilon()
            6*Lambda[0] - 4*Lambda[1]

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: y = C([[5,1,1], [4,2,2,1,1]])
            sage: y.weight()
            6*Lambda[0] - 4*Lambda[1] - 4*delta
            sage: y.Phi() - y.Epsilon()
            6*Lambda[0] - 4*Lambda[1]
        """
        WLR = self.parent().weight_lattice_realization()
        alpha = WLR.simple_roots()
        La = WLR.fundamental_weights()
        r = self.parent()._multicharge
        wt = WLR.sum(La[ZZ(x)] for x in r)
        return wt - WLR.sum(alpha[self.content(*c, multicharge=r)]
                            for c in self.cells())

class KleshchevPartitionCrystal(KleshchevPartition, KleshchevCrystalMixin):
    """
    Kleshchev partition with the crystal structure.
    """
    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, convention="left regular")
            sage: x = C([5,4,1])
            sage: x.e(0)
            sage: x.e(1)
            [5, 4]
        """
        P = self.parent()
        cell = self.good_cells(i)
        if cell is None:
            return None
        r,c = cell
        mu = list(self)
        mu[r] -= 1
        return type(self)(P, mu)

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, convention="left regular")
            sage: x = C([5,4,1])
            sage: x.f(0)
            [5, 5, 1]
            sage: x.f(1)
            sage: x.f(2)
            [5, 4, 2]
        """
        P = self.parent()
        cell = self.cogood_cells(i)
        if cell is None:
            return None
        r,c = cell
        mu = list(self)
        if c == 0:
            mu.append(1)
        else:
            mu[r] += 1
        return type(self)(P, mu)

class KleshchevPartitionTupleCrystal(KleshchevPartitionTuple, KleshchevCrystalMixin):
    """
    Kleshchev partition tuple with the crystal structure.
    """
    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.e(0)
            sage: x.e(1)
            ([5, 4, 1], [2, 2, 1, 1])
        """
        P = self.parent()
        cell = self.good_cells(i)
        if cell is None:
            return None
        k,r,c = cell
        mu = self.to_list()
        mu[k][r] -= 1
        return type(self)(P, mu)

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention="left regular")
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.f(0)
            ([5, 5, 1], [3, 2, 1, 1])
            sage: x.f(1)
            ([5, 4, 1], [3, 2, 2, 1])
            sage: x.f(2)
        """
        P = self.parent()
        cell = self.cogood_cells(i)
        if cell is None:
            return None
        k,r,c = cell
        mu = self.to_list()
        if c == 0:
            mu[k].append(1)
        else:
            mu[k][r] += 1
        return type(self)(P, mu)

#--------------------------------------------------
# Kleshchev partitions - parent classes
#--------------------------------------------------


class KleshchevPartitions(PartitionTuples):
    r"""
    Kleshchev partitions

    A partition (tuple) `\mu` is Kleshchev if it can be recursively
    obtained by adding a sequence of good nodes to the empty
    :class:`PartitionTuple` of the same :meth:`~PartitionTuple.level`
    and multicharge.

    There are four different conventions that are used in the literature for
    Kleshchev partitions, depending on whether we read partitions from top
    to bottom (regular) or bottom to top (restricted) and whether we read
    partition tuples from left to right or right to left. All of these
    conventions are supported::

        sage: KleshchevPartitions(2, [0,0], size=2, convention='left regular')[:]
        [([1], [1]), ([2], [])]
        sage: KleshchevPartitions(2, [0,0], size=2, convention='left restricted')[:]
        [([1], [1]), ([], [1, 1])]
        sage: KleshchevPartitions(2, [0,0], size=2, convention='right regular')[:]
        [([1], [1]), ([], [2])]
        sage: KleshchevPartitions(2, [0,0], size=2, convention='right restricted')[:]
        [([1], [1]), ([1, 1], [])]

    By default, the ``left restricted`` convention is used. As a shorthand,
    ``LG``, ``LS``, ``RG`` and ``RS``, respectively, can be used to specify
    the ``convention``. With the ``left`` convention the partition tuples
    should be ordered with the most dominant partitions in the partition
    tuple on the left and with the ``right`` convention the most dominant
    partition is on the right.

    The :class:`~KleshchevPartitions` class will automatically convert
    between these four different conventions::

        sage: KPlg = KleshchevPartitions(2, [0,0], size=2, convention='left regular')
        sage: KPls = KleshchevPartitions(2, [0,0], size=2, convention='left restricted')
        sage: [KPlg(mu) for mu in KPls]
        [([1], [1]), ([2], [])]

    EXAMPLES::

        sage: sorted(KleshchevPartitions(5,[3,2,1],1, convention='RS'))
        [([], [], [1]), ([], [1], []), ([1], [], [])]
        sage: sorted(KleshchevPartitions(5, [3,2,1], 1, convention='LS'))
        [([], [], [1]), ([], [1], []), ([1], [], [])]
        sage: sorted(KleshchevPartitions(5, [3,2,1], 3))
        [([], [], [1, 1, 1]),
         ([], [], [2, 1]),
         ([], [], [3]),
         ([], [1], [1, 1]),
         ([], [1], [2]),
         ([], [1, 1], [1]),
         ([], [2], [1]),
         ([], [3], []),
         ([1], [], [1, 1]),
         ([1], [], [2]),
         ([1], [1], [1]),
         ([1], [2], []),
         ([1, 1], [1], []),
         ([2], [], [1]),
         ([2], [1], []),
         ([3], [], [])]
        sage: sorted(KleshchevPartitions(5, [3,2,1], 3, convention="left regular"))
        [([], [], [1, 1, 1]),
         ([], [1], [1, 1]),
         ([], [1], [2]),
         ([], [1, 1], [1]),
         ([], [1, 1, 1], []),
         ([1], [], [1, 1]),
         ([1], [1], [1]),
         ([1], [1, 1], []),
         ([1], [2], []),
         ([1, 1], [], [1]),
         ([1, 1], [1], []),
         ([1, 1, 1], [], []),
         ([2], [], [1]),
         ([2], [1], []),
         ([2, 1], [], []),
         ([3], [], [])]

    REFERENCES:

    - [AM2000]_
    - [Ariki2001]_
    - [BK2009]_
    - [Kle2009]_
    """
    @staticmethod
    def __classcall_private__(cls, e, multicharge=(0,), size=None,
                              convention="left restricted"):
        r"""
        This is a factory class which returns the appropriate parent based on
        the values of `level` and `size`.

        EXAMPLES::

            sage: sorted(KleshchevPartitions(5, [3,2,1], 1, convention='RS'))
            [([], [], [1]), ([], [1], []), ([1], [], [])]
            sage: sorted(KleshchevPartitions(5, [3,2,1], 1, convention='LS'))
            [([], [], [1]), ([], [1], []), ([1], [], [])]
        """
        if size is None and multicharge in ZZ:
            size = ZZ(multicharge)
            multicharge = (0,)

        I = IntegerModRing(e)
        multicharge = tuple([I(x) for x in multicharge])

        convention = convention.upper()
        if 'S' in convention:
            convention = convention[0] + 'S'
        elif 'G' in convention:
            convention = convention[0] + 'G'
        if convention not in ['RG','LG', 'RS', 'LS']:
            raise ValueError('invalid convention')

        if size is None:
            return KleshchevPartitions_all(e, multicharge, convention)

        return KleshchevPartitions_size(e, multicharge, size, convention)

    def multicharge(self):
        """
        Return the multicharge of ``self``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(6, [2])
            sage: KP.multicharge()
            (2,)
            sage: KP = KleshchevPartitions(5, [3,0,1], 1, convention='LS')
            sage: KP.multicharge()
            (3, 0, 1)
            """
        return self._multicharge

    def convention(self):
        """
        Return the convention of ``self``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(4)
            sage: KP.convention()
            'restricted'
            sage: KP = KleshchevPartitions(6, [4], 3, convention="right regular")
            sage: KP.convention()
            'regular'
            sage: KP = KleshchevPartitions(5, [3,0,1], 1)
            sage: KP.convention()
            'left restricted'
            sage: KP = KleshchevPartitions(5, [3,0,1], 1, convention='right regular')
            sage: KP.convention()
            'right regular'
        """
        if self._convention[1] == 'S':
            convention = "restricted"
        else:
            convention = "regular"

        if self._level == 1:
            return convention

        if self._convention[0] == 'R':
            return "right " + convention

        return "left " + convention

    def _element_constructor_(self, mu):
        r"""
        Return ``mu`` as an element of :class:`~KleshchevPartitions`, or
        or raise an error if ``mu`` is not a Kleshchev partition (tuple).

        The main purpose of the element constructor code is to allow automatic
        conversion between the four possible conventions for Kleshchev
        partitions.

        EXAMPLES::

            sage: KPlg = KleshchevPartitions(2, [0,0], size=2, convention='left regular')
            sage: KPls = KleshchevPartitions(2, [0,0], size=2, convention='left restricted')
            sage: [KPlg(mu) for mu in KPls] # indirect doc test
            [([1], [1]), ([2], [])]

        """
        if isinstance(mu, (KleshchevPartition, KleshchevPartitionTuple)):
            KPmu = mu.parent()
            if KPmu == self:
                return mu

            if KPmu._level != self._level or KPmu._e != self._e:
                raise ValueError('%s is not an element of %s'%(mu, self))

            if KPmu._convention[1] != self._convention[1]:
                mu = [nu.conjugate() for nu in mu]
                if self._level>1 and KPmu._convention[0] == self._convention[0]:
                    mu = mu[::-1]

        return super(KleshchevPartitions, self)._element_constructor_(mu)

class KleshchevPartitions_all(KleshchevPartitions):
    r"""
    Class of all Kleshchev partitions.

    .. RUBRIC:: Crystal structure

    We consider type `A_{e-1}^{(1)}` crystals, and let `r = (r_i |
    r_i \in \ZZ / e \ZZ)` be a finite sequence of length `k`, which
    is the *level*, and `\lambda = \sum_i \Lambda_{r_i}`. We will
    model the highest weight `U_q(\mathfrak{g})`-crystal `B(\lambda)`
    by a particular subset of partition tuples of level `k`.

    Consider a partition tuple `\mu` with multicharge `r`.
    We define `e_i(\mu)` as the partition tuple obtained after the
    deletion of the `i`-:meth:`good cell
    <~sage.combinat.partition_kleshchev.KleshchevPartitionTuple.good_cell>`
    to `\mu` and `0` if there is no `i`-good cell. We define `f_i(\mu)` as
    the partition tuple obtained by the addition of the `i`-:meth:`cogood cell
    <~sage.combinat.partition_kleshchev.KleshchevPartitionTuple.cogood_cell>`
    to `\mu` and `0` if there is no `i`-good cell.

    The crystal `B(\lambda)` is the crystal generated by the empty
    partition tuple. We can compute the weight of an element `\mu` by taking
    `\lambda - \sum_{i=0}^n c_i \alpha_i` where `c_i` is the number of cells
    of `n`-residue `i` in `\mu`. Partition tuples in the crystal are known
    as *Kleshchev partitions*.

    .. NOTE::

        We can describe normal (not restricted) Kleshchev partition tuples
        in `B(\lambda)` as partition tuples `\mu` such that
        `\mu^{(t)}_{r_t - r_{t+1} + x} < \mu^{(t+1)}_x`
        for all `x \geq 1` and `1 \leq t \leq k - 1`.

    INPUT:

    - ``e`` -- for type `A_{e-1}^{(1)}` or `0`
    - ``multicharge`` -- the multicharge sequence `r`
    - ``convention`` -- (default: ``'LS'``) the reading convention

    EXAMPLES:

    We first do an example of a level 1 crystal::

        sage: C = crystals.KleshchevPartitions(3, [0], convention="left restricted")
        sage: C
        Kleshchev partitions with e=3
        sage: mg = C.highest_weight_vector()
        sage: mg
        []
        sage: mg.f(0)
        [1]
        sage: mg.f(1)
        sage: mg.f(2)
        sage: mg.f_string([0,2,1,0])
        [1, 1, 1, 1]
        sage: mg.f_string([0,1,2,0])
        [2, 2]
        sage: GC = C.subcrystal(max_depth=5).digraph()
        sage: B = crystals.LSPaths(['A',2,1], [1,0,0])
        sage: GB = B.subcrystal(max_depth=5).digraph()
        sage: GC.is_isomorphic(GB, edge_labels=True)
        True

    Now a higher level crystal::

        sage: C = crystals.KleshchevPartitions(3, [0,2], convention="right restricted")
        sage: mg = C.highest_weight_vector()
        sage: mg
        ([], [])
        sage: mg.f(0)
        ([1], [])
        sage: mg.f(2)
        ([], [1])
        sage: mg.f_string([0,1,2,0])
        ([2, 2], [])
        sage: mg.f_string([0,2,1,0])
        ([1, 1, 1, 1], [])
        sage: mg.f_string([2,0,1,0])
        ([2], [2])
        sage: GC = C.subcrystal(max_depth=5).digraph()
        sage: B = crystals.LSPaths(['A',2,1], [1,0,1])
        sage: GB = B.subcrystal(max_depth=5).digraph()
        sage: GC.is_isomorphic(GB, edge_labels=True)
        True

    The ordering of the residues gives a different representation of the
    higher level crystals (but it is still isomorphic)::

        sage: C2 = crystals.KleshchevPartitions(3, [2,0], convention="right restricted")
        sage: mg2 = C2.highest_weight_vector()
        sage: mg2.f_string([0,1,2,0])
        ([2], [2])
        sage: mg2.f_string([0,2,1,0])
        ([1, 1, 1], [1])
        sage: mg2.f_string([2,0,1,0])
        ([2, 1], [1])
        sage: GC2 = C2.subcrystal(max_depth=5).digraph()
        sage: GC.is_isomorphic(GC2, edge_labels=True)
        True

    TESTS:

    We check that all conventions give isomorphic crystals::

        sage: CLS = crystals.KleshchevPartitions(3, [2,0], convention="left restricted")
        sage: CRS = crystals.KleshchevPartitions(3, [2,0], convention="right restricted")
        sage: CLG = crystals.KleshchevPartitions(3, [2,0], convention="left regular")
        sage: CRG = crystals.KleshchevPartitions(3, [2,0], convention="right regular")
        sage: C = [CLS, CRS, CLG, CRG]
        sage: G = [B.subcrystal(max_depth=6).digraph() for B in C]
        sage: G[0].is_isomorphic(G[1], edge_labels=True)
        True
        sage: G[0].is_isomorphic(G[2], edge_labels=True)
        True
        sage: G[0].is_isomorphic(G[3], edge_labels=True)
        True

    REFERENCES:

    - [Ariki1996]_
    - [Ariki2001]_
    - [Tingley2007]_
    - [TingleyLN]_
    - [Vazirani2002]_
    """
    def __init__(self, e, multicharge, convention):
        r"""
        Initializes ``self``.

        EXAMPLES::

            sage: K = KleshchevPartitions(4, [2])
            sage: TestSuite(K).run()  # long time
            sage: K = KleshchevPartitions(4, [0,2,1])
            sage: TestSuite(K).run()  # long time

            sage: K = KleshchevPartitions(0, [2])
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(0, [0,2,1])
            sage: TestSuite(K).run()  # long time
        """
        if e not in NN or e == 1:
            raise ValueError('e must belong to {0,2,3,4,5,6,...}')
        if e > 0:
            from sage.combinat.root_system.cartan_type import CartanType
            from sage.categories.highest_weight_crystals import HighestWeightCrystals
            from sage.categories.regular_crystals import RegularCrystals
            self._cartan_type = CartanType(['A', e-1, 1])
            cat = (HighestWeightCrystals(), RegularCrystals().Infinite())
        else:
            cat = InfiniteEnumeratedSets()

        self._level = len(multicharge)
        if self._level == 1:
            self.Element = KleshchevPartitionCrystal
            self._element_constructor_ = getattr_from_other_class(self, Partitions, '_element_constructor_')
        else:
            self.Element = KleshchevPartitionTupleCrystal

        super(KleshchevPartitions_all, self).__init__(category=cat)
        self._e = e   # for printing
        self._index_set = IntegerModRing(e)
        self._multicharge = multicharge
        self._convention = convention
        if e > 0:
            if self._level == 1:
                self.module_generators = (self.element_class(self, []),)
            else:
                self.module_generators = (self.element_class(self, [[]]*self._level),)

    def _repr_(self):
        """
        EXAMPLES::

            sage: KleshchevPartitions(4, [2])
            Kleshchev partitions with e=4
            sage: KleshchevPartitions(3,[0,0,0])
            Kleshchev partitions with e=3 and multicharge=(0,0,0)
            sage: KleshchevPartitions(3,[0,0,1])
            Kleshchev partitions with e=3 and multicharge=(0,0,1)
        """
        if self._level == 1:
            return 'Kleshchev partitions with e=%s' % (self._e)

        return 'Kleshchev partitions with e=%s and multicharge=(%s)' % (
                        self._e,','.join('%s'%m for m in self._multicharge))

    def __contains__(self, mu):
        """
        Containment test for Kleshchev partitions.

        EXAMPLES::

            sage: PartitionTuple([[3,2],[2]]) in KleshchevPartitions(2, [0,0], 7)
            False
            sage: PartitionTuple([[],[2,1],[3,2]]) in KleshchevPartitions(5, [0,0,1], 7)
            False
            sage: PartitionTuple([[],[2,1],[3,2]]) in KleshchevPartitions(5, [0,1,1], 7)
            False
            sage: PartitionTuple([[],[2,1],[3,2]]) in KleshchevPartitions(5, [0,1,1], 8)
            True
            sage: all(mu in PartitionTuples(3,8) for mu in KleshchevPartitions(2, [0,0,0], 8))
            True
        """
        if isinstance(mu, (KleshchevPartition, KleshchevPartitionTuple)):
            if mu.level() != self._level:
                return False
            mu = self.element_class(self, list(mu))
            if self._convention[1] == 'G':
                return mu.is_regular()

            return mu.is_restricted()

        try:
            mu = self.element_class(self, mu)
        except ValueError:
            return False
        return mu in self

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: it = iter(KleshchevPartitions(2))
            sage: [next(it) for _ in range(10)]
            [[], [1], [1, 1], [2, 1], [1, 1, 1], [2, 1, 1],
             [1, 1, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: it = iter(KleshchevPartitions(2, convention='LG'))
            sage: [next(it) for _ in range(10)]
            [[], [1], [2], [3], [2, 1], [4], [3, 1], [5], [4, 1], [3, 2]]

            sage: it = iter(KleshchevPartitions(2, [0,1], convention='LS'))
            sage: [next(it) for _ in range(10)]
            [([], []),
             ([1], []),
             ([], [1]),
             ([1], [1]),
             ([], [1, 1]),
             ([1, 1], [1]),
             ([1], [1, 1]),
             ([], [2, 1]),
             ([], [1, 1, 1]),
             ([2, 1], [1])]
            sage: it = iter(KleshchevPartitions(2, [0,1], convention='RS'))
            sage: [next(it) for _ in range(10)]
            [([], []),
             ([1], []),
             ([], [1]),
             ([1, 1], []),
             ([1], [1]),
             ([2, 1], []),
             ([1, 1, 1], []),
             ([1, 1], [1]),
             ([1], [1, 1]),
             ([2, 1, 1], [])]
            sage: it = iter(KleshchevPartitions(2, [0,1], convention='LG'))
            sage: [next(it) for _ in range(10)]
            [([], []),
             ([1], []),
             ([], [1]),
             ([2], []),
             ([1], [1]),
             ([3], []),
             ([2, 1], []),
             ([2], [1]),
             ([1], [2]),
             ([4], [])]
            sage: it = iter(KleshchevPartitions(2, [0,1], convention='RG'))
            sage: [next(it) for _ in range(10)]
            [([], []),
             ([1], []),
             ([], [1]),
             ([1], [1]),
             ([], [2]),
             ([2], [1]),
             ([1], [2]),
             ([], [3]),
             ([], [2, 1]),
             ([2, 1], [1])]

            sage: it = iter(KleshchevPartitions(3, [0,1,2]))
            sage: [next(it) for _ in range(10)]
            [([], [], []), ([1], [], []), ([], [1], []), ([], [], [1]),
             ([1], [1], []), ([1], [], [1]), ([], [1, 1], []),
             ([], [1], [1]), ([], [], [2]), ([], [], [1, 1])]
        """
        # This is a modified form of what appears in the fixed size code
        if self._level == 1:
            if self._e == 0:
                P = Partitions()
            elif self._convention[1] == 'G':
                P = Partitions(regular=self._e)
            else:
                P = Partitions(restricted=self._e)

            for mu in P:
                yield self.element_class(self, list(mu))
        else:
            next_level = [self.element_class(self, [[]]*len(self._multicharge))]
            while True:
                cur = next_level
                next_level = []
                for mu in cur:
                    yield mu
                    mu_list = mu.to_list()
                    for cell in sorted(mu.cogood_cells().values()):
                        data = [list(p) for p in mu_list]
                        k, r, c = cell
                        if c == 0:
                            data[k].append(1)
                        else:
                            data[k][r] += 1
                        nu = self.element_class(self, data)
                        good_cells = nu.good_cells().values()
                        if self._convention[1] == "S":
                            if all(cell >= c for c in good_cells):
                                next_level.append(nu)
                        else:
                            if all(cell <= c for c in good_cells):
                                next_level.append(nu)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: KleshchevPartitions(3, [0,0,0,0], size=4).an_element()
            ([1], [1], [1], [1])
        """
        return self[12]


class KleshchevPartitions_size(KleshchevPartitions):
    """
    Kleshchev partitions of a fixed size.
    """
    def __init__(self, e, multicharge=(0,), size=0, convention='RS'):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: K = KleshchevPartitions(4, 2)
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(4, 4, convention='left regular')
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(4, 4, convention='left restricted')
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(4, [0,2,1], 4, convention='left regular')
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(0, 2, convention='right restricted')
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(0, [0,2,1], 4, convention='left restricted')
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(0, [0,2,1], 4, convention='left regular')
            sage: TestSuite(K).run()

        We verify that we obtain the same size for all conventions
        and that the result is equal to the number of elements in
        the crystal at the corresponding depth::

            sage: B = crystals.LSPaths(['A',2,1], [1,0,1])
            sage: nd4 = (B.subcrystal(max_depth=4).cardinality()
            ....:        - B.subcrystal(max_depth=3).cardinality())
            sage: K = KleshchevPartitions(3, [0,2], 4, convention='RS')
            sage: K.cardinality() == nd4
            True
            sage: K = KleshchevPartitions(3, [0,2], 4, convention='RG')
            sage: K.cardinality() == nd4
            True
            sage: K = KleshchevPartitions(3, [0,2], 4, convention='LS')
            sage: K.cardinality() == nd4
            True
            sage: K = KleshchevPartitions(3, [0,2], 4, convention='LG')
            sage: K.cardinality() == nd4
            True
        """
        self._level = len(multicharge)
        if self._level == 1:
            self.Element = KleshchevPartition
            self._element_constructor_ = getattr_from_other_class(self, Partitions, '_element_constructor_')
        else:
            self.Element = KleshchevPartitionTuple
        super(KleshchevPartitions_size, self).__init__(category=FiniteEnumeratedSets())
        self._size = size
        # As lists do not take negative indices the case e=0 needs to be handled
        # differently. Rather than doing this we set e equal to a "really big"
        # number. Mathematically, this is equivalent and it means that we don't
        # have an exception to cater for.
        self._e = e
        self._I = IntegerModRing(e)
        self._multicharge = tuple(self._I(m) for m in multicharge)
        self._convention = convention

    def _repr_(self):
        """
        EXAMPLES::

            sage: KleshchevPartitions(4, [0,0], 3)
            Kleshchev partitions with e=4 and multicharge=(0,0) and size 3
        """
        if self._level == 1:
            return 'Kleshchev partitions with e=%s and size %s' % (self._e, self._size)

        return 'Kleshchev partitions with e=%s and multicharge=(%s) and size %s' % (
            self._e,','.join('%s'%m for m in self._multicharge), self._size
        )

    def __contains__(self, mu):
        """
        Check if ``mu`` is in ``self``.

        TESTS::

            sage: PartitionTuple([[3,2],[2]]) in KleshchevPartitions(2,[0,0],7)
            False
            sage: PartitionTuple([[3,2],[],[],[],[2]]) in KleshchevPartitions(5,[0,0,0,0,0],7)
            False
            sage: PartitionTuple([[2,1],[],[1,1],[],[2]]) in KleshchevPartitions(5,[0,0,0,0,0],7, convention='RG')
            False
            sage: PartitionTuple([[2,1],[],[1,1],[],[3]]) in KleshchevPartitions(2,[0,0,0,0,0],9, convention='RS')
            False
            sage: all(mu in PartitionTuples(3,8) for mu in KleshchevPartitions(0,[0,0,0],8))
            True
        """
        if isinstance(mu, (KleshchevPartition, KleshchevPartitionTuple)):
            if not (mu.level() == self._level and mu.size() == self._size):
                return False
            mu = self.element_class(self, list(mu))
            if self._convention[1] == 'G':
                return mu.is_regular()

            return mu.is_restricted()

        try:
            mu = self.element_class(self, mu)
        except ValueError:
            return False
        return mu in self

    def __iter__level_one(self):
        r"""
        Iterate over all Kleshchev partitions of level one and a fixed size.

        EXAMPLES::

            sage: KleshchevPartitions(2,0)[:]  # indirect doctest
            [[]]
            sage: KleshchevPartitions(2,1)[:]  # indirect doctest
            [[1]]
            sage: KleshchevPartitions(2,2)[:]  # indirect doctest
            [[1, 1]]
            sage: KleshchevPartitions(3,2)[:]  # indirect doctest
            [[2], [1, 1]]
        """
        if self._size == 0:
            yield self.element_class(self, [])
        else:
            if self._e == 0:
                P = Partitions(self._size)
            elif self._convention[1] == 'G':
                P = Partitions(self._size, regular=self._e)
            else:
                P = Partitions(self._size, restricted=self._e)

            for mu in P:
                yield self.element_class(self, list(mu))

    def __iter__higher_levels(self):
        r"""
        Iterate over all KleshchevPartitions of a fixed level greater than 1
        and a fixed size.

        EXAMPLES::

            sage: KleshchevPartitions(2,[0,0],1)[:]  # indirect doctest
            [([], [1])]
            sage: KleshchevPartitions(2,[0,0],2)[:]  # indirect doctest
            [([1], [1]), ([], [1, 1])]
            sage: KleshchevPartitions(3,[0,0],2)[:]  # indirect doctest
            [([1], [1]), ([], [2]), ([], [1, 1])]
        """
        if self._size == 0:
            yield self.element_class(self, [[]]*len(self._multicharge))
            return

        # For higher levels we have to recursively construct the restricted partitions
        # by adding on co-good nodes to smaller restricted partition. To avoid over
        # counting we return a new restricted partition only if we added on its lowest
        # good node.
        for mu in KleshchevPartitions_size(self._e, self._multicharge,
                                           size=self._size-1,
                                           convention=self._convention):
            mu_list = mu.to_list()
            for cell in mu.cogood_cells().values():
                data = [list(p) for p in mu_list]
                k,r,c = cell
                if c == 0:
                    data[k].append(1)
                else:
                    data[k][r] += 1
                nu = self.element_class(self, data)
                good_cells = nu.good_cells().values()
                if self._convention[1] == "S":
                    if all(cell >= c for c in good_cells):
                        yield nu
                else:
                    if all(cell <= c for c in good_cells):
                        yield nu

    @lazy_attribute
    def __iter__(self):
        """
        Wrapper to return the correct iterator which is different for
        :class:`Partitions` (level 1) and for :class:PartitionTuples`
        (higher levels).

        EXAMPLES::

            sage: KleshchevPartitions(3, 3, convention='RS')[:]
            [[2, 1], [1, 1, 1]]
            sage: KleshchevPartitions(3, 3, convention='RG')[:]
            [[3], [2, 1]]
            sage: KleshchevPartitions(3, [0], 3)[:]
            [[2, 1], [1, 1, 1]]
            sage: KleshchevPartitions(3, [0,0], 3)[:]
            [([1], [2]), ([1], [1, 1]), ([], [2, 1]), ([], [1, 1, 1])]
            sage: KleshchevPartitions(2, [0,1], size=0)[:]
            [([], [])]
            sage: KleshchevPartitions(2, [0,1], size=1)[:]
            [([1], []), ([], [1])]
            sage: KleshchevPartitions(2, [0,1], size=2)[:]
            [([1], [1]), ([], [1, 1])]
            sage: KleshchevPartitions(3, [0,1,2], size=2)[:]
            [([1], [1], []), ([1], [], [1]), ([], [1, 1], []),
             ([], [1], [1]), ([], [], [2]), ([], [], [1, 1])]
        """
        if self.level() == 1:
            return self.__iter__level_one

        return self.__iter__higher_levels

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: KleshchevPartitions(4, [0,0,0,0], 4).an_element()
            ([1], [1], [1], [1])
            sage: KleshchevPartitions(4, [2], 4).an_element()
            [3, 1]
        """
        return self[0]

    Element = KleshchevPartitionTuple

#--------------------------------------------------
# helper functions
#--------------------------------------------------

def _a_good_cell(kpt, multicharge, convention):
    """
    Return a good cell from ``kpt`` considered as a Kleshchev partition
    with ``multicharge`` under ``convention``.

    EXAMPLES::

        sage: from sage.combinat.partition_kleshchev import _a_good_cell
        sage: I2 = IntegerModRing(2)
        sage: _a_good_cell([[3,2,1], [3,1,1]], [I2(0),I2(2)], 'RS')
        (1, 2, 0)
        sage: _a_good_cell([[3,1,1], [3,2]], [I2(0),I2(2)], 'LG')
        (1, 1, 1)
        sage: I3 = IntegerModRing(3)
        sage: _a_good_cell([[3,2,1], [3,1,1]], [I3(0),I3(2)], 'RS')
        (0, 2, 0)
        sage: _a_good_cell([[3,1,1], [3,2]], [I3(0),I3(2)], 'LG')
        (1, 0, 2)
        sage: _a_good_cell([[], []], [I3(0),I3(2)], 'RS') is None
        True
        sage: _a_good_cell([[1,1], []], [I3(0),I3(2)], 'LS') is None
        True
    """
    # We use a dictionary for the normal nodes as the indexing set is Z when e=0
    carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)
    ret = None

    if convention[0] == 'L':
        rows = [(k,r) for k,part in enumerate(kpt) for r in range(len(part)+1)]
    else:
        rows = [(k,r) for k,part in reversed(list(enumerate(kpt))) for r in range(len(part)+1)]
    if convention[1] == 'S':
        rows.reverse()

    for row in rows:
        k,r = row
        if r == len(kpt[k]): # addable cell at bottom of a component
            carry[multicharge[k]-r] += 1
        else:
            res = multicharge[k] + kpt[k][r] - r - 1
            if r == len(kpt[k])-1 or kpt[k][r] > kpt[k][r+1]: # removable cell
                if carry[res] == 0:
                    ret = (k, r, kpt[k][r]-1)
                else:
                    carry[res] -= 1
            if r == 0 or kpt[k][r-1] > kpt[k][r]:             # addable cell
                carry[res+1] += 1

    # finally return the result
    return ret

def _is_regular(kpt, multicharge, convention):
    """
    Return ``True`` if ``kpt`` is a ``multicharge``-regular
    Kleshchev partition.

    EXAMPLES::

        sage: from sage.combinat.partition_kleshchev import _is_regular
        sage: I2 = IntegerModRing(2)
        sage: _is_regular([[3,1,1], [3,2]], [I2(0),I2(2)], 'LS')
        False
        sage: I3 = IntegerModRing(3)
        sage: _is_regular([[3,1,1], [3,2]], [I3(0),I3(2)], 'LS')
        True
        sage: _is_regular([[], []], [I3(0),I3(2)], 'LS')
        True
    """
    if all(part == [] for part in kpt):
        return True
    convention = convention[0] + 'G'
    cell = _a_good_cell(kpt, multicharge, convention)
    while cell is not None:
        k,r,c = cell
        if kpt[k][r] == 1:
            kpt[k].pop()
        else:
            kpt[k][r] -= 1
        cell = _a_good_cell(kpt, multicharge, convention)
    return all(part == [] for part in kpt)

def _is_restricted(kpt, multicharge, convention):
    """
    Return ``True`` if ``kpt`` is an ``multicharge``-restricted
    Kleshchev partition.

    EXAMPLES::

        sage: from sage.combinat.partition_kleshchev import _is_restricted
        sage: I2 = IntegerModRing(2)
        sage: _is_restricted([[3,2,1], [3,1,1]], [I2(0),I2(2)], 'RG')
        False
        sage: I3 = IntegerModRing(3)
        sage: _is_restricted([[3,2,1], [3,1,1]], [I3(0),I3(2)], 'RG')
        True
        sage: _is_restricted([[], []], [I3(0),I3(2)], 'RG')
        True
    """
    if all(part == [] for part in kpt):
        return True
    convention = convention[0] + 'S'
    cell = _a_good_cell(kpt, multicharge, convention)
    while cell is not None:
        k,r,c = cell
        if kpt[k][r] == 1:
            kpt[k].pop()
        else:
            kpt[k][r] -= 1
        cell = _a_good_cell(kpt, multicharge, convention)
    return all(part == [] for part in kpt)
