# sage_setup: distribution = sagemath-combinat
r"""
Residue sequences of tableaux

A *residue sequence* for a :class:`~sage.combinat.tableau.StandardTableau`, or
:class:`~sage.combinat.tableau_tuple.StandardTableauTuple`, of size `n` is an
`n`-tuple `(i_1, i_2, \ldots, i_n)` of elements of `\ZZ / e\ZZ` for some
positive integer `e \ge 1`.  Such sequences arise in the representation
theory of the symmetric group and the closely related cyclotomic Hecke
algebras, and cyclotomic quiver Hecke algebras, where the residue sequences
play a similar role to weights in the representations of Lie groups and
Lie algebras. These Hecke algebras are semisimple when `e` is "large enough"
and in these cases residue sequences are essentially the same as content
sequences (see :meth:`sage.combinat.partition.Partition.content`) and it
is not difficult to see that residue sequences are in bijection with the
set of standard  tableaux. In the non-semisimple case, when `e` is "small",
different standard tableaux can have the same residue sequence. In this
case the residue sequences describe how to decompose modules into
generalised eigenspaces for the Jucys-Murphy elements for these algebras.

By definition, if `t` is a :class:`~sage.combinat.tableau.StandardTableau` of
size `n` then the residue sequence of `t` is the `n`-tuple `(i_1, \ldots, i_n)`
where `i_m = c - r + e\ZZ`, if `m` appears in row `r` and column `c` of `t`.
If `p` is prime then such sequence arise in the representation theory of the
symmetric group n characteristic `p`. More generally, `e`-residue sequences
arise in he representation theory of the Iwahori-Hecke algebra (see
:class:`~sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra`) the
symmetric group with Hecke parameter at an `e`-th root of unity.

More generally, the `e`-residue sequence of a
:class:`~sage.combinat.tableau.StandardTableau` of size `n` and level `l` is
the `n`-tuple `(i_1, \ldots, i_n)` determined by `e` and a *multicharge*
`\kappa = (\kappa_1, \ldots, \kappa_l)` by setting
`i_m = \kappa_k + c - r + e\ZZ`, if `m` appears in component `k`, row `r`
and column `c` of `t`.  These sequences arise in the representation theory
of the cyclotomic Hecke algebras of type A, which are also known
as Ariki-Koike algebras.

The residue classes are constructed from standard tableaux::

    sage: StandardTableau([[1,2],[3,4]]).residue_sequence(2)
    2-residue sequence (0,1,1,0) with multicharge (0)
    sage: StandardTableau([[1,2],[3,4]]).residue_sequence(3)
    3-residue sequence (0,1,2,0) with multicharge (0)

    sage: StandardTableauTuple([[[5]],[[1,2],[3,4]]]).residue_sequence(3,[0,0])
    3-residue sequence (0,1,2,0,0) with multicharge (0,0)
    sage: StandardTableauTuple([[[5]],[[1,2],[3,4]]]).residue_sequence(3,[0,1])
    3-residue sequence (1,2,0,1,0) with multicharge (0,1)
    sage: StandardTableauTuple([[[5]],[[1,2],[3,4]]]).residue_sequence(3,[0,2])
    3-residue sequence (2,0,1,2,0) with multicharge (0,2)

One of the most useful functions of a :class:`ResidueSequence` is that it can
return the :class:`~sage.combinat.tableau_tuple.StandardTableaux_residue` and
:class:`~sage.combinat.tableau_tuple.StandardTableaux_residue_shape` that
contain all of the tableaux with this residue sequence. Again, these are best
accessed via the standard tableaux classes::

    sage: res = StandardTableau([[1,2],[3,4]]).residue_sequence(2)
    sage: res.standard_tableaux()
    Standard tableaux with 2-residue sequence (0,1,1,0) and multicharge (0)
    sage: res.standard_tableaux()[:]
    [[[1, 2, 4], [3]],
    [[1, 2], [3, 4]],
    [[1, 2], [3], [4]],
    [[1, 3, 4], [2]],
    [[1, 3], [2, 4]],
    [[1, 3], [2], [4]]]
    sage: res.standard_tableaux(shape=[4])
    Standard (4)-tableaux with 2-residue sequence (0,1,1,0) and multicharge (0)
    sage: res.standard_tableaux(shape=[4])[:]
    []

    sage: res=StandardTableauTuple([[[5]],[[1,2],[3,4]]]).residue_sequence(3,[0,0])
    sage: res.standard_tableaux()
    Standard tableaux with 3-residue sequence (0,1,2,0,0) and multicharge (0,0)
    sage: res.standard_tableaux(shape=[[1],[2,2]])[:]
    [([[5]], [[1, 2], [3, 4]]), ([[4]], [[1, 2], [3, 5]])]

These residue sequences are particularly useful in the graded representation
theory of the cyclotomic KLR algebras and the cyclotomic Hecke algebras of type~A;
see [DJM1998]_ and [BK2009]_.

This module implements the following classes:

* :class:`ResidueSequence`
* :class:`ResidueSequences`

.. SEEALSO::

    * :class:`Partitions`
    * :class:`PartitionTuples`
    * :class:`~sage.combinat.tableau_tuple.StandardTableaux_residue`
    * :class:`~sage.combinat.tableau_tuple.StandardTableaux_residue_shape`
    * :class:`~sage.combinat.tableau_tuple.RowStandardTableauTuples_residue`
    * :class:`~sage.combinat.tableau_tuple.RowStandardTableauTuples_residue_shape`
    * :class:`StandardTableaux`
    * :class:`StandardTableauTuples`
    * :class:`Tableaux`
    * :class:`TableauTuples`

.. TODO::

    Strictly speaking this module implements residue sequences of
    type `A^{(1)}_e`. Residue sequences of other types also need
    to be implemented.

AUTHORS:

- Andrew Mathas (2016-07-01): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2012,2016 Andrew Mathas <andrew dot mathas at sydney dot edu dot au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.sets_cat import Sets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from .partition_tuple import PartitionTuple
from .tableau_tuple import (StandardTableaux_residue,
                            StandardTableaux_residue_shape,
                            RowStandardTableauTuples_residue,
                            RowStandardTableauTuples_residue_shape)

# -------------------------------------------------
# Residue sequences
# -------------------------------------------------


# needed for __classcall_private__
class ResidueSequence(ClonableArray,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A residue sequence.

    The *residue sequence* of a tableau `t` (of partition or partition tuple
    shape) is the sequence `(i_1, i_2, \ldots, i_n)` where `i_k` is the
    residue of `l` in `t`, for `k = 1, 2, \ldots, n`, where `n` is the
    size of `t`. Residue sequences are important in the representation
    theory of the cyclotomic Hecke algebras of type `G(r, 1, n)`, and
    of the cyclotomic quiver Hecke algebras, because they determine the
    eigenvalues of the Jucys-Murphy elements upon all modules. More precisely,
    they index and completely determine the irreducible representations
    of the (cyclotomic) Gelfand-Tsetlin algebras.

    Rather than being called directly, residue sequences are best accessed
    via the standard tableaux classes :class:`StandardTableau` and
    :class:`StandardTableauTuple`.

    INPUT:

    Can be of the form:

    - ``ResidueSequence(e, res)``,
    - ``ResidueSequence(e, multicharge, res)``,

    where ``e`` is a positive integer not equal to 1 and ``res`` is a
    sequence of integers (the residues).

    EXAMPLES::

        sage: res = StandardTableauTuple([[[1,3],[6]],[[2,7],[4],[5]]]).residue_sequence(3,(0,5))
        sage: res
        3-residue sequence (0,2,1,1,0,2,0) with multicharge (0,2)
        sage: res.quantum_characteristic()
        3
        sage: res.level()
        2
        sage: res.size()
        7
        sage: res.residues()
        [0, 2, 1, 1, 0, 2, 0]
        sage: res.restrict(2)
        3-residue sequence (0,2) with multicharge (0,2)
        sage: res.standard_tableaux([[2,1],[1],[2,1]])
        Standard (2,1|1|2,1)-tableaux with 3-residue sequence (0,2,1,1,0,2,0) and multicharge (0,2)
        sage: res.standard_tableaux([[2,2],[3]]).list()
        []
        sage: res.standard_tableaux([[2,2],[3]])[:]
        []
        sage: res.standard_tableaux()
        Standard tableaux with 3-residue sequence (0,2,1,1,0,2,0) and multicharge (0,2)
        sage: res.standard_tableaux()[:10]
        [([[1, 3, 6, 7], [2, 5], [4]], []),
         ([[1, 3, 6], [2, 5], [4], [7]], []),
         ([[1, 3], [2, 5], [4, 6], [7]], []),
         ([[1, 3], [2, 5], [4], [7]], [[6]]),
         ([[1, 3], [2, 5], [4]], [[6, 7]]),
         ([[1, 3, 6, 7], [2], [4], [5]], []),
         ([[1, 3, 6], [2, 7], [4], [5]], []),
         ([[1, 3], [2, 7], [4], [5], [6]], []),
         ([[1, 3], [2, 7], [4], [5]], [[6]]),
         ([[1, 3], [2], [4], [5]], [[6, 7]])]

    The TestSuite fails ``_test_pickling`` because ``__getitem__`` does
    not support slices, so we skip this.

    TESTS::

        sage: from sage.combinat.tableau_residues import ResidueSequence
        sage: TestSuite( ResidueSequence(3,(0,0,1), [0,1,2])).run(skip='_test_pickling')
    """
    @staticmethod
    def __classcall_private__(cls, e, multicharge, residues=None, check=True):
        r"""
        Magic to allow class to accept a list (which is not hashable) instead
        of a partition (which is). At the same time we ensure that every
        residue sequence is constructed as an ``element_class`` call of
        an appropriate parent.

        The ``residues`` must always be specified and, instead, it is the
        ``multicharge`` which is the optional argument with default ``[0]``.
        This means that we have to perform some tricks when ``residues``
        is ``None``.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, [0,0,1], [0,0,1,1,2,2,3,3])  # indirect doctest
            3-residue sequence (0,0,1,1,2,2,0,0) with multicharge (0,0,1)
        """
        # if the multicharge is omitted it defaults to (0,) in level 1
        if residues is None:
            residues = multicharge
            multicharge = (0,)
        multicharge = tuple(multicharge)
        return ResidueSequences(e, multicharge).element_class(ResidueSequences(e, multicharge), tuple(residues), check)

    def __init__(self, parent, residues, check):
        r"""
        Initialize ``self``.

        The ``multicharge`` is the optional argument which, if omitted,
        defaults to ``(0,)``. On the other hand, the ``residue`` must
        always be specified so, below, we check to see whether or note
        ``residues`` is ``None`` and adjust accordingly in this case.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, (0,0,1), [0,0,1,1,2,2,3,3])
            3-residue sequence (0,0,1,1,2,2,0,0) with multicharge (0,0,1)

        The TestSuite fails ``_test_pickling`` because ``__getitem__`` does
        not support slices, so we skip this.

        TESTS::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: TestSuite(ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3])).run(skip='_test_pickling')
            sage: TestSuite( ResidueSequence(3, [0,1,2])).run(skip='_test_pickling')
            sage: TestSuite( ResidueSequence(3, [0], [0,1,2])).run(skip='_test_pickling')
            sage: TestSuite( ResidueSequence(3, [0,0], [0,0,1,2])).run(skip='_test_pickling')
            sage: TestSuite( ResidueSequence(3, [0,0,1,2])).run(skip='_test_pickling')
        """
        residues = tuple(parent._base_ring(i) for i in residues)
        super().__init__(parent, residues, check)

    def check(self):
        r"""
        Raise a :exc:`ValueError` if ``self`` is not a residue sequence.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, [0,0,1], [0,0,1,1,2,2,3,3]).check()
            sage: ResidueSequence(3, [0,0,1], [2,0,1,1,2,2,3,3]).check()
        """
        self.parent().check_element(self)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3])
            3-residue sequence (0,0,1,1,2,2,0,0) with multicharge (0,0,1)
        """
        return self.__str__()

    def __str__(self, join='with'):
        r"""
        The string representation of a residue sequence is a comma separated
        tuple with no spaces.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]).__str__()
            '3-residue sequence (0,0,1,1,2,2,0,0) with multicharge (0,0,1)'
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]).__str__('and')
            '3-residue sequence (0,0,1,1,2,2,0,0) and multicharge (0,0,1)'
        """
        string = '{e}-residue sequence ({res}) {join} multicharge ({charge})'
        return string.format(e=self.quantum_characteristic(),
                             res=','.join('%s' % r for r in self), join=join,
                             charge=','.join('%s' % r for r in self.multicharge()))

    def __getitem__(self, k):
        r"""
        Return the ``k``-th residue.

        INPUT:

        - ``k`` --- an integer between 1 and the length of the residue
          sequence ``self``

        The ``k``-th residue is the ``e``-residue (see
        :meth:`sage.combinat.tableau.StandardTable.residue`) of the
        integer ``k`` in some standard tableaux. As the entries of standard
        tableaux are always between `1` and `n`, the size of the tableau,
        the integer ``k`` must also be in this range (that is, this
        is **not** 0-based!).

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3])[4]
            1
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3])[7]
            0
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3])[9]
            Traceback (most recent call last):
            ...
            IndexError: k must be in the range 1, 2, ..., 8
        """
        try:
            return ClonableArray.__getitem__(self, k - 1)
        except (IndexError, KeyError):
            raise IndexError('k must be in the range 1, 2, ..., {}'.format(len(self)))

    def residues(self) -> list:
        r"""
        Return a list of the residue sequence.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]).residues()
            [0, 0, 1, 1, 2, 2, 0, 0]
        """
        return list(self)

    def restrict(self, m):
        r"""
        Return the subsequence of this sequence of length `m`.

        The residue sequence ``self`` is of the form `(r_1, \ldots, r_n)`.
        The function returns the residue sequence `(r_1, \ldots, r_m)`, with
        the same  :meth:`quantum_characteristic` and :meth:`multicharge`.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]).restrict(7)
            3-residue sequence (0,0,1,1,2,2,0) with multicharge (0,0,1)
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]).restrict(6)
            3-residue sequence (0,0,1,1,2,2) with multicharge (0,0,1)
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]).restrict(4)
            3-residue sequence (0,0,1,1) with multicharge (0,0,1)
        """
        return ResidueSequence(self.quantum_characteristic(),
                               self.multicharge(), self.residues()[:m])

    def restrict_row(self, cell, row):
        r"""
        Return a residue sequence for the tableau obtained by swapping the row
        in ending in `cell` with the row that is `row` rows above it and which
        has the same length.

        The residue sequence ``self`` is of the form `(r_1, \ldots, r_n)`.
        The function returns the residue sequence `(r_1, \ldots, r_m)`, with
        the same  :meth:`quantum_characteristic` and :meth:`multicharge`.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, [0,1,2,2,0,1]).restrict_row((1,2),1)
            3-residue sequence (2,0,1,0,1) with multicharge (0)
            sage: ResidueSequence(3, [1,0], [0,1,2,2,0,1]).restrict_row((1,1,2),1)
            3-residue sequence (2,0,1,0,1) with multicharge (1,0)
        """
        residues = self.residues()  # residue sequence
        residues.reverse()          # reversed residue sequence

        if residues[0] + row == residues[0]:
            # if the residues in the two rows are the same we do not
            # need to do anything special
            return self.restrict(len(residues) - 1)

        # determine the sets of residues, one_res and two_res, that need to be
        # interchanged in order to swap the corresponding rows
        row_len = cell[-1]  # length of the row being swapped
        one_res = [0]       # last row of tableau will move
        two_res = [0]       # will prune this entry later
        try:
            for c in range(1, row_len + 1):
                # residues decrease by 1 from right to left in each row
                one_res.append(residues.index(residues[0] - c,
                                              one_res[c - 1] + 1))
            for c in range(row_len + 1):
                two_res.append(residues.index(residues[0] - c + row,
                                              two_res[c] + 1))
                while two_res[-1] in one_res:
                    # entries in one_res and two_res must be disjoint
                    two_res[-1] = residues.index(residues[0] - c + row,
                                                 two_res[-1] + 1)
        except ValueError:
            return None

        # now add row to the residues in two_res and subtract row from those in
        # one_res
        for c in range(row_len + 1):
            residues[one_res[c]] += row
            residues[two_res[c + 1]] -= row  # jump over two_res[0]

        # remove the first residue, reverse the order and return
        return ResidueSequence(self.quantum_characteristic(),
                               self.multicharge(),
                               residues[1:][::-1])

    def swap_residues(self, i, j):
        r"""
        Return the *new* residue sequence obtained by swapping the residues
        for ``i`` and ``j``.

        INPUT:

        - ``i``, ``j`` -- two integers between `1` and the length of
          the residue sequence

        If residue sequence ``self`` is of the form `(r_1, \ldots, r_n)`, and
        `i < j`, then the residue sequence
        `(r_1, \ldots, r_j, \ldots, r_i, \ldots, r_m)`, with the same
        :meth:`quantum_characteristic` and :meth:`multicharge`, is returned.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: res = ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]); res
            3-residue sequence (0,0,1,1,2,2,0,0) with multicharge (0,0,1)
            sage: ser = res.swap_residues(2,6); ser
            3-residue sequence (0,2,1,1,2,0,0,0) with multicharge (0,0,1)
            sage: res == ser
            False

        TESTS::

            sage: res.swap_residues(22,26)
            Traceback (most recent call last):
            ...
            IndexError: 22 and 26 must be between 1 and 8
        """
        with self.clone() as swap:
            try:
                # we have overridden __getitem__ so that indices are 1-based but
                # __setitem__ is still 0-based so we need to renormalise the LHS
                swap[i - 1], swap[j - 1] = self[j], self[i]
            except IndexError:
                raise IndexError('%s and %s must be between 1 and %s' % (i, j, self.size()))
        return swap

    def standard_tableaux(self, shape=None):
        r"""
        Return the residue-class of standard tableaux that have residue
        sequence ``self``.

        INPUT:

        - ``shape`` -- (optional) a partition or partition tuple of
          the correct level

        OUTPUT:

        An iterator for the standard tableaux with this residue sequence. If
        the ``shape`` is given then only tableaux of this shape are returned,
        otherwise all of the full residue-class of standard tableaux, or
        standard tableaux tuples, is returned. The residue sequence ``self``
        specifies the :meth:`multicharge` of the tableaux which, in turn,
        determines the :meth:`level` of the tableaux in the residue class.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,0),[0,1,2,0,1,2,0,1,2]).standard_tableaux()
            Standard tableaux with 3-residue sequence (0,1,2,0,1,2,0,1,2) and multicharge (0,0,0)
            sage: ResidueSequence(3,(0,0,0),[0,1,2,0,1,2,0,1,2]).standard_tableaux([[3],[3],[3]])
            Standard (3|3|3)-tableaux with 3-residue sequence (0,1,2,0,1,2,0,1,2) and multicharge (0,0,0)
        """
        if shape is None:
            return StandardTableaux_residue(residue=self)
        else:
            return StandardTableaux_residue_shape(residue=self,
                                                  shape=PartitionTuple(shape))

    def row_standard_tableaux(self, shape=None):
        r"""
        Return the residue-class of row standard tableaux that have residue
        sequence ``self``.

        INPUT:

        - ``shape`` -- (optional) a partition or partition tuple of
          the correct level

        OUTPUT:

        An iterator for the row standard tableaux with this residue sequence. If
        the ``shape`` is given then only tableaux of this shape are returned,
        otherwise all of the full residue-class of row standard tableaux, or row
        standard tableaux tuples, is returned. The residue sequence ``self``
        specifies the :meth:`multicharge` of the tableaux which, in turn,
        determines the :meth:`level` of the tableaux in the residue class.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,0),[0,1,2,0,1,2,0,1,2]).row_standard_tableaux()
            Row standard tableaux with 3-residue sequence (0,1,2,0,1,2,0,1,2) and multicharge (0,0,0)
            sage: ResidueSequence(3,(0,0,0),[0,1,2,0,1,2,0,1,2]).row_standard_tableaux([[3],[3],[3]])
            Row standard (3|3|3)-tableaux with 3-residue sequence (0,1,2,0,1,2,0,1,2) and multicharge (0,0,0)
        """
        if shape is None:
            return RowStandardTableauTuples_residue(residue=self)
        else:
            return RowStandardTableauTuples_residue_shape(residue=self, shape=PartitionTuple(shape))

    def negative(self):
        r"""
        Return the negative of the residue sequence ``self``.

        That is, if ``self`` is the residue sequence `(i_1, \ldots, i_n)`
        then return `(-i_1, \ldots, -i_n)`. Taking the negative residue
        sequences is a shadow of tensoring with the sign representation
        from the cyclotomic Hecke algebras of type `A`.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,[0,0,1],[0,0,1,1,2,2,3,3]).negative()
            3-residue sequence (0,0,2,2,1,1,0,0) with multicharge (0,0,1)
        """
        return ResidueSequence(self.quantum_characteristic(), self.multicharge(),
                               (self.base_ring()(-i) for i in self))

    def block(self):
        r"""
        Return a dictionary `\beta` that determines the block associated to
        the residue sequence ``self``.

        Two Specht modules for a cyclotomic Hecke algebra of type `A` belong to
        the same block, in this sense, if and only if the residue sequences of
        their standard tableaux have the same block in this sense.  The blocks
        of these algebras are actually indexed by positive roots in the root
        lattice of an affine special linear group. Instead of than constructing
        the root lattice, this method simply returns a dictionary `\beta` where
        the keys are residues `i` and where the value of the  key `i` is equal
        to the numbers of nodes in the residue sequence ``self`` that are equal
        to `i`. The dictionary `\beta` corresponds to the positive root:

        .. MATH::

            \sum_{i\in I} \beta_i \alpha_i \in Q^+,

        These positive roots also index the blocks of the cyclotomic KLR
        algebras of type `A`.

        We return a dictionary because when the :meth:`quantum_characteristic` is `0`,
        the Cartan type is `A_{\infty}`, in which case the simple roots are
        indexed by the integers, which is infinite.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, [0,0,0], [0,1,2,0,1,2,0,1,2]).block()
            {0: 3, 1: 3, 2: 3}
        """
        return {i: self.residues().count(i) for i in set(self.residues())}

    def base_ring(self):
        r"""
        Return the base ring for the residue sequence.

        If the :meth:`quantum_characteristic` of the residue sequence ``self``
        is `e` then the base ring for the sequence is `\ZZ / e\ZZ`,
        or `\ZZ` if `e=0`.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, (0,0,1), [0,0,1,1,2,2,3,3]).base_ring()
            Ring of integers modulo 3
        """
        return self.parent()._base_ring

    def quantum_characteristic(self):
        r"""
        Return the quantum characteristic of the residue sequence ``self``.

        The `e`-residue sequences are associated with a cyclotomic Hecke
        algebra that has a parameter `q` of *quantum characteristic* `e`.
        This is the smallest positive integer such that
        `1 + q + \cdots + q^{e-1} = 0`, or `e=0` if no such integer exists.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, (0,0,1), [0,0,1,1,2,2,3,3]).quantum_characteristic()
            3
        """
        return self.parent()._quantum_characteristic

    def multicharge(self):
        r"""
        Return the multicharge for the residue sequence ``self``.

        The `e`-residue sequences are associated with a cyclotomic Hecke
        algebra with Hecke parameter `q` of :meth:`quantum_characteristic` `e`
        and multicharge `(\kappa_1, \ldots, \kappa_l)`. This means that
        the cyclotomic parameters of the Hecke algebra are
        `q^{\kappa_1}, \ldots, q^{\kappa_l}`. Equivalently, the Hecke
        algebra is determined by the dominant weight

        .. MATH::

            \sum_{r \in \ZZ / e\ZZ} \kappa_r \Lambda_r \in P^+.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, (0,0,1), [0,0,1,1,2,2,3,3]).multicharge()
            (0, 0, 1)
        """
        return self.parent()._multicharge

    def level(self):
        r"""
        Return the level of the residue sequence. That is, the level of the
        corresponding (tuples of) standard tableaux.

        The *level* of a residue sequence is the length of its
        :meth:`multicharge`. This is the same as the level  of the
        :meth:`standard_tableaux` that belong to the residue class of tableaux
        determined by ``self``.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, (0,0,1), [0,0,1,1,2,2,3,3]).level()
            3
        """
        return len(self.multicharge())

    def size(self):
        r"""
        Return the size of the residue sequence.

        This is the size, or length, of the residue sequence, which is the
        same as  the size of the :meth:`standard_tableaux` that belong to
        the residue class of tableaux determined by ``self``.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3, (0,0,1), [0,0,1,1,2,2,3,3]).size()
            8
        """
        return len(self)


class ResidueSequences(UniqueRepresentation, Parent):
    r"""
    A parent class for :class:`ResidueSequence`.

    This class exists because :class:`ResidueSequence` needs to have a parent.
    Apart form being a parent the only useful method that it provides is
    :meth:`cell_residue`, which is a short-hand for computing the residue
    of a cell using the :meth:`ResidueSequence.quantum_characteristic`
    and :meth:`ResidueSequence.multicharge` for the residue class.

    EXAMPLES::

        sage: from sage.combinat.tableau_residues import ResidueSequences
        sage: ResidueSequences(e=0, multicharge=(0,1,2))
        0-residue sequences with multicharge (0, 1, 2)
        sage: ResidueSequences(e=0, multicharge=(0,1,2)) == ResidueSequences(e=0, multicharge=(0,1,2))
        True
        sage: ResidueSequences(e=0, multicharge=(0,1,2)) == ResidueSequences(e=3, multicharge=(0,1,2))
        False
        sage: ResidueSequences(e=0, multicharge=(0,1,2)).element_class
        <class 'sage.combinat.tableau_residues.ResidueSequences_with_category.element_class'>
    """

    Element = ResidueSequence

    def __init__(self, e, multicharge=(0,)):
        r"""
        Initialise the parent class for residue sequences.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequences
            sage: ResidueSequences(e=0, multicharge=(0,1,2))
            0-residue sequences with multicharge (0, 1, 2)
            sage: ResidueSequences(e=0, multicharge=(0,1,2)) == ResidueSequences(e=0, multicharge=(0,1,2))
            True
            sage: ResidueSequences(e=0, multicharge=(0,1,2)) == ResidueSequences(e=3, multicharge=(0,1,2))
            False

        The TestSuite fails ``_test_pickling`` because ``__getitem__`` does
        not support slices, so we skip this::

            sage: R = ResidueSequences(e=0, multicharge=(0,1,2))
            sage: TestSuite(R).run(skip='_test_elements')
        """
        self._quantum_characteristic = e
        self._base_ring = IntegerModRing(self._quantum_characteristic)
        self._multicharge = tuple(self._base_ring(i) for i in multicharge)
        super().__init__(category=Sets())

    def _repr_(self):
        r"""
        The string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequences
            sage: ResidueSequences(e=0, multicharge=(0,1,2))
            0-residue sequences with multicharge (0, 1, 2)
            sage: ResidueSequences(e=3)
            3-residue sequences with multicharge (0,)
            sage: ResidueSequences(2, (0,1,2,3))
            2-residue sequences with multicharge (0, 1, 0, 1)
        """
        return '{}-residue sequences with multicharge {}'.format(self._quantum_characteristic,
                                                                 self._multicharge)

    def an_element(self):
        r"""
        Return a particular element of ``self``.

        EXAMPLES::

            sage: TableauTuples().an_element()
            ([[1]], [[2]], [[3]], [[4]], [[5]], [[6]], [[7]])
        """
        return self.element_class(self, self._multicharge, check=True)

    def _cell_residue_level_one(self, r, c):
        r"""
        Return the residue a cell of level 1. It is called indirectly via
        :meth:`cell_residue`.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequences
            sage: ResidueSequences(3).cell_residue(1,0)  # indirect doctest
            2
        """
        return self._base_ring(c - r)

    def _cell_residue_higher_levels(self, k, r, c):
        r"""
        Return the residue a cell of level greater than 1.

        It is called indirectly via :meth:`cell_residue`.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequences
            sage: ResidueSequences(3,(0,0,1)).cell_residue(2,0,0) # indirect doctest
            1
        """
        return self._base_ring(self._multicharge[k] + c - r)

    @lazy_attribute
    def cell_residue(self, *args):
        r"""
        Return the residue a cell with respect to the quantum characteristic
        and the multicharge of the residue sequence.

        INPUT:

        - ``r``, ``c`` -- the row and column indices in level one
        - ``k``, ``r``, ``c`` -- the component, row and column indices
          in higher levels

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequences
            sage: ResidueSequences(3).cell_residue(1,1)
            0
            sage: ResidueSequences(3).cell_residue(2,1)
            2
            sage: ResidueSequences(3).cell_residue(3,1)
            1
            sage: ResidueSequences(3).cell_residue(3,2)
            2
            sage: ResidueSequences(3,(0,1,2)).cell_residue(0,0,0)
            0
            sage: ResidueSequences(3,(0,1,2)).cell_residue(0,1,0)
            2
            sage: ResidueSequences(3,(0,1,2)).cell_residue(0,1,2)
            1
            sage: ResidueSequences(3,(0,1,2)).cell_residue(1,0,0)
            1
            sage: ResidueSequences(3,(0,1,2)).cell_residue(1,1,0)
            0
            sage: ResidueSequences(3,(0,1,2)).cell_residue(1,0,1)
            2
            sage: ResidueSequences(3,(0,1,2)).cell_residue(2,0,0)
            2
            sage: ResidueSequences(3,(0,1,2)).cell_residue(2,1,0)
            1
            sage: ResidueSequences(3,(0,1,2)).cell_residue(2,0,1)
            0
        """
        # A shortcut for determining the residue of a cell, which depends on e
        # and the multicharge. The main advantage of this function is that it
        # automatically incorporates the level of this residue class. This is
        # used by the iterators for the corresponding standard tableaux classes.
        if len(self._multicharge) == 1:
            return self._cell_residue_level_one
        else:
            return self._cell_residue_higher_levels

    def check_element(self, element):
        r"""
        Check that ``element`` is a residue sequence with
        multicharge ``self.multicharge()``.

        This is weak criteria in that we only require that ``element`` is
        a tuple of elements in the underlying base ring of ``self``. Such
        a sequence is always a valid residue sequence, although there may
        be no tableaux with this residue sequence.

        EXAMPLES::

            sage: from sage.combinat.tableau_residues import ResidueSequence
            sage: ResidueSequence(3,(0,0,1),[0,0,1,1,2,2,3,3]) # indirect doctest
            3-residue sequence (0,0,1,1,2,2,0,0) with multicharge (0,0,1)
            sage: ResidueSequence(3,(0,0,1),[2,0,1,4,2,2,5,3]) # indirect doctest
            3-residue sequence (2,0,1,1,2,2,2,0) with multicharge (0,0,1)
            sage: ResidueSequence(3,(0,0,1),[2,0,1,1,2,2,3,3]) # indirect doctest
            3-residue sequence (2,0,1,1,2,2,0,0) with multicharge (0,0,1)
        """
        if any(r not in self._base_ring for r in element):
            raise ValueError('not a {}-residue sequence'.format(self._quantum_characteristic))
