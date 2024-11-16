r"""
Partition tuples

A :class:`PartitionTuple` is a tuple of partitions. That is, an ordered
`k`-tuple of partitions `\mu=(\mu^{(1)},\mu^{(2)},...,\mu^{(k)})`. If

.. MATH::

    n = \lvert \mu \rvert = \lvert \mu^{(1)} \rvert +
    \lvert \mu^{(2)} \rvert + \cdots + \lvert \mu^{(k)} \rvert

then we say that `\mu` is a `k`-partition of `n`.

In representation theory partition tuples arise as the natural indexing
set for the ordinary irreducible representations of:

- the wreath products of cyclic groups with symmetric groups,
- the Ariki-Koike algebras, or the cyclotomic Hecke algebras of
  the complex reflection groups of type `G(r,1,n)`,
- the degenerate cyclotomic Hecke algebras of type `G(r,1,n)`.

When these algebras are not semisimple, partition tuples index an important
class of modules for the algebras, which are generalisations of the Specht
modules of the symmetric groups.

Tuples of partitions also index the standard basis of the higher level
combinatorial Fock spaces. As a consequence, the combinatorics of partition
tuples encapsulates the canonical bases of crystal graphs for the irreducible
integrable highest weight modules of the (quantized) affine special linear
groups and the (quantized) affine general linear groups. By the
categorification theorems of Ariki, Varagnolo-Vasserot, Stroppel-Webster and
others, in characteristic zero the degenerate and non-degenerate cyclotomic
Hecke algebras, via their Khovanov-Lauda-Rouquier grading, categorify the
canonical bases of the quantum affine special and general linear groups.

Partitions are naturally in bijection with 1-tuples of partitions. Most of the
combinatorial operations defined on partitions extend to partition tuples in
a meaningful way. For example, the semisimple branching rules for the Specht
modules are described by adding and removing cells from partition tuples and
the modular branching rules correspond to adding and removing good and
cogood nodes, which is the underlying combinatorics for the associated
crystal graphs.

A :class:`PartitionTuple` belongs to :class:`PartitionTuples` and its derived
classes. :class:`PartitionTuples` is the parent class for all partitions
tuples. Four different classes of tuples of partitions are currently supported:

- ``PartitionTuples(level=k,size=n)`` are `k`-tuple of partitions of `n`.
- ``PartitionTuples(level=k)`` are `k`-tuple of partitions.
- ``PartitionTuples(size=n)`` are tuples of partitions of `n`.
- ``PartitionTuples()`` are tuples of partitions.

.. NOTE::

    As with :class:`Partitions`, in sage the cells, or nodes, of partition
    tuples are 0-based. For example, the (lexicographically) first cell in
    any non-empty partition tuple is `[0,0,0]`.

    EXAMPLES::

        sage: PartitionTuple([[2,2],[1,1],[2]]).cells()
        [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 1, 0), (2, 0, 0), (2, 0, 1)]

.. NOTE::

    Many :class:`PartitionTuple` methods take the individual coordinates `(k,r,c)`
    as their arguments, here `k` is the component, `r` is the row index and `c` is
    the column index.  If your coordinates are in the form ``(k,r,c)`` then use
    Python's \*-operator.

    EXAMPLES::

        sage: mu=PartitionTuple([[1,1],[2],[2,1]])
        sage: [ mu.arm_length(*c) for c in mu.cells()]
        [0, 0, 1, 0, 1, 0, 0]

.. WARNING::

    In sage, if ``mu`` is a partition tuple then ``mu[k]`` most naturally refers
    to the `k`-th component of ``mu``, so we use the convention of the
    `(k,r,c)`-th cell in a partition tuple refers to the cell in component `k`,
    row `r`, and column `c`. In the literature, the cells of a partition tuple
    are usually written in the form `(r,c,k)`, where `r` is the row index, `c`
    is the column index, and `k` is the component index.

REFERENCES:

- [DJM1998]_
- [BK2009]_

AUTHORS:

- Andrew Mathas (2012-06-01): Initial classes.

EXAMPLES:

First is a finite enumerated set and the remaining classes are infinite
enumerated sets::

    sage: PartitionTuples().an_element()
    ([1, 1, 1, 1], [2, 1, 1], [3, 1], [4])
    sage: PartitionTuples(4).an_element()
    ([], [1], [2], [3])
    sage: PartitionTuples(size=5).an_element()
    ([1], [1], [1], [1], [1])
    sage: PartitionTuples(4,5).an_element()
    ([1], [], [], [4])
    sage: PartitionTuples(3,2)[:]                                                       # needs sage.libs.flint
    [([2], [], []),
     ([1, 1], [], []),
     ([1], [1], []),
     ([1], [], [1]),
     ([], [2], []),
     ([], [1, 1], []),
     ([], [1], [1]),
     ([], [], [2]),
     ([], [], [1, 1])]
    sage: PartitionTuples(2,3).list()                                                   # needs sage.libs.flint
    [([3], []),
     ([2, 1], []),
     ([1, 1, 1], []),
     ([2], [1]),
     ([1, 1], [1]),
     ([1], [2]),
     ([1], [1, 1]),
     ([], [3]),
     ([], [2, 1]),
     ([], [1, 1, 1])]

One tuples of partitions are naturally in bijection with partitions and, as far
as possible, partition tuples attempts to identify one tuples with partitions::

    sage: Partition([4,3]) == PartitionTuple([[4,3]])
    True
    sage: Partition([4,3]) == PartitionTuple([4,3])
    True
    sage: PartitionTuple([4,3])
    [4, 3]
    sage: Partition([4,3]) in PartitionTuples()
    True

Partition tuples come equipped with many of the corresponding methods for
partitions. For example, it is possible to add and remove cells, to conjugate
partition tuples, to work with their diagrams, compare partition tuples in
dominance and so::

    sage: PartitionTuple([[4,1],[],[2,2,1],[3]]).pp()
       ****   -   **   ***
       *          **
                  *
    sage: PartitionTuple([[4,1],[],[2,2,1],[3]]).conjugate()
    ([1, 1, 1], [3, 2], [], [2, 1, 1, 1])
    sage: PartitionTuple([[4,1],[],[2,2,1],[3]]).conjugate().pp()
       *   ***   -   **
       *   **        *
       *             *
                     *
    sage: lam = PartitionTuples(3)([[3,2],[],[1,1,1,1]]); lam
    ([3, 2], [], [1, 1, 1, 1])
    sage: lam.level()
    3
    sage: lam.size()
    9
    sage: lam.category()
    Category of elements of Partition tuples of level 3
    sage: lam.parent()
    Partition tuples of level 3
    sage: lam[0]
    [3, 2]
    sage: lam[1]
    []
    sage: lam[2]
    [1, 1, 1, 1]
    sage: lam.pp()
       ***   -   *
       **        *
                 *
                 *
    sage: lam.removable_cells()
    [(0, 0, 2), (0, 1, 1), (2, 3, 0)]
    sage: lam.down_list()
    [([2, 2], [], [1, 1, 1, 1]),
     ([3, 1], [], [1, 1, 1, 1]),
     ([3, 2], [], [1, 1, 1])]
    sage: lam.addable_cells()
    [(0, 0, 3), (0, 1, 2), (0, 2, 0), (1, 0, 0), (2, 0, 1), (2, 4, 0)]
    sage: lam.up_list()
    [([4, 2], [], [1, 1, 1, 1]),
     ([3, 3], [], [1, 1, 1, 1]),
     ([3, 2, 1], [], [1, 1, 1, 1]),
     ([3, 2], [1], [1, 1, 1, 1]),
     ([3, 2], [], [2, 1, 1, 1]),
     ([3, 2], [], [1, 1, 1, 1, 1])]
    sage: lam.conjugate()
    ([4], [], [2, 2, 1])
    sage: lam.dominates( PartitionTuple([[3],[1],[2,2,1]]) )
    False
    sage: lam.dominates( PartitionTuple([[3],[2],[1,1,1]]))
    True

Every partition tuple behaves every much like a tuple of partitions::

    sage: mu = PartitionTuple([[4,1],[],[2,2,1],[3]])
    sage: [ nu for nu in mu ]
    [[4, 1], [], [2, 2, 1], [3]]
    sage: Set([ type(nu) for nu in mu ])
    {<class 'sage.combinat.partition.Partitions_all_with_category.element_class'>}
    sage: mu[2][2]
    1
    sage: mu[3]
    [3]
    sage: mu.components()
    [[4, 1], [], [2, 2, 1], [3]]
    sage: mu.components() == [ nu for nu in mu ]
    True
    sage: mu[0]
    [4, 1]
    sage: mu[1]
    []
    sage: mu[2]
    [2, 2, 1]
    sage: mu[2][0]
    2
    sage: mu[2][1]
    2
    sage: mu.level()
    4
    sage: len(mu)
    4
    sage: mu.cells()
    [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 1, 0), (2, 0, 0), (2, 0, 1),
     (2, 1, 0), (2, 1, 1), (2, 2, 0), (3, 0, 0), (3, 0, 1), (3, 0, 2)]
    sage: mu.addable_cells()
    [(0, 0, 4), (0, 1, 1), (0, 2, 0), (1, 0, 0), (2, 0, 2), (2, 2, 1),
     (2, 3, 0), (3, 0, 3), (3, 1, 0)]
    sage: mu.removable_cells()
    [(0, 0, 3), (0, 1, 0), (2, 1, 1), (2, 2, 0), (3, 0, 2)]

Attached to a partition tuple is the corresponding Young, or parabolic,
subgroup::

    sage: mu.young_subgroup()                                                           # needs sage.groups
    Permutation Group with generators
     [(), (12,13), (11,12), (8,9), (6,7), (3,4), (2,3), (1,2)]
    sage: mu.young_subgroup_generators()                                                # needs sage.groups
    [1, 2, 3, 6, 8, 11, 12]
"""

# ****************************************************************************
#       Copyright (C) 2012 Andrew Mathas <andrew.mathas@sydney.edu.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools

from .combinat import CombinatorialElement
from .integer_vector import IntegerVectors
from .partition import (Partition, Partitions, Partitions_n, _Partitions,
                        RegularPartitions_all, RegularPartitions_n)
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer_ring import ZZ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.rings.integer import Integer
from sage.sets.positive_integers import PositiveIntegers
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.groups.perm_gps.permgroup', 'PermutationGroup')
lazy_import('sage.libs.pari.all', 'pari')

# -------------------------------------------------
# Partition tuple - element class
# -------------------------------------------------


class PartitionTuple(CombinatorialElement):
    r"""
    A tuple of :class:`Partition`.

    A tuple of partition comes equipped with many of methods available to
    partitions. The ``level`` of the PartitionTuple is the length of the tuple.

    This is an ordered `k`-tuple of partitions
    `\mu=(\mu^{(1)},\mu^{(2)},...,\mu^{(k)})`. If

    .. MATH::

        n = \lvert \mu \rvert = \lvert \mu^{(1)} \rvert +
        \lvert \mu^{(2)} \rvert + \cdots + \lvert \mu^{(k)} \rvert

    then `\mu` is a `k`-partition of `n`.

    In representation theory PartitionTuples arise as the natural indexing
    set for the ordinary irreducible representations of:

    - the wreath products of cyclic groups with symmetric groups
    - the Ariki-Koike algebras, or the cyclotomic Hecke algebras of
      the complex reflection groups of type `G(r,1,n)`
    - the degenerate cyclotomic Hecke algebras of type `G(r,1,n)`

    When these algebras are not semisimple, partition tuples index an important
    class of modules for the algebras which are generalisations of the Specht
    modules of the symmetric groups.

    Tuples of partitions also index the standard basis of the higher level
    combinatorial Fock spaces. As a consequence, the combinatorics of partition
    tuples encapsulates the canonical bases of crystal graphs for the irreducible
    integrable highest weight modules of the (quantized) affine special linear
    groups and the (quantized) affine general linear groups. By the
    categorification theorems of Ariki, Varagnolo-Vasserot, Stroppel-Webster and
    others, in characteristic zero the degenerate and non-degenerate cyclotomic
    Hecke algebras, via their Khovanov-Lauda-Rouquier grading, categorify the
    canonical bases of the quantum affine special and general linear groups.

    Partitions are naturally in bijection with 1-tuples of partitions. Most of the
    combinatorial operations defined on partitions extend to PartitionTuples in
    a meaningful way. For example, the semisimple branching rules for the Specht
    modules are described by adding and removing cells from partition tuples and
    the modular branching rules correspond to adding and removing good and
    cogood nodes, which is the underlying combinatorics for the associated
    crystal graphs.

    .. WARNING::

        In the literature, the cells of a partition tuple are usually written
        in the form `(r,c,k)`, where `r` is the row index, `c` is the column
        index, and `k` is the component index. In sage, if ``mu`` is a
        partition tuple then ``mu[k]`` most naturally refers to the `k`-th
        component of ``mu``, so we use the convention of the `(k,r,c)`-th cell
        in a partition tuple refers to the cell in component `k`, row `r`, and
        column `c`.

    INPUT:

    Anything which can reasonably be interpreted as a tuple of partitions.
    That is, a list or tuple of partitions or valid input to
    :class:`Partition`.

    EXAMPLES::

        sage: mu=PartitionTuple( [[3,2],[2,1],[],[1,1,1,1]] ); mu
        ([3, 2], [2, 1], [], [1, 1, 1, 1])
        sage: nu=PartitionTuple( ([3,2],[2,1],[],[1,1,1,1]) ); nu
        ([3, 2], [2, 1], [], [1, 1, 1, 1])
        sage: mu == nu
        True
        sage: mu is nu
        False
        sage: mu in PartitionTuples()
        True
        sage: mu.parent()
        Partition tuples

        sage: lam=PartitionTuples(3)([[3,2],[],[1,1,1,1]]); lam
        ([3, 2], [], [1, 1, 1, 1])
        sage: lam.level()
        3
        sage: lam.size()
        9
        sage: lam.category()
        Category of elements of Partition tuples of level 3
        sage: lam.parent()
        Partition tuples of level 3
        sage: lam[0]
        [3, 2]
        sage: lam[1]
        []
        sage: lam[2]
        [1, 1, 1, 1]
        sage: lam.pp()
           ***   -   *
           **        *
                     *
                     *
        sage: lam.removable_cells()
        [(0, 0, 2), (0, 1, 1), (2, 3, 0)]
        sage: lam.down_list()
        [([2, 2], [], [1, 1, 1, 1]), ([3, 1], [], [1, 1, 1, 1]), ([3, 2], [], [1, 1, 1])]
        sage: lam.addable_cells()
        [(0, 0, 3), (0, 1, 2), (0, 2, 0), (1, 0, 0), (2, 0, 1), (2, 4, 0)]
        sage: lam.up_list()
        [([4, 2], [], [1, 1, 1, 1]), ([3, 3], [], [1, 1, 1, 1]), ([3, 2, 1], [], [1, 1, 1, 1]), ([3, 2], [1], [1, 1, 1, 1]), ([3, 2], [], [2, 1, 1, 1]), ([3, 2], [], [1, 1, 1, 1, 1])]
        sage: lam.conjugate()
        ([4], [], [2, 2, 1])
        sage: lam.dominates( PartitionTuple([[3],[1],[2,2,1]]) )
        False
        sage: lam.dominates( PartitionTuple([[3],[2],[1,1,1]]))
        True

    TESTS::

        sage: TestSuite( PartitionTuple([4,3,2]) ).run()
        sage: TestSuite( PartitionTuple([[4,3,2],[],[],[3,2,1]]) ).run()

    .. SEEALSO::

        - :class:`PartitionTuples`
        - :class:`Partitions`
    """
    Element = Partition

    @staticmethod
    def __classcall_private__(self, mu):
        """
        This delegates the construction of a :class:`PartitionTuple` to the
        ``element_class()`` call of the appropriate
        :class:`PartitionTuples_level`.

        TESTS::

            sage: mu=PartitionTuple([[1,1],[1]])
            sage: mu.category()
            Category of elements of Partition tuples
            sage: type(mu)
            <class 'sage.combinat.partition_tuple.PartitionTuples_all_with_category.element_class'>
        """
        if isinstance(mu, (Partition, PartitionTuple)):
            return mu

        # one way or another these two cases need to be treated separately
        if mu == [] or mu == [[]]:
            return _Partitions([])

        # We must check mu is a list of partitions
        try:
            mu = [_Partitions(mu)]
        except ValueError:
            try:
                mu = [_Partitions(nu) for nu in mu]
            except ValueError:
                raise ValueError('%s is not a tuple of Partitions' % mu)

        if len(mu) == 1:
            return _Partitions(mu[0])
        else:
            return PartitionTuples_all().element_class(PartitionTuples_all(), mu)

    def __init__(self, parent, mu):
        """
        Initialize ``self`` and checks that the input determines a tuple of
        partitions.

        EXAMPLES::

            sage: PartitionTuple([])
            []
            sage: P = PartitionTuple([[2,1,1,0],[2,1]]); P
            ([2, 1, 1], [2, 1])
            sage: TestSuite(P).run()
            sage: PartitionTuple([[],[],[2,1,2,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[], [], [2, 1, 2, 1]] is not a tuple of Partitions
        """
        mu = [_Partitions(nu) for nu in mu]
        CombinatorialElement.__init__(self, parent, mu)

    def level(self):
        """
        Return the level of this partition tuple.

        The level is the length of the tuple.

        EXAMPLES::

            sage: PartitionTuple([[2,1,1,0],[2,1]]).level()
            2
            sage: PartitionTuple([[],[],[2,1,1]]).level()
            3
        """
        return len(self._list)

    def __len__(self):
        """
        Return the length of this partition tuple.

        The length is also known as the level.

        EXAMPLES::

            sage: len( PartitionTuple([[2,1],[3,2],[1,1,1]]) )
            3
        """
        return self.level()

    def _repr_(self, compact=None):
        """
        Return a string representation of ``self`` depending on
        :meth:`PartitionTuples.options`.

        EXAMPLES::

            sage: mu=PartitionTuple(([2,1],[3,2],[1,1,1]))      # indirect doctest

            sage: PartitionTuples.options(display='list'); mu
            ([2, 1], [3, 2], [1, 1, 1])
            sage: PartitionTuples.options(display='diagram'); mu
            **   ***   *
            *    **    *
            *
            sage: PartitionTuples.options(display='compact_low'); mu
            1,2|2,3|1^3
            sage: PartitionTuples.options(display='compact_high'); mu
            2,1|3,2|1^3
            sage: PartitionTuples.options(display='exp_low'); mu
            1, 2 | 2, 3 | 1^3
            sage: PartitionTuples.options(display='exp_high'); mu
            2, 1 | 3, 2 | 1^3
            sage: PartitionTuples.options._reset()

            sage: Partitions.options(convention='French')
            sage: PartitionTuples.options(display='diagram'); mu
            *
            *    **    *
            **   ***   *
            sage: PartitionTuples.options(display='list'); mu
            ([2, 1], [3, 2], [1, 1, 1])
            sage: PartitionTuples.options(display='compact_low'); mu
            1,2|2,3|1^3
            sage: PartitionTuples.options(display='compact_high'); mu
            2,1|3,2|1^3
            sage: PartitionTuples.options(display='exp_low'); mu
            1, 2 | 2, 3 | 1^3
            sage: PartitionTuples.options(display='exp_high'); mu
            2, 1 | 3, 2 | 1^3
            sage: PartitionTuples.options._reset()
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as a Ferrers diagram.

        EXAMPLES::

            sage: print(PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_diagram())
               **   ***   *
               *    **    *
                          *
        """
        return self.diagram()

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list.

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_list()
            '([2, 1], [3, 2], [1, 1, 1])'
        """
        return '(' + ', '.join(nu._repr_() for nu in self) + ')'

    def _repr_exp_low(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_exp_low()
            '1, 2 | 2, 3 | 1^3'
            sage: PartitionTuple(([],[3,2],[1,1,1]))._repr_exp_low()
            '- | 2, 3 | 1^3'
        """
        return ' | '.join(nu._repr_exp_low() for nu in self)

    def _repr_exp_high(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1,1,1,1,1,1,1,1]))._repr_exp_high()
            '2, 1 | 3, 2 | 1^10'
            sage: PartitionTuple(([],[3,2],[1,1,1]))._repr_exp_high()
            '- | 3, 2 | 1^3'
        """
        return ' | '.join(nu._repr_exp_high() for nu in self)

    def _repr_compact_low(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_compact_low()
            '1,2|2,3|1^3'
            sage: PartitionTuple(([],[3,2],[1,1,1]))._repr_compact_low()
            '-|2,3|1^3'
        """
        return '%s' % '|'.join(mu._repr_compact_low() for mu in self)

    def _repr_compact_high(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_compact_high()
            '2,1|3,2|1^3'
            sage: PartitionTuple(([],[3,2],[1,1,1]))._repr_compact_high()
            '-|3,2|1^3'
        """
        return '|'.join(mu._repr_compact_high() for mu in self)

    # override default string representation which is str(self._list)
    __str__ = lambda self: self._repr_()  # type:ignore

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        For more on the latex options, see :meth:`Partitions.option`.

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1]])
            sage: PartitionTuples.options(latex='diagram'); latex(mu)       # indirect doctest
            {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\\
            \lr{\ast}&\lr{\ast}\\
            \lr{\ast}\\
            \end{array}$},\raisebox{-.6ex}{$\begin{array}[b]{*{1}c}\\
            \lr{\ast}\\
            \lr{\ast}\\
            \lr{\ast}\\
            \end{array}$}
            }
            sage: PartitionTuples.options(latex='exp_high'); latex(mu)      # indirect doctest
            (2,1|1^{3})
            sage: PartitionTuples.options(latex='exp_low'); latex(mu)       # indirect doctest
            (1,2|1^{3})
            sage: PartitionTuples.options(latex='list'); latex(mu)          # indirect doctest
            [[2, 1], [1, 1, 1]]
            sage: PartitionTuples.options(latex='young_diagram'); latex(mu) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$},\raisebox{-.6ex}{$\begin{array}[b]{*{1}c}\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }

            sage: PartitionTuples.options(latex='young_diagram', convention='french')
            sage: PartitionTuples.options(latex='exp_high'); latex(mu)      # indirect doctest
            (2,1|1^{3})
            sage: PartitionTuples.options(latex='exp_low'); latex(mu)       # indirect doctest
            (1,2|1^{3})
            sage: PartitionTuples.options(latex='list'); latex(mu)          # indirect doctest
            [[2, 1], [1, 1, 1]]
            sage: PartitionTuples.options(latex='young_diagram'); latex(mu) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[t]{*{2}c}\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \end{array}$},\raisebox{-.6ex}{$\begin{array}[t]{*{1}c}\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }

            sage: PartitionTuples.options._reset()
        """
        return self.parent().options._dispatch(self, '_latex_', 'latex')

    def _latex_young_diagram(self):
        """
        LaTeX output as a Young diagram.

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1]])._latex_young_diagram()
        """
        from sage.combinat.output import tex_from_array_tuple
        return tex_from_array_tuple([[["\\phantom{x}"] * row for row in mu]
                                     for mu in self._list])

    def _latex_diagram(self):
        """
        LaTeX output as a Ferrers' diagram.

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1]])._latex_diagram()
        """
        entry = self.parent().options("latex_diagram_str")
        from sage.combinat.output import tex_from_array_tuple
        return tex_from_array_tuple([[[entry] * row for row in mu]
                                     for mu in self._list], with_lines=False)

    def _latex_list(self):
        """
        LaTeX output as a list.

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1]])._latex_list()
        """
        return repr(self._list)

    def _latex_exp_low(self):
        """
        LaTeX output in exponential notation (lowest first).

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1,1,1,1,1,1,1,1]])._latex_exp_low()
        """
        txt = '|'.join(','.join('%s%s' % (a + 1, '' if e == 1 else '^{%s}' % e)
                                for a, e in enumerate(mu))
                       for mu in self.to_exp())
        return '(' + txt + ')'

    def _latex_exp_high(self):
        """
        LaTeX output in exponential notation (highest first).

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1,1,1,1,1,1,1,1]])._latex_exp_high()
        """
        txt = '|'.join(','.join(['%s%s' % (a + 1, '' if e == 1 else '^{%s}' % e)
                                 for a, e in enumerate(mu)][::-1])
                       for mu in self.to_exp())
        return '(' + txt + ')'

    def components(self):
        r"""
        Return a list containing the shape of this partition.

        This function exists in order to give a uniform way of iterating over
        the \"components\" of partition tuples of level 1 (partitions) and for
        higher levels.

        EXAMPLES::

            sage: for t in PartitionTuple([[2,1],[3,2],[3]]).components():
            ....:     print('%s\n' % t.ferrers_diagram())
            **
            *
            <BLANKLINE>
            ***
            **
            <BLANKLINE>
            ***
            <BLANKLINE>
            sage: for t in PartitionTuple([3,2]).components():
            ....:     print('%s\n' % t.ferrers_diagram())
            ***
            **
        """
        return list(self)

    def diagram(self):
        r"""
        Return a string for the Ferrers diagram of ``self``.

        EXAMPLES::

            sage: print(PartitionTuple([[2,1],[3,2],[1,1,1]]).diagram())
               **   ***   *
               *    **    *
                          *
            sage: print(PartitionTuple([[3,2],[2,1],[],[1,1,1,1]]).diagram())
               ***   **   -   *
               **    *        *
                              *
                              *
            sage: PartitionTuples.options(convention='french')
            sage: print(PartitionTuple([[3,2],[2,1],[],[1,1,1,1]]).diagram())
                              *
                              *
               **    *        *
               ***   **   -   *
            sage: PartitionTuples.options._reset()
        """
        col_len = [mu and mu[0] or 1 for mu in self]  # columns per component
        row_max = max(len(mu) for mu in self)                # maximum row length
        # There should be a fancier list compression for this but I couldn't get
        # one to work in the cases where a component was the empty partition
        diag = []
        diag_str = PartitionTuples.options('diagram_str')
        for row in range(row_max):
            line = ''
            for c in range(len(self)):
                if row == 0 and self[c] == []:
                    line += '   -'
                elif row < len(self[c]):
                    line += '   {:{}}'.format(diag_str*self[c][row],col_len[c])
                else:
                    line += '   {:{}}'.format('',col_len[c])
            diag.append(line.rstrip())
        if PartitionTuples.options('convention') == "English":
            return '\n'.join(map(str, diag))
        else:
            return '\n'.join(map(str, diag[::-1]))

    ferrers_diagram = diagram

    def pp(self):
        r"""
        Pretty print this partition tuple. See :meth:`diagram`.

        EXAMPLES::

            sage: PartitionTuple([[5,5,2,1],[3,2]]).pp()
            *****   ***
            *****   **
            **
            *
        """
        print(self.diagram())

    def size(self):
        """
        Return the size of a partition tuple.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[],[2,2]]).size()
            7
            sage: PartitionTuple([[],[],[1],[3,2,1]]).size()
            7
        """
        return sum(mu.size() for mu in self)

    def row_standard_tableaux(self):
        """
        Return the :class:`row standard tableau tuples
        <sage.combinat.tableau_tuple.RowStandardTableauTuples>`
        of shape ``self``.

        EXAMPLES::

            sage: PartitionTuple([[],[3,2,2,1],[2,2,1],[3]]).row_standard_tableaux()
            Row standard tableau tuples of shape ([], [3, 2, 2, 1], [2, 2, 1], [3])
        """
        from .tableau_tuple import RowStandardTableauTuples
        return RowStandardTableauTuples(shape=self)

    def standard_tableaux(self):
        """
        Return the :class:`standard tableau tuples<StandardTableauTuples>`
        of shape ``self``.

        EXAMPLES::

            sage: PartitionTuple([[],[3,2,2,1],[2,2,1],[3]]).standard_tableaux()
            Standard tableau tuples of shape ([], [3, 2, 2, 1], [2, 2, 1], [3])
        """
        from .tableau_tuple import StandardTableauTuples
        return StandardTableauTuples(shape=self)

    def up(self):
        r"""
        Generator (iterator) for the partition tuples that are obtained from
        ``self`` by adding a cell.

        EXAMPLES::

            sage: [mu for mu in PartitionTuple([[],[3,1],[1,1]]).up()]
            [([1], [3, 1], [1, 1]), ([], [4, 1], [1, 1]), ([], [3, 2], [1, 1]), ([], [3, 1, 1], [1, 1]), ([], [3, 1], [2, 1]), ([], [3, 1], [1, 1, 1])]
            sage: [mu for mu in PartitionTuple([[],[],[],[]]).up()]
            [([1], [], [], []), ([], [1], [], []), ([], [], [1], []), ([], [], [], [1])]
        """
        for c in range(len(self)):
            for nu in self[c].up():
                up = list(self)
                up[c] = nu
                yield PartitionTuple(up)

    def up_list(self):
        """
        Return a list of the partition tuples that can be formed from ``self``
        by adding a cell.

        EXAMPLES::

            sage: PartitionTuple([[],[3,1],[1,1]]).up_list()
            [([1], [3, 1], [1, 1]), ([], [4, 1], [1, 1]), ([], [3, 2], [1, 1]), ([], [3, 1, 1], [1, 1]), ([], [3, 1], [2, 1]), ([], [3, 1], [1, 1, 1])]
            sage: PartitionTuple([[],[],[],[]]).up_list()
            [([1], [], [], []), ([], [1], [], []), ([], [], [1], []), ([], [], [], [1])]
        """
        return list(self.up())

    def down(self):
        r"""
        Generator (iterator) for the partition tuples that are obtained from
        ``self`` by removing a cell.

        EXAMPLES::

            sage: [mu for mu in PartitionTuple([[],[3,1],[1,1]]).down()]
            [([], [2, 1], [1, 1]), ([], [3], [1, 1]), ([], [3, 1], [1])]
            sage: [mu for mu in PartitionTuple([[],[],[]]).down()]
            []
        """
        for c in range(len(self)):
            for nu in self[c].down():
                down = list(self)
                down[c] = nu
                yield PartitionTuple(down)

    def down_list(self):
        """
        Return a list of the partition tuples that can be formed from ``self``
        by removing a cell.

        EXAMPLES::

            sage: PartitionTuple([[],[3,1],[1,1]]).down_list()
            [([], [2, 1], [1, 1]), ([], [3], [1, 1]), ([], [3, 1], [1])]
            sage: PartitionTuple([[],[],[]]).down_list()
            []
        """
        return list(self.down())

    def cells(self):
        """
        Return the coordinates of the cells of ``self``. Coordinates are given
        as (component index, row index, column index) and are 0 based.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[1],[1,1,1]]).cells()
            [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 2, 0)]
        """
        return [(c,a,b) for c in range(len(self)) for (a,b) in self[c].cells()]

    def content(self, k, r, c, multicharge):
        r"""
        Return the content of the cell.

        Let `m_k =` ``multicharge[k]``, then the content of a cell is
        `m_k + c - r`.

        If the ``multicharge`` is a list of integers then it simply offsets the
        values of the contents in each component. On the other hand, if the
        ``multicharge`` belongs to `\ZZ/e\ZZ` then the corresponding
        `e`-residue is returned (that is, the content mod `e`).

        As with the content method for partitions, the content of a cell does
        not technically depend on the partition tuple, but this method is
        included because it is often useful.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[2],[1,1,1]]).content(0,1,0, [0,0,0])
            -1
            sage: PartitionTuple([[2,1],[2],[1,1,1]]).content(0,1,0, [1,0,0])
            0
            sage: PartitionTuple([[2,1],[2],[1,1,1]]).content(2,1,0, [0,0,0])
            -1

        and now we return the 3-residue of a cell::

            sage: multicharge = [IntegerModRing(3)(c) for c in [0,0,0]]
            sage: PartitionTuple([[2,1],[2],[1,1,1]]).content(0,1,0, multicharge)
            2
        """
        return multicharge[k]-r+c

    def content_tableau(self, multicharge):
        """
        Return the tableau which has (k,r,c)th entry equal to the content
        ``multicharge[k]-r+c`` of this cell.

        As with the content function, by setting the ``multicharge``
        appropriately the tableau containing the residues is returned.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[2],[1,1,1]]).content_tableau([0,0,0])
            ([[0, 1], [-1]], [[0, 1]], [[0], [-1], [-2]])
            sage: PartitionTuple([[2,1],[2],[1,1,1]]).content_tableau([0,0,1]).pp()
                0  1     0  1     1
               -1                 0
                                 -1

        as with the content function the multicharge can be used to return the
        tableau containing the residues of the cells::

            sage: multicharge=[ IntegerModRing(3)(c) for c in [0,0,1] ]
            sage: PartitionTuple([[2,1],[2],[1,1,1]]).content_tableau(multicharge).pp()
                0  1     0  1     1
                2                 0
                                  2
        """
        from sage.combinat.tableau_tuple import TableauTuple
        return TableauTuple([[[multicharge[k] - r + c
                               for c in range(self[k][r])]
                              for r in range(len(self[k]))]
                             for k in range(len(self))])

    def conjugate(self):
        """
        Return the conjugate partition tuple of ``self``.

        The conjugate partition tuple is obtained by reversing the order of the
        components and then swapping the rows and columns in each component.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[1],[1,1,1]]).conjugate()
            ([3], [1], [2, 1])
        """
        return PartitionTuple([nu.conjugate() for nu in self[::-1]])

    def dominates(self, mu):
        r"""
        Return ``True`` if the PartitionTuple dominates or equals `\mu` and
        ``False`` otherwise.

        Given partition tuples `\mu=(\mu^{(1)},...,\mu^{(m)})` and `\nu=(\nu^{(1)},...,\nu^{(n)})`
        then `\mu` dominates `\nu` if

        .. MATH::

            \sum_{k=1}^{l-1} |\mu^{(k)}| +\sum_{r \geq 1} \mu^{(l)}_r
               \geq \sum_{k=1}^{l-1} |\nu^{(k)}| + \sum_{r \geq 1} \nu^{(l)}_r

        EXAMPLES::

            sage: mu=PartitionTuple([[1,1],[2],[2,1]])
            sage: nu=PartitionTuple([[1,1],[1,1],[2,1]])
            sage: mu.dominates(mu)
            True
            sage: mu.dominates(nu)
            True
            sage: nu.dominates(mu)
            False
            sage: tau=PartitionTuple([[],[2,1],[]])
            sage: tau.dominates([[2,1],[],[]])
            False
            sage: tau.dominates([[],[],[2,1]])
            True
        """
        try:
            mu = PartitionTuple(mu)
        except ValueError:
            raise ValueError('%s must be a PartitionTuple' % mu)

        if mu == self:
            return True
        level = 0
        ssum = 0  # sum of successive rows in self
        musum = 0  # sum of successive rows in self
        while level < self.level() and level < mu.level():
            row = 0
            while row < len(self[level]) and row < len(mu[level]):
                ssum += self[level][row]
                musum += mu[level][row]
                if musum > ssum:
                    return False
                row += 1
            if row < len(self[level]):
                ssum += sum(self[level][row:])
            elif row < len(mu[level]):
                musum += sum(mu[level][row:])
                if musum > ssum:
                    return False
            level += 1
        return True

    @cached_method
    def initial_tableau(self):
        r"""
        Return the :class:`StandardTableauTuple` which has the numbers
        `1, 2, \ldots, n`, where `n` is the :meth:`size` of ``self``,
        entered in order from left to right along the rows of each component,
        where the components are ordered from left to right.

        EXAMPLES::

            sage: PartitionTuple([ [2,1],[3,2] ]).initial_tableau()
            ([[1, 2], [3]], [[4, 5, 6], [7, 8]])
        """
        from .tableau_tuple import StandardTableauTuples
        return StandardTableauTuples(self).first()

    @cached_method
    def initial_column_tableau(self):
        r"""
        Return the initial column tableau of shape ``self``.

        The initial column tableau of shape `\lambda` is the standard tableau
        that has the numbers `1` to `n`, where `n` is the :meth:`size`
        of `\lambda`, entered in order from top to bottom, and then left
        to right, down the columns of each component, starting from the
        rightmost component and working to the left.

        EXAMPLES::

            sage: PartitionTuple([ [3,1],[3,2] ]).initial_column_tableau()
            ([[6, 8, 9], [7]], [[1, 3, 5], [2, 4]])
        """
        return self.conjugate().initial_tableau().conjugate()

    def garnir_tableau(self, *cell):
        r"""
        Return the Garnir tableau of shape ``self`` corresponding to the cell
        ``cell``.

        If ``cell`` `= (k,a,c)` then `(k,a+1,c)` must belong to the diagram of
        the :class:`PartitionTuple`. If this is not the case then we return
        ``False``.

        .. NOTE::

            The function also sets ``g._garnir_cell`` equal to ``cell``
            which is used by some other functions.

        The Garnir tableaux play an important role in integral and
        non-semisimple representation theory because they determine the
        "straightening" rules for the Specht modules over an arbitrary ring.

        The Garnir tableau are the "first" non-standard tableaux which arise
        when you act by simple transpositions. If `(k,a,c)` is a cell in the
        Young diagram of a partition, which is not at the bottom of its
        column, then the corresponding Garnir tableau has the integers
        `1, 2, \ldots, n` entered in order from left to right along the rows
        of the diagram up to the cell `(k,a,c-1)`, then along the cells
        `(k,a+1,1)` to `(k,a+1,c)`, then `(k,a,c)` until the end of row `a`
        and then continuing from left to right in the remaining positions.
        The examples below probably make this clearer!

        EXAMPLES::

            sage: PartitionTuple([[5,3],[2,2],[4,3]]).garnir_tableau((0,0,2)).pp()
                 1  2  6  7  8     9 10    13 14 15 16
                 3  4  5          11 12    17 18 19
            sage: PartitionTuple([[5,3,3],[2,2],[4,3]]).garnir_tableau((0,0,2)).pp()
                 1  2  6  7  8    12 13    16 17 18 19
                 3  4  5          14 15    20 21 22
                 9 10 11
            sage: PartitionTuple([[5,3,3],[2,2],[4,3]]).garnir_tableau((0,1,2)).pp()
                 1  2  3  4  5    12 13    16 17 18 19
                 6  7 11          14 15    20 21 22
                 8  9 10
            sage: PartitionTuple([[5,3,3],[2,2],[4,3]]).garnir_tableau((1,0,0)).pp()
                 1  2  3  4  5    13 14    16 17 18 19
                 6  7  8          12 15    20 21 22
                 9 10 11
            sage: PartitionTuple([[5,3,3],[2,2],[4,3]]).garnir_tableau((1,0,1)).pp()
                 1  2  3  4  5    12 15    16 17 18 19
                 6  7  8          13 14    20 21 22
                 9 10 11
            sage: PartitionTuple([[5,3,3],[2,2],[4,3]]).garnir_tableau((2,0,1)).pp()
                 1  2  3  4  5    12 13    16 19 20 21
                 6  7  8          14 15    17 18 22
                 9 10 11
            sage: PartitionTuple([[5,3,3],[2,2],[4,3]]).garnir_tableau((2,1,1)).pp()
            Traceback (most recent call last):
            ...
            ValueError: (comp, row+1, col) must be inside the diagram

        .. SEEALSO::

            - :meth:`top_garnir_tableau`
        """
        try:
            (comp, row,col) = cell
        except ValueError:
            (comp, row,col) = cell[0]

        if comp >= len(self) or row+1 >= len(self[comp]) or col >= self[comp][row+1]:
            raise ValueError('(comp, row+1, col) must be inside the diagram')
        g = self.initial_tableau().to_list()
        a = g[comp][row][col]
        g[comp][row][col:] = list(range(a+col+1, g[comp][row+1][col]+1))
        g[comp][row+1][:col+1] = list(range(a, a+col+1))
        from .tableau_tuple import TableauTuple
        g = TableauTuple(g)
        g._garnir_cell = (comp,row,col)
        return g

    def top_garnir_tableau(self, e, cell):
        r"""
        Return the most dominant *standard* tableau which dominates the
        corresponding Garnir tableau and has the same residue that has shape
        ``self`` and is determined by ``e`` and ``cell``.

        The Garnir tableau play an important role in integral and
        non-semisimple representation theory because they determine the
        "straightening" rules for the Specht modules over an arbitrary ring.
        The *top Garnir tableaux* arise in the graded representation theory of
        the symmetric groups and higher level Hecke algebras. They were
        introduced in [KMR2012]_.

        If the Garnir node is ``cell=(k,r,c)`` and `m` and `M` are the entries
        in the cells ``(k,r,c)`` and ``(k,r+1,c)``, respectively, in the
        initial tableau then the top ``e``-Garnir tableau is obtained by
        inserting the numbers `m, m+1, \ldots, M` in order from left to right
        first in the cells in row ``r+1`` which are not in the ``e``-Garnir
        belt, then in the cell in rows ``r`` and ``r+1`` which are in the
        Garnir belt and then, finally, in the remaining cells in row ``r``
        which are not in the Garnir belt. All other entries in the tableau
        remain unchanged.

        If ``e = 0``, or if there are no ``e``-bricks in either row ``r`` or
        ``r+1``, then the top Garnir tableau is the corresponding Garnir
        tableau.

        EXAMPLES::

            sage: PartitionTuple([[3,3,2],[5,4,3,2]]).top_garnir_tableau(2,(1,0,2)).pp()
                1  2  3     9 10 12 13 16
                4  5  6    11 14 15 17
                7  8       18 19 20
                           21 22
            sage: PartitionTuple([[3,3,2],[5,4,3,2]]).top_garnir_tableau(2,(1,0,1)).pp()
                1  2  3     9 10 11 12 13
                4  5  6    14 15 16 17
                7  8       18 19 20
                           21 22
            sage: PartitionTuple([[3,3,2],[5,4,3,2]]).top_garnir_tableau(3,(1,0,1)).pp()
                1  2  3     9 12 13 14 15
                4  5  6    10 11 16 17
                7  8       18 19 20
                           21 22

            sage: PartitionTuple([[3,3,2],[5,4,3,2]]).top_garnir_tableau(3,(3,0,1)).pp()
            Traceback (most recent call last):
            ...
            ValueError: (comp, row+1, col) must be inside the diagram

        .. SEEALSO::

            - :meth:`~sage.combinat.partition.Partition_tuple.garnir_tableau`
        """
        (comp,row,col) = cell
        if comp >= len(self) or row+1 >= len(self[comp]) or col >= self[comp][row+1]:
            raise ValueError('(comp, row+1, col) must be inside the diagram')

        g = self.garnir_tableau(cell)

        if e == 0:
            return      # no more dominant tableau of the same residue

        a = e*int((self[comp][row]-col)/e)    # number of cells in the e-bricks in row `row`
        b = e*int((col+1)/e)            # number of cells in the e-bricks in row `row+1`

        if a == 0 or b == 0:
            return self.garnir_tableau(cell)

        t = g.to_list()
        m = t[comp][row+1][0]            # smallest number of 0-Garnir belt
        # now we will put the number m,m+1,...,t[row+1][col] in order into t
        t[comp][row][col:a+col] = [m+col-b+1+i for i in range(a)]
        t[comp][row+1][col-b+1:col+1] = [m+a+col-b+1+i for i in range(b)]
        from .tableau_tuple import StandardTableauTuple
        return StandardTableauTuple(t)

    def arm_length(self, k, r, c):
        """
        Return the length of the arm of cell ``(k, r, c)`` in ``self``.

        INPUT:

        - ``k`` -- the component
        - ``r`` -- the row
        - ``c`` -- the cell

        OUTPUT: the arm length as an integer

        The arm of cell ``(k, r, c)`` is the number of cells in the ``k``-th
        component which are to the right of the cell in row ``r`` and column
        ``c``.

        EXAMPLES::

            sage: PartitionTuple([[],[2,1],[2,2,1],[3]]).arm_length(2,0,0)
            1
            sage: PartitionTuple([[],[2,1],[2,2,1],[3]]).arm_length(2,0,1)
            0
            sage: PartitionTuple([[],[2,1],[2,2,1],[3]]).arm_length(2,2,0)
            0
        """
        try:
            return self[k][r]-(c+1)
        except IndexError:
            raise ValueError("The cell %s is not in the diagram" % ((k,r,c),))

    def leg_length(self, k, r, c):
        """
        Return the length of the leg of cell ``(k, r, c)`` in ``self``.

        INPUT:

        - ``k`` -- the component
        - ``r`` -- the row
        - ``c`` -- the cell

        OUTPUT: the leg length as an integer

        The leg of cell ``(k, r, c)`` is the number of cells in the ``k``-th
        component which are below the node in row ``r`` and column ``c``.

        EXAMPLES::

            sage: PartitionTuple([[],[2,1],[2,2,1],[3]]).leg_length(2,0,0)
            2
            sage: PartitionTuple([[],[2,1],[2,2,1],[3]]).leg_length(2,0,1)
            1
            sage: PartitionTuple([[],[2,1],[2,2,1],[3]]).leg_length(2,2,0)
            0
        """
        try:
            return self[k].leg_length(r,c)
        except IndexError:
            raise ValueError("The cell is not in the diagram")

    def contains(self, mu):
        r"""
        Return ``True`` if this partition tuple contains `\mu`.

        If `\lambda=(\lambda^{(1)}, \ldots, \lambda^{(l)})` and
        `\mu=(\mu^{(1)}, \ldots, \mu^{(m)})` are two partition tuples then
        `\lambda` contains `\mu` if `m \leq l` and
        `\mu^{(i)}_r \leq \lambda^{(i)}_r` for `1 \leq i \leq m` and `r \geq 0`.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[2],[2,1]]).contains( PartitionTuple([[1,1],[2],[2,1]]) )
            True
        """
        return mu.level() <= self.level() and all(self[c].contains(mu[c]) for c in range(len(mu)))

    def hook_length(self, k, r, c):
        r"""
        Return the length of the hook of cell ``(k, r, c)`` in the partition.

        The hook of cell ``(k, r, c)`` is defined as the cells to the right or
        below (in the English convention). If your coordinates are in the
        form ``(k,r,c)``, use Python's \*-operator.

        EXAMPLES::

            sage: mu=PartitionTuple([[1,1],[2],[2,1]])
            sage: [ mu.hook_length(*c) for c in mu.cells()]
            [2, 1, 2, 1, 3, 1, 1]
        """
        try:
            return self[k].hook_length(r,c)
        except IndexError:
            raise ValueError("The cell is not in the diagram")

    def to_exp(self, k=0):
        """
        Return a tuple of the multiplicities of the parts of a partition.

        Use the optional parameter ``k`` to get a return list of length at
        least ``k``.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[2],[2,1]]).to_exp()
            ([2], [0, 1], [1, 1])
            sage: PartitionTuple([[1,1],[2,2,2,2],[2,1]]).to_exp()
            ([2], [0, 4], [1, 1])
        """
        return tuple(self[c].to_exp(k) for c in range(len(self)))

    def removable_cells(self):
        """
        Return a list of the removable cells of this partition tuple.

        All indices are of the form ``(k, r, c)``, where ``r`` is the
        row-index, ``c`` is the column index and ``k`` is the component.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[2],[2,1]]).removable_cells()
            [(0, 1, 0), (1, 0, 1), (2, 0, 1), (2, 1, 0)]
            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).removable_cells()
            [(0, 1, 0), (1, 0, 3), (1, 1, 2), (2, 0, 1), (2, 2, 0)]
        """
        return [(k,r,c) for k in range(len(self)) for (r,c) in self[k].removable_cells()]

    corners = removable_cells  # for compatibility with partitions

    def addable_cells(self):
        """
        Return a list of the removable cells of this partition tuple.

        All indices are of the form ``(k, r, c)``, where ``r`` is the
        row-index, ``c`` is the column index and ``k`` is the component.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[2],[2,1]]).addable_cells()
            [(0, 0, 1), (0, 2, 0), (1, 0, 2), (1, 1, 0), (2, 0, 2), (2, 1, 1), (2, 2, 0)]
            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).addable_cells()
            [(0, 0, 1), (0, 2, 0), (1, 0, 4), (1, 1, 3), (1, 2, 0), (2, 0, 2), (2, 1, 1), (2, 3, 0)]
        """
        return [(k,r,c) for k in range(len(self)) for (r,c) in self[k].addable_cells()]

    outside_corners = addable_cells  # for compatibility with partitions

    def add_cell(self, k, r, c):
        r"""
        Return the partition tuple obtained by adding a cell in row ``r``,
        column ``c``, and component ``k``.

        This does not change ``self``.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).add_cell(0,0,1)
            ([2, 1], [4, 3], [2, 1, 1])
        """
        if (k, r, c) in self.addable_cells():  # an addable cell
            mu = self.to_list()
            if c == 0:
                mu[k].append(1)
            else:
                mu[k][r] += 1
            return PartitionTuple(mu)
        else:
            raise ValueError("%s is not an addable cell" % ((k, r, c),))

    def remove_cell(self, k, r, c):
        """
        Return the partition tuple obtained by removing a cell in row ``r``,
        column ``c``, and component ``k``.

        This does not change ``self``.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).remove_cell(0,1,0)
            ([1], [4, 3], [2, 1, 1])
        """
        if (k, r, c) in self.removable_cells():  # a removable cell
            mu = self.to_list()
            mu[k][r] -= 1
            return PartitionTuple(mu)
        else:
            raise ValueError("%s is not a removable cell" % ((k, r, c),))

    def to_list(self):
        r"""
        Return ``self`` as a list of lists.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).to_list()
            [[1, 1], [4, 3], [2, 1, 1]]

        TESTS::

            sage: all(mu==PartitionTuple(mu.to_list()) for mu in PartitionTuples(4,4))  # needs sage.libs.flint
            True
        """
        return [mu.to_list() for mu in self]

    def young_subgroup(self):
        """
        Return the corresponding Young, or parabolic, subgroup of the
        symmetric group.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[4,2],[1]]).young_subgroup()                    # needs sage.groups
            Permutation Group with generators [(), (8,9), (6,7), (5,6), (4,5), (1,2)]
        """
        gens = []
        m = 0
        for comp in self:
            for row in comp:
                gens.extend([(c, c+1) for c in range(m+1, m+row)])
                m += row
        gens.append(list(range(1, self.size()+1)))  # to ensure we get a subgroup of Sym_n
        return PermutationGroup(gens)

    def young_subgroup_generators(self):
        """
        Return an indexing set for the generators of the corresponding Young
        subgroup.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[4,2],[1]]).young_subgroup_generators()
            [1, 4, 5, 6, 8]
        """
        gens = []
        m = 0
        for comp in self:
            for row in comp:
                gens.extend(list(range(m + 1, m + row)))
                m += row
        return gens

    @cached_method
    def _initial_degree(self, e, multicharge):
        r"""
        Return the Brundan-Kleshchev-Wang degree of the initial tableau
        of shape ``self``.

        This degree depends only the shape of the tableau and it is
        used as the base case for computing the degrees of all tableau of
        shape ``self``, which is why this method is cached. See
        :meth:`sage.combinat.tableau.Tableau.degree` for more information.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[2,2]])._initial_degree(0,(0,0))
            1
            sage: PartitionTuple([[2,1],[2,2]])._initial_degree(2,(0,0))
            4
            sage: PartitionTuple([[2,1],[2,2]])._initial_degree(3,(0,0))
            1
            sage: PartitionTuple([[2,1],[2,2]])._initial_degree(4,(0,0))
            1
        """
        if e == 0:
            deg = 0
        else:
            deg = sum(mu._initial_degree(e) for mu in self)
        I = IntegerModRing(e)
        multires = [I(k) for k in multicharge]
        for (k,r,c) in self.cells():
            res = I(multicharge[k]-r+c)
            for l in range(k+1, self.level()):
                if res == multires[l]:
                    deg += 1
        return deg

    def degree(self, e):
        r"""
        Return the ``e``-th degree of ``self``.

        The `e`-th degree is the sum of the degrees of the standard
        tableaux of shape `\lambda`. The `e`-th degree is the exponent
        of `\Phi_e(q)` in the Gram determinant of the Specht module for a
        semisimple cyclotomic Hecke algebra of type `A` with parameter `q`.

        For this calculation the multicharge `(\kappa_1, \ldots, \kappa_l)`
        is chosen so that `\kappa_{r+1} - \kappa_r > n`, where `n` is
        the :meth:`size` of `\lambda` as this ensures that the Hecke algebra
        is semisimple.

        INPUT:

        - ``e`` -- integer `e > 1`

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: PartitionTuple([[2,1],[2,2]]).degree(2)
            532
            sage: PartitionTuple([[2,1],[2,2]]).degree(3)
            259
            sage: PartitionTuple([[2,1],[2,2]]).degree(4)
            196
            sage: PartitionTuple([[2,1],[2,2]]).degree(5)
            105
            sage: PartitionTuple([[2,1],[2,2]]).degree(6)
            105
            sage: PartitionTuple([[2,1],[2,2]]).degree(7)
            0

        Therefore,  the Gram determinant of `S(2,1|2,2)` when the Hecke parameter
        `q` is "generic" is

        .. MATH::

            q^N \Phi_2(q)^{532}\Phi_3(q)^{259}\Phi_4(q)^{196}\Phi_5(q)^{105}\Phi_6(q)^{105}

        for some integer `N`.  Compare with :meth:`prime_degree`.
        """
        multicharge = tuple([i*self.size() for i in range(self.size())])
        return sum(t.degree(e, multicharge) for t in self.standard_tableaux())

    def prime_degree(self, p):
        r"""
        Return the ``p``-th prime degree of ``self``.

        The degree of a partition `\lambda` is the sum of the `e`-degrees`
        of the standard tableaux of shape `\lambda` (see :meth:`degree`),
        for `e` a power of the prime `p`. The prime degree gives the
        exponent of `p` in the Gram determinant of the integral Specht
        module of the symmetric group.

        The `p`-th degree is the sum of the degrees of the standard tableaux
        of shape `\lambda`. The `p`-th degree is the exponent of `p` in the
        Gram determinant of a semisimple cyclotomic Hecke algebra of type `A`
        with parameter `q = 1`.

        As with :meth:`degree`, for this calculation the multicharge
        `(\kappa_1, \ldots, \kappa_l)` is chosen so that
        `\kappa_{r+1} - \kappa_r > n`, where `n` is the :meth:`size`
        of `\lambda` as this ensures that the Hecke algebra is semisimple.

        INPUT:

        - ``e`` -- an  integer `e > 1`
        - ``multicharge`` -- an `l`-tuple of integers, where `l` is
          the :meth:`level` of ``self``

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: PartitionTuple([[2,1],[2,2]]).prime_degree(2)
            728
            sage: PartitionTuple([[2,1],[2,2]]).prime_degree(3)
            259
            sage: PartitionTuple([[2,1],[2,2]]).prime_degree(5)
            105
            sage: PartitionTuple([[2,1],[2,2]]).prime_degree(7)
            0

        Therefore, the Gram determinant of `S(2,1|2,2)` when `q=1` is
        `2^{728} 3^{259}5^{105}`. Compare with :meth:`degree`.
        """
        ps = [p]

        while ps[-1]*p < self.size():
            ps.append(ps[-1] * p)
        multicharge = tuple([i*self.size() for i in range(self.size())])
        return sum(t.degree(pk, multicharge) for pk in ps for t in self.standard_tableaux())

    @cached_method
    def block(self, e, multicharge):
        r"""
        Return a dictionary `\beta` that determines the block associated to
        the partition ``self`` and the
        :meth:`~sage.combinat.tableau_residues.ResidueSequence.quantum_characteristic` ``e``.

        INPUT:

        - ``e`` -- the quantum characteristic

        - ``multicharge`` -- the multicharge (default: `(0,)`)

        OUTPUT:

        - a dictionary giving the multiplicities of the residues in the
          partition tuple ``self``

        In more detail, the value ``beta[i]`` is equal to the
        number of nodes of residue ``i``. This corresponds to
        the positive root

        .. MATH::

            \sum_{i\in I} \beta_i \alpha_i \in Q^+,

        a element of the positive root lattice of the corresponding
        Kac-Moody algebra. See [DJM1998]_ and [BK2009]_ for more details.

        This is a useful statistics because two Specht modules for a cyclotomic
        Hecke algebra of type `A` belong to the same block if and only if they
        correspond to same element `\beta` of the root lattice, given above.

        We return a dictionary because when the quantum characteristic is `0`,
        the Cartan type is `A_{\infty}`, in which case the simple roots are
        indexed by the integers.

        EXAMPLES::

            sage: PartitionTuple([[2,2],[2,2]]).block(0,(0,0))
            {-1: 2, 0: 4, 1: 2}
            sage: PartitionTuple([[2,2],[2,2]]).block(2,(0,0))
            {0: 4, 1: 4}
            sage: PartitionTuple([[2,2],[2,2]]).block(2,(0,1))
            {0: 4, 1: 4}
            sage: PartitionTuple([[2,2],[2,2]]).block(3,(0,2))
            {0: 3, 1: 2, 2: 3}
            sage: PartitionTuple([[2,2],[2,2]]).block(3,(0,2))
            {0: 3, 1: 2, 2: 3}
            sage: PartitionTuple([[2,2],[2,2]]).block(3,(3,2))
            {0: 3, 1: 2, 2: 3}
            sage: PartitionTuple([[2,2],[2,2]]).block(4,(0,0))
            {0: 4, 1: 2, 3: 2}
        """
        block = {}
        Ie = IntegerModRing(e)
        for (k,r,c) in self.cells():
            i = Ie(multicharge[k] + c - r)
            block[i] = block.get(i, 0) + 1
        return block

    def defect(self, e, multicharge):
        r"""
        Return the ``e``-defect or the ``e``-weight ``self``.

        The `e`-defect is the number of (connected) `e`-rim hooks
        that can be removed from the partition.

        The defect of a partition tuple is given by

        .. MATH::

            \text{defect}(\beta) = (\Lambda, \beta) - \tfrac12(\beta, \beta),

        where `\Lambda = \sum_r \Lambda_{\kappa_r}` for the multicharge
        `(\kappa_1, \ldots, \kappa_{\ell})` and
        `\beta = \sum_{(r,c)} \alpha_{(c-r) \pmod e}`, with the sum
        being over the cells in the partition.

        INPUT:

        - ``e`` -- the quantum characteristic

        - ``multicharge`` -- the multicharge (default: `(0,)`)

        OUTPUT: a nonnegative integer, which is the defect of the block
        containing the partition tuple ``self``

        EXAMPLES::

            sage: PartitionTuple([[2,2],[2,2]]).defect(0,(0,0))
            0
            sage: PartitionTuple([[2,2],[2,2]]).defect(2,(0,0))
            8
            sage: PartitionTuple([[2,2],[2,2]]).defect(2,(0,1))
            8
            sage: PartitionTuple([[2,2],[2,2]]).defect(3,(0,2))
            5
            sage: PartitionTuple([[2,2],[2,2]]).defect(3,(0,2))
            5
            sage: PartitionTuple([[2,2],[2,2]]).defect(3,(3,2))
            2
            sage: PartitionTuple([[2,2],[2,2]]).defect(4,(0,0))
            0
        """
        # Will correspond to an element of the positive root lattice
        #   corresponding to the block.
        # We use a dictionary to cover the case when e = 0.
        beta = self.block(e, multicharge)
        Ie = IntegerModRing(e)
        return (sum(beta.get(r, 0) for r in multicharge)
                - sum(beta[r]**2 - beta[r] * beta.get(Ie(r+1), 0) for r in beta))

# -------------------------------------------------
# Partition tuples - parent classes
# -------------------------------------------------


class PartitionTuples(UniqueRepresentation, Parent):
    r"""
    Class of all partition tuples.

    For more information about partition tuples, see :class:`PartitionTuple`.

    This is a factory class which returns the appropriate parent based on
    the values of ``level``, ``size``, and ``regular``

    INPUT:

    - ``level`` -- the length of the tuple

    - ``size`` -- the total number of cells

    - ``regular`` -- positive integer or a tuple of nonnegative
      integers; if an integer, the highest multiplicity an entry may
      have in a component plus `1`

    If a level `k` is specified and ``regular`` is a tuple of integers
    `\ell_1, \ldots, \ell_k`, then this specifies partition tuples `\mu`
    such that `\mu_i` is `\ell_i`-regular, where `0` here
    represents `\infty`-regular partitions (equivalently, partitions
    without restrictions). If ``regular`` is an integer `\ell`, then
    we set `\ell_i = \ell` for all `i`.

    TESTS::

        sage: [ [2,1],[],[3] ] in PartitionTuples()
        True
        sage: ( [2,1],[],[3] ) in PartitionTuples()
        True
        sage: ( [] ) in PartitionTuples()
        True
        sage: PartitionTuples(level=1, regular=(0,))
        Partitions
        sage: PartitionTuples(level=1, size=3, regular=(0,))
        Partitions of the integer 3

    Check that :issue:`14145` has been fixed::

        sage: 1 in PartitionTuples()
        False
    """

    @staticmethod
    def __classcall_private__(klass, level=None, size=None, regular=None):
        r"""
        Return the correct parent object based upon the input.

        TESTS::

            sage: PartitionTuples()
            Partition tuples
            sage: PartitionTuples(3)
            Partition tuples of level 3
            sage: PartitionTuples(size=3)
            Partition tuples of size 3
            sage: PartitionTuples(3,8)
            Partition tuples of level 3 and size 8
            sage: PartitionTuples(level=3, regular=(0,2,4))
            (0, 2, 4)-Regular partition tuples of level 3
            sage: PartitionTuples(level=1,regular=(4,)) is PartitionTuples(level=1, regular=4)
            True
        """
        # sanity testing
        if level is not None and (not isinstance(level, (int, Integer)) or level < 1):
            raise ValueError('the level must be a positive integer')

        if size is not None and (not isinstance(size, (int, Integer)) or size < 0):
            raise ValueError('the size must be a nonnegative integer')

        if isinstance(regular, (list, tuple)):
            if level is None:
                raise ValueError("When no level is specified, regular must be "
                                 "a positive integer")
            if len(regular) != level:
                raise ValueError("regular must be a list of length {}, got {}".format(
                                 level, regular))
        if regular == 0:
            raise ValueError("regular must be a positive integer or a tuple "
                             "of nonnegative integers")
        if level is None:
            if size is None:
                if regular is None:
                    return PartitionTuples_all()
                return RegularPartitionTuples_all(regular)

            if regular is None:
                return PartitionTuples_size(size)
            return RegularPartitionTuples_size(size, regular)

        elif level == 1:
            if isinstance(regular, (list, tuple)):
                regular = regular[0]
            if size is None:
                if regular is None or regular == 0:
                    return _Partitions
                return RegularPartitions_all(regular)

            if regular is None or regular == 0:
                return Partitions_n(size)
            return RegularPartitions_n(size, regular)

        # Higher level
        if regular is not None:
            if not isinstance(regular, (list, tuple)):
                regular = (regular,) * level
            else:
                regular = tuple(regular)
        if size is None:
            if regular is None:
                return PartitionTuples_level(level)
            return RegularPartitionTuples_level(level, regular)
        if regular is None:
            return PartitionTuples_level_size(level, size)
        return RegularPartitionTuples_level_size(level, size, regular)

    Element = PartitionTuple
    options = Partitions.options

    # default for level
    _level = None
    _size = None

    def _element_constructor_(self, mu):
        r"""
        Construct an element of :class:`PartitionTuple`.

        INPUT:

        - ``mu`` -- tuple of partitions

        OUTPUT: the corresponding :class:`PartitionTuple` object

        TESTS::

            sage: PartitionTuple([[2],[2],[]]).parent()
            Partition tuples
            sage: parts = PartitionTuples(3)
            sage: parts([[2,1],[1],[2,2,2]]).parent() is parts
            True
            sage: PartitionTuples._element_constructor_(PartitionTuples(), [[2,1],[3,2],[1,1,1]])
            ([2, 1], [3, 2], [1, 1, 1])
            sage: parts([[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2]] is not a Partition tuples of level 3
        """
        # one way or another these two cases need to be treated separately
        if mu == [] or mu == () or mu == [[]]:
            if mu not in self:
                raise ValueError('{} is not a {}'.format(mu, self))
            return self.element_class(self, [_Partitions([])])

        # As partitions are 1-tuples of partitions we need to treat them separately
        try:
            mu = [_Partitions(mu)]
        except ValueError:
            try:
                mu = [_Partitions(nu) for nu in mu]
            except ValueError:
                raise ValueError('{} is not a {}'.format(mu, self))

        if mu not in self:
            raise ValueError('{} is not a {}'.format(mu, self))
        return self.element_class(self, mu)

    def __contains__(self, mu):
        r"""
        Return ``True`` if `\mu` is in ``self``.

        TESTS::

            sage: PartitionTuple([[3,2],[2]]) in PartitionTuples()
            True
            sage: PartitionTuple([[3,2],[],[],[],[2]]) in PartitionTuples()
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[2]]) in PartitionTuples()
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[3]]) in PartitionTuples()
            True
            sage: all(mu in PartitionTuples() for mu in PartitionTuples(3,8))           # needs sage.libs.flint
            True
            sage: [5,1,1] in PartitionTuples()
            True
            sage: [[5,1,1]] in PartitionTuples()
            True
            sage: la = Partition([3,3,1])
            sage: PT = PartitionTuples()
            sage: la in PT
            True
            sage: PT(la)
            ([3, 3, 1])

        Check that :issue:`14145` is fixed::

            sage: 1 in PartitionTuples()
            False
        """
        if isinstance(mu, (PartitionTuple, Partition)):
            return True
        if isinstance(mu, (tuple, list)):
            if not mu:
                return True
            if mu[0] in ZZ:
                return mu in _Partitions
            return all(m in _Partitions for m in mu)
        return False

    def __getitem__(self, r):
        r"""
        The default implementation of ``__getitem__()`` for enumerated sets
        does not allow slices, so we override it.

        EXAMPLES::

            sage: PartitionTuples()[10:20]                                              # needs sage.libs.flint
            [([1, 1, 1]),
             ([2], []),
             ([1, 1], []),
             ([1], [1]),
             ([], [2]),
             ([], [1, 1]),
             ([1], [], []),
             ([], [1], []),
             ([], [], [1]),
             ([], [], [], [])]
        """
        if isinstance(r,(int,Integer)):
            return self.unrank(r)
        elif isinstance(r,slice):
            start = 0 if r.start is None else r.start
            stop = r.stop
            if stop is None and not self.is_finite():
                raise ValueError('infinite set')
        else:
            raise ValueError('r must be an integer or a slice')
        count = 0
        parts = []
        for t in self:
            if count == stop:
                break
            if count >= start:
                parts.append(t)
            count += 1

        # this is to cope with empty slice endpoints like [6:] or [:]
        if count == stop or stop is None:
            return parts
        raise IndexError('value out of range')

    def level(self):
        """
        Return the level or ``None`` if it is not defined.

        EXAMPLES::

            sage: PartitionTuples().level() is None
            True
            sage: PartitionTuples(7).level()
            7
        """
        return self._level

    def size(self):
        """
        Return the size or ``None`` if it is not defined.

        EXAMPLES::

            sage: PartitionTuples().size() is None
            True
            sage: PartitionTuples(size=7).size()
            7
        """
        return self._size

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples().an_element()
            ([1, 1, 1, 1], [2, 1, 1], [3, 1], [4])
        """
        return PartitionTuple(([1, 1, 1, 1], [2, 1, 1], [3, 1], [4]))


class PartitionTuples_all(PartitionTuples):
    """
    Class of partition tuples of a arbitrary level and arbitrary sum.
    """

    def __init__(self):
        r"""
        Initialize the class.

        EXAMPLES::

            sage: TestSuite( PartitionTuples() ).run()                                  # needs sage.libs.flint
        """
        super().__init__(category=InfiniteEnumeratedSets())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples()
            Partition tuples
        """
        return 'Partition tuples'

    def __iter__(self):
        r"""
        Iterate through the infinite class of partition tuples of arbitrary
        level and size.

        EXAMPLES::

            sage: PartitionTuples()[:20]                                                # needs sage.libs.flint
            [([]),
             ([1]),
             ([], []),
             ([2]),
             ([1, 1]),
             ([1], []),
             ([], [1]),
             ([], [], []),
             ([3]),
             ([2, 1]),
             ([1, 1, 1]),
             ([2], []),
             ([1, 1], []),
             ([1], [1]),
             ([], [2]),
             ([], [1, 1]),
             ([1], [], []),
             ([], [1], []),
             ([], [], [1]),
             ([], [], [], [])]
        """
        for size in NN:
            for level in range(size+1):
                for mu in PartitionTuples_level_size(level+1,size-level):
                    yield self._element_constructor_(mu)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples().an_element()
            ([1, 1, 1, 1], [2, 1, 1], [3, 1], [4])
        """
        return self.element_class(self,([1,1,1,1],[2,1,1],[3,1],[4]))


class PartitionTuples_level(PartitionTuples):
    """
    Class of partition tuples of a fixed level, but summing to an arbitrary
    integer.
    """

    def __init__(self, level, category=None):
        r"""
        Initialize this class.

        EXAMPLES::

            sage: PartitionTuples(4)
            Partition tuples of level 4
            sage: PartitionTuples(level=6)
            Partition tuples of level 6
            sage: TestSuite( PartitionTuples(level=4) ).run()                           # needs sage.libs.flint
        """
        if level not in NN:
            raise ValueError('level must be a nonnegative integer')
        if category is None:
            category = InfiniteEnumeratedSets()
        super().__init__(category=category)
        self._level = level

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples(2)
            Partition tuples of level 2
        """
        return 'Partition tuples of level {}'.format(self._level)

    def __contains__(self, mu):
        r"""
        Return ``True`` if `\mu` is in ``self``.

        TESTS::

            sage: PartitionTuple([[3,2],[2]]) in PartitionTuples(2)
            True
            sage: PartitionTuple([[3,2],[2]]) in PartitionTuples(level=2)
            True
            sage: PartitionTuple([[2,2,1],[2]]) in PartitionTuples(level=2)
            True
            sage: PartitionTuple([[2,2,1],[],[2]]) in PartitionTuples(level=2)
            False
            sage: all(mu in PartitionTuples(3) for mu in PartitionTuples(3,8))          # needs sage.libs.flint
            True

        Check that :issue:`14145` is fixed::

            sage: 1 in PartitionTuples(level=2)
            False
        """
        # Note that self._level > 1
        return PartitionTuples.__contains__(self, mu) and len(mu) == self._level

    def __iter__(self):
        r"""
        Iterate through the infinite class of partition tuples of fixed level.

        EXAMPLES::

            sage: parts = PartitionTuples(3)
            sage: [parts[k] for k in range(20)]                                         # needs sage.libs.flint
            [([], [], []),
             ([1], [], []),
             ([], [1], []),
             ([], [], [1]),
             ([2], [], []),
             ([1, 1], [], []),
             ([1], [1], []),
             ([1], [], [1]),
             ([], [2], []),
             ([], [1, 1], []),
             ([], [1], [1]),
             ([], [], [2]),
             ([], [], [1, 1]),
             ([3], [], []),
             ([2, 1], [], []),
             ([1, 1, 1], [], []),
             ([2], [1], []),
             ([1, 1], [1], []),
             ([2], [], [1]),
             ([1, 1], [], [1])]
        """
        for size in NN:
            for mu in PartitionTuples_level_size(self._level, size):
                yield self.element_class(self, list(mu))

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples(level=4).an_element()
            ([], [1], [2], [3])
        """
        return self.element_class(self, tuple([l] for l in range(self.level())))


class PartitionTuples_size(PartitionTuples):
    """
    Class of partition tuples of a fixed size, but arbitrary level.
    """

    def __init__(self, size):
        r"""
        Initialize this class.

        EXAMPLES::

            sage: PartitionTuples(size=4)
            Partition tuples of size 4
            sage: PartitionTuples(size=6)
            Partition tuples of size 6

            sage: TestSuite( PartitionTuples(size=6) ).run()                            # needs sage.libs.flint
        """
        if size not in NN:
            raise ValueError('size must be a nonnegative integer')
        super().__init__(category=InfiniteEnumeratedSets())
        self._size = size

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples(size=4)    # indirect doctest
            Partition tuples of size 4
        """
        return 'Partition tuples of size {}'.format(self._size)

    def __contains__(self, mu):
        r"""
        Return ``True`` if `\mu` is in ``self``.

        TESTS::

            sage: PartitionTuple([[3,2],[2]]) in PartitionTuples(size=7)
            True
            sage: PartitionTuple([[3,2],[],[],[],[2]]) in PartitionTuples(size=7)
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[2]]) in PartitionTuples(size=7)
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[3]]) in PartitionTuples(size=7)
            False
            sage: all(mu in PartitionTuples(size=8) for mu in PartitionTuples(3,8))     # needs sage.libs.flint
            True
            sage: [3, 2, 1] in PartitionTuples(size=7)
            False

        Check that :issue:`14145` is fixed::

            sage: 1 in PartitionTuples(size=7)
            False
        """
        if mu in _Partitions:
            return self._size == sum(mu)
        return PartitionTuples.__contains__(self, mu) and self._size == sum(map(sum, mu))

    def __iter__(self):
        r"""
        Iterate through the infinite class of partition tuples of a fixed size.

        EXAMPLES::

            sage: PartitionTuples(size=3)[:20]                                          # needs sage.libs.flint
            [([3]),
             ([2, 1]),
             ([1, 1, 1]),
             ([3], []),
             ([2, 1], []),
             ([1, 1, 1], []),
             ([2], [1]),
             ([1, 1], [1]),
             ([1], [2]),
             ([1], [1, 1]),
             ([], [3]),
             ([], [2, 1]),
             ([], [1, 1, 1]),
             ([3], [], []),
             ([2, 1], [], []),
             ([1, 1, 1], [], []),
             ([2], [1], []),
             ([1, 1], [1], []),
             ([2], [], [1]),
             ([1, 1], [], [1])]
        """
        for level in NN:
            for mu in PartitionTuples_level_size(level, self._size):
                yield self.element_class(self, list(mu))

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples(size=4).an_element()
            ([1], [1], [1], [1])
        """
        return self.element_class(self, tuple([1] for l in range(self._size)))


class PartitionTuples_level_size(PartitionTuples):
    """
    Class of partition tuples with a fixed level and a fixed size.
    """

    def __init__(self, level, size):
        r"""
        Initialize this class.

        EXAMPLES::

            sage: TestSuite( PartitionTuples(4,2) ).run()                               # needs sage.libs.flint sage.libs.pari
            sage: TestSuite( PartitionTuples(level=4, size=5) ).run()                   # needs sage.libs.flint sage.libs.pari
        """
        if not (level in NN and size in NN):
            raise ValueError('n and level must be nonnegative integers')
        super().__init__(category=FiniteEnumeratedSets())
        self._level = level
        self._size = size

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples(4,2)
            Partition tuples of level 4 and size 2
            sage: PartitionTuples(size=2, level=4)
            Partition tuples of level 4 and size 2
        """
        return 'Partition tuples of level {} and size {}'.format(self._level, self._size)

    def __contains__(self, mu):
        r"""
        Return ``True`` if ``mu`` is in ``self``.

        TESTS::

            sage: PartitionTuple([[3,2],[2]]) in PartitionTuples(2,7)
            True
            sage: PartitionTuple([[3,2],[],[],[],[2]]) in PartitionTuples(5,7)
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[2]]) in PartitionTuples(5,7)
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[3]]) in PartitionTuples(2,8)
            False
            sage: all(mu in PartitionTuples(3,8) for mu in PartitionTuples(3,8))        # needs sage.libs.flint
            True

        Check that :issue:`14145` is fixed::

            sage: 1 in PartitionTuples(5,7)
            False
        """
        if self._level == 1 and mu in _Partitions:
            return self._size == sum(mu)
        return (PartitionTuples.__contains__(self, mu)
                and self._level == len(mu)
                and self._size == sum(map(sum,mu)))

    def __iter__(self):
        r"""
        Iterate through the finite class of partition tuples of a fixed level
        and a fixed size.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: PartitionTuples(2,0).list() #indirect doctest
            [([], [])]
            sage: PartitionTuples(2,1).list() #indirect doctest
            [([1], []), ([], [1])]
            sage: PartitionTuples(2,2).list() #indirect doctest
            [([2], []), ([1, 1], []), ([1], [1]), ([], [2]), ([], [1, 1])]
            sage: PartitionTuples(3,2).list() #indirect doctest
            [([2], [], []),
             ([1, 1], [], []),
             ([1], [1], []),
             ([1], [], [1]),
             ([], [2], []),
             ([], [1, 1], []),
             ([], [1], [1]),
             ([], [], [2]),
             ([], [], [1, 1])]
        """
        p = [Partitions_n(i) for i in range(self._size+1)]
        for iv in IntegerVectors(self._size, self._level):
            for cp in itertools.product(*[p[i] for i in iv]):
                yield self._element_constructor_(cp)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples(level=4,size=4).an_element()
            ([1], [], [], [3])
        """
        mu = [[] for _ in itertools.repeat(None, self._level)]
        if self._size > 0:
            if self._level == 1:
                mu = [self._size-1,1]
            else:
                mu[0] = [1]
                mu[-1] = [self._size-1]
        return self.element_class(self, mu)

    def cardinality(self):
        r"""
        Return the number of ``level``-tuples of partitions of size ``n``.

        Wraps a pari function call using :pari:`eta`.

        EXAMPLES::

            sage: PartitionTuples(2,3).cardinality()                                    # needs sage.libs.pari
            10
            sage: PartitionTuples(2,8).cardinality()                                    # needs sage.libs.pari
            185

        TESTS:

        The following calls used to fail (:issue:`11476`)::

            sage: # needs sage.libs.pari
            sage: PartitionTuples(17,2).cardinality()
            170
            sage: PartitionTuples(2,17).cardinality()
            8470
            sage: PartitionTuples(100,13).cardinality()
            110320020147886800
            sage: PartitionTuples(13,90).cardinality()
            91506473741200186152352843611

        These answers were checked against Gap4 (the last of which takes an
        awful long time for gap to compute).
        """
        eta = pari(f'Ser(x,x,{self.size()})').eta()
        return ZZ((1 / eta**self.level()).polcoef(self.size(), pari('x')))

    def __setstate__(self, state):
        r"""
        In order to maintain backwards compatibility and be able to unpickle a
        old pickle from PartitionTuples_nk we have to override the default
        ``__setstate__``.

        TESTS::

            sage: loads(b"x\x9cM\x90\xcdN\xc30\x0c\x80\xd5\xc1\x06\xeb\x80\xf1{\xe0\r\xe0\xd2\x0b\x07\x1e\x02)B\x88\x9c-7\xb5\xba\xa8MR')\x12\x07$8p\xe0\xadq\x996q\xb1b\xfb\xb3\xf59\x9f3\x93\xb0\xa5\xca\x04W[\x8f\xb9\x1a0f\x9bm\xf0\xe5\xf3\xee\xf5:\x0e=%\xf0]\xc9\xc5\xfd\x17\xcf>\xf8\xe0N_\x83\xf5\xd2\xc5\x1e\xd0L\x10\xf46e>T\xba\x04r55\x8d\xf5-\xcf\x95p&\xf87\x8a\x19\x1c\xe5Mh\xc0\xa3#^(\xbd\x00\xd3`F>Rz\t\x063\xb5!\xbe\xf3\xf1\xd4\x98\x90\xc4K\xa5\x0b\xbf\xb5\x8b\xb2,U\xd6\x0bD\xb1t\xd8\x11\xec\x12.u\xf1\xf0\xfd\xc2+\xbd\x82\x96<E\xcc!&>Qz\x0e5&\xe2S\xa5\xd70X\xd3\xf5\x04\xe2\x91\xc4\x95\xcf\x9e\n\x11\xa3\x9e\x1c\xf9<\t\xa6\x1cG#\x83\xbcV\xfaf\x7f\xd9\xce\xfc\xef\xb4s\xa5o\xf7#\x13\x01\x03\xa6$!J\x81/~t\xd1m\xc4\xe5Q\\.\xff\xfd\x8e\t\x14\rmW\\\xa9\xb1\xae~\x01/\x8f\x85\x02")
            Partition tuples of level 7 and size 3
            sage: loads(dumps( PartitionTuples(7,3) ))  # indirect doctest for unpickling a Tableau element
            Partition tuples of level 7 and size 3
        """
        if isinstance(state, dict):   # for old pickles from Tableau_class
            parts = PartitionTuples(state['k'], state['n'])
            self.__class__ = parts.__class__
            self.__dict__ = parts.__dict__
        else:
            super().__setstate__(state)

###############################################################################
# Regular partition tuples


class RegularPartitionTuples(PartitionTuples):
    r"""
    Abstract base class for `\ell`-regular partition tuples.
    """

    def __init__(self, regular, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: RPT = PartitionTuples(regular=3)
            sage: TestSuite(RPT).run()                                                  # needs sage.libs.flint
        """
        if regular not in ZZ or regular < 1:
            raise ValueError("regular must be an integer greater than 1")
        self._ell = regular
        PartitionTuples.__init__(self, **kwds)

    def __contains__(self, mu):
        r"""
        Check if ``mu`` is an `\ell`-regular partition tuple.

        TESTS::

            sage: RPT = PartitionTuples(regular=2)
            sage: [[11,1], [2]] in RPT
            True
            sage: Partition([4,1]) in RPT
            True
            sage: [5,4,3,2,1] in RPT
            True
            sage: [[6,3,1], [], [], [3,1], [1], [1], [1]] in RPT
            True
            sage: [[10], [1], [1,1], [4,2]] in RPT
            False
            sage: [[5,2], [17, 1], [], [3,3,1], [1,1]] in RPT
            False
            sage: RPT = PartitionTuples(4,2,3)
            sage: elt = RPT([[1], [], [], [1]])
            sage: elt in RPT
            True
        """
        if not PartitionTuples.__contains__(self, mu):
            return False
        if isinstance(mu, Partition):
            return max(mu.to_exp() + [0]) < self._ell
        if isinstance(mu, PartitionTuple):
            return all(max(nu.to_exp() + [0]) < self._ell for nu in mu)
        if not mu:
            return True
        if mu in _Partitions:
            return all(mu.count(i) < self._ell for i in set(mu) if i > 0)
        return all(list(nu).count(i) < self._ell for nu in mu for i in set(nu) if i > 0)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples(regular=2).an_element()
            ([1], [], [], [2])
        """
        if self._level is None:
            lvl = 4
        else:
            lvl = self._level
        if self._size is None:
            size = 3
        else:
            size = self._size
        elt = RegularPartitionTuples_level_size(lvl, size, self._ell).an_element()
        return self.element_class(self, list(elt))


class RegularPartitionTuples_all(RegularPartitionTuples):
    r"""
    Class of `\ell`-regular partition tuples.
    """

    def __init__(self, regular):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: RPT = PartitionTuples(regular=3)
            sage: TestSuite(RPT).run()                                                  # needs sage.libs.flint
        """
        RegularPartitionTuples.__init__(self, regular, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples(regular=3)
            3-Regular partition tuples
        """
        return '{}-Regular partition tuples'.format(self._ell)

    def __iter__(self):
        r"""
        Iterate through the class of `\ell`-regular partition tuples.

        EXAMPLES::

            sage: PartitionTuples(regular=2)[:20]                                       # needs sage.libs.flint
            [([]),
             ([], []),
             ([1]),
             ([], [], []),
             ([1], []),
             ([], [1]),
             ([2]),
             ([], [], [], []),
             ([1], [], []),
             ([], [1], []),
             ([], [], [1]),
             ([2], []),
             ([1], [1]),
             ([], [2]),
             ([3]),
             ([2, 1]),
             ([], [], [], [], []),
             ([1], [], [], []),
             ([], [1], [], []),
             ([], [], [1], [])]
        """
        for N in NN:
            for size in range(N+1):
                for mu in RegularPartitionTuples_level_size(N-size+1, size, self._ell):
                    yield self.element_class(self, list(mu))


class RegularPartitionTuples_level(PartitionTuples_level):
    r"""
    Regular Partition tuples of a fixed level.

    INPUT:

    - ``level`` -- nonnegative integer; the level
    - ``regular`` -- positive integer or a tuple of nonnegative
      integers; if an integer, the highest multiplicity an entry may
      have in a component plus `1` with `0` representing `\infty`-regular
      (equivalently, partitions without restrictions)

    ``regular`` is a tuple of integers `(\ell_1, \ldots, \ell_k)` that
    specifies partition tuples `\mu` such that `\mu_i` is `\ell_i`-regular.
    If ``regular`` is an integer `\ell`, then we set `\ell_i = \ell` for
    all `i`.

    EXAMPLES::

        sage: RPT = PartitionTuples(level=4, regular=(2,3,0,2))
        sage: RPT[:24]                                                                  # needs sage.libs.flint
        [([], [], [], []),
         ([1], [], [], []),
         ([], [1], [], []),
         ([], [], [1], []),
         ([], [], [], [1]),
         ([2], [], [], []),
         ([1], [1], [], []),
         ([1], [], [1], []),
         ([1], [], [], [1]),
         ([], [2], [], []),
         ([], [1, 1], [], []),
         ([], [1], [1], []),
         ([], [1], [], [1]),
         ([], [], [2], []),
         ([], [], [1, 1], []),
         ([], [], [1], [1]),
         ([], [], [], [2]),
         ([3], [], [], []),
         ([2, 1], [], [], []),
         ([2], [1], [], []),
         ([2], [], [1], []),
         ([2], [], [], [1]),
         ([1], [2], [], []),
         ([1], [1, 1], [], [])]
        sage: [[1,1],[3],[5,5,5],[7,2]] in RPT
        False
        sage: [[3,1],[3],[5,5,5],[7,2]] in RPT
        True
        sage: [[3,1],[3],[5,5,5]] in RPT
        False
    """

    def __init__(self, level, regular):
        r"""
        Initialize ``self``.

        TESTS::

            sage: RPT = PartitionTuples(level=2, regular=(0,0))
            sage: RPT.category()
            Category of infinite enumerated sets
            sage: RPT = PartitionTuples(level=4, regular=3)
            sage: TestSuite(RPT).run()                                                  # needs sage.libs.flint
        """
        if level not in NN:
            raise ValueError('level must be a nonnegative integer')
        if not isinstance(regular, tuple):
            # This should not happen if called from RegularPartitionTuples
            regular = (regular,) * level
        if any(r != 1 for r in regular):
            category = InfiniteEnumeratedSets()
        else:
            category = FiniteEnumeratedSets()
        if any(r not in NN for r in regular):
            raise ValueError('regular must be a tuple of nonnegative integers')
        if len(regular) != level:
            raise ValueError("regular must be a tuple with length {}".format(level))
        PartitionTuples_level.__init__(self, level, category=category)
        self._ell = regular

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples(level=4, regular=3)
            3-Regular partition tuples of level 4
            sage: PartitionTuples(level=4, regular=(2,3,0,2))
            (2, 3, 0, 2)-Regular partition tuples of level 4
        """
        if self._ell[1:] == self._ell[:-1]:
            return '{}-Regular partition tuples of level {}'.format(self._ell[0],
                                                                    self._level)
        return '{}-Regular partition tuples of level {}'.format(self._ell,
                                                                self._level)

    def __contains__(self, mu):
        r"""
        Return ``True`` if ``mu`` is in ``self``.

        TESTS::

            sage: RPT = PartitionTuples(level=4, regular=2)
            sage: [[4,2,1], [], [2], [2]] in RPT
            True
            sage: [[10], [1], [1,1], [4,2]] in RPT
            False
            sage: [[5,2], [], [3,3,1], [1,1]] in RPT
            False
            sage: [4, 3, 2] in RPT
            False

            sage: RPT = PartitionTuples(level=3, regular=(2,1,4))
            sage: [[4], [2], [5]] in RPT
            False
            sage: [[4], [], [5]] in RPT
            True
            sage: [[4,3], [], [5]] in RPT
            True
            sage: [[4,4], [], [5]] in RPT
            False
            sage: [[4,3], [5]] in RPT
            False
            sage: [5, 4, 3] in RPT
            False
            sage: [] in RPT
            False
            sage: [[], [], []] in RPT
            True
            sage: [[], [], [], [2]] in RPT
            False

            sage: from sage.combinat.partition_tuple import RegularPartitionTuples_level
            sage: RPT = RegularPartitionTuples_level(1, (3,)); RPT
            3-Regular partition tuples of level 1
            sage: [[2,2]] in RPT
            True
            sage: [[2,2,2]] in RPT
            False
        """
        if self._level == 1:
            try:
                if mu[0] in ZZ:
                    return mu in RegularPartitions_all(self._ell[0])
            except (TypeError, ValueError):
                return False
            return mu[0] in RegularPartitions_all(self._ell[0])
        if mu not in PartitionTuples_level(self._level):
            return False
        if isinstance(mu, Partition):  # it is level 1
            return False
        if isinstance(mu, PartitionTuple):
            return all(max(p.to_exp() + [0]) < ell for p, ell in zip(mu, self._ell)
                       if ell > 0)
        return all(p in RegularPartitions_all(ell) for p, ell in zip(mu, self._ell)
                   if ell > 0)

    def __iter__(self):
        r"""
        Iterate through the class of `\ell`-regular partition tuples
        of a fixed level.

        EXAMPLES::

            sage: PartitionTuples(level=3, regular=(2,1,4))[:24]                        # needs sage.libs.flint
            [([], [], []),
             ([1], [], []),
             ([], [], [1]),
             ([2], [], []),
             ([1], [], [1]),
             ([], [], [2]),
             ([], [], [1, 1]),
             ([3], [], []),
             ([2, 1], [], []),
             ([2], [], [1]),
             ([1], [], [2]),
             ([1], [], [1, 1]),
             ([], [], [3]),
             ([], [], [2, 1]),
             ([], [], [1, 1, 1]),
             ([4], [], []),
             ([3, 1], [], []),
             ([3], [], [1]),
             ([2, 1], [], [1]),
             ([2], [], [2]),
             ([2], [], [1, 1]),
             ([1], [], [3]),
             ([1], [], [2, 1]),
             ([1], [], [1, 1, 1])]
            sage: PartitionTuples(level=4, regular=2)[:20]                              # needs sage.libs.flint
            [([], [], [], []),
             ([1], [], [], []),
             ([], [1], [], []),
             ([], [], [1], []),
             ([], [], [], [1]),
             ([2], [], [], []),
             ([1], [1], [], []),
             ([1], [], [1], []),
             ([1], [], [], [1]),
             ([], [2], [], []),
             ([], [1], [1], []),
             ([], [1], [], [1]),
             ([], [], [2], []),
             ([], [], [1], [1]),
             ([], [], [], [2]),
             ([3], [], [], []),
             ([2, 1], [], [], []),
             ([2], [1], [], []),
             ([2], [], [1], []),
             ([2], [], [], [1])]
        """
        for size in NN:
            for mu in RegularPartitionTuples_level_size(self._level, size, self._ell):
                yield self.element_class(self, list(mu))


class RegularPartitionTuples_size(RegularPartitionTuples):
    r"""
    Class of `\ell`-regular partition tuples with a fixed size.
    """

    def __init__(self, size, regular):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: RPT = PartitionTuples(size=4, regular=3)
            sage: TestSuite(RPT).run()                                                  # needs sage.libs.flint
        """
        if size not in NN:
            raise ValueError('size must be a nonnegative integer')
        RegularPartitionTuples.__init__(self, regular, category=InfiniteEnumeratedSets())
        self._size = size

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples(size=4, regular=3)
            3-Regular partition tuples of size 4
        """
        return '{}-Regular partition tuples of size {}'.format(self._ell, self._size)

    def __contains__(self, mu):
        r"""
        Return ``True`` if ``mu`` is in ``self``.

        TESTS::

            sage: RPT = PartitionTuples(size=4, regular=2)
            sage: [[2, 1], [1]] in RPT
            True
            sage: [3, 1] in RPT
            True
            sage: [[1], [], [], [2,1]] in RPT
            True
            sage: [[1], [1], [1], [1]] in RPT
            True
            sage: [[1], [1,1,1]] in RPT
            False
            sage: [[2,1,1]] in RPT
            False
            sage: [2,1,1] in RPT
            False
            sage: RPT = PartitionTuples(size=7, regular=2)
            sage: [[], [3,2,2,1], [1], [1]] in RPT
            False
            sage: RPT = PartitionTuples(size=9, regular=2)
            sage: [4, 3, 2] in RPT
            True
        """
        return ((mu in RegularPartitions_all(self._ell)
                 and self._size == sum(mu))
                or (RegularPartitionTuples.__contains__(self, mu)
                    and self._size == sum(map(sum, mu))))

    def __iter__(self):
        r"""
        Iterate through the class of `\ell`-regular partition tuples
        of a fixed size.

        EXAMPLES::

            sage: PartitionTuples(size=4, regular=2)[:10]                               # needs sage.libs.flint
            [([4]),
             ([3, 1]),
             ([4], []),
             ([3, 1], []),
             ([3], [1]),
             ([2, 1], [1]),
             ([2], [2]),
             ([1], [3]),
             ([1], [2, 1]),
             ([], [4])]
        """
        for level in PositiveIntegers():
            for mu in RegularPartitionTuples_level_size(level, self._size, self._ell):
                yield self.element_class(self, list(mu))


class RegularPartitionTuples_level_size(PartitionTuples_level_size):
    r"""
    Class of `\ell`-regular partition tuples with a fixed level and
    a fixed size.

    INPUT:

    - ``level`` -- nonnegative integer; the level
    - ``size`` -- nonnegative integer; the size
    - ``regular`` -- positive integer or a tuple of nonnegative
      integers; if an integer, the highest multiplicity an entry may
      have in a component plus `1` with `0` representing `\infty`-regular
      (equivalently, partitions without restrictions)

    ``regular`` is a tuple of integers `(\ell_1, \ldots, \ell_k)` that
    specifies partition tuples `\mu` such that `\mu_i` is `\ell_i`-regular.
    If ``regular`` is an integer `\ell`, then we set `\ell_i = \ell` for
    all `i`.

    EXAMPLES::

        sage: PartitionTuples(level=3, size=7, regular=(2,1,3))[0:24]                   # needs sage.libs.flint
        [([7], [], []),
         ([6, 1], [], []),
         ([5, 2], [], []),
         ([4, 3], [], []),
         ([4, 2, 1], [], []),
         ([6], [], [1]),
         ([5, 1], [], [1]),
         ([4, 2], [], [1]),
         ([3, 2, 1], [], [1]),
         ([5], [], [2]),
         ([5], [], [1, 1]),
         ([4, 1], [], [2]),
         ([4, 1], [], [1, 1]),
         ([3, 2], [], [2]),
         ([3, 2], [], [1, 1]),
         ([4], [], [3]),
         ([4], [], [2, 1]),
         ([3, 1], [], [3]),
         ([3, 1], [], [2, 1]),
         ([3], [], [4]),
         ([3], [], [3, 1]),
         ([3], [], [2, 2]),
         ([3], [], [2, 1, 1]),
         ([2, 1], [], [4])]
    """

    def __init__(self, level, size, regular):
        r"""
        Initialize ``self``.

        TESTS::

            sage: RPT = PartitionTuples(4,2,3)
            sage: TestSuite(RPT).run()                                                  # needs sage.libs.flint sage.libs.pari
        """
        if size not in NN:
            raise ValueError('size must be a nonnegative integer')
        if not (level in ZZ and level > 0):
            raise ValueError('level must be a positive integer')
        if not isinstance(regular, tuple):
            # This should not happen if called from RegularPartitionTuples
            regular = (regular,) * level
        if len(regular) != level:
            raise ValueError(f'regular must be a list with length {level}')
        if any(i not in NN for i in regular):
            raise ValueError('regular must be a list of nonnegative integers')
        PartitionTuples_level_size.__init__(self, level, size)
        self._ell = regular

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionTuples(level=3, size=7, regular=(2,1,4))
            (2, 1, 4)-Regular partition tuples of level 3 and size 7
            sage: PartitionTuples(4,2,3)
            3-Regular partition tuples of level 4 and size 2
            sage: PartitionTuples(size=2, level=4, regular=3)
            3-Regular partition tuples of level 4 and size 2
        """
        if self._ell[1:] == self._ell[:-1]:
            return '{}-Regular partition tuples of level {} and size {}'.format(
                                         self._ell[0], self._level, self._size)
        return '{}-Regular partition tuples of level {} and size {}'.format(
                                            self._ell, self._level, self._size)

    def __contains__(self, mu):
        r"""
        Return ``True`` if `\mu` is in ``self``.

        TESTS::

            sage: RPT = PartitionTuples(level=3, size=7, regular=(2,1,4))
            sage: RPT
            (2, 1, 4)-Regular partition tuples of level 3 and size 7
            sage: [[3,1],[],[3]] in RPT
            True
            sage: [[3],[1],[3]] in RPT
            False
            sage: [[3,2],[],[3]] in RPT
            False
            sage: [[3,3],[],[1]] in RPT
            False
            sage: RPT = PartitionTuples(4,3,2)
            sage: [[], [], [2], [1]] in RPT
            True
            sage: [[1], [1], [], [1]] in RPT
            True
            sage: [[1,1,1], [], [], []] in RPT
            False
            sage: RPT = PartitionTuples(9, 3, 2)
            sage: [4, 3, 2] in RPT
            False
        """
        if mu not in RegularPartitionTuples_level(self._level, self._ell):
            return False
        return self._size == sum(map(sum, mu))

    def __iter__(self):
        r"""
        Iterate through the finite class of `\ell`-regular partition tuples
        of a fixed level and a fixed size.

        EXAMPLES::

            sage: list(PartitionTuples(3,3,2))                                          # needs sage.libs.pari
            [([3], [], []),
             ([2, 1], [], []),
             ([2], [1], []),
             ([2], [], [1]),
             ([1], [2], []),
             ([1], [1], [1]),
             ([1], [], [2]),
             ([], [3], []),
             ([], [2, 1], []),
             ([], [2], [1]),
             ([], [1], [2]),
             ([], [], [3]),
             ([], [], [2, 1])]
        """
        for iv in IntegerVectors(self._size, self._level):
            p = [RegularPartitions_n(v, ell) if ell > 0 else Partitions_n(v)
                 for v, ell in zip(iv, self._ell)]
            for cp in itertools.product(*[p[i] for i in range(self._level)]):
                yield self._element_constructor_(cp)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples(level=4, size=4, regular=3).an_element()
            ([1], [], [], [3])
        """
        mu = [[] for _ in itertools.repeat(None, self._level)]
        if self._size > 0:
            if self._level == 1:
                mu = [[self._size - 1, 1]]
            else:
                mu[0] = [1]
                mu[-1] = [self._size - 1]
        return self.element_class(self, mu)
