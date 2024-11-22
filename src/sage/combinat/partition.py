r"""
Integer partitions

A partition `p` of a nonnegative integer `n` is a
non-increasing list of positive integers (the *parts* of the
partition) with total sum `n`.

A partition can be depicted by a diagram made of rows of cells,
where the number of cells in the `i`-th row starting from
the top is the `i`-th part of the partition.

The coordinate system related to a partition applies from the top
to the bottom and from left to right. So, the corners of the
partition `[5, 3, 1]` are `[[0,4], [1,2], [2,0]]`.

For display options, see :obj:`Partitions.options`.

.. NOTE::

    - Boxes is a synonym for cells. All methods will use 'cell' and 'cells'
      instead of 'box' and 'boxes'.

    - Partitions are 0 based with coordinates in the form of (row-index,
      column-index).

    - If given coordinates of the form ``(r, c)``, then use Python's
      \*-operator.

    - Throughout this documentation, for a partition `\lambda` we will denote
      its conjugate partition by `\lambda^{\prime}`. For more on conjugate
      partitions, see :meth:`Partition.conjugate()`.

    - The comparisons on partitions use lexicographic order.

.. NOTE::

    We use the convention that lexicographic ordering is read from
    left-to-right. That is to say `[1, 3, 7]` is smaller than `[2, 3, 4]`.

AUTHORS:

- Mike Hansen (2007): initial version

- Dan Drake (2009-03-28): deprecate RestrictedPartitions and implement
  Partitions_parts_in

- Travis Scrimshaw (2012-01-12): Implemented latex function to Partition_class

- Travis Scrimshaw (2012-05-09): Fixed Partitions(-1).list() infinite recursion
  loop by saying Partitions_n is the empty set.

- Travis Scrimshaw (2012-05-11): Fixed bug in inner where if the length was
  longer than the length of the inner partition, it would include 0's.

- Andrew Mathas (2012-06-01): Removed deprecated functions and added
  compatibility with the PartitionTuple classes.  See :issue:`13072`

- Travis Scrimshaw (2012-10-12): Added options. Made
  ``Partition_class`` to the element ``Partition``. ``Partitions*`` are now
  all in the category framework except ``PartitionsRestricted`` (which will
  eventually be removed). Cleaned up documentation.

- Matthew Lancellotti (2018-09-14): Added a bunch of "k" methods to Partition.

EXAMPLES:

There are `5` partitions of the integer `4`::

    sage: Partitions(4).cardinality()                                                   # needs sage.libs.flint
    5
    sage: Partitions(4).list()
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

We can use the method ``.first()`` to get the 'first' partition of a
number::

    sage: Partitions(4).first()
    [4]

Using the method ``.next(p)``, we can calculate the 'next' partition
after `p`. When we are at the last partition, ``None`` will be returned::

    sage: Partitions(4).next([4])
    [3, 1]
    sage: Partitions(4).next([1,1,1,1]) is None
    True

We can use ``iter`` to get an object which iterates over the partitions
one by one to save memory.  Note that when we do something like
``for part in Partitions(4)`` this iterator is used in the background::

    sage: g = iter(Partitions(4))
    sage: next(g)
    [4]
    sage: next(g)
    [3, 1]
    sage: next(g)
    [2, 2]
    sage: for p in Partitions(4): print(p)
    [4]
    [3, 1]
    [2, 2]
    [2, 1, 1]
    [1, 1, 1, 1]

We can add constraints to the type of partitions we want. For
example, to get all of the partitions of `4` of length `2`, we'd do the
following::

    sage: Partitions(4, length=2).list()
    [[3, 1], [2, 2]]

Here is the list of partitions of length at least `2` and the list of
ones with length at most `2`::

    sage: Partitions(4, min_length=2).list()
    [[3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    sage: Partitions(4, max_length=2).list()
    [[4], [3, 1], [2, 2]]

The options ``min_part`` and ``max_part`` can be used to set constraints
on the sizes of all parts. Using ``max_part``, we can select
partitions having only 'small' entries. The following is the list
of the partitions of `4` with parts at most `2`::

    sage: Partitions(4, max_part=2).list()
    [[2, 2], [2, 1, 1], [1, 1, 1, 1]]

The ``min_part`` options is complementary to ``max_part`` and selects
partitions having only 'large' parts. Here is the list of all
partitions of `4` with each part at least `2`::

    sage: Partitions(4, min_part=2).list()
    [[4], [2, 2]]

The options ``inner`` and ``outer`` can be used to set part-by-part
constraints. This is the list of partitions of `4` with ``[3, 1, 1]`` as
an outer bound (that is, partitions of `4` contained in the partition
``[3, 1, 1]``)::

    sage: Partitions(4, outer=[3,1,1]).list()
    [[3, 1], [2, 1, 1]]

``outer`` sets ``max_length`` to the length of its argument. Moreover, the
parts of ``outer`` may be infinite to clear constraints on specific
parts. Here is the list of the partitions of `4` of length at most `3`
such that the second and third part are `1` when they exist::

    sage: Partitions(4, outer=[oo,1,1]).list()
    [[4], [3, 1], [2, 1, 1]]

Finally, here are the partitions of `4` with ``[1,1,1]`` as an inner
bound (i. e., the partitions of `4` containing the partition ``[1,1,1]``).
Note that ``inner`` sets ``min_length`` to the length of its argument::

    sage: Partitions(4, inner=[1,1,1]).list()
    [[2, 1, 1], [1, 1, 1, 1]]

The options ``min_slope`` and ``max_slope`` can be used to set
constraints on the slope, that is on the difference ``p[i+1]-p[i]`` of
two consecutive parts. Here is the list of the strictly decreasing
partitions of `4`::

    sage: Partitions(4, max_slope=-1).list()
    [[4], [3, 1]]

The constraints can be combined together in all reasonable ways.
Here are all the partitions of `11` of length between `2` and `4` such
that the difference between two consecutive parts is between `-3` and
`-1`::

    sage: Partitions(11, min_slope=-3, max_slope=-1, min_length=2, max_length=4).list()
    [[7, 4], [6, 5], [6, 4, 1], [6, 3, 2], [5, 4, 2], [5, 3, 2, 1]]

Partition objects can also be created individually with :class:`Partition`::

    sage: Partition([2,1])
    [2, 1]

Once we have a partition object, then there are a variety of
methods that we can use. For example, we can get the conjugate of a
partition. Geometrically, the conjugate of a partition is the
reflection of that partition through its main diagonal. Of course,
this operation is an involution::

    sage: Partition([4,1]).conjugate()
    [2, 1, 1, 1]
    sage: Partition([4,1]).conjugate().conjugate()
    [4, 1]

If we create a partition with extra zeros at the end, they will be dropped::

    sage: Partition([4,1,0,0])
    [4, 1]
    sage: Partition([0])
    []
    sage: Partition([0,0])
    []

The idea of a partition being followed by infinitely many parts of size
`0` is consistent with the ``get_part`` method::

    sage: p = Partition([5, 2])
    sage: p.get_part(0)
    5
    sage: p.get_part(10)
    0

We can go back and forth between the standard and the exponential
notations of a partition. The exponential notation can be padded with
extra zeros::

    sage: Partition([6,4,4,2,1]).to_exp()
    [1, 1, 0, 2, 0, 1]
    sage: Partition(exp=[1,1,0,2,0,1])
    [6, 4, 4, 2, 1]
    sage: Partition([6,4,4,2,1]).to_exp(5)
    [1, 1, 0, 2, 0, 1]
    sage: Partition([6,4,4,2,1]).to_exp(7)
    [1, 1, 0, 2, 0, 1, 0]
    sage: Partition([6,4,4,2,1]).to_exp(10)
    [1, 1, 0, 2, 0, 1, 0, 0, 0, 0]

We can get the (zero-based!) coordinates of the corners of a
partition::

    sage: Partition([4,3,1]).corners()
    [(0, 3), (1, 2), (2, 0)]

We can compute the core and quotient of a partition and build
the partition back up from them::

    sage: Partition([6,3,2,2]).core(3)
    [2, 1, 1]
    sage: Partition([7,7,5,3,3,3,1]).quotient(3)
    ([2], [1], [2, 2, 2])
    sage: p = Partition([11,5,5,3,2,2,2])
    sage: p.core(3)
    []
    sage: p.quotient(3)
    ([2, 1], [4], [1, 1, 1])
    sage: Partition(core=[],quotient=([2, 1], [4], [1, 1, 1]))
    [11, 5, 5, 3, 2, 2, 2]

We can compute the `0-1` sequence and go back and forth::

    sage: Partitions().from_zero_one([1, 1, 1, 1, 0, 1, 0])
    [5, 4]
    sage: all(Partitions().from_zero_one(mu.zero_one_sequence())
    ....:     == mu for n in range(5) for mu in Partitions(n))
    True

We can compute the Frobenius coordinates and go back and forth::

    sage: Partition([7,3,1]).frobenius_coordinates()
    ([6, 1], [2, 0])
    sage: Partition(frobenius_coordinates=([6,1],[2,0]))
    [7, 3, 1]
    sage: all(mu == Partition(frobenius_coordinates=mu.frobenius_coordinates())
    ....:     for n in range(12) for mu in Partitions(n))
    True

We use the lexicographic ordering::

    sage: pl = Partition([4,1,1])
    sage: ql = Partitions()([3,3])
    sage: pl > ql
    True
    sage: PL = Partitions()
    sage: pl = PL([4,1,1])
    sage: ql = PL([3,3])
    sage: pl > ql
    True
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from itertools import accumulate

from sage.arith.misc import binomial, factorial, gcd, multinomial
from sage.structure.global_options import GlobalOptions
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.misc.prandom import randrange
from sage.misc.cachefunc import cached_method, cached_function

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.infinity import infinity

from .combinat import CombinatorialElement
from . import tableau
from . import permutation
from . import composition
from sage.combinat.partitions import ZS1_iterator, ZS1_iterator_nk, ZS1_next, ZS2_next
from sage.combinat.integer_lists import IntegerListsLex
from sage.combinat.integer_lists.invlex import IntegerListsBackend_invlex
from sage.combinat.integer_vector_weighted import iterator_fast as weighted_iterator_fast
from sage.combinat.combinat_cython import conjugate
from sage.combinat.combinatorial_map import combinatorial_map

lazy_import('sage.combinat.skew_partition', 'SkewPartition')
lazy_import('sage.combinat.partition_tuple', 'PartitionTuple')
lazy_import('sage.combinat.root_system.weyl_group', 'WeylGroup')
lazy_import('sage.libs.pari.all', 'pari')
lazy_import('sage.groups.perm_gps.permgroup', 'PermutationGroup')
lazy_import("sage.symbolic.ring", "var")


class Partition(CombinatorialElement):
    r"""
    A partition `p` of a nonnegative integer `n` is a
    non-increasing list of positive integers (the *parts* of the
    partition) with total sum `n`.

    A partition is often represented as a diagram consisting of **cells**,
    or **boxes**, placed in rows on top of each other such that the number of
    cells in the `i`-th row, reading from top to bottom, is the `i`-th
    part of the partition. The rows are left-justified (and become shorter
    and shorter the farther down one goes). This diagram is called the
    **Young diagram** of the partition, or more precisely its Young diagram
    in English notation. (French and Russian notations are variations on this
    representation.)

    The coordinate system related to a partition applies from the top
    to the bottom and from left to right. So, the corners of the
    partition ``[5, 3, 1]`` are ``[[0,4], [1,2], [2,0]]``.

    For display options, see :meth:`Partitions.options`.

    .. NOTE::

        Partitions are 0 based with coordinates in the form of (row-index,
        column-index). For example consider the partition
        ``mu=Partition([4,3,2,2])``, the first part is ``mu[0]`` (which is 4),
        the second is ``mu[1]``, and so on, and the upper-left cell in English
        convention is ``(0, 0)``.

    A partition can be specified in one of the following ways:

    - a list (the default)
    - using exponential notation
    - by Frobenius coordinates
    - specifying its `0-1` sequence
    - specifying the core and the quotient

    See the examples below.

    EXAMPLES:

    Creating partitions though parents::

        sage: mu = Partitions(8)([3,2,1,1,1]); mu
        [3, 2, 1, 1, 1]
        sage: nu = Partition([3,2,1,1,1]); nu
        [3, 2, 1, 1, 1]
        sage: mu == nu
        True
        sage: mu is nu
        False
        sage: mu in Partitions()
        True
        sage: mu.parent()
        Partitions of the integer 8
        sage: mu.size()
        8
        sage: mu.category()
        Category of elements of Partitions of the integer 8
        sage: nu.parent()
        Partitions
        sage: nu.category()
        Category of elements of Partitions
        sage: mu[0]
        3
        sage: mu[1]
        2
        sage: mu[2]
        1
        sage: mu.pp()
        ***
        **
        *
        *
        *
        sage: mu.removable_cells()
        [(0, 2), (1, 1), (4, 0)]
        sage: mu.down_list()
        [[2, 2, 1, 1, 1], [3, 1, 1, 1, 1], [3, 2, 1, 1]]
        sage: mu.addable_cells()
        [(0, 3), (1, 2), (2, 1), (5, 0)]
        sage: mu.up_list()
        [[4, 2, 1, 1, 1], [3, 3, 1, 1, 1], [3, 2, 2, 1, 1], [3, 2, 1, 1, 1, 1]]
        sage: mu.conjugate()
        [5, 2, 1]
        sage: mu.dominates(nu)
        True
        sage: nu.dominates(mu)
        True

    Creating partitions using ``Partition``::

        sage: Partition([3,2,1])
        [3, 2, 1]
        sage: Partition(exp=[2,1,1])
        [3, 2, 1, 1]
        sage: Partition(core=[2,1], quotient=[[2,1],[3],[1,1,1]])
        [11, 5, 5, 3, 2, 2, 2]
        sage: Partition(frobenius_coordinates=([3,2],[4,0]))
        [4, 4, 1, 1, 1]
        sage: Partitions().from_zero_one([1, 1, 1, 1, 0, 1, 0])
        [5, 4]
        sage: [2,1] in Partitions()
        True
        sage: [2,1,0] in Partitions()
        True
        sage: Partition([1,2,3])
        Traceback (most recent call last):
        ...
        ValueError: [1, 2, 3] is not an element of Partitions

    Sage ignores trailing zeros at the end of partitions::

        sage: Partition([3,2,1,0])
        [3, 2, 1]
        sage: Partitions()([3,2,1,0])
        [3, 2, 1]
        sage: Partitions(6)([3,2,1,0])
        [3, 2, 1]

    TESTS:

    Check that only trailing zeros are stripped::

        sage: TestSuite( Partition([]) ).run()
        sage: TestSuite( Partition([4,3,2,2,2,1]) ).run()
        sage: Partition([3,2,2,2,1,0,0,0])
        [3, 2, 2, 2, 1]
        sage: Partition([3,0,2,2,2,1,0])
        Traceback (most recent call last):
        ...
        ValueError: [3, 0, 2, 2, 2, 1, 0] is not an element of Partitions
        sage: Partition([0,7,3])
        Traceback (most recent call last):
        ...
        ValueError: [0, 7, 3] is not an element of Partitions
    """
    @staticmethod
    def __classcall_private__(cls, mu=None, **keyword):
        """
        This constructs a list from optional arguments and delegates the
        construction of a :class:`Partition` to the ``element_class()`` call
        of the appropriate parent.

        EXAMPLES::

            sage: Partition([3,2,1])
            [3, 2, 1]
            sage: Partition(exp=[2,1,1])
            [3, 2, 1, 1]
            sage: Partition(core=[2,1], quotient=[[2,1],[3],[1,1,1]])
            [11, 5, 5, 3, 2, 2, 2]
        """
        l = len(keyword)
        if l == 0:
            if mu is not None:
                if isinstance(mu, Partition):
                    return mu
                return _Partitions(list(mu))
        if l == 1:
            if 'beta_numbers' in keyword:
                return _Partitions.from_beta_numbers(keyword['beta_numbers'])
            elif 'exp' in keyword:
                return _Partitions.from_exp(keyword['exp'])
            elif 'frobenius_coordinates' in keyword:
                return _Partitions.from_frobenius_coordinates(keyword['frobenius_coordinates'])
            elif 'zero_one' in keyword:
                return _Partitions.from_zero_one(keyword['zero_one'])

        if l == 2 and 'core' in keyword and 'quotient' in keyword:
            return _Partitions.from_core_and_quotient(keyword['core'], keyword['quotient'])
        raise ValueError('incorrect syntax for Partition()')

    def __setstate__(self, state):
        r"""
        In order to maintain backwards compatibility and be able to unpickle a
        old pickle from ``Partition_class`` we have to override the default
        ``__setstate__``.

        EXAMPLES::

            sage: loads(b'x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+H,*\xc9,\xc9\xcc\xcf\xe3\n\x80\xb1\xe2\x93s\x12\x8b\x8b\xb9\n\x195\x1b\x0b\x99j\x0b\x995BY\xe33\x12\x8b3\nY\xfc\x80\xac\x9c\xcc\xe2\x92B\xd6\xd8B6\r\x88IE\x99y\xe9\xc5z\x99y%\xa9\xe9\xa9E\\\xb9\x89\xd9\xa9\xf10N!{(\xa3qkP!G\x06\x90a\x04dp\x82\x18\x86@\x06Wji\x92\x1e\x00x0.\xb5')
            [3, 2, 1]
            sage: loads(dumps( Partition([3,2,1]) ))  # indirect doctest
            [3, 2, 1]
        """
        if isinstance(state, dict):   # for old pickles from Partition_class
            self._set_parent(_Partitions)
            self.__dict__ = state
        else:
            self._set_parent(state[0])
            self.__dict__ = state[1]

    def __init__(self, parent, mu):
        """
        Initialize ``self``.

        We assume that ``mu`` is a weakly decreasing list of
        nonnegative elements in ``ZZ``.

        EXAMPLES::

            sage: p = Partition([3,1])
            sage: TestSuite(p).run()

        TESTS:

        Fix that tuples raise the correct error::

            sage: Partition((3,1,7))
            Traceback (most recent call last):
            ...
            ValueError: [3, 1, 7] is not an element of Partitions
        """
        if isinstance(mu, Partition):
            # since we are (suppose to be) immutable, we can share the underlying data
            CombinatorialElement.__init__(self, parent, mu._list)
        else:
            if mu and not mu[-1]:
                # direct callers might assume that mu is not modified
                mu = mu[:-1]
                while mu and not mu[-1]:
                    mu.pop()
            CombinatorialElement.__init__(self, parent, mu)

    @cached_method
    def __hash__(self):
        r"""
        Return the hash of ``self``.

        TESTS::

            sage: P = Partition([4,2,2,1])
            sage: hash(P) == hash(P)
            True
        """
        return hash(tuple(self._list))

    def _repr_(self):
        r"""
        Return a string representation of ``self`` depending on
        :meth:`Partitions.options`.

        EXAMPLES::

            sage: mu=Partition([7,7,7,3,3,2,1,1,1,1,1,1,1]); mu # indirect doctest
            [7, 7, 7, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1]
            sage: Partitions.options.display="diagram"; mu
            *******
            *******
            *******
            ***
            ***
            **
            *
            *
            *
            *
            *
            *
            *
            sage: Partitions.options.display="list"; mu
            [7, 7, 7, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1]
            sage: Partitions.options.display="compact_low"; mu
            1^7,2,3^2,7^3
            sage: Partitions.options.display="compact_high"; mu
            7^3,3^2,2,1^7
            sage: Partitions.options.display="exp_low"; mu
            1^7, 2, 3^2, 7^3
            sage: Partitions.options.display="exp_high"; mu
            7^3, 3^2, 2, 1^7

            sage: Partitions.options.convention="French"
            sage: mu = Partition([7,7,7,3,3,2,1,1,1,1,1,1,1]); mu # indirect doctest
            7^3, 3^2, 2, 1^7
            sage: Partitions.options.display="diagram"; mu
            *
            *
            *
            *
            *
            *
            *
            **
            ***
            ***
            *******
            *******
            *******
            sage: Partitions.options.display="list"; mu
            [7, 7, 7, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1]
            sage: Partitions.options.display="compact_low"; mu
            1^7,2,3^2,7^3
            sage: Partitions.options.display="compact_high"; mu
            7^3,3^2,2,1^7
            sage: Partitions.options.display="exp_low"; mu
            1^7, 2, 3^2, 7^3
            sage: Partitions.options.display="exp_high"; mu
            7^3, 3^2, 2, 1^7

            sage: Partitions.options._reset()
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _ascii_art_(self):
        """
        TESTS::

            sage: ascii_art(Partitions(5).list())
            [                                * ]
            [                            **  * ]
            [                   ***  **  *   * ]
            [        ****  ***  *    **  *   * ]
            [ *****, *   , ** , *  , * , * , * ]
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._repr_diagram().splitlines(), baseline=0)

    def _unicode_art_(self):
        """
        TESTS::

            sage: unicode_art(Partitions(5).list())
            ⎡                                      ┌┐ ⎤
            ⎢                                 ┌┬┐  ├┤ ⎥
            ⎢                      ┌┬┬┐  ┌┬┐  ├┼┘  ├┤ ⎥
            ⎢         ┌┬┬┬┐  ┌┬┬┐  ├┼┴┘  ├┼┤  ├┤   ├┤ ⎥
            ⎢ ┌┬┬┬┬┐  ├┼┴┴┘  ├┼┼┘  ├┤    ├┼┘  ├┤   ├┤ ⎥
            ⎣ └┴┴┴┴┘, └┘   , └┴┘ , └┘  , └┘ , └┘ , └┘ ⎦
            sage: Partitions.options.convention = "French"
            sage: unicode_art(Partitions(5).list())
            ⎡                                      ┌┐ ⎤
            ⎢                                 ┌┐   ├┤ ⎥
            ⎢                      ┌┐    ┌┐   ├┤   ├┤ ⎥
            ⎢         ┌┐     ┌┬┐   ├┤    ├┼┐  ├┤   ├┤ ⎥
            ⎢ ┌┬┬┬┬┐  ├┼┬┬┐  ├┼┼┐  ├┼┬┐  ├┼┤  ├┼┐  ├┤ ⎥
            ⎣ └┴┴┴┴┘, └┴┴┴┘, └┴┴┘, └┴┴┘, └┴┘, └┴┘, └┘ ⎦
            sage: Partitions.options._reset()
        """
        from sage.typeset.unicode_art import UnicodeArt

        if not self._list:
            return UnicodeArt('∅', baseline=0)
        if self.parent().options.convention == "English":
            data = list(self)
        else:
            data = list(reversed(self))

        txt = ['┌' + '┬' * (data[0] - 1) + '┐']
        for i in range(len(data) - 1):
            p = data[i]
            q = data[i + 1]
            if p < q:
                txt += ['├' + '┼' * p + '┬' * (q - p - 1) + '┐']
            elif p == q:
                txt += ['├' + '┼' * (p - 1) + '┤']
            else:
                txt += ['├' + '┼' * q + '┴' * (p - q - 1) + '┘']
        txt += ['└' + '┴' * (data[-1] - 1) + '┘']

        return UnicodeArt(txt, baseline=0)

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list.

        EXAMPLES::

            sage: print(Partition([7,7,7,3,3,2,1,1,1,1,1,1,1])._repr_list())
            [7, 7, 7, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1]
        """
        return '[%s]' % ', '.join('%s' % m for m in self)

    def _repr_exp_low(self):
        """
        Return a string representation of ``self`` in exponential form (lowest
        first).

        EXAMPLES::

            sage: print(Partition([7,7,7,3,3,2,1,1,1,1,1,1,1])._repr_exp_low())
            1^7, 2, 3^2, 7^3
            sage: print(Partition([])._repr_exp_low())
            -
        """
        if not self._list:
            return '-'
        exp = self.to_exp()
        return ', '.join('{}{}'.format(m, '' if e == 1 else '^%s' % e)
                         for m, e in enumerate(exp, start=1) if e)

    def _repr_exp_high(self):
        """
        Return a string representation of ``self`` in exponential form (highest
        first).

        EXAMPLES::

            sage: print(Partition([7,7,7,3,3,2,1,1,1,1,1,1,1])._repr_exp_high())
            7^3, 3^2, 2, 1^7

            sage: print(Partition([])._repr_exp_high())
            -
        """
        if not self._list:
            return '-'
        exp = self.to_exp()[::-1]         # reversed list of exponents
        M = max(self)
        return ', '.join('{}{}'.format(M - m, '' if e == 1 else '^%s' % e)
                         for m, e in enumerate(exp) if e)

    def _repr_compact_low(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with lowest first).

        EXAMPLES::

            sage: print(Partition([7,7,7,3,3,2,1,1,1,1,1,1,1])._repr_compact_low())
            1^7,2,3^2,7^3
            sage: print(Partition([])._repr_compact_low())
            -
        """
        if not self._list:
            return '-'
        exp = self.to_exp()
        return ','.join('{}{}'.format(m, '' if e == 1 else '^%s' % e)
                        for m, e in enumerate(exp, start=1) if e)

    def _repr_compact_high(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: print(Partition([7,7,7,3,3,2,1,1,1,1,1,1,1])._repr_compact_high())
            7^3,3^2,2,1^7
            sage: print(Partition([])._repr_compact_low())
            -
        """
        if not self._list:
            return '-'
        exp = self.to_exp()[::-1]         # reversed list of exponents
        M = max(self)
        return ','.join('{}{}'.format(M - m, '' if e == 1 else '^%s' % e)
                        for m, e in enumerate(exp) if e)

    def _repr_diagram(self):
        r"""
        Return a representation of ``self`` as a Ferrers diagram.

        EXAMPLES::

            sage: print(Partition([7,7,7,3,3,2,1,1,1,1,1,1,1])._repr_diagram())
            *******
            *******
            *******
            ***
            ***
            **
            *
            *
            *
            *
            *
            *
            *
        """
        return self.ferrers_diagram()

    def level(self) -> int:
        """
        Return the level of ``self``, which is always 1.

        This method exists only for compatibility with
        :class:`PartitionTuples`.

        EXAMPLES::

            sage: Partition([4,3,2]).level()
            1
        """
        return 1

    def components(self) -> list:
        """
        Return a list containing the shape of ``self``.

        This method exists only for compatibility with
        :class:`PartitionTuples`.

        EXAMPLES::

            sage: Partition([3,2]).components()
            [[3, 2]]
        """
        return [self]

    def _latex_(self) -> str:
        r"""
        Return a LaTeX version of ``self``.

        For more on the latex options, see :meth:`Partitions.options`.

        EXAMPLES::

            sage: mu = Partition([2, 1])
            sage: Partitions.options.latex='diagram'; latex(mu)       # indirect doctest
            {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\\
            \lr{\ast}&\lr{\ast}\\
            \lr{\ast}\\
            \end{array}$}
            }
            sage: Partitions.options.latex='exp_high'; latex(mu)      # indirect doctest
            2,1
            sage: Partitions.options.latex='exp_low'; latex(mu)       # indirect doctest
            1,2
            sage: Partitions.options.latex='list'; latex(mu)          # indirect doctest
            [2, 1]
            sage: Partitions.options.latex='young_diagram'; latex(mu) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }

            sage: Partitions.options(latex='young_diagram', convention='french')
            sage: Partitions.options.latex='exp_high'; latex(mu)      # indirect doctest
            2,1
            sage: Partitions.options.latex='exp_low'; latex(mu)       # indirect doctest
            1,2
            sage: Partitions.options.latex='list'; latex(mu)          # indirect doctest
            [2, 1]
            sage: Partitions.options.latex='young_diagram'; latex(mu) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[t]{*{2}c}\cline{1-1}
            \lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \end{array}$}
            }

            sage: Partitions.options._reset()
        """
        return self.parent().options._dispatch(self, '_latex_', 'latex')

    def _latex_young_diagram(self) -> str:
        r"""
        LaTeX output as a Young diagram.

        EXAMPLES::

            sage: print(Partition([2, 1])._latex_young_diagram())
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }
            sage: print(Partition([])._latex_young_diagram())
            {\emptyset}
        """
        if not self._list:
            return "{\\emptyset}"

        from sage.combinat.output import tex_from_array
        return tex_from_array([["\\phantom{x}"] * row_size
                               for row_size in self._list])

    def _latex_diagram(self) -> str:
        r"""
        LaTeX output as a Ferrers' diagram.

        EXAMPLES::

            sage: print(Partition([2, 1])._latex_diagram())
            {\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\\
            \lr{\ast}&\lr{\ast}\\
            \lr{\ast}\\
            \end{array}$}
            }
            sage: print(Partition([])._latex_diagram())
            {\emptyset}
        """
        if not self._list:
            return "{\\emptyset}"

        entry = self.parent().options("latex_diagram_str")
        from sage.combinat.output import tex_from_array
        return tex_from_array([[entry] * row_size
                               for row_size in self._list], False)

    def _latex_list(self) -> str:
        r"""
        LaTeX output as a list.

        EXAMPLES::

            sage: print(Partition([2, 1])._latex_list())
            [2, 1]
            sage: print(Partition([])._latex_list())
            []
        """
        return repr(self._list)

    def _latex_exp_low(self) -> str:
        r"""
        LaTeX output in exponential notation (lowest first).

        EXAMPLES::

            sage: print(Partition([2,2,1])._latex_exp_low())
            1,2^{2}
            sage: print(Partition([])._latex_exp_low())
            {\emptyset}
        """
        if not self._list:
            return "{\\emptyset}"
        exp = self.to_exp()
        return '%s' % ','.join('{}{}'.format(m + 1, '' if e == 1 else '^{%s}' % e)
                               for (m, e) in enumerate(exp) if e > 0)

    def _latex_exp_high(self):
        r"""
        LaTeX output in exponential notation (highest first).

        EXAMPLES::

            sage: print(Partition([2,2,1])._latex_exp_high())
            2^{2},1
            sage: print(Partition([])._latex_exp_high())
            {\emptyset}
        """
        if not self._list:
            return "{\\emptyset}"
        exp = self.to_exp()[::-1]  # reversed list of exponents
        M = max(self)
        return ','.join('{}{}'.format(M - m, '' if e == 1 else '^{%s}' % e)
                        for m, e in enumerate(exp) if e)

    def ferrers_diagram(self) -> str:
        r"""
        Return the Ferrers diagram of ``self``.

        EXAMPLES::

            sage: mu = Partition([5,5,2,1])
            sage: Partitions.options(diagram_str='*', convention='english')
            sage: print(mu.ferrers_diagram())
            *****
            *****
            **
            *
            sage: Partitions.options(diagram_str='▉')
            sage: print(mu.ferrers_diagram())
            ▉▉▉▉▉
            ▉▉▉▉▉
            ▉▉
            ▉
            sage: Partitions.options.convention="french"
            sage: print(mu.ferrers_diagram())
            ▉
            ▉▉
            ▉▉▉▉▉
            ▉▉▉▉▉
            sage: print(Partition([]).ferrers_diagram())
            -
            sage: Partitions.options(diagram_str='-')
            sage: print(Partition([]).ferrers_diagram())
            (/)
            sage: Partitions.options._reset()
        """
        diag_str = self.parent().options.diagram_str
        if not self._list:
            return '-' if diag_str != '-' else "(/)"
        if self.parent().options.convention == "English":
            return '\n'.join(diag_str * p for p in self)
        return '\n'.join(diag_str * p for p in reversed(self))

    def pp(self) -> None:
        r"""
        Print the Ferrers diagram.

        See :meth:`ferrers_diagram` for more on the Ferrers diagram.

        EXAMPLES::

            sage: Partition([5,5,2,1]).pp()
            *****
            *****
            **
            *
            sage: Partitions.options.convention='French'
            sage: Partition([5,5,2,1]).pp()
            *
            **
            *****
            *****
            sage: Partitions.options._reset()
        """
        print(self.ferrers_diagram())

    def __truediv__(self, p):
        """
        Return the skew partition ``self / p``.

        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p/[1,1]
            [3, 2, 1] / [1, 1]
            sage: p/[3,2,1]
            [3, 2, 1] / [3, 2, 1]
            sage: p/Partition([1,1])
            [3, 2, 1] / [1, 1]
            sage: p/[2,2,2]
            Traceback (most recent call last):
            ...
            ValueError: to form a skew partition p/q, q must be contained in p
        """
        if not self.contains(p):
            raise ValueError("to form a skew partition p/q, q must be contained in p")

        return SkewPartition([self[:], p])

    def stretch(self, k):
        """
        Return the partition obtained by multiplying each part with the
        given number.

        EXAMPLES::

            sage: p = Partition([4,2,2,1,1])
            sage: p.stretch(3)
            [12, 6, 6, 3, 3]
        """
        return _Partitions([k * p for p in self])

    def power(self, k):
        r"""
        Return the cycle type of the `k`-th power of any permutation
        with cycle type ``self`` (thus describes the powermap of
        symmetric groups).

        Equivalent to GAP's ``PowerPartition``.

        EXAMPLES::

            sage: p = Partition([5,3])
            sage: p.power(1)
            [5, 3]
            sage: p.power(2)
            [5, 3]
            sage: p.power(3)
            [5, 1, 1, 1]
            sage: p.power(4)
            [5, 3]

        Now let us compare this to the power map on `S_8`::

            sage: # needs sage.groups
            sage: G = SymmetricGroup(8)
            sage: g = G([(1,2,3,4,5),(6,7,8)]); g
            (1,2,3,4,5)(6,7,8)
            sage: g^2
            (1,3,5,2,4)(6,8,7)
            sage: g^3
            (1,4,2,5,3)
            sage: g^4
            (1,5,4,3,2)(6,7,8)

        ::

            sage: Partition([3,2,1]).power(3)
            [2, 1, 1, 1, 1]
        """
        res = []
        for i in self:
            g = int(gcd(i, k))
            res.extend([i // g] * g)
        res.sort(reverse=True)
        return Partition(res)

    def __next__(self):
        """
        Return the partition that lexicographically follows ``self``, of the
        same size. If ``self`` is the last partition, then return ``False``.

        EXAMPLES::

            sage: next(Partition([4]))
            [3, 1]
            sage: next(Partition([1,1,1,1]))
            False
        """
        p = self
        n = 0
        m = 0
        for i in p:
            n += i
            m += 1

        next_p = p[:] + [1]*(n - len(p))

        # Check to see if we are at the last (all ones) partition
        if p == [1] * n:
            return False

        #
        # If we are not, then run the ZS1 algorithm.
        #

        # Let h be the number of non-one  entries in the
        # partition
        h = 0
        for i in next_p:
            if i != 1:
                h += 1

        if next_p[h-1] == 2:
            m += 1
            next_p[h-1] = 1
            h -= 1
        else:
            r = next_p[h-1] - 1
            t = m - h + 1
            next_p[h-1] = r

            while t >= r:
                h += 1
                next_p[h-1] = r
                t -= r

            if t == 0:
                m = h
            else:
                m = h + 1
                if t > 1:
                    h += 1
                    next_p[h-1] = t

        return self.parent()(next_p[:m])

    next = __next__

    def size(self):
        """
        Return the size of ``self``.

        EXAMPLES::

            sage: Partition([2,2]).size()
            4
            sage: Partition([3,2,1]).size()
            6
        """
        return sum(self)

    def sign(self):
        r"""
        Return the sign of any permutation with cycle type ``self``.

        This function corresponds to a homomorphism from the symmetric
        group `S_n` into the cyclic group of order 2, whose kernel
        is exactly the alternating group `A_n`. Partitions of sign
        `1` are called even partitions while partitions of sign
        `-1` are called odd.

        EXAMPLES::

            sage: Partition([5,3]).sign()
            1
            sage: Partition([5,2]).sign()
            -1

        Zolotarev's lemma states that the Legendre symbol
        `\left(\frac{a}{p}\right)` for an integer
        `a \pmod p` (`p` a prime number), can be computed
        as sign(p_a), where sign denotes the sign of a permutation and
        p_a the permutation of the residue classes `\pmod p`
        induced by modular multiplication by `a`, provided
        `p` does not divide `a`.

        We verify this in some examples.

        ::

            sage: F = GF(11)                                                            # needs sage.rings.finite_rings
            sage: a = F.multiplicative_generator();a                                    # needs sage.rings.finite_rings
            2
            sage: plist = [int(a*F(x)) for x in range(1,11)]; plist                     # needs sage.rings.finite_rings
            [2, 4, 6, 8, 10, 1, 3, 5, 7, 9]

        This corresponds to the permutation (1, 2, 4, 8, 5, 10, 9, 7, 3, 6)
        (acting the set `\{1,2,...,10\}`) and to the partition
        [10].

        ::

            sage: p = PermutationGroupElement('(1, 2, 4, 8, 5, 10, 9, 7, 3, 6)')        # needs sage.groups
            sage: p.sign()                                                              # needs sage.groups
            -1
            sage: Partition([10]).sign()
            -1
            sage: kronecker_symbol(11,2)
            -1

        Now replace `2` by `3`::

            sage: plist = [int(F(3*x)) for x in range(1,11)]; plist                     # needs sage.rings.finite_rings
            [3, 6, 9, 1, 4, 7, 10, 2, 5, 8]
            sage: list(range(1, 11))
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: p = PermutationGroupElement('(3,4,8,7,9)')                            # needs sage.groups
            sage: p.sign()                                                              # needs sage.groups
            1
            sage: kronecker_symbol(3,11)
            1
            sage: Partition([5,1,1,1,1,1]).sign()
            1

        In both cases, Zolotarev holds.

        REFERENCES:

        - :wikipedia:`Zolotarev%27s_lemma`
        """
        return (-1)**(self.size() - self.length())

    def k_size(self, k):
        r"""
        Given a partition ``self`` and a ``k``, return the size of the
        `k`-boundary.

        This is the same as the length method
        :meth:`sage.combinat.core.Core.length` of the
        :class:`sage.combinat.core.Core` object, with the exception that here we
        don't require ``self`` to be a `k+1`-core.

        EXAMPLES::

            sage: Partition([2, 1, 1]).k_size(1)
            2
            sage: Partition([2, 1, 1]).k_size(2)
            3
            sage: Partition([2, 1, 1]).k_size(3)
            3
            sage: Partition([2, 1, 1]).k_size(4)
            4

        .. SEEALSO::

            :meth:`k_boundary`, :meth:`SkewPartition.size`
        """
        return self.k_boundary(k).size()

    def boundary(self):
        r"""
        Return the integer coordinates of points on the boundary of ``self``.

        For the following description, picture the Ferrer's diagram of ``self``
        using the French convention.  Recall that the French convention puts
        the longest row on the bottom and the shortest row on the top.  In
        addition, interpret the Ferrer's diagram as 1 x 1 cells in the Euclidean
        plane.  So if ``self`` was the partition [3, 1], the lower-left vertices
        of the 1 x 1 cells in the Ferrer's diagram would be (0, 0), (1, 0),
        (2, 0), and (0, 1).

        The boundary of a partition is the set `\{ \text{NE}(d) \mid \forall
        d\:\text{diagonal} \}`.  That is, for every diagonal line `y = x + b`
        where `b \in \mathbb{Z}`, we find the northeasternmost (NE) point on
        that diagonal which is also in the Ferrer's diagram.

        The boundary will go from bottom-right to top-left.

        EXAMPLES:

        Consider the partition (1) depicted as a square on a cartesian plane
        with vertices (0, 0), (1, 0), (1, 1), and (0, 1).  Three of those
        vertices in the appropriate order form the boundary::

            sage: Partition([1]).boundary()
            [(1, 0), (1, 1), (0, 1)]

        The partition (3, 1) can be visualized as three squares on a cartesian
        plane. The coordinates of the appropriate vertices form the boundary::

            sage: Partition([3, 1]).boundary()
            [(3, 0), (3, 1), (2, 1), (1, 1), (1, 2), (0, 2)]

        TESTS::

            sage: Partition([1]).boundary()
            [(1, 0), (1, 1), (0, 1)]
            sage: Partition([2, 1]).boundary()
            [(2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]
            sage: Partition([3, 1]).boundary()
            [(3, 0), (3, 1), (2, 1), (1, 1), (1, 2), (0, 2)]
            sage: Partition([2, 1, 1]).boundary()
            [(2, 0), (2, 1), (1, 1), (1, 2), (1, 3), (0, 3)]

        .. SEEALSO::

            :meth:`k_rim`.  You might have been looking for :meth:`k_boundary`
            instead.
        """
        def horizontal_piece(xy, bdy):
            (start_x, start_y) = xy
            if not bdy:
                h_piece = [(start_x, start_y)]
            else:
                stop_x = bdy[-1][0]
                y = start_y  # y never changes
                h_piece = [(x, y) for x in range(start_x, stop_x)]
            return list(reversed(h_piece))
        bdy = []
        for i, part in enumerate(self):
            (cell_x, cell_y) = (part - 1, i)
            (x, y) = (cell_x + 1, cell_y + 1)
            bdy += horizontal_piece((x, y - 1), bdy)
            bdy.append((x, y))
        # add final "top-left" horizontal piece
        (top_left_x, top_left_y) = (0, len(self))
        bdy += horizontal_piece((top_left_x, top_left_y), bdy)
        return bdy

    def k_rim(self, k):
        r"""
        Return the ``k``-rim of ``self`` as a list of integer coordinates.

        The `k`-rim of a partition is the "line between" (or "intersection of")
        the `k`-boundary and the `k`-interior.  (Section 2.3 of [HM2011]_)

        It will be output as an ordered list of integer coordinates, where the
        origin is `(0, 0)`.  It will start at the top-left of the `k`-rim (using
        French convention) and end at the bottom-right.

        EXAMPLES:

        Consider the partition (3, 1) split up into its 1-interior and
        1-boundary:

        .. image:: ../../media/k-rim.JPG
            :height: 180px
            :align: center

        The line shown in bold is the 1-rim, and that information is equivalent
        to the integer coordinates of the points that occur along that line::

            sage: Partition([3, 1]).k_rim(1)
            [(3, 0), (2, 0), (2, 1), (1, 1), (0, 1), (0, 2)]

        TESTS::

            sage: Partition([1]).k_rim(0)
            [(1, 0),  (1, 1),  (0, 1)]
            sage: Partition([3,  1]).k_rim(0)
            [(3, 0),  (3, 1),  (2, 1),  (1, 1),  (1, 2),  (0, 2)]
            sage: Partition([3,  1]).k_rim(1)
            [(3, 0),  (2, 0),  (2, 1),  (1, 1),  (0, 1),  (0, 2)]
            sage: Partition([3,  1]).k_rim(2)
            [(3, 0),  (2, 0),  (1, 0),  (1, 1),  (0, 1),  (0, 2)]
            sage: Partition([3,  1]).k_rim(3)
            [(3, 0),  (2, 0),  (1, 0),  (1, 1),  (0, 1),  (0, 2)]

        .. SEEALSO::

            :meth:`k_interior`, :meth:`k_boundary`, :meth:`boundary`
        """
        interior_rim = self.k_interior(k).boundary()
        # get leftmost vertical line
        interior_top_left_y = interior_rim[-1][1]
        v_piece = [(0, y) for y in range(interior_top_left_y+1, len(self)+1)]
        # get bottommost horizontal line
        interior_bottom_right_x = interior_rim[0][0]
        if self:
            ptn_bottom_right_x = self[0]
        else:
            ptn_bottom_right_x = 0
        h_piece = [(x, 0) for x in
                   range(ptn_bottom_right_x, interior_bottom_right_x, -1)]
        # glue together with boundary
        rim = h_piece + interior_rim + v_piece
        return rim

    def k_row_lengths(self, k):
        r"""
        Return the ``k``-row-shape of the partition ``self``.

        This is equivalent to taking the `k`-boundary of the partition and then
        returning the row-shape of that.  We do *not* discard rows of length 0.
        (Section 2.2 of [LLMS2013]_)

        EXAMPLES::

            sage: Partition([6, 1]).k_row_lengths(2)
            [2, 1]

            sage: Partition([4, 4, 4, 3, 2]).k_row_lengths(2)
            [0, 1, 1, 1, 2]

        .. SEEALSO::

            :meth:`k_column_lengths`, :meth:`k_boundary`,
            :meth:`SkewPartition.row_lengths`,
            :meth:`SkewPartition.column_lengths`
        """
        return self.k_boundary(k).row_lengths()

    def k_column_lengths(self, k):
        r"""
        Return the ``k``-column-shape of the partition ``self``.

        This is the 'column' analog of :meth:`k_row_lengths`.

        EXAMPLES::

            sage: Partition([6, 1]).k_column_lengths(2)
            [1, 0, 0, 0, 1, 1]

            sage: Partition([4, 4, 4, 3, 2]).k_column_lengths(2)
            [1, 1, 1, 2]

        .. SEEALSO::

            :meth:`k_row_lengths`, :meth:`k_boundary`,
            :meth:`SkewPartition.row_lengths`,
            :meth:`SkewPartition.column_lengths`
        """
        return self.k_boundary(k).column_lengths()

    def has_rectangle(self, h, w) -> bool:
        r"""
        Return ``True`` if the Ferrer's diagram of ``self`` has ``h``
        (*or more*) rows of length ``w`` (*exactly*).

        INPUT:

        - ``h`` -- integer `h \geq 1`;  the (*minimum*) height of the
          rectangle

        - ``w`` -- integer `w \geq 1`;  the width of the rectangle

        EXAMPLES::

            sage: Partition([3, 3, 3, 3]).has_rectangle(2, 3)
            True
            sage: Partition([3, 3]).has_rectangle(2, 3)
            True
            sage: Partition([4, 3]).has_rectangle(2, 3)
            False
            sage: Partition([3]).has_rectangle(2, 3)
            False

        TESTS::

            sage: Partition([1, 1, 1]).has_rectangle(4, 1)
            False
            sage: Partition([1, 1, 1]).has_rectangle(3, 1)
            True
            sage: Partition([1, 1, 1]).has_rectangle(2, 1)
            True
            sage: Partition([1, 1, 1]).has_rectangle(1, 2)
            False
            sage: Partition([3]).has_rectangle(1, 3)
            True
            sage: Partition([3]).has_rectangle(1, 2)
            False
            sage: Partition([3]).has_rectangle(2, 3)
            False

        .. SEEALSO::

            :meth:`has_k_rectangle`
        """
        assert h >= 1
        assert w >= 1
        return self.to_exp(w)[w - 1] >= h

    def has_k_rectangle(self, k) -> bool:
        r"""
        Return ``True`` if the Ferrer's diagram of ``self`` contains `k-i+1`
        rows (*or more*) of length `i` (*exactly*) for any `i` in `[1, k]`.

        This is mainly a helper function for :meth:`is_k_reducible` and
        :meth:`is_k_irreducible`, the only difference between this function and
        :meth:`is_k_reducible` being that this function allows any partition as
        input while :meth:`is_k_reducible` requires the input to be `k`-bounded.

        EXAMPLES:

        The partition [1, 1, 1] has at least 2 rows of length 1::

            sage: Partition([1, 1, 1]).has_k_rectangle(2)
            True

        The partition [1, 1, 1] does *not* have 4 rows of length 1, 3 rows of
        length 2, 2 rows of length 3, nor 1 row of length 4::

            sage: Partition([1, 1, 1]).has_k_rectangle(4)
            False

        TESTS::

            sage: Partition([1]).has_k_rectangle(1)
            True
            sage: Partition([1]).has_k_rectangle(2)
            False
            sage: Partition([1, 1, 1]).has_k_rectangle(3)
            True
            sage: Partition([1, 1, 1]).has_k_rectangle(2)
            True
            sage: Partition([1, 1, 1]).has_k_rectangle(4)
            False
            sage: Partition([3]).has_k_rectangle(3)
            True
            sage: Partition([3]).has_k_rectangle(2)
            False
            sage: Partition([3]).has_k_rectangle(4)
            False

        .. SEEALSO::

            :meth:`is_k_irreducible`, :meth:`is_k_reducible`,
            :meth:`has_rectangle`
        """
        return any(self.has_rectangle(k - i + 1, i)
                   for i in range(1, k + 1))

    def is_k_bounded(self, k) -> bool:
        r"""
        Return ``True`` if the partition ``self`` is bounded by ``k``.

        EXAMPLES::

            sage: Partition([4, 3, 1]).is_k_bounded(4)
            True
            sage: Partition([4, 3, 1]).is_k_bounded(7)
            True
            sage: Partition([4, 3, 1]).is_k_bounded(3)
            False
        """
        assert k >= 0
        if self.is_empty():
            return True
        else:
            return self[0] <= k

    def is_k_reducible(self, k):
        r"""
        Return ``True`` if the partition ``self`` is ``k``-reducible.

        A `k`-bounded partition is `k`-*reducible* if its Ferrer's diagram
        contains `k-i+1` rows (or more) of length `i` (exactly) for some
        `i \in [1, k]`.

        (Also, a `k`-bounded partition is `k`-reducible if and only if it is not `k`-irreducible.)

        EXAMPLES:

        The partition [1, 1, 1] has at least 2 rows of length 1::

            sage: Partition([1, 1, 1]).is_k_reducible(2)
            True

        The partition [1, 1, 1] does *not* have 4 rows of length 1, 3 rows of
        length 2, 2 rows of length 3, nor 1 row of length 4::

            sage: Partition([1, 1, 1]).is_k_reducible(4)
            False

        .. SEEALSO::

            :meth:`is_k_irreducible`, :meth:`has_k_rectangle`
        """
        if not self.is_k_bounded(k):
            raise ValueError('we only talk about k-reducible / k-irreducible for k-bounded partitions')
        return self.has_k_rectangle(k)

    def is_k_irreducible(self, k):
        r"""
        Return ``True`` if the partition ``self`` is ``k``-irreducible.

        A `k`-bounded partition is `k`-*irreducible* if its Ferrer's diagram
        does *not* contain `k-i+1` rows (or more) of length `i` (exactly) for
        every `i \in [1, k]`.

        (Also, a `k`-bounded partition is `k`-irreducible if and only if it is
        not `k`-reducible.)

        EXAMPLES:

        The partition [1, 1, 1] has at least 2 rows of length 1::

            sage: Partition([1, 1, 1]).is_k_irreducible(2)
            False

        The partition [1, 1, 1] does *not* have 4 rows of length 1, 3 rows of
        length 2, 2 rows of length 3, nor 1 row of length 4::

            sage: Partition([1, 1, 1]).is_k_irreducible(4)
            True

        .. SEEALSO::

            :meth:`is_k_reducible`, :meth:`has_k_rectangle`
        """
        return not self.is_k_reducible(k)

    def is_symmetric(self):
        r"""
        Return ``True`` if the partition ``self`` equals its own transpose.

        EXAMPLES::

            sage: Partition([2, 1]).is_symmetric()
            True
            sage: Partition([3, 1]).is_symmetric()
            False
        """
        return self == self.conjugate()

    def next_within_bounds(self, min=[], max=None, partition_type=None):
        r"""
        Get the next partition lexicographically that contains ``min`` and is
        contained in ``max``.

        INPUT:

        - ``min`` -- (default: ``[]``, the empty partition) the
          'minimum partition' that ``next_within_bounds(self)`` must contain

        - ``max`` -- (default: ``None``) the 'maximum partition' that
          ``next_within_bounds(self)`` must be contained in;  if set to ``None``,
          then there is no restriction

        - ``partition_type`` -- (default: ``None``) the type of partitions
          allowed;  for example, 'strict' for strictly decreasing partitions, or
          ``None`` to allow any valid partition

        EXAMPLES::

            sage: m = [1, 1]
            sage: M = [3, 2, 1]
            sage: Partition([1, 1]).next_within_bounds(min=m, max=M)
            [1, 1, 1]
            sage: Partition([1, 1, 1]).next_within_bounds(min=m, max=M)
            [2, 1]
            sage: Partition([2, 1]).next_within_bounds(min=m, max=M)
            [2, 1, 1]
            sage: Partition([2, 1, 1]).next_within_bounds(min=m, max=M)
            [2, 2]
            sage: Partition([2, 2]).next_within_bounds(min=m, max=M)
            [2, 2, 1]
            sage: Partition([2, 2, 1]).next_within_bounds(min=m, max=M)
            [3, 1]
            sage: Partition([3, 1]).next_within_bounds(min=m, max=M)
            [3, 1, 1]
            sage: Partition([3, 1, 1]).next_within_bounds(min=m, max=M)
            [3, 2]
            sage: Partition([3, 2]).next_within_bounds(min=m, max=M)
            [3, 2, 1]
            sage: Partition([3, 2, 1]).next_within_bounds(min=m, max=M) == None
            True

        .. SEEALSO::

            :meth:`next`
        """
        # make sure min <= self <= max
        if max is not None:
            assert _Partitions(max).contains(_Partitions(self))
        assert _Partitions(self).contains(_Partitions(min))
        # check for empty max
        if max is not None and _Partitions(max).is_empty():
            return None
        # convert partitions to lists to make them mutable
        p = list(self)
        min = list(min)
        # if there is no max, the next partition just tacks a '1' on to the end!
        if max is None:
            return _Partitions(p + [1])
        # extend p and min to include 0s at the end
        p = p + [0] * (len(max) - len(p))
        min = min + [0] * (len(max) - len(min))
        # finally, run the algo to find next_p
        next_p = copy(p)

        def condition(a, b):
            if partition_type in ('strict', 'strictly decreasing'):
                return a < b - 1
            if partition_type in (None, 'weak', 'weakly decreasing'):
                return a < b
            raise ValueError('unrecognized partition type')

        for r in range(len(p) - 1, -1, -1):
            if not r:
                if max is None or p[r] < max[r]:
                    next_p[r] += 1
                    break
                return None
            elif (max is None or p[r] < max[r]) and condition(p[r], p[r-1]):
                next_p[r] += 1
                break
            next_p[r] = min[r]
            continue

        return _Partitions(next_p)

    def row_standard_tableaux(self):
        """
        Return the :class:`row standard tableaux
        <sage.combinat.tableau.RowStandardTableaux>` of shape ``self``.

        EXAMPLES::

            sage: Partition([3,2,2,1]).row_standard_tableaux()
            Row standard tableaux of shape [3, 2, 2, 1]
        """
        return tableau.RowStandardTableaux(self)

    def standard_tableaux(self):
        """
        Return the :class:`standard tableaux<StandardTableaux>`
        of shape ``self``.

        EXAMPLES::

            sage: Partition([3,2,2,1]).standard_tableaux()
            Standard tableaux of shape [3, 2, 2, 1]
        """
        return tableau.StandardTableaux(self)

    def up(self):
        r"""
        Return a generator for partitions that can be obtained from ``self``
        by adding a cell.

        EXAMPLES::

            sage: list(Partition([2,1,1]).up())
            [[3, 1, 1], [2, 2, 1], [2, 1, 1, 1]]
            sage: list(Partition([3,2]).up())
            [[4, 2], [3, 3], [3, 2, 1]]
            sage: [p for p in Partition([]).up()]
            [[1]]
        """
        p = self
        previous = p.get_part(0) + 1
        for i, current in enumerate(p):
            if current < previous:
                yield Partition(p[:i] + [current + 1] + p[i + 1:])
            previous = current
        yield Partition(p + [1])

    def up_list(self):
        """
        Return a list of the partitions that can be formed from ``self`` by
        adding a cell.

        EXAMPLES::

            sage: Partition([2,1,1]).up_list()
            [[3, 1, 1], [2, 2, 1], [2, 1, 1, 1]]
            sage: Partition([3,2]).up_list()
            [[4, 2], [3, 3], [3, 2, 1]]
            sage: Partition([]).up_list()
            [[1]]
        """
        return list(self.up())

    def down(self):
        r"""
        Return a generator for partitions that can be obtained from ``self``
        by removing a cell.

        EXAMPLES::

            sage: [p for p in Partition([2,1,1]).down()]
            [[1, 1, 1], [2, 1]]
            sage: [p for p in Partition([3,2]).down()]
            [[2, 2], [3, 1]]
            sage: [p for p in Partition([3,2,1]).down()]
            [[2, 2, 1], [3, 1, 1], [3, 2]]

        TESTS:

        We check that :issue:`11435` is fixed::

            sage: Partition([]).down_list() #indirect doctest
            []
        """
        p = self
        l = len(p)
        for i in range(l-1):
            if p[i] > p[i+1]:
                yield Partition(p[:i] + [p[i]-1] + p[i+1:])
        if l >= 1:
            last = p[-1]
            if last == 1:
                yield Partition(p[:-1])
            else:
                yield Partition(p[:-1] + [p[-1] - 1])

    def down_list(self):
        """
        Return a list of the partitions that can be obtained from ``self``
        by removing a cell.

        EXAMPLES::

            sage: Partition([2,1,1]).down_list()
            [[1, 1, 1], [2, 1]]
            sage: Partition([3,2]).down_list()
            [[2, 2], [3, 1]]
            sage: Partition([3,2,1]).down_list()
            [[2, 2, 1], [3, 1, 1], [3, 2]]
            sage: Partition([]).down_list()  #checks :issue:`11435`
            []
        """
        return list(self.down())

    @combinatorial_map(name="cell poset")
    def cell_poset(self, orientation='SE'):
        """
        Return the Young diagram of ``self`` as a poset. The optional
        keyword variable ``orientation`` determines the order relation
        of the poset.

        The poset always uses the set of cells of the Young diagram
        of ``self`` as its ground set. The order relation of the poset
        depends on the ``orientation`` variable (which defaults to
        ``'SE'``). Concretely, ``orientation`` has to be specified to
        one of the strings ``'NW'``, ``'NE'``, ``'SW'``, and ``'SE'``,
        standing for "northwest", "northeast", "southwest" and
        "southeast", respectively. If ``orientation`` is ``'SE'``, then
        the order relation of the poset is such that a cell `u` is
        greater or equal to a cell `v` in the poset if and only if `u`
        lies weakly southeast of `v` (this means that `u` can be
        reached from `v` by a sequence of south and east steps; the
        sequence is allowed to consist of south steps only, or of east
        steps only, or even be empty). Similarly the order relation is
        defined for the other three orientations. The Young diagram is
        supposed to be drawn in English notation.

        The elements of the poset are the cells of the Young diagram
        of ``self``, written as tuples of zero-based coordinates (so
        that `(3, 7)` stands for the `8`-th cell of the `4`-th row,
        etc.).

        EXAMPLES::

            sage: # needs sage.graphs
            sage: p = Partition([3,3,1])
            sage: Q = p.cell_poset(); Q
            Finite poset containing 7 elements
            sage: sorted(Q)
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0)]
            sage: sorted(Q.maximal_elements())
            [(1, 2), (2, 0)]
            sage: Q.minimal_elements()
            [(0, 0)]
            sage: sorted(Q.upper_covers((1, 0)))
            [(1, 1), (2, 0)]
            sage: Q.upper_covers((1, 1))
            [(1, 2)]

            sage: # needs sage.graphs
            sage: P = p.cell_poset(orientation="NW"); P
            Finite poset containing 7 elements
            sage: sorted(P)
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0)]
            sage: sorted(P.minimal_elements())
            [(1, 2), (2, 0)]
            sage: P.maximal_elements()
            [(0, 0)]
            sage: P.upper_covers((2, 0))
            [(1, 0)]
            sage: sorted(P.upper_covers((1, 2)))
            [(0, 2), (1, 1)]
            sage: sorted(P.upper_covers((1, 1)))
            [(0, 1), (1, 0)]
            sage: sorted([len(P.upper_covers(v)) for v in P])
            [0, 1, 1, 1, 1, 2, 2]

            sage: # needs sage.graphs
            sage: R = p.cell_poset(orientation="NE"); R
            Finite poset containing 7 elements
            sage: sorted(R)
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0)]
            sage: R.maximal_elements()
            [(0, 2)]
            sage: R.minimal_elements()
            [(2, 0)]
            sage: sorted([len(R.upper_covers(v)) for v in R])
            [0, 1, 1, 1, 1, 2, 2]
            sage: R.is_isomorphic(P)
            False
            sage: R.is_isomorphic(P.dual())
            False

        Linear extensions of ``p.cell_poset()`` are in 1-to-1 correspondence
        with standard Young tableaux of shape `p`::

            sage: all( len(p.cell_poset().linear_extensions())                          # needs sage.graphs
            ....:      == len(p.standard_tableaux())
            ....:      for n in range(8) for p in Partitions(n) )
            True

        This is not the case for northeast orientation::

            sage: q = Partition([3, 1])
            sage: q.cell_poset(orientation="NE").is_chain()                             # needs sage.graphs
            True

        TESTS:

        We check that the posets are really what they should be for size
        up to `7`::

            sage: def check_NW(n):
            ....:     for p in Partitions(n):
            ....:         P = p.cell_poset(orientation="NW")
            ....:         for c in p.cells():
            ....:             for d in p.cells():
            ....:                 if P.le(c, d) != (c[0] >= d[0]
            ....:                                   and c[1] >= d[1]):
            ....:                     return False
            ....:     return True
            sage: all( check_NW(n) for n in range(8) )                                  # needs sage.graphs
            True

            sage: def check_NE(n):
            ....:     for p in Partitions(n):
            ....:         P = p.cell_poset(orientation="NE")
            ....:         for c in p.cells():
            ....:             for d in p.cells():
            ....:                 if P.le(c, d) != (c[0] >= d[0]
            ....:                                   and c[1] <= d[1]):
            ....:                     return False
            ....:     return True
            sage: all( check_NE(n) for n in range(8) )                                  # needs sage.graphs
            True

            sage: def test_duality(n, ori1, ori2):
            ....:     for p in Partitions(n):
            ....:         P = p.cell_poset(orientation=ori1)
            ....:         Q = p.cell_poset(orientation=ori2)
            ....:         for c in p.cells():
            ....:             for d in p.cells():
            ....:                 if P.lt(c, d) != Q.lt(d, c):
            ....:                     return False
            ....:     return True
            sage: all( test_duality(n, "NW", "SE") for n in range(8) )                  # needs sage.graphs
            True
            sage: all( test_duality(n, "NE", "SW") for n in range(8) )                  # needs sage.graphs
            True
            sage: all( test_duality(n, "NE", "SE") for n in range(4) )                  # needs sage.graphs
            False
        """
        from sage.combinat.posets.posets import Poset
        covers = {}
        if orientation == "NW":
            for i, row in enumerate(self):
                if i == 0:
                    covers[(0, 0)] = []
                    for j in range(1, row):
                        covers[(0, j)] = [(0, j - 1)]
                else:
                    covers[(i, 0)] = [(i - 1, 0)]
                    for j in range(1, row):
                        covers[(i, j)] = [(i - 1, j), (i, j - 1)]
        elif orientation == "NE":
            for i, row in enumerate(self):
                if i == 0:
                    covers[(0, row - 1)] = []
                    for j in range(row - 1):
                        covers[(0, j)] = [(0, j + 1)]
                else:
                    covers[(i, row - 1)] = [(i - 1, row - 1)]
                    for j in range(row - 1):
                        covers[(i, j)] = [(i - 1, j), (i, j + 1)]
        elif orientation == "SE":
            l = len(self) - 1
            for i, row in enumerate(self):
                if i == l:
                    covers[(i, row - 1)] = []
                    for j in range(row - 1):
                        covers[(i, j)] = [(i, j + 1)]
                else:
                    next_row = self[i + 1]
                    if row == next_row:
                        covers[(i, row - 1)] = [(i + 1, row - 1)]
                        for j in range(row - 1):
                            covers[(i, j)] = [(i + 1, j), (i, j + 1)]
                    else:
                        covers[(i, row - 1)] = []
                        for j in range(next_row):
                            covers[(i, j)] = [(i + 1, j), (i, j + 1)]
                        for j in range(next_row, row - 1):
                            covers[(i, j)] = [(i, j + 1)]
        elif orientation == "SW":
            l = len(self) - 1
            for i, row in enumerate(self):
                if i == l:
                    covers[(i, 0)] = []
                    for j in range(1, row):
                        covers[(i, j)] = [(i, j - 1)]
                else:
                    covers[(i, 0)] = [(i + 1, 0)]
                    next_row = self[i + 1]
                    for j in range(1, next_row):
                        covers[(i, j)] = [(i + 1, j), (i, j - 1)]
                    for j in range(next_row, row):
                        covers[(i, j)] = [(i, j - 1)]
        return Poset(covers)

    def frobenius_coordinates(self):
        """
        Return a pair of sequences of Frobenius coordinates aka beta numbers
        of the partition.

        These are two strictly decreasing sequences of nonnegative integers
        of the same length.

        EXAMPLES::

            sage: Partition([]).frobenius_coordinates()
            ([], [])
            sage: Partition([1]).frobenius_coordinates()
            ([0], [0])
            sage: Partition([3,3,3]).frobenius_coordinates()
            ([2, 1, 0], [2, 1, 0])
            sage: Partition([9,1,1,1,1,1,1]).frobenius_coordinates()
            ([8], [6])
        """
        mu = self
        muconj = mu.conjugate()     # Naive implementation
        if len(mu) <= len(muconj):
            a = [x for x in (val-i-1 for i, val in enumerate(mu)) if x >= 0]
            b = [x for x in (muconj[i]-i-1 for i in range(len(a))) if x >= 0]
        else:
            b = [x for x in (val-i-1 for i, val in enumerate(muconj)) if x >= 0]
            a = [x for x in (mu[i]-i-1 for i in range(len(b))) if x >= 0]
        return (a, b)

    def frobenius_rank(self):
        r"""
        Return the Frobenius rank of the partition ``self``.

        The Frobenius rank of a partition
        `\lambda = (\lambda_1, \lambda_2, \lambda_3, \cdots)` is
        defined to be the largest `i` such that `\lambda_i \geq i`.
        In other words, it is the number of cells on the main diagonal
        of `\lambda`. In yet other words, it is the size of the largest
        square fitting into the Young diagram of `\lambda`.

        EXAMPLES::

            sage: Partition([]).frobenius_rank()
            0
            sage: Partition([1]).frobenius_rank()
            1
            sage: Partition([3,3,3]).frobenius_rank()
            3
            sage: Partition([9,1,1,1,1,1]).frobenius_rank()
            1
            sage: Partition([2,1,1,1,1,1]).frobenius_rank()
            1
            sage: Partition([2,2,1,1,1,1]).frobenius_rank()
            2
            sage: Partition([3,2]).frobenius_rank()
            2
            sage: Partition([3,2,2]).frobenius_rank()
            2
            sage: Partition([8,4,4,4,4]).frobenius_rank()
            4
            sage: Partition([8,4,1]).frobenius_rank()
            2
            sage: Partition([3,3,1]).frobenius_rank()
            2
        """
        for i, x in enumerate(self):
            if x <= i:
                return i
        return len(self)

    def beta_numbers(self, length=None):
        """
        Return the set of beta numbers corresponding to ``self``.

        The optional argument ``length`` specifies the length of the beta set
        (which must be at least the length of ``self``).

        For more on beta numbers, see :meth:`frobenius_coordinates`.

        EXAMPLES::

            sage: Partition([4,3,2]).beta_numbers()
            [6, 4, 2]
            sage: Partition([4,3,2]).beta_numbers(5)
            [8, 6, 4, 1, 0]
            sage: Partition([]).beta_numbers()
            []
            sage: Partition([]).beta_numbers(3)
            [2, 1, 0]
            sage: Partition([6,4,1,1]).beta_numbers()
            [9, 6, 2, 1]
            sage: Partition([6,4,1,1]).beta_numbers(6)
            [11, 8, 4, 3, 1, 0]
            sage: Partition([1,1,1]).beta_numbers()
            [3, 2, 1]
            sage: Partition([1,1,1]).beta_numbers(4)
            [4, 3, 2, 0]
        """
        true_length = len(self)
        if length is None:
            length = true_length
        elif length < true_length:
            raise ValueError("length must be at least the length of the partition")
        beta = [l + length - i for i, l in enumerate(self, start=1)]
        if length > true_length:
            beta.extend(range(length - true_length - 1, -1, -1))
        return beta

    def crank(self):
        r"""
        Return the Dyson crank of ``self``.

        The Dyson crank of a partition `\lambda` is defined as follows:
        If `\lambda` contains at least one `1`, then the crank is
        `\mu(\lambda) - \omega(\lambda)`, where `\omega(\lambda)` is the
        number of `1`s in `\lambda`, and `\mu(\lambda)` is the number of
        parts of `\lambda` larger than `\omega(\lambda)`. If `\lambda`
        contains no `1`, then the crank is simply the largest part of
        `\lambda`.

        REFERENCES:

        - [AG1988]_

        EXAMPLES::

            sage: Partition([]).crank()
            0
            sage: Partition([3,2,2]).crank()
            3
            sage: Partition([5,4,2,1,1]).crank()
            0
            sage: Partition([1,1,1]).crank()
            -3
            sage: Partition([6,4,4,3]).crank()
            6
            sage: Partition([6,3,3,1,1]).crank()
            1
            sage: Partition([6]).crank()
            6
            sage: Partition([5,1]).crank()
            0
            sage: Partition([4,2]).crank()
            4
            sage: Partition([4,1,1]).crank()
            -1
            sage: Partition([3,3]).crank()
            3
            sage: Partition([3,2,1]).crank()
            1
            sage: Partition([3,1,1,1]).crank()
            -3
            sage: Partition([2,2,2]).crank()
            2
            sage: Partition([2,2,1,1]).crank()
            -2
            sage: Partition([2,1,1,1,1]).crank()
            -4
            sage: Partition([1,1,1,1,1,1]).crank()
            -6
        """
        l = len(self)
        if l == 0:
            return 0
        if self[-1] > 1:
            return self[0]
        ind_1 = self.index(1)
        w = l - ind_1      # w is omega(self).
        m = len([x for x in self if x > w])
        return m - w

    def t_completion(self, t):
        r"""
        Return the ``t``-completion of the partition ``self``.

        If `\lambda = (\lambda_1, \lambda_2, \lambda_3, \ldots)` is a
        partition and `t` is an integer greater or equal to
        `\left\lvert \lambda \right\rvert + \lambda_1`, then the
        `t`-*completion of* `\lambda` is defined as the partition
        `(t - \left\lvert \lambda \right\rvert, \lambda_1, \lambda_2,
        \lambda_3, \ldots)` of `t`. This partition is denoted by `\lambda[t]`
        in [BOR2009]_, by `\lambda_{[t]}` in [BdVO2012]_, and by `\lambda(t)`
        in [CO2010]_.

        EXAMPLES::

            sage: Partition([]).t_completion(0)
            []
            sage: Partition([]).t_completion(1)
            [1]
            sage: Partition([]).t_completion(2)
            [2]
            sage: Partition([]).t_completion(3)
            [3]
            sage: Partition([2, 1]).t_completion(5)
            [2, 2, 1]
            sage: Partition([2, 1]).t_completion(6)
            [3, 2, 1]
            sage: Partition([4, 2, 2, 1]).t_completion(13)
            [4, 4, 2, 2, 1]
            sage: Partition([4, 2, 2, 1]).t_completion(19)
            [10, 4, 2, 2, 1]
            sage: Partition([4, 2, 2, 1]).t_completion(10)
            Traceback (most recent call last):
            ...
            ValueError: 10-completion is not defined
            sage: Partition([4, 2, 2, 1]).t_completion(5)
            Traceback (most recent call last):
            ...
            ValueError: 5-completion is not defined
        """
        if self._list and t < self.size() + self._list[0]:
            raise ValueError(f"{t}-completion is not defined")
        return Partition([t - self.size()] + self._list)

    def larger_lex(self, rhs):
        """
        Return ``True`` if ``self`` is larger than ``rhs`` in lexicographic
        order. Otherwise return ``False``.

        EXAMPLES::

            sage: p = Partition([3,2])
            sage: p.larger_lex([3,1])
            True
            sage: p.larger_lex([1,4])
            True
            sage: p.larger_lex([3,2,1])
            False
            sage: p.larger_lex([3])
            True
            sage: p.larger_lex([5])
            False
            sage: p.larger_lex([3,1,1,1,1,1,1,1])
            True
        """
        return CombinatorialElement.__gt__(self, rhs)

    def dominates(self, p2):
        r"""
        Return ``True`` if ``self`` dominates the partition ``p2``. Otherwise
        it returns ``False``.

        EXAMPLES::

            sage: p = Partition([3,2])
            sage: p.dominates([3,1])
            True
            sage: p.dominates([2,2])
            True
            sage: p.dominates([2,1,1])
            True
            sage: p.dominates([3,3])
            False
            sage: p.dominates([4])
            False
            sage: Partition([4]).dominates(p)
            False
            sage: Partition([]).dominates([1])
            False
            sage: Partition([]).dominates([])
            True
            sage: Partition([1]).dominates([])
            True
        """
        p1 = self
        sum1 = 0
        sum2 = 0
        min_length = min(len(p1), len(p2))
        if min_length == 0:
            return not p2  # equivalent to len(p1) >= len(p2) = 0

        for i in range(min_length):
            sum1 += p1[i]
            sum2 += p2[i]
            if sum2 > sum1:
                return False
        return sum(p1) >= sum(p2)

    def cells(self):
        """
        Return the coordinates of the cells of ``self``.

        EXAMPLES::

            sage: Partition([2,2]).cells()
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: Partition([3,2]).cells()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
        """
        return [(i, j) for i, si in enumerate(self) for j in range(si)]

    def generalized_pochhammer_symbol(self, a, alpha):
        r"""
        Return the generalized Pochhammer symbol
        `(a)_{self}^{(\alpha)}`. This is the product over all
        cells `(i,j)` in ``self`` of `a - (i-1) / \alpha + j - 1`.

        EXAMPLES::

            sage: Partition([2,2]).generalized_pochhammer_symbol(2,1)
            12
        """
        res = 1
        for (i, j) in self.cells():
            res *= (a - (i-1)/alpha + j-1)
        return res

    def get_part(self, i, default=Integer(0)):
        r"""
        Return the `i`-th part of ``self``, or ``default`` if it does
        not exist.

        EXAMPLES::

            sage: p = Partition([2,1])
            sage: p.get_part(0), p.get_part(1), p.get_part(2)
            (2, 1, 0)
            sage: p.get_part(10,-1)
            -1
            sage: Partition([]).get_part(0)
            0
        """
        if i < len(self._list):
            return self._list[i]
        else:
            return default

    @combinatorial_map(name="partition to minimal Dyck word")
    def to_dyck_word(self, n=None):
        r"""
        Return the ``n``-Dyck word whose corresponding partition is
        ``self`` (or, if ``n`` is not specified, the `n`-Dyck word with
        smallest `n` to satisfy this property).

        If `w` is an `n`-Dyck word (that is, a Dyck word with `n` open
        symbols and `n` close symbols), then the Dyck path corresponding
        to `w` can be regarded as a lattice path in the northeastern
        half of an `n \times n`-square. The region to the northeast of
        this Dyck path can be regarded as a partition. It is called the
        partition corresponding to the Dyck word `w`. (See
        :meth:`~sage.combinat.dyck_word.DyckWord.to_partition`.)

        For every partition `\lambda` and every nonnegative integer `n`,
        there exists at most one `n`-Dyck word `w` such that the
        partition corresponding to `w` is `\lambda` (in fact, such `w`
        exists if and only if `\lambda_i + i \leq n` for every `i`,
        where `\lambda` is written in the form
        `(\lambda_1, \lambda_2, \ldots, \lambda_k)` with `\lambda_k > 0`).
        This method computes this `w` for a given `\lambda` and `n`.
        If `n` is not specified, this method computes the `w` for the
        smallest possible `n` for which such an `w` exists.
        (The minimality of `n` means that the partition demarcated by the
        Dyck path touches the diagonal.)

        EXAMPLES::

            sage: Partition([2,2]).to_dyck_word()
            [1, 1, 0, 0, 1, 1, 0, 0]
            sage: Partition([2,2]).to_dyck_word(4)
            [1, 1, 0, 0, 1, 1, 0, 0]
            sage: Partition([2,2]).to_dyck_word(5)
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 0]
            sage: Partition([6,3,1]).to_dyck_word()
            [1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0]
            sage: Partition([]).to_dyck_word()
            []
            sage: Partition([]).to_dyck_word(3)
            [1, 1, 1, 0, 0, 0]

        The partition corresponding to ``self.dyck_word()`` is ``self``
        indeed::

            sage: all( p.to_dyck_word().to_partition() == p
            ....:      for p in Partitions(5) )
            True
        """
        from sage.combinat.dyck_word import DyckWord
        if not self._list:
            if n is None:
                return DyckWord([])
            return DyckWord([1]*n + [0]*n)
        list_of_word = []
        if n is None:
            n = max(i + l + 1 for (i, l) in enumerate(self))
            # This n is also max(i+j for (i,j) in self.cells()) + 2.
        list_of_word.extend([1]*(n-self.length()))
        copy_part = list(self)
        while copy_part:
            c = copy_part.pop()
            list_of_word.extend([0]*c)
            for i in range(len(copy_part)):
                copy_part[i] -= c
            list_of_word.append(1)
        list_of_word.extend([0]*(n-self[0]))
        return DyckWord(list_of_word)

    @combinatorial_map(order=2, name="conjugate partition")
    def conjugate(self):
        """
        Return the conjugate partition of the partition ``self``. This
        is also called the associated partition or the transpose in the
        literature.

        EXAMPLES::

            sage: Partition([2,2]).conjugate()
            [2, 2]
            sage: Partition([6,3,1]).conjugate()
            [3, 2, 2, 1, 1, 1]

        The conjugate partition is obtained by transposing the Ferrers
        diagram of the partition (see :meth:`.ferrers_diagram`)::

            sage: print(Partition([6,3,1]).ferrers_diagram())
            ******
            ***
            *
            sage: print(Partition([6,3,1]).conjugate().ferrers_diagram())
            ***
            **
            **
            *
            *
            *
        """
        if not self:
            par = Partitions_n(0)
            return par.element_class(par, [])
        par = Partitions_n(sum(self))
        return par.element_class(par, conjugate(self))

    def glaisher_franklin(self, s):
        r"""
        Apply the Glaisher-Franklin bijection to ``self``.

        The Franklin-Glaisher bijection, with parameter `s`, returns
        a partition whose set of parts that are repeated at least `s`
        times equals the set of parts divisible by `s` in ``self``,
        after dividing each part by `s`.

        INPUT:

        - ``s`` -- positive integer

        EXAMPLES::

            sage: Partition([4, 3, 2, 2, 1]).glaisher_franklin(2)
            [3, 2, 2, 1, 1, 1, 1, 1]

        TESTS:

        The map preserves the size::

            sage: all(mu.glaisher_franklin(s).size() == n
            ....:     for n in range(20) for mu in Partitions(n)
            ....:     for s in range(1, 5))
            True

        The map is bijective::

            sage: l = [[mu.glaisher_franklin(s)
            ....:      for n in range(20) for mu in Partitions(n)]
            ....:     for s in range(1, 5)]
            sage: all(len(set(ls)) == len(ls) for ls in l)
            True

        The map transports the statistics::

            sage: d = lambda la, s: set(p / s for p in la if p % s == 0)
            sage: r = lambda la, s: set(p for p in la if list(la).count(p) >= s)
            sage: all(d(mu, s) == r(mu.glaisher_franklin(s), s)
            ....:     for n in range(20) for mu in Partitions(n)
            ....:     for s in range(1, 5))
            True

        For `s=2`, the map is known to findstat::

            sage: findmap(Partitions, lambda mu: mu.glaisher_franklin(2))       # optional - internet
            0: Mp00312 (quality [100])
        """
        s = ZZ(s)
        if s.is_one():
            return self
        mu = []
        for p, m in enumerate(self.to_exp(), 1):
            if not p % s:
                mu.extend([p // s] * (m*s))
            else:
                for i, v in enumerate(m.digits(s)):
                    mu.extend([p * s**i]*v)

        P = self.parent()
        return P.element_class(P, sorted(mu, reverse=True))

    def glaisher_franklin_inverse(self, s):
        r"""
        Apply the inverse of the Glaisher-Franklin bijection to ``self``.

        The inverse of the Franklin-Glaisher bijection, with
        parameter `s`, returns a partition whose set of parts that
        are divisible by `s`, after dividing each by `s`, equals the
        equals the set of parts repeated at least `s` times in
        ``self``.

        INPUT:

        - ``s`` -- positive integer

        EXAMPLES::

            sage: Partition([4, 3, 2, 2, 1]).glaisher_franklin(2)
            [3, 2, 2, 1, 1, 1, 1, 1]
            sage: Partition([3, 2, 2, 1, 1, 1, 1, 1]).glaisher_franklin_inverse(2)
            [4, 3, 2, 2, 1]

        TESTS:

        The map is inverse to :meth:`glaisher_franklin`::

            sage: all(mu.glaisher_franklin(s).glaisher_franklin_inverse(s) == mu
            ....:     and mu.glaisher_franklin_inverse(s).glaisher_franklin(s) == mu
            ....:     for n in range(20) for mu in Partitions(n)
            ....:     for s in range(1, 5))
            True

        For `s=2`, the map is known to findstat::

            sage: findmap(Partitions, lambda mu: mu.glaisher_franklin_inverse(2))         # optional - internet
            0: Mp00313 (quality [100])
        """
        s = ZZ(s)
        if s.is_one():
            return self
        mu = []
        for p, m in enumerate(self.to_exp(), 1):
            p = ZZ(p)
            mu.extend([p * s] * (m // s))
            m1, p1 = p.val_unit(s)
            mu.extend([p1] * ((m % s) * s**m1))

        P = self.parent()
        return P.element_class(P, sorted(mu, reverse=True))

    def suter_diagonal_slide(self, n, exp=1):
        r"""
        Return the image of ``self`` in `Y_n` under Suter's diagonal slide
        `\sigma_n`, where the notations used are those defined in [Sut2002]_.

        The set `Y_n` is defined as the set of all partitions
        `\lambda` such that the hook length of the `(0, 0)`-cell (i.e. the
        northwestern most cell in English notation) of `\lambda` is less
        than `n`, including the empty partition.

        The map `\sigma_n` sends a partition (with nonzero entries)
        `(\lambda_1, \lambda_2, \ldots, \lambda_m) \in Y_n` to the partition
        `(\lambda_2 + 1, \lambda_3 + 1, \ldots, \lambda_m + 1,
        \underbrace{1, 1, \ldots, 1}_{n - m - \lambda_1\text{ ones}})`.
        In other words, it pads the partition with trailing zeroes
        until it has length `n - \lambda_1`, then removes its first
        part, and finally adds `1` to each part.

        By Theorem 2.1 of [Sut2002]_, the dihedral group `D_n` with
        `2n` elements acts on `Y_n` by letting the primitive rotation
        act as `\sigma_n` and the reflection act as conjugation of
        partitions (:meth:`conjugate()`). This action is faithful if
        `n \geq 3`.

        INPUT:

        - ``n`` -- nonnegative integer

        - ``exp`` -- (default: 1) how many times `\sigma_n` should be applied

        OUTPUT:

        The result of applying Suter's diagonal slide `\sigma_n` to
        ``self``, assuming that ``self`` lies in `Y_n`. If the
        optional argument ``exp`` is set, then the slide
        `\sigma_n` is applied not just once, but ``exp`` times
        (note that ``exp`` is allowed to be negative, since
        the slide has finite order).

        EXAMPLES::

            sage: Partition([5,4,1]).suter_diagonal_slide(8)
            [5, 2]
            sage: Partition([5,4,1]).suter_diagonal_slide(9)
            [5, 2, 1]
            sage: Partition([]).suter_diagonal_slide(7)
            [1, 1, 1, 1, 1, 1]
            sage: Partition([]).suter_diagonal_slide(1)
            []
            sage: Partition([]).suter_diagonal_slide(7, exp=-1)
            [6]
            sage: Partition([]).suter_diagonal_slide(1, exp=-1)
            []
            sage: P7 = Partitions(7)
            sage: all( p == p.suter_diagonal_slide(9, exp=-1).suter_diagonal_slide(9)
            ....:      for p in P7 )
            True
            sage: all( p == p.suter_diagonal_slide(9, exp=3)
            ....:            .suter_diagonal_slide(9, exp=3)
            ....:            .suter_diagonal_slide(9, exp=3)
            ....:      for p in P7 )
            True
            sage: all( p == p.suter_diagonal_slide(9, exp=6)
            ....:            .suter_diagonal_slide(9, exp=6)
            ....:            .suter_diagonal_slide(9, exp=6)
            ....:      for p in P7 )
            True
            sage: all( p == p.suter_diagonal_slide(9, exp=-1)
            ....:            .suter_diagonal_slide(9, exp=1)
            ....:      for p in P7 )
            True

        Check of the assertion in [Sut2002]_ that `\sigma_n\bigl( \sigma_n(
        \lambda^{\prime})^{\prime} \bigr) = \lambda`::

            sage: all( p.suter_diagonal_slide(8).conjugate()
            ....:      == p.conjugate().suter_diagonal_slide(8, exp=-1)
            ....:      for p in P7 )
            True

        Check of Claim 1 in [Sut2002]_::

            sage: P5 = Partitions(5)
            sage: all( all( (p.suter_diagonal_slide(6) in q.suter_diagonal_slide(6).down())
            ....:           or (q.suter_diagonal_slide(6) in p.suter_diagonal_slide(6).down())
            ....:           for p in q.down() )
            ....:      for q in P5 )
            True

        TESTS:

        Check for ``exp = 0``::

            sage: P = Partitions(4)
            sage: all(p == p.suter_diagonal_slide(7, 0) for p in P)
            True

        Check for invalid input::

            sage: p = Partition([2,1])
            sage: p.hook_length(0, 0)
            3
            sage: p.suter_diagonal_slide(2)
            Traceback (most recent call last):
            ...
            ValueError: the hook length must be less than n
        """
        # Check for valid input
        if len(self) > 0 and len(self) + self._list[0] > n:  # >, not >=, since we double count the (0,0) cell
            raise ValueError("the hook length must be less than n")
        ret = self
        # Arbitrary exp
        exp = exp % n  # It is at most order n
        if exp > n / 2:
            exp -= n
        while exp != 0:
            leng = len(ret)
            if exp > 0:
                # Suter's map \sigma_n
                if leng == 0:   # Taking extra care about the empty partition.
                    ret = Partition([1] * (n - 1))
                    exp -= 1
                    continue
                res = [i + 1 for i in ret._list[1:]]
                res += [1] * (n - leng - ret._list[0])
                ret = Partition(res)
                exp -= 1
            else:  # exp < 0 since if exp == 0, we would exit the while loop
                # inverse map \sigma_n^{-1}
                if leng == 0:   # Taking extra care about the empty partition.
                    ret = Partition([n - 1])
                    exp += 1
                    continue
                res = [n - leng - 1]
                res.extend([i - 1 for i in ret._list if i > 1])
                ret = Partition(res)
                exp += 1
        return ret

    @combinatorial_map(name="reading tableau")
    def reading_tableau(self):
        r"""
        Return the RSK recording tableau of the reading word of the
        (standard) tableau `T` labeled down (in English convention)
        each column to the shape of ``self``.

        For an example of the tableau `T`, consider the partition
        `\lambda = (3,2,1)`, then we have::

            1 4 6
            2 5
            3

        For more, see :func:`~sage.combinat.rsk.RSK()`.

        EXAMPLES::

            sage: Partition([3,2,1]).reading_tableau()
            [[1, 3, 6], [2, 5], [4]]
        """
        st = tableau.StandardTableaux(self).first()
        return st.reading_word_permutation().right_tableau()

    @combinatorial_map(name="initial tableau")
    def initial_tableau(self):
        r"""
        Return the :class:`standard tableau<StandardTableau>` which has the
        numbers `1, 2, \ldots, n` where `n` is the :meth:`size` of ``self``
        entered in order from left to right along the rows of each component,
        where the components are ordered from left to right.

        EXAMPLES::

            sage: Partition([3,2,2]).initial_tableau()
            [[1, 2, 3], [4, 5], [6, 7]]
        """
        sigma = list(accumulate([1] + self._list))
        tab = [list(range(sigma[i], sigma[i + 1]))
               for i in range(len(sigma) - 1)]
        return tableau.StandardTableau(tab)

    def initial_column_tableau(self):
        r"""
        Return the initial column tableau of shape ``self``.

        The initial column tableau of shape ``self`` is the standard tableau
        that has the numbers `1` to `n`, where `n` is the :meth:`size` of ``self``,
        entered in order from top to bottom and then left to right down the
        columns of ``self``.

        EXAMPLES::

            sage: Partition([3,2]).initial_column_tableau()
            [[1, 3, 5], [2, 4]]
        """
        return self.conjugate().initial_tableau().conjugate()

    def garnir_tableau(self, *cell):
        r"""
        Return the Garnir tableau of shape ``self`` corresponding to the cell
        ``cell``. If ``cell`` `= (a,c)` then `(a+1,c)` must belong to the
        diagram of ``self``.

        The Garnir tableaux play an important role in integral and
        non-semisimple representation theory because they determine the
        "straightening" rules for the Specht modules over an arbitrary ring.

        The Garnir tableaux are the "first" non-standard tableaux which arise
        when you act by simple transpositions. If `(a,c)` is a cell in the
        Young diagram of a partition, which is not at the bottom of its
        column, then the corresponding Garnir tableau has the integers
        `1, 2, \ldots, n` entered in order from left to right along the rows
        of the diagram up to the cell `(a,c-1)`, then along the cells
        `(a+1,1)` to `(a+1,c)`, then `(a,c)` until the end of row `a` and
        then continuing from left to right in the remaining positions. The
        examples below probably make this clearer!

        .. NOTE::

            The function also sets ``g._garnir_cell``, where ``g`` is the
            resulting Garnir tableau, equal to ``cell`` which is used by
            some other functions.

        EXAMPLES::

            sage: g = Partition([5,3,3,2]).garnir_tableau((0,2)); g.pp()
              1  2  6  7  8
              3  4  5
              9 10 11
             12 13
            sage: g.is_row_strict(); g.is_column_strict()
            True
            False

            sage: Partition([5,3,3,2]).garnir_tableau(0,2).pp()
              1  2  6  7  8
              3  4  5
              9 10 11
             12 13
            sage: Partition([5,3,3,2]).garnir_tableau(2,1).pp()
              1  2  3  4  5
              6  7  8
              9 12 13
             10 11
            sage: Partition([5,3,3,2]).garnir_tableau(2,2).pp()
            Traceback (most recent call last):
            ...
            ValueError: (row+1, col) must be inside the diagram

        .. SEEALSO::

            - :meth:`top_garnir_tableau`
        """
        try:
            (row, col) = cell
        except ValueError:
            (row, col) = cell[0]

        if row + 1 >= len(self) or col >= self[row+1]:
            raise ValueError('(row+1, col) must be inside the diagram')
        g = self.initial_tableau().to_list()
        a = g[row][col]
        g[row][col:] = list(range(a+col+1, g[row+1][col]+1))
        g[row+1][:col+1] = list(range(a, a+col+1))
        g = tableau.Tableau(g)
        g._garnir_cell = (row, col)
        return g

    def top_garnir_tableau(self, e, cell):
        r"""
        Return the most dominant *standard* tableau which dominates the
        corresponding Garnir tableau and has the same ``e``-residue.

        The Garnir tableau play an important role in integral and non-semisimple
        representation theory because they determine the "straightening" rules
        for the Specht modules. The *top Garnir tableaux* arise in the graded
        representation theory of the symmetric groups and higher level Hecke
        algebras. They were introduced in [KMR2012]_.

        If the Garnir node is ``cell=(r,c)`` and `m` and `M` are the entries
        in the cells ``(r,c)`` and ``(r+1,c)``, respectively, in the initial
        tableau then the top ``e``-Garnir tableau is obtained by inserting the
        numbers `m, m+1, \ldots, M` in order from left to right first in the
        cells in row ``r+1`` which are not in the ``e``-Garnir belt, then in
        the cell in rows ``r`` and ``r+1`` which are in the Garnir belt and
        then, finally, in the remaining cells in row ``r`` which are not in
        the Garnir belt. All other entries in the tableau remain unchanged.

        If ``e = 0``, or if there are no ``e``-bricks in either row ``r``
        or ``r+1``, then the top Garnir tableau is the corresponding Garnir
        tableau.

        EXAMPLES::

            sage: Partition([5,4,3,2]).top_garnir_tableau(2,(0,2)).pp()
               1  2  4  5  8
               3  6  7  9
              10 11 12
              13 14
            sage: Partition([5,4,3,2]).top_garnir_tableau(3,(0,2)).pp()
               1  2  3  4  5
               6  7  8  9
              10 11 12
              13 14
            sage: Partition([5,4,3,2]).top_garnir_tableau(4,(0,2)).pp()
               1  2  6  7  8
               3  4  5  9
              10 11 12
              13 14
            sage: Partition([5,4,3,2]).top_garnir_tableau(0,(0,2)).pp()
               1  2  6  7  8
               3  4  5  9
              10 11 12
              13 14

        TESTS::

            sage: Partition([5,4,3,2]).top_garnir_tableau(0,(3,2)).pp()
            Traceback (most recent call last):
            ...
            ValueError: (4,2)=(row+1,col) must be inside the diagram

        REFERENCES:

        - [KMR2012]_
        """
        (row, col) = cell
        if row+1 >= len(self) or col >= self[row+1]:
            raise ValueError(f'({row+1},{col})=(row+1,col) must be inside the diagram')

        g = self.garnir_tableau(cell)   # start with the Garnir tableau and modify

        if e == 0:
            return g             # no more dominant tableau of the same residue

        a = e*int((self[row]-col)/e)    # number of cells in the e-bricks in row `row`
        b = e*int((col+1)/e)            # number of cells in the e-bricks in row `row+1`

        if a == 0 or b == 0:
            return g

        t = g.to_list()
        m = g[row+1][0]                 # smallest  number in 0-Garnir belt
        # now we will put the number m,m+1,...,t[row+1][col] in order into t
        t[row][col:a+col] = [m+col-b+1+i for i in range(a)]
        t[row+1][col-b+1:col+1] = [m+a+col-b+1+i for i in range(b)]
        return tableau.StandardTableau(t)

    def ladder_tableau(self, e, ladder_lengths=False):
        r"""
        Return the ladder tableau of shape ``self``.

        The `e`-*ladder tableau* is the standard Young tableau obtained
        by reading the *ladders*, the set of cells `(i, j)` that differ
        from `(i+e-1, j-1)`, of the partition `\lambda` from left-to-right.

        INPUT:

        - ``e`` -- nonnegative integer; `0` is considered as `\infty`
          (analogous to the characteristic of a ring)
        - ``ladder_sizes`` -- boolean (default: ``False``); if ``True``, also
          return the sizes of the ladders

        .. SEEALSO::

            :meth:`ladders`

        EXAMPLES::

            sage: la = Partition([6, 5, 3, 1])
            sage: ascii_art(la.ladder_tableau(3))
              1  2  3  5  7 10
              4  6  8 11 13
              9 12 14
             15
            sage: la.ladder_tableau(3, ladder_lengths=True)[1]
            [1, 1, 2, 2, 3, 3, 3]

            sage: ascii_art(la.ladder_tableau(0))
              1  2  3  4  5  6
              7  8  9 10 11
             12 13 14
             15
            sage: all(ll == 1 for ll in la.ladder_tableau(0, ladder_lengths=True)[1])
            True
        """
        Tlad = [[None] * val for val in self]
        counter = 0
        start = 0
        n = sum(self)
        sizes = []
        e = e - 1 if e > 0 else n  # change to the slope
        while counter < n:
            cur = start
            size = 0
            for i, val in enumerate(self):
                if cur < 0:
                    break
                if cur < val:
                    counter += 1
                    Tlad[i][cur] = counter
                    size += 1
                cur -= e
            if ladder_lengths and size:
                sizes.append(size)
            start += 1
        ret = tableau.StandardTableaux(self)(Tlad)
        if ladder_lengths:
            return (ret, sizes)
        return ret

    def ladders(self, e):
        r"""
        Return a dictionary containing the ladders in the diagram of ``self``.

        For `e > 0`, a node `(i, j)` in a partition belongs to the `l`-th
        `e`-ladder if `l = (e - 1) r + c`.

        INPUT:

        - ``e`` -- nonnegative integer; if ``0``, then we
          set ``e = self.size() + 1``

        EXAMPLES::

            sage: Partition([3, 2]).ladders(3)
            {0: [(0, 0)], 1: [(0, 1)], 2: [(0, 2), (1, 0)], 3: [(1, 1)]}

        When ``e`` is ``0``, the cells are in bijection with the ladders,
        but the index of the ladder depends on the size of the partition::

            sage: Partition([3, 2]).ladders(0)
            {0: [(0, 0)], 1: [(0, 1)], 2: [(0, 2)], 5: [(1, 0)], 6: [(1, 1)]}
            sage: Partition([3, 2, 1]).ladders(0)
            {0: [(0, 0)], 1: [(0, 1)], 2: [(0, 2)], 6: [(1, 0)], 7: [(1, 1)],
             12: [(2, 0)]}
            sage: Partition([3, 1, 1]).ladders(0)
            {0: [(0, 0)], 1: [(0, 1)], 2: [(0, 2)], 5: [(1, 0)], 10: [(2, 0)]}
            sage: Partition([1, 1, 1]).ladders(0)
            {0: [(0, 0)], 3: [(1, 0)], 6: [(2, 0)]}
        """
        if e == 0:
            e = sum(self) + 1
        ladders = {}
        for row, val in enumerate(self):
            for col in range(val):
                ell = col + row * (e - 1)
                if ell not in ladders:
                    ladders[ell] = []
                ladders[ell].append((row, col))
        return ladders

    @cached_method
    def young_subgroup(self):
        r"""
        Return the corresponding Young, or parabolic, subgroup of the symmetric
        group.

        The Young subgroup of a partition
        `\lambda = (\lambda_1, \lambda_2, \ldots, \lambda_{\ell})` of `n` is
        the group:

        .. MATH::

            S_{\lambda_1} \times S_{\lambda_2} \times \cdots \times
            S_{\lambda_{\ell}}

        embedded into `S_n` in the standard way (i.e.,
        the `S_{\lambda_i}` factor acts on the numbers from
        `\lambda_1 + \lambda_2 + \cdots + \lambda_{i-1} + 1` to
        `\lambda_1 + \lambda_2 + \cdots + \lambda_i`).

        EXAMPLES::

            sage: Partition([4,2]).young_subgroup()                                     # needs sage.groups
            Permutation Group with generators [(), (5,6), (3,4), (2,3), (1,2)]
        """
        gens = []
        m = 0
        for row in self:
            gens.extend([(c, c+1) for c in range(m+1, m+row)])
            m += row
        gens.append(list(range(1, self.size() + 1)))  # to ensure we get a subgroup of Sym_n
        return PermutationGroup(gens)

    def young_subgroup_generators(self):
        r"""
        Return an indexing set for the generators of the corresponding Young
        subgroup. Here the generators correspond to the simple adjacent
        transpositions `s_i = (i \; i+1)`.

        EXAMPLES::

            sage: Partition([4,2]).young_subgroup_generators()
            [1, 2, 3, 5]
            sage: Partition([1,1,1]).young_subgroup_generators()
            []
            sage: Partition([2,2]).young_subgroup_generators()
            [1, 3]

        .. SEEALSO::

            :meth:`young_subgroup`
        """
        gens = []
        m = 0
        for row in self:
            gens.extend(range(m + 1, m + row))
            m += row
        return gens

    @cached_method
    def _initial_degree(self, e, multicharge=(0,)):
        r"""
        Return the Brundan-Kleshchev-Wang degree of the initial row tableau
        of shape ``self``.

        This degree depends only the shape of the tableau and it is
        used as the base case for computing the degrees of all tableau
        of shape ``self``, which is why this method is cached. See
        :meth:`sage.combinat.tableau.Tableau.degree` for more information.

        EXAMPLES::

            sage: Partition([5,3,2])._initial_degree(0)
            0
            sage: Partition([5,3,2])._initial_degree(2)
            4
            sage: Partition([5,3,2])._initial_degree(3)
            2
            sage: Partition([5,3,2])._initial_degree(4)
            1
        """
        if e == 0:
            return ZZ.zero()
        else:
            return sum(m // e for m in self)

    def degree(self, e):
        r"""
        Return the ``e``-th degree of ``self``.

        The `e`-th degree of a partition `\lambda` is the sum of the `e`-th
        degrees of the standard tableaux of shape `\lambda`. The `e`-th degree
        is the exponent of `\Phi_e(q)` in the Gram determinant of the Specht
        module for a semisimple Iwahori-Hecke algebra of type `A` with
        parameter `q`.

        INPUT:

        - ``e`` -- an  integer  `e > 1`

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: Partition([4,3]).degree(2)
            28
            sage: Partition([4,3]).degree(3)
            15
            sage: Partition([4,3]).degree(4)
            8
            sage: Partition([4,3]).degree(5)
            13
            sage: Partition([4,3]).degree(6)
            0
            sage: Partition([4,3]).degree(7)
            0

        Therefore, the Gram determinant of `S(5,3)` when the Hecke parameter
        `q` is "generic" is

        .. MATH::

            q^N \Phi_2(q)^{28} \Phi_3(q)^{15} \Phi_4(q)^8 \Phi_5(q)^{13}

        for some integer `N`. Compare with :meth:`prime_degree`.
        """
        return sum(t.degree(e) for t in self.standard_tableaux())

    def prime_degree(self, p):
        r"""
        Return the prime degree for the prime integer``p`` for ``self``.

        INPUT:

        - ``p`` -- prime integer

        OUTPUT: nonnegative integer

        The degree of a partition `\lambda` is the sum of the
        `e`-:meth:`degree` of the standard tableaux of shape `\lambda`, for
        `e` a power of the prime `p`. The prime degree gives the exponent of
        `p` in the Gram determinant of the integral Specht module of the
        symmetric group.

        EXAMPLES::

            sage: Partition([4,3]).prime_degree(2)
            36
            sage: Partition([4,3]).prime_degree(3)
            15
            sage: Partition([4,3]).prime_degree(5)
            13
            sage: Partition([4,3]).prime_degree(7)
            0

        Therefore, the Gram determinant of `S(5,3)` when `q = 1` is
        `2^{36} 3^{15} 5^{13}`.  Compare with :meth:`degree`.
        """
        ps = [p]

        while ps[-1] * p < self.size():
            ps.append(ps[-1] * p)
        return sum(t.degree(pk) for pk in ps for t in self.standard_tableaux())

    def arm_length(self, i, j):
        r"""
        Return the length of the arm of cell `(i,j)` in ``self``.

        The arm of cell `(i,j)` is the cells that appear to the right of
        cell `(i,j)`.

        The cell coordinates are zero-based, i. e., the northwesternmost
        cell is `(0,0)`.

        INPUT:

        - ``i``, ``j`` -- two integers

        OUTPUT: integer or a :exc:`ValueError`

        EXAMPLES::

            sage: p = Partition([2,2,1])
            sage: p.arm_length(0, 0)
            1
            sage: p.arm_length(0, 1)
            0
            sage: p.arm_length(2, 0)
            0
            sage: Partition([3,3]).arm_length(0, 0)
            2
            sage: Partition([3,3]).arm_length(*[0,0])
            2
        """
        p = self
        if i < len(p) and j < p[i]:
            return p[i]-(j+1)
        raise ValueError("the cell is not in the diagram")

    def arm_lengths(self, flat=False):
        """
        Return a tableau of shape ``self`` where each cell is filled with
        its arm length.

        The optional boolean parameter ``flat`` provides the option of
        returning a flat list.

        EXAMPLES::

            sage: Partition([2,2,1]).arm_lengths()
            [[1, 0], [1, 0], [0]]
            sage: Partition([2,2,1]).arm_lengths(flat=True)
            [1, 0, 1, 0, 0]
            sage: Partition([3,3]).arm_lengths()
            [[2, 1, 0], [2, 1, 0]]
            sage: Partition([3,3]).arm_lengths(flat=True)
            [2, 1, 0, 2, 1, 0]
        """
        p = self
        if not flat:
            return [[pi - (j + 1) for j in range(pi)] for pi in p]
        return [pi - (j + 1) for pi in p for j in range(pi)]

    def arm_cells(self, i, j):
        r"""
        Return the list of the cells of the arm of cell `(i,j)` in ``self``.

        The arm of cell `c = (i,j)` is the boxes that appear to the right of
        `c`.

        The cell coordinates are zero-based, i. e., the northwesternmost
        cell is `(0,0)`.

        INPUT:

        - ``i``, ``j`` -- two integers

        OUTPUT: list of pairs of integers

        EXAMPLES::

            sage: Partition([4,4,3,1]).arm_cells(1,1)
            [(1, 2), (1, 3)]

            sage: Partition([]).arm_cells(0,0)
            Traceback (most recent call last):
            ...
            ValueError: the cell is not in the diagram
        """
        p = self
        if i < len(p) and j < p[i]:
            return [(i, x) for x in range(j + 1, p[i])]
        raise ValueError("the cell is not in the diagram")

    def leg_length(self, i, j):
        """
        Return the length of the leg of cell `(i,j)` in ``self``.

        The leg of cell `c = (i,j)` is defined to be the cells below `c`
        (in English convention).

        The cell coordinates are zero-based, i. e., the northwesternmost
        cell is `(0,0)`.

        INPUT:

        - ``i``, ``j`` -- two integers

        OUTPUT: integer or a :exc:`ValueError`

        EXAMPLES::

            sage: p = Partition([2,2,1])
            sage: p.leg_length(0, 0)
            2
            sage: p.leg_length(0,1)
            1
            sage: p.leg_length(2,0)
            0
            sage: Partition([3,3]).leg_length(0, 0)
            1
            sage: cell = [0,0]; Partition([3,3]).leg_length(*cell)
            1
        """
        conj = self.conjugate()
        if j < len(conj) and i < conj[j]:
            return conj[j] - (i + 1)
        raise ValueError("the cell is not in the diagram")

    def leg_lengths(self, flat=False):
        """
        Return a tableau of shape ``self`` with each cell filled in with
        its leg length.  The optional boolean parameter ``flat`` provides
        the option of returning a flat list.

        EXAMPLES::

            sage: Partition([2,2,1]).leg_lengths()
            [[2, 1], [1, 0], [0]]
            sage: Partition([2,2,1]).leg_lengths(flat=True)
            [2, 1, 1, 0, 0]
            sage: Partition([3,3]).leg_lengths()
            [[1, 1, 1], [0, 0, 0]]
            sage: Partition([3,3]).leg_lengths(flat=True)
            [1, 1, 1, 0, 0, 0]
        """
        p = self
        conj = p.conjugate()
        if not flat:
            return [[conj[j] - (i + 1) for j in range(pi)]
                    for i, pi in enumerate(p)]
        return [conj[j] - (i + 1) for i, pi in enumerate(p)
                for j in range(pi)]

    def leg_cells(self, i, j):
        r"""
        Return the list of the cells of the leg of cell `(i,j)` in ``self``.

        The leg of cell `c = (i,j)` is defined to be the cells below `c` (in
        English convention).

        The cell coordinates are zero-based, i. e., the northwesternmost
        cell is `(0,0)`.

        INPUT:

        - ``i``, ``j`` -- two integers

        OUTPUT: list of pairs of integers

        EXAMPLES::

            sage: Partition([4,4,3,1]).leg_cells(1,1)
            [(2, 1)]
            sage: Partition([4,4,3,1]).leg_cells(0,1)
            [(1, 1), (2, 1)]

            sage: Partition([]).leg_cells(0,0)
            Traceback (most recent call last):
            ...
            ValueError: the cell is not in the diagram
        """
        l = self.leg_length(i, j)
        return [(x, j) for x in range(i + 1, i + l + 1)]

    def attacking_pairs(self):
        """
        Return a list of the attacking pairs of the Young diagram of
        ``self``.

        A pair of cells `(c, d)` of a Young diagram (in English notation) is
        said to be attacking if one of the following conditions holds:

        1. `c` and `d` lie in the same row with `c` strictly to the west
           of `d`.

        2. `c` is in the row immediately to the south of `d`, and `c`
           lies strictly east of `d`.

        This particular method returns each pair `(c, d)` as a tuple,
        where each of `c` and `d` is given as a tuple `(i, j)` with
        `i` and `j` zero-based (so `i = 0` means that the cell lies
        in the topmost row).

        EXAMPLES::

            sage: p = Partition([3, 2])
            sage: p.attacking_pairs()
            [((0, 0), (0, 1)),
             ((0, 0), (0, 2)),
             ((0, 1), (0, 2)),
             ((1, 0), (1, 1)),
             ((1, 1), (0, 0))]
            sage: Partition([]).attacking_pairs()
            []
        """
        attacking_pairs = []
        for i, r in enumerate(self):
            for j in range(r):
                # c is in position (i,j)
                # Find the d that satisfy condition 1
                attacking_pairs.extend(((i, j), (i, k))
                                       for k in range(j + 1, r))

                # Find the d that satisfy condition 2
                if i == 0:
                    continue
                attacking_pairs.extend(((i, j), (i - 1, k))
                                       for k in range(j))

        return attacking_pairs

    def dominated_partitions(self, rows=None):
        """
        Return a list of the partitions dominated by `n`. If ``rows`` is
        specified, then it only returns the ones whose number of rows
        is at most ``rows``.

        EXAMPLES::

            sage: Partition([3,2,1]).dominated_partitions()
            [[3, 2, 1], [3, 1, 1, 1], [2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]
            sage: Partition([3,2,1]).dominated_partitions(rows=3)
            [[3, 2, 1], [2, 2, 2]]
        """
        # Naive implementation because iteration is so fast
        n = sum(self)
        P = Partitions_n(n)
        if rows:
            return [P(x) for x in ZS1_iterator_nk(n, rows) if self.dominates(x)]
        else:
            return [P(x) for x in ZS1_iterator(n) if self.dominates(x)]

    def contains(self, x):
        """
        Return ``True`` if ``x`` is a partition whose Ferrers diagram is
        contained in the Ferrers diagram of ``self``.

        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p.contains([2,1])
            True
            sage: all(p.contains(mu) for mu in Partitions(3))
            True
            sage: all(p.contains(mu) for mu in Partitions(4))
            False
        """
        return len(self) >= len(x) and all(self[i] >= x[i] for i in range(len(x)))

    def hook_product(self, a):
        """
        Return the Jack hook-product.

        EXAMPLES::

            sage: Partition([3,2,1]).hook_product(x)                                    # needs sage.symbolic
            (2*x + 3)*(x + 2)^2
            sage: Partition([2,2]).hook_product(x)                                      # needs sage.symbolic
            2*(x + 2)*(x + 1)
        """

        nu = self.conjugate()
        res = 1
        for i in range(len(self)):
            for j in range(self[i]):
                res *= a*(self[i]-j-1)+nu[j]-i
        return res

    def hook_polynomial(self, q, t):
        """
        Return the two-variable hook polynomial.

        EXAMPLES::

            sage: R.<q,t> = PolynomialRing(QQ)
            sage: a = Partition([2,2]).hook_polynomial(q,t)
            sage: a == (1 - t)*(1 - q*t)*(1 - t^2)*(1 - q*t^2)
            True
            sage: a = Partition([3,2,1]).hook_polynomial(q,t)
            sage: a == (1 - t)^3*(1 - q*t^2)^2*(1 - q^2*t^3)
            True
        """
        nu = self.conjugate()
        res = 1
        for i in range(len(self)):
            for j in range(self[i]):
                res *= 1-q**(self[i]-j-1)*t**(nu[j]-i)
        return res

    def hook_length(self, i, j):
        r"""
        Return the length of the hook of cell `(i,j)` in ``self``.

        The (length of the) hook of cell `(i,j)` of a partition `\lambda`
        is

        .. MATH::

            \lambda_i + \lambda^{\prime}_j - i - j + 1

        where `\lambda^{\prime}` is the conjugate partition. In English
        convention, the hook length is the number of cells horizontally
        to the right and vertically below the cell `(i,j)` (including
        that cell).

        EXAMPLES::

            sage: p = Partition([2,2,1])
            sage: p.hook_length(0, 0)
            4
            sage: p.hook_length(0, 1)
            2
            sage: p.hook_length(2, 0)
            1
            sage: Partition([3,3]).hook_length(0, 0)
            4
            sage: cell = [0,0]; Partition([3,3]).hook_length(*cell)
            4
        """
        return self.leg_length(i, j) + self.arm_length(i, j) + 1

    def hooks(self):
        """
        Return a sorted list of the hook lengths in ``self``.

        EXAMPLES::

            sage: Partition([3,2,1]).hooks()
            [5, 3, 3, 1, 1, 1]
        """
        res = []
        for row in self.hook_lengths():
            res += row
        res.sort(reverse=True)
        return res

    def hook_lengths(self):
        r"""
        Return a tableau of shape ``self`` with the cells filled in with the
        hook lengths.

        In each cell, put the sum of one plus the number of cells
        horizontally to the right and vertically below the cell (the
        hook length).

        For example, consider the partition ``[3,2,1]`` of 6 with Ferrers
        diagram::

            # # #
            # #
            #

        When we fill in the cells with the hook lengths, we obtain::

            5 3 1
            3 1
            1

        EXAMPLES::

            sage: Partition([2,2,1]).hook_lengths()
            [[4, 2], [3, 1], [1]]
            sage: Partition([3,3]).hook_lengths()
            [[4, 3, 2], [3, 2, 1]]
            sage: Partition([3,2,1]).hook_lengths()
            [[5, 3, 1], [3, 1], [1]]
            sage: Partition([2,2]).hook_lengths()
            [[3, 2], [2, 1]]
            sage: Partition([5]).hook_lengths()
            [[5, 4, 3, 2, 1]]

        REFERENCES:

        - http://mathworld.wolfram.com/HookLengthFormula.html
        """
        p = self
        conj = p.conjugate()
        return [[p[i]-(i+1)+conj[j]-(j+1)+1 for j in range(p[i])] for i in range(len(p))]

    def upper_hook(self, i, j, alpha):
        r"""
        Return the upper hook length of the cell `(i,j)` in ``self``.
        When ``alpha = 1``, this is just the normal hook length.

        The upper hook length of a cell `(i,j)` in a partition
        `\kappa` is defined by

        .. MATH::

            h^*_\kappa(i,j) = \kappa^\prime_j - i + \alpha(\kappa_i - j + 1).

        EXAMPLES::

            sage: p = Partition([2,1])
            sage: p.upper_hook(0,0,1)
            3
            sage: p.hook_length(0,0)
            3
            sage: [ p.upper_hook(i,j,x) for i,j in p.cells() ]                          # needs sage.symbolic
            [2*x + 1, x, x]
        """
        p = self
        conj = self.conjugate()
        return conj[j] - (i+1) + alpha*(p[i]-j)

    def upper_hook_lengths(self, alpha):
        r"""
        Return a tableau of shape ``self`` with the cells filled in with the
        upper hook lengths. When ``alpha = 1``, these are just the normal hook
        lengths.

        The upper hook length of a cell `(i,j)` in a partition
        `\kappa` is defined by

        .. MATH::

            h^*_\kappa(i,j) = \kappa^\prime_j - i + \alpha(\kappa_i - j + 1).

        EXAMPLES::

            sage: Partition([3,2,1]).upper_hook_lengths(x)                              # needs sage.symbolic
            [[3*x + 2, 2*x + 1, x], [2*x + 1, x], [x]]
            sage: Partition([3,2,1]).upper_hook_lengths(1)
            [[5, 3, 1], [3, 1], [1]]
            sage: Partition([3,2,1]).hook_lengths()
            [[5, 3, 1], [3, 1], [1]]
        """
        p = self
        conj = p.conjugate()
        return [[conj[j] - (i+1) + alpha*(p[i]-j) for j in range(p[i])] for i in range(len(p))]

    def lower_hook(self, i, j, alpha):
        r"""
        Return the lower hook length of the cell `(i,j)` in ``self``.
        When ``alpha = 1``, this is just the normal hook length.

        The lower hook length of a cell `(i,j)` in a partition
        `\kappa` is defined by

        .. MATH::

            h_*^\kappa(i,j) = \kappa^\prime_j - i + 1 + \alpha(\kappa_i - j).

        EXAMPLES::

            sage: p = Partition([2,1])
            sage: p.lower_hook(0,0,1)
            3
            sage: p.hook_length(0,0)
            3
            sage: [ p.lower_hook(i,j,x) for i,j in p.cells() ]                          # needs sage.symbolic
            [x + 2, 1, 1]
        """
        p = self
        conj = self.conjugate()
        return conj[j] - i + alpha*(p[i] - (j+1))

    def lower_hook_lengths(self, alpha):
        r"""
        Return a tableau of shape ``self`` with the cells filled in with the
        lower hook lengths. When ``alpha = 1``, these are just the normal hook
        lengths.

        The lower hook length of a cell `(i,j)` in a partition
        `\kappa` is defined by

        .. MATH::

            h_*^\kappa(i,j) = \kappa^\prime_j - i + 1 + \alpha(\kappa_i - j).

        EXAMPLES::

            sage: Partition([3,2,1]).lower_hook_lengths(x)                              # needs sage.symbolic
            [[2*x + 3, x + 2, 1], [x + 2, 1], [1]]
            sage: Partition([3,2,1]).lower_hook_lengths(1)
            [[5, 3, 1], [3, 1], [1]]
            sage: Partition([3,2,1]).hook_lengths()
            [[5, 3, 1], [3, 1], [1]]
        """
        p = self
        conj = p.conjugate()
        return [[conj[j] - i + alpha*(p[i] - (j + 1)) for j in range(p[i])]
                for i in range(len(p))]

    def weighted_size(self):
        r"""
        Return the weighted size of ``self``.

        The weighted size of a partition `\lambda` is

        .. MATH::

            \sum_i i \cdot \lambda_i,

        where `\lambda = (\lambda_0, \lambda_1, \lambda_2, \cdots )`.

        This also the sum of the leg length of every cell in `\lambda`, or

        .. MATH::

            \sum_i \binom{\lambda^{\prime}_i}{2}

        where `\lambda^{\prime}` is the conjugate partition of `\lambda`.

        EXAMPLES::

            sage: Partition([2,2]).weighted_size()
            2
            sage: Partition([3,3,3]).weighted_size()
            9
            sage: Partition([5,2]).weighted_size()
            2
            sage: Partition([]).weighted_size()
            0
        """
        p = self
        return sum([i*p[i] for i in range(len(p))])

    def is_empty(self):
        """
        Return ``True`` if ``self`` is the empty partition.

        EXAMPLES::

            sage: Partition([]).is_empty()
            True
            sage: Partition([2,1,1]).is_empty()
            False
        """
        return len(self) == 0

    def length(self):
        """
        Return the number of parts in ``self``.

        EXAMPLES::

            sage: Partition([3,2]).length()
            2
            sage: Partition([2,2,1]).length()
            3
            sage: Partition([]).length()
            0
        """
        return len(self)

    def to_exp(self, k=0):
        """
        Return a list of the multiplicities of the parts of a partition.
        Use the optional parameter ``k`` to get a return list of length at
        least ``k``.

        EXAMPLES::

            sage: Partition([3,2,2,1]).to_exp()
            [1, 2, 1]
            sage: Partition([3,2,2,1]).to_exp(5)
            [1, 2, 1, 0, 0]

        TESTS::

            sage: [parent(x) for x in Partition([3,2,2,1]).to_exp(5)]
            [Integer Ring, Integer Ring, Integer Ring, Integer Ring, Integer Ring]
        """
        p = self
        if len(p) > 0:
            k = max(k, p[0])
        a = [ZZ.zero()] * k
        for i in p:
            a[i-1] += 1
        return a

    def evaluation(self):
        r"""
        Return the evaluation of ``self``.

        The **commutative evaluation**, often shortened to **evaluation**, of
        a word (we think of a partition as a word in `\{1, 2, 3, \ldots\}`)
        is its image in the free commutative monoid. In other words,
        this counts how many occurrences there are of each letter.

        This is also is known as **Parikh vector** and **abelianization** and
        has the same output as :meth:`to_exp()`.

        EXAMPLES::

            sage: Partition([4,3,1,1]).evaluation()
            [2, 0, 1, 1]
        """
        return self.to_exp()

    def to_exp_dict(self):
        """
        Return a dictionary containing the multiplicities of the parts of
        ``self``.

        EXAMPLES::

            sage: p = Partition([4,2,2,1])
            sage: d = p.to_exp_dict()
            sage: d[4]
            1
            sage: d[2]
            2
            sage: d[1]
            1
            sage: 5 in d
            False
        """
        d = {}
        for part in self:
            d[part] = d.get(part, 0) + 1
        return d

    def centralizer_size(self, t=0, q=0):
        r"""
        Return the size of the centralizer of any permutation of cycle type
        ``self``.

        If `m_i` is the multiplicity of `i` as a part of `p`, this is given by

        .. MATH::

           \prod_i m_i! i^{m_i}.

        Including the optional parameters `t` and `q` gives the `q,t` analog,
        which is the former product times

        .. MATH::

           \prod_{i=1}^{\mathrm{length}(p)} \frac{1 - q^{p_i}}{1 - t^{p_i}}.

        See Section 1.3, p. 24, in [Ke1991]_.

        EXAMPLES::

            sage: Partition([2,2,1]).centralizer_size()
            8
            sage: Partition([2,2,2]).centralizer_size()
            48
            sage: Partition([2,2,1]).centralizer_size(q=2, t=3)
            9/16
            sage: Partition([]).centralizer_size()
            1
            sage: Partition([]).centralizer_size(q=2, t=4)
            1

        TESTS::

            sage: Partition([2,2,2]).aut()
            48
        """
        size = prod(i**mi * factorial(mi)
                    for i, mi in self.to_exp_dict().items())
        if t or q:
            size *= prod((ZZ.one() - q ** j) / (ZZ.one() - t ** j)
                         for j in self)
        return size

    aut = centralizer_size

    def content(self, r, c, multicharge=(0,)):
        r"""
        Return the content of the cell at row `r` and column `c`.

        The content of a cell is `c - r`.

        For consistency with partition tuples there is also an optional
        ``multicharge`` argument which is an offset to the usual content. By
        setting the ``multicharge`` equal to the 0-element of the ring
        `\ZZ/e\ZZ`, the corresponding `e`-residue will be returned. This is
        the content modulo `e`.

        The content (and residue) do not strictly depend on the partition,
        however, this method is included because it is often useful in the
        context of partitions.

        EXAMPLES::

            sage: Partition([2,1]).content(1,0)
            -1
            sage: p = Partition([3,2])
            sage: sum([p.content(*c) for c in p.cells()])
            2

        and now we return the 3-residue of a cell::

            sage: Partition([2,1]).content(1,0, multicharge=[IntegerModRing(3)(0)])
            2
        """
        return c - r + multicharge[0]

    def residue(self, r, c, l):
        r"""
        Return the ``l``-residue of the cell at row ``r`` and column ``c``.

        The `\ell`-residue of a cell is `c - r` modulo `\ell`.

        This does not strictly depend upon the partition, however, this method
        is included because it is often useful in the context of partitions.

        EXAMPLES::

            sage: Partition([2,1]).residue(1, 0, 3)
            2
        """
        return (c - r) % l

    @cached_method
    def block(self, e, multicharge=(0,)):
        r"""
        Return a dictionary `\beta` that determines the block associated to
        the partition ``self`` and the
        :meth:`~sage.combinat.tableau_residues.ResidueSequence.quantum_characteristic` ``e``.

        INPUT:

        - ``e`` -- the quantum characteristic

        - ``multicharge`` -- the multicharge (default: `(0,)`)

        OUTPUT:

        - A dictionary giving the multiplicities of the residues in the
          partition tuple ``self``

        In more detail, the value ``beta[i]`` is equal to the
        number of nodes of residue ``i``. This corresponds to
        the positive root

        .. MATH::

            \sum_{i\in I} \beta_i \alpha_i \in Q^+,

        a element of the positive root lattice of the corresponding
        Kac-Moody algebra. See [DJM1998]_ and [BK2009]_ for more details.

        This is a useful statistics because two Specht modules for a
        Hecke algebra of type `A` belong to the same block if and only if they
        correspond to same element `\beta` of the root lattice, given above.

        We return a dictionary because when the quantum characteristic is `0`,
        the Cartan type is `A_{\infty}`, in which case the simple roots are
        indexed by the integers.

        EXAMPLES::

            sage: Partition([4,3,2]).block(0)
            {-2: 1, -1: 2, 0: 2, 1: 2, 2: 1, 3: 1}
            sage: Partition([4,3,2]).block(2)
            {0: 4, 1: 5}
            sage: Partition([4,3,2]).block(2, multicharge=(1,))
            {0: 5, 1: 4}
            sage: Partition([4,3,2]).block(3)
            {0: 3, 1: 3, 2: 3}
            sage: Partition([4,3,2]).block(4)
            {0: 2, 1: 2, 2: 2, 3: 3}
        """
        block = {}
        Ie = IntegerModRing(e)
        for (r, c) in self.cells():
            i = Ie(multicharge[0] + c - r)
            block[i] = block.get(i, 0) + 1
        return block

    def defect(self, e, multicharge=(0,)):
        r"""
        Return the ``e``-defect or the ``e``-weight of ``self``.

        The `e`-defect is the number of (connected) `e`-rim hooks that
        can be removed from the partition.

        The defect of a partition is given by

        .. MATH::

            \text{defect}(\beta) = (\Lambda, \beta) - \tfrac12(\beta, \beta),

        where `\Lambda = \sum_r \Lambda_{\kappa_r}` for the multicharge
        `(\kappa_1, \ldots, \kappa_{\ell})` and
        `\beta = \sum_{(r,c)} \alpha_{(c-r) \pmod e}`, with the sum
        being over the cells in the partition.

        INPUT:

        - ``e`` -- the quantum characteristic

        - ``multicharge`` -- the multicharge (default: `(0,)`)

        OUTPUT: nonnegative integer, which is the defect of the block
        containing the partition ``self``

        EXAMPLES::

            sage: Partition([4,3,2]).defect(2)
            3
            sage: Partition([0]).defect(2)
            0
            sage: Partition([3]).defect(2)
            1
            sage: Partition([6]).defect(2)
            3
            sage: Partition([9]).defect(2)
            4
            sage: Partition([12]).defect(2)
            6
            sage: Partition([4,3,2]).defect(3)
            3
            sage: Partition([0]).defect(3)
            0
            sage: Partition([3]).defect(3)
            1
            sage: Partition([6]).defect(3)
            2
            sage: Partition([9]).defect(3)
            3
            sage: Partition([12]).defect(3)
            4

        TESTS::

            sage: all(mu.core(e).size() + e * mu.defect(e) == 9
            ....:     for mu in Partitions(9) for e in [2,3,4])
            True
        """
        beta = self.block(e, multicharge)
        Ie = IntegerModRing(e)
        return beta.get(multicharge[0], 0) - sum(beta[r]**2 - beta[r] * beta.get(Ie(r+1), 0)
                                                 for r in beta)

    def contents_tableau(self, multicharge=(0,)):
        """
        Return the tableau which has ``(k,r,c)``-th cell equal to the
        content ``multicharge[k] - r + c`` of the cell.

        EXAMPLES::

            sage: Partition([2,1]).contents_tableau()
            [[0, 1], [-1]]
            sage: Partition([3,2,1,1]).contents_tableau().pp()
                0  1  2
                -1  0
                -2
                -3
            sage: Partition([3,2,1,1]).contents_tableau([ IntegerModRing(3)(0)] ).pp()
                0  1  2
                2  0
                1
                0
        """
        return tableau.Tableau([[multicharge[0]-r+c for c in range(self[r])]
                                for r in range(len(self))])

    def is_restricted(self, e, multicharge=(0,)):
        """
        Return ``True`` is this is an ``e``-restricted partition.

        An `e`-restricted partition is a partition such that the
        difference of consecutive parts is always strictly less
        than `e`, where partitions are considered to have an infinite
        number of `0` parts. I.e., the last part must be strictly
        less than `e`.

        EXAMPLES::

          sage: Partition([4,3,3,2]).is_restricted(2)
          False
          sage: Partition([4,3,3,2]).is_restricted(3)
          True
          sage: Partition([4,3,3,2]).is_restricted(4)
          True
          sage: Partition([4]).is_restricted(4)
          False
        """
        return (not self
                or (self[-1] < e and all(self[r] - self[r+1] < e for r in range(len(self) - 1))))

    def is_regular(self, e, multicharge=(0,)) -> bool:
        """
        Return ``True`` is this is an ``e``-regular partition.

        A partition is `e`-regular if it does not have `e` equal
        nonzero parts.

        EXAMPLES::

          sage: Partition([4,3,3,3]).is_regular(2)
          False
          sage: Partition([4,3,3,3]).is_regular(3)
          False
          sage: Partition([4,3,3,3]).is_regular(4)
          True
        """
        return all(self[r] > self[r+e-1] for r in range(len(self)-e+1))

    def conjugacy_class_size(self):
        """
        Return the size of the conjugacy class of the symmetric group
        indexed by ``self``.

        EXAMPLES::

            sage: Partition([2,2,2]).conjugacy_class_size()
            15
            sage: Partition([2,2,1]).conjugacy_class_size()
            15
            sage: Partition([2,1,1]).conjugacy_class_size()
            6
        """
        return factorial(sum(self)) / self.centralizer_size()

    def corners(self) -> list:
        r"""
        Return a list of the corners of the partition ``self``.

        A corner of a partition `\lambda` is a cell of the Young diagram
        of `\lambda` which can be removed from the Young diagram while
        still leaving a straight shape behind.

        The entries of the list returned are pairs of the form `(i,j)`,
        where `i` and `j` are the coordinates of the respective corner.
        The coordinates are counted from `0`.

        .. NOTE::

            This is referred to as an "inner corner" in [Sag2001]_.

        EXAMPLES::

            sage: Partition([3,2,1]).corners()
            [(0, 2), (1, 1), (2, 0)]
            sage: Partition([3,3,1]).corners()
            [(1, 2), (2, 0)]
            sage: Partition([]).corners()
            []
        """
        p = self
        if p.is_empty():
            return []

        lcors = [[0, p[0]-1]]
        nn = len(p)
        if nn == 1:
            return [tuple(c) for c in lcors]

        lcors_index = 0
        for i in range(1, nn):
            if p[i] == p[i-1]:
                lcors[lcors_index][0] += 1
            else:
                lcors.append([i, p[i]-1])
                lcors_index += 1

        return [tuple(c) for c in lcors]

    inside_corners = corners
    removable_cells = corners     # for compatibility with partition tuples

    def corners_residue(self, i, l):
        r"""
        Return a list of the corners of the partition ``self`` having
        ``l``-residue ``i``.

        A corner of a partition `\lambda` is a cell of the Young diagram
        of `\lambda` which can be removed from the Young diagram while
        still leaving a straight shape behind. See :meth:`residue` for
        the definition of the ``l``-residue.

        The entries of the list returned are pairs of the form `(i,j)`,
        where `i` and `j` are the coordinates of the respective corner.
        The coordinates are counted from `0`.

        EXAMPLES::

            sage: Partition([3,2,1]).corners_residue(0, 3)
            [(1, 1)]
            sage: Partition([3,2,1]).corners_residue(1, 3)
            [(2, 0)]
            sage: Partition([3,2,1]).corners_residue(2, 3)
            [(0, 2)]
        """
        return [x for x in self.corners() if self.residue(*x, l=l) == i]

    inside_corners_residue = corners_residue
    removable_cells_residue = corners_residue

    def outside_corners(self):
        r"""
        Return a list of the outside corners of the partition ``self``.

        An outside corner (also called a cocorner) of a partition
        `\lambda` is a cell on `\ZZ^2` which does not belong to
        the Young diagram of `\lambda` but can be added to this Young
        diagram to still form a straight-shape Young diagram.

        The entries of the list returned are pairs of the form `(i,j)`,
        where `i` and `j` are the coordinates of the respective corner.
        The coordinates are counted from `0`.

        .. NOTE::

            These are called "outer corners" in [Sag2001]_.

        EXAMPLES::

            sage: Partition([2,2,1]).outside_corners()
            [(0, 2), (2, 1), (3, 0)]
            sage: Partition([2,2]).outside_corners()
            [(0, 2), (2, 0)]
            sage: Partition([6,3,3,1,1,1]).outside_corners()
            [(0, 6), (1, 3), (3, 1), (6, 0)]
            sage: Partition([]).outside_corners()
            [(0, 0)]
        """
        p = self._list
        if not p:
            return [(0, 0)]
        res = [(0, p[0])]
        res.extend((n, j) for n, (i, j) in enumerate(zip(p[:-1], p[1:]), start=1) if i != j)
        res.append((len(p), 0))
        return res

    addable_cells = outside_corners   # for compatibility with partition tuples

    def outside_corners_residue(self, i, l):
        r"""
        Return a list of the outside corners of the partition ``self``
        having ``l``-residue ``i``.

        An outside corner (also called a cocorner) of a partition
        `\lambda` is a cell on `\ZZ^2` which does not belong to
        the Young diagram of `\lambda` but can be added to this Young
        diagram to still form a straight-shape Young diagram. See
        :meth:`residue` for the definition of the ``l``-residue.

        The entries of the list returned are pairs of the form `(i,j)`,
        where `i` and `j` are the coordinates of the respective corner.
        The coordinates are counted from `0`.

        EXAMPLES::

            sage: Partition([3,2,1]).outside_corners_residue(0, 3)
            [(0, 3), (3, 0)]
            sage: Partition([3,2,1]).outside_corners_residue(1, 3)
            [(1, 2)]
            sage: Partition([3,2,1]).outside_corners_residue(2, 3)
            [(2, 1)]
        """
        return [x for x in self.outside_corners() if self.residue(*x, l=l) == i]

    addable_cells_residue = outside_corners_residue

    def rim(self):
        r"""
        Return the rim of ``self``.

        The rim of a partition `\lambda` is defined as the cells which belong
        to `\lambda` and which are adjacent to cells not in `\lambda`.

        EXAMPLES:

        The rim of the partition `[5,5,2,1]` consists of the cells marked with
        ``#`` below::

            ****#
            *####
            ##
            #

            sage: Partition([5,5,2,1]).rim()
            [(3, 0), (2, 0), (2, 1), (1, 1), (1, 2), (1, 3), (1, 4), (0, 4)]

            sage: Partition([2,2,1]).rim()
            [(2, 0), (1, 0), (1, 1), (0, 1)]
            sage: Partition([2,2]).rim()
            [(1, 0), (1, 1), (0, 1)]
            sage: Partition([6,3,3,1,1]).rim()
            [(4, 0), (3, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2), (0, 3), (0, 4), (0, 5)]
            sage: Partition([]).rim()
            []
        """
        p = self
        res = []
        prevLen = 1
        for i in range(len(p) - 1, -1, -1):
            res.extend((i, c) for c in range(prevLen - 1, p[i]))
            prevLen = p[i]
        return res

    def outer_rim(self):
        r"""
        Return the outer rim of ``self``.

        The outer rim of a partition `\lambda` is defined as the cells which do
        not belong to `\lambda` and which are adjacent to cells in `\lambda`.

        EXAMPLES:

        The outer rim of the partition `[4,1]` consists of the cells marked
        with ``#`` below::

            ****#
            *####
            ##

        ::

            sage: Partition([4,1]).outer_rim()
            [(2, 0), (2, 1), (1, 1), (1, 2), (1, 3), (1, 4), (0, 4)]

            sage: Partition([2,2,1]).outer_rim()
            [(3, 0), (3, 1), (2, 1), (2, 2), (1, 2), (0, 2)]
            sage: Partition([2,2]).outer_rim()
            [(2, 0), (2, 1), (2, 2), (1, 2), (0, 2)]
            sage: Partition([6,3,3,1,1]).outer_rim()
            [(5, 0), (5, 1), (4, 1), (3, 1), (3, 2), (3, 3), (2, 3), (1, 3), (1, 4), (1, 5), (1, 6), (0, 6)]
            sage: Partition([]).outer_rim()
            [(0, 0)]
        """
        p = self
        res = []
        prevLen = 0
        for i in range(len(p) - 1, -1, -1):
            res.extend((i + 1, c) for c in range(prevLen, p[i] + 1))
            prevLen = p[i]
        res.append((0, prevLen))
        return res

    def zero_one_sequence(self):
        r"""
        Compute the finite `0-1` sequence of the partition.

        The full `0-1` sequence is the sequence (infinite in both
        directions) indicating the steps taken when following the
        outer rim of the diagram of the partition. We use the convention
        that in English convention, a 1 corresponds to an East step, and
        a 0 corresponds to a North step.

        Note that every full `0-1` sequence starts with infinitely many 0s and
        ends with infinitely many 1s.

        One place where these arise is in the affine symmetric group where
        one takes an affine permutation `w` and every `i` such that
        `w(i) \leq 0` corresponds to a 1 and `w(i) > 0` corresponds to a 0.
        See pages 24-25 of [LLMSSZ2013]_ for connections to affine Grassmannian
        elements (note there they use the French convention for their
        partitions).

        These are also known as **path sequences**, **Maya diagrams**,
        **plus-minus diagrams**, **Comet code** [Sta-EC2]_, among others.

        OUTPUT:

        The finite `0-1` sequence is obtained from the full `0-1`
        sequence by omitting all heading 0s and trailing 1s. The
        output sequence is finite, starts with a 1 and ends with a
        0 (unless it is empty, for the empty partition). Its length
        is the sum of the first part of the partition with the
        length of the partition.

        EXAMPLES::

            sage: Partition([5,4]).zero_one_sequence()
            [1, 1, 1, 1, 0, 1, 0]
            sage: Partition([]).zero_one_sequence()
            []
            sage: Partition([2]).zero_one_sequence()
            [1, 1, 0]

        TESTS::

            sage: all(Partitions().from_zero_one(mu.zero_one_sequence()) == mu for n in range(10) for mu in Partitions(n))
            True
        """
        tmp = {si - i for i, si in enumerate(self)}
        return [Integer(i not in tmp)
                for i in range(-len(self) + 1, self.get_part(0) + 1)]

    def core(self, length):
        r"""
        Return the ``length``-core of the partition -- in the literature
        the core is commonly referred to as the `k`-core, `p`-core,
        `r`-core, ... .

        The `r`-core of a partition `\lambda` can be obtained by
        repeatedly removing rim hooks of size `r` from (the Young diagram
        of) `\lambda` until this is no longer possible. The remaining
        partition is the core.

        EXAMPLES::

            sage: Partition([6,3,2,2]).core(3)
            [2, 1, 1]
            sage: Partition([]).core(3)
            []
            sage: Partition([8,7,7,4,1,1,1,1,1]).core(3)
            [2, 1, 1]

        TESTS::

            sage: Partition([3,3,3,2,1]).core(3)
            []
            sage: Partition([10,8,7,7]).core(4)
            []
            sage: Partition([21,15,15,9,6,6,6,3,3]).core(3)
            []
        """
        p = self
        # Normalize the length
        remainder = len(p) % length
        part = p[:] + [0]*remainder

        # Add the canonical vector to the partition
        part = [part[i-1] + len(part)-i for i in range(1, len(part)+1)]

        for e in range(length):
            k = e
            for i in reversed(range(1, len(part)+1)):
                if part[i-1] % length == e:
                    part[i-1] = k
                    k += length
        part.sort()
        part.reverse()

        # Remove the canonical vector
        part = [part[i-1]-len(part)+i for i in range(1, len(part)+1)]
        # Select the r-core
        return Partition([x for x in part if x != 0])

    def quotient(self, length):
        r"""
        Return the quotient of the partition  -- in the literature the
        quotient is commonly referred to as the `k`-quotient, `p`-quotient,
        `r`-quotient, ... .

        The `r`-quotient of a partition `\lambda` is a list of `r`
        partitions (labelled from `0` to `r-1`), constructed in the following
        way. Label each cell in the Young diagram of `\lambda` with its
        content modulo `r`. Let `R_i` be the set of rows ending in a cell
        labelled `i`, and `C_i` be the set of columns ending in a cell
        labelled `i`. Then the `j`-th component of the quotient of
        `\lambda` is the partition defined by intersecting `R_j` with
        `C_{j+1}`. (See Theorem 2.7.37 in [JK1981]_.)

        EXAMPLES::

            sage: Partition([7,7,5,3,3,3,1]).quotient(3)
            ([2], [1], [2, 2, 2])

        TESTS::

            sage: Partition([8,7,7,4,1,1,1,1,1]).quotient(3)
            ([2, 1], [2, 2], [2])
            sage: Partition([10,8,7,7]).quotient(4)
            ([2], [3], [2], [1])
            sage: Partition([6,3,3]).quotient(3)
            ([1], [1], [2])
            sage: Partition([3,3,3,2,1]).quotient(3)
            ([1], [1, 1], [1])
            sage: Partition([6,6,6,3,3,3]).quotient(3)
            ([2, 1], [2, 1], [2, 1])
            sage: Partition([21,15,15,9,6,6,6,3,3]).quotient(3)
            ([5, 2, 1], [5, 2, 1], [7, 3, 2])
            sage: Partition([21,15,15,9,6,6,3,3]).quotient(3)
            ([5, 2], [5, 2, 1], [7, 3, 1])
            sage: Partition([14,12,11,10,10,10,10,9,6,4,3,3,2,1]).quotient(5)
            ([3, 3], [2, 2, 1], [], [3, 3, 3], [1])

            sage: all(p == Partition(core=p.core(k), quotient=p.quotient(k))
            ....:     for i in range(10) for p in Partitions(i)
            ....:     for k in range(1,6))
            True
        """
        p = self
        # Normalize the length
        remainder = len(p) % length
        part = p[:] + [0]*(length-remainder)

        # Add the canonical vector to the partition
        part = [part[i-1] + len(part)-i for i in range(1, len(part)+1)]
        result = [None]*length

        # Reducing vector
        for e in range(length):
            k = e
            tmp = []
            for i in reversed(range(len(part))):
                if part[i] % length == e:
                    tmp.append(ZZ((part[i]-k)//length))
                    k += length

            a = [i for i in tmp if i != 0]
            a.reverse()
            result[e] = a

        from .partition_tuple import PartitionTuple
        return PartitionTuple(result)  # tuple(map(Partition, result))

    def is_core(self, k):
        r"""
        Return ``True`` if the Partition ``self`` is a ``k``-core.

        A partition is said to be a *`k`-core* if it has no hooks of length
        `k`. Equivalently, a partition is said to be a `k`-core if it is its
        own `k`-core (where the latter is defined as in :meth:`core`).

        Visually, this can be checked by trying to remove border strips of size
        `k` from ``self``.  If this is not possible, then ``self`` is a
        `k`-core.

        EXAMPLES:

        In the partition (2, 1), a hook length of 2 does not occur, but a hook
        length of 3 does::

            sage: p = Partition([2, 1])
            sage: p.is_core(2)
            True
            sage: p.is_core(3)
            False

            sage: q = Partition([12, 8, 5, 5, 2, 2, 1])
            sage: q.is_core(4)
            False
            sage: q.is_core(5)
            True
            sage: q.is_core(0)
            True

        .. SEEALSO::

            :meth:`core`, :class:`Core`
        """
        return k not in self.hooks()

    def k_interior(self, k):
        r"""
        Return the partition consisting of the cells of ``self`` whose hook
        lengths are greater than ``k``.

        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p.hook_lengths()
            [[5, 3, 1], [3, 1], [1]]
            sage: p.k_interior(2)
            [2, 1]
            sage: p.k_interior(3)
            [1]

            sage: p = Partition([])
            sage: p.k_interior(3)
            []
        """
        return Partition([len([i for i in row if i > k])
                          for row in self.hook_lengths()])

    def k_boundary(self, k):
        r"""
        Return the skew partition formed by removing the cells of the
        ``k``-interior, see :meth:`k_interior`.

        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p.k_boundary(2)
            [3, 2, 1] / [2, 1]
            sage: p.k_boundary(3)
            [3, 2, 1] / [1]

            sage: p = Partition([12,8,5,5,2,2,1])
            sage: p.k_boundary(4)
            [12, 8, 5, 5, 2, 2, 1] / [8, 5, 2, 2]
        """
        return SkewPartition([self, self.k_interior(k)])

    def add_cell(self, i, j=None):
        r"""
        Return a partition corresponding to ``self`` with a cell added in
        row ``i``. (This does not change ``self``.)

        EXAMPLES::

            sage: Partition([3, 2, 1, 1]).add_cell(0)
            [4, 2, 1, 1]
            sage: cell = [4, 0]; Partition([3, 2, 1, 1]).add_cell(*cell)
            [3, 2, 1, 1, 1]
        """

        if j is None:
            if i >= len(self):
                j = 0
            else:
                j = self[i]

        if (i, j) in self.outside_corners():
            pl = self.to_list()
            if i == len(pl):
                pl.append(1)
            else:
                pl[i] += 1
            return Partition(pl)

        raise ValueError(f"[{i}, {j}] is not an addable cell")

    def remove_cell(self, i, j=None):
        """
        Return the partition obtained by removing a cell at the end of row
        ``i`` of ``self``.

        EXAMPLES::

            sage: Partition([2,2]).remove_cell(1)
            [2, 1]
            sage: Partition([2,2,1]).remove_cell(2)
            [2, 2]
            sage: #Partition([2,2]).remove_cell(0)

        ::

            sage: Partition([2,2]).remove_cell(1,1)
            [2, 1]
            sage: #Partition([2,2]).remove_cell(1,0)
        """

        if i >= len(self):
            raise ValueError("i must be less than the length of the partition")

        if j is None:
            j = self[i] - 1

        if (i, j) not in self.corners():
            raise ValueError("[%d,%d] is not a corner of the partition" % (i, j))

        if self[i] == 1:
            return Partition(self[:-1])
        else:
            return Partition(self[:i] + [self[i:i+1][0] - 1] + self[i+1:])

    def k_irreducible(self, k):
        r"""
        Return the partition with all `r \times (k+1-r)` rectangles removed.

        If ``self`` is a `k`-bounded partition, then this method will return the partition
        where all rectangles of dimension `r \times (k+1-r)` for `1 \leq r \leq k`
        have been deleted.

        If ``self`` is not a `k`-bounded partition then the method will raise an error.

        INPUT:

        - ``k`` -- nonnegative integer

        OUTPUT: a partition

        EXAMPLES::

            sage: Partition([3,2,2,1,1,1]).k_irreducible(4)
            [3, 2, 2, 1, 1, 1]
            sage: Partition([3,2,2,1,1,1]).k_irreducible(3)
            []
            sage: Partition([3,3,3,2,2,2,2,2,1,1,1,1]).k_irreducible(3)
            [2, 1]
        """
        pexp = self.to_exp()
        return Partition(sum(([r+1] for r in range(len(pexp)-1, -1, -1) for m in range(pexp[r] % (k-r))), []))

    def k_skew(self, k):
        r"""
        Return the `k`-skew partition.

        The `k`-skew diagram of a `k`-bounded partition is the skew diagram
        denoted `\lambda/^k` satisfying the conditions:

        1. row `i` of `\lambda/^k` has length `\lambda_i`,

        2. no cell in `\lambda/^k` has hook-length exceeding `k`,

        3. every square above the diagram of `\lambda/^k` has hook
           length exceeding `k`.

        REFERENCES:

        - [LM2004]_

        EXAMPLES::

            sage: p = Partition([4,3,2,2,1,1])
            sage: p.k_skew(4)
            [9, 5, 3, 2, 1, 1] / [5, 2, 1]
        """

        if len(self) == 0:
            return SkewPartition([[], []])

        if self[0] > k:
            raise ValueError(f"the partition must be {k}-bounded")

        # Find the k-skew diagram of the partition formed
        # by removing the first row
        s = Partition(self[1:]).k_skew(k)

        s_inner = list(s.inner())
        s_outer = list(s.outer())
        s_conj_rl = s.conjugate().row_lengths()

        # Find the leftmost column with less than
        # or equal to kdiff cells
        kdiff = k - self[0]

        if s_outer == []:
            spot = 0
        else:
            spot = s_outer[0]

        for i in range(len(s_conj_rl)):
            if s_conj_rl[i] <= kdiff:
                spot = i
                break

        outer = [self[0] + spot] + s_outer[:]
        if spot > 0:
            inner = [spot] + s_inner[:]
        else:
            inner = s_inner[:]

        return SkewPartition([outer, inner])

    def to_core(self, k):
        r"""
        Map the `k`-bounded partition ``self`` to its corresponding `k+1`-core.

        See also :meth:`k_skew`.

        EXAMPLES::

            sage: p = Partition([4,3,2,2,1,1])
            sage: c = p.to_core(4); c
            [9, 5, 3, 2, 1, 1]
            sage: type(c)
            <class 'sage.combinat.core.Cores_length_with_category.element_class'>
            sage: c.to_bounded_partition() == p
            True
        """
        from sage.combinat.core import Core
        return Core(self.k_skew(k)[0], k+1)

    def from_kbounded_to_reduced_word(self, k):
        r"""
        Map a `k`-bounded partition to a reduced word for an element in
        the affine permutation group.

        This uses the fact that there is a bijection between `k`-bounded
        partitions and `(k+1)`-cores and an action of the affine nilCoxeter
        algebra of type `A_k^{(1)}` on `(k+1)`-cores as described in [LM2006b]_.

        EXAMPLES::

            sage: p = Partition([2,1,1])
            sage: p.from_kbounded_to_reduced_word(2)
            [2, 1, 2, 0]
            sage: p = Partition([3,1])
            sage: p.from_kbounded_to_reduced_word(3)
            [3, 2, 1, 0]
            sage: p.from_kbounded_to_reduced_word(2)
            Traceback (most recent call last):
            ...
            ValueError: the partition must be 2-bounded
            sage: p = Partition([])
            sage: p.from_kbounded_to_reduced_word(2)
            []
        """
        p = self.k_skew(k)[0]
        result = []
        while not p.is_empty():
            corners = p.corners()
            c = p.content(corners[0][0], corners[0][1]) % (k+1)
            result.append(Integer(c))
            list = [x for x in corners if p.content(x[0], x[1]) % (k+1) == c]
            for x in list:
                p = p.remove_cell(x[0])
        return result

    def from_kbounded_to_grassmannian(self, k):
        r"""
        Map a `k`-bounded partition to a Grassmannian element in
        the affine Weyl group of type `A_k^{(1)}`.

        For details, see the documentation of the method
        :meth:`from_kbounded_to_reduced_word` .

        EXAMPLES::

            sage: p = Partition([2,1,1])
            sage: p.from_kbounded_to_grassmannian(2)                                    # needs sage.modules
            [-1  1  1]
            [-2  2  1]
            [-2  1  2]
            sage: p = Partition([])
            sage: p.from_kbounded_to_grassmannian(2)                                    # needs sage.modules
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return WeylGroup(['A', k, 1]).from_reduced_word(self.from_kbounded_to_reduced_word(k))

    def to_list(self):
        r"""
        Return ``self`` as a list.

        EXAMPLES::

            sage: p = Partition([2,1]).to_list(); p
            [2, 1]
            sage: type(p)
            <class 'list'>

        TESTS::

            sage: p = Partition([2,1])
            sage: pl = p.to_list()
            sage: pl[0] = 0; p
            [2, 1]
        """
        return self._list[:]

    def add_vertical_border_strip(self, k):
        """
        Return a list of all the partitions that can be obtained by adding
        a vertical border strip of length ``k`` to ``self``.

        EXAMPLES::

            sage: Partition([]).add_vertical_border_strip(0)
            [[]]
            sage: Partition([3,2,1]).add_vertical_border_strip(0)
            [[3, 2, 1]]
            sage: Partition([]).add_vertical_border_strip(2)
            [[1, 1]]
            sage: Partition([2,2]).add_vertical_border_strip(2)
            [[3, 3], [3, 2, 1], [2, 2, 1, 1]]
            sage: Partition([3,2,2]).add_vertical_border_strip(2)
            [[4, 3, 2], [4, 2, 2, 1], [3, 3, 3], [3, 3, 2, 1], [3, 2, 2, 1, 1]]
        """
        if k == 0:
            return [self]

        shelf = []
        res = []
        i = 0
        ell = len(self._list)
        while i < ell:
            tmp = 1
            while i+1 < ell and self._list[i] == self._list[i+1]:
                tmp += 1
                i += 1
            if i == ell-1 and i > 0 and self._list[i] != self._list[i-1]:
                tmp = 1
            shelf.append(tmp)
            i += 1

        # added the last shelf on the right side of
        # the first line
        shelf.append(k)

        # list all of the positions for cells
        # filling each self from the left to the right
        for iv in IntegerListsBackend_invlex(k, length=len(shelf), ceiling=shelf, check=False)._iter():
            tmp = self._list + [0]*k
            j = 0
            for t in range(len(iv)):
                for _ in range(iv[t]):
                    tmp[j] += 1
                    j += 1
                j = sum(shelf[:t+1])
            # This should never return the empty partition.
            # So tmp should never be [0, ..., 0].
            while not tmp[-1]:
                tmp.pop()
            res.append(_Partitions(tmp))
        return res

    def add_horizontal_border_strip(self, k):
        """
        Return a list of all the partitions that can be obtained by adding
        a horizontal border strip of length ``k`` to ``self``.

        EXAMPLES::

            sage: Partition([]).add_horizontal_border_strip(0)
            [[]]
            sage: Partition([3,2,1]).add_horizontal_border_strip(0)
            [[3, 2, 1]]
            sage: Partition([]).add_horizontal_border_strip(2)
            [[2]]
            sage: Partition([2,2]).add_horizontal_border_strip(2)
            [[4, 2], [3, 2, 1], [2, 2, 2]]
            sage: Partition([3,2,2]).add_horizontal_border_strip(2)
            [[5, 2, 2], [4, 3, 2], [4, 2, 2, 1], [3, 3, 2, 1], [3, 2, 2, 2]]
        """
        if k == 0:
            return [self]

        L = self._list
        res = []
        mapping = [0]
        shelf = [k]
        for i in range(len(L)-1):
            val = L[i] - L[i+1]
            if not val:
                continue
            mapping.append(i+1)
            shelf.append(val)

        # add the last shelf
        if L:
            mapping.append(len(L))
            shelf.append(L[-1])

        # list all of the positions for cells
        # filling each self from the top to bottom
        for iv in IntegerListsBackend_invlex(k, length=len(shelf), ceiling=shelf, check=False)._iter():
            tmp = self._list + [0]
            for i, val in enumerate(iv):
                if val:
                    tmp[mapping[i]] += val
            # Only the last row is possibly empty
            if not tmp[-1]:
                tmp.pop()
            res.append(_Partitions(tmp))
        return res

    def vertical_border_strip_cells(self, k):
        """
        Return a list of all the vertical border strips of length ``k``
        which can be added to ``self``, where each horizontal border strip is
        a ``generator`` of cells.

        EXAMPLES::

            sage: list(Partition([]).vertical_border_strip_cells(0))
            []
            sage: list(Partition([3,2,1]).vertical_border_strip_cells(0))
            []
            sage: list(Partition([]).vertical_border_strip_cells(2))
            [[(0, 0), (1, 0)]]
            sage: list(Partition([2,2]).vertical_border_strip_cells(2))
            [[(0, 2), (1, 2)],
             [(0, 2), (2, 0)],
             [(2, 0), (3, 0)]]
            sage: list(Partition([3,2,2]).vertical_border_strip_cells(2))
            [[(0, 3), (1, 2)],
             [(0, 3), (3, 0)],
             [(1, 2), (2, 2)],
             [(1, 2), (3, 0)],
             [(3, 0), (4, 0)]]
        """
        if k == 0:
            return []

        shelf = []
        i = 0
        ell = len(self._list)
        while i < ell:
            tmp = 1
            while i+1 < ell and self._list[i] == self._list[i+1]:
                tmp += 1
                i += 1
            if i == ell-1 and i > 0 and self._list[i] != self._list[i-1]:
                tmp = 1
            shelf.append(tmp)
            i += 1

        # added the last shelf on the right side of
        # the first line
        shelf.append(k)
        # list all of the positions for cells
        tmp = self._list + [0]*k
        for iv in IntegerListsBackend_invlex(k, length=len(shelf),
                                             ceiling=shelf,
                                             check=False)._iter():
            j = 0
            current_strip = []
            for t in range(len(iv)):
                for _ in range(iv[t]):
                    current_strip.append((j, tmp[j]))
                    j += 1
                j = sum(shelf[:t+1])
            yield current_strip

    def horizontal_border_strip_cells(self, k):
        """
        Return a list of all the horizontal border strips of length ``k``
        which can be added to ``self``, where each horizontal border strip is
        a ``generator`` of cells.

        EXAMPLES::

            sage: list(Partition([]).horizontal_border_strip_cells(0))
            []
            sage: list(Partition([3,2,1]).horizontal_border_strip_cells(0))
            []
            sage: list(Partition([]).horizontal_border_strip_cells(2))
            [[(0, 0), (0, 1)]]
            sage: list(Partition([2,2]).horizontal_border_strip_cells(2))
            [[(0, 2), (0, 3)], [(0, 2), (2, 0)], [(2, 0), (2, 1)]]
            sage: list(Partition([3,2,2]).horizontal_border_strip_cells(2))
            [[(0, 3), (0, 4)],
             [(0, 3), (1, 2)],
             [(0, 3), (3, 0)],
             [(1, 2), (3, 0)],
             [(3, 0), (3, 1)]]
        """
        if k == 0:
            return []

        L = self._list
        shelf = [k]  # the number of boxes which will fit in a row
        mapping = [0]  # a record of the rows
        for i in range(len(L)-1):
            val = L[i] - L[i+1]
            if not val:
                continue
            mapping.append(i+1)
            shelf.append(val)

        # add the last shelf
        if L:
            mapping.append(len(L))
            shelf.append(L[-1])

        L.append(0)  # add room on the bottom
        # list all of the positions for cells
        # filling each self from the top to bottom
        for iv in IntegerListsBackend_invlex(k, length=len(shelf), ceiling=shelf, check=False)._iter():
            tmp = []
            # mapping[i] is the row index, val is the number of cells added to the row.
            for i, val in enumerate(iv):
                tmp.extend((mapping[i], L[mapping[i]] + j) for j in range(val))
            yield tmp

    def remove_horizontal_border_strip(self, k):
        """
        Return the partitions obtained from ``self`` by removing an
        horizontal border strip of length ``k``.

        EXAMPLES::

            sage: Partition([5,3,1]).remove_horizontal_border_strip(0).list()
            [[5, 3, 1]]
            sage: Partition([5,3,1]).remove_horizontal_border_strip(1).list()
            [[5, 3], [5, 2, 1], [4, 3, 1]]
            sage: Partition([5,3,1]).remove_horizontal_border_strip(2).list()
            [[5, 2], [5, 1, 1], [4, 3], [4, 2, 1], [3, 3, 1]]
            sage: Partition([5,3,1]).remove_horizontal_border_strip(3).list()
            [[5, 1], [4, 2], [4, 1, 1], [3, 3], [3, 2, 1]]
            sage: Partition([5,3,1]).remove_horizontal_border_strip(4).list()
            [[4, 1], [3, 2], [3, 1, 1]]
            sage: Partition([5,3,1]).remove_horizontal_border_strip(5).list()
            [[3, 1]]
            sage: Partition([5,3,1]).remove_horizontal_border_strip(6).list()
            []

        The result is returned as an instance of
        :class:`Partitions_with_constraints`::

            sage: Partition([5,3,1]).remove_horizontal_border_strip(5)
            The subpartitions of [5, 3, 1] obtained by removing a horizontal border strip of length 5

        TESTS::

            sage: Partition([3,2,2]).remove_horizontal_border_strip(2).list()
            [[3, 2], [2, 2, 1]]
            sage: Partition([3,2,2]).remove_horizontal_border_strip(2).first().parent()
            The subpartitions of [3, 2, 2] obtained by removing a horizontal border strip of length 2
            sage: Partition([]).remove_horizontal_border_strip(0).list()
            [[]]
            sage: Partition([]).remove_horizontal_border_strip(6).list()
            []
        """
        return Partitions_with_constraints(n=self.size() - k,
                                           min_length=len(self) - 1,
                                           max_length=len(self),
                                           floor=self[1:] + [0],
                                           ceiling=self[:],
                                           max_slope=0,
                                           name=f"The subpartitions of {self} obtained by removing a horizontal border strip of length {k}")

    def k_conjugate(self, k):
        r"""
        Return the ``k``-conjugate of ``self``.

        The `k`-conjugate is the partition that is given by the columns of
        the `k`-skew diagram of the partition.

        We can also define the `k`-conjugate in the following way. Let `P`
        denote the bijection from `(k+1)`-cores to `k`-bounded partitions. The
        `k`-conjugate of a `(k+1)`-core `\lambda` is

        .. MATH::

            \lambda^{(k)} = P^{-1}\left( (P(\lambda))^{\prime} \right).

        EXAMPLES::

            sage: p = Partition([4,3,2,2,1,1])
            sage: p.k_conjugate(4)
            [3, 2, 2, 1, 1, 1, 1, 1, 1]
        """
        return Partition(self.k_skew(k).conjugate().row_lengths())

    def arms_legs_coeff(self, i, j):
        r"""
        This is a statistic on a cell `c = (i,j)` in the diagram of partition
        `p` given by

        .. MATH::

            \frac{ 1 - q^a \cdot t^{\ell + 1} }{ 1 - q^{a + 1} \cdot t^{\ell} }

        where `a` is the arm length of `c` and `\ell` is the leg length of `c`.

        The coordinates ``i`` and ``j`` of the cell are understood to be
        `0`-based, so that ``(0, 0)`` is the northwesternmost cell (in
        English notation).

        EXAMPLES::

            sage: Partition([3,2,1]).arms_legs_coeff(1,1)
            (-t + 1)/(-q + 1)
            sage: Partition([3,2,1]).arms_legs_coeff(0,0)
            (-q^2*t^3 + 1)/(-q^3*t^2 + 1)
            sage: Partition([3,2,1]).arms_legs_coeff(*[0,0])
            (-q^2*t^3 + 1)/(-q^3*t^2 + 1)
        """
        QQqt = PolynomialRing(QQ, ['q', 't'])
        (q, t) = QQqt.gens()
        if i < len(self) and j < self[i]:
            res = 1 - q**self.arm_length(i, j) * t**(self.leg_length(i, j)+1)
            res /= 1 - q**(self.arm_length(i, j)+1) * t**self.leg_length(i, j)
            return res
        return ZZ.one()

    def atom(self):
        """
        Return a list of the standard tableaux of size ``self.size()`` whose
        atom is equal to ``self``.

        EXAMPLES::

            sage: Partition([2,1]).atom()
            [[[1, 2], [3]]]
            sage: Partition([3,2,1]).atom()
            [[[1, 2, 3, 6], [4, 5]], [[1, 2, 3], [4, 5], [6]]]
        """
        return [tab for tab in tableau.StandardTableaux_size(self.size())
                if tab.atom() == self]

    def k_atom(self, k):
        r"""
        Return a list of the standard tableaux of size ``self.size()`` whose
        ``k``-atom is equal to ``self``.

        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p.k_atom(1)
            []
            sage: p.k_atom(3)
            [[[1, 1, 1, 2, 3], [2]],
             [[1, 1, 1, 3], [2, 2]],
             [[1, 1, 1, 2], [2], [3]],
             [[1, 1, 1], [2, 2], [3]]]
            sage: Partition([3,2,1]).k_atom(4)
            [[[1, 1, 1, 3], [2, 2]], [[1, 1, 1], [2, 2], [3]]]

        TESTS::

            sage: Partition([1]).k_atom(1)
            [[[1]]]
            sage: Partition([1]).k_atom(2)
            [[[1]]]
            sage: Partition([]).k_atom(1)
            [[]]
        """
        res = [tableau.Tableau([])]
        for i in range(len(self)):
            res = (x.promotion_operator(self[-i - 1]) for x in res)
            res = sum(res, [])
            res = (y.catabolism_projector(Partition(self[-i - 1:]).k_split(k))
                   for y in res)
            res = [i for i in res if i]
        return res

    def k_split(self, k):
        """
        Return the ``k``-split of ``self``.

        EXAMPLES::

            sage: Partition([4,3,2,1]).k_split(3)
            []
            sage: Partition([4,3,2,1]).k_split(4)
            [[4], [3, 2], [1]]
            sage: Partition([4,3,2,1]).k_split(5)
            [[4, 3], [2, 1]]
            sage: Partition([4,3,2,1]).k_split(6)
            [[4, 3, 2], [1]]
            sage: Partition([4,3,2,1]).k_split(7)
            [[4, 3, 2, 1]]
            sage: Partition([4,3,2,1]).k_split(8)
            [[4, 3, 2, 1]]
        """
        if self == []:
            return []
        elif k < self[0]:
            return []
        else:
            res = []
            part = list(self)
            while part and part[0] + len(part) - 1 >= k:
                p = k - part[0]
                res.append(part[:p + 1])
                part = part[p + 1:]
            if part:
                res.append(part)
        return res

    def jacobi_trudi(self):
        """
        Return the Jacobi-Trudi matrix of ``self`` thought of as a skew
        partition. See :meth:`SkewPartition.jacobi_trudi()
        <sage.combinat.skew_partition.SkewPartition.jacobi_trudi>`.

        EXAMPLES::

            sage: # needs sage.modules
            sage: part = Partition([3,2,1])
            sage: jt = part.jacobi_trudi(); jt
            [h[3] h[1]    0]
            [h[4] h[2]  h[]]
            [h[5] h[3] h[1]]
            sage: s = SymmetricFunctions(QQ).schur()
            sage: h = SymmetricFunctions(QQ).homogeneous()
            sage: h( s(part) )
            h[3, 2, 1] - h[3, 3] - h[4, 1, 1] + h[5, 1]
            sage: jt.det()
            h[3, 2, 1] - h[3, 3] - h[4, 1, 1] + h[5, 1]
        """
        return SkewPartition([self, []]).jacobi_trudi()

    def character_polynomial(self):
        r"""
        Return the character polynomial associated to the partition ``self``.

        The character polynomial `q_\mu` associated to a partition `\mu`
        is defined by

        .. MATH::

            q_\mu(x_1, x_2, \ldots, x_k) = \downarrow \sum_{\alpha \vdash k}
            \frac{ \chi^\mu_\alpha }{1^{a_1}2^{a_2}\cdots k^{a_k}a_1!a_2!\cdots
            a_k!} \prod_{i=1}^{k} (ix_i-1)^{a_i}

        where `k` is the size of `\mu`, and `a_i` is the multiplicity of
        `i` in `\alpha`.

        It is computed in the following manner:

        1. Expand the Schur function `s_\mu` in the power-sum basis,

        2. Replace each `p_i` with `ix_i-1`,

        3. Apply the umbral operator `\downarrow` to the resulting polynomial.

        EXAMPLES::

            sage: Partition([1]).character_polynomial()                                 # needs sage.modules
            x - 1
            sage: Partition([1,1]).character_polynomial()                               # needs sage.modules
            1/2*x0^2 - 3/2*x0 - x1 + 1
            sage: Partition([2,1]).character_polynomial()                               # needs sage.modules
            1/3*x0^3 - 2*x0^2 + 8/3*x0 - x2
        """
        # Create the polynomial ring we will use
        k = self.size()
        P = PolynomialRing(QQ, k, 'x')
        x = P.gens()

        # Expand s_mu in the power sum basis
        from sage.combinat.sf.sf import SymmetricFunctions
        Sym = SymmetricFunctions(QQ)
        s = Sym.schur()
        p = Sym.power()
        ps_mu = p(s(self))

        # Replace each p_i by i*x_i-1
        items = ps_mu.monomial_coefficients().items()  # items contains a list of (partition, coeff) pairs

        def partition_to_monomial(part):
            return prod([i*x[i-1] - 1 for i in part])

        res = [[partition_to_monomial(mc[0]), mc[1]] for mc in items]

        # Write things in the monomial basis
        res = [prod(pair) for pair in res]
        res = sum(res)

        # Apply the umbral operator and return the result
        from sage.combinat.misc import umbral_operation
        return umbral_operation(res)

    def dimension(self, smaller=None, k=1):
        r"""
        Return the number of paths from the ``smaller`` partition to
        the partition ``self``, where each step consists of adding a
        `k`-ribbon while keeping a partition.

        Note that a 1-ribbon is just a single cell, so this counts paths
        in the Young graph when `k = 1`.

        Note also that the default case (`k = 1` and ``smaller = []``)
        gives the dimension of the irreducible representation of the
        symmetric group corresponding to ``self``.

        INPUT:

        - ``smaller`` -- a partition (default: an empty list ``[]``)

        - ``k`` -- positive integer (default: 1)

        OUTPUT: the number of such paths

        EXAMPLES:

        Looks at the number of ways of getting from ``[5,4]`` to the empty
        partition, removing one cell at a time::

            sage: mu = Partition([5,4])
            sage: mu.dimension()
            42

        Same, but removing one 3-ribbon at a time. Note that the 3-core of
        ``mu`` is empty::

            sage: mu.dimension(k=3)
            3

        The 2-core of ``mu`` is not the empty partition::

            sage: mu.dimension(k=2)
            0

        Indeed, the 2-core of ``mu`` is ``[1]``::

            sage: mu.dimension(Partition([1]),k=2)
            2

        TESTS:

        Checks that the sum of squares of dimensions of characters of the
        symmetric group is the order of the group::

            sage: all(sum(mu.dimension()^2 for mu in Partitions(i)) == factorial(i)
            ....:     for i in range(10))
            True

        A check coming from the theory of `k`-differentiable posets::

            sage: k = 2; core = Partition([2,1])
            sage: all(sum(mu.dimension(core,k=2)^2
            ....:         for mu in Partitions(3+i*2) if mu.core(2) == core)
            ....:     == 2^i*factorial(i) for i in range(10))
            True

        Checks that the dimension satisfies the obvious recursion relation::

            sage: test = lambda larger, smaller: larger.dimension(smaller) == sum(mu.dimension(smaller) for mu in larger.down())
            sage: all(test(larger,smaller) for l in range(1,8) for s in range(8)
            ....:     for larger in Partitions(l) for smaller in Partitions(s) if smaller != larger)
            True

        ALGORITHM:

        Depending on the parameters given, different simplifications
        occur. When `k=1` and ``smaller`` is empty, this function uses
        the hook formula. When `k=1` and ``smaller`` is not empty, it
        uses a formula from [ORV]_.

        When `k \neq 1`, we first check that both ``self`` and
        ``smaller`` have the same `k`-core, then use the `k`-quotients
        and the same algorithm on each of the `k`-quotients.

        AUTHORS:

        - Paul-Olivier Dehaye (2011-06-07)
        """
        larger = self
        if smaller is None:
            smaller = Partition([])
        if k == 1:
            if smaller == Partition([]):        # In this case, use the hook dimension formula
                return factorial(larger.size()) / prod(larger.hooks())
            if not larger.contains(smaller):    # easy case
                return 0

            # relative dimension
            # Uses a formula of Olshanski, Regev, Vershik (see reference)
            def inv_factorial(i):
                if i < 0:
                    return 0
                return 1 / factorial(i)

            len_range = range(larger.length())
            from sage.matrix.constructor import matrix
            M = matrix(QQ, [[inv_factorial(larger.get_part(i) - smaller.get_part(j) - i + j)
                             for i in len_range] for j in len_range])
            return factorial(larger.size() - smaller.size()) * M.determinant()

        larger_core = larger.core(k)
        smaller_core = smaller.core(k)
        if smaller_core != larger_core:  # easy case
            return 0
        larger_quotients = larger.quotient(k)
        smaller_quotients = smaller.quotient(k)

        def multinomial_with_partitions(sizes, path_counts):
            # count the number of ways of performing the k paths in parallel,
            # if we know the total length allotted for each of the paths (sizes), and the number
            # of paths for each component. A multinomial picks the ordering of the components where
            # each step is taken.
            return prod(path_counts) * multinomial(sizes)

        sizes = [larger_quotients[i].size() - smaller_quotients[i].size() for i in range(k)]
        path_counts = [larger_quotients[i].dimension(smaller_quotients[i]) for i in range(k)]
        return multinomial_with_partitions(sizes, path_counts)

    def plancherel_measure(self):
        r"""
        Return the probability of ``self`` under the Plancherel probability
        measure on partitions of the same size.

        This probability distribution comes from the uniform distribution
        on permutations via the Robinson-Schensted correspondence.

        See :wikipedia:`Plancherel\_measure`
        and :meth:`Partitions_n.random_element_plancherel`.

        EXAMPLES::

            sage: Partition([]).plancherel_measure()
            1
            sage: Partition([1]).plancherel_measure()
            1
            sage: Partition([2]).plancherel_measure()
            1/2
            sage: [mu.plancherel_measure() for mu in Partitions(3)]
            [1/6, 2/3, 1/6]
            sage: Partition([5,4]).plancherel_measure()
            7/1440

        TESTS::

            sage: all(sum(mu.plancherel_measure() for mu in Partitions(n))==1 for n in range(10))
            True
        """
        return self.dimension()**2 / factorial(self.size())

    def outline(self, variable=None):
        r"""
        Return the outline of the partition ``self``.

        This is a piecewise linear function, normalized so that the area
        under the partition ``[1]`` is 2.

        INPUT:

        - ``variable`` -- a variable (default: ``'x'`` in the symbolic ring)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: [Partition([5,4]).outline()(x=i) for i in range(-10, 11)]
            [10, 9, 8, 7, 6, 5, 6, 5, 6, 5, 4, 3, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: Partition([]).outline()
            abs(x)
            sage: Partition([1]).outline()
            abs(x + 1) + abs(x - 1) - abs(x)
            sage: y = SR.var("y")
            sage: Partition([6,5,1]).outline(variable=y)
            abs(y + 6) - abs(y + 5) + abs(y + 4) - abs(y + 3)
             + abs(y - 1) - abs(y - 2) + abs(y - 3)

        TESTS::

            sage: integrate(Partition([1]).outline() - abs(x), (x, -10, 10))              # needs sage.symbolic
            2
        """
        if variable is None:
            variable = var('x')
        outside_contents = [self.content(*c) for c in self.outside_corners()]
        inside_contents = [self.content(*c) for c in self.corners()]
        return sum(abs(variable+c) for c in outside_contents)\
            - sum(abs(variable+c) for c in inside_contents)

    def dual_equivalence_graph(self, directed=False, coloring=None):
        r"""
        Return the dual equivalence graph of ``self``.

        Two permutations `p` and `q` in the symmetric group `S_n`
        differ by an `i`-*elementary dual equivalence (or dual Knuth)
        relation* (where `i` is an integer with `1 < i < n`) when the
        following two conditions are satisfied:

        - In the one-line notation of the permutation `p`, the letter
          `i` does not appear between `i-1` and `i+1`.

        - The permutation `q` is obtained from `p` by switching two
          of the three letters `i-1, i, i+1` (in its one-line
          notation) -- namely, the leftmost and the rightmost one
          in order of their appearance in `p`.

        Notice that this is equivalent to the statement that the
        permutations `p^{-1}` and `q^{-1}` differ by an elementary
        Knuth equivalence at positions `i-1, i, i+1`.

        Two standard Young tableaux of shape `\lambda` differ by an
        `i`-elementary dual equivalence relation (of color `i`), if
        their reading words differ by an `i`-elementary dual
        equivalence relation.

        The *dual equivalence graph* of the partition `\lambda` is the
        edge-colored graph whose vertices are the standard Young
        tableaux of shape `\lambda`, and whose edges colored by `i` are
        given by the `i`-elementary dual equivalences.

        INPUT:

        - ``directed`` -- boolean (default: ``False``); whether to have the
          dual equivalence graph be directed (where we have a directed edge
          `S \to T` if `i` appears to the left of `i+1` in the reading word of
          `T`; otherwise we have the directed edge `T \to S`)

        - ``coloring`` -- (optional) a function which sends each
          integer `i > 1` to a color (as a string, e.g., ``'red'`` or
          ``'black'``) to be used when visually representing the
          resulting graph using dot2tex; the default choice is
          ``2 -> 'red', 3 -> 'blue', 4 -> 'green', 5 -> 'purple',
          6 -> 'brown', 7 -> 'orange', 8 -> 'yellow', anything greater
          than 8 -> 'black'``.

        REFERENCES:

        - [As2008b]_

        EXAMPLES::

            sage: # needs sage.graphs
            sage: P = Partition([3,1,1])
            sage: G = P.dual_equivalence_graph()
            sage: G.edges(sort=True)
            [([[1, 2, 3], [4], [5]], [[1, 2, 4], [3], [5]], 3),
             ([[1, 2, 4], [3], [5]], [[1, 2, 5], [3], [4]], 4),
             ([[1, 2, 4], [3], [5]], [[1, 3, 4], [2], [5]], 2),
             ([[1, 2, 5], [3], [4]], [[1, 3, 5], [2], [4]], 2),
             ([[1, 3, 4], [2], [5]], [[1, 3, 5], [2], [4]], 4),
             ([[1, 3, 5], [2], [4]], [[1, 4, 5], [2], [3]], 3)]
            sage: G = P.dual_equivalence_graph(directed=True)
            sage: G.edges(sort=True)
            [([[1, 2, 4], [3], [5]], [[1, 2, 3], [4], [5]], 3),
             ([[1, 2, 5], [3], [4]], [[1, 2, 4], [3], [5]], 4),
             ([[1, 3, 4], [2], [5]], [[1, 2, 4], [3], [5]], 2),
             ([[1, 3, 5], [2], [4]], [[1, 2, 5], [3], [4]], 2),
             ([[1, 3, 5], [2], [4]], [[1, 3, 4], [2], [5]], 4),
             ([[1, 4, 5], [2], [3]], [[1, 3, 5], [2], [4]], 3)]

        TESTS::

            sage: # needs sage.graphs
            sage: G = Partition([1]).dual_equivalence_graph()
            sage: G.vertices(sort=False)
            [[[1]]]
            sage: G = Partition([]).dual_equivalence_graph()
            sage: G.vertices(sort=False)
            [[]]
            sage: P = Partition([3,1,1])
            sage: G = P.dual_equivalence_graph(coloring=lambda x: 'red')
            sage: G2 = P.dual_equivalence_graph(coloring={2: 'black', 3: 'blue',
            ....:                                         4: 'cyan', 5: 'grey'})
            sage: G is G2
            False
            sage: G == G2
            True
        """
        # We do some custom caching to not recreate the graph, but to make
        #   copies with the desired coloring (i.e., act like a factory).
        try:
            if directed:
                G = self._DDEG.copy(immutable=False)
            else:
                G = self._DEG.copy(immutable=False)

            try:
                from sage.graphs.dot2tex_utils import have_dot2tex
                have = have_dot2tex()
            except ImportError:
                have = False

            if have:
                if coloring is None:
                    d = {2: 'red', 3: 'blue', 4: 'green', 5: 'purple',
                         6: 'brown', 7: 'orange', 8: 'yellow'}

                    def coloring(i):
                        if i in d:
                            return d[i]
                        return 'black'
                elif isinstance(coloring, dict):
                    d = coloring
                    coloring = lambda x: d[x]
                G.set_latex_options(format='dot2tex',
                                    edge_labels=True,
                                    color_by_label=coloring)
            return G
        except AttributeError:
            pass

        T = list(tableau.StandardTableaux(self))
        n = sum(self)
        edges = []
        to_perms = {t: t.reading_word_permutation() for t in T}
        to_tab = {to_perms[k]: k for k in to_perms}
        Perm = permutation.Permutations()
        for t in T:
            pt = list(to_perms[t])
            for i in range(2, n):
                ii = pt.index(i)
                iip = pt.index(i+1)
                iim = pt.index(i-1)
                l = sorted([iim, ii, iip])
                if l[0] != ii:
                    continue
                x = pt[:]
                x[l[0]], x[l[2]] = x[l[2]], x[l[0]]
                if ii < iip:
                    e = [t, to_tab[Perm(x)], i]
                    edges.append(e)
                else:
                    e = [to_tab[Perm(x)], t, i]
                    edges.append(e)

        if directed:
            from sage.graphs.digraph import DiGraph
            self._DDEG = DiGraph([T, edges], format='vertices_and_edges',
                                 immutable=True, multiedges=True)
        else:
            from sage.graphs.graph import Graph
            self._DEG = Graph([T, edges], format='vertices_and_edges',
                              immutable=True, multiedges=True)
        return self.dual_equivalence_graph(directed, coloring)

    def specht_module(self, base_ring=None):
        r"""
        Return the Specht module corresponding to ``self``.

        EXAMPLES::

            sage: SM = Partition([2,2,1]).specht_module(QQ); SM                         # needs sage.modules
            Specht module of [2, 2, 1] over Rational Field
            sage: SM.frobenius_image()                                                  # needs sage.modules
            s[2, 2, 1]
        """
        from sage.combinat.specht_module import SpechtModule
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ
        R = SymmetricGroupAlgebra(base_ring, sum(self))
        return SpechtModule(R, self)

    def garsia_procesi_module(self, base_ring=None):
        r"""
        Return the :class:`Garsia-Procesi module
        <sage.combinat.symmetric_group_representations.GarsiaProcesiModule>`
        corresponding to ``self``.

        INPUT:

        - ``base_ring`` -- (default: `\QQ`) the base ring

        EXAMPLES::

            sage: GP = Partition([3,2,1]).garsia_procesi_module(QQ); GP
            Garsia-Procesi module of shape [3, 2, 1] over Rational Field
            sage: GP.graded_frobenius_image()
            q^4*s[3, 2, 1] + q^3*s[3, 3] + q^3*s[4, 1, 1] + (q^3+q^2)*s[4, 2]
             + (q^2+q)*s[5, 1] + s[6]

            sage: Partition([3,2,1]).garsia_procesi_module(GF(3))
            Garsia-Procesi module of shape [3, 2, 1] over Finite Field of size 3
        """
        from sage.combinat.symmetric_group_representations import GarsiaProcesiModule
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ
        R = SymmetricGroupAlgebra(base_ring, sum(self))
        return GarsiaProcesiModule(R, self)

    def specht_module_dimension(self, base_ring=None):
        r"""
        Return the dimension of the Specht module corresponding to ``self``.

        This is equal to the number of standard tableaux of shape ``self`` when
        over a field of characteristic `0`.

        INPUT:

        - ``base_ring`` -- (default: `\QQ`) the base ring

        EXAMPLES::

            sage: Partition([2,2,1]).specht_module_dimension()
            5
            sage: Partition([2,2,1]).specht_module_dimension(GF(2))                     # needs sage.rings.finite_rings
            5
        """
        from sage.categories.fields import Fields
        if base_ring is None or (base_ring in Fields() and base_ring.characteristic() == 0):
            from sage.combinat.tableau import StandardTableaux
            return StandardTableaux(self).cardinality()
        from sage.combinat.specht_module import specht_module_rank
        return specht_module_rank(self, base_ring)

    def simple_module_dimension(self, base_ring=None):
        r"""
        Return the dimension of the simple module corresponding to ``self``.

        When the base ring is a field of characteristic `0`, this is equal
        to the dimension of the Specht module.

        INPUT:

        - ``base_ring`` -- (default: `\QQ`) the base ring

        EXAMPLES::

            sage: Partition([2,2,1]).simple_module_dimension()
            5
            sage: Partition([2,2,1]).specht_module_dimension(GF(3))                     # needs sage.rings.finite_rings
            5
            sage: Partition([2,2,1]).simple_module_dimension(GF(3))                     # needs sage.rings.finite_rings
            4

            sage: for la in Partitions(6, regular=3):
            ....:     print(la, la.specht_module_dimension(), la.simple_module_dimension(GF(3)))
            [6] 1 1
            [5, 1] 5 4
            [4, 2] 9 9
            [4, 1, 1] 10 6
            [3, 3] 5 1
            [3, 2, 1] 16 4
            [2, 2, 1, 1] 9 9
        """
        from sage.categories.fields import Fields
        if base_ring is None or (base_ring in Fields() and base_ring.characteristic() == 0):
            from sage.combinat.tableau import StandardTableaux
            return StandardTableaux(self).cardinality()
        from sage.combinat.specht_module import simple_module_rank
        return simple_module_rank(self, base_ring)

    def tabloid_module(self, base_ring=None):
        r"""
        Return the tabloid module corresponding to ``self``.

        EXAMPLES::

            sage: TM = Partition([2,2,1]).tabloid_module(QQ); TM
            Tabloid module of [2, 2, 1] over Rational Field
            sage: TM.frobenius_image()
            s[2, 2, 1] + s[3, 1, 1] + 2*s[3, 2] + 2*s[4, 1] + s[5]
        """
        from sage.combinat.specht_module import TabloidModule
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ
        R = SymmetricGroupAlgebra(base_ring, sum(self))
        return TabloidModule(R, self)


##############
# Partitions #
##############

class Partitions(UniqueRepresentation, Parent):
    r"""
    ``Partitions(n, **kwargs)`` returns the combinatorial class of
    integer partitions of `n` subject to the constraints given by the
    keywords.

    Valid keywords are: ``starting``, ``ending``, ``min_part``,
    ``max_part``, ``max_length``, ``min_length``, ``length``,
    ``max_slope``, ``min_slope``, ``inner``, ``outer``, ``parts_in``,
    ``regular``, and ``restricted``. They have the following meanings:

    - ``starting=p`` specifies that the partitions should all be less
      than or equal to `p` in lex order. This argument cannot be combined
      with any other (see :issue:`15467`).

    - ``ending=p`` specifies that the partitions should all be greater than
      or equal to `p` in lex order. This argument cannot be combined with any
      other (see :issue:`15467`).

    - ``length=k`` specifies that the partitions have
      exactly `k` parts.

    - ``min_length=k`` specifies that the partitions have
      at least `k` parts.

    - ``min_part=k`` specifies that all parts of the
      partitions are at least `k`.

    - ``inner=p`` specifies that the partitions must contain the
      partition `p`.

    - ``outer=p`` specifies that the partitions
      be contained inside the partition `p`.

    - ``min_slope=k`` specifies that the partitions have slope at least
      `k`; the slope at position `i` is the difference between the
      `(i+1)`-th part and the `i`-th part.

    - ``parts_in=S`` specifies that the partitions have parts in the set
      `S`, which can be any sequence of pairwise distinct positive
      integers. This argument cannot be combined with any other
      (see :issue:`15467`).

    - ``regular=ell`` specifies that the partitions are `\ell`-regular,
      and can only be combined with the ``max_length`` or ``max_part``, but
      not both, keywords if `n` is not specified

    - ``restricted=ell`` specifies that the partitions are `\ell`-restricted,
      and cannot be combined with any other keywords

    The ``max_*`` versions, along with ``inner`` and ``ending``, work
    analogously.

    Right now, the ``parts_in``, ``starting``, ``ending``, ``regular``, and
    ``restricted`` keyword arguments are mutually exclusive, both of each
    other and of other keyword arguments. If you specify, say, ``parts_in``,
    all other keyword arguments will be ignored; ``starting``, ``ending``,
    ``regular``, and ``restricted`` work the same way.

    EXAMPLES:

    If no arguments are passed, then the combinatorial class
    of all integer partitions is returned::

        sage: Partitions()
        Partitions
        sage: [2,1] in Partitions()
        True

    If an integer `n` is passed, then the combinatorial class of integer
    partitions of `n` is returned::

        sage: Partitions(3)
        Partitions of the integer 3
        sage: Partitions(3).list()
        [[3], [2, 1], [1, 1, 1]]

    If ``starting=p`` is passed, then the combinatorial class of partitions
    greater than or equal to `p` in lexicographic order is returned::

        sage: Partitions(3, starting=[2,1])
        Partitions of the integer 3 starting with [2, 1]
        sage: Partitions(3, starting=[2,1]).list()
        [[2, 1], [1, 1, 1]]

    If ``ending=p`` is passed, then the combinatorial class of
    partitions at most `p` in lexicographic order is returned::

        sage: Partitions(3, ending=[2,1])
        Partitions of the integer 3 ending with [2, 1]
        sage: Partitions(3, ending=[2,1]).list()
        [[3], [2, 1]]

    Using ``max_slope=-1`` yields partitions into distinct parts -- each
    part differs from the next by at least 1. Use a different
    ``max_slope`` to get parts that differ by, say, 2::

        sage: Partitions(7, max_slope=-1).list()
        [[7], [6, 1], [5, 2], [4, 3], [4, 2, 1]]
        sage: Partitions(15, max_slope=-1).cardinality()
        27

    The number of partitions of `n` into odd parts equals the number of
    partitions into distinct parts. Let's test that for `n` from 10 to 20::

        sage: def test(n):
        ....:     return (Partitions(n, max_slope=-1).cardinality()
        ....:              == Partitions(n, parts_in=[1,3..n]).cardinality())
        sage: all(test(n) for n in [10..20])                                            # needs sage.libs.gap
        True

    The number of partitions of `n` into distinct parts that differ by
    at least 2 equals the number of partitions into parts that equal 1
    or 4 modulo 5; this is one of the Rogers-Ramanujan identities::

        sage: def test(n):
        ....:     return (Partitions(n, max_slope=-2).cardinality()
        ....:              == Partitions(n, parts_in=([1,6..n] + [4,9..n])).cardinality())
        sage: all(test(n) for n in [10..20])                                            # needs sage.libs.gap
        True

    Here are some more examples illustrating ``min_part``, ``max_part``,
    and ``length``::

        sage: Partitions(5, min_part=2)
        Partitions of 5 whose parts are at least 2
        sage: Partitions(5, min_part=2).list()
        [[5], [3, 2]]

    ::

        sage: Partitions(3, max_length=2).list()
        [[3], [2, 1]]

    ::

        sage: Partitions(10, min_part=2, length=3).list()
        [[6, 2, 2], [5, 3, 2], [4, 4, 2], [4, 3, 3]]

    Some examples using the ``regular`` keyword::

        sage: Partitions(regular=4)
        4-Regular Partitions
        sage: Partitions(regular=4, max_length=3)
        4-Regular Partitions with max length 3
        sage: Partitions(regular=4, max_part=3)
        4-Regular 3-Bounded Partitions
        sage: Partitions(3, regular=4)
        4-Regular Partitions of the integer 3

    Some examples using the ``restricted`` keyword::

        sage: Partitions(restricted=4)
        4-Restricted Partitions
        sage: Partitions(3, restricted=4)
        4-Restricted Partitions of the integer 3

    Here are some further examples using various constraints::

        sage: [x for x in Partitions(4)]
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        sage: [x for x in Partitions(4, length=2)]
        [[3, 1], [2, 2]]
        sage: [x for x in Partitions(4, min_length=2)]
        [[3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        sage: [x for x in Partitions(4, max_length=2)]
        [[4], [3, 1], [2, 2]]
        sage: [x for x in Partitions(4, min_length=2, max_length=2)]
        [[3, 1], [2, 2]]
        sage: [x for x in Partitions(4, max_part=2)]
        [[2, 2], [2, 1, 1], [1, 1, 1, 1]]
        sage: [x for x in Partitions(4, min_part=2)]
        [[4], [2, 2]]
        sage: [x for x in Partitions(4, outer=[3,1,1])]
        [[3, 1], [2, 1, 1]]
        sage: [x for x in Partitions(4, outer=[infinity, 1, 1])]
        [[4], [3, 1], [2, 1, 1]]
        sage: [x for x in Partitions(4, inner=[1,1,1])]
        [[2, 1, 1], [1, 1, 1, 1]]
        sage: [x for x in Partitions(4, max_slope=-1)]
        [[4], [3, 1]]
        sage: [x for x in Partitions(4, min_slope=-1)]
        [[4], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        sage: [x for x in Partitions(11, max_slope=-1, min_slope=-3, min_length=2, max_length=4)]
        [[7, 4], [6, 5], [6, 4, 1], [6, 3, 2], [5, 4, 2], [5, 3, 2, 1]]
        sage: [x for x in Partitions(11, max_slope=-1, min_slope=-3, min_length=2, max_length=4, outer=[6,5,2])]
        [[6, 5], [6, 4, 1], [6, 3, 2], [5, 4, 2]]

    Note that if you specify ``min_part=0``, then it will treat the minimum
    part as being 1 (see :issue:`13605`)::

        sage: [x for x in Partitions(4, length=3, min_part=0)]
        [[2, 1, 1]]
        sage: [x for x in Partitions(4, min_length=3, min_part=0)]
        [[2, 1, 1], [1, 1, 1, 1]]

    Except for very special cases, counting is done by brute force iteration
    through all the partitions. However the iteration itself has a reasonable
    complexity (see :class:`IntegerListsLex`), which allows for
    manipulating large partitions::

        sage: Partitions(1000, max_length=1).list()
        [[1000]]

    In particular, getting the first element is also constant time::

        sage: Partitions(30, max_part=29).first()
        [29, 1]

    TESTS::

        sage: TestSuite(Partitions(0)).run()                                            # needs sage.libs.flint
        sage: TestSuite(Partitions(5)).run()                                            # needs sage.libs.flint
        sage: TestSuite(Partitions(5, min_part=2)).run()                                # needs sage.libs.flint

        sage: repr(Partitions(5, min_part=2))
        'Partitions of 5 whose parts are at least 2'

        sage: P = Partitions(5, min_part=2)
        sage: P.first().parent()
        Partitions...
        sage: [2,1] in P
        False
        sage: [2,2,1] in P
        False
        sage: [3,2] in P
        True

        sage: Partitions(5, inner=[2,1], min_length=3).list()
        [[3, 1, 1], [2, 2, 1], [2, 1, 1, 1]]
        sage: Partitions(5, inner=Partition([2,2]), min_length=3).list()
        [[2, 2, 1]]
        sage: Partitions(7, inner=(2, 2), min_length=3).list()
        [[4, 2, 1], [3, 3, 1], [3, 2, 2], [3, 2, 1, 1], [2, 2, 2, 1], [2, 2, 1, 1, 1]]
        sage: Partitions(5, inner=[2,0,0,0,0,0]).list()
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1]]
        sage: Partitions(6, length=2, max_slope=-1).list()
        [[5, 1], [4, 2]]

        sage: Partitions(length=2, max_slope=-1).list()
        Traceback (most recent call last):
        ...
        ValueError: the size must be specified with any keyword argument

        sage: Partitions(max_part=3)
        3-Bounded Partitions

    Check that :issue:`14145` has been fixed::

        sage: 1 in Partitions()
        False

    Check :issue:`15467`::

        sage: Partitions(5, parts_in=[1,2,3,4], length=4)
        Traceback (most recent call last):
        ...
        ValueError: the parameters 'parts_in', 'starting', 'ending', 'regular' and 'restricted' cannot be combined with anything else
        sage: Partitions(5, starting=[3,2], length=2)
        Traceback (most recent call last):
        ...
        ValueError: the parameters 'parts_in', 'starting', 'ending', 'regular' and 'restricted' cannot be combined with anything else
        sage: Partitions(5, ending=[3,2], length=2)
        Traceback (most recent call last):
        ...
        ValueError: the parameters 'parts_in', 'starting', 'ending', 'regular' and 'restricted' cannot be combined with anything else
        sage: Partitions(5, restricted=2, length=2)
        Traceback (most recent call last):
        ...
        ValueError: the parameters 'parts_in', 'starting', 'ending', 'regular' and 'restricted' cannot be combined with anything else
        sage: Partitions(5, regular=5, length=2)
        Traceback (most recent call last):
        ...
        ValueError: the parameters 'parts_in', 'starting', 'ending', 'regular' and 'restricted' cannot be combined with anything else
        sage: Partitions(NN, length=2)
        Traceback (most recent call last):
        ...
        ValueError: the size must be specified with any keyword argument
        sage: Partitions(('la','la','laaaa'), max_part=8)
        Traceback (most recent call last):
        ...
        ValueError: n must be an integer or be equal to one of None, NN, NonNegativeIntegers()

    Check that calling ``Partitions`` with ``outer=a`` no longer
    mutates ``a`` (:issue:`16234`)::

        sage: a = [4,3,2,1,1,1,1]
        sage: for p in Partitions(8, outer=a, min_slope=-1):
        ....:    print(p)
        [3, 3, 2]
        [3, 2, 2, 1]
        [3, 2, 1, 1, 1]
        [2, 2, 2, 1, 1]
        [2, 2, 1, 1, 1, 1]
        [2, 1, 1, 1, 1, 1, 1]
        sage: a
        [4, 3, 2, 1, 1, 1, 1]

    Check that ``inner`` and ``outer`` indeed accept a partition as
    argument (:issue:`18423`)::

        sage: P = Partitions(5, inner=Partition([2,1]), outer=Partition([3,2])); P
        Partitions of the integer 5 satisfying constraints inner=[2, 1], outer=[3, 2]
        sage: P.list()
        [[3, 2]]

    Check that contradictory length requirements are handled correctly::

        sage: list(Partitions(5, max_length=1, length=3))
        Traceback (most recent call last):
        ...
        ValueError: do not specify the length together with the minimal or maximal length

        sage: list(Partitions(5, min_length=2, max_length=1))
        []

    Check that :issue:`38897` is fixed::

        sage: Partitions(40, min_length=10).cardinality()
        24000

        sage: Partitions(40, max_length=10).cardinality()
        16928
    """
    @staticmethod
    def __classcall_private__(cls, n=None, **kwargs):
        """
        Return the correct parent based upon the input.

        TESTS::

            sage: P = Partitions()
            sage: P2 = Partitions(NN)
            sage: P is P2
            True
            sage: P2 = Partitions(NonNegativeIntegers())
            sage: P is P2
            True
            sage: P = Partitions(4)
            sage: P2 = Partitions(int(4))
            sage: P is P2
            True

        Check that :issue:`17898` is fixed::

            sage: P = Partitions(5, min_slope=0)
            sage: list(P)
            [[5], [1, 1, 1, 1, 1]]
        """
        if n == infinity:
            raise ValueError("n cannot be infinite")
        if isinstance(n, (int, Integer)):
            if not kwargs:
                return Partitions_n(n)

            if len(kwargs) == 1:
                if 'max_part' in kwargs:
                    return Partitions_parts_length_restricted(n, ZZ.one(), kwargs['max_part'],
                                                              ZZ.zero(), n)
                if 'min_part' in kwargs:
                    return Partitions_parts_length_restricted(n, kwargs['min_part'], n,
                                                              ZZ.zero(), n)
                if 'length' in kwargs:
                    return Partitions_nk(n, kwargs['length'])
                if 'parts_in' in kwargs:
                    return Partitions_parts_in(n, kwargs['parts_in'])
                if 'starting' in kwargs:
                    return Partitions_starting(n, kwargs['starting'])
                if 'ending' in kwargs:
                    return Partitions_ending(n, kwargs['ending'])
                if 'regular' in kwargs:
                    return RegularPartitions_n(n, kwargs['regular'])
                if 'restricted' in kwargs:
                    return RestrictedPartitions_n(n, kwargs['restricted'])

            else:
                if ('parts_in' in kwargs or
                    'starting' in kwargs or
                    'ending' in kwargs or
                    'regular' in kwargs or
                    'restricted' in kwargs):
                    raise ValueError("the parameters 'parts_in', 'starting', "
                                     + "'ending', 'regular' and 'restricted' "
                                     + "cannot be combined with anything else")

                if 'length' in kwargs and ('min_length' in kwargs or 'max_length' in kwargs):
                    raise ValueError("do not specify the length together with the minimal or maximal length")

            if set(kwargs).issubset(['length', 'min_part', 'max_part',
                                     'min_length', 'max_length']):
                min_part = max(kwargs.get('min_part', ZZ.one()), ZZ.one())
                max_part = max(min(kwargs.get('max_part', n), n), ZZ.zero())
                if 'length' in kwargs:
                    k = ZZ(kwargs['length'])
                    return Partitions_parts_length_restricted(n, min_part, max_part, k, k)

                min_length = max(kwargs.get('min_length', ZZ.zero()), ZZ.zero())
                max_length = min(kwargs.get('max_length', n), n)
                return Partitions_parts_length_restricted(n, min_part, max_part, min_length, max_length)

            # FIXME: should inherit from IntegerListLex, and implement repr, or _name as a lazy attribute
            kwargs['name'] = "Partitions of the integer {} satisfying constraints {}".format(n, ", ".join(["{}={}".format(key, kwargs[key]) for key in sorted(kwargs)]))

            # min_part is at least 1, and it is 1 by default
            kwargs['min_part'] = max(1, kwargs.get('min_part', 1))

            # max_slope is at most 0, and it is 0 by default
            kwargs['max_slope'] = min(0, kwargs.get('max_slope', 0))

            if kwargs.get('min_slope', -float('inf')) > 0:
                raise ValueError("the minimum slope must be nonnegative")

            if 'outer' in kwargs:
                kwargs['max_length'] = min(len(kwargs['outer']),
                                           kwargs.get('max_length', infinity))

                kwargs['ceiling'] = tuple(kwargs['outer'])
                del kwargs['outer']

            if 'inner' in kwargs:
                inner = [x for x in kwargs['inner'] if x > 0]
                kwargs['floor'] = inner
                kwargs['min_length'] = max(len(inner),
                                           kwargs.get('min_length', 0))
                del kwargs['inner']
            return Partitions_with_constraints(n, **kwargs)

        if n is None or n is NN or n is NonNegativeIntegers():
            if not kwargs:
                return Partitions_all()

            if len(kwargs) == 1:
                if 'max_part' in kwargs:
                    return Partitions_all_bounded(kwargs['max_part'])
                if 'regular' in kwargs:
                    return RegularPartitions_all(kwargs['regular'])
                if 'restricted' in kwargs:
                    return RestrictedPartitions_all(kwargs['restricted'])
            elif len(kwargs) == 2:
                if 'regular' in kwargs:
                    if kwargs['regular'] < 1 or kwargs['regular'] not in ZZ:
                        raise ValueError("the regularity must be a positive integer")
                    if 'max_part' in kwargs:
                        return RegularPartitions_bounded(kwargs['regular'], kwargs['max_part'])
                    if 'max_length' in kwargs:
                        return RegularPartitions_truncated(kwargs['regular'], kwargs['max_length'])
            raise ValueError("the size must be specified with any keyword argument")

        raise ValueError("n must be an integer or be equal to one of "
                         "None, NN, NonNegativeIntegers()")

    def __init__(self, is_infinite=False):
        """
        Initialize ``self``.

        INPUT:

        - ``is_infinite`` -- boolean (default: ``False``); if ``True``, then
          the number of partitions in this set is infinite

        EXAMPLES::

            sage: Partitions()
            Partitions
            sage: Partitions(2)
            Partitions of the integer 2
        """
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

    Element = Partition

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the global options for elements of the partition,
        skew partition, and partition tuple classes.  If no parameters are
        set, then the function returns a copy of the options dictionary.

        The ``options`` to partitions can be accessed as the method
        :obj:`Partitions.options` of :class:`Partitions` and
        related parent classes.

        @OPTIONS@

        EXAMPLES::

            sage: P = Partition([4,2,2,1])
            sage: P
            [4, 2, 2, 1]
            sage: Partitions.options.display="exp"
            sage: P
            1, 2^2, 4
            sage: Partitions.options.display="exp_high"
            sage: P
            4, 2^2, 1

        It is also possible to use user defined functions for the ``display`` and
        ``latex`` options::

            sage: Partitions.options(display=lambda mu: '<%s>' % ','.join('%s'%m for m in mu._list)); P
            <4,2,2,1>
            sage: Partitions.options(latex=lambda mu: '\\Diagram{%s}' % ','.join('%s'%m for m in mu._list)); latex(P)
            \Diagram{4,2,2,1}
            sage: Partitions.options(display='diagram', diagram_str='#')
            sage: P
            ####
            ##
            ##
            #
            sage: Partitions.options(diagram_str='*', convention='french')
            sage: print(P.ferrers_diagram())
            *
            **
            **
            ****

        Changing the ``convention`` for partitions also changes the ``convention``
        option for tableaux and vice versa::

            sage: T = Tableau([[1,2,3],[4,5]])
            sage: T.pp()
              4  5
              1  2  3
            sage: Tableaux.options.convention="english"
            sage: print(P.ferrers_diagram())
            ****
            **
            **
            *
            sage: T.pp()
              1  2  3
              4  5
            sage: Partitions.options._reset()
        """
        NAME = 'Partitions'
        module = 'sage.combinat.partition'
        display = {'default': "list",
                   'description': 'Specifies how partitions should be printed',
                   'values': {'list': 'displayed as a list',
                              'exp_low': 'in exponential form (lowest first)',
                              'exp_high': 'in exponential form (highest first)',
                              'diagram': 'as a Ferrers diagram',
                              'compact_low': 'compact form of ``exp_low``',
                              'compact_high': 'compact form of ``exp_high``'},
                   'alias': {'exp': "exp_low", 'compact': "compact_low", 'array': "diagram",
                             'ferrers_diagram': "diagram", 'young_diagram': "diagram"},
                   'case_sensitive': False}
        latex = {'default': "young_diagram",
                 'description': 'Specifies how partitions should be latexed',
                 'values': {'diagram': 'latex as a Ferrers diagram',
                            'young_diagram': 'latex as a Young diagram',
                            'list': 'latex as a list',
                            'exp_high': 'latex as a list in exponential notation (highest first)',
                            'exp_low': 'as a list latex in exponential notation (lowest first)'},
                 'alias': {'exp': "exp_low", 'array': "diagram", 'ferrers_diagram': "diagram"},
                 'case_sensitive': False}
        diagram_str = {'default': "*",
                       'description': 'The character used for the cells when printing Ferrers diagrams',
                       'checker': lambda char: isinstance(char, str)}
        latex_diagram_str = {'default': "\\ast",
                             'description': 'The character used for the cells when latexing Ferrers diagrams',
                             'checker': lambda char: isinstance(char, str)}
        convention = {'link_to': (tableau.Tableaux.options, 'convention')}
        notation = {'alt_name': 'convention'}

    def __reversed__(self):
        """
        A reversed iterator.

        EXAMPLES::

            sage: [x for x in reversed(Partitions(4))]
            [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
        """
        if not self.is_finite():
            raise NotImplementedError("the set is infinite, so this needs a custom reverse iterator")

        for i in reversed(range(self.cardinality())):
            yield self[i]

    def _element_constructor_(self, lst):
        """
        Construct an element with ``self`` as parent.

        EXAMPLES::

            sage: P = Partitions()
            sage: p = P([3,3,1]); p
            [3, 3, 1]
            sage: P(p) is p
            True
            sage: P([3, 2, 1, 0])
            [3, 2, 1]

            sage: PT = PartitionTuples()
            sage: elt = PT([[4,4,2,2,1]]); elt
            ([4, 4, 2, 2, 1])
            sage: P(elt)
            [4, 4, 2, 2, 1]

        TESTS::

            sage: Partition([3/2])
            Traceback (most recent call last):
            ...
            ValueError: all parts of [3/2] should be nonnegative integers
        """
        if isinstance(lst, PartitionTuple):
            if lst.level() != 1:
                raise ValueError(f'{lst} is not an element of {self}')
            lst = lst[0]
            if lst.parent() is self:
                return lst
        try:
            lst = list(map(ZZ, lst))
        except TypeError:
            raise ValueError('all parts of %s should be nonnegative integers' % repr(lst))

        if lst in self:
            # trailing zeros are removed in Partition.__init__
            return self.element_class(self, lst)

        raise ValueError(f'{lst} is not an element of {self}')

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        TESTS::

            sage: P = Partitions()
            sage: Partition([2,1]) in P
            True
            sage: [2,1] in P
            True
            sage: [3,2,1] in P
            True
            sage: [1,2] in P
            False
            sage: [] in P
            True
            sage: [0] in P
            True

        Check that types that represent integers are not excluded::

            sage: P = Partitions()
            sage: [3/1, 2/2] in P
            True
            sage: Partition([3/1, 2]) in P
            True

        Check that non-integers and non-lists are excluded::

            sage: P = Partitions()
            sage: [2,1.5] in P
            False

            sage: 0 in P
            False
        """
        if isinstance(x, Partition):
            return True
        if isinstance(x, (list, tuple)):
            return not x or (all((a in ZZ) and (a >= b) for a, b in zip(x, x[1:]))
                             and (x[-1] in ZZ) and (x[-1] >= 0))
        return False

    def subset(self, *args, **kwargs):
        r"""
        Return ``self`` if no arguments are given.

        Otherwise, it raises a :exc:`ValueError`.

        EXAMPLES::

            sage: P = Partitions(5, starting=[3,1]); P
            Partitions of the integer 5 starting with [3, 1]
            sage: P.subset()
            Partitions of the integer 5 starting with [3, 1]
            sage: P.subset(ending=[3,1])
            Traceback (most recent call last):
            ...
            ValueError: invalid combination of arguments
        """
        if len(args) != 0 or len(kwargs) != 0:
            raise ValueError("invalid combination of arguments")
        return self


class Partitions_all(Partitions):
    """
    Class of all partitions.

    TESTS::

        sage: TestSuite( sage.combinat.partition.Partitions_all() ).run()
    """

    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: P = Partitions()
            sage: P.category()
            Category of infinite enumerated sets
            sage: Partitions().cardinality()
            +Infinity
            sage: TestSuite(P).run()
        """
        Partitions.__init__(self, is_infinite=True)

    def subset(self, size=None, **kwargs):
        """
        Return the subset of partitions of a given size and additional
        keyword arguments.

        EXAMPLES::

            sage: P = Partitions()
            sage: P.subset(4)
            Partitions of the integer 4
        """
        if size is None:
            return self
        return Partitions(size, **kwargs)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: Partitions() # indirect doctest
            Partitions
        """
        return "Partitions"

    def __iter__(self):
        """
        An iterator for all partitions.

        EXAMPLES::

            sage: p = Partitions()
            sage: it = p.__iter__()
            sage: [next(it) for i in range(10)]
            [[], [1], [2], [1, 1], [3], [2, 1], [1, 1, 1], [4], [3, 1], [2, 2]]
        """
        n = 0
        while True:
            for p in ZS1_iterator(n):
                yield self.element_class(self, p)
            n += 1

    def __reversed__(self):
        """
        A reversed iterator for all partitions.

        This reverse iterates through partitions of fixed `n` and incrementing
        `n` after reaching the end.

        EXAMPLES::

            sage: p = Partitions()
            sage: revit = p.__reversed__()
            sage: [next(revit) for i in range(10)]
            [[], [1], [1, 1], [2], [1, 1, 1], [2, 1], [3], [1, 1, 1, 1], [2, 1, 1], [2, 2]]
        """
        n = 0
        while True:
            for p in reversed(list(ZS1_iterator(n))):
                yield self.element_class(self, p)
            n += 1

    def from_frobenius_coordinates(self, frobenius_coordinates):
        """
        Return a partition from a pair of sequences of Frobenius coordinates.

        EXAMPLES::

            sage: Partitions().from_frobenius_coordinates(([],[]))
            []
            sage: Partitions().from_frobenius_coordinates(([0],[0]))
            [1]
            sage: Partitions().from_frobenius_coordinates(([1],[1]))
            [2, 1]
            sage: Partitions().from_frobenius_coordinates(([6,3,2],[4,1,0]))
            [7, 5, 5, 1, 1]
        """
        if len(frobenius_coordinates) != 2:
            raise ValueError('%s is not a valid partition, two sequences of coordinates are needed' % str(frobenius_coordinates))
        else:
            a = frobenius_coordinates[0]
            b = frobenius_coordinates[1]
            if len(a) != len(b):
                raise ValueError('%s is not a valid partition, the sequences of coordinates need to be the same length' % str(frobenius_coordinates))
                # should add tests to see if a and b are sorted down, nonnegative and strictly decreasing
        r = len(a)
        if r == 0:
            return self.element_class(self, [])
        tmp = [a[i]+i+1 for i in range(r)]
        # should check that a is strictly decreasing
        if a[-1] < 0:
            raise ValueError('%s is not a partition, no coordinate can be negative' % str(frobenius_coordinates))
        if b[-1] >= 0:
            tmp.extend([r]*b[r-1])
        else:
            raise ValueError('%s is not a partition, no coordinate can be negative' % str(frobenius_coordinates))
        for i in range(r - 1, 0, -1):
            if b[i-1]-b[i] > 0:
                tmp.extend([i]*(b[i-1]-b[i]-1))
            else:
                raise ValueError('%s is not a partition, the coordinates need to be strictly decreasing' % str(frobenius_coordinates))
        return self.element_class(self, tmp)

    def from_beta_numbers(self, beta):
        r"""
        Return a partition corresponding to a sequence of beta numbers.

        A sequence of beta numbers is a strictly increasing sequence
        `0 \leq b_1 < \cdots < b_k` of nonnegative integers. The
        corresponding partition `\mu = (\mu_k, \ldots, \mu_1)` is
        given by `\mu_i = [1,i) \setminus \{ b_1, \ldots, b_i \}`. This gives
        a bijection from the set of partitions with at most `k` nonzero parts
        to the set of strictly increasing sequences of nonnegative integers
        of length `k`.

        EXAMPLES::

            sage: Partitions().from_beta_numbers([0,1,2,4,5,8])
            [3, 1, 1]
            sage: Partitions().from_beta_numbers([0,2,3,6])
            [3, 1, 1]
        """
        beta.sort()  # put them into increasing order just in case
        offset = 0
        while offset < len(beta)-1 and beta[offset] == offset:
            offset += 1
        beta = beta[offset:]
        mu = [beta[i]-offset-i for i in range(len(beta))]
        return self.element_class(self, list(reversed(mu)))

    def from_exp(self, exp):
        """
        Return a partition from its list of multiplicities.

        EXAMPLES::

            sage: Partitions().from_exp([2,2,1])
            [3, 2, 2, 1, 1]
        """
        p = []
        for i in reversed(range(len(exp))):
            p += [i+1]*exp[i]
        return self.element_class(self, p)

    def from_zero_one(self, seq):
        r"""
        Return a partition from its `0-1` sequence.

        The full `0-1` sequence is the sequence (infinite in both
        directions) indicating the steps taken when following the
        outer rim of the diagram of the partition. We use the convention
        that in English convention, a 1 corresponds to an East step, and
        a 0 corresponds to a North step.

        Note that every full `0-1` sequence starts with infinitely many 0s and
        ends with infinitely many 1s.

        .. SEEALSO::

            :meth:`Partition.zero_one_sequence()`

        INPUT:

        The input should be a finite sequence of 0s and 1s. The
        heading 0s and trailing 1s will be discarded.

        EXAMPLES::

            sage: Partitions().from_zero_one([])
            []
            sage: Partitions().from_zero_one([1,0])
            [1]
            sage: Partitions().from_zero_one([1, 1, 1, 1, 0, 1, 0])
            [5, 4]

        Heading 0s and trailing 1s are correctly handled::

            sage: Partitions().from_zero_one([0,0,1,1,1,1,0,1,0,1,1,1])
            [5, 4]

        TESTS::

            sage: all(Partitions().from_zero_one(mu.zero_one_sequence()) == mu for n in range(10) for mu in Partitions(n))
            True
        """
        tmp = [i for i in range(len(seq)) if seq[i] == 0]
        return self.element_class(self, [tmp[i]-i for i in range(len(tmp)-1, -1, -1)])

    def from_core_and_quotient(self, core, quotient):
        """
        Return a partition from its core and quotient.

        Algorithm from mupad-combinat.

        EXAMPLES::

            sage: Partitions().from_core_and_quotient([2,1], [[2,1],[3],[1,1,1]])
            [11, 5, 5, 3, 2, 2, 2]

        TESTS::

            sage: Partitions().from_core_and_quotient([2,1], [[2,1],[2,3,1],[1,1,1]])
            Traceback (most recent call last):
            ...
            ValueError: the quotient [[2, 1], [2, 3, 1], [1, 1, 1]] must be a tuple of partitions

        We check that :issue:`11412` is actually fixed::

            sage: test = lambda x, k: x == Partition(core=x.core(k),
            ....:                                    quotient=x.quotient(k))
            sage: all(test(mu,k) for k in range(1,5)
            ....:     for n in range(10) for mu in Partitions(n))
            True
            sage: test2 = lambda core, mus: (
            ....:     Partition(core=core, quotient=mus).core(mus.level()) == core
            ....:     and
            ....:     Partition(core=core, quotient=mus).quotient(mus.level()) == mus)
            sage: all(test2(core,mus)  # long time (5s on sage.math, 2011)
            ....:     for k in range(1,10)
            ....:     for n_core in range(10-k)
            ....:     for core in Partitions(n_core)
            ....:     if core.core(k) == core
            ....:     for n_mus in range(10-k)
            ....:     for mus in PartitionTuples(k,n_mus))
            True
        """
        from .partition_tuple import PartitionTuple, PartitionTuples
        if quotient not in PartitionTuples():
            raise ValueError('the quotient %s must be a tuple of partitions' % quotient)
        components = PartitionTuple(quotient).components()
        length = len(components)
        k = length*max(len(q) for q in components) + len(core)
        # k needs to be large enough. this seems to me like the smallest it can be
        v = [core[i]-i for i in range(len(core))] + [-i for i in range(len(core), k)]
        w = [[x for x in v if (x-i) % length == 0] for i in range(1, length+1)]
        new_w = []
        for i in range(length):
            lw = len(w[i])
            lq = len(components[i])
            # k needs to be chosen so lw >= lq
            new_w += [w[i][j] + length*components[i][j] for j in range(lq)]
            new_w += [w[i][j] for j in range(lq, lw)]
        new_w.sort(reverse=True)
        return self.element_class(self, [new_w[i]+i for i in range(len(new_w))])


class Partitions_all_bounded(Partitions):
    def __init__(self, k):
        """
        TESTS::

            sage: TestSuite( sage.combinat.partition.Partitions_all_bounded(3) ).run() # long time
        """
        self.k = k
        Partitions.__init__(self, is_infinite=True)

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions(max_part=3)
            sage: Partition([2,1]) in P
            True
            sage: [2,1] in P
            True
            sage: [3,2,1] in P
            True
            sage: [1,2] in P
            False
            sage: [5,1] in P
            False
            sage: [0] in P
            True
            sage: [] in P
            True
        """
        return not x or (x[0] <= self.k and x in _Partitions)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.partition import Partitions_all_bounded
            sage: Partitions_all_bounded(3)
            3-Bounded Partitions
        """
        return "%d-Bounded Partitions" % self.k

    def __iter__(self):
        """
        An iterator for all `k`-bounded partitions.

        EXAMPLES::

            sage: p = Partitions(max_part=3)
            sage: it = p.__iter__()
            sage: [next(it) for i in range(10)]
            [[], [1], [2], [1, 1], [3], [2, 1], [1, 1, 1], [3, 1], [2, 2], [2, 1, 1]]
        """
        n = 0
        while True:
            for p in Partitions(n, max_part=self.k):
                yield self.element_class(self, p)
            n += 1


class Partitions_n(Partitions):
    """
    Partitions of the integer `n`.

    TESTS::

        sage: TestSuite( sage.combinat.partition.Partitions_n(0) ).run()                # needs sage.libs.flint
        sage: TestSuite( sage.combinat.partition.Partitions_n(0) ).run()                # needs sage.libs.flint
    """

    def __init__(self, n):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(  Partitions(5) ).run()
        """
        Partitions.__init__(self)
        self.n = n

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        TESTS::

            sage: p = Partitions(5)
            sage: [2,1] in p
            False
            sage: [2,2,1] in p
            True
            sage: [3,2] in p
            True
            sage: [2,3] in p
            False
        """
        return x in _Partitions and sum(x) == self.n

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        TESTS::

            sage: Partitions(5) # indirect doctest
            Partitions of the integer 5
        """
        return "Partitions of the integer %s" % self.n

    def _an_element_(self):
        """
        Return a partition in ``self``.

        EXAMPLES::

            sage: Partitions(4).an_element()  # indirect doctest
            [3, 1]
            sage: Partitions(0).an_element()
            []
            sage: Partitions(1).an_element()
            [1]
        """
        if self.n == 0:
            lst = []
        elif self.n == 1:
            lst = [1]
        else:
            lst = [self.n - 1, 1]
        return self.element_class(self, lst)

    def cardinality(self, algorithm='flint'):
        r"""
        Return the number of partitions of the specified size.

        INPUT:

        - ``algorithm`` -- (default: ``'flint'``)

          - ``'flint'`` -- use FLINT (currently the fastest)
          - ``'gap'`` -- use GAP (VERY *slow*)
          - ``'pari'`` -- use PARI. Speed seems the same as GAP until
            `n` is in the thousands, in which case PARI is faster

        It is possible to associate with every partition of the integer `n` a
        conjugacy class of permutations in the symmetric group on `n` points
        and vice versa. Therefore the number of partitions `p_n` is the number
        of conjugacy classes of the symmetric group on `n` points.

        EXAMPLES::

            sage: v = Partitions(5).list(); v
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: len(v)
            7
            sage: Partitions(5).cardinality(algorithm='gap')                            # needs sage.libs.gap
            7

        ::

            sage: # needs sage.libs.flint
            sage: Partitions(3).cardinality()
            3
            sage: number_of_partitions(5, algorithm='flint')
            7
            sage: Partitions(10).cardinality()
            42
            sage: Partitions(40).cardinality()
            37338
            sage: Partitions(100).cardinality()
            190569292

        ::

            sage: # needs sage.libs.pari
            sage: Partitions(3).cardinality(algorithm='pari')
            3
            sage: Partitions(5).cardinality(algorithm='pari')
            7
            sage: Partitions(10).cardinality(algorithm='pari')
            42

        A generating function for `p_n` is given by the reciprocal of
        Euler's function:

        .. MATH::

           \sum_{n=0}^{\infty} p_n x^n = \prod_{k=1}^{\infty} \frac{1}{1-x^k}.

        We use Sage to verify that the first several coefficients do
        indeed agree::

            sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
            sage: prod([(1-q^k)^(-1) for k in range(1,9)])  # partial product of
            1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
            sage: [Partitions(k).cardinality() for k in range(2,10)]                    # needs sage.libs.flint
            [2, 3, 5, 7, 11, 15, 22, 30]

        Another consistency test for ``n`` up to 500::

            sage: len([n for n in [1..500]                                              # needs sage.libs.flint sage.libs.pari
            ....:     if Partitions(n).cardinality() != Partitions(n).cardinality(algorithm='pari')])
            0

        For negative inputs, the result is zero (the algorithm is ignored)::

            sage: Partitions(-5).cardinality()
            0

        REFERENCES:

        - :wikipedia:`Partition\_(number\_theory)`
        """
        if self.n < 0:
            return ZZ.zero()

        if algorithm == 'flint':
            return cached_number_of_partitions(self.n)

        elif algorithm == 'gap':
            from sage.libs.gap.libgap import libgap
            return ZZ(libgap.NrPartitions(ZZ(self.n)))

        elif algorithm == 'pari':
            return ZZ(pari(ZZ(self.n)).numbpart())

        raise ValueError("unknown algorithm '%s'" % algorithm)

    def random_element(self, measure='uniform'):
        """
        Return a random partitions of `n` for the specified measure.

        INPUT:

        - ``measure`` -- ``'uniform'`` or ``'Plancherel'``
          (default: ``'uniform'``)

        .. SEEALSO::

            - :meth:`random_element_uniform`
            - :meth:`random_element_plancherel`

        EXAMPLES::

            sage: Partitions(5).random_element()  # random                              # needs sage.libs.flint
            [2, 1, 1, 1]
            sage: Partitions(5).random_element(measure='Plancherel')  # random          # needs sage.libs.flint
            [2, 1, 1, 1]
        """
        if measure == 'uniform':
            return self.random_element_uniform()
        elif measure == 'Plancherel':
            return self.random_element_plancherel()
        else:
            raise ValueError("Unknown measure: %s" % measure)

    def random_element_uniform(self):
        """
        Return a random partition of `n` with uniform probability.

        EXAMPLES::

            sage: Partitions(5).random_element_uniform()  # random                      # needs sage.libs.flint
            [2, 1, 1, 1]
            sage: Partitions(20).random_element_uniform()  # random                     # needs sage.libs.flint
            [9, 3, 3, 2, 2, 1]

        TESTS::

            sage: all(Part.random_element_uniform() in Part                             # needs sage.libs.flint
            ....:     for Part in map(Partitions, range(10)))
            True

        Check that :issue:`18752` is fixed::

            sage: P = Partitions(5)
            sage: la = P.random_element_uniform()                                       # needs sage.libs.flint
            sage: la.parent() is P                                                      # needs sage.libs.flint
            True

        ALGORITHM:

        - It is a python Implementation of RANDPAR, see [NW1978]_.  The
          complexity is unknown, there may be better algorithms.

           .. TODO::

               Check in Knuth AOCP4.

        - There is also certainly a lot of room for optimizations, see
          comments in the code.

        AUTHOR:

        - Florent Hivert (2009-11-23)
        """
        n = self.n
        res = []  # A dictionary of multiplicities could be faster.
        while n > 0:
            # Choose a pair d,j = 1,2..., with d*j <= n with probability
            #        d*numpart(n-d*j) / n / numpart(n)
            # and add d^j to the result partition. The resulting partitions is
            # equiprobable.

            # The following could be made faster by a clever use of floats
            rand = randrange(0, n*cached_number_of_partitions(n))  # cached number_of_partition

            # It is better to start by the j = 1 pairs because they are the
            # most probable. Maybe there is an even more clever order.
            for j in range(1, n+1):
                d = 1
                r = n-j        # n - d*j
                while r >= 0:
                    rand -= d * cached_number_of_partitions(r)
                    if rand < 0:
                        break
                    d += 1
                    r -= j
                else:
                    continue
                break
            res.extend([d]*j)
            n = r
        res.sort(reverse=True)
        return self.element_class(self, res)

    def random_element_plancherel(self):
        r"""
        Return a random partition of `n` (for the Plancherel measure).

        This probability distribution comes from the uniform distribution
        on permutations via the Robinson-Schensted correspondence.

        See :wikipedia:`Plancherel\_measure`
        and :meth:`Partition.plancherel_measure`.

        EXAMPLES::

            sage: Partitions(5).random_element_plancherel()   # random
            [2, 1, 1, 1]
            sage: Partitions(20).random_element_plancherel()  # random
            [9, 3, 3, 2, 2, 1]

        TESTS::

            sage: all(Part.random_element_plancherel() in Part
            ....:     for Part in map(Partitions, range(10)))
            True

        Check that :issue:`18752` is fixed::

            sage: P = Partitions(5)
            sage: la = P.random_element_plancherel()
            sage: la.parent() is P
            True

        ALGORITHM:

        - insert by Robinson-Schensted a uniform random permutations of n and
          returns the shape of the resulting tableau. The complexity is
          `O(n\ln(n))` which is likely optimal. However, the implementation
          could be optimized.

        AUTHOR:

        - Florent Hivert (2009-11-23)
        """
        T = permutation.Permutations(self.n).random_element().left_tableau()
        return self.element_class(self, [len(row) for row in T])

    def first(self):
        """
        Return the lexicographically first partition of a positive integer
        `n`. This is the partition ``[n]``.

        EXAMPLES::

            sage: Partitions(4).first()
            [4]
        """
        return self.element_class(self, [self.n])

    def next(self, p):
        """
        Return the lexicographically next partition after the partition ``p``.

        EXAMPLES::

            sage: Partitions(4).next([4])
            [3, 1]
            sage: Partitions(4).next([1,1,1,1]) is None
            True
        """
        if p is None:
            return None
        q = ZS1_next(list(p))
        if q:
            return self.element_class(self, q)
        return None

    def prev(self, p):
        r"""
        Return the lexicographically previous partition before partition ``p``.

        EXAMPLES::

            sage: Partitions(4).prev([3, 1])
            [4]
            sage: Partitions(4).prev([4]) is None
            True
        """
        if p is None:
            return None
        q = ZS2_next(p)
        if q:
            return self.element_class(self, q)
        return None

    def last(self):
        """
        Return the lexicographically last partition of the positive
        integer `n`. This is the all-ones partition.

        EXAMPLES::

            sage: Partitions(4).last()
            [1, 1, 1, 1]
        """
        return self.element_class(self, [1]*self.n)

    def __iter__(self):
        """
        An iterator for the partitions of `n`.

        EXAMPLES::

            sage: [x for x in Partitions(4)]
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

        TESTS::

            sage: all(isinstance(i, Integer) for p in Partitions(4) for i in p)
            True
        """
        for p in ZS1_iterator(self.n):
            yield self.element_class(self, [Integer(i) for i in p])

    def subset(self, **kwargs):
        r"""
        Return a subset of ``self`` with the additional optional arguments.

        EXAMPLES::

            sage: P = Partitions(5); P
            Partitions of the integer 5
            sage: P.subset(starting=[3,1])
            Partitions of the integer 5 starting with [3, 1]
        """
        return Partitions(self.n, **kwargs)


class Partitions_nk(Partitions):
    """
    Partitions of the integer `n` of length equal to `k`.

    TESTS::

        sage: TestSuite( sage.combinat.partition.Partitions_nk(0,0) ).run()
        sage: TestSuite( sage.combinat.partition.Partitions_nk(0,0) ).run()
    """

    def __init__(self, n, k):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(  Partitions(5, length=2) ).run()
        """
        Partitions.__init__(self)
        self.n = n
        self.k = k

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        TESTS::

            sage: p = Partitions(5, length=2)
            sage: [2,1] in p
            False
            sage: [2,2,1] in p
            False
            sage: [3,2] in p
            True
            sage: [2,3] in p
            False
            sage: [4,1] in p
            True
            sage: [1,1,1,1,1] in p
            False
            sage: [5] in p
            False
        """
        return x in _Partitions and sum(x) == self.n and len(x) == self.k

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: Partitions(5, length=2) # indirect doctest
            Partitions of the integer 5 of length 2
        """
        return f"Partitions of the integer {self.n} of length {self.k}"

    def _an_element_(self):
        """
        Return a partition in ``self``.

        EXAMPLES::

            sage: Partitions(4, length=1).an_element()  # indirect doctest
            [4]
            sage: Partitions(4, length=2).an_element()
            [3, 1]
            sage: Partitions(4, length=3).an_element()
            [2, 1, 1]
            sage: Partitions(4, length=4).an_element()
            [1, 1, 1, 1]

            sage: Partitions(1, length=1).an_element()
            [1]

            sage: Partitions(0, length=0).an_element()
            []
        """
        if self.n == 0:
            if self.k == 0:
                lst = []
            else:
                from sage.categories.sets_cat import EmptySetError
                raise EmptySetError
        elif self.n >= self.k > 0:
            lst = [self.n - self.k + 1] + [1] * (self.k-1)
        else:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError
        return self.element_class(self, lst)

    def __iter__(self):
        """
        An iterator for all partitions of `n` of length `k`.

        EXAMPLES::

            sage: p = Partitions(9, length=3)
            sage: it = p.__iter__()
            sage: list(it)
            [[7, 1, 1], [6, 2, 1], [5, 3, 1], [5, 2, 2], [4, 4, 1], [4, 3, 2], [3, 3, 3]]

            sage: p = Partitions(9, length=10)
            sage: list(p.__iter__())
            []

            sage: p = Partitions(0, length=0)
            sage: list(p.__iter__())
            [[]]

            sage: from sage.combinat.partition import number_of_partitions_length
            sage: all( len(Partitions(n, length=k).list())                              # needs sage.libs.flint
            ....:      == number_of_partitions_length(n, k)
            ....:      for n in range(9) for k in range(n+2) )
            True

        TESTS::

            sage: partitions = Partitions(9, length=3)
            sage: all(isinstance(i, Integer) for p in partitions for i in p)
            True
        """
        for p in ZS1_iterator_nk(self.n - self.k, self.k):
            v = [Integer(i + 1) for i in p]
            adds = [Integer(1)] * (self.k - len(v))
            yield self.element_class(self, v + adds)

    def cardinality(self, algorithm='hybrid'):
        r"""
        Return the number of partitions of the specified size with the
        specified length.

        INPUT:

        - ``algorithm`` -- (default: ``'hybrid'``) the algorithm to compute
          the cardinality and can be one of the following:

          * ``'hybrid'`` -- use a hybrid algorithm which uses heuristics to
            reduce the complexity
          * ``'gap'`` -- use GAP

        EXAMPLES::

            sage: v = Partitions(5, length=2).list(); v
            [[4, 1], [3, 2]]
            sage: len(v)
            2
            sage: Partitions(5, length=2).cardinality()
            2

        More generally, the number of partitions of `n` of length `2`
        is `\left\lfloor \frac{n}{2} \right\rfloor`::

            sage: all( Partitions(n, length=2).cardinality()
            ....:      == n // 2 for n in range(10) )
            True

        The number of partitions of `n` of length `1` is `1` for `n`
        positive::

            sage: all( Partitions(n, length=1).cardinality() == 1
            ....:      for n in range(1, 10) )
            True

        Further examples::

            sage: # needs sage.libs.flint
            sage: Partitions(5, length=3).cardinality()
            2
            sage: Partitions(6, length=3).cardinality()
            3
            sage: Partitions(8, length=4).cardinality()
            5
            sage: Partitions(8, length=5).cardinality()
            3
            sage: Partitions(15, length=6).cardinality()
            26
            sage: Partitions(0, length=0).cardinality()
            1
            sage: Partitions(0, length=1).cardinality()
            0
            sage: Partitions(1, length=0).cardinality()
            0
            sage: Partitions(1, length=4).cardinality()
            0

        TESTS:

        We check the hybrid approach gives the same results as GAP::

            sage: N = [0, 1, 2, 3, 5, 10, 20, 500, 850]
            sage: K = [0, 1, 2, 3, 5, 10, 11, 20, 21, 250, 499, 500]
            sage: all(Partitions(n, length=k).cardinality()                             # needs sage.libs.flint
            ....:       == Partitions(n,length=k).cardinality('gap')
            ....:     for n in N for k in K)
            True
            sage: P = Partitions(4562, length=2800)
            sage: P.cardinality() == P.cardinality('gap')                               # needs sage.libs.flint
            True
        """
        return number_of_partitions_length(self.n, self.k, algorithm)

    def subset(self, **kwargs):
        r"""
        Return a subset of ``self`` with the additional optional arguments.

        EXAMPLES::

            sage: P = Partitions(5, length=2); P
            Partitions of the integer 5 of length 2
            sage: P.subset(max_part=3)
            Partitions of 5 having length 2 and whose parts are at most 3
        """
        return Partitions(self.n, length=self.k, **kwargs)


class Partitions_parts_in(Partitions):
    """
    Partitions of `n` with parts in a given set `S`.

    This is invoked indirectly when calling
    ``Partitions(n, parts_in=parts)``, where ``parts`` is a list of
    pairwise distinct integers.

    TESTS::

        sage: TestSuite( sage.combinat.partition.Partitions_parts_in(6, parts=[2,1]) ).run()        # needs sage.libs.gap
    """

    @staticmethod
    def __classcall_private__(cls, n, parts):
        """
        Normalize the input to ensure a unique representation.

        TESTS::

            sage: P = Partitions(4, parts_in=[2,1])
            sage: P2 = Partitions(4, parts_in=(1,2))
            sage: P is P2
            True
        """
        parts = tuple(sorted(parts))
        return super().__classcall__(cls, Integer(n), parts)

    def __init__(self, n, parts):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(Partitions(5, parts_in=[1,2,3])).run()                      # needs sage.libs.gap
        """
        Partitions.__init__(self)
        self.n = n
        self.parts = list(parts)

    def __contains__(self, x):
        """
        TESTS::

            sage: p = Partitions(5, parts_in=[1,2])
            sage: [2,1,1,1] in p
            True
            sage: [4,1] in p
            False
        """
        return (x in _Partitions and sum(x) == self.n and
                all(p in self.parts for p in x))

    def _repr_(self):
        """
        TESTS::

            sage: Partitions(5, parts_in=[1,2,3]) # indirect doctest
            Partitions of the integer 5 with parts in [1, 2, 3]
        """
        return f"Partitions of the integer {self.n} with parts in {self.parts}"

    def cardinality(self):
        r"""
        Return the number of partitions with parts in ``self``. Wraps GAP's
        ``NrRestrictedPartitions``.

        EXAMPLES::

            sage: Partitions(15, parts_in=[2,3,7]).cardinality()                        # needs sage.libs.gap
            5

        If you can use all parts 1 through `n`, we'd better get `p(n)`::

            sage: (Partitions(20, parts_in=[1..20]).cardinality()                       # needs sage.libs.gap
            ....:   == Partitions(20).cardinality())
            True

        TESTS:

        Let's check the consistency of GAP's function and our own
        algorithm that actually generates the partitions::

            sage: # needs sage.libs.gap
            sage: ps = Partitions(15, parts_in=[1,2,3])
            sage: ps.cardinality() == len(ps.list())
            True
            sage: ps = Partitions(15, parts_in=[])
            sage: ps.cardinality() == len(ps.list())
            True
            sage: ps = Partitions(3000, parts_in=[50,100,500,1000])
            sage: ps.cardinality() == len(ps.list())
            True
            sage: ps = Partitions(10, parts_in=[3,6,9])
            sage: ps.cardinality() == len(ps.list())
            True
            sage: ps = Partitions(0, parts_in=[1,2])
            sage: ps.cardinality() == len(ps.list())
            True
        """
        # GAP complains if you give it an empty list
        if self.parts:
            from sage.libs.gap.libgap import libgap
            return ZZ(libgap.NrRestrictedPartitions(ZZ(self.n), self.parts))
        return Integer(self.n == 0)

    def first(self):
        """
        Return the lexicographically first partition of a positive
        integer `n` with the specified parts, or ``None`` if no such
        partition exists.

        EXAMPLES::

            sage: Partitions(9, parts_in=[3,4]).first()
            [3, 3, 3]
            sage: Partitions(6, parts_in=[1..6]).first()
            [6]
            sage: Partitions(30, parts_in=[4,7,8,10,11]).first()
            [11, 11, 8]
        """
        try:
            return self.element_class(self, self._findfirst(self.n, self.parts[:]))
        except TypeError:
            return None

    def _findfirst(self, n, parts):
        """
        TESTS::

            sage: p = Partitions(9, parts_in=[3,4])
            sage: p._findfirst(p.n, p.parts[:])
            [3, 3, 3]
            sage: p._findfirst(0, p.parts[:])
            []
            sage: p._findfirst(p.n, [10])
        """
        if n == 0:
            return []
        else:
            while parts:
                p = parts.pop()
                for k in range(n.quo_rem(p)[0], 0, -1):
                    try:
                        return k * [p] + self._findfirst(n - k * p, parts[:])
                    except TypeError:
                        pass

    def last(self):
        """
        Return the lexicographically last partition of the positive
        integer `n` with the specified parts, or ``None`` if no such
        partition exists.

        EXAMPLES::

            sage: Partitions(15, parts_in=[2,3]).last()
            [3, 2, 2, 2, 2, 2, 2]
            sage: Partitions(30, parts_in=[4,7,8,10,11]).last()
            [7, 7, 4, 4, 4, 4]
            sage: Partitions(10, parts_in=[3,6]).last() is None
            True
            sage: Partitions(50, parts_in=[11,12,13]).last()
            [13, 13, 12, 12]
            sage: Partitions(30, parts_in=[4,7,8,10,11]).last()
            [7, 7, 4, 4, 4, 4]

        TESTS::

            sage: Partitions(6, parts_in=[1..6]).last()
            [1, 1, 1, 1, 1, 1]
            sage: Partitions(0, parts_in=[]).last()
            []
            sage: Partitions(50, parts_in=[11,12]).last() is None
            True
        """
        try:
            return self.element_class(self, self._findlast(self.n, self.parts))
        except TypeError:
            return None

    def _findlast(self, n, parts):
        """
        Return the lexicographically largest partition of `n` using the
        given parts, or ``None`` if no such partition exists. This function
        is not intended to be called directly.

        INPUT:

        - ``n`` -- nonnegative integer

        - ``parts`` -- a sorted list of positive integers

        OUTPUT:

        A list of integers in weakly decreasing order, or ``None``. The
        output is just a list, not a partition object.

        EXAMPLES::

            sage: ps = Partitions(1, parts_in=[1])
            sage: ps._findlast(15, [2,3])
            [3, 2, 2, 2, 2, 2, 2]
            sage: ps._findlast(9, [2,4]) is None
            True
            sage: ps._findlast(0, [])
            []
            sage: ps._findlast(100, [9,17,31])
            [31, 17, 17, 17, 9, 9]
        """
        if n < 0:
            return None
        elif n == 0:
            return []
        elif parts:
            p = parts[0]
            q, r = n.quo_rem(p)
            if r == 0:
                return [p] * q
            # If the smallest part doesn't divide n, try using the next
            # largest part
            else:
                for i, p in enumerate(parts[1:]):
                    rest = self._findlast(n - p, parts[:i + 2])
                    if rest is not None:
                        return [p] + rest
        # If we get to here, nothing ever worked, so there's no such
        # partitions, and we return None.
        return None

    def __iter__(self):
        """
        An iterator through the partitions of `n` with all parts belonging
        to a particular set.

        EXAMPLES::

            sage: [x for x in Partitions(5, parts_in=[1,2,3])]
            [[3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        """
        for p in self._other_iterator(self.n, self.parts):
            yield self.element_class(self, p)

    def _fast_iterator(self, n, parts):
        """
        A fast iterator for the partitions of `n` which returns lists and
        not partition types. This function is not intended to be called
        directly.

        INPUT:

        - ``n`` -- nonnegative integer

        - ``parts`` -- list of parts to use. This list will be
          destroyed, so pass things here with ``foo[:]`` (or something
          equivalent) if you want to preserve your list. In particular,
          the ``__iter__`` method needs to use ``self.parts[:]``, or else we
          forget which parts we're using!

        OUTPUT:

        A generator object for partitions of `n` with parts in
        ``parts``.

        If the parts in ``parts`` are sorted in increasing order, this
        function returns weakly decreasing lists. If ``parts`` is not
        sorted, your lists won't be, either.

        EXAMPLES::

            sage: P = Partitions(4, parts_in=[2,4])
            sage: it = P._fast_iterator(4, [2,4])
            sage: next(it)
            [4]
            sage: type(_)
            <class 'list'>
        """
        if n == 0:
            yield []
        else:
            while parts:
                p = parts.pop()
                for k in range(n.quo_rem(p)[0], 0, -1):
                    for q in self._fast_iterator(n - k * p, parts[:]):
                        yield k * [p] + q

    def _other_iterator(self, n, parts):
        """
        A fast iterator for the partitions of `n` which returns lists and
        not partition types. This function is not intended to be called
        directly.

        INPUT:

        - ``n`` -- nonnegative integer

        - ``parts`` -- list of parts to use

        OUTPUT: a generator object for partitions of `n` with parts in
        ``parts``

        EXAMPLES::

            sage: P = Partitions(4, parts_in=[2,4])
            sage: it = P._other_iterator(4, [2,4])
            sage: next(it)
            [4]
            sage: type(_)
            <class 'list'>
        """
        sorted_parts = sorted(parts, reverse=True)
        for vec in weighted_iterator_fast(n, sorted_parts):
            yield sum(([pi] * multi
                       for pi, multi in zip(sorted_parts, vec)), [])


class Partitions_starting(Partitions):
    """
    All partitions with a given start.
    """

    @staticmethod
    def __classcall_private__(cls, n, starting_partition):
        """
        Normalize the input to ensure a unique representation.

        TESTS::

            sage: P = Partitions(4, starting=[2,1])
            sage: P2 = Partitions(4, starting=[2,1])
            sage: P is P2
            True
        """
        starting_partition = Partition(starting_partition)
        return super().__classcall__(cls, Integer(n),
                                     starting_partition)

    def __init__(self, n, starting_partition):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Partitions(3, starting=[2,1])
            Partitions of the integer 3 starting with [2, 1]
            sage: Partitions(3, starting=[2,1]).list()
            [[2, 1], [1, 1, 1]]

            sage: Partitions(7, starting=[2,2,1]).list()
            [[2, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]]

            sage: Partitions(7, starting=[3,2]).list()
            [[3, 1, 1, 1, 1],
             [2, 2, 2, 1],
             [2, 2, 1, 1, 1],
             [2, 1, 1, 1, 1, 1],
             [1, 1, 1, 1, 1, 1, 1]]

            sage: Partitions(4, starting=[3,2]).list()
            [[3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

            sage: Partitions(3, starting=[1,1]).list()
            []

        TESTS::

            sage: p = Partitions(3, starting=[2,1])
            sage: TestSuite(p).run()
        """
        Partitions.__init__(self)
        self.n = n
        self._starting = starting_partition

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Partitions(3, starting=[2,1]) # indirect doctest
            Partitions of the integer 3 starting with [2, 1]
        """
        return f"Partitions of the integer {self.n} starting with {self._starting}"

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        EXAMPLES::

            sage: p = Partitions(3, starting=[2,1])
            sage: [1,1] in p
            False
            sage: [2,1] in p
            True
            sage: [1,1,1] in p
            True
            sage: [3] in p
            False
        """
        return x in Partitions_n(self.n) and x <= self._starting

    def first(self):
        """
        Return the first partition in ``self``.

        EXAMPLES::

            sage: Partitions(3, starting=[2,1]).first()
            [2, 1]
            sage: Partitions(3, starting=[1,1,1]).first()
            [1, 1, 1]
            sage: Partitions(3, starting=[1,1]).first()
            False
            sage: Partitions(3, starting=[3,1]).first()
            [3]
            sage: Partitions(3, starting=[2,2]).first()
            [2, 1]
        """
        if sum(self._starting) == self.n:
            return self._starting

        if (k := self._starting.size()) < self.n:
            mu = list(self._starting) + [1] * (self.n - k)
            return next(Partition(mu))

        # if self._starting.size() > self.n:
        return self.element_class(self, Partitions(self.n, outer=self._starting).first())

    def next(self, part):
        """
        Return the next partition after ``part`` in ``self``.

        EXAMPLES::

            sage: Partitions(3, starting=[2,1]).next(Partition([2,1]))
            [1, 1, 1]
        """
        return next(part)


class Partitions_ending(Partitions):
    """
    All partitions with a given ending.
    """

    @staticmethod
    def __classcall_private__(cls, n, ending_partition):
        """
        Normalize the input to ensure a unique representation.

        TESTS::

            sage: P = Partitions(4)
            sage: P2 = Partitions(4)
            sage: P is P2
            True
        """
        ending_partition = Partition(ending_partition)
        return super().__classcall__(cls, Integer(n),
                                     ending_partition)

    def __init__(self, n, ending_partition):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Partitions(4, ending=[1,1,1,1]).list()
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
            sage: Partitions(4, ending=[2,2]).list()
            [[4], [3, 1], [2, 2]]
            sage: Partitions(4, ending=[4]).list()
            [[4]]
            sage: Partitions(4, ending=[5]).list()
            []

        TESTS::

            sage: p = Partitions(4, ending=[1,1,1,1])
            sage: TestSuite(p).run()
        """
        Partitions.__init__(self)
        self.n = n
        self._ending = ending_partition
        self._ending_size_is_not_same = (n != sum(self._ending))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Partitions(4, ending=[1,1,1,1]) # indirect doctest
            Partitions of the integer 4 ending with [1, 1, 1, 1]
        """
        return f"Partitions of the integer {self.n} ending with {self._ending}"

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        EXAMPLES::

            sage: p = Partitions(4, ending=[2,2])
            sage: [4] in p
            True
            sage: [2,1,1] in p
            False
            sage: [2,1] in p
            False
        """
        return x in Partitions_n(self.n) and x >= self._ending

    def first(self):
        """
        Return the first partition in ``self``.

        EXAMPLES::

            sage: Partitions(4, ending=[1,1,1,1]).first()
            [4]
            sage: Partitions(4, ending=[5]).first() is None
            True
        """
        if self._ending and self.n <= self._ending[0] and not (self.n == self._ending[0] and len(self._ending) == 1):
            return None
        return self.element_class(self, [self.n])

    def next(self, part):
        """
        Return the next partition after ``part`` in ``self``.

        EXAMPLES::
            sage: Partitions(4, ending=[1,1,1,1,1]).next(Partition([4]))
            [3, 1]
            sage: Partitions(4, ending=[3,2]).next(Partition([3,1])) is None
            True
            sage: Partitions(4, ending=[1,1,1,1]).next(Partition([4]))
            [3, 1]
            sage: Partitions(4, ending=[1,1,1,1]).next(Partition([1,1,1,1])) is None
            True
            sage: Partitions(4, ending=[3]).next(Partition([3,1])) is None
            True
        """
        # if we have passed the last Partition, there is no next partition
        if part == self._ending:
            return None

        # if self._ending is a different size, we should make the comparison
        mu = next(part)
        if self._ending_size_is_not_same and mu < self._ending:
            return None
        return mu


class PartitionsInBox(Partitions):
    r"""
    All partitions which fit in an `h \times w` box.

    EXAMPLES::

        sage: PartitionsInBox(2, 2)
        Integer partitions which fit in a 2 x 2 box
        sage: PartitionsInBox(2, 2).list()
        [[], [1], [1, 1], [2], [2, 1], [2, 2]]
    """

    def __init__(self, h, w):
        """
        Initialize ``self``.

        TESTS::

            sage: p = PartitionsInBox(2,2)
            sage: TestSuite(p).run()
        """
        Partitions.__init__(self)
        self.h = h
        self.w = w

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionsInBox(2,2) # indirect doctest
            Integer partitions which fit in a 2 x 2 box
        """
        return f"Integer partitions which fit in a {self.h} x {self.w} box"

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        EXAMPLES::

            sage: [] in PartitionsInBox(2,2)
            True
            sage: [2,1] in PartitionsInBox(2,2)
            True
            sage: [3,1] in PartitionsInBox(2,2)
            False
            sage: [2,1,1] in PartitionsInBox(2,2)
            False
            sage: [3,1] in PartitionsInBox(3, 2)
            False
            sage: [3,1] in PartitionsInBox(2, 3)
            True
        """
        return x in _Partitions and len(x) <= self.h \
            and (len(x) == 0 or x[0] <= self.w)

    def list(self):
        """
        Return a list of all the partitions inside a box of height `h` and
        width `w`.

        EXAMPLES::

            sage: PartitionsInBox(2,2).list()
            [[], [1], [1, 1], [2], [2, 1], [2, 2]]
            sage: PartitionsInBox(2,3).list()
            [[], [1], [1, 1], [2], [2, 1], [2, 2], [3], [3, 1], [3, 2], [3, 3]]

        TESTS:

        Check :issue:`10890`::

            sage: type(PartitionsInBox(0,0)[0])
            <class 'sage.combinat.partition.PartitionsInBox_with_category.element_class'>
        """
        h = self.h
        w = self.w
        if h == 0:
            return [self.element_class(self, [])]
        else:
            l = [[i] for i in range(w + 1)]

            def add(x):
                return [x + [i] for i in range(x[-1] + 1)]

            for i in range(h-1):
                new_list = []
                for element in l:
                    new_list += add(element)
                l = new_list

            return [self.element_class(self, [x for x in p if x != 0]) for p in l]

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: PartitionsInBox(2, 3).cardinality()
            10

        TESTS:

        Check the corner case::

            sage: PartitionsInBox(0, 0).cardinality()
            1

            sage: PartitionsInBox(0, 1).cardinality()
            1

            sage: all(PartitionsInBox(a, b).cardinality() ==
            ....:     len(PartitionsInBox(a, b).list())
            ....:     for a in range(6) for b in range(6))
            True
        """
        return binomial(self.h + self.w, self.w)


class Partitions_constraints(IntegerListsLex):
    """
    For unpickling old constrained ``Partitions_constraints`` objects created
    with sage <= 3.4.1. See :class:`Partitions`.
    """

    def __setstate__(self, data):
        r"""
        TESTS::

            sage: dmp = b'x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+H,*\xc9,\xc9\xcc\xcf\xe3\n\x80\xb1\x8a\xe3\x93\x81DIQbf^I1W!\xa3fc!Sm!\xb3F(7\x92x!Km!k(GnbE<\xc8\x88B6\x88\xb9E\x99y\xe9\xc5z@\x05\xa9\xe9\xa9E\\\xb9\x89\xd9\xa9\xf10N!{(\xa3QkP!Gq(c^\x06\x90c\x0c\xe4p\x96&\xe9\x01\x00\xc2\xe53\xfd'
            sage: sp = loads(dmp); sp
            Integer lists of sum 3 satisfying certain constraints
            sage: sp.list()
            [[2, 1], [1, 1, 1]]
        """
        n = data['n']
        self.__class__ = Partitions_with_constraints
        constraints = {'max_slope': 0,
                       'min_part': 1}
        constraints.update(data['constraints'])
        self.__init__(n, **constraints)


class Partitions_with_constraints(IntegerListsLex):
    """
    Partitions which satisfy a set of constraints.

    EXAMPLES::

        sage: P = Partitions(6, inner=[1,1], max_slope=-1)
        sage: list(P)
        [[5, 1], [4, 2], [3, 2, 1]]

    TESTS::

        sage: P = Partitions(6, min_part=2, max_slope=-1)
        sage: TestSuite(P).run()

    Test that :issue:`15525` is fixed::

        sage: loads(dumps(P)) == P
        True
    """
#    def __init__(self, n, **kwargs):
#        """
#        Initialize ``self``.
#        """
#        IntegerListsLex.__init__(self, n, **kwargs)

    Element = Partition
    options = Partitions.options


######################
# Regular Partitions #
######################

class RegularPartitions(Partitions):
    r"""
    Base class for `\ell`-regular partitions.

    Let `\ell` be a positive integer. A partition `\lambda` is
    `\ell`-*regular* if `m_i < \ell` for all `i`, where `m_i` is the
    multiplicity of `i` in `\lambda`.

    .. NOTE::

        This is conjugate to the notion of `\ell`-*restricted* partitions,
        where the difference between any two consecutive
        parts is `< \ell`.

    INPUT:

    - ``ell`` -- the positive integer `\ell`
    - ``is_infinite`` -- boolean; if the subset of `\ell`-regular
      partitions is infinite
    """

    def __init__(self, ell, is_infinite=False):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=2)
            sage: TestSuite(P).run()
        """
        self._ell = ell
        Partitions.__init__(self, is_infinite)

    def ell(self):
        r"""
        Return the value `\ell`.

        EXAMPLES::

            sage: P = Partitions(regular=2)
            sage: P.ell()
            2
        """
        return self._ell

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions(regular=3)
            sage: [5] in P
            True
            sage: [] in P
            True
            sage: [3, 3, 2, 2] in P
            True
            sage: [3, 3, 3, 1] in P
            False
            sage: [4, 0, 0, 0, 0, 0] in P
            True
            sage: Partition([4,2,2,1]) in P
            True
            sage: Partition([4,2,2,2]) in P
            False
            sage: Partition([10,1]) in P
            True
        """
        if not Partitions.__contains__(self, x):
            return False
        if isinstance(x, Partition):
            return max(x.to_exp() + [0]) < self._ell
        return all(x.count(i) < self._ell for i in set(x) if i > 0)

    def _fast_iterator(self, n, max_part):
        """
        A fast (recursive) iterator which returns a list.

        EXAMPLES::

            sage: P = Partitions(regular=3)
            sage: list(P._fast_iterator(5, 5))
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]
            sage: list(P._fast_iterator(5, 3))
            [[3, 2], [3, 1, 1], [2, 2, 1]]
            sage: list(P._fast_iterator(5, 6))
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]
        """
        if n == 0:
            yield []
            return

        max_part = min(n, max_part)
        bdry = self._ell - 1

        for i in reversed(range(1, max_part + 1)):
            for p in self._fast_iterator(n - i, i):
                if p.count(i) < bdry:
                    yield [i] + p


class RegularPartitions_all(RegularPartitions):
    r"""
    The class of all `\ell`-regular partitions.

    INPUT:

    - ``ell`` -- the positive integer `\ell`

    .. SEEALSO::

        :class:`~sage.combinat.partition.RegularPartitions`
    """

    def __init__(self, ell):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=4)
            sage: TestSuite(P).run()

        1-regular partitions::

            sage: P = Partitions(regular=1)
            sage: P in FiniteEnumeratedSets()
            True
            sage: TestSuite(P).run()
        """
        RegularPartitions.__init__(self, ell, bool(ell > 1))

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.partition import RegularPartitions_all
            sage: RegularPartitions_all(3)
            3-Regular Partitions
        """
        return f"{self._ell}-Regular Partitions"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=3)
            sage: it = P.__iter__()
            sage: [next(it) for x in range(10)]
            [[], [1], [2], [1, 1], [3], [2, 1], [4], [3, 1], [2, 2], [2, 1, 1]]

        Check that 1-regular partitions works (:issue:`20584`)::

            sage: P = Partitions(regular=1)
            sage: list(P)
            [[]]
        """
        if self._ell == 1:
            yield self.element_class(self, [])
            return

        n = 0
        while True:
            for p in self._fast_iterator(n, n):
                yield self.element_class(self, p)
            n += 1


class RegularPartitions_truncated(RegularPartitions):
    r"""
    The class of `\ell`-regular partitions with max length `k`.

    INPUT:

    - ``ell`` -- the integer `\ell`
    - ``max_len`` -- integer; the maximum length

    .. SEEALSO::

        :class:`~sage.combinat.partition.RegularPartitions`
    """

    def __init__(self, ell, max_len):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=4, max_length=3)
            sage: TestSuite(P).run()
        """
        self._max_len = max_len
        RegularPartitions.__init__(self, ell, bool(ell > 1))

    def max_length(self):
        """
        Return the maximum length of the partitions of ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=4, max_length=3)
            sage: P.max_length()
            3
        """
        return self._max_len

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions(regular=4, max_length=3)
            sage: [3, 3, 3] in P
            True
            sage: [] in P
            True
            sage: [4, 2, 1, 1] in P
            False
        """
        return len(x) <= self._max_len and RegularPartitions.__contains__(self, x)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.partition import RegularPartitions_truncated
            sage: RegularPartitions_truncated(4, 3)
            4-Regular Partitions with max length 3
        """
        return f"{self._ell}-Regular Partitions with max length {self._max_len}"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=3, max_length=2)
            sage: it = P.__iter__()
            sage: [next(it) for x in range(10)]
            [[], [1], [2], [1, 1], [3], [2, 1], [4], [3, 1], [2, 2], [5]]

        Check that 1-regular partitions works (:issue:`20584`)::

            sage: P = Partitions(regular=1, max_length=2)
            sage: list(P)
            [[]]
        """
        if self._ell == 1:
            yield self.element_class(self, [])
            return

        n = 0
        while True:
            for p in self._fast_iterator(n, n):
                yield self.element_class(self, p)
            n += 1

    def _fast_iterator(self, n, max_part, depth=0):
        """
        A fast (recursive) iterator which returns a list.

        EXAMPLES::

            sage: P = Partitions(regular=2, max_length=2)
            sage: list(P._fast_iterator(5, 5))
            [[5], [4, 1], [3, 2]]
            sage: list(P._fast_iterator(5, 3))
            [[3, 2]]
            sage: list(P._fast_iterator(5, 6))
            [[5], [4, 1], [3, 2]]
        """
        if n == 0 or depth >= self._max_len:
            yield []
            return

        # Special case
        if depth + 1 == self._max_len:
            if max_part >= n:
                yield [n]
            return

        max_part = min(n, max_part)
        bdry = self._ell - 1

        for i in reversed(range(1, max_part + 1)):
            for p in self._fast_iterator(n - i, i, depth + 1):
                if p.count(i) < bdry:
                    yield [i] + p


class RegularPartitions_bounded(RegularPartitions):
    r"""
    The class of `\ell`-regular `k`-bounded partitions.

    INPUT:

    - ``ell`` -- the integer `\ell`
    - ``k`` -- integer; the value `k`

    .. SEEALSO::

        :class:`~sage.combinat.partition.RegularPartitions`
    """

    def __init__(self, ell, k):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=4, max_part=3)
            sage: TestSuite(P).run()

        1-regular partitions::

            sage: P = Partitions(regular=1, max_part=3)
            sage: P in FiniteEnumeratedSets()
            True
            sage: TestSuite(P).run()
        """
        self.k = k
        RegularPartitions.__init__(self, ell, False)

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions(regular=4, max_part=3)
            sage: [3, 3, 3] in P
            True
            sage: [] in P
            True
            sage: [4, 2, 1] in P
            False
        """
        return len(x) == 0 or (x[0] <= self.k and RegularPartitions.__contains__(self, x))

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.partition import RegularPartitions_bounded
            sage: RegularPartitions_bounded(4, 3)
            4-Regular 3-Bounded  Partitions
        """
        return f"{self._ell}-Regular {self.k}-Bounded Partitions"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: P = Partitions(regular=2, max_part=3)
            sage: list(P)
            [[3, 2, 1], [3, 2], [3, 1], [3], [2, 1], [2], [1], []]

        Check that 1-regular partitions works (:issue:`20584`)::

            sage: P = Partitions(regular=1, max_part=3)
            sage: list(P)
            [[]]
        """
        k = self.k
        for n in reversed(range(k*(k+1)/2 * self._ell)):
            for p in self._fast_iterator(n, k):
                yield self.element_class(self, p)


class RegularPartitions_n(RegularPartitions, Partitions_n):
    r"""
    The class of `\ell`-regular partitions of `n`.

    INPUT:

    - ``n`` -- the integer `n` to partition
    - ``ell`` -- the integer `\ell`

    .. SEEALSO::

        :class:`~sage.combinat.partition.RegularPartitions`
    """

    def __init__(self, n, ell):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(5, regular=3)
            sage: TestSuite(P).run()

        1-regular partitions::

            sage: P = Partitions(5, regular=1)
            sage: TestSuite(P).run()
        """
        RegularPartitions.__init__(self, ell)
        Partitions_n.__init__(self, n)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.partition import RegularPartitions_n
            sage: RegularPartitions_n(3, 5)
            5-Regular Partitions of the integer 3
        """
        return f"{self._ell}-Regular Partitions of the integer {self.n}"

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions(5, regular=3)
            sage: [3, 1, 1] in P
            True
            sage: [3, 2, 1] in P
            False
        """
        return RegularPartitions.__contains__(self, x) and sum(x) == self.n

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: P = Partitions(5, regular=3)
            sage: list(P)
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]
        """
        for p in self._fast_iterator(self.n, self.n):
            yield self.element_class(self, p)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: P = Partitions(5, regular=3)
            sage: P.cardinality()                                                       # needs sage.libs.flint
            5
            sage: P = Partitions(5, regular=6)
            sage: P.cardinality()                                                       # needs sage.libs.flint
            7
            sage: P.cardinality() == Partitions(5).cardinality()                        # needs sage.libs.flint
            True

        TESTS:

        Check the corner case::

            sage: P = Partitions(0, regular=3)
            sage: P.cardinality()                                                       # needs sage.libs.flint
            1

        Check for 1-regular partitions::

            sage: P = Partitions(0, regular=1)
            sage: P.cardinality()                                                       # needs sage.libs.flint
            1
            sage: P = Partitions(5, regular=1)
            sage: P.cardinality()                                                       # needs sage.libs.flint
            0
        """
        if self._ell > self.n:
            return Partitions_n.cardinality(self)
        return ZZ.sum(1 for x in self)

    def _an_element_(self):
        """
        Return a partition in ``self``.

        EXAMPLES::

            sage: P = Partitions(5, regular=2)
            sage: P._an_element_()
            [4, 1]

            sage: P = Partitions(0, regular=1)
            sage: P._an_element_()
            []

            sage: P = Partitions(5, regular=1)
            sage: P._an_element_()
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        if self._ell == 1 and self.n > 0:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError
        return Partitions_n._an_element_(self)


######################
# Ordered Partitions #
######################

class OrderedPartitions(Partitions):
    """
    The class of ordered partitions of `n`. If `k` is specified, then this
    contains only the ordered partitions of length `k`.

    An *ordered partition* of a nonnegative integer `n` means a list of
    positive integers whose sum is `n`. This is the same as a composition
    of `n`.

    .. NOTE::

       It is recommended that you use :meth:`Compositions` instead as
       :meth:`OrderedPartitions` wraps GAP.

    EXAMPLES::

        sage: OrderedPartitions(3)
        Ordered partitions of 3
        sage: OrderedPartitions(3).list()                                               # needs sage.libs.gap
        [[3], [2, 1], [1, 2], [1, 1, 1]]
        sage: OrderedPartitions(3,2)
        Ordered partitions of 3 of length 2
        sage: OrderedPartitions(3,2).list()                                             # needs sage.libs.gap
        [[2, 1], [1, 2]]

        sage: OrderedPartitions(10, k=2).list()                                         # needs sage.libs.gap
        [[9, 1], [8, 2], [7, 3], [6, 4], [5, 5], [4, 6], [3, 7], [2, 8], [1, 9]]
        sage: OrderedPartitions(4).list()                                               # needs sage.libs.gap
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
    """
    @staticmethod
    def __classcall_private__(cls, n, k=None):
        """
        Normalize the input to ensure a unique representation.

        TESTS::

            sage: P = OrderedPartitions(3,2)
            sage: P2 = OrderedPartitions(3,2)
            sage: P is P2
            True
        """
        if k is not None:
            k = Integer(k)
        return super().__classcall__(cls, Integer(n), k)

    def __init__(self, n, k):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: o = OrderedPartitions(4,2)

        TESTS::

            sage: TestSuite( OrderedPartitions(5,3) ).run()                             # needs sage.libs.gap
        """
        Partitions.__init__(self)
        self.n = n
        self.k = k

    def __contains__(self, x):
        """
        Check to see if ``x`` is an element of ``self``.

        EXAMPLES::

            sage: o = OrderedPartitions(4,2)
            sage: [2,1] in o
            False
            sage: [2,2] in o
            True
            sage: [1,2,1] in o
            False
        """
        C = composition.Compositions(self.n, length=self.k)
        return C(x) in composition.Compositions(self.n, length=self.k)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: OrderedPartitions(3) # indirect doctest
            Ordered partitions of 3
            sage: OrderedPartitions(3,2) # indirect doctest
            Ordered partitions of 3 of length 2
        """
        string = "Ordered partitions of %s" % self.n
        if self.k is not None:
            string += " of length %s" % self.k
        return string

    def list(self):
        """
        Return a list of partitions in ``self``.

        EXAMPLES::

            sage: OrderedPartitions(3).list()                                           # needs sage.libs.gap
            [[3], [2, 1], [1, 2], [1, 1, 1]]
            sage: OrderedPartitions(3,2).list()                                         # needs sage.libs.gap
            [[2, 1], [1, 2]]
        """
        from sage.libs.gap.libgap import libgap
        n = self.n
        k = self.k
        if k is None:
            ans = libgap.OrderedPartitions(ZZ(n))
        else:
            ans = libgap.OrderedPartitions(ZZ(n), ZZ(k))
        result = ans.sage()
        result.reverse()
        return result

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: OrderedPartitions(3).cardinality()
            4
            sage: OrderedPartitions(3,2).cardinality()
            2
            sage: OrderedPartitions(10,2).cardinality()
            9
            sage: OrderedPartitions(15).cardinality()
            16384
        """
        from sage.libs.gap.libgap import libgap
        n = self.n
        k = self.k
        if k is None:
            ans = libgap.NrOrderedPartitions(n)
        else:
            ans = libgap.NrOrderedPartitions(n, k)
        return ZZ(ans)


######################################
# Partitions_parts_length_restricted #
######################################

class Partitions_parts_length_restricted(Partitions):
    r"""
    The class of all integer partitions having parts and length in a
    given range.

    This class is strictly more general than
    :class:`PartitionsGreatestLE`.

    INPUT:

    - ``n`` -- the size of the partition
    - ``min_part`` -- the bound on the smallest part
    - ``max_part`` -- the bound on the largest part
    - ``min_length`` -- the lower bound on the number of parts
    - ``max_length`` -- the upper bound on the number of parts

    EXAMPLES::

        sage: from sage.combinat.partition import Partitions_parts_length_restricted
        sage: Partitions_parts_length_restricted(10, 2, 5, 0, 10)
        Partitions of 10 whose parts are between 2 and 5
        sage: list(Partitions_parts_length_restricted(9, 2, 4, 3, 4))
        [[4, 3, 2], [3, 3, 3], [3, 2, 2, 2]]

        sage: [4,3,2,1] in Partitions_parts_length_restricted(10, 2, 10, 0, 10)
        False
        sage: [2,2,2,2,2] in Partitions_parts_length_restricted(10, 2, 10, 0, 10)
        True
    """
    def __init__(self, n, min_part, max_part, min_length, max_length):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.partition import Partitions_parts_length_restricted
            sage: p = Partitions_parts_length_restricted(10, 2, 5, 3, 4)
            sage: TestSuite(p).run()
        """
        Partitions.__init__(self)
        self._n = n
        self._min_part = min_part
        self._max_part = max_part
        self._min_length = min_length
        self._max_length = max_length

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.partition import Partitions_parts_length_restricted
            sage: Partitions_parts_length_restricted(9, 2, 9, 0, 9)
            Partitions of 9 whose parts are at least 2
            sage: Partitions_parts_length_restricted(9, 2, 9, 3, 5)
            Partitions of 9 having length between 3 and 5 and whose parts are at least 2
        """
        if not self._min_length and self._max_length == self._n:
            length_str = ""
        elif self._min_length == self._max_length:
            length_str = f"having length {self._min_length}"
        elif not self._min_length:
            length_str = f"having length at most {self._max_length}"
        elif self._max_length == self._n:
            length_str = f"having length at least {self._min_length}"
        else:
            length_str = f"having length between {self._min_length} and {self._max_length}"

        if self._min_part == ZZ.one() and self._max_part == self._n:
            parts_str = ""
        elif self._min_part == self._max_part:
            parts_str = f"having parts equal to {self._min_part}"
        elif self._min_part == ZZ.one():
            parts_str = f"whose parts are at most {self._max_part}"
        elif self._max_part == self._n:
            parts_str = f"whose parts are at least {self._min_part}"
        else:
            parts_str = f"whose parts are between {self._min_part} and {self._max_part}"

        if length_str:
            if parts_str:
                return f"Partitions of {self._n} " + length_str + " and " + parts_str
            return f"Partitions of {self._n} " + length_str
        if parts_str:
            return f"Partitions of {self._n} " + parts_str
        return f"Partitions of {self._n}"

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        TESTS::

            sage: P = Partitions(10, min_part=2, max_part=5, min_length=2, max_length=4)
            sage: Partition([]) in P
            False

            sage: Partition([3]) in P
            False

            sage: Partition([5, 3, 2]) in P
            True
        """
        try:
            mu = Partition(x)
        except ValueError:
            return False
        return (mu.size() == self._n
                and (not mu
                     or (min(mu) >= self._min_part
                         and max(mu) <= self._max_part
                         and self._min_length <= len(mu) <= self._max_length)))

    def __iter__(self):
        """
        Iterator over the set of partitions in ``self``.

        EXAMPLES::

            sage: list(Partitions(9, min_part=2, max_part=4, min_length=3, max_length=4))
            [[4, 3, 2], [3, 3, 3], [3, 2, 2, 2]]
        """
        yield from IntegerListsLex(self._n, max_slope=0,
                                   min_part=self._min_part,
                                   max_part=self._max_part,
                                   min_length=self._min_length,
                                   max_length=self._max_length,
                                   element_constructor=Partition)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: from sage.combinat.partition import Partitions_parts_length_restricted
            sage: list(Partitions_parts_length_restricted(9, 3, 9, 0, 2))
            [[9], [6, 3], [5, 4]]
            sage: Partitions_parts_length_restricted(9, 3, 9, 0, 2).cardinality()
            3

        TESTS::

            sage: from itertools import product
            sage: P = Partitions
            sage: all(P(n, min_part=a, max_part=b, min_length=k, max_length=m).cardinality()
            ....:     == len(list(P(n, min_part=a, max_part=b, min_length=k, max_length=m)))
            ....:     for n, a, b, k, m in product(range(-1, 5), repeat=5))
            True
        """
        n = self._n
        a = self._min_part - 1
        if not self._min_length and self._max_length == n and not a:
            # unrestricted length, parts smaller max_part
            return ZZ.sum(number_of_partitions_length(n, i)
                          for i in range(self._max_part + 1))

        m = self._max_part - a
        return ZZ.sum(number_of_partitions_length_max_part(n - a * ell, ell, m)
                      for ell in range(self._min_length, self._max_length + 1))


##########################
# Partitions Greatest LE #
##########################

class PartitionsGreatestLE(UniqueRepresentation, IntegerListsLex):
    r"""
    The class of all (unordered) "restricted" partitions of the
    integer `n` having parts less than or equal to the integer `k`.

    EXAMPLES::

        sage: PartitionsGreatestLE(10, 2)
        Partitions of 10 having parts less than or equal to 2
        sage: PartitionsGreatestLE(10, 2).list()
        [[2, 2, 2, 2, 2],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

        sage: [4,3,2,1] in PartitionsGreatestLE(10, 2)
        False
        sage: [2,2,2,2,2] in PartitionsGreatestLE(10, 2)
        True
        sage: PartitionsGreatestLE(10, 2).first().parent()
        Partitions...
    """
    def __init__(self, n, k):
        """
        Initialize ``self``.

        TESTS::

            sage: p = PartitionsGreatestLE(10, 2)
            sage: p.n, p.k
            (10, 2)
            sage: TestSuite(p).run()
        """
        IntegerListsLex.__init__(self, n, max_slope=0, min_part=1, max_part=k)
        self.n = n
        self.k = k

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: PartitionsGreatestLE(10, 2) # indirect doctest
            Partitions of 10 having parts less than or equal to 2
        """
        return f"Partitions of {self.n} having parts less than or equal to {self.k}"

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: PartitionsGreatestLE(9, 5).cardinality()                              # needs sage.libs.gap
            23

        TESTS::

            sage: all(PartitionsGreatestLE(n, a).cardinality() ==
            ....:     len(list(PartitionsGreatestLE(n, a)))
            ....:     for n in range(20) for a in range(6))
            True
        """
        return sum(number_of_partitions_length(self.n, i) for i in range(self.k + 1))

    Element = Partition
    options = Partitions.options


##########################
# Partitions Greatest EQ #
##########################

class PartitionsGreatestEQ(UniqueRepresentation, IntegerListsLex):
    """
    The class of all (unordered) "restricted" partitions of the integer `n`
    having all its greatest parts equal to the integer `k`.

    EXAMPLES::

        sage: PartitionsGreatestEQ(10, 2)
        Partitions of 10 having greatest part equal to 2
        sage: PartitionsGreatestEQ(10, 2).list()
        [[2, 2, 2, 2, 2],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 1, 1, 1, 1, 1, 1, 1, 1]]

        sage: [4,3,2,1] in PartitionsGreatestEQ(10, 2)
        False
        sage: [2,2,2,2,2] in PartitionsGreatestEQ(10, 2)
        True

    The empty partition has no maximal part, but it is contained in
    the set of partitions with any specified maximal part::

        sage: PartitionsGreatestEQ(0, 2).list()
        [[]]

    TESTS::

        sage: [1]*10 in PartitionsGreatestEQ(10, 2)
        False

        sage: PartitionsGreatestEQ(10, 2).first().parent()
        Partitions...
    """

    def __init__(self, n, k):
        """
        Initialize ``self``.

        TESTS::

            sage: p = PartitionsGreatestEQ(10, 2)
            sage: p.n, p.k
            (10, 2)
            sage: TestSuite(p).run()
        """
        IntegerListsLex.__init__(self, n, max_slope=0, max_part=k, floor=[k])
        self.n = n
        self.k = k

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: PartitionsGreatestEQ(10, 2) # indirect doctest
            Partitions of 10 having greatest part equal to 2
        """
        return f"Partitions of {self.n} having greatest part equal to {self.k}"

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: PartitionsGreatestEQ(10, 2).cardinality()
            5

        TESTS::

            sage: all(PartitionsGreatestEQ(n, a).cardinality() ==                       # needs sage.libs.flint
            ....:     len(PartitionsGreatestEQ(n, a).list())
            ....:     for n in range(20) for a in range(6))
            True
        """
        if not self.n:
            return 1
        return number_of_partitions_length(self.n, self.k)

    Element = Partition
    options = Partitions.options


#########################
# Restricted Partitions #
#########################

class RestrictedPartitions_generic(Partitions):
    r"""
    Base class for `\ell`-restricted partitions.

    Let `\ell` be a positive integer. A partition `\lambda` is
    `\ell`-*restricted* if `\lambda_i - \lambda_{i+1} < \ell` for all `i`,
    including rows of length 0.

    .. NOTE::

        This is conjugate to the notion of `\ell`-*regular* partitions,
        where the multiplicity of any parts is at most `\ell`.

    INPUT:

    - ``ell`` -- the positive integer `\ell`
    - ``is_infinite`` -- boolean; if the subset of `\ell`-restricted
      partitions is infinite
    """

    def __init__(self, ell, is_infinite=False):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(restricted=2)
            sage: TestSuite(P).run()
        """
        self._ell = ell
        Partitions.__init__(self, is_infinite)

    def ell(self):
        r"""
        Return the value `\ell`.

        EXAMPLES::

            sage: P = Partitions(restricted=2)
            sage: P.ell()
            2
        """
        return self._ell

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions(restricted=3)
            sage: [5] in P
            False
            sage: [2] in P
            True
            sage: [] in P
            True
            sage: [3, 3, 3, 3, 2, 2] in P
            True
            sage: [3, 3, 3, 1] in P
            True
            sage: [8, 3, 3, 1] in P
            False
            sage: [2, 0, 0, 0, 0, 0] in P
            True
            sage: Partition([4,2,2,1]) in P
            True
            sage: Partition([4,2,2,2]) in P
            True
            sage: Partition([6,6,6,6,4,3,2]) in P
            True
            sage: Partition([7,6,6,2]) in P
            False
            sage: Partition([6,5]) in P
            False
            sage: Partition([10,1]) in P
            False
            sage: Partition([3,3] + [1]*10) in P
            True
        """
        if not Partitions.__contains__(self, x):
            return False
        if x == []:
            return True
        return (all(x[i] - x[i+1] < self._ell for i in range(len(x)-1))
                and x[-1] < self._ell)

    def _fast_iterator(self, n, max_part):
        """
        A fast (recursive) iterator which returns a list.

        EXAMPLES::

            sage: P = Partitions(restricted=3)
            sage: list(P._fast_iterator(5, 5))
            [[3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: list(P._fast_iterator(5, 2))
            [[2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]

        TESTS::

            sage: for n in range(10):
            ....:     for ell in range(2, n):
            ....:         P_res = Partitions(n, restricted=ell)
            ....:         P_reg = Partitions(n, regular=ell)
            ....:         assert set(P_res) == set(p.conjugate() for p in P_reg)
        """
        if n == 0:
            yield []
            return

        max_part = min(n, max_part)

        for i in range(max_part, 0, -1):
            for p in self._fast_iterator(n-i, i):
                if (p and i - p[0] >= self._ell) or (not p and i >= self._ell):
                    break
                yield [i] + p


class RestrictedPartitions_all(RestrictedPartitions_generic):
    r"""
    The class of all `\ell`-restricted partitions.

    INPUT:

    - ``ell`` -- the positive integer `\ell`

    .. SEEALSO::

        :class:`~sage.combinat.partition.RestrictedPartitions_generic`
    """

    def __init__(self, ell):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(restricted=4)
            sage: TestSuite(P).run()
        """
        RestrictedPartitions_generic.__init__(self, ell, True)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.partition import RestrictedPartitions_all
            sage: RestrictedPartitions_all(3)
            3-Restricted Partitions
        """
        return f"{self._ell}-Restricted Partitions"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: P = Partitions(restricted=3)
            sage: it = P.__iter__()
            sage: [next(it) for x in range(10)]
            [[], [1], [2], [1, 1], [2, 1], [1, 1, 1],
             [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        """
        n = 0
        while True:
            for p in self._fast_iterator(n, n):
                yield self.element_class(self, p)
            n += 1


class RestrictedPartitions_n(RestrictedPartitions_generic, Partitions_n):
    r"""
    The class of `\ell`-restricted partitions of `n`.

    INPUT:

    - ``n`` -- the integer `n` to partition
    - ``ell`` -- the integer `\ell`

    .. SEEALSO::

        :class:`~sage.combinat.partition.RestrictedPartitions_generic`
    """

    def __init__(self, n, ell):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: P = Partitions(5, restricted=3)
            sage: TestSuite(P).run()
        """
        RestrictedPartitions_generic.__init__(self, ell)
        Partitions_n.__init__(self, n)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.partition import RestrictedPartitions_n
            sage: RestrictedPartitions_n(3, 5)
            5-Restricted Partitions of the integer 3
        """
        return f"{self._ell}-Restricted Partitions of the integer {self.n}"

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions(5, regular=3)
            sage: [3, 1, 1] in P
            True
            sage: [3, 2, 1] in P
            False
        """
        return RestrictedPartitions_generic.__contains__(self, x) and sum(x) == self.n

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: P = Partitions(5, restricted=3)
            sage: list(P)
            [[3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        """
        for p in self._fast_iterator(self.n, self.n):
            yield self.element_class(self, p)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: P = Partitions(5, restricted=3)
            sage: P.cardinality()                                                       # needs sage.libs.flint
            5
            sage: P = Partitions(5, restricted=6)
            sage: P.cardinality()                                                       # needs sage.libs.flint
            7
            sage: P.cardinality() == Partitions(5).cardinality()                        # needs sage.libs.flint
            True
        """
        if self._ell > self.n:
            return Partitions_n.cardinality(self)
        return ZZ.sum(ZZ.one() for x in self)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: P = Partitions(5, restricted=3)
            sage: P.an_element()
            [2, 1, 1, 1]

            sage: Partitions(0, restricted=3).an_element()
            []
            sage: Partitions(1, restricted=3).an_element()
            [1]
        """
        return self.element_class(self, Partitions_n._an_element_(self).conjugate())


#########################################################################

# partitions

def number_of_partitions(n, algorithm='default'):
    r"""
    Return the number of partitions of `n` with, optionally, at most `k`
    parts.

    The options of :meth:`number_of_partitions()` are being deprecated
    :issue:`13072` in favour of :meth:`Partitions_n.cardinality()` so that
    :meth:`number_of_partitions()` can become a stripped down version of
    the fastest algorithm available (currently this is using FLINT).

    INPUT:

    - ``n`` -- integer

    - ``algorithm`` -- (default: ``'default'``)
       [Will be deprecated except in Partition().cardinality() ]

       - ``'default'`` -- if ``k`` is not ``None``, then use Gap (very slow);
          if  ``k`` is ``None``, use FLINT

       - ``'flint'`` -- use FLINT

    EXAMPLES::

        sage: v = Partitions(5).list(); v
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        sage: len(v)
        7

    The input must be a nonnegative integer or a :exc:`ValueError` is raised.

    ::

        sage: number_of_partitions(-5)
        Traceback (most recent call last):
        ...
        ValueError: n (=-5) must be a nonnegative integer

    ::

        sage: # needs sage.libs.flint
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(3)
        3
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(40)
        37338
        sage: number_of_partitions(100)
        190569292
        sage: number_of_partitions(100000)
        27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519

    A generating function for the number of partitions `p_n` is given by the
    reciprocal of Euler's function:

    .. MATH::

        \sum_{n=0}^{\infty} p_n x^n = \prod_{k=1}^{\infty} \left(
        \frac{1}{1-x^k} \right).

    We use Sage to verify that the first several coefficients do
    instead agree::

        sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
        sage: prod([(1-q^k)^(-1) for k in range(1,9)])  # partial product of
        1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
        sage: [number_of_partitions(k) for k in range(2,10)]                            # needs sage.libs.flint
        [2, 3, 5, 7, 11, 15, 22, 30]

    REFERENCES:

    - :wikipedia:`Partition\_(number\_theory)`

    TESTS::

        sage: # needs sage.libs.flint
        sage: n = 500 + randint(0,500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1500 + randint(0,1500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 100000000 + randint(0,100000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0     # long time (4s on sage.math, 2011)
        True
    """
    n = ZZ(n)
    if n < 0:
        raise ValueError("n (=%s) must be a nonnegative integer" % n)
    elif n == 0:
        return ZZ.one()

    if algorithm == 'default':
        algorithm = 'flint'

    if algorithm == 'flint':
        return cached_number_of_partitions(n)

    raise ValueError("unknown algorithm '%s'" % algorithm)


def number_of_partitions_length(n, k, algorithm='hybrid'):
    r"""
    Return the number of partitions of `n` with length `k`.

    This is a wrapper for GAP's ``NrPartitions`` function.

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: from sage.combinat.partition import number_of_partitions_length
        sage: number_of_partitions_length(5, 2)
        2
        sage: number_of_partitions_length(10, 2)
        5
        sage: number_of_partitions_length(10, 4)
        9
        sage: number_of_partitions_length(10, 0)
        0
        sage: number_of_partitions_length(10, 1)
        1
        sage: number_of_partitions_length(0, 0)
        1
        sage: number_of_partitions_length(0, 1)
        0
    """
    if algorithm == 'hybrid':
        # Do the hybrid algorithm

        # Special relations between n and k
        if n < k:
            return ZZ.zero()
        if n == k and n >= 0:
            return ZZ.one()

        # Special case of n
        if n <= 0:
            # Note: we've already checked the case when n == k == 0
            return ZZ.zero()

        # Small values of k
        if k <= 0:
            return ZZ.zero()
        if k == 1:
            return ZZ.one()
        if k == 2:
            return n // 2

        # We have one column of length `k` and all (inner) partitions of
        #    size `n-k` can't have length more than `k`
        if n <= k*2:
            return number_of_partitions(n - k)

        # Fall back to GAP
    from sage.libs.gap.libgap import libgap
    return ZZ(libgap.NrPartitions(ZZ(n), ZZ(k)))


@cached_function
def number_of_partitions_length_max_part(n, k, b):
    r"""
    Return the number of partitions of `n` with exactly `k` parts and
    the largest part at most `b`.

    EXAMPLES::

        sage: from sage.combinat.partition import number_of_partitions_length_max_part
        sage: number_of_partitions_length_max_part(10, 5, 3)
        3
        sage: list(Partitions(10, length=5, max_part=3))
        [[3, 3, 2, 1, 1], [3, 2, 2, 2, 1], [2, 2, 2, 2, 2]]
    """
    if not n:
        if not k:
            return ZZ.one()
        return ZZ.zero()
    if not k or k > n or n > b * k:
        return ZZ.zero()
    if b >= n:
        return number_of_partitions_length(n, k)

    return ZZ.sum(number_of_partitions_length_max_part(n - m, k - 1, m)
                  for m in range(1, b+1))


##########
# issue 14225: Partitions() is frequently used, but only weakly cached.
# Hence, establish a strong reference to it.

_Partitions = Partitions()

# Rather than caching an under-used function I have cached the default
# number_of_partitions functions which is currently using FLINT.
# AM issue #13072
try:
    from sage.libs.flint.arith_sage import number_of_partitions as flint_number_of_partitions
    cached_number_of_partitions = cached_function(flint_number_of_partitions)
except ImportError:
    pass

# October 2012: fixing outdated pickles which use classes being deprecated
from sage.misc.persist import register_unpickle_override
from sage.combinat.partition_tuple import PartitionTuples_level_size
register_unpickle_override('sage.combinat.partition', 'PartitionTuples_nk', PartitionTuples_level_size)
register_unpickle_override('sage.combinat.partition', 'Partition_class', Partition)
register_unpickle_override('sage.combinat.partition', 'OrderedPartitions_nk', OrderedPartitions)
register_unpickle_override('sage.combinat.partition', 'PartitionsInBox_hw', PartitionsInBox)
register_unpickle_override('sage.combinat.partition', 'PartitionsGreatestLE_nk', PartitionsGreatestLE)
register_unpickle_override('sage.combinat.partition', 'PartitionsGreatestEQ_nk', PartitionsGreatestEQ)
