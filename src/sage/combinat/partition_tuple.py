r"""
Partition tuples

A :class:`PartitionTuple` is a tuple of partitions. That is, an ordered
`k`-tuple of partitions `\mu=(\mu^{(1)},\mu^{(2)},...,\mu^{(k)})`. If

.. MATH::

    n = \lvert \mu \rvert = \lvert \mu^{(1)} \rvert +
    \lvert \mu^{(2)} \rvert + \cdots + \lvert \mu^{(k)} \rvert

then `\mu` is a `k`-partition of `n`.

In representation theory partition tuples arise as the natural indexing
set for the ordinary irreducible representations of:

- the wreath products of cyclic groups with symmetric groups
- the Ariki-Koike algebras, or the cyclotomic Hecke algebras of
  the complex reflection groups of type `G(r,1,n)`
- the degenerate cyclotomic Hecke algebras of type `G(r,1,n)`

When these algebras are not semisimple partition, tuples index an important
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

.. [DJM99] R. Dipper, G. James and A. Mathas "The cyclotomic q-Schur algebra", Math. Z,
      229 (1999), 385-416.
.. [BK09] J. Brundan and A. Kleshchev "Graded decomposition numbers for cyclotomic Hecke algebras",
      Adv. Math., 222 (2009), 1883-1942"

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
    sage: PartitionTuples(3,2)[:]
    [([2], [], []),
     ([1, 1], [], []),
     ([1], [1], []),
     ([1], [], [1]),
     ([], [2], []),
     ([], [1, 1], []),
     ([], [1], [1]),
     ([], [], [2]),
     ([], [], [1, 1])]

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

    sage: mu=PartitionTuple([[4,1],[],[2,2,1],[3]])
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
    [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 1, 0), (2, 0, 0), (2, 0, 1), (2, 1, 0), (2, 1, 1), (2, 2, 0), (3, 0, 0), (3, 0, 1), (3, 0, 2)]
    sage: mu.addable_cells()
    [(0, 0, 4), (0, 1, 1), (0, 2, 0), (1, 0, 0), (2, 0, 2), (2, 2, 1), (2, 3, 0), (3, 0, 3), (3, 1, 0)]
    sage: mu.removable_cells()
    [(0, 0, 3), (0, 1, 0), (2, 1, 1), (2, 2, 0), (3, 0, 2)]

Attached to a partition tuple is the corresponding Young, or parabolic,
subgroup::

    sage: mu.young_subgroup()
    Permutation Group with generators [(), (12,13), (11,12), (8,9), (6,7), (3,4), (2,3), (1,2)]
    sage: mu.young_subgroup_generators()
    [1, 2, 3, 6, 8, 11, 12]

"""
#*****************************************************************************
#       Copyright (C) 2012 Andrew Mathas <andrew.mathas@sydney.edu.au>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartesian_product import CartesianProduct
from combinat import CombinatorialObject
from integer_vector import IntegerVectors
from partition import Partition, Partitions, Partitions_n, _Partitions
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.interfaces.all import gp
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.all import NN, ZZ
from sage.rings.integer import Integer
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

#--------------------------------------------------
# Partition tuple - element class
#--------------------------------------------------
class PartitionTuple(CombinatorialObject,Element):
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
    __metaclass__ = ClasscallMetaclass
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
        if mu==[] or mu==[[]]:
            return Partition([])

        # We must check mu is a list of partitions
        try:
            mu=[Partition(mu)]
        except ValueError:
            try:
                mu=[Partition(nu) for nu in mu]
            except ValueError:
                raise ValueError('%s is not a tuple of Partitions' % mu)

        if len(mu)==1:
            return Partition(mu[0])
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
        Element.__init__(self, parent)
        mu=[Partition(nu) for nu in mu]
        CombinatorialObject.__init__(self, mu)

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

        EXAMPLE::

            sage: len( PartitionTuple([[2,1],[3,2],[1,1,1]]) )
            3

        """
        return self.level()

    def _repr_(self, compact=False):
        """
        Return a string representation of ``self`` depending on
        :meth:`PartitionTuples.global_options`.

        EXAMPLES::

            sage: mu=PartitionTuple(([2,1],[3,2],[1,1,1]))      # indirect doctest

            sage: PartitionTuples.global_options(display="list"); mu
            ([2, 1], [3, 2], [1, 1, 1])
            sage: PartitionTuples.global_options(display="diagram"); mu
            **   ***   *
            *    **    *
            *
            sage: PartitionTuples.global_options(display="compact_low"); mu
            1,2|2,3|1^3
            sage: PartitionTuples.global_options(display="compact_high"); mu
            2,1|3,2|1^3
            sage: PartitionTuples.global_options(display="exp_low"); mu
            1, 2 | 2, 3 | 1^3
            sage: PartitionTuples.global_options(display="exp_high"); mu
            2, 1 | 3, 2 | 1^3
            sage: PartitionTuples.global_options.reset()

            sage: Partitions.global_options(convention="French");
            sage: PartitionTuples.global_options(display="diagram"); mu
            *
            *    **    *
            **   ***   *
            sage: PartitionTuples.global_options(display="list"); mu
            ([2, 1], [3, 2], [1, 1, 1])
            sage: PartitionTuples.global_options(display="compact_low"); mu
            1,2|2,3|1^3
            sage: PartitionTuples.global_options(display="compact_high"); mu
            2,1|3,2|1^3
            sage: PartitionTuples.global_options(display="exp_low"); mu
            1, 2 | 2, 3 | 1^3
            sage: PartitionTuples.global_options(display="exp_high"); mu
            2, 1 | 3, 2 | 1^3
            sage: PartitionTuples.global_options.reset()
        """
        if compact:
            from sage.misc.superseded import deprecation
            deprecation(13605, 'compact option is deprecated. Use PartitionTuples.global_options instead.')
            return self._repr_compact_high()

        return self.parent().global_options.dispatch(self, '_repr_', 'display')

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as a Ferrers diagram.

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_diagram()
            '   **   ***   *\n   *    **    *\n              *'
        """
        return self.diagram()

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list.

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_list()
            '([2, 1], [3, 2], [1, 1, 1])'
        """
        return '('+', '.join(nu._repr_() for nu in self)+')'

    def _repr_exp_low(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_exp_low()
            '1, 2 | 2, 3 | 1^3'
        """
        return ' | '.join(nu._repr_exp_low() for nu in self)

    def _repr_exp_high(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1,1,1,1,1,1,1,1]))._repr_exp_high()
            '2, 1 | 3, 2 | 1^10'
        """
        return ' | '.join(nu._repr_exp_high() for nu in self)

    def _repr_compact_low(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_compact_low()
            '1,2|2,3|1^3'
        """
        return '%s' % '|'.join('-' if mu==[] else mu._repr_compact_low() for mu in self)

    def _repr_compact_high(self):
        """
        Return a string representation of ``self`` in compact form (exponential
        form with highest first).

        EXAMPLES::

            sage: PartitionTuple(([2,1],[3,2],[1,1,1]))._repr_compact_high()
            '2,1|3,2|1^3'
        """
        return '%s' % '|'.join('-' if mu==[] else mu._repr_compact_high() for mu in self)

    # override default string representation which is str(self._list)
    __str__=lambda self: self._repr_()


    def _latex_(self):
        r"""
        Returns a LaTeX version of ``self``.

        For more on the latex options, see :meth:`Partitions.option`.

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1]])
            sage: PartitionTuples.global_options(latex='diagram'); latex(mu)       # indirect doctest
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
            sage: PartitionTuples.global_options(latex='exp_high'); latex(mu)      # indirect doctest
            (2,1|1^{3})
            sage: PartitionTuples.global_options(latex='exp_low'); latex(mu)       # indirect doctest
            (1,2|1^{3})
            sage: PartitionTuples.global_options(latex='list'); latex(mu)          # indirect doctest
            [[2, 1], [1, 1, 1]]
            sage: PartitionTuples.global_options(latex='young_diagram'); latex(mu) # indirect doctest
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

            sage: PartitionTuples.global_options(latex="young_diagram", convention="french")
            sage: PartitionTuples.global_options(latex='exp_high'); latex(mu)      # indirect doctest
            (2,1|1^{3})
            sage: PartitionTuples.global_options(latex='exp_low'); latex(mu)       # indirect doctest
            (1,2|1^{3})
            sage: PartitionTuples.global_options(latex='list'); latex(mu)          # indirect doctest
            [[2, 1], [1, 1, 1]]
            sage: PartitionTuples.global_options(latex='young_diagram'); latex(mu) # indirect doctest
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

            sage: PartitionTuples.global_options.reset()
        """
        return self.parent().global_options.dispatch(self, '_latex_', 'latex')

    def _latex_young_diagram(self):
        """
        LaTeX output as a Young diagram.

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1]])._latex_young_diagram()
        """
        from sage.combinat.output import tex_from_array_tuple
        return tex_from_array_tuple([ [["\\phantom{x}"]*row for row in mu] for mu in self._list ])

    def _latex_diagram(self):
        """
        LaTeX output as a Ferrers' diagram.

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1]])._latex_diagram()
        """
        entry = self.parent().global_options("latex_diagram_str")
        from sage.combinat.output import tex_from_array_tuple
        return tex_from_array_tuple([ [[entry]*row for row in mu] for mu in self._list ], with_lines=False)

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
        return '(%s)' % '|'.join(','.join('%s%s'%(a+1,'' if e==1 else '^{%s}'%e)
                        for (a,e) in enumerate(mu)) for mu in self.to_exp())

    def _latex_exp_high(self):
        """
        LaTeX output in exponential notation (highest first).

        EXAMPLES::

            sage: mu = PartitionTuple([[2, 1],[1,1,1,1,1,1,1,1,1,1]])._latex_exp_high()
        """
        return '(%s)' % '|'.join(','.join(['%s%s'%(a+1,'' if e==1 else '^{%s}'%e)
            for (a,e) in enumerate(mu)][::-1]) for mu in self.to_exp())


    def components(self):
        r"""
        Return a list containing the shape of this partition.

        This function exists in order to give a uniform way of iterating over
        the \"components\" of partition tuples of level 1 (partitions) and for
        higher levels.

        EXAMPLE::

            sage: for t in PartitionTuple([[2,1],[3,2],[3]]).components():  print '%s\n' % t.ferrers_diagram()
            ...
            **
            *
            <BLANKLINE>
            ***
            **
            <BLANKLINE>
            ***
            <BLANKLINE>
            sage: for t in PartitionTuple([3,2]).components(): print '%s\n'% t.ferrers_diagram()
            ...
            ***
            **
        """
        return [ t for t in self ]


    def diagram(self):
        r"""
        Return a string for the Ferrers diagram of ``self``.

        EXAMPLES::

            sage: print PartitionTuple([[2,1],[3,2],[1,1,1]]).diagram()
               **   ***   *
               *    **    *
                          *
            sage: print PartitionTuple([[3,2],[2,1],[],[1,1,1,1]]).diagram()
               ***   **   -   *
               **    *        *
                              *
                              *
            sage: PartitionTuples.global_options(convention="french")
            sage: print PartitionTuple([[3,2],[2,1],[],[1,1,1,1]]).diagram()
                              *
                              *
               **    *        *
               ***   **   -   *
            sage: PartitionTuples.global_options.reset()
        """
        col_len = [len(mu)>0 and mu[0] or 1 for mu in self]  # columns per component
        row_max = max(len(mu) for mu in self)                # maximum row length
        # There should be a fancier list compression for this but I couldn't get
        # one to work in the cases where a component was the empty partition
        diag = []
        diag_str=PartitionTuples.global_options('diagram_str')
        for row in range(row_max):
            line=''
            for c in range(len(self)):
                if row == 0 and self[c] == []:
                    line += '   -'
                elif row < len(self[c]):
                    line += '   {:{}}'.format(diag_str*self[c][row],col_len[c])
                else:
                    line += '   {:{}}'.format('',col_len[c])
            diag.append(line.rstrip())
        if PartitionTuples.global_options('convention') == "English":
            return '\n'.join(map(str, diag))
        else:
            return '\n'.join(map(str, diag[::-1]))

    ferrers_diagram = diagram


    def pp(self):
        r"""
        Pretty prints this partition tuple. See :meth:`diagram`.

        EXAMPLES::

            sage: PartitionTuple([[5,5,2,1],[3,2]]).pp()
            *****   ***
            *****   **
            **
            *
        """
        print self.diagram()


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


    def standard_tableaux(self):
        """
        Return the :class:`standard tableau tuples<StandardTableauTuples>` of
        this shape.

        EXAMPLE::

            sage: PartitionTuple([[],[3,2,2,1],[2,2,1],[3]]).standard_tableaux()
            Standard tableau tuples of shape ([], [3, 2, 2, 1], [2, 2, 1], [3])
        """
        from tableau_tuple import StandardTableauTuples
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
                up=[tau for tau in self]
                up[c]=nu
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
        return [mu for mu in self.up()]

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
                down=[tau for tau in self]
                down[c]=nu
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
        return [mu for mu in self.down()]

    def cells(self):
        """
        Return the coordinates of the cells of ``self``. Coordinates are given
        as (component index, row index, column index) and are 0 based.

        EXAMPLES::

            sage: PartitionTuple([[2,1],[1],[1,1,1]]).cells()
            [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 2, 0)]
        """
        return [(c,a,b) for c in range(len(self)) for (a,b) in self[c].cells()]

    def content(self, k,r,c, multicharge):
        r"""
        Returns the content of the cell.

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
            mu=PartitionTuple(mu)
        except ValueError:
            raise ValueError('%s must be a PartitionTuple' % mu)

        if mu==self: return True
        level=0
        ssum=0  # sum of successive rows in self
        musum=0 # sum of successive rows in self
        while level<self.level() and level<mu.level():
            row=0
            while row<len(self[level]) and row<len(mu[level]):
                ssum+=self[level][row]
                musum+=mu[level][row]
                if musum>ssum: return False
                row+=1
            if row<len(self[level]):
                ssum+=sum(self[level][row:])
            elif row<len(mu[level]):
                musum+=sum(mu[level][row:])
                if musum>ssum: return False
            level+=1
        return True

    def initial_tableau(self):
        r"""
        Return the :class:`StandardTableauTuple` which has the numbers
        `1, 2, \ldots, n`, where `n` is the :meth:`size` of ``self``,
        entered in order from left to right along the rows of each component,
        where the components are ordered from left to right.

        EXAMPLE::

            sage: PartitionTuple([ [2,1],[3,2] ]).initial_tableau()
            ([[1, 2], [3]], [[4, 5, 6], [7, 8]])
        """
        from tableau_tuple import StandardTableauTuples
        return StandardTableauTuples(self).first()

    def garnir_tableau(self, *cell):
        r"""
        Return the Garnir tableau of shape ``self`` corresponding to the cell
        ``cell``.

        If ``cell`` `= (k,a,c)` then `(k,a+1,c)` must belong to the diagram of
        the :class:`PartitionTuple`. If this is not the case when we return
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
            (comp, row,col)=cell
        except ValueError:
            (comp, row,col)=cell[0]

        if comp>=len(self) or row+1>=len(self[comp]) or col>=self[comp][row+1]:
            raise ValueError('(comp, row+1, col) must be inside the diagram')
        from tableau_tuple import TableauTuple
        g = TableauTuple(self.initial_tableau().to_list())
        a=g[comp][row][col]
        g[comp][row][col:]=range(a+col+1,g[comp][row+1][col]+1)
        g[comp][row+1][:col+1]=range(a,a+col+1)
        g._garnir_cell=(comp,row,col)
        return g

    def top_garnir_tableau(self,e,cell):
        r"""
        Return the most dominant *standard* tableau which dominates the
        corresponding Garnir tableau and has the same residue that has shape
        ``self`` and is determined by ``e`` and ``cell``.

        The Garnir tableau play an important role in integral and
        non-semisimple representation theory because they determine the
        "straightening" rules for the Specht modules over an arbitrary ring.
        The *top Garnir tableaux* arise in the graded representation theory of
        the symmetric groups and higher level Hecke algebras. They were
        introduced in [KMR]_.

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

        REFERENCE:

        .. [KMR] A. Kleshchev, A. Mathas, and A. Ram, *Universal Specht
           modules for cyclotomic Hecke algebras*,
           Proc. London Math. Soc. (2012) 105 (6): 1245-1289.
           :arxiv:`1102.3519v1`

        """
        (comp,row,col)=cell
        if comp>=len(self) or row+1>=len(self[comp]) or col>=self[comp][row+1]:
            raise ValueError('(comp, row+1, col) must be inside the diagram')

        g=self.garnir_tableau(cell)

        if e==0: return      # no more dominant tableau of the same residue

        a=e*int((self[comp][row]-col)/e)    # number of cells in the e-bricks in row `row`
        b=e*int((col+1)/e)            # number of cells in the e-bricks in row `row+1`

        if a==0 or b==0: return self.garnir_tableau(cell)

        t=g.to_list()
        m=t[comp][row+1][0]            # smallest number of 0-Garnir belt
        # now we will put the number m,m+1,...,t[row+1][col] in order into t
        t[comp][row][col:a+col]=[m+col-b+1+i for i in range(a)]
        t[comp][row+1][col-b+1:col+1]=[m+a+col-b+1+i for i in range(b)]
        from tableau_tuple import StandardTableauTuple
        return StandardTableauTuple(t)

    def arm_length(self, k,r,c):
        """
        Return the length of the arm of cell ``(k, r, c)`` in ``self``.

        INPUT:

        - ``k`` -- The compoenent
        - ``r`` -- The row
        - ``c`` -- The cell

        OUTPUT:

        - The arm length as an integer

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
            raise ValueError("The cell %s is not in the diagram" %((k,r,c),))

    def leg_length(self, k,r,c):
        """
        Return the length of the leg of cell ``(k, r, c)`` in ``self``.

        INPUT:

        - ``k`` -- The compoenent
        - ``r`` -- The row
        - ``c`` -- The cell

        OUTPUT:

        - The leg length as an integer

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
        Returns ``True`` if this partition tuple contains `\mu`.

        If `\lambda=(\lambda^{(1)}, \ldots, \lambda^{(l)})` and
        `\mu=(\mu^{(1)}, \ldots, \mu^{(m)})` are two partition tuples then
        `\lambda` contains `\mu` if `m \leq l` and
        `\mu^{(i)}_r \leq \lambda^{(i)}_r` for `1 \leq i \leq m` and `r \geq 0`.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[2],[2,1]]).contains( PartitionTuple([[1,1],[2],[2,1]]) )
            True
        """
        return mu.level()<=self.level() and all(self[c].contains(mu[c]) for c in range(len(mu)))

    def hook_length(self, k,r,c):
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
        Returns a list of the removable cells of this partition tuple.

        All indice are of the form ``(k, r, c)``, where ``r`` is the
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

    def add_cell(self, k,r,c):
        r"""
        Return the partition tuple obtained by adding a cell in row ``r``,
        column ``c``, and component ``k``.

        This does not change ``self``.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).add_cell(0,0,1)
            ([2, 1], [4, 3], [2, 1, 1])

        """
        if (k,r,c) in self.addable_cells(): # an addable cell
            mu=self.to_list()
            if c==0: mu[k].append(1)
            else: mu[k][r]+=1
            return PartitionTuple(mu)
        else:
            raise ValueError("%s is not an addable cell"%((k,r,c),))

    def remove_cell(self, k,r,c):
        """
        Return the partition tuple obtained by removing a cell in row ``r``,
        column ``c``, and component ``k``.

        This does not change ``self``.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).remove_cell(0,1,0)
            ([1], [4, 3], [2, 1, 1])

        """
        if (k,r,c) in self.removable_cells():  # a removable cell
            mu=self.to_list()
            mu[k][r]-=1
            return PartitionTuple(mu)
        else:
            raise ValueError("%s is not a removable cell"%((k,r,c),))

    def to_list(self):
        r"""
        Return ``self`` as a list of lists.

        EXAMPLES::

            sage: PartitionTuple([[1,1],[4,3],[2,1,1]]).to_list()
            [[1, 1], [4, 3], [2, 1, 1]]

        TESTS::

            sage: all(mu==PartitionTuple(mu.to_list()) for mu in PartitionTuples(4,4))
            True
        """
        return [ mu.to_list() for mu in self]

    def young_subgroup(self):
        """
        Return the corresponding Young, or parabolic, subgroup of the
        symmetric group.

        EXAMPLE::

            sage: PartitionTuple([[2,1],[4,2],[1]]).young_subgroup()
            Permutation Group with generators [(), (8,9), (6,7), (5,6), (4,5), (1,2)]
        """
        gens=[]
        m=0
        for comp in self:
            for row in comp:
                gens.extend([(c,c+1) for c in range(m+1,m+row)])
                m+=row
        gens.append( range(1,self.size()+1) )  # to ensure we get a subgroup of Sym_n
        return PermutationGroup( gens )

    def young_subgroup_generators(self):
        """
        Return an indexing set for the generators of the corresponding Young
        subgroup.

        EXAMPLE::

            sage: PartitionTuple([[2,1],[4,2],[1]]).young_subgroup_generators()
            [1, 4, 5, 6, 8]
        """
        gens=[]
        m=0
        for comp in self:
            for row in comp:
                gens.extend([c for c in range(m+1,m+row)])
                m+=row
        return gens

#--------------------------------------------------
# Partition tuples - parent classes
#--------------------------------------------------
class PartitionTuples(UniqueRepresentation, Parent):
    """
    List of all partition tuples.

    For more information about partition tuples, see :class:`PartitionTuple`.

    This is a factory class which returns the appropriate parent based on
    the values of ``level`` and ``size``.

    INPUT:

    - ``level`` -- The length of the tuple

    - ``size``  -- The total number of cells

    TESTS::

        sage: TestSuite( PartitionTuples() ).run()
        sage: TestSuite( PartitionTuples(level=4) ).run()
        sage: TestSuite( PartitionTuples(size=6) ).run()
        sage: TestSuite( PartitionTuples(level=4, size=5) ).run()

        sage: [ [2,1],[],[3] ] in PartitionTuples()
        True
        sage: ( [2,1],[],[3] ) in PartitionTuples()
        True
        sage: ( [] ) in PartitionTuples()
        True

    Check that trac:`14145` has been fixed::

        sage: 1 in PartitionTuples()
        False
    """

    @staticmethod
    def __classcall_private__(klass, level=None, size=None):
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
        """

        # sanity testing
        if level is not None and (not isinstance(level,(int,Integer)) or level<1):
            raise ValueError('the level must be a positive integer')

        if size is not None and (not isinstance(size,(int,Integer)) or size<0):
            raise ValueError('the size must be a non-negative integer')

        if level==1:
            if size is not None: return Partitions_n(size)
            else: return Partitions()
        elif level is not None and size is not None:
            return PartitionTuples_level_size(level=level, size=size)
        elif level is not None:
            return PartitionTuples_level(level=level)
        elif size is not None:
            return PartitionTuples_size(size=size)
        else:
            return PartitionTuples_all()

    Element = PartitionTuple
    global_options=Partitions.global_options

    # default for level
    _level=None
    _size=None

    def _element_constructor_(self, mu):
        r"""
        Constructs an element of :class:`PartitionTuple`.

        INPUT:

        - ``mu`` -- A tuple of partitions

        OUTPUT:

        - The corresponding :class:`PartitionTuple` object

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
        if mu==[] or mu==[[]]:
            return Partition([])

        # As partitions are 1-tuples of partitions we need to treat them separately
        try:
            mu=[Partition(mu)]
        except ValueError:
            try:
                mu=[Partition(nu) for nu in mu]
            except ValueError:
                raise ValueError('%s is not a %s'%(mu, self))

        if mu in self:
            if len(mu)==1:
                return Partition(mu[0])
            else:
                return self.element_class(self, mu)

    def __contains__(self,mu):
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
            sage: all(mu in PartitionTuples() for mu in PartitionTuples(3,8))
            True

        Check that :trac:`14145` is fixed::

            sage: 1 in PartitionTuples()
            False
        """
        if isinstance(mu, PartitionTuple) or isinstance(mu, Partition):
            return True
        if isinstance(mu, (tuple, list)):
            return all(m in _Partitions for m in mu)
        return False

    def __getitem__(self, r):
        r"""
        The default implementation of ``__getitem__()`` for enumerated sets
        does not allow slices, so we override it.

        EXAMPLES::

            sage: PartitionTuples()[10:20]
            [[1, 1, 1],
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
            start=0 if r.start is None else r.start
            stop=r.stop
            if stop is None and not self.is_finite():
                raise ValueError('infinite set')
        else:
            raise ValueError('r must be an integer or a slice')
        count=0
        parts=[]
        for t in self:
            if count==stop:
                break
            if count>=start:
                parts.append(t)
            count+=1

        # this is to cope with empty slice endpoints like [6:] or [:]
        if count==stop or stop is None:
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
            sage: PartitionTuples().an_element()  # indirect doctest
            ([1, 1, 1, 1], [2, 1, 1], [3, 1], [4])

        """
        return PartitionTuple( ([1,1,1,1],[2,1,1],[3,1],[4]) )


class PartitionTuples_all(PartitionTuples):
    """
    Class of partition tuples of a arbitrary level and arbitrary sum.
    """

    def __init__(self):
        r"""
        Initializes the class.

        EXAMPLE::

            sage: PartitionTuples()
            Partition tuples

        """
        super(PartitionTuples_all, self).__init__(category=InfiniteEnumeratedSets())


    def _repr_(self):
        r"""
        Return a string represention of ``self``.

        EXAMPLES::

            sage: PartitionTuples() # indirect doctest
            Partition tuples
        """
        return 'Partition tuples'

    def __iter__(self):
        r"""
        Iterate through the infinite class of partition tuples of arbitrary
        level and size.

        EXAMPLES::

            sage: PartitionTuples()[:20]
            [[],
             [1],
             ([], []),
             [2],
             [1, 1],
             ([1], []),
             ([], [1]),
             ([], [], []),
             [3],
             [2, 1],
             [1, 1, 1],
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
            for level in xrange(size+1):
                for mu in PartitionTuples_level_size(level+1,size-level):
                    yield self._element_constructor_(mu)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::
            sage: PartitionTuples().an_element()   # indirect doctest
            ([1, 1, 1, 1], [2, 1, 1], [3, 1], [4])

        """
        return self.element_class(self,([1,1,1,1],[2,1,1],[3,1],[4]))


class PartitionTuples_level(PartitionTuples):
    """
    Class of partition tuples of a fixed level, but summing to an arbitrary
    integer.
    """

    def __init__(self, level):
        r"""
        Initializes this class.

        EXAMPLES::

            sage: PartitionTuples(4)
            Partition tuples of level 4
            sage: PartitionTuples(level=6)
            Partition tuples of level 6
        """
        if not level in NN:
            raise ValueError('level must be a non-negative integer')
        super(PartitionTuples_level, self).__init__(category=InfiniteEnumeratedSets())
        self._level=level

    def _repr_(self):
        """
        Return a string represention of ``self``.

        EXAMPLES::

            sage: PartitionTuples(2)    # indirect doctest
            Partition tuples of level 2
        """
        return 'Partition tuples of level %s' % self.level()

    def __contains__(self,mu):
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
            sage: all(mu in PartitionTuples(3) for mu in PartitionTuples(3,8))
            True

        Check that :trac:`14145` is fixed::

            sage: 1 in PartitionTuples(level=2)
            False
        """
        return PartitionTuples.__contains__(self, mu) and len(mu) == self.level()

    def __iter__(self):
        r"""
        Iterate through the infinite class of partition tuples of fixed level.

        EXAMPLES::

        sage: parts=PartitionTuples(3)
        sage: [parts[k] for k in range(20)]
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
            for mu in PartitionTuples_level_size(self.level(),size):
                yield self._element_constructor_(mu)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::
            sage: PartitionTuples(level=4).an_element()  # indirect doctest
            ([], [1], [2], [3])

        """
        return self.element_class(self, tuple([l] for l in range(self.level()) ))


class PartitionTuples_size(PartitionTuples):
    """
    Class of partition tuples of a fixed size, but arbitrary level.
    """
    def __init__(self, size):
        r"""
        Initializes this class.

        EXAMPLES::
            sage: PartitionTuples(size=4)
            Partition tuples of size 4
            sage: PartitionTuples(size=6)
            Partition tuples of size 6
        """
        if not size in NN:
            raise ValueError('size must be a non-negative integer')
        super(PartitionTuples_size, self).__init__(category=InfiniteEnumeratedSets())
        self._size=size

    def _repr_(self):
        """
        Return a string represention of ``self``.

        EXAMPLES::

            sage: PartitionTuples(size=4)    # indirect doctest
            Partition tuples of size 4
        """
        return 'Partition tuples of size %d' % self.size()

    def __contains__(self,mu):
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
            sage: all(mu in PartitionTuples(size=8) for mu in PartitionTuples(3,8))
            True

        Check that :trac:`14145` is fixed::

            sage: 1 in PartitionTuples(size=7)
            False
        """
        return (mu in _Partitions and self.size() == sum(mu)) or \
            (PartitionTuples.__contains__(self, mu) and self.size() == sum(map(sum, mu)))

    def __iter__(self):
        r"""
        Iterates through the infinite class of partition tuples of a fixed size.

        EXAMPLES::

            sage: PartitionTuples(size=3)[:20]
            [[3],
             [2, 1],
             [1, 1, 1],
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
            for mu in PartitionTuples_level_size(level,self.size()):
                yield self._element_constructor_(mu)

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::
            sage: PartitionTuples(size=4).an_element()  # indirect doctest
            ([1], [1], [1], [1])

        """
        return self.element_class(self, tuple([1] for l in range(self.size()) ))


class PartitionTuples_level_size(PartitionTuples):
    """
    Class of partition tuples with a fixed level and a fixed size.
    """

    def __init__(self, level, size):
        r"""
        Initializes this class.

        EXAMPLES::

            sage: PartitionTuples(4,2)
            Partition tuples of level 4 and size 2
        """
        if not (level in NN and size in NN):
            raise ValueError('n and level must be non-negative integers (or None)')
        super(PartitionTuples_level_size, self).__init__(category=FiniteEnumeratedSets())
        self._level=level
        self._size=size

    def _repr_(self):
        """
        Return a string represention of ``self``.

        EXAMPLES::

            sage: PartitionTuples(4,2)    # indirect doctest
            Partition tuples of level 4 and size 2
            sage: PartitionTuples(size=2,level=4)
            Partition tuples of level 4 and size 2
        """
        return 'Partition tuples of level %s and size %s' % (self.level(), self.size())

    def __contains__(self,mu):
        r"""
        Return ``True`` if `\mu` is in ``self``.

        TESTS::

            sage: PartitionTuple([[3,2],[2]]) in PartitionTuples(2,7)
            True
            sage: PartitionTuple([[3,2],[],[],[],[2]]) in PartitionTuples(5,7)
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[2]]) in PartitionTuples(5,7)
            True
            sage: PartitionTuple([[2,1],[],[1,1],[],[3]]) in PartitionTuples(2,8)
            False
            sage: all(mu in PartitionTuples(3,8) for mu in PartitionTuples(3,8))
            True

        Check that :trac:`14145` is fixed::

            sage: 1 in PartitionTuples(5,7)
            False
        """
        return (   (mu in Partitions() and self.size()==sum(mu) and self.level()==1)
                or (mu in PartitionTuples() and self.size()==sum(map(sum,mu))) and self.level()==len(mu))

    def __iter__(self):
        r"""
        Iterates through the finite class of partition tuples of a fixed level
        and a fixed size.

        EXAMPLES::

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
        p = [Partitions(i) for i in range(self.size()+1)]
        for iv in IntegerVectors(self.size(),self.level()):
            for cp in CartesianProduct(*[p[i] for i in iv]):
                yield self._element_constructor_(cp)


    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: PartitionTuples(level=4,size=4).an_element()   # indirect doctest
            ([1], [], [], [3])

        """
        mu=[[] for l in range(self.level())]
        if self.size()>0:
            if self.level()==1: mu=[self.size()-1,1]
            else:
                mu[0]=[1]
                mu[-1]=[self.size()-1]
        return self.element_class(self, mu)

    def cardinality(self):
        r"""
        Returns the number of ``level``-tuples of partitions of size ``n``.

        Wraps a pari function call.

        EXAMPLES::

            sage: PartitionTuples(2,3).cardinality()
            10
            sage: PartitionTuples(2,8).cardinality()
            185

        TESTS:

        The following calls used to fail (:trac:`11476`)::

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
        return ZZ(gp.eval('polcoeff((1/eta(x+O(x^%s)))^%s, %s, x)'%(self.size()+1,self.level(), self.size())))

    def __setstate__(self, state):
        """
        In order to maintain backwards compatibility and be able to unpickle a
        old pickle from PartitionTuples_nk we have to override the default
        __setstate__.

        TESTS::

            sage: loads("x\x9cM\x90\xcdN\xc30\x0c\x80\xd5\xc1\x06\xeb\x80\xf1{\xe0\r\xe0\xd2\x0b\x07\x1e\x02)B\x88\x9c-7\xb5\xba\xa8MR')\x12\x07$8p\xe0\xadq\x996q\xb1b\xfb\xb3\xf59\x9f3\x93\xb0\xa5\xca\x04W[\x8f\xb9\x1a0f\x9bm\xf0\xe5\xf3\xee\xf5:\x0e=%\xf0]\xc9\xc5\xfd\x17\xcf>\xf8\xe0N_\x83\xf5\xd2\xc5\x1e\xd0L\x10\xf46e>T\xba\x04r55\x8d\xf5-\xcf\x95p&\xf87\x8a\x19\x1c\xe5Mh\xc0\xa3#^(\xbd\x00\xd3`F>Rz\t\x063\xb5!\xbe\xf3\xf1\xd4\x98\x90\xc4K\xa5\x0b\xbf\xb5\x8b\xb2,U\xd6\x0bD\xb1t\xd8\x11\xec\x12.u\xf1\xf0\xfd\xc2+\xbd\x82\x96<E\xcc!&>Qz\x0e5&\xe2S\xa5\xd70X\xd3\xf5\x04\xe2\x91\xc4\x95\xcf\x9e\n\x11\xa3\x9e\x1c\xf9<\t\xa6\x1cG#\x83\xbcV\xfaf\x7f\xd9\xce\xfc\xef\xb4s\xa5o\xf7#\x13\x01\x03\xa6$!J\x81/~t\xd1m\xc4\xe5Q\\.\xff\xfd\x8e\t\x14\rmW\\\xa9\xb1\xae~\x01/\x8f\x85\x02")
            Partition tuples of level 7 and size 3
            sage: loads(dumps( PartitionTuples(7,3) ))  # indirect doctest for unpickling a Tableau element
            Partition tuples of level 7 and size 3
        """
        if isinstance(state, dict):   # for old pickles from Tableau_class
            parts=PartitionTuples(state['k'], state['n'])
            self.__class__=parts.__class__
            self.__dict__=parts.__dict__
        else:
            super(PartitionTuples, self).__setstate__(state)
