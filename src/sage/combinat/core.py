"""
Cores

A `k`-core is a partition from which no rim hook of size `k` can be removed.
Alternatively, a `k`-core is an integer partition such that the Ferrers
diagram for the partition contains no cells with a hook of size (a
multiple of) `k`.

Authors:

- Anne Schilling and Mike Zabrocki (2011): initial version
- Travis Scrimshaw (2012): Added latex output for Core class
"""
# ****************************************************************************
#       Copyright (C) 2011 Anne Schilling <anne at math.ucdavis.edu>
#                          Mike Zabrocki  <zabrocki at mathstat.yorku.ca>
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
# ***************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.combinat.partition import Partitions, Partition
from sage.combinat.combinat import CombinatorialElement
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.combinatorial_map import combinatorial_map


class Core(CombinatorialElement):
    r"""
    A `k`-core is an integer partition from which no rim hook of size `k`
    can be removed.

    EXAMPLES::

        sage: c = Core([2,1],4); c
        [2, 1]
        sage: c = Core([3,1],4); c
        Traceback (most recent call last):
        ...
        ValueError: [3, 1] is not a 4-core
    """
    @staticmethod
    def __classcall_private__(cls, part, k):
        r"""
        Implement the shortcut ``Core(part, k)`` to ``Cores(k,l)(part)``
        where `l` is the length of the core.

        TESTS::

            sage: c = Core([2,1],4); c
            [2, 1]
            sage: c.parent()
            4-Cores of length 3
            sage: type(c)
            <class 'sage.combinat.core.Cores_length_with_category.element_class'>

            sage: Core([2,1],3)
            Traceback (most recent call last):
            ...
            ValueError: [2, 1] is not a 3-core
        """
        if isinstance(part, cls):
            return part
        part = Partition(part)
        if not part.is_core(k):
            raise ValueError("%s is not a %s-core" % (part, k))
        l = sum(part.k_boundary(k).row_lengths())
        return Cores(k, l)(part)

    def __init__(self, parent, core):
        """
        TESTS::

            sage: C = Cores(4,3)
            sage: c = C([2,1]); c
            [2, 1]
            sage: type(c)
            <class 'sage.combinat.core.Cores_length_with_category.element_class'>
            sage: c.parent()
            4-Cores of length 3
            sage: TestSuite(c).run()

            sage: C = Cores(3,3)
            sage: C([2,1])
            Traceback (most recent call last):
            ...
            ValueError: [2, 1] is not a 3-core
        """
        k = parent.k
        part = Partition(core)
        if not part.is_core(k):
            raise ValueError("%s is not a %s-core" % (part, k))
        CombinatorialElement.__init__(self, parent, core)

    def __eq__(self, other):
        """
        Test for equality.

        EXAMPLES::

            sage: c = Core([4,2,1,1],5)
            sage: d = Core([4,2,1,1],5)
            sage: e = Core([4,2,1,1],6)
            sage: c == [4,2,1,1]
            False
            sage: c == d
            True
            sage: c == e
            False
        """
        if isinstance(other, Core):
            return (self._list == other._list and
                    self.parent().k == other.parent().k)
        return False

    def __ne__(self, other):
        """
        Test for un-equality.

        EXAMPLES::

            sage: c = Core([4,2,1,1],5)
            sage: d = Core([4,2,1,1],5)
            sage: e = Core([4,2,1,1],6)
            sage: c != [4,2,1,1]
            True
            sage: c != d
            False
            sage: c != e
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Compute the hash of ``self`` by computing the hash of the
        underlying list and of the additional parameter.

        The hash is cached and stored in ``self._hash``.

        EXAMPLES::

            sage: c = Core([4,2,1,1],3)
            sage: c._hash is None
            True
            sage: hash(c) #random
            1335416675971793195
            sage: c._hash #random
            1335416675971793195

        TESTS::

            sage: c = Core([4,2,1,1],5)
            sage: d = Core([4,2,1,1],6)
            sage: hash(c) == hash(d)
            False
        """
        if self._hash is None:
            self._hash = hash(tuple(self._list)) + hash(self.parent().k)
        return self._hash

    def _latex_(self):
        r"""
        Output the LaTeX representation of this core as a partition.

        See the ``_latex_`` method of :class:`Partition`.

        EXAMPLES::

            sage: c = Core([2,1],4)
            sage: latex(c)
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-2}
            \lr{\phantom{x}}\\\cline{1-1}
            \end{array}$}
            }
        """
        return self.to_partition()._latex_()

    def k(self):
        r"""
        Return `k` of the `k`-core ``self``.

        EXAMPLES::

            sage: c = Core([2,1],4)
            sage: c.k()
            4
        """
        return self.parent().k

    @combinatorial_map(name="to partition")
    def to_partition(self):
        r"""
        Turn the core ``self`` into the partition identical to ``self``.

        EXAMPLES::

            sage: mu = Core([2,1,1],3)
            sage: mu.to_partition()
            [2, 1, 1]
        """
        return Partition(self)

    @combinatorial_map(name="to bounded partition")
    def to_bounded_partition(self):
        r"""
        Bijection between `k`-cores and `(k-1)`-bounded partitions.

        This maps the `k`-core ``self`` to the corresponding `(k-1)`-bounded partition.
        This bijection is achieved by deleting all cells in ``self`` of hook length
        greater than `k`.

        EXAMPLES::

            sage: gamma = Core([9,5,3,2,1,1], 5)
            sage: gamma.to_bounded_partition()
            [4, 3, 2, 2, 1, 1]
        """
        k_boundary = self.to_partition().k_boundary(self.k())
        return Partition(k_boundary.row_lengths())

    def size(self):
        r"""
        Return the size of ``self`` as a partition.

        EXAMPLES::

            sage: Core([2,1],4).size()
            3
            sage: Core([4,2],3).size()
            6
        """
        return self.to_partition().size()

    def length(self):
        r"""
        Return the length of ``self``.

        The length of a `k`-core is the size of the corresponding `(k-1)`-bounded partition
        which agrees with the length of the corresponding Grassmannian element,
        see :meth:`to_grassmannian`.

        EXAMPLES::

            sage: c = Core([4,2],3); c.length()
            4
            sage: c.to_grassmannian().length()                                          # needs sage.modules
            4

            sage: Core([9,5,3,2,1,1], 5).length()
            13
        """
        return self.to_bounded_partition().size()

    def to_grassmannian(self):
        r"""
        Bijection between `k`-cores and Grassmannian elements in the affine Weyl group of type `A_{k-1}^{(1)}`.

        For further details, see the documentation of the method
        :meth:`~sage.combinat.partition.Partition.from_kbounded_to_reduced_word` and
        :meth:`~sage.combinat.partition.Partition.from_kbounded_to_grassmannian`.

        EXAMPLES::

            sage: c = Core([3,1,1],3)
            sage: w = c.to_grassmannian(); w                                            # needs sage.modules
            [-1  1  1]
            [-2  2  1]
            [-2  1  2]
            sage: c.parent()
            3-Cores of length 4
            sage: w.parent()                                                            # needs sage.modules
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root space)

            sage: c = Core([],3)
            sage: c.to_grassmannian()                                                   # needs sage.modules
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        bp = self.to_bounded_partition()
        return bp.from_kbounded_to_grassmannian(self.k() - 1)

    def affine_symmetric_group_simple_action(self, i):
        r"""
        Return the action of the simple transposition `s_i` of the affine symmetric group on ``self``.

        This gives the action of the affine symmetric group of type `A_k^{(1)}` on the `k`-core
        ``self``. If ``self`` has outside (resp. inside) corners of content `i` modulo `k`, then
        these corners are added (resp. removed). Otherwise the action is trivial.

        EXAMPLES::

            sage: c = Core([4,2],3)
            sage: c.affine_symmetric_group_simple_action(0)                             # needs sage.modules
            [3, 1]
            sage: c.affine_symmetric_group_simple_action(1)                             # needs sage.modules
            [5, 3, 1]
            sage: c.affine_symmetric_group_simple_action(2)                             # needs sage.modules
            [4, 2]

        This action corresponds to the left action by the `i`-th simple reflection in the affine
        symmetric group::

            sage: c = Core([4,2],3)
            sage: W = c.to_grassmannian().parent()                                      # needs sage.modules
            sage: i = 0
            sage: (c.affine_symmetric_group_simple_action(i).to_grassmannian()          # needs sage.modules
            ....:     == W.simple_reflection(i)*c.to_grassmannian())
            True
            sage: i = 1
            sage: (c.affine_symmetric_group_simple_action(i).to_grassmannian()          # needs sage.modules
            ....:     == W.simple_reflection(i)*c.to_grassmannian())
            True
        """
        mu = self.to_partition()
        corners = [p for p in mu.outside_corners()
                   if mu.content(p[0], p[1]) % self.k() == i]
        if not corners:
            corners = [p for p in mu.corners()
                       if mu.content(p[0], p[1]) % self.k() == i]
            if not corners:
                return self
            for p in corners:
                mu = mu.remove_cell(p[0])
        else:
            for p in corners:
                mu = mu.add_cell(p[0])
        return Core(mu, self.k())

    def affine_symmetric_group_action(self, w, transposition=False):
        r"""
        Return the (left) action of the affine symmetric group on ``self``.

        INPUT:

        - ``w`` -- tuple of integers `[w_1,\ldots,w_m]` with `0\le w_j<k`.
          If transposition is set to be ``True``, then `w = [w_0,w_1]` is
          interpreted as a transposition `t_{w_0, w_1}`
          (see :meth:`_transposition_to_reduced_word`).

        - ``transposition`` -- boolean (default: ``False``)

        The output is the (left) action of the product of the corresponding simple transpositions
        on ``self``, that is `s_{w_1} \cdots s_{w_m}(self)`. See :meth:`affine_symmetric_group_simple_action`.

        EXAMPLES::

            sage: c = Core([4,2],3)
            sage: c.affine_symmetric_group_action([0,1,0,2,1])
            [8, 6, 4, 2]
            sage: c.affine_symmetric_group_action([0,2], transposition=True)
            [4, 2, 1, 1]

            sage: c = Core([11,8,5,5,3,3,1,1,1],4)
            sage: c.affine_symmetric_group_action([2,5],transposition=True)
            [11, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        c = self
        if transposition:
            w = self._transposition_to_reduced_word(w)
        w.reverse()
        for i in w:
            c = c.affine_symmetric_group_simple_action(i)
        return c

    def _transposition_to_reduced_word(self, t):
        r"""
        Convert the transposition `t = [r,s]` to a reduced word.

        INPUT:

        - ``t`` -- a tuple `[r,s]` such that `r` and `s` are not equivalent mod `k`

        OUTPUT: list of integers in `\{0,1,\ldots,k-1\}` representing a
        reduced word for the transposition `t`

        EXAMPLES::

            sage: c = Core([],4)
            sage: c._transposition_to_reduced_word([2, 5])
            [2, 3, 0, 3, 2]
            sage: c._transposition_to_reduced_word([2, 5]) == c._transposition_to_reduced_word([5,2])
            True
            sage: c._transposition_to_reduced_word([2, 2])
            Traceback (most recent call last):
            ...
            ValueError: t_0 and t_1 cannot be equal mod k

            sage: c = Core([],30)
            sage: c._transposition_to_reduced_word([4, 12])
            [4, 5, 6, 7, 8, 9, 10, 11, 10, 9, 8, 7, 6, 5, 4]

            sage: c = Core([],3)
            sage: c._transposition_to_reduced_word([4, 12])
            [1, 2, 0, 1, 2, 0, 2, 1, 0, 2, 1]
        """
        k = self.k()
        if (t[0] - t[1]) % k == 0:
            raise ValueError("t_0 and t_1 cannot be equal mod k")
        if t[0] > t[1]:
            return self._transposition_to_reduced_word([t[1], t[0]])
        else:
            resu = [i % k for i in range(t[0], t[1] - (t[1] - t[0]) // k)]
            resu += [(t[1] - (t[1] - t[0]) // k - 2 - i) % k
                     for i in range(t[1] - (t[1] - t[0]) // k - t[0] - 1)]
            return resu

    def weak_le(self, other):
        r"""
        Weak order comparison on cores.

        INPUT:

        - ``other`` -- another `k`-core

        OUTPUT: boolean

        This returns whether ``self`` <= ``other`` in weak order.

        EXAMPLES::

            sage: c = Core([4,2],3)
            sage: x = Core([5,3,1],3)
            sage: c.weak_le(x)                                                          # needs sage.modules
            True
            sage: c.weak_le([5,3,1])                                                    # needs sage.modules
            True

            sage: x = Core([4,2,2,1,1],3)
            sage: c.weak_le(x)                                                          # needs sage.modules
            False

            sage: x = Core([5,3,1],6)
            sage: c.weak_le(x)
            Traceback (most recent call last):
            ...
            ValueError: the two cores do not have the same k
        """
        if type(self) is type(other):
            if self.k() != other.k():
                raise ValueError("the two cores do not have the same k")
        else:
            other = Core(other, self.k())
        w = self.to_grassmannian()
        v = other.to_grassmannian()
        return w.weak_le(v, side='left')

    def weak_covers(self):
        r"""
        Return a list of all elements that cover ``self`` in weak order.

        EXAMPLES::

            sage: c = Core([1],3)
            sage: c.weak_covers()                                                       # needs sage.modules
            [[1, 1], [2]]

            sage: c = Core([4,2],3)
            sage: c.weak_covers()                                                       # needs sage.modules
            [[5, 3, 1]]
        """
        w = self.to_grassmannian()
        S = w.upper_covers(side='left')
        S = (x for x in S if x.is_affine_grassmannian())
        return [x.affine_grassmannian_to_core() for x in set(S)]

    def strong_le(self, other):
        r"""
        Strong order (Bruhat) comparison on cores.

        INPUT:

        - ``other`` -- another `k`-core

        OUTPUT: boolean

        This returns whether ``self`` <= ``other`` in Bruhat (or strong) order.

        EXAMPLES::

            sage: c = Core([4,2],3)
            sage: x = Core([4,2,2,1,1],3)
            sage: c.strong_le(x)
            True
            sage: c.strong_le([4,2,2,1,1])
            True

            sage: x = Core([4,1],4)
            sage: c.strong_le(x)
            Traceback (most recent call last):
            ...
            ValueError: the two cores do not have the same k
        """
        if type(self) is type(other):
            if self.k() != other.k():
                raise ValueError("the two cores do not have the same k")
        else:
            other = Core(other, self.k())
        return other.contains(self)

    def contains(self, other):
        r"""
        Check whether ``self`` contains ``other``.

        INPUT:

        - ``other`` -- another `k`-core or a list

        OUTPUT: boolean

        This returns ``True`` if the Ferrers diagram of ``self`` contains the
        Ferrers diagram of ``other``.

        EXAMPLES::

            sage: c = Core([4,2],3)
            sage: x = Core([4,2,2,1,1],3)
            sage: x.contains(c)
            True
            sage: c.contains(x)
            False
        """
        la = self.to_partition()
        mu = Core(other, self.k()).to_partition()
        return la.contains(mu)

    def strong_covers(self):
        r"""
        Return a list of all elements that cover ``self`` in strong order.

        EXAMPLES::

            sage: c = Core([1],3)
            sage: c.strong_covers()
            [[2], [1, 1]]
            sage: c = Core([4,2],3)
            sage: c.strong_covers()
            [[5, 3, 1], [4, 2, 1, 1]]
        """
        S = Cores(self.k(), length=self.length() + 1)
        return [ga for ga in S if ga.contains(self)]

    def strong_down_list(self):
        r"""
        Return a list of all elements that are covered by ``self`` in strong order.

        EXAMPLES::

            sage: c = Core([1],3)
            sage: c.strong_down_list()
            [[]]
            sage: c = Core([5,3,1],3)
            sage: c.strong_down_list()
            [[4, 2], [3, 1, 1]]
        """
        if not self:
            return []
        return [ga for ga in Cores(self.k(), length=self.length() - 1)
                if self.contains(ga)]


def Cores(k, length=None, **kwargs):
    r"""
    A `k`-core is a partition from which no rim hook of size `k` can be removed.
    Alternatively, a `k`-core is an integer partition such that the Ferrers
    diagram for the partition contains no cells with a hook of size (a multiple of) `k`.

    The `k`-cores generally have two notions of size which are
    useful for different applications. One is the number of cells in the
    Ferrers diagram with hook less than `k`, the other is the total
    number of cells of the Ferrers diagram.  In the implementation in
    Sage, the first of notion is referred to as the ``length`` of the `k`-core
    and the second is the ``size`` of the `k`-core.  The class
    of Cores requires that either the size or the length of the elements in
    the class is specified.

    EXAMPLES:

    We create the set of the `4`-cores of length `6`. Here the length of a `k`-core is the size
    of the corresponding `(k-1)`-bounded partition, see also :meth:`~sage.combinat.core.Core.length`::

        sage: C = Cores(4, 6); C
        4-Cores of length 6
        sage: C.list()
        [[6, 3], [5, 2, 1], [4, 1, 1, 1], [4, 2, 2], [3, 3, 1, 1], [3, 2, 1, 1, 1], [2, 2, 2, 1, 1, 1]]
        sage: C.cardinality()
        7
        sage: C.an_element()
        [6, 3]

    We may also list the set of `4`-cores of size `6`, where the size is the number of boxes in the
    core, see also :meth:`~sage.combinat.core.Core.size`::

        sage: C = Cores(4, size=6); C
        4-Cores of size 6
        sage: C.list()
        [[4, 1, 1], [3, 2, 1], [3, 1, 1, 1]]
        sage: C.cardinality()
        3
        sage: C.an_element()
        [4, 1, 1]
    """
    if length is None and 'size' in kwargs:
        return Cores_size(k, kwargs['size'])
    if length is not None:
        return Cores_length(k, length)
    raise ValueError("you need to either specify the length or size of the cores considered")


class Cores_length(UniqueRepresentation, Parent):
    r"""
    The class of `k`-cores of length `n`.
    """

    def __init__(self, k, n):
        """
        TESTS::

            sage: C = Cores(3, 4)
            sage: TestSuite(C).run()
        """
        self.k = k
        self.n = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: repr(Cores(4, 3))  #indirect doctest
            '4-Cores of length 3'
        """
        return "%s-Cores of length %s" % (self.k, self.n)

    def list(self):
        r"""
        Return the list of all `k`-cores of length `n`.

        EXAMPLES::

            sage: C = Cores(3,4)
            sage: C.list()
            [[4, 2], [3, 1, 1], [2, 2, 1, 1]]
        """
        return [la.to_core(self.k - 1)
                for la in Partitions(self.n, max_part=self.k - 1)]

    def from_partition(self, part):
        r"""
        Convert the partition ``part`` into a core (as the identity map).

        This is the inverse method to :meth:`~sage.combinat.core.Core.to_partition`.

        EXAMPLES::

            sage: C = Cores(3,4)
            sage: c = C.from_partition([4,2]); c
            [4, 2]

            sage: mu = Partition([2,1,1])
            sage: C = Cores(3,3)
            sage: C.from_partition(mu).to_partition() == mu
            True

            sage: mu = Partition([])
            sage: C = Cores(3,0)
            sage: C.from_partition(mu).to_partition() == mu
            True
        """
        return Core(part, self.k)

    Element = Core


class Cores_size(UniqueRepresentation, Parent):
    r"""
    The class of `k`-cores of size `n`.
    """

    def __init__(self, k, n):
        """
        TESTS::

            sage: C = Cores(3, size = 4)
            sage: TestSuite(C).run()
        """
        self.k = k
        self.n = n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: repr(Cores(4, size = 3))  #indirect doctest
            '4-Cores of size 3'
        """
        return "%s-Cores of size %s" % (self.k, self.n)

    def list(self):
        r"""
        Return the list of all `k`-cores of size `n`.

        EXAMPLES::

            sage: C = Cores(3, size = 4)
            sage: C.list()
            [[3, 1], [2, 1, 1]]
        """
        return [Core(x, self.k) for x in Partitions(self.n)
                if x.is_core(self.k)]

    def from_partition(self, part):
        r"""
        Convert the partition ``part`` into a core (as the identity map).

        This is the inverse method to :meth:`to_partition`.

        EXAMPLES::

            sage: C = Cores(3,size=4)
            sage: c = C.from_partition([2,1,1]); c
            [2, 1, 1]

            sage: mu = Partition([2,1,1])
            sage: C = Cores(3,size=4)
            sage: C.from_partition(mu).to_partition() == mu
            True

            sage: mu = Partition([])
            sage: C = Cores(3,size=0)
            sage: C.from_partition(mu).to_partition() == mu
            True
        """
        return Core(part, self.k)

    Element = Core
