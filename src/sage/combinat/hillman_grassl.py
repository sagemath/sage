r"""
The Hillman-Grassl correspondence

This module implements weak reverse plane partitions and four
correspondences on them: the Hillman-Grassl correspondence and
its inverse, as well as the Sulzgruber correspondence and its
inverse (the Pak correspondence).

Fix a partition `\lambda`
(see :meth:`~sage.combinat.partition.Partition`).
We draw all partitions and tableaux in English notation.

A `\lambda`-*array* will mean a tableau of shape `\lambda` whose
entries are nonnegative integers. (No conditions on the order of
these entries are made. Note that `0` is allowed.)

A *weak reverse plane partition of shape* `\lambda` (short:
`\lambda`-*rpp*) will mean a `\lambda`-array whose entries weakly
increase along each row and weakly increase along each column.
(The name "weak reverse plane partition" comes from Stanley in
[EnumComb2]_ Section 7.22; other authors -- such as Pak
[Sulzgr2017]_, or Hillman and Grassl in [HilGra1976]_ -- just
call it a reverse plane partition.)

The Hillman-Grassl correspondence is a bijection from the
set of `\lambda`-arrays to the set of `\lambda`-rpps.
For its definition, see
:meth:`~sage.combinat.tableau.Tableau.hillman_grassl`;
for its inverse, see
:meth:`~sage.combinat.hillman_grassl.WeakReversePlanePartition.hillman_grassl_inverse`.

The Sulzgruber correspondence `\Phi_\lambda` and the Pak
correspondence `\xi_\lambda` are two further mutually
inverse bijections between the set of
`\lambda`-arrays and the set of `\lambda`-rpps.
They appear (sometimes with different definitions, but
defining the same maps) in [Pak2002]_, [Hopkins2017]_ and
[Sulzgr2017]_. For their definitions, see
:meth:`~sage.combinat.tableau.Tableau.sulzgruber_correspondence` and
:meth:`~sage.combinat.hillman_grassl.WeakReversePlanePartition.pak_correspondence`.

EXAMPLES:

We construct a `\lambda`-rpp for `\lambda = (3, 3, 1)`
(note that `\lambda` needs not be specified explicitly)::

    sage: p = WeakReversePlanePartition([[0, 1, 3], [2, 4, 4], [3]])
    sage: p.parent()
    Weak Reverse Plane Partitions

(This is the example in Section 7.22 of [EnumComb2]_.)

Next, we apply the inverse of the Hillman-Grassl correspondence
to it::

    sage: HGp = p.hillman_grassl_inverse(); HGp
    [[1, 2, 0], [1, 0, 1], [1]]
    sage: HGp.parent()
    Tableaux

This is a `\lambda`-array, encoded as a tableau. We can
recover our original `\lambda`-rpp from it using the
Hillman-Grassl correspondence::

    sage: HGp.hillman_grassl() == p
    True

We can also apply the Pak correspondence to our rpp::

    sage: Pp = p.pak_correspondence(); Pp
    [[2, 0, 1], [0, 2, 0], [1]]
    sage: Pp.parent()
    Tableaux

This is undone by the Sulzgruber correspondence::

    sage: Pp.sulzgruber_correspondence() == p
    True

These four correspondences can also be accessed as standalone
functions (:meth:`hillman_grassl_inverse`, :meth:`hillman_grassl`,
:meth:`pak_correspondence` and :meth:`sulzgruber_correspondence`)
that transform lists of lists into lists of lists;
this may be more efficient. For example, the above computation
of ``HGp`` can also be obtained as follows::

    sage: from sage.combinat.hillman_grassl import hillman_grassl_inverse
    sage: HGp_bare = hillman_grassl_inverse([[0, 1, 3], [2, 4, 4], [3]])
    sage: HGp_bare
    [[1, 2, 0], [1, 0, 1], [1]]
    sage: isinstance(HGp_bare, list)
    True

REFERENCES:

- [Gans1981]_
- [HilGra1976]_
- [EnumComb2]_
- [Pak2002]_
- [Sulzgr2017]_
- [Hopkins2017]_

AUTHORS:

- Darij Grinberg and Tom Roby (2018): Initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2018 Darij Grinberg <darijgrinberg@gmail.com>,
#                     2018 Tom Roby <tomrobyuconn@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.tableau import Tableau, Tableaux
from sage.categories.sets_cat import Sets
from sage.combinat.combinatorial_map import combinatorial_map


class WeakReversePlanePartition(Tableau):
    r"""
    A weak reverse plane partition (short: rpp).

    A weak reverse plane partition is a tableau with nonnegative
    entries that are weakly increasing in each row and weakly
    increasing in each column.

    EXAMPLES::

        sage: x = WeakReversePlanePartition([[0, 1, 1], [0, 1, 3], [1, 2, 2], [1, 2, 3], [2]]); x
        [[0, 1, 1], [0, 1, 3], [1, 2, 2], [1, 2, 3], [2]]
        sage: x.pp()
          0  1  1
          0  1  3
          1  2  2
          1  2  3
          2
        sage: x.shape()
        [3, 3, 3, 3, 1]
    """
    @staticmethod
    def __classcall_private__(cls, r):
        r"""
        Return an rpp object.

        EXAMPLES::

            sage: WeakReversePlanePartition([[1, 2], [1, 3], [1]])
            [[1, 2], [1, 3], [1]]

        TESTS::

            sage: a1 = [[1, 2], [1, 3], [1]]
            sage: a2 = [(1, 2), (1, 3), (1,)]
            sage: A1 = WeakReversePlanePartition(a1)
            sage: A2 = WeakReversePlanePartition(a2)
            sage: A3 = WeakReversePlanePartition(A1)
            sage: A4 = Tableau(A1)
            sage: A1 == A2 == A3 == A4
            True
        """
        try:
            r = list(map(tuple, r))
        except TypeError:
            raise TypeError("r must be a list of positive integers")
        return WeakReversePlanePartitions()(r)

    def __init__(self, parent, t):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R = WeakReversePlanePartition([[0, 1, 2], [0, 2]])
            sage: TestSuite(R).run()
        """
        if not isinstance(t, Tableau):
            t = [list(row) for row in t]
        else:
            t = list(t)

        Tableau.__init__(self, parent, t)

    @combinatorial_map(order=2, name='conjugate')
    def conjugate(self):
        """
        Return the conjugate of ``self``.

        EXAMPLES::

            sage: c = WeakReversePlanePartition([[1,1],[1,3],[2]]).conjugate(); c
            [[1, 1, 2], [1, 3]]
            sage: c.parent()
            Weak Reverse Plane Partitions
        """
        C = super().conjugate()
        return WeakReversePlanePartition(C)

    def hillman_grassl_inverse(self):
        r"""
        Return the image of the `\lambda`-rpp ``self`` under the
        inverse of the Hillman-Grassl correspondence (as a
        :class:`~sage.combinat.tableau.Tableau`).

        Fix a partition `\lambda`
        (see :meth:`~sage.combinat.partition.Partition`).
        We draw all partitions and tableaux in English notation.

        A `\lambda`-*array* will mean a tableau of shape `\lambda` whose
        entries are nonnegative integers. (No conditions on the order of
        these entries are made. Note that `0` is allowed.)

        A *weak reverse plane partition of shape* `\lambda` (short:
        `\lambda`-*rpp*) will mean a `\lambda`-array whose entries weakly
        increase along each row and weakly increase along each column.

        The inverse `H^{-1}` of the Hillman-Grassl correspondence (see
        (:meth:`~sage.combinat.tableau.Tableau.hillman_grassl` for the
        latter) sends a `\lambda`-rpp `\pi` to a `\lambda`-array
        `H^{-1}(\pi)` constructed recursively as follows:

        * If all entries of `\pi` are `0`, then `H^{-1}(\pi) = \pi`.

        * Otherwise, let `s` be the index of the leftmost column of `\pi`
          containing a nonzero entry. Write the `\lambda`-array `M`
          as `(m_{i, j})`.

        * Define a sequence `((i_1, j_1), (i_2, j_2), \ldots,
          (i_n, j_n))` of boxes in the diagram of `\lambda` (actually a
          lattice path made of northward and eastward steps) as follows:
          Let `(i_1, j_1)` be the bottommost box in the `s`-th column
          of `\pi`.
          If `(i_k, j_k)` is defined for some `k \geq 1`, then
          `(i_{k+1}, j_{k+1})` is constructed as follows:
          If `q_{i_k - 1, j_k}` is well-defined and equals `q_{i_k, j_k}`,
          then we set `(i_{k+1}, j_{k+1}) = (i_k - 1, j_k)`. Otherwise,
          we set `(i_{k+1}, j_{k+1}) = (i_k, j_k + 1)` if this is still
          a box of `\lambda`. Otherwise, the sequence ends here.

        * Let `\pi'` be the `\lambda`-rpp obtained from `\pi` by
          subtracting `1` from the `(i_k, j_k)`-th entry of `\pi` for each
          `k \in \{1, 2, \ldots, n\}`.

        * Let `N'` be the image `H^{-1}(\pi')` (which is already
          constructed by recursion).
          Then, `H^{-1}(\pi)` is obtained from `N'` by adding `1` to the
          `(i_n, s)`-th entry of `N'`.

        This construction appears in [HilGra1976]_ Section 6 (where
        `\lambda`-arrays are re-encoded as sequences of "hook number
        multiplicities") and [EnumComb2]_ Section 7.22.

        .. SEEALSO::

            :meth:`~sage.combinat.hillman_grassl.hillman_grassl_inverse`
            for the inverse of the Hillman-Grassl correspondence as a
            standalone function.

            :meth:`~sage.combinat.tableau.Tableau.hillman_grassl`
            for the inverse map.

        EXAMPLES::

            sage: a = WeakReversePlanePartition([[2, 2, 4], [2, 3, 4], [3, 5]])
            sage: a.hillman_grassl_inverse()
            [[2, 1, 1], [0, 2, 0], [1, 1]]
            sage: b = WeakReversePlanePartition([[1, 1, 2, 2], [1, 1, 2, 2], [2, 2, 3, 3], [2, 2, 3, 3]])
            sage: B = b.hillman_grassl_inverse(); B
            [[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: b.parent(), B.parent()
            (Weak Reverse Plane Partitions, Tableaux)

        Applying the inverse of the Hillman-Grassl correspondence
        to the transpose of a `\lambda`-rpp `M` yields the same
        result as applying it to `M` and then transposing the
        result ([Gans1981]_ Corollary 3.4)::

            sage: a = WeakReversePlanePartition([[1,3,5],[2,4]])
            sage: aic = a.hillman_grassl_inverse().conjugate()
            sage: aic == a.conjugate().hillman_grassl_inverse()
            True
        """
        return Tableau(hillman_grassl_inverse(list(self)))

    def pak_correspondence(self):
        r"""
        Return the image of the `\lambda`-rpp ``self`` under the Pak
        correspondence (as a :class:`~sage.combinat.tableau.Tableau`).

        See :mod:`~sage.combinat.hillman_grassl`.

        The Pak correspondence is the map `\xi_\lambda`
        from [Sulzgr2017]_ Section 7, and is the map
        `\xi_\lambda` from [Pak2002]_ Section 4.
        It is the inverse of the Sulzgruber correspondence
        (:meth:`sulzgruber_correspondence`).
        The following description of the Pak correspondence follows
        [Hopkins2017]_ (which denotes it by `\mathcal{RSK}^{-1}`):

        Fix a partition `\lambda`
        (see :meth:`~sage.combinat.partition.Partition`).
        We draw all partitions and tableaux in English notation.

        A `\lambda`-*array* will mean a tableau of shape `\lambda` whose
        entries are nonnegative integers. (No conditions on the order of
        these entries are made. Note that `0` is allowed.)

        A *weak reverse plane partition of shape* `\lambda` (short:
        `\lambda`-*rpp*) will mean a `\lambda`-array whose entries weakly
        increase along each row and weakly increase along each column.

        We shall also use the following notation:
        If `(u, v)` is a cell of `\lambda`, and if `\pi` is a
        `\lambda`-rpp, then:

        * the *lower bound* of `\pi` at `(u, v)` (denoted by
          `\pi_{<(u, v)}`) is defined to be
          `\max \{ \pi_{u-1, v} , \pi_{u, v-1} \}`
          (where `\pi_{0, v}` and `\pi_{u, 0}` are understood to
          mean `0`).

        * the *upper bound* of `\pi` at `(u, v)` (denoted by
          `\pi_{>(u, v)}`) is defined to be
          `\min \{ \pi_{u+1, v} , \pi_{u, v+1} \}`
          (where `\pi_{i, j}` is understood to mean `+ \infty`
          if `(i, j)` is not in `\lambda`; thus, the upper
          bound at a corner cell is `+ \infty`).

        * *toggling* `\pi` at `(u, v)` means replacing the entry
          `\pi_{u, v}` of `\pi` at `(u, v)` by
          `\pi_{<(u, v)} + \pi_{>(u, v)} - \pi_{u, v}`
          (this is well-defined as long as `(u, v)` is not a
          corner of `\lambda`).

        Note that every `\lambda`-rpp `\pi` and every cell
        `(u, v)` of `\lambda` satisfy
        `\pi_{<(u, v)} \leq \pi_{u, v} \leq \pi_{>(u, v)}`.
        Note that toggling a `\lambda`-rpp (at a cell that is not
        a corner) always results in a `\lambda`-rpp. Also,
        toggling is an involution).

        Note also that the lower bound of `\pi` at `(u, v)` is
        defined (and finite) even when `(u, v)` is not a cell of
        `\lambda`, as long as both `(u-1, v)` and `(u, v-1)` are
        cells of `\lambda`.

        The Pak correspondence `\Phi_\lambda` sends a `\lambda`-array
        `M = (m_{i, j})` to a `\lambda`-rpp `\Phi_\lambda(M)`. It
        is defined by recursion on `\lambda` (that is, we assume that
        `\Phi_\mu` is already defined for every partition `\mu`
        smaller than `\lambda`), and its definition proceeds as
        follows:

        * If `\lambda = \varnothing`, then `\Phi_\lambda` is the
          obvious bijection sending the only `\varnothing`-array
          to the only `\varnothing`-rpp.

        * Pick any corner `c = (i, j)` of `\lambda`, and let `\mu`
          be the result of removing this corner `c` from the partition
          `\lambda`. (The exact choice of `c` is immaterial.)

        * Let `M'` be what remains of `M` when the corner cell `c`
          is removed.

        * Let `\pi' = \Phi_\mu(M')`.

        * For each positive integer `k` such that `(i-k, j-k)` is a
          cell of `\lambda`, toggle `\pi'` at `(i-k, j-k)`.
          (All these togglings commute, so the order in which they
          are made is immaterial.)

        * Extend the `\mu`-rpp `\pi'` to a `\lambda`-rpp `\pi` by
          adding the cell `c` and writing the number
          `m_{i, j} - \pi'_{<(i, j)}` into this cell.

        * Set `\Phi_\lambda(M) = \pi`.

        .. SEEALSO::

            :meth:`~sage.combinat.hillman_grassl.pak_correspondence`
            for the Pak correspondence as a standalone function.

            :meth:`~sage.combinat.tableau.Tableau.sulzgruber_correspondence`
            for the inverse map.

        EXAMPLES::

            sage: a = WeakReversePlanePartition([[1, 2, 3], [1, 2, 3], [2, 4, 4]])
            sage: A = a.pak_correspondence(); A
            [[1, 0, 2], [0, 2, 0], [1, 1, 0]]
            sage: a.parent(), A.parent()
            (Weak Reverse Plane Partitions, Tableaux)

        Applying the Pak correspondence to the transpose of a
        `\lambda`-rpp `M` yields the same result as applying it to
        `M` and then transposing the result::

            sage: a = WeakReversePlanePartition([[1,3,5],[2,4]])
            sage: acc = a.pak_correspondence().conjugate()
            sage: acc == a.conjugate().pak_correspondence()
            True
        """
        return Tableau(pak_correspondence(list(self)))


class WeakReversePlanePartitions(Tableaux):
    r"""
    The set of all weak reverse plane partitions.
    """
    @staticmethod
    def __classcall_private__(cls, shape=None, **kwds):
        """
        Normalize input to ensure a unique representation and
        return the correct class based on input.

        The ``shape`` parameter is currently not implemented.

        EXAMPLES::

            sage: S1 = WeakReversePlanePartitions()
            sage: S2 = WeakReversePlanePartitions()
            sage: S1 is S2
            True

            sage: S1 = WeakReversePlanePartitions([4, 2, 2, 1])  # not tested (not implemented)
            sage: S2 = WeakReversePlanePartitions((4, 2, 2, 1))  # not tested (not implemented)
            sage: S1 is S2  # not tested (not implemented)
            True
        """
        if shape is not None:
            raise NotImplementedError("shape cannot be specified")
            # from sage.combinat.partition import Partition
            # return RibbonShapedTableaux_shape(Partition(shape))

        return super().__classcall__(cls, **kwds)

    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = WeakReversePlanePartitions()
            sage: TestSuite(S).run()
        """
        Tableaux.__init__(self, category=Sets())

    def _repr_(self):
        """
        TESTS::

            sage: WeakReversePlanePartitions()
            Weak Reverse Plane Partitions
        """
        return "Weak Reverse Plane Partitions"

    Element = WeakReversePlanePartition

    def an_element(self):
        r"""
        Return a particular element of the class.

        TESTS::

            sage: T = WeakReversePlanePartitions()
            sage: T.an_element()
            [[0, 0, 1, 2], [0, 1, 1], [0], [2]]
        """
        return self.element_class(self, [[0, 0, 1, 2], [0, 1, 1], [0], [2]])


def transpose(M):
    r"""
    Return the transpose of a `\lambda`-array.

    The transpose of a `\lambda`-array `(m_{i, j})` is the
    `\lambda^t`-array `(m_{j, i})`
    (where `\lambda^t` is the conjugate of the partition
    `\lambda`).

    EXAMPLES::

        sage: from sage.combinat.hillman_grassl import transpose
        sage: transpose([[1, 2, 3], [4, 5]])
        [[1, 4], [2, 5], [3]]
        sage: transpose([[5, 0, 3], [4, 1, 0], [7]])
        [[5, 4, 7], [0, 1], [3, 0]]

    TESTS::

        sage: transpose(((2, 1), (3,)))
        [[2, 3], [1]]
        sage: transpose([])
        []
        sage: transpose(WeakReversePlanePartition([[1, 2, 3], [4, 5]]))
        [[1, 4], [2, 5], [3]]
        sage: transpose(WeakReversePlanePartition([]))
        []
    """
    if not M:  # empty array
        return []
    l = len(M[0])
    res = []
    for j in range(l):
        lj = []
        for r in M:
            if len(r) <= j:
                break
            lj.append(r[j])
        res.append(lj)
    return res


def hillman_grassl(M):
    r"""
    Return the image of the `\lambda`-array ``M``
    under the Hillman-Grassl correspondence.

    The Hillman-Grassl correspondence is a bijection
    between the tableaux with nonnegative entries
    (otherwise arbitrary) and the weak reverse plane
    partitions with nonnegative entries.
    This bijection preserves the shape of the
    tableau. See :mod:`~sage.combinat.hillman_grassl`.

    See :meth:`~sage.combinat.tableau.Tableau.hillman_grassl`
    for a description of this map.

    .. SEEALSO::

        :meth:`hillman_grassl_inverse`

    EXAMPLES::

        sage: from sage.combinat.hillman_grassl import hillman_grassl
        sage: hillman_grassl([[2, 1, 1], [0, 2, 0], [1, 1]])
        [[2, 2, 4], [2, 3, 4], [3, 5]]
        sage: hillman_grassl([[1, 2, 0], [1, 0, 1], [1]])
        [[0, 1, 3], [2, 4, 4], [3]]
        sage: hillman_grassl([])
        []
        sage: hillman_grassl([[3, 1, 2]])
        [[3, 4, 6]]
        sage: hillman_grassl([[2, 2, 0], [1, 1, 1], [1]])
        [[1, 2, 4], [3, 5, 5], [4]]
        sage: hillman_grassl([[1, 1, 1, 1]]*3)
        [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]]

    TESTS::

        sage: hillman_grassl(((2, 2, 0), (1, 1, 1), (1,)))
        [[1, 2, 4], [3, 5, 5], [4]]
    """
    lam = [len(row) for row in M]
    l = len(lam)
    Mt = transpose(M)
    # We won't touch M from now on; we'll only modify Mt
    # in place.
    hook_mults = []
    for j, col_j in enumerate(Mt):
        col_j_hook_mults = []
        for r, entry in enumerate(col_j):
            if entry != 0:
                col_j_hook_mults += [(r, j)] * entry
        hook_mults += reversed(col_j_hook_mults)
    res = [[0] * rowlen for rowlen in lam]
    for (r, s) in reversed(hook_mults):
        i = r
        j = lam[r] - 1
        while True:
            old = res[i][j]
            res[i][j] += 1
            if i + 1 != l and j < lam[i + 1] and old == res[i + 1][j]:
                i += 1
            else:
                if j == s:
                    break
                j -= 1
    return res


def hillman_grassl_inverse(M):
    r"""
    Return the image of the `\lambda`-rpp ``M`` under the
    inverse of the Hillman-Grassl correspondence.

    See :mod:`~sage.combinat.hillman_grassl`.

    See
    :meth:`~sage.combinat.hillman_grassl.WeakReversePlanePartition.hillman_grassl_inverse`
    for a description of this map.

    .. SEEALSO::

        :meth:`hillman_grassl`

    EXAMPLES::

        sage: from sage.combinat.hillman_grassl import hillman_grassl_inverse
        sage: hillman_grassl_inverse([[2, 2, 4], [2, 3, 4], [3, 5]])
        [[2, 1, 1], [0, 2, 0], [1, 1]]
        sage: hillman_grassl_inverse([[0, 1, 3], [2, 4, 4], [3]])
        [[1, 2, 0], [1, 0, 1], [1]]

    Applying the inverse of the Hillman-Grassl correspondence
    to the transpose of a `\lambda`-rpp `M` yields the same
    result as applying it to `M` and then transposing the
    result ([Gans1981]_ Corollary 3.4)::

        sage: hillman_grassl_inverse([[1,3,5],[2,4]])
        [[1, 2, 2], [1, 1]]
        sage: hillman_grassl_inverse([[1,2],[3,4],[5]])
        [[1, 1], [2, 1], [2]]
        sage: hillman_grassl_inverse([[1, 2, 3], [1, 2, 3], [2, 4, 4]])
        [[1, 2, 0], [0, 1, 1], [1, 0, 1]]
        sage: hillman_grassl_inverse([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]])
        [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]

    TESTS::

        sage: a = [[3], [4], [6]]
        sage: hillman_grassl_inverse(a)
        [[3], [1], [2]]
        sage: a
        [[3], [4], [6]]

        sage: hillman_grassl_inverse(((1,2),(3,4),(5,)))
        [[1, 1], [2, 1], [2]]
    """
    lam = [len(row) for row in M]
    res = [[0] * rowlen for rowlen in lam]
    Mt = transpose(M)
    # We won't touch M from now on; we'll only modify Mt
    # in place.
    while True:
        for j, col_j in enumerate(Mt):
            if all(entry == 0 for entry in col_j):
                continue
            else:
                break
        else:  # all entries of Mt are 0.
            break
        # Now, j is the index of the leftmost nonzero column of
        # the array.
        s = j
        i = len(col_j) - 1  # We already have j = s.
        while True:
            old = col_j[i]
            col_j[i] -= 1
            if i > 0 and old == col_j[i - 1]:
                i -= 1
            else:
                j += 1
                if j == lam[i]:
                    break
                col_j = Mt[j]
        res[i][s] += 1
    return res


def sulzgruber_correspondence(M):
    r"""
    Return the image of a `\lambda`-array ``M``
    under the Sulzgruber correspondence.

    The Sulzgruber correspondence is the map `\Phi_\lambda`
    from [Sulzgr2017]_ Section 7, and is the map
    `\xi_\lambda^{-1}` from [Pak2002]_ Section 5.
    It is denoted by `\mathcal{RSK}` in [Hopkins2017]_.
    It is the inverse of the Pak correspondence
    (:meth:`pak_correspondence`).

    See :meth:`~sage.combinat.tableau.Tableau.sulzgruber_correspondence`
    for a description of this map.

    EXAMPLES::

        sage: from sage.combinat.hillman_grassl import sulzgruber_correspondence
        sage: sulzgruber_correspondence([[1, 0, 2], [0, 2, 0], [1, 1, 0]])
        [[1, 2, 3], [1, 2, 3], [2, 4, 4]]
        sage: sulzgruber_correspondence([[1, 1, 2], [0, 1, 0], [3, 0, 0]])
        [[1, 1, 4], [2, 3, 4], [4, 4, 4]]
        sage: sulzgruber_correspondence([[1, 0, 2], [0, 2, 0], [1, 1]])
        [[0, 2, 3], [1, 3, 3], [2, 4]]
        sage: sulzgruber_correspondence([[0, 2, 2], [1, 1], [2]])
        [[1, 2, 4], [1, 3], [3]]
        sage: sulzgruber_correspondence([[1, 1, 1, 1]]*3)
        [[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]]

    The Sulzgruber correspondence can actually be
    extended (by the same definition) to arrays
    of nonnegative reals rather than nonnegative
    integers. This implementation supports this::

        sage: sulzgruber_correspondence([[1/2, 0, 1], [0, 1, 0], [1/2, 1/2]])
        [[0, 1, 3/2], [1/2, 3/2, 3/2], [1, 2]]

    TESTS::

        sage: sulzgruber_correspondence([])
        []
        sage: sulzgruber_correspondence(((0, 2, 2), (1, 1), (2,)))
        [[1, 2, 4], [1, 3], [3]]
    """
    lam = [len(row) for row in M]
    l = len(lam)
    if l == 0:
        return []
    # Finding a corner of lam...
    lam_0 = lam[0]
    for i, lam_i in enumerate(lam):
        if lam_i != lam_0:
            i -= 1
            break
    j = lam_0 - 1
    # Now, i is the index of the last row of ``M`` that
    # has the same length as the first row; hence, (i, j)
    # is a corner of lam.
    x = M[i][j]
    N = [list(row) for row in M]

    # remove the corner (i, j):
    N[i].pop()
    if not N[i]:
        N.pop()

    N = sulzgruber_correspondence(N)

    # toggling and inserting the new entry:
    for k in range(min(i, j) + 1):
        u = i - k
        v = j - k
        if u > 0 and v > 0:
            lower_bound = max(N[u - 1][v], N[u][v - 1])
        elif u > 0:
            lower_bound = N[u - 1][v]
        elif v > 0:
            lower_bound = N[u][v - 1]
        else:
            lower_bound = 0
        if k > 0:
            val = N[u][v]
            upper_bound = min(N[u + 1][v], N[u][v + 1])
            N[u][v] = lower_bound + upper_bound - val
        else:
            if len(N) <= u:
                N.append([])
            N[u].append(lower_bound + x)

    return N


def pak_correspondence(M, copy=True):
    r"""
    Return the image of a `\lambda`-rpp ``M``
    under the Pak correspondence.

    The Pak correspondence is the map `\xi_\lambda`
    from [Sulzgr2017]_ Section 7, and is the map
    `\xi_\lambda` from [Pak2002]_ Section 4.
    It is the inverse of the Sulzgruber correspondence
    (:meth:`sulzgruber_correspondence`).

    See
    :meth:`~sage.combinat.hillman_grassl.WeakReversePlanePartition.pak_correspondence`
    for a description of this map.

    INPUT:

    - ``copy`` -- boolean (default: ``True``); if set to ``False``, the
      algorithm will mutate the input (but be more efficient)

    EXAMPLES::

        sage: from sage.combinat.hillman_grassl import pak_correspondence
        sage: pak_correspondence([[1, 2, 3], [1, 2, 3], [2, 4, 4]])
        [[1, 0, 2], [0, 2, 0], [1, 1, 0]]
        sage: pak_correspondence([[1, 1, 4], [2, 3, 4], [4, 4, 4]])
        [[1, 1, 2], [0, 1, 0], [3, 0, 0]]
        sage: pak_correspondence([[0, 2, 3], [1, 3, 3], [2, 4]])
        [[1, 0, 2], [0, 2, 0], [1, 1]]
        sage: pak_correspondence([[1, 2, 4], [1, 3], [3]])
        [[0, 2, 2], [1, 1], [2]]
        sage: pak_correspondence([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]])
        [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]

    The Pak correspondence can actually be
    extended (by the same definition) to "rpps"
    of nonnegative reals rather than nonnegative
    integers. This implementation supports this::

        sage: pak_correspondence([[0, 1, 3/2], [1/2, 3/2, 3/2], [1, 2]])
        [[1/2, 0, 1], [0, 1, 0], [1/2, 1/2]]

    TESTS::

        sage: pak_correspondence([])
        []
        sage: pak_correspondence(((1, 2, 4), (1, 3), (3,)))
        [[0, 2, 2], [1, 1], [2]]

        sage: a = [[0, 2, 3], [1, 3, 3], [2, 4]]
        sage: pak_correspondence(a)
        [[1, 0, 2], [0, 2, 0], [1, 1]]
        sage: a
        [[0, 2, 3], [1, 3, 3], [2, 4]]
        sage: pak_correspondence(a, copy=False)
        [[1, 0, 2], [0, 2, 0], [1, 1]]
        sage: a
        []
    """
    if not M:
        return []
    lam = [len(row) for row in M]
    # Finding a corner of lam...
    lam_0 = lam[0]
    for i, lam_i in enumerate(lam):
        if lam_i != lam_0:
            i -= 1
            break
    j = lam_0 - 1
    # Now, i is the index of the last row of ``M`` that
    # has the same length as the first row; hence, (i, j)
    # is a corner of lam.
    x = M[i][j]

    if copy:
        N = [list(row) for row in M]
        # make a deep copy of M to avoid vandalizing M
    else:
        N = M

    # toggling and obtaining x (the corner entry of the output):
    for k in range(min(i, j) + 1):
        u = i - k
        v = j - k
        if u > 0 and v > 0:
            lower_bound = max(N[u - 1][v], N[u][v - 1])
        elif u > 0:
            lower_bound = N[u - 1][v]
        elif v > 0:
            lower_bound = N[u][v - 1]
        else:
            lower_bound = 0
        if k > 0:
            val = N[u][v]
            upper_bound = min(N[u + 1][v], N[u][v + 1])
            N[u][v] = lower_bound + upper_bound - val
        else:
            x -= lower_bound

    # remove the corner (i, j):
    N[i].pop()
    if not N[i]:
        N.pop()

    N = pak_correspondence(N, copy=False)
    if len(N) <= i:
        N.append([])
    N[i].append(x)
    return N
