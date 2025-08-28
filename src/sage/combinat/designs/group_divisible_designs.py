# sage.doctest: needs sage.rings.finite_rings
r"""
Group-divisible designs (GDD)

This module gathers everything related to Group-Divisible Designs. The
constructions defined here can be accessed through ``designs.<tab>``::

    sage: designs.group_divisible_design(14,{4},{2})
    Group Divisible Design on 14 points of type 2^7

The main function implemented here is :meth:`group_divisible_design` (which
calls all others) and the main class is :class:`GroupDivisibleDesign`. The
following functions are available:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`group_divisible_design` | Return a `(v,K,G)`-Group Divisible Design.
    :func:`GDD_4_2` | Return a `(2q,\{4\},\{2\})`-GDD for `q` a prime power with `q\equiv 1\pmod{6}`.

Functions
---------
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import is_prime_power
from sage.misc.unknown import Unknown
from .incidence_structures import IncidenceStructure


def group_divisible_design(v, K, G, existence=False, check=False):
    r"""
    Return a `(v,K,G)`-Group Divisible Design.

    A `(v,K,G)`-GDD is a pair `(\mathcal G, \mathcal B)` where:

    - `\mathcal G` is a partition of `X=\bigcup \mathcal G` where `|X|=v`

    - `\forall S\in \mathcal G, |S| \in G`

    - `\forall S\in \mathcal B, |S| \in K`

    - `\mathcal G\cup \mathcal B` is a `(v,K\cup G)`-PBD

    For more information, see the documentation of
    :class:`~sage.combinat.designs.incidence_structures.GroupDivisibleDesign` or
    :class:`~sage.combinat.designs.bibd.PairwiseBalancedDesign`.

    INPUT:

    - ``v`` -- integer

    - ``K``, ``G`` -- sets of integers

    - ``existence`` -- boolean; instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`)

        - ``False`` -- meaning that the design does not exist

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    .. NOTE::

        The GDD returned by this function are defined on ``range(v)``, and its
        groups are sets of consecutive integers.

    EXAMPLES::

        sage: designs.group_divisible_design(14,{4},{2})
        Group Divisible Design on 14 points of type 2^7
    """
    G = list(set(G))
    K = list(set(K))

    blocks = None

    # from a (v+1,k,1)-BIBD
    if len(G) == 1 == len(K) and G[0] + 1 in K:
        from .bibd import balanced_incomplete_block_design
        k = K[0]
        if existence:
            return balanced_incomplete_block_design(v + 1, k, existence=True)
        BIBD = balanced_incomplete_block_design(v + 1, k)
        groups = [[x for x in S if x != v] for S in BIBD if v in S]
        d = {p: i for i, p in enumerate(sum(groups, []))}
        d[v] = v
        BIBD.relabel(d)
        groups = [list(range((k - 1) * i, (k - 1) * (i + 1)))
                  for i in range(v // (k - 1))]
        blocks = [S for S in BIBD if v not in S]

    # (v,{4},{2})-GDD
    elif (v % 2 == 0 and K == [4] and
          G == [2] and GDD_4_2(v // 2, existence=True)):
        if existence:
            return True
        return GDD_4_2(v // 2, check=check)

    # From a TD(k,g)
    elif len(G) == 1 == len(K) and K[0] * G[0] == v:
        from .orthogonal_arrays import transversal_design
        return transversal_design(k=K[0], n=G[0], existence=existence)

    if blocks:
        return GroupDivisibleDesign(v,
                                    groups=groups,
                                    blocks=blocks,
                                    G=G,
                                    K=K,
                                    check=check,
                                    copy=True)

    if existence:
        return Unknown
    raise NotImplementedError


def GDD_4_2(q, existence=False, check=True):
    r"""
    Return a `(2q,\{4\},\{2\})`-GDD for `q` a prime power with `q\equiv 1\pmod{6}`.

    This method implements Lemma VII.5.17 from [BJL99] (p.495).

    INPUT:

    - ``q`` -- integer

    - ``existence`` -- boolean; instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`)

        - ``False`` -- meaning that the design does not exist

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    EXAMPLES::

        sage: from sage.combinat.designs.group_divisible_designs import GDD_4_2
        sage: GDD_4_2(7,existence=True)
        True
        sage: GDD_4_2(7)
        Group Divisible Design on 14 points of type 2^7
        sage: GDD_4_2(8,existence=True)
        Unknown
        sage: GDD_4_2(8)
        Traceback (most recent call last):
        ...
        NotImplementedError
    """
    if q <= 1 or q % 6 != 1 or not is_prime_power(q):
        if existence:
            return Unknown
        raise NotImplementedError
    if existence:
        return True

    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    G = FiniteField(q, 'x')
    w = G.primitive_element()
    e = w**((q - 1) // 3)

    # A first parallel class is defined. G acts on it, which yields all others.
    first_class = [[(0, 0), (1, w**i), (1, e * w**i), (1, e * e * w**i)]
                   for i in range((q - 1) // 6)]

    label = {p: i for i, p in enumerate(G)}
    classes = [[[2 * label[x[1] + g] + (x[0] + j) % 2 for x in S]
                for S in first_class]
               for g in G for j in range(2)]

    return GroupDivisibleDesign(2 * q,
        groups=[[i, i + 1] for i in range(0, 2 * q, 2)],
        blocks=sum(classes, []),
        K=[4],
        G=[2],
        check=check,
        copy=False)


class GroupDivisibleDesign(IncidenceStructure):
    r"""
    Group Divisible Design (GDD).

    Let `K` and `G` be sets of positive integers and let `\lambda` be a positive
    integer. A Group Divisible Design of index `\lambda` and order `v` is a
    triple `(V,\mathcal G,\mathcal B)` where:

    - `V` is a set of cardinality `v`

    - `\mathcal G` is a partition of `V` into groups whose size belongs to `G`

    - `\mathcal B` is a family of subsets of `V` whose size belongs to `K` such
      that any two points `p_1,p_2\in V` from different groups appear
      simultaneously in exactly `\lambda` elements of `\mathcal B`. Besides, a
      group and a block intersect on at most one point.

    If `K=\{k_1,...,k_l\}` and `G` has exactly `m_i` groups of cardinality `k_i`
    then `G` is said to have type `k_1^{m_1}...k_l^{m_l}`.

    INPUT:

    - ``points`` -- the underlying set. If ``points`` is an integer `v`, then
      the set is considered to be `\{0, ..., v-1\}`

    - ``groups`` -- the groups of the design. Set to ``None`` for an automatic
      guess (this triggers ``check=True`` and can thus cost some time)

    - ``blocks`` -- collection of blocks

    - ``G`` -- list of integers of which the sizes of the groups must be
      elements. Set to ``None`` (automatic guess) by default

    - ``K`` -- list of integers of which the sizes of the blocks must be
      elements. Set to ``None`` (automatic guess) by default

    - ``lambd`` -- integer (default: `1`); value of `\lambda`

    - ``check`` -- boolean (default: ``True``); whether to check that the design is indeed a `GDD`
      with the right parameters

    - ``copy`` -- (use with caution) if set to ``False`` then ``blocks`` must be
      a list of lists of integers. The list will not be copied but will be
      modified in place (each block is sorted, and the whole list is
      sorted). Your ``blocks`` object will become the instance's internal data.

    EXAMPLES::

        sage: from sage.combinat.designs.group_divisible_designs import GroupDivisibleDesign
        sage: TD = designs.transversal_design(4,10)
        sage: groups = [list(range(i*10,(i+1)*10)) for i in range(4)]
        sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
        Group Divisible Design on 40 points of type 10^4

    With unspecified groups::

        sage: # needs sage.schemes
        sage: D = designs.transversal_design(4,3).relabel(list('abcdefghiklm'),inplace=False).blocks()
        sage: GDD = GroupDivisibleDesign('abcdefghiklm',None,D)
        sage: sorted(GDD.groups())
        [['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h', 'i'], ['k', 'l', 'm']]
    """
    def __init__(self, points, groups, blocks, G=None, K=None, lambd=1,
                 check=True, copy=True, **kwds):
        r"""
        Constructor function.

        EXAMPLES::

            sage: from sage.combinat.designs.group_divisible_designs import GroupDivisibleDesign
            sage: TD = designs.transversal_design(4,10)
            sage: groups = [list(range(i*10,(i+1)*10)) for i in range(4)]
            sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
            Group Divisible Design on 40 points of type 10^4
        """
        from .designs_pyx import is_group_divisible_design

        self._lambd = lambd

        IncidenceStructure.__init__(self,
                                    points,
                                    blocks,
                                    copy=copy,
                                    check=False,
                                    **kwds)

        if groups is None or (copy is False and self._point_to_index is None):
            self._groups = groups
        elif self._point_to_index is None:
            self._groups = [g[:] for g in groups]
        else:
            self._groups = [[self._point_to_index[x] for x in g] for g in groups]

        if check or groups is None:
            is_gdd = is_group_divisible_design(self._groups, self._blocks,
                                               self.num_points(), G, K,
                                               lambd, verbose=1)
            assert is_gdd
            if groups is None:
                self._groups = is_gdd[1]

    def groups(self):
        r"""
        Return the groups of the Group-Divisible Design.

        EXAMPLES::

            sage: from sage.combinat.designs.group_divisible_designs import GroupDivisibleDesign
            sage: TD = designs.transversal_design(4,10)
            sage: groups = [list(range(i*10,(i+1)*10)) for i in range(4)]
            sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
            Group Divisible Design on 40 points of type 10^4
            sage: GDD.groups()
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
             [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
             [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
             [30, 31, 32, 33, 34, 35, 36, 37, 38, 39]]

        TESTS:

        Non-integer ground set::

            sage: # needs sage.schemes
            sage: TD = designs.transversal_design(5,5)
            sage: TD.relabel({i: chr(97+i) for i in range(25)})
            sage: TD.groups()
            [['a', 'b', 'c', 'd', 'e'],
             ['f', 'g', 'h', 'i', 'j'],
             ['k', 'l', 'm', 'n', 'o'],
             ['p', 'q', 'r', 's', 't'],
             ['u', 'v', 'w', 'x', 'y']]
        """
        if self._point_to_index is None:
            return [list(g) for g in self._groups]
        else:
            return [[self._points[i] for i in g] for g in self._groups]

    def __repr__(self):
        r"""
        Return a string that describes ``self``.

        EXAMPLES::

            sage: from sage.combinat.designs.group_divisible_designs import GroupDivisibleDesign
            sage: TD = designs.transversal_design(4,10)
            sage: groups = [list(range(i*10,(i+1)*10)) for i in range(4)]
            sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
            Group Divisible Design on 40 points of type 10^4
        """
        group_sizes = [len(g) for g in self._groups]

        gdd_type = ("{}^{}".format(s, group_sizes.count(s))
                    for s in sorted(set(group_sizes)))
        gdd_type = ".".join(gdd_type)

        if not gdd_type:
            gdd_type = "1^0"

        v = self.num_points()

        return "Group Divisible Design on {} points of type {}".format(v, gdd_type)
