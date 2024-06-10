"""
Examples of group cycle indices
"""
#*****************************************************************************
#       Copyright (C) 2013 Andrew Gainer-Dewar <andrew.gainer.dewar@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from .group_cycle_index_series import GroupCycleIndexSeriesRing
from sage.misc.cachefunc import cached_function
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.rings.rational_field import RationalField
from .generating_series import _integers_from, CycleIndexSeriesRing
from sage.combinat.sf.sf import SymmetricFunctions

@cached_function
def LinearOrderWithReversalGroupCycleIndex():
    r"""
    Returns the $\mathfrak{S}_{2}$-cycle index of the species $L$
    of linearly-ordered sets with the order-reversing action of
    $\mathfrak{S}_{2}$.

    The quotient $L / \mathfrak{S}_{2}$ is exactly the species of chains, as in [BLL, Table 4.1].

    EXAMPLES::

        sage: from sage.combinat.species.group_cycle_index_series_library import LinearOrderWithReversalGroupCycleIndex
        sage: L = LinearOrderWithReversalGroupCycleIndex()
        sage: e,t = L.parent().basis().keys()
        sage: L[e].coefficients(6) == species.LinearOrderSpecies().cycle_index_series().coefficients(6)
        True
        sage: L[t].coefficients(6)
        [p[], p[1], p[2], p[2, 1], p[2, 2], p[2, 2, 1]]
        sage: L.quotient().coefficients(6)
        [p[],
         p[1],
         1/2*p[1, 1] + 1/2*p[2],
         1/2*p[1, 1, 1] + 1/2*p[2, 1],
         1/2*p[1, 1, 1, 1] + 1/2*p[2, 2],
         1/2*p[1, 1, 1, 1, 1] + 1/2*p[2, 2, 1]]
    """
    from sage.combinat.species.linear_order_species import LinearOrderSpecies

    QQ = RationalField()
    CISR = CycleIndexSeriesRing(QQ)

    S2 = SymmetricGroup(2)
    GCISR = GroupCycleIndexSeriesRing(S2)
    G = GCISR.basis()
    e,t = G.keys()
    ge = G[e]
    gt = G[t]

    Le = LinearOrderSpecies().cycle_index_series()

    def Ltgen():
        p = SymmetricFunctions(QQ).power()
        yield p[0]
        yield p[1]
        for n in _integers_from(1):
            yield p([2]*n)
            yield p([2]*n+[1])

    Lt = CISR(Ltgen())

    L = ge*Le + gt*Lt
    return L

@cached_function
def CyclicOrderWithReversalGroupCycleIndex():
    r"""
    Returns the $\mathfrak{S}_{2}$-cycle index of the species $C$
    of cyclically-ordered sets with the order-reversing action of
    $\mathfrak{S}_{2}$.

    The quotient $C / \mathfrak{S}_{2}$ is exactly the species of polygons, as in [BLL, Sec. 2.6 ex. 9].

    EXAMPLES::

        sage: from sage.combinat.species.group_cycle_index_series_library import CyclicOrderWithReversalGroupCycleIndex
        sage: C = CyclicOrderWithReversalGroupCycleIndex()
        sage: e,t = C.parent().basis().keys()
        sage: C[e].coefficients(6) == species.CycleSpecies().cycle_index_series().coefficients(6)
        True
        sage: C[t].coefficients(6)
        [0, p[1], 1/2*p[1, 1] + 1/2*p[2], p[2, 1], 1/2*p[2, 1, 1] + 1/2*p[2, 2], p[2, 2, 1]]
        sage: C.quotient().generating_series().counts(9)
        [0, 1, 1, 1, 3, 12, 60, 360, 2520]
        sage: C.quotient().coefficients(5)
        [0,
         p[1],
         1/2*p[1, 1] + 1/2*p[2],
         1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3],
         1/8*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] + 3/8*p[2, 2] + 1/4*p[4]]
    """
    from sage.combinat.species.cycle_species import CycleSpecies

    QQ = RationalField()
    CISR = CycleIndexSeriesRing(QQ)
    p = SymmetricFunctions(QQ).power()

    S2 = SymmetricGroup(2)
    GCISR = GroupCycleIndexSeriesRing(S2)
    G = GCISR.basis()
    e,t = G.keys()
    ge = G[e]
    gt = G[t]

    Ce = CycleSpecies().cycle_index_series()

    def Ctgen():
        yield 0
        yield p[1]
        for n in _integers_from(1):
            yield QQ(1)/QQ(2) * (p([2]*n) + p([2]*(n-1) + [1]*2))
            yield p([2]*n+[1])

    Ct = CISR(Ctgen())

    C = ge*Ce + gt*Ct
    return C
