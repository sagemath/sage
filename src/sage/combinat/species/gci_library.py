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
from group_cycle_index_series import GroupCycleIndexSeriesRing
from sage.misc.cachefunc import cached_function
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.rings.rational_field import RationalField
from generating_series import _integers_from, CycleIndexSeriesRing
from sage.combinat.sf.sf import SymmetricFunctions

@cached_function
def LinearOrderWithReversalGroupCycleIndex():
    """
    Returns the $\mathfrak{S}_{2}$-cycle index of the species $L$
    of linearly-ordered sets with the order-reversing action of
    $\mathfrak{S}_{2}$.

    EXAMPLES::

        sage: L = gci.LinearOrderWithReversalGroupCycleIndex()
        sage: e,t = L.parent().basis().keys()
        sage: L[e].coefficients(6) == species.LinearOrderSpecies().cycle_index_series().coefficients(6)
        True
        sage: L[t].coefficients(6)
        [0, 0, p[2], p[2, 1], p[2, 2], p[2, 2, 1]]
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
        yield 0
        yield 0
        for n in _integers_from(1):
            yield p([2]*n)
            yield p([2]*n+[1])

    Lt = CISR(Ltgen())

    L = ge*Le + gt*Lt
    return L

@cached_function
def CyclicOrderWithReversalGroupCycleIndex():
    """
    Returns the $\mathfrak{S}_{2}$-cycle index of the species $L$
    of linearly-ordered sets with the order-reversing action of
    $\mathfrak{S}_{2}$.

    EXAMPLES::

        sage: C = gci.CyclicOrderWithReversalGroupCycleIndex()
        sage: e,t = C.parent().basis().keys()
        sage: C[e].coefficients(6) == species.CycleSpecies().cycle_index_series().coefficients(6)
        True
        sage: C[t].coefficients(6)
        [0, 0, 0, p[2, 1], 0, p[2, 2, 1]]
    """
    from sage.combinat.species.cycle_species import CycleSpecies
    
    QQ = RationalField()
    CISR = CycleIndexSeriesRing(QQ)

    S2 = SymmetricGroup(2)
    GCISR = GroupCycleIndexSeriesRing(S2)
    G = GCISR.basis()
    e,t = G.keys()
    ge = G[e]
    gt = G[t]

    Ce = CycleSpecies().cycle_index_series()
    
    def Ctgen():
        p = SymmetricFunctions(QQ).power()
        yield 0
        yield 0
        for n in _integers_from(1):
            yield 1/2*p([2]*n)+1/2*p([2]*(n-1)+[1]*2)
            yield p([2]*n+[1])

    Ct = CISR(Ctgen())

    C = ge*Ce + gt*Ct
    return C
