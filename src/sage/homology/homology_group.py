# sage_setup: distribution = sagemath-modules
"""
Homology Groups

This module defines a :meth:`HomologyGroup` class which is an abelian
group that prints itself in a way that is suitable for homology
groups.
"""
########################################################################
#       Copyright (C) 2013 John H. Palmieri <palmieri@math.washington.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
########################################################################

from sage.modules.free_module import VectorSpace
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup_fixed_gens
from sage.rings.integer_ring import ZZ


class HomologyGroup_class(AdditiveAbelianGroup_fixed_gens):
    """
    Discrete Abelian group on `n` generators. This class inherits from
    :class:`~sage.groups.additive_abelian.additive_abelian_group.AdditiveAbelianGroup_fixed_gens`;
    see :mod:`sage.groups.additive_abelian.additive_abelian_group` for more
    documentation. The main difference between the classes is in the print
    representation.

    EXAMPLES::

        sage: from sage.homology.homology_group import HomologyGroup
        sage: G = AbelianGroup(5, [5,5,7,8,9]); G                                       # needs sage.groups
        Multiplicative Abelian group isomorphic to C5 x C5 x C7 x C8 x C9
        sage: H = HomologyGroup(5, ZZ, [5,5,7,8,9]); H
        C5 x C5 x C7 x C8 x C9
        sage: G == loads(dumps(G))                                                      # needs sage.groups
        True
        sage: AbelianGroup(4)                                                           # needs sage.groups
        Multiplicative Abelian group isomorphic to Z x Z x Z x Z
        sage: HomologyGroup(4, ZZ)
        Z x Z x Z x Z
        sage: HomologyGroup(100, ZZ)
        Z^100
    """
    def __init__(self, n, invfac):
        """
        See :func:`HomologyGroup` for full documentation.

        EXAMPLES::

            sage: from sage.homology.homology_group import HomologyGroup
            sage: H = HomologyGroup(5, ZZ, [5,5,7,8,9]); H
            C5 x C5 x C7 x C8 x C9
        """
        n = len(invfac)
        A = ZZ ** n
        B = A.span([A.gen(i) * invfac[i] for i in range(n)])

        AdditiveAbelianGroup_fixed_gens.__init__(self, A, B, A.gens())
        self._original_invts = invfac

    def _repr_(self):
        """
        Print representation of ``self``.

        EXAMPLES::

            sage: from sage.homology.homology_group import HomologyGroup
            sage: H = HomologyGroup(7, ZZ, [4,4,4,4,4,7,7])
            sage: H._repr_()
            'C4^5 x C7 x C7'
            sage: HomologyGroup(6, ZZ)
            Z^6
        """
        eldv = self._original_invts
        if len(eldv) == 0:
            return "0"
        rank = len([x for x in eldv if x == 0])
        torsion = sorted(x for x in eldv if x)
        if rank > 4:
            g = ["Z^%s" % rank]
        else:
            g = ["Z"] * rank
        if len(torsion) != 0:
            printed = []
            for t in torsion:
                numfac = torsion.count(t)
                too_many = (numfac > 4)
                if too_many:
                    if t not in printed:
                        g.append("C{}^{}".format(t, numfac))
                        printed.append(t)
                else:
                    g.append("C%s" % t)
        times = " x "
        return times.join(g)

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.homology.homology_group import HomologyGroup
            sage: H = HomologyGroup(7, ZZ, [4,4,4,4,4,7,7])
            sage: H._latex_()
            'C_{4}^{5} \\times C_{7} \\times C_{7}'
            sage: latex(HomologyGroup(6, ZZ))
            \ZZ^{6}
        """
        eldv = self._original_invts
        if len(eldv) == 0:
            return "0"
        rank = len([x for x in eldv if x == 0])
        torsion = sorted(x for x in eldv if x)
        if rank > 4:
            g = ["\\ZZ^{{{}}}".format(rank)]
        else:
            g = ["\\ZZ"] * rank
        if torsion:
            printed = []
            for t in torsion:
                numfac = torsion.count(t)
                too_many = (numfac > 4)
                if too_many:
                    if t not in printed:
                        g.append("C_{{{}}}^{{{}}}".format(t, numfac))
                        printed.append(t)
                else:
                    g.append("C_{{{}}}".format(t))
        times = " \\times "
        return times.join(g)


def HomologyGroup(n, base_ring, invfac=None):
    """
    Abelian group on `n` generators which represents a homology group in a
    fixed degree.

    INPUT:

    - ``n`` -- integer; the number of generators

    - ``base_ring`` -- ring; the base ring over which the homology is computed

    - ``inv_fac`` -- list of integers; the invariant factors -- ignored
      if the base ring is a field

    OUTPUT:

    A class that can represent the homology group in a fixed
    homological degree.

    EXAMPLES::

        sage: from sage.homology.homology_group import HomologyGroup
        sage: G = AbelianGroup(5, [5,5,7,8,9]); G                                       # needs sage.groups
        Multiplicative Abelian group isomorphic to C5 x C5 x C7 x C8 x C9
        sage: H = HomologyGroup(5, ZZ, [5,5,7,8,9]); H
        C5 x C5 x C7 x C8 x C9
        sage: AbelianGroup(4)                                                           # needs sage.groups
        Multiplicative Abelian group isomorphic to Z x Z x Z x Z
        sage: HomologyGroup(4, ZZ)
        Z x Z x Z x Z

        sage: # needs sage.libs.flint (otherwise timeout)
        sage: HomologyGroup(100, ZZ)
        Z^100
    """
    if base_ring.is_field():
        return VectorSpace(base_ring, n)

    # copied from AbelianGroup:
    if invfac is None:
        if isinstance(n, (list, tuple)):
            invfac = n
            n = len(n)
        else:
            invfac = []
    if len(invfac) < n:
        invfac = [0] * (n - len(invfac)) + invfac
    elif len(invfac) > n:
        raise ValueError("invfac (={}) must have length n (={})".format(invfac, n))
    return HomologyGroup_class(n, invfac)
