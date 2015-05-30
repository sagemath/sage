# -*- coding: utf-8 -*-
r"""
Torsion subgroups of elliptic curves over number fields (including `\QQ`)

AUTHORS:

- Nick Alexander: original implementation over `\QQ`
- Chris Wuthrich: original implementation over number fields
- John Cremona: rewrote p-primary part to use division
    polynomials, added some features, unified Number Field and `\QQ` code.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.all import (Integer, RationalField, ZZ)
import sage.groups.additive_abelian.additive_abelian_wrapper as groups

class EllipticCurveTorsionSubgroup(groups.AdditiveAbelianGroupWrapper):
    r"""
    The torsion subgroup of an elliptic curve over a number field.

    EXAMPLES:

    Examples over `\QQ`::

        sage: E = EllipticCurve([-4, 0]); E
        Elliptic Curve defined by y^2  = x^3 - 4*x over Rational Field
        sage: G = E.torsion_subgroup(); G
        Torsion Subgroup isomorphic to Z/2 + Z/2 associated to the Elliptic Curve defined by y^2  = x^3 - 4*x over Rational Field
        sage: G.order()
        4
        sage: G.gen(0)
        (-2 : 0 : 1)
        sage: G.gen(1)
        (0 : 0 : 1)
        sage: G.ngens()
        2

    ::

        sage: E = EllipticCurve([17, -120, -60, 0, 0]); E
        Elliptic Curve defined by y^2 + 17*x*y - 60*y = x^3 - 120*x^2 over Rational Field
        sage: G = E.torsion_subgroup(); G
        Torsion Subgroup isomorphic to Trivial group associated to the Elliptic Curve defined by y^2 + 17*x*y - 60*y = x^3 - 120*x^2 over Rational Field
        sage: G.gens()
        ()
        sage: e = EllipticCurve([0, 33076156654533652066609946884,0,\
        347897536144342179642120321790729023127716119338758604800,\
        1141128154369274295519023032806804247788154621049857648870032370285851781352816640000])
        sage: e.torsion_order()
        16

    Constructing points from the torsion subgroup::

        sage: E = EllipticCurve('14a1')
        sage: T = E.torsion_subgroup()
        sage: [E(t) for t in T]
        [(0 : 1 : 0),
        (9 : 23 : 1),
        (2 : 2 : 1),
        (1 : -1 : 1),
        (2 : -5 : 1),
        (9 : -33 : 1)]

    An example where the torsion subgroup is not cyclic::

        sage: E = EllipticCurve([0,0,0,-49,0])
        sage: T = E.torsion_subgroup()
        sage: [E(t) for t in T]
        [(0 : 1 : 0), (-7 : 0 : 1), (0 : 0 : 1), (7 : 0 : 1)]

    An example where the torsion subgroup is trivial::

        sage: E = EllipticCurve('37a1')
        sage: T = E.torsion_subgroup()
        sage: T
        Torsion Subgroup isomorphic to Trivial group associated to the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: [E(t) for t in T]
        [(0 : 1 : 0)]

    Examples over other Number Fields::

        sage: E=EllipticCurve('11a1')
        sage: K.<i>=NumberField(x^2+1)
        sage: EK=E.change_ring(K)
        sage: from sage.schemes.elliptic_curves.ell_torsion import EllipticCurveTorsionSubgroup
        sage: EllipticCurveTorsionSubgroup(EK)
        Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1

        sage: E=EllipticCurve('11a1')
        sage: K.<i>=NumberField(x^2+1)
        sage: EK=E.change_ring(K)
        sage: T = EK.torsion_subgroup()
        sage: T.ngens()
        1
        sage: T.gen(0)
        (5 : -6 : 1)

    Note: this class is normally constructed indirectly as follows::

        sage: T = EK.torsion_subgroup(); T
        Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1
        sage: type(T)
        <class 'sage.schemes.elliptic_curves.ell_torsion.EllipticCurveTorsionSubgroup_with_category'>


    AUTHORS:

    - Nick Alexander - initial implementation over `\QQ`.
    - Chris Wuthrich - initial implementation over number fields.
    - John Cremona - additional features and unification.
    """
    def __init__(self, E, algorithm=None):
        r"""
        Initialization function for EllipticCurveTorsionSubgroup class

        INPUT:

        - ``E`` - An elliptic curve defined over a number field (including `\Q`)

        - ``algorithm`` - (string, default None): If not None, must be one
                     of 'PARI', 'doud', 'lutz_nagell'.  For curves
                     defined over `\QQ`, PARI is then used with the
                     appropriate flag passed to PARI's ``elltors()``
                     function; this parameter is ignored for curves
                     whose base field is not `\QQ`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_torsion import EllipticCurveTorsionSubgroup
            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: EllipticCurveTorsionSubgroup(EK)
            Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1

        Note: this class is normally constructed indirectly as follows::

            sage: T = EK.torsion_subgroup(); T
            Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1
            sage: type(T)
            <class 'sage.schemes.elliptic_curves.ell_torsion.EllipticCurveTorsionSubgroup_with_category'>

            sage: T == loads(dumps(T))  # known bug, see http://trac.sagemath.org/sage_trac/ticket/11599#comment:7
            True
        """
        self.__E = E
        self.__K = E.base_field()

        pari_torsion_algorithms = ["pari","doud","lutz_nagell"]

        if self.__K is RationalField() and algorithm in pari_torsion_algorithms:
            flag = pari_torsion_algorithms.index(algorithm)

            G = self.__E.pari_curve().elltors(flag)
            order = G[0].python()
            structure = G[1].python()
            gens = G[2].python()

            self.__torsion_gens = [ self.__E(P) for P in gens ]
            from sage.groups.additive_abelian.additive_abelian_group import cover_and_relations_from_invariants
            groups.AdditiveAbelianGroupWrapper.__init__(self, self.__E(0).parent(), self.__torsion_gens, structure)
            return

        T1 = E(0) # these will be the two generators
        T2 = E(0)
        k1 = 1    # with their order
        k2 = 1

        # find a multiple of the order of the torsion group
        bound = E._torsion_bound(number_of_places=20)

        # now do prime by prime
        for p,e in bound.factor():
            ptor = E._p_primary_torsion_basis(p,e)
            # print p,'-primary part is ',ptor
            if len(ptor)>0:
                T1 += ptor[0][0]
                k1 *= p**(ptor[0][1])
            if len(ptor)>1:
                T2 += ptor[1][0]
                k2 *= p**(ptor[1][1])

        order = k1*k2
        if k1 == 1:
            structure = []
            gens = []
        elif k2 == 1:
            structure = [k1]
            gens = [T1]
        else:
            structure = [k1,k2]
            gens = [T1,T2]

        #self.__torsion_gens = gens
        self._structure = structure
        groups.AdditiveAbelianGroupWrapper.__init__(self, T1.parent(), [T1, T2], structure)


    def _repr_(self):
        """
        String representation of an instance of the EllipticCurveTorsionSubgroup class.

        EXAMPLES::

            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: T = EK.torsion_subgroup(); T._repr_()
            'Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1'
        """
        return "Torsion Subgroup isomorphic to %s associated to the %s" % (self.short_name(), self.__E)

    def __cmp__(self,other):
        r"""
        Compares two torsion groups by simply comparing the elliptic curves.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: tor  = E.torsion_subgroup()
            sage: tor == tor
            True
        """
        c = cmp(type(self), type(other))
        if c:
            return c
        return cmp(self.__E, other.__E)

    def curve(self):
        """
        Return the curve of this torsion subgroup.

        EXAMPLES::

            sage: E=EllipticCurve('11a1')
            sage: K.<i>=NumberField(x^2+1)
            sage: EK=E.change_ring(K)
            sage: T = EK.torsion_subgroup()
            sage: T.curve() is EK
            True
        """
        return self.__E

    @cached_method
    def points(self):
        """
        Return a list of all the points in this torsion subgroup.
        The list is cached.

        EXAMPLES::

            sage: K.<i>=NumberField(x^2 + 1)
            sage: E = EllipticCurve(K,[0,0,0,1,0])
            sage: tor = E.torsion_subgroup()
            sage: tor.points()
            [(0 : 1 : 0), (-i : 0 : 1), (0 : 0 : 1), (i : 0 : 1)]
        """
        return [x.element() for x in self]
