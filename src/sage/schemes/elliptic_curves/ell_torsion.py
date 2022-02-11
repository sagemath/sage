# -*- coding: utf-8 -*-
r"""
Torsion subgroups of elliptic curves over number fields (including `\QQ`)

AUTHORS:

- Nick Alexander: original implementation over `\QQ`
- Chris Wuthrich: original implementation over number fields
- John Cremona: rewrote p-primary part to use division
    polynomials, added some features, unified Number Field and `\QQ` code.
"""

# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.all import RationalField
import sage.groups.additive_abelian.additive_abelian_wrapper as groups
from sage.structure.richcmp import richcmp_method, richcmp


@richcmp_method
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
        [(0 : 1 : 0), (0 : 0 : 1), (-7 : 0 : 1), (7 : 0 : 1)]

    An example where the torsion subgroup is trivial::

        sage: E = EllipticCurve('37a1')
        sage: T = E.torsion_subgroup()
        sage: T
        Torsion Subgroup isomorphic to Trivial group associated to the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: [E(t) for t in T]
        [(0 : 1 : 0)]

    Examples over other Number Fields::

        sage: E = EllipticCurve('11a1')
        sage: K.<i> = NumberField(x^2+1)
        sage: EK = E.change_ring(K)
        sage: from sage.schemes.elliptic_curves.ell_torsion import EllipticCurveTorsionSubgroup
        sage: EllipticCurveTorsionSubgroup(EK)
        Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1

        sage: E = EllipticCurve('11a1')
        sage: K.<i> = NumberField(x^2+1)
        sage: EK = E.change_ring(K)
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
    def __init__(self, E):
        r"""
        Initialization function for EllipticCurveTorsionSubgroup class

        INPUT:

        - ``E`` - An elliptic curve defined over a number field (including `\Q`)

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_torsion import EllipticCurveTorsionSubgroup
            sage: E = EllipticCurve('11a1')
            sage: K.<i> = NumberField(x^2+1)
            sage: EK = E.change_ring(K)
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

        if self.__K is RationalField():
            G = self.__E.pari_curve().elltors()
            structure = G[1].sage()
            gens = G[2].sage()

            self.__torsion_gens = [self.__E(P) for P in gens]
            groups.AdditiveAbelianGroupWrapper.__init__(self, self.__E(0).parent(), self.__torsion_gens, structure)
            return

        T1 = E(0) # these will be the two generators
        T2 = E(0)
        k1 = 1    # with their order
        k2 = 1

        # find a multiple of the order of the torsion group
        bound = torsion_bound(E, number_of_places=20)

        # now do prime by prime
        for p, e in bound.factor():
            ptor = E._p_primary_torsion_basis(p, e)
            if ptor:
                T1 += ptor[0][0]
                k1 *= p**(ptor[0][1])
            if len(ptor) > 1:
                T2 += ptor[1][0]
                k2 *= p**(ptor[1][1])

        if k1 == 1:
            structure = []
            gens = []
        elif k2 == 1:
            structure = [k1]
            gens = [T1]
        else:
            structure = [k1, k2]
            gens = [T1, T2]

        #self.__torsion_gens = gens
        self._structure = structure
        groups.AdditiveAbelianGroupWrapper.__init__(self, T1.parent(),
                                                    [T1, T2], structure)

    def _repr_(self):
        """
        String representation of an instance of the EllipticCurveTorsionSubgroup class.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: K.<i> = NumberField(x^2+1)
            sage: EK = E.change_ring(K)
            sage: T = EK.torsion_subgroup(); T._repr_()
            'Torsion Subgroup isomorphic to Z/5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in i with defining polynomial x^2 + 1'
        """
        return "Torsion Subgroup isomorphic to %s associated to the %s" % (self.short_name(), self.__E)

    def __richcmp__(self, other, op):
        r"""
        Compare two torsion groups by simply comparing the elliptic curves.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: tor = E.torsion_subgroup()
            sage: tor == tor
            True
        """
        if type(self) != type(other):
            return NotImplemented
        return richcmp(self.__E, other.__E, op)

    def curve(self):
        """
        Return the curve of this torsion subgroup.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: K.<i> = NumberField(x^2+1)
            sage: EK = E.change_ring(K)
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

            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve(K,[0,0,0,1,0])
            sage: tor = E.torsion_subgroup()
            sage: tor.points()
            [(0 : 1 : 0), (0 : 0 : 1), (-i : 0 : 1), (i : 0 : 1)]
        """
        return [x.element() for x in self]


def torsion_bound(E, number_of_places=20):
    r"""
    Return an upper bound on the order of the torsion subgroup.

    INPUT:

    - ``E`` -- an elliptic curve over `\QQ` or a number field

    - ``number_of_places`` (positive integer, default = 20) -- the
        number of places that will be used to find the bound

    OUTPUT:

    (integer) An upper bound on the torsion order.

    ALGORITHM:

    An upper bound on the order of the torsion group of the elliptic
    curve is obtained by counting points modulo several primes of good
    reduction. Note that the upper bound returned by this function is
    a multiple of the order of the torsion group, and in general will
    be greater than the order.

    To avoid nontrivial arithmetic in the base field (in particular,
    to avoid having to compute the maximal order) we only use prime
    `P` above rational primes `p` which do not divide the discriminant
    of the equation order.

    EXAMPLES::

        sage: CDB = CremonaDatabase()
        sage: from sage.schemes.elliptic_curves.ell_torsion import torsion_bound
        sage: [torsion_bound(E) for E in CDB.iter([14])]
        [6, 6, 6, 6, 6, 6]
        sage: [E.torsion_order() for E in CDB.iter([14])]
        [6, 6, 2, 6, 2, 6]

    An example over a relative number field (see :trac:`16011`)::

        sage: R.<x> = QQ[]
        sage: F.<a> = QuadraticField(5)
        sage: K.<b> = F.extension(x^2-3)
        sage: E = EllipticCurve(K,[0,0,0,b,1])
        sage: E.torsion_subgroup().order()
        1

    An example of a base-change curve from `\QQ` to a degree 16 field::

        sage: from sage.schemes.elliptic_curves.ell_torsion import torsion_bound
        sage: f = PolynomialRing(QQ,'x')([5643417737593488384,0,
        ....:     -11114515801179776,0,-455989850911004,0,379781901872,
        ....:     0,14339154953,0,-1564048,0,-194542,0,-32,0,1])
        sage: K = NumberField(f,'a')
        sage: E = EllipticCurve(K, [1, -1, 1, 824579, 245512517])
        sage: torsion_bound(E)
        16
        sage: E.torsion_subgroup().invariants()
        (4, 4)
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.finite_field_constructor import GF
    from sage.schemes.elliptic_curves.constructor import EllipticCurve

    K = E.base_field()

    # Special case K = QQ

    if K is RationalField():
        bound = ZZ.zero()
        k = 0
        p = ZZ(2)  # so we start with 3
        E = E.integral_model()
        disc_E = E.discriminant()

        while k < number_of_places:
            p = p.next_prime()
            if p.divides(disc_E):
                continue
            k += 1
            Fp = GF(p)
            new_bound = E.reduction(p).cardinality()
            bound = bound.gcd(new_bound)
            if bound == 1:
                return bound
        return bound

    # In case K is a relative extension we absolutize:

    absK = K.absolute_field('a_')
    f = absK.defining_polynomial()
    abs_map = absK.structure()[1]

    # Ensure f is monic and in ZZ[x]

    f = f.monic()
    den = f.denominator()
    if den != 1:
        x = f.parent().gen()
        n = f.degree()
        f = den**n * f(x/den)
    disc_f = f.discriminant()
    d = K.absolute_degree()

    # Now f is monic in ZZ[x] of degree d and defines the extension K = Q(a)

    # Make sure that we have a model for E with coefficients in ZZ[a]

    E = E.integral_model()
    disc_E = E.discriminant().norm()
    ainvs = [abs_map(c) for c in E.a_invariants()]

    bound = ZZ.zero()
    k = 0
    p = ZZ(2)  # so we start with 3

    try:  # special case, useful for base-changes from QQ
        ainvs = [ZZ(ai)  for ai in ainvs]
        while k < number_of_places:
            p = p.next_prime()
            if p.divides(disc_E) or p.divides(disc_f):
                continue
            k += 1
            for fi, ei in f.factor_mod(p):
                di = fi.degree()
                Fp = GF(p)
                new_bound = EllipticCurve(Fp, ainvs).cardinality(extension_degree=di)
                bound = bound.gcd(new_bound)
                if bound == 1:
                    return bound
        return bound
    except (ValueError, TypeError):
        pass

    # General case

    while k < number_of_places:
        p = p.next_prime()
        if p.divides(disc_E) or p.divides(disc_f):
            continue
        k += 1
        for fi, ei in f.factor_mod(p):
            di = fi.degree()
            Fq = GF((p, di))
            ai = fi.roots(Fq, multiplicities=False)[0]

            def red(c):
                return Fq.sum(Fq(c[j]) * ai**j for j in range(d))
            new_bound = EllipticCurve([red(c) for c in ainvs]).cardinality()
            bound = bound.gcd(new_bound)
            if bound == 1:
                return bound
    return bound
