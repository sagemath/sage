r"""
Elliptic curves over a general field

This module defines the class :class:`EllipticCurve_field`, based on
:class:`EllipticCurve_generic`, for elliptic curves over general fields.
"""
# *****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import sage.rings.abc
from sage.categories.number_fields import NumberFields
from sage.categories.finite_fields import FiniteFields
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.rational_field import QQ
from sage.misc.misc_c import prod
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_field

from .constructor import EllipticCurve
from .ell_curve_isogeny import EllipticCurveIsogeny, isogeny_codomain_from_kernel
from . import ell_generic


class EllipticCurve_field(ell_generic.EllipticCurve_generic, ProjectivePlaneCurve_field):

    def __init__(self, R, data, category=None):
        r"""
        Constructor for elliptic curves over fields.

        Identical to the constructor for elliptic curves over
        general rings, except for setting the default category
        to :class:`AbelianVarieties`.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [1,1])
            sage: E.category()
            Category of abelian varieties over Rational Field
            sage: E = EllipticCurve(GF(101), [1,1])
            sage: E.category()
            Category of abelian varieties over Finite Field of size 101
        """
        from sage.categories.schemes import AbelianVarieties
        if category is None:
            category = AbelianVarieties(R)
        super().__init__(R, data, category=category)

    base_field = ell_generic.EllipticCurve_generic.base_ring

    _point = EllipticCurvePoint_field

    # Twists: rewritten by John Cremona as follows:
    #
    # Quadratic twist allowed except when char=2, j=0
    # Quartic twist allowed only if j=1728!=0 (so char!=2,3)
    # Sextic  twist allowed only if j=0!=1728 (so char!=2,3)
    #
    # More complicated twists exist in theory for char=2,3 and
    # j=0=1728, but I have never worked them out or seen them used!
    #

    def genus(self):
        """
        Return 1 for elliptic curves.

        EXAMPLES::

            sage: E = EllipticCurve(GF(3), [0, -1, 0, -346, 2652])
            sage: E.genus()
            1

            sage: R = FractionField(QQ['z'])
            sage: E = EllipticCurve(R, [0, -1, 0, -346, 2652])
            sage: E.genus()
            1
        """
        return ZZ.one()

    r"""
    Twists: rewritten by John Cremona as follows:

    The following twists are implemented:

    - Quadratic twist:  except when char=2 and `j=0`.
    - Quartic twist: only if `j=1728\not=0` (so not if char=2,3).
    - Sextic  twist: only if `j=0\not=1728` (so not if char=2,3).

    More complicated twists exist in theory for char=2,3 and j=0=1728,
    but are not implemented.
    """

    def quadratic_twist(self, D=None):
        r"""
        Return the quadratic twist of this curve by ``D``.

        INPUT:

        - ``D`` -- (default: ``None``) the twisting parameter (see below)

        In characteristics other than 2, `D` must be nonzero, and the
        twist is isomorphic to ``self`` after adjoining `\sqrt(D)` to the
        base.

        In characteristic 2, `D` is arbitrary, and the twist is
        isomorphic to ``self`` after adjoining a root of `x^2+x+D` to the
        base.

        In characteristic 2 when `j=0`, this is not implemented.

        If the base field `F` is finite, `D` need not be specified,
        and the curve returned is the unique curve (up to isomorphism)
        defined over `F` isomorphic to the original curve over the
        quadratic extension of `F` but not over `F` itself.  Over
        infinite fields, an error is raised if `D` is not given.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve([GF(1103)(1), 0, 0, 107, 340]); E
            Elliptic Curve defined by y^2 + x*y  = x^3 + 107*x + 340
             over Finite Field of size 1103
            sage: F = E.quadratic_twist(-1); F
            Elliptic Curve defined by y^2  = x^3 + 1102*x^2 + 609*x + 300
             over Finite Field of size 1103
            sage: E.is_isomorphic(F)
            False
            sage: E.is_isomorphic(F, GF(1103^2,'a'))
            True

        A characteristic 2 example::

            sage: E = EllipticCurve(GF(2), [1,0,1,1,1])
            sage: E1 = E.quadratic_twist(1)
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1, GF(4,'a'))
            True

        Over finite fields, the twisting parameter may be omitted::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(2^10)
            sage: E = EllipticCurve(k, [a^2,a,1,a+1,1])
            sage: Et = E.quadratic_twist()
            sage: Et  # random (only determined up to isomorphism)
            Elliptic Curve defined
             by y^2 + x*y  = x^3 + (a^7+a^4+a^3+a^2+a+1)*x^2 + (a^8+a^6+a^4+1)
             over Finite Field in a of size 2^10
            sage: E.is_isomorphic(Et)
            False
            sage: E.j_invariant() == Et.j_invariant()
            True

            sage: # needs sage.rings.finite_rings
            sage: p = next_prime(10^10)
            sage: k = GF(p)
            sage: E = EllipticCurve(k, [1,2,3,4,5])
            sage: Et = E.quadratic_twist()
            sage: Et  # random (only determined up to isomorphism)
            Elliptic Curve defined
             by y^2  = x^3 + 7860088097*x^2 + 9495240877*x + 3048660957
             over Finite Field of size 10000000019
            sage: E.is_isomorphic(Et)
            False
            sage: k2 = GF(p^2,'a')
            sage: E.change_ring(k2).is_isomorphic(Et.change_ring(k2))
            True
        """
        K = self.base_ring()
        char = K.characteristic()

        if D is None:
            if K.is_finite():
                x = polygen(K)
                if char == 2:
                    # We find D such that x^2+x+D is irreducible. If the
                    # degree is odd we can take D=1; otherwise it suffices to
                    # consider odd powers of a generator.
                    D = K(1)
                    if K.degree() % 2 == 0:
                        D = K.gen()
                        a = D**2
                        while (x**2 + x + D).roots():
                            D *= a
                else:
                    # We could take a multiplicative generator but
                    # that might be expensive to compute; otherwise
                    # half the elements will do, and testing squares
                    # is very fast.
                    D = K.random_element()
                    while D.is_square():
                        D = K.random_element()
            else:
                raise ValueError("twisting parameter D must be specified over infinite fields.")
        else:
            try:
                D = K(D)
            except ValueError:
                raise ValueError("twisting parameter D must be in the base field.")

            if char != 2 and D.is_zero():
                raise ValueError("twisting parameter D must be nonzero when characteristic is not 2")

        if char != 2:
            b2,b4,b6,b8 = self.b_invariants()
            # E is isomorphic to  [0,b2,0,8*b4,16*b6]
            return EllipticCurve(K,[0,b2*D,0,8*b4*D**2,16*b6*D**3])

        # now char==2
        if self.j_invariant() != 0: # iff a1!=0
            a1,a2,a3,a4,a6 = self.ainvs()
            E0 = self.change_weierstrass_model(a1,a3/a1,0,(a1**2*a4+a3**2)/a1**3)
            # which has the form = [1,A2,0,0,A6]
            assert E0.a1() == K(1)
            assert E0.a3() == K(0)
            assert E0.a4() == K(0)
            return EllipticCurve(K,[1,E0.a2()+D,0,0,E0.a6()])
        else:
            raise ValueError("Quadratic twist not implemented in char 2 when j=0")

    def two_torsion_rank(self):
        r"""
        Return the dimension of the 2-torsion subgroup of
        `E(K)`.

        This will be 0, 1 or 2.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: E.two_torsion_rank()
            0
            sage: K.<alpha> = QQ.extension(E.division_polynomial(2).monic())            # needs sage.rings.number_field
            sage: E.base_extend(K).two_torsion_rank()                                   # needs sage.rings.number_field
            1
            sage: E.reduction(53).two_torsion_rank()
            2

        ::

            sage: E = EllipticCurve('14a1')
            sage: E.two_torsion_rank()
            1
            sage: f = E.division_polynomial(2).monic().factor()[1][0]
            sage: K.<alpha> = QQ.extension(f)                                           # needs sage.rings.number_field
            sage: E.base_extend(K).two_torsion_rank()                                   # needs sage.rings.number_field
            2

        ::

            sage: EllipticCurve('15a1').two_torsion_rank()
            2
        """
        f = self.division_polynomial(Integer(2))
        n = len(f.roots())+1
        return Integer(n).ord(Integer(2))

    def quartic_twist(self, D):
        r"""
        Return the quartic twist of this curve by `D`.

        INPUT:

        - ``D`` -- (must be nonzero) the twisting parameter

        .. NOTE::

            The characteristic must not be 2 or 3, and the `j`-invariant must be 1728.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve_from_j(GF(13)(1728)); E
            Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 13
            sage: E1 = E.quartic_twist(2); E1
            Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 13
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1, GF(13^2,'a'))
            False
            sage: E.is_isomorphic(E1, GF(13^4,'a'))
            True
        """
        K = self.base_ring()
        char = K.characteristic()
        D = K(D)

        if char == 2 or char == 3:
            raise ValueError("Quartic twist not defined in chars 2,3")

        if self.j_invariant() != K(1728):
            raise ValueError("Quartic twist not defined when j!=1728")

        if D.is_zero():
            raise ValueError("quartic twist requires a nonzero argument")

        c4,c6 = self.c_invariants()
        # E is isomorphic to  [0,0,0,-27*c4,0]
        assert c6 == 0
        return EllipticCurve(K,[0,0,0,-27*c4*D,0])

    def sextic_twist(self, D):
        r"""
        Return the sextic twist of this curve by `D`.

        INPUT:

        - ``D`` -- (must be nonzero) the twisting parameter

        .. NOTE::

            The characteristic must not be 2 or 3, and the `j`-invariant must be 0.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve_from_j(GF(13)(0)); E
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 13
            sage: E1 = E.sextic_twist(2); E1
            Elliptic Curve defined by y^2 = x^3 + 11 over Finite Field of size 13
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1, GF(13^2,'a'))
            False
            sage: E.is_isomorphic(E1, GF(13^4,'a'))
            False
            sage: E.is_isomorphic(E1, GF(13^6,'a'))
            True
        """
        K = self.base_ring()
        char = K.characteristic()
        D = K(D)

        if char == 2 or char == 3:
            raise ValueError("Sextic twist not defined in chars 2,3")

        if self.j_invariant() != K(0):
            raise ValueError("Sextic twist not defined when j!=0")

        if D.is_zero():
            raise ValueError("Sextic twist requires a nonzero argument")

        c4,c6 = self.c_invariants()
        # E is isomorphic to  [0,0,0,0,-54*c6]
        assert c4 == 0
        return EllipticCurve(K,[0,0,0,0,-54*c6*D])

    def is_quadratic_twist(self, other):
        r"""
        Determine whether this curve is a quadratic twist of another.

        INPUT:

        - ``other`` -- an elliptic curve with the same base field as ``self``

        OUTPUT:

        Either 0, if the curves are not quadratic twists, or `D` if
        ``other`` is ``self.quadratic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        If the curves are defined over `\QQ`, the output `D` is
        a squarefree integer.

        .. NOTE::

            Not fully implemented in characteristic 2, or in
            characteristic 3 when both `j`-invariants are 0.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Et = E.quadratic_twist(-24)
            sage: E.is_quadratic_twist(Et)
            -6

            sage: E1 = EllipticCurve([0,0,1,0,0])
            sage: E1.j_invariant()
            0
            sage: E2 = EllipticCurve([0,0,0,0,2])
            sage: E1.is_quadratic_twist(E2)
            2
            sage: E1.is_quadratic_twist(E1)
            1
            sage: type(E1.is_quadratic_twist(E1)) == type(E1.is_quadratic_twist(E2))   # Issue #6574
            True

        ::

            sage: E1 = EllipticCurve([0,0,0,1,0])
            sage: E1.j_invariant()
            1728
            sage: E2 = EllipticCurve([0,0,0,2,0])
            sage: E1.is_quadratic_twist(E2)
            0
            sage: E2 = EllipticCurve([0,0,0,25,0])
            sage: E1.is_quadratic_twist(E2)
            5

        ::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(101)
            sage: E1 = EllipticCurve(F, [4,7])
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2); D != 0
            True
            sage: F = GF(101)
            sage: E1 = EllipticCurve(F, [4,7])
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2)
            sage: E1.quadratic_twist(D).is_isomorphic(E2)
            True
            sage: E1.is_isomorphic(E2)
            False
            sage: F2 = GF(101^2,'a')
            sage: E1.change_ring(F2).is_isomorphic(E2.change_ring(F2))
            True

        A characteristic 3 example::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(3^5,'a')
            sage: E1 = EllipticCurve_from_j(F(1))
            sage: E2 = E1.quadratic_twist(-1)
            sage: D = E1.is_quadratic_twist(E2); D != 0
            True
            sage: E1.quadratic_twist(D).is_isomorphic(E2)
            True

        ::

            sage: # needs sage.rings.finite_rings
            sage: E1 = EllipticCurve_from_j(F(0))
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2); D
            1
            sage: E1.is_isomorphic(E2)
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
        E = self
        F = other
        if not isinstance(E, EllipticCurve_generic) or not isinstance(F, EllipticCurve_generic):
            raise ValueError("arguments are not elliptic curves")
        K = E.base_ring()
        zero = K.zero()
        if not K == F.base_ring():
            return zero
        j = E.j_invariant()
        if j != F.j_invariant():
            return zero

        if E.is_isomorphic(F):
            if K is QQ:
                return ZZ(1)
            return K.one()

        char = K.characteristic()

        if char == 2:
            raise NotImplementedError("not implemented in characteristic 2")
        elif char == 3:
            if j == 0:
                raise NotImplementedError("not implemented in characteristic 3 for curves of j-invariant 0")
            D = E.b2()/F.b2()

        else:
            # now char!=2,3:
            c4E,c6E = E.c_invariants()
            c4F,c6F = F.c_invariants()

            if j == 0:
                um = c6E/c6F
                x = polygen(K)
                ulist = (x**3-um).roots(multiplicities=False)
                if not ulist:
                    D = zero
                else:
                    D = ulist[0]
            elif j == 1728:
                um = c4E/c4F
                x = polygen(K)
                ulist = (x**2-um).roots(multiplicities=False)
                if not ulist:
                    D = zero
                else:
                    D = ulist[0]
            else:
                D = (c6E*c4F)/(c6F*c4E)

        # Normalization of output:

        if D.is_zero():
            return D

        if K is QQ:
            D = D.squarefree_part()

        assert E.quadratic_twist(D).is_isomorphic(F)

        return D

    def is_quartic_twist(self, other):
        r"""
        Determine whether this curve is a quartic twist of another.

        INPUT:

        - ``other`` -- an elliptic curves with the same base field as ``self``

        OUTPUT:

        Either 0, if the curves are not quartic twists, or `D` if
        ``other`` is ``self.quartic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        .. NOTE::

            Not fully implemented in characteristics 2 or 3.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(1728))
            sage: E1 = E.quartic_twist(2)
            sage: D = E.is_quartic_twist(E1); D!=0
            True
            sage: E.quartic_twist(D).is_isomorphic(E1)
            True

        ::

            sage: E = EllipticCurve_from_j(1728)
            sage: E1 = E.quartic_twist(12345)
            sage: D = E.is_quartic_twist(E1); D
            15999120
            sage: (D/12345).is_perfect_power(4)
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
        E = self
        F = other
        if not isinstance(E, EllipticCurve_generic) or not isinstance(F, EllipticCurve_generic):
            raise ValueError("arguments are not elliptic curves")
        K = E.base_ring()
        zero = K.zero()
        if not K == F.base_ring():
            return zero
        j = E.j_invariant()
        if j != F.j_invariant() or j != K(1728):
            return zero

        if E.is_isomorphic(F):
            return K.one()

        char = K.characteristic()

        if char == 2:
            raise NotImplementedError("not implemented in characteristic 2")
        elif char == 3:
            raise NotImplementedError("not implemented in characteristic 3")
        else:
            # now char!=2,3:
            D = F.c4()/E.c4()

        if D.is_zero():
            return D

        assert E.quartic_twist(D).is_isomorphic(F)

        return D

    def is_sextic_twist(self, other):
        r"""
        Determine whether this curve is a sextic twist of another.

        INPUT:

        - ``other`` -- an elliptic curves with the same base field as ``self``

        OUTPUT:

        Either 0, if the curves are not sextic twists, or `D` if
        ``other`` is ``self.sextic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        .. NOTE::

            Not fully implemented in characteristics 2 or 3.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(0))
            sage: E1 = E.sextic_twist(2)
            sage: D = E.is_sextic_twist(E1); D != 0
            True
            sage: E.sextic_twist(D).is_isomorphic(E1)
            True

        ::

            sage: E = EllipticCurve_from_j(0)
            sage: E1 = E.sextic_twist(12345)
            sage: D = E.is_sextic_twist(E1); D
            575968320
            sage: (D/12345).is_perfect_power(6)
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
        E = self
        F = other
        if not isinstance(E, EllipticCurve_generic) or not isinstance(F, EllipticCurve_generic):
            raise ValueError("arguments are not elliptic curves")
        K = E.base_ring()
        zero = K.zero()
        if not K == F.base_ring():
            return zero
        j = E.j_invariant()
        if j != F.j_invariant() or not j.is_zero():
            return zero

        if E.is_isomorphic(F):
            return K.one()

        char = K.characteristic()

        if char == 2:
            raise NotImplementedError("not implemented in characteristic 2")
        elif char == 3:
            raise NotImplementedError("not implemented in characteristic 3")
        else:
            # now char!=2,3:
            D = F.c6()/E.c6()

        if D.is_zero():
            return D

        assert E.sextic_twist(D).is_isomorphic(F)

        return D

    def descend_to(self, K, f=None):
        r"""
        Given an elliptic curve ``self`` defined over a field `L` and a
        subfield `K` of `L`, return all elliptic curves over `K` which
        are isomorphic over `L` to ``self``.

        INPUT:

        - ``K`` -- a field which embeds into the base field `L` of ``self``

        - ``f`` -- (optional) an embedding of `K` into `L`; ignored if
          `K` is `\QQ`

        OUTPUT:

        A list (possibly empty) of elliptic curves defined over `K`
        which are isomorphic to ``self`` over `L`, up to isomorphism over `K`.

        .. NOTE::

            Currently only implemented over number fields.  To extend
            to other fields of characteristic not 2 or 3, what is
            needed is a method giving the preimages in `K^*/(K^*)^m` of
            an element of the base field, for `m=2,4,6`.

        EXAMPLES::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.descend_to(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Input must be a field.

        ::

            sage: # needs sage.rings.number_field
            sage: F.<b> = QuadraticField(23)
            sage: x = polygen(ZZ, 'x')
            sage: G.<a> = F.extension(x^3 + 5)
            sage: E = EllipticCurve(j=1728*b).change_ring(G)
            sage: EF = E.descend_to(F); EF
            [Elliptic Curve defined by y^2 = x^3 + (27*b-621)*x + (-1296*b+2484)
              over Number Field in b with defining polynomial x^2 - 23
              with b = 4.795831523312720?]
            sage: all(Ei.change_ring(G).is_isomorphic(E) for Ei in EF)
            True

        ::

            sage: # needs sage.rings.number_field
            sage: L.<a> = NumberField(x^4 - 7)
            sage: K.<b> = NumberField(x^2 - 7, embedding=a^2)
            sage: E = EllipticCurve([a^6, 0])
            sage: EK = E.descend_to(K); EK
            [Elliptic Curve defined by y^2 = x^3 + b*x over Number Field in b
              with defining polynomial x^2 - 7 with b = a^2,
             Elliptic Curve defined by y^2 = x^3 + 7*b*x over Number Field in b
              with defining polynomial x^2 - 7 with b = a^2]
            sage: all(Ei.change_ring(L).is_isomorphic(E) for Ei in EK)
            True

        ::

            sage: K.<a> = QuadraticField(17)                                            # needs sage.rings.number_field
            sage: E = EllipticCurve(j=2*a)                                              # needs sage.rings.number_field
            sage: E.descend_to(QQ)                                                      # needs sage.rings.number_field
            []

        TESTS:

        Check that :issue:`16456` is fixed::

            sage: # needs sage.rings.number_field
            sage: K.<a> = NumberField(x^3 - 2)
            sage: E = EllipticCurve('11a1').quadratic_twist(2)
            sage: EK = E.change_ring(K)
            sage: EK2 = EK.change_weierstrass_model((a,a,a,a+1))
            sage: EK2.descend_to(QQ)
            [Elliptic Curve defined by y^2 = x^3 + x^2 - 41*x - 199 over Rational Field]

            sage: k.<i> = QuadraticField(-1)                                            # needs sage.rings.number_field
            sage: E = EllipticCurve(k,[0,0,0,1,0])                                      # needs sage.rings.number_field
            sage: E.descend_to(QQ)                                                      # needs sage.rings.number_field
            [Elliptic Curve defined by y^2 = x^3 + x over Rational Field,
             Elliptic Curve defined by y^2 = x^3 - 4*x over Rational Field]
        """
        if not K.is_field():
            raise TypeError("Input must be a field.")
        L = self.base_field()
        if L is K:
            return self
        elif L == K:  # number fields can be equal but not identical
            return self.base_extend(K)

        # Construct an embedding f of K in L, and check that the
        # j-invariant is in the image, otherwise return an empty list:

        j = self.j_invariant()
        if K == QQ:
            try:
                jK = QQ(j)
            except (ValueError, TypeError):
                return []
        elif f is None:
            embeddings = K.embeddings(L)
            if not embeddings:
                raise TypeError("Input must be a subfield of the base field of the curve.")
            for g in embeddings:
                try:
                    jK = g.preimage(j)
                    f = g
                    break
                except Exception:
                    pass
            if f is None:
                return []
        else:
            try:
                if f.domain() != K:
                    raise ValueError("embedding has wrong domain")
                if f.codomain() != L:
                    raise ValueError("embedding has wrong codomain")
            except AttributeError:
                raise ValueError("invalid embedding: {}".format(f))
            try:
                jK = f.preimage(j)
            except Exception:
                return []

        # Now we have the j-invariant in K and must find all twists
        # which work, separating the cases of j=0 and j=1728.

        if L.characteristic():
            raise NotImplementedError("Not implemented in positive characteristic")

        if jK == 0:
            t = -54*self.c6()
            try:
                dlist = t.descend_mod_power(K,6)
                # list of d in K such that t/d is in L*^6
            except AttributeError:
                raise NotImplementedError("Not implemented over %s" % L)
            Elist = [EllipticCurve([0,0,0,0,d]) for d in dlist]
        elif jK == 1728:
            t = -27*self.c4()
            try:
                dlist = t.descend_mod_power(K,4)
                # list of d in K such that t/d is in L*^4
            except AttributeError:
                raise NotImplementedError("Not implemented over %s" % L)
            Elist = [EllipticCurve([0,0,0,d,0]) for d in dlist]
        else:
            c4, c6 = self.c_invariants()
            t = c6/c4
            try:
                dlist = t.descend_mod_power(K,2)
                # list of d in K such that t/d is in L*^2
            except AttributeError:
                raise NotImplementedError("Not implemented over %s" % L)
            c = -27*jK/(jK-1728) # =-27c4^3/c6^2
            a4list = [c*d**2 for d in dlist]
            a6list = [2*a4*d for a4,d in zip(a4list,dlist)]
            Elist = [EllipticCurve([0,0,0,a4,a6]) for a4,a6 in zip(a4list,a6list)]

        if K is QQ:
            Elist = [E.minimal_model() for E in Elist]
        return Elist

    def division_field(self, n, names='t', map=False, **kwds):
        r"""
        Given an elliptic curve over a number field or finite field `F` and
        a positive integer `n`, construct the `n`-division field `F(E[n])`.

        The `n`-division field is the smallest extension of `F` over which
        all `n`-torsion points of `E` are defined.

        INPUT:

        - ``n`` -- positive integer
        - ``names`` -- (default: ``'t'``) a variable name for the division field
        - ``map`` -- boolean (default: ``False``); also return an embedding of the
          :meth:`base_field` into the resulting field
        - ``kwds`` -- additional keyword arguments passed to
          :func:`~sage.rings.polynomial.polynomial_element.Polynomial.splitting_field`

        OUTPUT:

        If ``map`` is ``False``, the division field `K` as an absolute
        number field or a finite field.
        If ``map`` is ``True``, a tuple `(K, \phi)` where `\phi` is an
        embedding of the base field in the division field `K`.

        .. WARNING::

            This can take a very long time when the degree of the division
            field is large (e.g. when `n` is large or when the Galois
            representation is surjective).  The ``simplify`` flag also
            has a big influence on the running time over number fields:
            sometimes ``simplify=False`` is faster, sometimes the default
            ``simplify=True`` is faster.

        EXAMPLES:

        The 2-division field is the same as the splitting field of
        the 2-division polynomial (therefore, it has degree 1, 2, 3 or 6)::

            sage: # needs sage.rings.number_field
            sage: E = EllipticCurve('15a1')
            sage: K.<b> = E.division_field(2); K
            Number Field in b with defining polynomial x
            sage: E = EllipticCurve('14a1')
            sage: K.<b> = E.division_field(2); K
            Number Field in b with defining polynomial x^2 + 5*x + 92
            sage: E = EllipticCurve('196b1')
            sage: K.<b> = E.division_field(2); K
            Number Field in b with defining polynomial x^3 + x^2 - 114*x - 127
            sage: E = EllipticCurve('19a1')
            sage: K.<b> = E.division_field(2); K
            Number Field in b with defining polynomial
             x^6 + 10*x^5 + 24*x^4 - 212*x^3 + 1364*x^2 + 24072*x + 104292

        For odd primes `n`, the division field is either the splitting
        field of the `n`-division polynomial, or a quadratic extension
        of it. ::

            sage: # needs sage.rings.number_field
            sage: E = EllipticCurve('50a1')
            sage: F.<a> = E.division_polynomial(3).splitting_field(simplify_all=True); F
            Number Field in a
             with defining polynomial x^6 - 3*x^5 + 4*x^4 - 3*x^3 - 2*x^2 + 3*x + 3
            sage: K.<b> = E.division_field(3, simplify_all=True); K
            Number Field in b
             with defining polynomial x^6 - 3*x^5 + 4*x^4 - 3*x^3 - 2*x^2 + 3*x + 3

        If we take any quadratic twist, the splitting field of the
        3-division polynomial remains the same, but the 3-division field
        becomes a quadratic extension::

            sage: # needs sage.rings.number_field
            sage: E = E.quadratic_twist(5)  # 50b3
            sage: F.<a> = E.division_polynomial(3).splitting_field(simplify_all=True); F
            Number Field in a
             with defining polynomial x^6 - 3*x^5 + 4*x^4 - 3*x^3 - 2*x^2 + 3*x + 3
            sage: K.<b> = E.division_field(3, simplify_all=True); K
            Number Field in b with defining polynomial x^12 - 3*x^11 + 8*x^10 - 15*x^9
             + 30*x^8 - 63*x^7 + 109*x^6 - 144*x^5 + 150*x^4 - 120*x^3 + 68*x^2 - 24*x + 4

        Try another quadratic twist, this time over a subfield of `F`::

            sage: # needs sage.rings.number_field
            sage: G.<c>,_,_ = F.subfields(3)[0]
            sage: E = E.base_extend(G).quadratic_twist(c); E
            Elliptic Curve defined
             by y^2 = x^3 + 5*a0*x^2 + (-200*a0^2)*x + (-42000*a0^2+42000*a0+126000)
             over Number Field in a0 with defining polynomial x^3 - 3*x^2 + 3*x + 9
            sage: K.<b> = E.division_field(3, simplify_all=True); K
            Number Field in b with defining polynomial x^12 - 25*x^10 + 130*x^8 + 645*x^6 + 1050*x^4 + 675*x^2 + 225

        Some higher-degree examples::

            sage: # needs sage.rings.number_field
            sage: E = EllipticCurve('11a1')
            sage: K.<b> = E.division_field(2); K
            Number Field in b with defining polynomial
             x^6 + 2*x^5 - 48*x^4 - 436*x^3 + 1668*x^2 + 28792*x + 73844
            sage: K.<b> = E.division_field(3); K        # long time
            Number Field in b with defining polynomial x^48 ...
            sage: K.<b> = E.division_field(5); K
            Number Field in b with defining polynomial x^4 - x^3 + x^2 - x + 1
            sage: E.division_field(5, 'b', simplify=False)
            Number Field in b with defining polynomial x^4 + x^3 + 11*x^2 + 41*x + 101
            sage: E.base_extend(K).torsion_subgroup()   # long time
            Torsion Subgroup isomorphic to Z/5 + Z/5 associated to the Elliptic Curve
             defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
             over Number Field in b with defining polynomial x^4 - x^3 + x^2 - x + 1

            sage: # needs sage.rings.number_field
            sage: E = EllipticCurve('27a1')
            sage: K.<b> = E.division_field(3); K
            Number Field in b with defining polynomial x^2 + 3*x + 9
            sage: K.<b> = E.division_field(2); K
            Number Field in b with defining polynomial
             x^6 + 6*x^5 + 24*x^4 - 52*x^3 - 228*x^2 + 744*x + 3844
            sage: K.<b> = E.division_field(2, simplify_all=True); K
            Number Field in b with defining polynomial x^6 - 3*x^5 + 5*x^3 - 3*x + 1
            sage: K.<b> = E.division_field(5); K        # long time
            Number Field in b with defining polynomial x^48 ...
            sage: K.<b> = E.division_field(7); K        # long time
            Number Field in b with defining polynomial x^72 ...

        Over a number field::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([0,0,0,0,i])
            sage: L.<b> = E.division_field(2); L
            Number Field in b with defining polynomial x^4 - x^2 + 1
            sage: L.<b>, phi = E.division_field(2, map=True); phi
            Ring morphism:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Number Field in b with defining polynomial x^4 - x^2 + 1
              Defn: i |--> -b^3
            sage: L.<b>, phi = E.division_field(3, map=True)
            sage: L
            Number Field in b with defining polynomial x^24 - 6*x^22 - 12*x^21
             - 21*x^20 + 216*x^19 + 48*x^18 + 804*x^17 + 1194*x^16 - 13488*x^15
             + 21222*x^14 + 44196*x^13 - 47977*x^12 - 102888*x^11 + 173424*x^10
             - 172308*x^9 + 302046*x^8 + 252864*x^7 - 931182*x^6 + 180300*x^5
             + 879567*x^4 - 415896*x^3 + 1941012*x^2 + 650220*x + 443089
            sage: phi
            Ring morphism:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Number Field in b with defining polynomial x^24 ...
              Defn: i |--> -215621657062634529/183360797284413355040732*b^23 ...

        Over a finite field::

            sage: E = EllipticCurve(GF(431^2), [1,0])                                   # needs sage.rings.finite_rings
            sage: E.division_field(5, map=True)                                         # needs sage.rings.finite_rings
            (Finite Field in t of size 431^4,
             Ring morphism:
               From: Finite Field in z2 of size 431^2
               To:   Finite Field in t of size 431^4
               Defn: z2 |--> 52*t^3 + 222*t^2 + 78*t + 105)

        ::

            sage: E = EllipticCurve(GF(433^2), [1,0])                                   # needs sage.rings.finite_rings
            sage: K.<v> = E.division_field(7); K                                        # needs sage.rings.finite_rings
            Finite Field in v of size 433^16

        It also works for composite orders::

            sage: E = EllipticCurve(GF(11), [5,5])
            sage: E.change_ring(E.division_field(8)).abelian_group().torsion_subgroup(8).invariants()
            (8, 8)
            sage: E.change_ring(E.division_field(9)).abelian_group().torsion_subgroup(9).invariants()
            (9, 9)
            sage: E.change_ring(E.division_field(10)).abelian_group().torsion_subgroup(10).invariants()
            (10, 10)
            sage: E.change_ring(E.division_field(36)).abelian_group().torsion_subgroup(36).invariants()
            (36, 36)
            sage: E.change_ring(E.division_field(11)).abelian_group().torsion_subgroup(11).invariants()
            (11,)
            sage: E.change_ring(E.division_field(66)).abelian_group().torsion_subgroup(66).invariants()
            (6, 66)

        ...also over number fields::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)
            sage: E = EllipticCurve([0,0,0,0,i])
            sage: L,emb = E.division_field(6, names='b', map=True); L
            Number Field in b with defining polynomial x^24 + 12*x^23 + ...
            sage: E.change_ring(emb).torsion_subgroup().invariants()
            (6, 6)

        .. SEEALSO::

            To compute a basis of the `n`-torsion once the base field
            has been extended, you may use
            :meth:`sage.schemes.elliptic_curves.ell_number_field.EllipticCurve_number_field.torsion_subgroup`
            or
            :meth:`sage.schemes.elliptic_curves.ell_finite_field.EllipticCurve_finite_field.torsion_basis`.

        TESTS:

        Some random for prime orders::

            sage: # needs sage.rings.finite_rings
            sage: def check(E, l, K):
            ....:     EE = E.change_ring(K)
            ....:     cof = EE.order().prime_to_m_part(l)
            ....:     pts = (cof * EE.random_point() for _ in iter(int, 1))
            ....:     mul = lambda P: P if not l*P else mul(l*P)
            ....:     pts = map(mul, filter(bool, pts))
            ....:     if l == EE.base_field().characteristic():
            ....:         if EE.is_supersingular():
            ....:             Ps = ()
            ....:         else:
            ....:             assert l.divides(EE.order())
            ....:             Ps = (next(pts),)
            ....:     else:
            ....:         assert l.divides(EE.order())
            ....:         for _ in range(9999):
            ....:             P,Q = next(pts), next(pts)
            ....:             if P.weil_pairing(Q,l) != 1:
            ....:                 Ps = (P,Q)
            ....:                 break
            ....:         else:
            ....:             assert False
            ....:     deg = lcm(el.minpoly().degree() for el in sum(map(list,Ps),[]))
            ....:     assert max(deg, E.base_field().degree()) == K.degree()
            sage: q = next_prime_power(randrange(1, 10^9))
            sage: F.<a> = GF(q)
            sage: while True:
            ....:     try:
            ....:         E = EllipticCurve([F.random_element() for _ in range(5)])
            ....:     except ArithmeticError:
            ....:         continue
            ....:     break
            sage: l = random_prime(8)
            sage: K = E.division_field(l)
            sage: n = E.cardinality(extension_degree=K.degree()//F.degree())
            sage: (l^2 if q%l else 0 + E.is_ordinary()).divides(n)
            True
            sage: check(E, l, K)                # long time

        AUTHORS:

        - Jeroen Demeyer (2014-01-06): :issue:`11905`, use
          ``splitting_field`` method, moved from ``gal_reps.py``, make
          it work over number fields.
        - Lorenz Panny (2022): extend to finite fields
        - Lorenz Panny (2023): extend to composite `n`.
        """
        from sage.misc.verbose import verbose

        n = Integer(n)
        if n <= 0:
            raise ValueError("n must be a positive integer")

        verbose("Adjoining X-coordinates of %s-torsion points" % n)

        F = self.base_ring()
        f = self.division_polynomial(n).radical()

        if n == 2 or f.is_constant():
            # For n = 2, the division field is the splitting field of
            # the division polynomial.
            # If f is a nonzero constant, the n-torsion is trivial:
            # This means the curve must be supersingular and n == p.
            return f.splitting_field(names, map=map, **kwds)

        # We divide out the part defining points of non-maximal order.
        # Clearly all points of non-maximal order are multiples of points
        # of maximal order, so they cannot be defined over a larger field.
        if not n.is_prime():
            for d in n.prime_divisors():
                g = self.division_polynomial(n // d)
                f //= f.gcd(g)

        # Compute splitting field of X-coordinates.
        # The Galois group of the division field is a subgroup of GL(2,n).
        # The Galois group of the X-coordinates is a subgroup of GL(2,n)/{-1,+1}.
        if F in NumberFields():
            from sage.misc.misc_c import prod
            deg_mult = F.degree() * prod(l * (l+1) * (l-1)**2 * l**(4*(e-1)) for l,e in n.factor()) // 2
            K, F_to_K = f.splitting_field(names, degree_multiple=deg_mult, map=True, **kwds)
        elif F in FiniteFields():
            K, F_to_K = f.splitting_field('u', map=True, **kwds)
        else:
            raise NotImplementedError('only number fields and finite fields are currently supported')

        verbose("Adjoining Y-coordinates of %s-torsion points" % n)

        # THEOREM
        # (Cremona, https://github.com/sagemath/sage/issues/11905#comment:21)
        # (Later generalized to composite n by Lorenz Panny)
        #
        # Let K be a field, E an elliptic curve over K and n a positive
        # integer. Assume that K contains all roots of the n-division
        # polynomial of E, and that at least one point P of full order n
        # is defined over K. Then K contains all n-torsion points on E.
        #
        # PROOF. Let G be the absolute Galois group of K (every element
        # in it fixes all elements of K). For any n-torsion point Q
        # over the algebraic closure and any sigma in G, we must have
        # either sigma(Q) = Q or sigma(Q) = -Q (since K contains the
        # X-coordinate of Q). Similarly, sigma(P+Q) must equal either
        # P+Q or -(P+Q). However, since sigma is a group homomorphism,
        # we have sigma(P+Q) = sigma(P) + sigma(Q) = P + sigma(Q),
        # so either P + sigma(Q) = P+Q, which implies sigma(Q) = Q,
        # or P + sigma(Q) = -(P+Q), which implies sigma(Q) = -2P-Q.
        # The latter is impossible except for the easier case n = 2.
        # Hence, sigma(Q) = Q in all cases.
        #
        # This implies that it suffices to adjoin the Y-coordinate
        # of just one full-order point.

        x = f.change_ring(F_to_K).any_root(assume_squarefree=True)
        h = self.defining_polynomial().change_ring(F_to_K)(x, polygen(K), 1)
        L = h.splitting_field(names, map=map, **kwds)

        if map:
            L, K_to_L = L
            L = L, F_to_K.post_compose(K_to_L)
        return L

    def _Hom_(self, other, category=None):
        r"""
        Hook to make :class:`~sage.categories.homset.Hom`
        set the correct parent
        :class:`~sage.schemes.elliptic_curves.homset.EllipticCurveHomset`
        for
        :class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom`
        objects.

        EXAMPLES::

            sage: E = EllipticCurve(GF(19), [1,0])
            sage: type(E._Hom_(E))
            <class 'sage.schemes.elliptic_curves.homset.EllipticCurveHomset_with_category'>
        """
        if isinstance(other, ell_generic.EllipticCurve_generic) and self.base_ring() == other.base_ring():
            from . import homset
            return homset.EllipticCurveHomset(self, other, category=category)
        from sage.schemes.generic.homset import SchemeHomset_generic
        return SchemeHomset_generic(self, other, category=category)

    def isogeny(self, kernel, codomain=None, degree=None, model=None, check=True, algorithm=None, velu_sqrt_bound=None):
        r"""
        Return an elliptic-curve isogeny from this elliptic curve.

        The isogeny can be specified in two ways, by passing either a
        polynomial or a set of torsion points.  The methods used are:

        - Factored Isogenies (see
          :mod:`~sage.schemes.elliptic_curves.hom_composite`):
          Given a point, or a list of points which generate a
          composite-order subgroup, decomposes the isogeny into
          prime-degree steps. This can be used to construct isogenies
          of extremely large, smooth degree. When applicable, this
          algorithm is selected as default (see below). After factoring
          the degree single isogenies are computed using the other
          methods.
          This algorithm is selected using ``algorithm="factored"``.

        - Vélu's Formulas: Vélu's original formulas for computing
          isogenies.  This algorithm is selected by giving as the
          ``kernel`` parameter a single point generating a finite
          subgroup.

        - Kohel's Formulas: Kohel's original formulas for computing
          isogenies.  This algorithm is selected by giving as the
          ``kernel`` parameter a monic polynomial (or a coefficient list
          in little endian) which will define the kernel of the isogeny.
          Kohel's algorithm is currently only implemented for cyclic
          isogenies, with the exception of `[2]`.

        - √élu Algorithm (see
          :mod:`~sage.schemes.elliptic_curves.hom_velusqrt`):
          A variant of Vélu's formulas with essentially square-root
          instead of linear complexity (in the degree). Currently only
          available over finite fields. The input must be a single
          kernel point of odd order `\geq 5`.
          This algorithm is selected using ``algorithm="velusqrt"``.

        INPUT:

        - ``kernel`` -- a kernel; either a point on this curve, a list of
          points on this curve, a monic kernel polynomial, or ``None``.
          If initializing from a codomain, this must be ``None``.

        - ``codomain`` -- an elliptic curve (default: ``None``).

          - If ``kernel`` is ``None``, then ``degree`` must be given as well
            and the given ``codomain`` must be the codomain of a cyclic,
            separable, normalized isogeny of the given degree.

          - If ``kernel`` is not ``None``, then this must be isomorphic to
            the codomain of the separable isogeny defined by ``kernel``; in
            this case, the isogeny is post-composed with an isomorphism so
            that the codomain equals the given curve.

        - ``degree`` -- integer (default: ``None``).

          - If ``kernel`` is ``None``, then this is the degree of the isogeny
            from this curve to ``codomain``.

          - If ``kernel`` is not ``None``, then this is used to determine
            whether or not to skip a `\gcd` of the given kernel polynomial
            with the two-torsion polynomial of this curve.

        - ``model`` -- string (default: ``None``); supported values
          (cf. :func:`~sage.schemes.elliptic_curves.ell_field.compute_model`):

          - ``'minimal'``: if ``self`` is a curve over the rationals or
            over a number field, then the codomain is a global minimal
            model where this exists.

          - ``'short_weierstrass'``: the codomain is a short Weierstrass curve,
            assuming one exists.

          - ``'montgomery'``: the codomain is an (untwisted) Montgomery
            curve, assuming one exists over this field.

        - ``check`` -- boolean (default: ``True``); check whether the input is valid.
          Setting this to ``False`` can lead to significant speedups.

        - ``algorithm`` -- string (optional); the possible choices are:

          - ``'velusqrt'``: Use
            :class:`~sage.schemes.elliptic_curves.hom_velusqrt.EllipticCurveHom_velusqrt`.

          - ``'factored'``: Use
            :class:`~sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite`
            to decompose the isogeny into prime-degree steps.

          - ``'traditional'``: Use
            :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`.

          When ``algorithm`` is not specified, and ``kernel`` is not ``None``, an
          algorithm is selected using the following criteria:

            - if ``kernel`` is a list of multiple points, ``'factored'`` is selected.

            - If ``kernel`` is a single point, or a list containing a single point:

              - if the order of the point is unknown, ``'traditional'`` is selected.

              - If the order is known and composite, ``'factored'`` is selected.

              - If the order is known and prime, a choice between ``'velusqrt'`` and
                ``'traditional'`` is done according to the ``velu_sqrt_bound``
                parameter (see below).

          If none of the previous apply, ``'traditional'`` is selected.

        - ``velu_sqrt_bound`` -- integer (default: ``None``); establish the highest
          (prime) degree for which the ``'traditional'`` algorithm should be selected
          instead of ``'velusqrt'``. If ``None``, the default value from
          :class:`~sage.schemes.elliptic_curves.hom_velusqrt._VeluBoundObj` is used.
          This value is initially set to 1000, but can be modified by the user.
          If an integer is supplied and the isogeny computation goes through the
          ``'factored'`` algorithm, the same integer is supplied to each factor.

        The ``degree`` parameter is not supported when an ``algorithm`` is
        specified.

        OUTPUT:

        An isogeny between elliptic curves. This is a morphism of curves.
        (In all cases, the returned object will be an instance of
        :class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom`.)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(2^5, 'alpha'); alpha = F.gen()
            sage: E = EllipticCurve(F, [1,0,1,1,1])
            sage: R.<x> = F[]
            sage: phi = E.isogeny(x + 1)
            sage: phi.rational_maps()
            ((x^2 + x + 1)/(x + 1), (x^2*y + x)/(x^2 + 1))

        ::

            sage: E = EllipticCurve('11a1')
            sage: P = E.torsion_points()[1]
            sage: E.isogeny(P)
            Isogeny of degree 5
             from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20
                  over Rational Field
               to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580
                  over Rational Field

        ::

            sage: E = EllipticCurve(GF(19),[1,1])
            sage: P = E(15,3); Q = E(2,12)
            sage: (P.order(), Q.order())
            (7, 3)
            sage: phi = E.isogeny([P,Q]); phi
            Composite morphism of degree 21 = 7*3:
              From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19
              To:   Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19
            sage: phi(E.random_point())  # all points defined over GF(19) are in the kernel
            (0 : 1 : 0)

        ::

            sage: E = EllipticCurve(GF(2^32 - 5), [170246996, 2036646110])              # needs sage.rings.finite_rings
            sage: P = E.lift_x(2)                                                       # needs sage.rings.finite_rings
            sage: E.isogeny(P, algorithm='factored')                                    # needs sage.rings.finite_rings
            Composite morphism of degree 1073721825 = 3^4*5^2*11*19*43*59:
              From: Elliptic Curve defined by y^2 = x^3 + 170246996*x + 2036646110
                     over Finite Field of size 4294967291
              To:   Elliptic Curve defined by y^2 = x^3 + 272790262*x + 1903695400
                     over Finite Field of size 4294967291

        Not all polynomials define a finite subgroup (:issue:`6384`)::

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: phi = E.isogeny([14,27,4,1])
            Traceback (most recent call last):
            ...
            ValueError: the polynomial x^3 + 4*x^2 + 27*x + 14 does not define a finite
            subgroup of Elliptic Curve defined by y^2 + x*y = x^3 + x + 2
            over Finite Field of size 31

        Order of the point known and composite::

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: P = E(26, 4)
            sage: assert P.order() == 12
            sage: print(P._order)
            12
            sage: E.isogeny(P)
            Composite morphism of degree 12 = 2^2*3:
              From: Elliptic Curve defined by y^2 + x*y = x^3 + x + 2 over Finite Field of size 31
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + 26*x + 8 over Finite Field of size 31

        ``kernel`` is a list of points::

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: P = E(21,2)
            sage: Q = E(7, 12)
            sage: print(P.order())
            6
            sage: print(Q.order())
            2
            sage: E.isogeny([P, Q])
            Composite morphism of degree 12 = 2*3*2:
              From: Elliptic Curve defined by y^2 + x*y = x^3 + x + 2 over Finite Field of size 31
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + 2*x + 26 over Finite Field of size 31

        Multiple ways to set the `velu_sqrt_bound`::

            sage: E = EllipticCurve_from_j(GF(97)(42))
            sage: P = E.gens()[0]*4
            sage: print(P.order())
            23
            sage: E.isogeny(P)
            Isogeny of degree 23 from Elliptic Curve defined by y^2 = x^3 + 6*x + 46 over Finite Field of size 97 to Elliptic Curve defined by y^2 = x^3 + 72*x + 29 over Finite Field of size 97
            sage: E.isogeny(P, velu_sqrt_bound=10)
            Elliptic-curve isogeny (using square-root Vélu) of degree 23:
              From: Elliptic Curve defined by y^2 = x^3 + 6*x + 46 over Finite Field of size 97
              To:   Elliptic Curve defined by y^2 = x^3 + 95*x + 68 over Finite Field of size 97
            sage: from sage.schemes.elliptic_curves.hom_velusqrt import _velu_sqrt_bound
            sage: _velu_sqrt_bound.set(10)
            sage: E.isogeny(P)
            Elliptic-curve isogeny (using square-root Vélu) of degree 23:
              From: Elliptic Curve defined by y^2 = x^3 + 6*x + 46 over Finite Field of size 97
              To:   Elliptic Curve defined by y^2 = x^3 + 95*x + 68 over Finite Field of size 97
            sage: _velu_sqrt_bound.set(1000) # Reset bound

        If the order of the point is unknown, fall back to ``'traditional'``::

            sage: E = EllipticCurve_from_j(GF(97)(42))
            sage: P = E(2, 39)
            sage: from sage.schemes.elliptic_curves.hom_velusqrt import _velu_sqrt_bound
            sage: _velu_sqrt_bound.set(1)
            sage: E.isogeny(P)
            Isogeny of degree 46 from Elliptic Curve defined by y^2 = x^3 + 6*x + 46 over Finite Field of size 97 to Elliptic Curve defined by y^2 = x^3 + 87*x + 47 over Finite Field of size 97
            sage: _velu_sqrt_bound.set(1000) # Reset bound

        .. SEEALSO::

            - :class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom`
            - :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`
            - :class:`~sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite`

        TESTS:

        Until the checking of kernel polynomials was implemented in
        :issue:`23222`, the following raised no error but returned an
        invalid morphism.  See also :issue:`11578`::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2 - x - 1)
            sage: E = EllipticCurve(K, [-13392, -1080432])
            sage: R.<x> = K[]
            sage: phi = E.isogeny( (x-564)*(x - 396/5*a + 348/5) )
            Traceback (most recent call last):
            ...
            ValueError: the polynomial x^2 + (-396/5*a - 2472/5)*x + 223344/5*a - 196272/5 does not
            define a finite subgroup of Elliptic Curve defined by y^2 = x^3 + (-13392)*x + (-1080432)
            over Number Field in a with defining polynomial x^2 - x - 1

        We check that the cached order is correctly copied over::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(2^127 - 1), [1,2,3,4,5])
            sage: E.set_order(170141183460469231746191640949390434666)
            sage: phi = E.isogeny(E.lift_x(77347718128277853096420969229987528666))
            sage: phi.codomain()._order
            170141183460469231746191640949390434666

        Check that ``'factored'`` recursively apply `velu_sqrt_bound`::

            sage: from sage.schemes.elliptic_curves.hom_velusqrt import _velu_sqrt_bound
            sage: _velu_sqrt_bound.get()
            1000
            sage: _velu_sqrt_bound.set(50)
            sage: _velu_sqrt_bound.get()
            50
            sage: from sage.schemes.elliptic_curves import hom_composite
            sage: p = 3217
            sage: E = EllipticCurve_from_j(GF(p)(42))
            sage: P = E.gens()[0]
            sage: phis = hom_composite._compute_factored_isogeny_single_generator(P, velu_sqrt_bound=50)
            sage: for phi in phis:
            ....:     print(phi)
            ....:
            Isogeny of degree 31 from Elliptic Curve defined by y^2 = x^3 + 114*x + 544 over Finite Field of size 3217 to Elliptic Curve defined by y^2 = x^3 + 277*x + 1710 over Finite Field of size 3217
            Elliptic-curve isogeny (using square-root Vélu) of degree 103:
              From: Elliptic Curve defined by y^2 = x^3 + 277*x + 1710 over Finite Field of size 3217
              To:   Elliptic Curve defined by y^2 = x^3 + 2979*x + 1951 over Finite Field of size 3217
        """
        if algorithm is not None and degree is not None:
            raise TypeError('cannot pass "degree" and "algorithm" parameters simultaneously')
        if algorithm == "velusqrt":
            from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            return EllipticCurveHom_velusqrt(self, kernel, codomain=codomain, model=model)
        if algorithm == "factored":
            from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            return EllipticCurveHom_composite(self, kernel, codomain=codomain, model=model, velu_sqrt_bound=velu_sqrt_bound)
        if algorithm == "traditional":
            return EllipticCurveIsogeny(self, kernel, codomain, degree, model, check=check)

        if kernel is not None:
            # Check for multiple points or point of known order
            kernel_is_list = isinstance(kernel, (list, tuple))
            if kernel_is_list and kernel[0] in self and len(kernel) > 1:
                from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
                return EllipticCurveHom_composite(self, kernel, codomain=codomain, model=model, velu_sqrt_bound=velu_sqrt_bound)

            if not kernel_is_list or (len(kernel) == 1 and kernel[0] in self):
                # Single point on the curve; unpack the list for compatibility with velusqrt
                if kernel_is_list:
                    kernel = kernel[0]

                known_order = hasattr(kernel, "_order")

                if known_order and kernel._order.is_pseudoprime():
                    if not velu_sqrt_bound:
                        from sage.schemes.elliptic_curves.hom_velusqrt import _velu_sqrt_bound
                        velu_sqrt_bound = _velu_sqrt_bound.get()

                    if kernel._order > velu_sqrt_bound:
                        from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
                        return EllipticCurveHom_velusqrt(self, kernel, codomain=codomain, model=model)
                    # Otherwise fall back to the standard case
                elif known_order:
                    from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
                    return EllipticCurveHom_composite(self, kernel, codomain=codomain, model=model, velu_sqrt_bound=velu_sqrt_bound)
        try:
            return EllipticCurveIsogeny(self, kernel, codomain, degree, model, check=check)
        except AttributeError as e:
            raise RuntimeError("Unable to construct isogeny: %s" % e)

    def isogeny_codomain(self, kernel):
        r"""
        Return the codomain of the isogeny from ``self`` with given kernel.

        INPUT:

        - ``kernel`` -- either a list of points in the kernel of the isogeny,
          or a kernel polynomial (specified as either a univariate polynomial
          or a coefficient list)

        OUTPUT:

        An elliptic curve, the codomain of the separable normalized
        isogeny defined by this kernel.

        EXAMPLES::

            sage: E = EllipticCurve('17a1')
            sage: R.<x> = QQ[]
            sage: E2 = E.isogeny_codomain(x - 11/4); E2
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 1461/16*x - 19681/64
             over Rational Field

        TESTS:

        We check that the cached order is correctly copied over::

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(2^127 - 1), [1,2,3,4,5])
            sage: E.set_order(170141183460469231746191640949390434666)
            sage: E2 = E.isogeny_codomain(E.lift_x(77347718128277853096420969229987528666))
            sage: E2._order
            170141183460469231746191640949390434666
        """
        E = isogeny_codomain_from_kernel(self, kernel)
        if self.base_field().is_finite():
            E._fetch_cached_order(self)
        return E

    def period_lattice(self):
        r"""
        Return the period lattice of the elliptic curve for the given
        embedding of its base field with respect to the differential
        `dx/(2y + a_1x + a_3)`.

        Only supported for some base rings.

        EXAMPLES::

            sage: EllipticCurve(RR, [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 6.00000000000000 over Real Field with 53 bits of precision

        TESTS::

            sage: EllipticCurve(QQ, [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x + 6 over Rational Field
            sage: EllipticCurve(RR, [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 6.00000000000000 over Real Field with 53 bits of precision
            sage: EllipticCurve(RealField(100), [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + 1.0000000000000000000000000000*x + 6.0000000000000000000000000000 over Real Field with 100 bits of precision
            sage: EllipticCurve(CC, [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 6.00000000000000 over Complex Field with 53 bits of precision
            sage: EllipticCurve(ComplexField(100), [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + 1.0000000000000000000000000000*x + 6.0000000000000000000000000000 over Complex Field with 100 bits of precision
            sage: EllipticCurve(AA, [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x + 6 over Algebraic Real Field
            sage: EllipticCurve(QQbar, [1, 6]).period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x + 6 over Algebraic Field

        Unsupported cases::

            sage: EllipticCurve(ZZ, [1, 6]).period_lattice()
            Traceback (most recent call last):
            ...
            AttributeError: 'EllipticCurve_generic_with_category' object has no attribute 'period_lattice'
            sage: QQt.<t> = QQ[]
            sage: EllipticCurve(QQt.fraction_field(), [1, 6]).period_lattice()
            Traceback (most recent call last):
            ...
            AttributeError: 'FractionField_1poly_field_with_category' object has no attribute 'embeddings'
            sage: EllipticCurve(GF(7), [1, 6]).period_lattice()
            Traceback (most recent call last):
            ...
            AttributeError: 'FiniteField_prime_modn_with_category' object has no attribute 'embeddings'
        """
        from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
        return PeriodLattice_ell(self)

    def kernel_polynomial_from_point(self, P, *, algorithm=None):
        r"""
        Given a point `P` on this curve which generates a rational subgroup,
        return the kernel polynomial of that subgroup as a polynomial over
        the base field of the curve.
        (The point `P` itself may be defined over an extension.)

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [1,1])
            sage: F = GF(101^3)
            sage: EE = E.change_ring(F)
            sage: xK = F([77, 28, 8]); xK
            8*z3^2 + 28*z3 + 77
            sage: K = EE.lift_x(xK); K.order()
            43
            sage: E.kernel_polynomial_from_point(K)
            x^21 + 7*x^20 + 22*x^19 + 4*x^18 + 7*x^17 + 81*x^16 + 41*x^15 + 68*x^14 + 18*x^13 + 58*x^12 + 31*x^11 + 26*x^10 + 62*x^9 + 20*x^8 + 73*x^7 + 23*x^6 + 66*x^5 + 79*x^4 + 12*x^3 + 40*x^2 + 50*x + 93

        The ``'minpoly'`` algorithm is often much faster than the
        ``'basic'`` algorithm::

            sage: from sage.schemes.elliptic_curves.ell_field import EllipticCurve_field, point_of_order
            sage: p = 2^127 - 1
            sage: E = EllipticCurve(GF(p), [1,0])
            sage: P = point_of_order(E, 31)                                             # long time (8.5s)
            sage: %timeit E.kernel_polynomial_from_point(P, algorithm='basic')          # not tested
            4.38 ms ± 13.7 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
            sage: %timeit E.kernel_polynomial_from_point(P, algorithm='minpoly')        # not tested
            854 µs ± 1.56 µs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)

        Example of finding all the rational isogenies using this method::

            sage: E = EllipticCurve(GF(71), [1,2,3,4,5])
            sage: F = E.division_field(11)
            sage: EE = E.change_ring(F)
            sage: fs = set()
            sage: for K in EE(0).division_points(11):
            ....:     if not K:
            ....:         continue
            ....:     Kp = EE.frobenius_isogeny()(K)
            ....:     if Kp.weil_pairing(K, 11) == 1:
            ....:         fs.add(E.kernel_polynomial_from_point(K))
            sage: fs = sorted(fs); fs
            [x^5 + 10*x^4 + 18*x^3 + 10*x^2 + 43*x + 46,
             x^5 + 65*x^4 + 39*x^2 + 20*x + 63]
            sage: from sage.schemes.elliptic_curves.isogeny_small_degree import is_kernel_polynomial
            sage: {is_kernel_polynomial(E, 11, f) for f in fs}
            {True}
            sage: isogs = [E.isogeny(f) for f in fs]
            sage: isogs[0]
            Isogeny of degree 11 from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 71 to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 34*x + 42 over Finite Field of size 71
            sage: isogs[1]
            Isogeny of degree 11 from Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 71 to Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 12*x + 40 over Finite Field of size 71
            sage: set(isogs) == set(E.isogenies_prime_degree(11))
            True

        ALGORITHM:

        - The ``'basic'`` algorithm is to multiply together all the linear
          factors `(X - x([i]P))` of the kernel polynomial using a product
          tree, then converting the result to the base field of the curve.
          Its complexity is `\widetilde O(\ell k)` where `k` is the
          extension degree.

        - The ``'minpoly'`` algorithm is
          [EPSV2023]_, Algorithm 4 (``KernelPolynomialFromIrrationalX``).
          Over finite fields, its complexity is `O(\ell k) + \widetilde O(\ell)`
          where `k` is the extension degree.
        """
        R = self.base_ring()

        if not P:
            return R['x'].one()

        S = P.base_ring()
        if not S.has_coerce_map_from(R):
            raise TypeError(f'{R} does not coerce into {S}')

        EE = self.change_ring(S)
        if P.curve() is not EE:
            raise TypeError(f'{P} is not a point on {EE}')

        l = P.order()

        if algorithm is None:
            if R in FiniteFields():
                # In this case the minpoly approach is likely to be faster.
                if l & 1 and l.is_prime_power():
                    algorithm = 'minpoly'
            if algorithm is None:
                algorithm = 'basic'

        if algorithm == 'basic':
            from sage.groups.generic import multiples
            Qs = multiples(P, l//2, P)
            x = polygen(S)
            f = prod(x - Q.xy()[0] for Q in Qs)
            return f.change_ring(R)

        if algorithm == 'minpoly':
            if not l & 1 or not l.is_prime_power():
                raise ValueError('algorithm "minpoly" only supports odd prime-power degrees')

            xx = P.xy()[0]
            ext = xx.parent().over(self.base_ring())
            mu = ext(xx).minpoly()
            assert mu.base_ring() == self.base_ring()

            return self.kernel_polynomial_from_divisor(mu, P.order(), check=False)

        raise ValueError('unknown algorithm')

    def kernel_polynomial_from_divisor(self, f, l, *, check=True):
        r"""
        Given an irreducible divisor `f` of the `l`-division polynomial
        on this curve, return the kernel polynomial defining the subgroup
        defined by `f`.

        If the given polynomial does not define a rational subgroup, a
        :exc:`ValueError` is raised.

        This method is currently only implemented for prime `l`.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101^2), [0,1])
            sage: f,_ = E.division_polynomial(5).factor()[0]
            sage: ker = E.kernel_polynomial_from_divisor(f, 5); ker
            x^2 + (49*z2 + 10)*x + 30*z2 + 80
            sage: E.isogeny(ker)
            Isogeny of degree 5
             from Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field in z2 of size 101^2
             to Elliptic Curve defined by y^2 = x^3 + (6*z2+16)*x + 18 over Finite Field in z2 of size 101^2

        The method detects invalid inputs::

            sage: E = EllipticCurve(GF(101), [0,1])
            sage: f,_ = E.division_polynomial(5).factor()[-1]
            sage: E.kernel_polynomial_from_divisor(f, 5)
            Traceback (most recent call last):
            ...
            ValueError: given polynomial does not define a rational 5-isogeny

        ::

            sage: E = EllipticCurve(GF(101), [1,1])
            sage: f,_ = E.division_polynomial(7).factor()[-1]
            sage: E.kernel_polynomial_from_divisor(f, 7)
            Traceback (most recent call last):
            ...
            ValueError: given polynomial does not define a rational 7-isogeny

        ::

            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^12 - 2*x^10 + 3*x^8 + 228/13*x^6 + 235/13*x^4 + 22/13*x^2 + 1/13)
            sage: E = EllipticCurve(K, [1,0])
            sage: ker = E.kernel_polynomial_from_divisor(x - t, 13); ker
            x^6 + (-169/64*t^10 + 169/32*t^8 - 247/32*t^6 - 377/8*t^4 - 2977/64*t^2 - 105/32)*x^4 + (-169/32*t^10 + 169/16*t^8 - 247/16*t^6 - 377/4*t^4 - 2977/32*t^2 - 89/16)*x^2 - 13/64*t^10 + 13/32*t^8 - 19/32*t^6 - 29/8*t^4 - 229/64*t^2 - 13/32
            sage: phi = E.isogeny(ker, check=True); phi
            Isogeny of degree 13
             from Elliptic Curve defined by y^2 = x^3 + x
              over Number Field in t with defining polynomial x^12 - 2*x^10 + 3*x^8 + 228/13*x^6 + 235/13*x^4 + 22/13*x^2 + 1/13
             to Elliptic Curve defined by y^2 = x^3 + (-2535/16*t^10+2535/8*t^8-3705/8*t^6-5655/2*t^4-44655/16*t^2-2047/8)*x
              over Number Field in t with defining polynomial x^12 - 2*x^10 + 3*x^8 + 228/13*x^6 + 235/13*x^4 + 22/13*x^2 + 1/13

        ALGORITHM: [EPSV2023]_, Algorithm 3 (``KernelPolynomialFromDivisor``).
        """
        l = ZZ(l)
        if check:
            if not l.is_prime():
                raise NotImplementedError('currently, kernel_polynomial_from_divisor() only supports prime orders')
            if not f.is_irreducible():
                raise NotImplementedError('currently, kernel_polynomial_from_divisor() only supports irreducible polynomials')
            if f.parent().base_ring() != self.base_ring():
                raise TypeError(f'given polynomial is not defined over the base ring of the curve')
            if self.division_polynomial(l, x=f.parent().quotient_ring(f).gen()):
                raise ValueError(f'given polynomial does not divide the {l}-division polynomial')

        if l == 2:
            return f

        if not f.degree().divides(l//2):
            raise ValueError(f'given polynomial does not define a rational {l}-isogeny')

        from sage.schemes.elliptic_curves.isogeny_small_degree import _least_semi_primitive
        a = _least_semi_primitive(l)
        mul_a = lambda x: self._multiple_x_numerator(a, x=x) / self._multiple_x_denominator(a, x=x)
        x_mod = lambda g: g.parent().quotient(g).gen()

        fs = [f]
        m = l//2//f.degree()

        for i in range(1, m):
            fs.append(mul_a(x_mod(fs[-1])).minpoly())

        if fs[0](mul_a(x_mod(fs[-1]))):
            raise ValueError(f'given polynomial does not define a rational {l}-isogeny')

        return prod(fs)

    def isogenies_prime_degree(self, l=None, max_l=31):
        """
        Return a list of all separable isogenies (up to post-composition with
        isomorphisms) of given prime degree(s) with domain equal to ``self``,
        which are defined over the base field of ``self``.

        INPUT:

        - ``l`` -- a prime or a list of primes

        - ``max_l`` -- (default: 31) a bound on the primes to be tested.
          This is only used if ``l`` is ``None``.

        OUTPUT:

        (list) All separable `l`-isogenies for the given `l` with domain ``self``.

        ALGORITHM:

        Calls the generic function :func:`isogenies_prime_degree()`.
        This is generic code, valid for all fields. It requires that
        certain operations have been implemented over the base field,
        such as root-finding for univariate polynomials.

        EXAMPLES:

        Examples over finite fields::

            sage: # needs sage.libs.pari
            sage: E = EllipticCurve(GF(next_prime(1000000)), [7,8])
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(3)
            []
            sage: E.isogenies_prime_degree(5)
            []
            sage: E.isogenies_prime_degree(7)
            []
            sage: E.isogenies_prime_degree(11)
            []
            sage: E.isogenies_prime_degree(13)
            [Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
             Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(max_l=13)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003,
             Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
             Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree()  # Default limit of 31
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003,
             Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
             Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003,
             Isogeny of degree 17
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 347438*x + 594729 over Finite Field of size 1000003,
             Isogeny of degree 17
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 674846*x + 7392 over Finite Field of size 1000003,
             Isogeny of degree 23
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 390065*x + 605596 over Finite Field of size 1000003]

            sage: E = EllipticCurve(GF(17), [2,0])
            sage: E.isogenies_prime_degree(3)
            []
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17
              to Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 17,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17
              to Elliptic Curve defined by y^2 = x^3 + 5*x + 9 over Finite Field of size 17,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17
              to Elliptic Curve defined by y^2 = x^3 + 5*x + 8 over Finite Field of size 17]

        The base field matters, over a field extension we find more
        isogenies::

            sage: E = EllipticCurve(GF(13), [2,8])
            sage: E.isogenies_prime_degree(max_l=3)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13
                to Elliptic Curve defined by y^2 = x^3 + 7*x + 4 over Finite Field of size 13,
             Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13
                to Elliptic Curve defined by y^2 = x^3 + 9*x + 11 over Finite Field of size 13]

            sage: # needs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(13^6), [2,8])
            sage: E.isogenies_prime_degree(max_l=3)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6
                to Elliptic Curve defined by y^2 = x^3 + 7*x + 4 over Finite Field in z6 of size 13^6,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6
                to Elliptic Curve defined by y^2 = x^3 + (2*z6^5+6*z6^4+9*z6^3+8*z6+7)*x + (3*z6^5+9*z6^4+7*z6^3+12*z6+7) over Finite Field in z6 of size 13^6,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6
                to Elliptic Curve defined by y^2 = x^3 + (11*z6^5+7*z6^4+4*z6^3+5*z6+9)*x + (10*z6^5+4*z6^4+6*z6^3+z6+10) over Finite Field in z6 of size 13^6,
             Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6
                to Elliptic Curve defined by y^2 = x^3 + 9*x + 11 over Finite Field in z6 of size 13^6,
             Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6
                to Elliptic Curve defined by y^2 = x^3 + (3*z6^5+5*z6^4+8*z6^3+11*z6^2+5*z6+12)*x + (12*z6^5+6*z6^4+8*z6^3+4*z6^2+7*z6+6) over Finite Field in z6 of size 13^6,
             Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6
                to Elliptic Curve defined by y^2 = x^3 + (7*z6^4+12*z6^3+7*z6^2+4)*x + (6*z6^5+10*z6^3+12*z6^2+10*z6+8) over Finite Field in z6 of size 13^6,
             Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6
                to Elliptic Curve defined by y^2 = x^3 + (10*z6^5+z6^4+6*z6^3+8*z6^2+8*z6)*x + (8*z6^5+7*z6^4+8*z6^3+10*z6^2+9*z6+7) over Finite Field in z6 of size 13^6]

        If the degree equals the characteristic, we find only separable
        isogenies::

            sage: E = EllipticCurve(GF(13), [2,8])
            sage: E.isogenies_prime_degree(13)
            [Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13
                to Elliptic Curve defined by y^2 = x^3 + 6*x + 5 over Finite Field of size 13]
            sage: E = EllipticCurve(GF(5), [1,1])
            sage: E.isogenies_prime_degree(5)
            [Isogeny of degree 5
              from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 5
                to Elliptic Curve defined by y^2 = x^3 + x + 4 over Finite Field of size 5]

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(3^4)
            sage: E = EllipticCurve(k, [0,1,0,0,a])
            sage: E.isogenies_prime_degree(3)
            [Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + x^2 + a
                   over Finite Field in a of size 3^4
                to Elliptic Curve defined by y^2 = x^3 + x^2 + (2*a^3+a^2+2)*x + (a^2+2)
                   over Finite Field in a of size 3^4]

        In the supersingular case, there are no separable isogenies of
        degree equal to the characteristic::

            sage: E = EllipticCurve(GF(5), [0,1])
            sage: E.isogenies_prime_degree(5)
            []

        An example over a rational function field::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: E = EllipticCurve(K, [1, t^5])
            sage: E.isogenies_prime_degree(5)
            [Isogeny of degree 5
              from Elliptic Curve defined by y^2 = x^3 + x + t^5 over Fraction Field
                   of Univariate Polynomial Ring in t over Finite Field of size 5
                to Elliptic Curve defined by y^2 = x^3 + x + 4*t over Fraction Field
                   of Univariate Polynomial Ring in t over Finite Field of size 5]

        Examples over number fields (other than QQ)::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: QQroot2.<e> = NumberField(x^2 - 2)
            sage: E = EllipticCurve(QQroot2, j=8000)
            sage: E.isogenies_prime_degree()
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000)
                   over Number Field in e with defining polynomial x^2 - 2
                to Elliptic Curve defined by y^2 = x^3 + (-36750)*x + 2401000
                   over Number Field in e with defining polynomial x^2 - 2,
            Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000)
                   over Number Field in e with defining polynomial x^2 - 2
                to Elliptic Curve defined by y^2 = x^3 + (220500*e-257250)*x + (54022500*e-88837000)
                   over Number Field in e with defining polynomial x^2 - 2,
            Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000)
                   over Number Field in e with defining polynomial x^2 - 2
                to Elliptic Curve defined by y^2 = x^3 + (-220500*e-257250)*x + (-54022500*e-88837000)
                   over Number Field in e with defining polynomial x^2 - 2]
            sage: E = EllipticCurve(QQroot2, [1,0,1,4, -6]); E
            Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6)
             over Number Field in e with defining polynomial x^2 - 2
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6)
                   over Number Field in e with defining polynomial x^2 - 2
                to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-36)*x + (-70)
                   over Number Field in e with defining polynomial x^2 - 2]
            sage: E.isogenies_prime_degree(3)
            [Isogeny of degree 3
              from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6)
                   over Number Field in e with defining polynomial x^2 - 2
                to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-1)*x
                   over Number Field in e with defining polynomial x^2 - 2,
             Isogeny of degree 3
              from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6)
                   over Number Field in e with defining polynomial x^2 - 2
                to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-171)*x + (-874)
                   over Number Field in e with defining polynomial x^2 - 2]

        These are not implemented yet::

            sage: E = EllipticCurve(QQbar, [1,18]); E                                   # needs sage.rings.number_field
            Elliptic Curve defined by y^2 = x^3 + x + 18 over Algebraic Field
            sage: E.isogenies_prime_degree()                                            # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            NotImplementedError: This code could be implemented for QQbar, but has not been yet.

            sage: E = EllipticCurve(CC, [1,18]); E
            Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 18.0000000000000
            over Complex Field with 53 bits of precision
            sage: E.isogenies_prime_degree(11)
            Traceback (most recent call last):
            ...
            NotImplementedError: This code could be implemented for general complex fields,
            but has not been yet.

        TESTS::

            sage: E = EllipticCurve(QQ, [1,1])
            sage: E.isogenies_prime_degree([2, 4])
            Traceback (most recent call last):
            ...
            ValueError: 4 is not prime.
            sage: E.isogenies_prime_degree(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not prime.
        """
        F = self.base_ring()
        if isinstance(F, sage.rings.abc.RealField):
            raise NotImplementedError("This code could be implemented for general real fields, but has not been yet.")
        if isinstance(F, sage.rings.abc.ComplexField):
            raise NotImplementedError("This code could be implemented for general complex fields, but has not been yet.")
        if isinstance(F, sage.rings.abc.AlgebraicField):
            raise NotImplementedError("This code could be implemented for QQbar, but has not been yet.")

        if l is None:
            from sage.rings.fast_arith import prime_range
            L = prime_range(max_l + 1)
        else:
            try:
                l = list(l)
            except TypeError:
                L = [ZZ(l)]
            else:
                L = [ZZ(d) for d in l]

        from .isogeny_small_degree import isogenies_prime_degree
        return sum([isogenies_prime_degree(self, d) for d in L], [])

    def isogenies_degree(self, n, *, _intermediate=False):
        r"""
        Return an iterator of all separable isogenies of given degree (up to
        post-composition with isomorphisms) with domain equal to ``self``,
        which are defined over the base field of ``self``.

        ALGORITHM:

        The prime factors `p` of `n` are processed one by one in decreasing
        order, each time "branching" out by taking isogenies of degree `p`.

        INPUT:

        - ``n`` -- integer, or its
          :class:`~sage.structure.factorization.Factorization`.

        - ``_intermediate`` -- (bool, default: False): If set, the isogenies
          from this curve to the curves traversed within the depth-first search
          are returned. This is for internal use only.

        EXAMPLES::

            sage: E = EllipticCurve(GF(11), [1, 1])
            sage: list(E.isogenies_degree(23 * 19))
            [Composite morphism of degree 437 = 23*19:
               From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11
               To:   Elliptic Curve defined by y^2 = x^3 + 8*x + 7 over Finite Field of size 11,
             Composite morphism of degree 437 = 23*19:
               From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11
               To:   Elliptic Curve defined by y^2 = x^3 + 6*x + 2 over Finite Field of size 11,
             Composite morphism of degree 437 = 23*19:
               From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11
               To:   Elliptic Curve defined by y^2 = x^3 + 2*x + 6 over Finite Field of size 11,
             Composite morphism of degree 437 = 23*19:
               From: Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 11
               To:   Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 11]

        ::

            sage: E = EllipticCurve(GF(next_prime(2^32)), j=1728)
            sage: sorted([phi.codomain().j_invariant() for phi in E.isogenies_degree(11 * 17 * 19^2)])
            [1348157279, 1348157279, 1713365879, 1713365879, 3153894341, 3153894341,
             3225140514, 3225140514, 3673460198, 3673460198, 3994312564, 3994312564]
            sage: it = E.isogenies_degree(2^2); it
            <generator object EllipticCurve_field.isogenies_degree at 0x...>
            sage: all(phi.degree() == 2^2 for phi in it)
            True

        We verify that the isogenies outputted are distinct. Note that we do
        not use a ``set`` or any hash-based data structure, as hashing
        isogenies is slow::

            sage: import itertools
            sage: all_distinct = lambda arr: all(x != y for x, y in itertools.combinations(arr, 2))
            sage: K.<z> = GF((19, 2))
            sage: E = EllipticCurve(K, [11*z+5, 14*z+3])
            sage: S = list(E.isogenies_degree(5^2)); len(S), all_distinct(S)
            (3, True)
            sage: S = list(E.isogenies_degree(5^2*11)); len(S), all_distinct(S)
            (6, True)
            sage: S = list(E.isogenies_degree(5^2*11^4)); len(S), all_distinct(S)       # long time (2s)
            (15, True)

        For curves over number fields, the number of distinct isogenies will usually be small::

            sage: E = EllipticCurve(QQ, [0, 1, 0, -2, 0])
            sage: len(list(E.isogenies_degree(2**1)))
            3
            sage: len(list(E.isogenies_degree(2**5)))
            3
            sage: len(list(E.isogenies_degree(2**8)))                                   # long time (8s)
            1

        ::

            sage: pol = PolynomialRing(QQ, 'x')([529, 782, 1])
            sage: L.<a> = NumberField(pol)
            sage: E = EllipticCurve(j=-7072/1127*a + 2016)
            sage: len(list(E.isogenies_degree(2)))
            3
            sage: len(list(E.isogenies_degree(2**5)))
            3

        ::

            sage: pol = PolynomialRing(QQ, 'x')([1, -3, 5, -5, 5, -3, 1])
            sage: L.<a> = NumberField(pol)
            sage: js = hilbert_class_polynomial(-23).roots(L, multiplicities=False)
            sage: E = EllipticCurve(j=choice(js))
            sage: len(list(E.isogenies_degree(2^3)))                                    # long time (9s)
            10

        TESTS::

            sage: E = EllipticCurve(GF(next_prime(2^32)), j=1728)
            sage: list(E.isogenies_degree(2^2, _intermediate=True))
            [Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 4294967311
               Via:  (u,r,s,t) = (1, 0, 0, 0),
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 4294967311
              to Elliptic Curve defined by y^2 = x^3 + 4294967307*x over Finite Field of size 4294967311,
             Composite morphism of degree 4 = 2^2:
               From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 4294967311
               To:   Elliptic Curve defined by y^2 = x^3 + 16*x over Finite Field of size 4294967311,
             Composite morphism of degree 4 = 2^2:
               From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 4294967311
               To:   Elliptic Curve defined by y^2 = x^3 + 4294967267*x + 4294967199 over Finite Field of size 4294967311,
             Composite morphism of degree 4 = 2^2:
               From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 4294967311
               To:   Elliptic Curve defined by y^2 = x^3 + 4294967267*x + 112 over Finite Field of size 4294967311]
            sage: all(isog.domain() is E for isog in _)
            True
            sage: all(isog.domain() is E for isog in E.isogenies_degree(2^5, _intermediate=True))
            True

        The following curve has no degree-`53` isogenies, so the code is quick::

            sage: E = EllipticCurve(GF(103), [3, 5])
            sage: E.isogenies_prime_degree(53)
            []
            sage: list(E.isogenies_degree(product(prime_range(3, 53)) * 53))
            []
        """
        def compute_key(phi):
            """
            Data used in ``hash(phi)`` excluding the expensive `.kernel_polynomial`.
            """
            return (phi.domain(), phi.codomain(), phi.degree(), phi.scaling_factor())

        from sage.schemes.elliptic_curves.weierstrass_morphism import identity_morphism
        from sage.structure.factorization import Factorization

        if not isinstance(n, Factorization):
            n = Integer(n).factor()

        if n.value() == 1:
            yield identity_morphism(self)
            return

        p = n[-1][0]
        seen = {}

        def insert_seen(phi) -> bool:
            key = compute_key(phi)
            if key not in seen:
                seen[key] = [phi]
                return True
            for psi in seen[key]:
                if psi == phi:
                    return False
            seen[key].append(phi)
            return True

        if _intermediate:
            yield identity_morphism(self)

        # isog: self -> E1
        for isog in self.isogenies_prime_degree(p):
            if _intermediate:
                if insert_seen(isog):
                    # self -> E1
                    yield isog

            Eiso = isog.codomain()
            # next_isog : E1 -> E2
            for next_isog in Eiso.isogenies_degree(n / p, _intermediate=_intermediate):
                # psi: self -> E2
                psi = next_isog * isog
                if insert_seen(psi):
                    # self -> E2
                    yield psi

    def is_isogenous(self, other, field=None):
        """
        Return whether or not ``self`` is isogenous to ``other``.

        INPUT:

        - ``other`` -- another elliptic curve

        - ``field`` -- (default: ``None``) currently not implemented. A
          field containing the base fields of the two elliptic curves
          onto which the two curves may be extended to test if they
          are isogenous over this field. By default ``is_isogenous`` will
          not try to find this field unless one of the curves can be
          be extended into the base field of the ``other``, in which case
          it will test over the larger base field.

        OUTPUT: boolean; ``True`` if there is an isogeny from curve ``self`` to
        curve ``other`` defined over ``field``

        METHOD:

        Over general fields this is only implemented in trivial cases.

        EXAMPLES::

            sage: E1 = EllipticCurve(CC, [1,18]); E1
            Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 18.0000000000000
             over Complex Field with 53 bits of precision
            sage: E2 = EllipticCurve(CC, [2,7]); E2
            Elliptic Curve defined by y^2 = x^3 + 2.00000000000000*x + 7.00000000000000
             over Complex Field with 53 bits of precision
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented for isomorphic curves over general fields.

            sage: E1 = EllipticCurve(Frac(PolynomialRing(ZZ,'t')), [2,19]); E1
            Elliptic Curve defined by y^2 = x^3 + 2*x + 19
             over Fraction Field of Univariate Polynomial Ring in t over Integer Ring
            sage: E2 = EllipticCurve(CC, [23,4]); E2
            Elliptic Curve defined by y^2 = x^3 + 23.0000000000000*x + 4.00000000000000
             over Complex Field with 53 bits of precision
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented for isomorphic curves over general fields.
        """
        from .ell_generic import EllipticCurve_generic
        if not isinstance(other, EllipticCurve_generic):
            raise ValueError("Second argument is not an Elliptic Curve.")
        if self.is_isomorphic(other):
            return True
        else:
            raise NotImplementedError("Only implemented for isomorphic curves over general fields.")

    def weierstrass_p(self, prec=20, algorithm=None):
        r"""
        Compute the Weierstrass `\wp`-function of this elliptic curve.

        ALGORITHM: :func:`sage.schemes.elliptic_curves.ell_wp.weierstrass_p`

        INPUT:

        - ``prec`` -- precision

        - ``algorithm`` -- string or ``None`` (default: ``None``);
          a choice of algorithm among ``'pari'``, ``'fast'``, ``'quadratic'``
          or ``None`` to let this function determine the best algorithm to use

        OUTPUT:

        A Laurent series in one variable `z` with coefficients in the
        base field `k` of `E`.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: E.weierstrass_p(prec=10)
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + O(z^10)
            sage: E.weierstrass_p(prec=8)
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + O(z^8)
            sage: Esh = E.short_weierstrass_model()
            sage: Esh.weierstrass_p(prec=8)
            z^-2 + 13392/5*z^2 + 1080432/7*z^4 + 59781888/25*z^6 + O(z^8)
            sage: E.weierstrass_p(prec=20, algorithm='fast')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8
            + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14
            + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
            sage: E.weierstrass_p(prec=20, algorithm='pari')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8
            + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14
            + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
            sage: E.weierstrass_p(prec=20, algorithm='quadratic')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8
            + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14
            + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
        """
        from .ell_wp import weierstrass_p
        return weierstrass_p(self, prec=prec, algorithm=algorithm)

    def hasse_invariant(self):
        r"""
        Return the Hasse invariant of this elliptic curve.

        OUTPUT:

        The Hasse invariant of this elliptic curve, as an element of
        the base field.  This is only defined over fields of positive
        characteristic, and is an element of the field which is zero
        if and only if the curve is supersingular.  Over a field of
        characteristic zero, where the Hasse invariant is undefined,
        a :exc:`ValueError` is raised.

        EXAMPLES::

            sage: E = EllipticCurve([Mod(1,2), Mod(1,2), 0, 0, Mod(1,2)])
            sage: E.hasse_invariant()
            1
            sage: E = EllipticCurve([0, 0, Mod(1,3), Mod(1,3), Mod(1,3)])
            sage: E.hasse_invariant()
            0
            sage: E = EllipticCurve([0, 0, Mod(1,5), 0, Mod(2,5)])
            sage: E.hasse_invariant()
            0
            sage: E = EllipticCurve([0, 0, Mod(1,5), Mod(1,5), Mod(2,5)])
            sage: E.hasse_invariant()
            2

        Some examples over larger fields::

            sage: # needs sage.rings.finite_rings
            sage: EllipticCurve(GF(101), [0,0,0,0,1]).hasse_invariant()
            0
            sage: EllipticCurve(GF(101), [0,0,0,1,1]).hasse_invariant()
            98
            sage: EllipticCurve(GF(103), [0,0,0,0,1]).hasse_invariant()
            20
            sage: EllipticCurve(GF(103), [0,0,0,1,1]).hasse_invariant()
            17
            sage: F.<a> = GF(107^2)
            sage: EllipticCurve(F, [0,0,0,a,1]).hasse_invariant()
            62*a + 75
            sage: EllipticCurve(F, [0,0,0,0,a]).hasse_invariant()
            0

        Over fields of characteristic zero, the Hasse invariant is
        undefined::

            sage: E = EllipticCurve([0,0,0,0,1])
            sage: E.hasse_invariant()
            Traceback (most recent call last):
            ...
            ValueError: Hasse invariant only defined in positive characteristic
        """
        k = self.base_field()
        p = k.characteristic()
        if p == 0:
            raise ValueError('Hasse invariant only defined in positive characteristic')
        elif p == 2:
            return self.a1()
        elif p == 3:
            return self.b2()
        elif p == 5:
            return self.c4()
        elif p == 7:
            return -self.c6()
        else:
            R = k['x']
            x = R.gen()
            E = self.short_weierstrass_model()
            f = (x**3+E.a4()*x+E.a6())**((p-1)//2)
            return f.coefficients(sparse=False)[p-1]

    def isogeny_ell_graph(self, l, directed=True, label_by_j=False):
        """
        Return a graph representing the ``l``-degree ``K``-isogenies between
        ``K``-isomorphism classes of elliptic curves for ``K =
        self.base_field()``.

        INPUT:

        - ``l`` -- prime degree of isogenies

        - ``directed`` -- boolean (default: ``True``); whether to return a
          directed or undirected graph.  In the undirected case, the in-degrees
          and out-degrees of the vertices must be balanced and therefore the
          number of out-edges from the vertices corresponding to j-invariants 0
          and 1728 (if they are part of the graph) are reduced to match the
          number of in-edges.

        - ``label_by_j`` -- boolean (default: ``False``); whether to label
          graph vertices by the j-invariant corresponding to the isomorphism
          class of curves.  If the j-invariant is not unique in the isogeny
          class, append ``*`` to it to indicate a twist.  Otherwise, if
          ``False`` label vertices by the equation of a representative curve.

        OUTPUT: a :class:`Graph` or :class:`DiGraph`

        EXAMPLES:

        Ordinary curve over finite extension field of degree 2::

            sage: # needs sage.graphs sage.rings.finite_rings
            sage: x = polygen(ZZ, 'x')
            sage: E = EllipticCurve(GF(59^2, "i", x^2 + 1), j=5)
            sage: G = E.isogeny_ell_graph(5, directed=False, label_by_j=True); G
            Graph on 20 vertices
            sage: G.vertices(sort=True)
            ['1',
             '12',
             ...
             'i + 55']
            sage: G.edges(sort=True)
            [('1', '28*i + 11', None),
             ('1', '31*i + 11', None),
             ...
             ('8', 'i + 1', None)]

        Supersingular curve over prime field::

            sage: # needs sage.graphs sage.rings.finite_rings
            sage: E = EllipticCurve(GF(419), j=1728)
            sage: G3 = E.isogeny_ell_graph(3, directed=False, label_by_j=True); G3
            Graph on 27 vertices
            sage: G3.vertices(sort=True)
            ['0',
             '0*',
             ...
             '98*']
            sage: G3.edges(sort=True)
            [('0', '0*', None),
             ('0', '13', None),
             ...
             ('48*', '98*', None)]
             sage: G5 = E.isogeny_ell_graph(5, directed=False, label_by_j=True); G5
             Graph on 9 vertices
             sage: G5.vertices(sort=True)
             ['13', '13*', '407', '407*', '52', '62', '62*', '98', '98*']
             sage: G5.edges(sort=True)
             [('13', '52', None),
              ('13', '98', None),
              ...
              ('62*', '98*', None)]

        Supersingular curve over finite extension field of degree 2::

            sage: # needs sage.graphs sage.rings.finite_rings
            sage: K = GF(431^2, "i", x^2 + 1)
            sage: E = EllipticCurve(K, j=0)
            sage: E.is_supersingular()
            True
            sage: G = E.isogeny_ell_graph(2, directed=True, label_by_j=True); G
            Looped multi-digraph on 37 vertices
            sage: G.vertices(sort=True)
            ['0',
             '102',
             ...
             '87*i + 190']
            sage: G.edges(sort=True)
            [('0', '125', None),
             ('0', '125', None),
             ...
             '81*i + 65', None)]
            sage: H = E.isogeny_ell_graph(2, directed=False, label_by_j=True); H
            Looped multi-graph on 37 vertices
            sage: H.vertices(sort=True)
            ['0',
             '102',
             ...
             '87*i + 190']
            sage: H.edges(sort=True)
            [('0', '125', None),
             ('102', '125', None),
             ...
             ('81*i + 65', '87*i + 190', None)]

        Curve over a quadratic number field::

            sage: # needs sage.graphs sage.rings.finite_rings sage.rings.number_field
            sage: K.<e> = NumberField(x^2 - 2)
            sage: E = EllipticCurve(K, [1, 0, 1, 4, -6])
            sage: G2 = E.isogeny_ell_graph(2, directed=False)
            sage: G2.vertices(sort=True)
            ['y^2 + x*y + y = x^3 + (-130*e-356)*x + (-2000*e-2038)',
             'y^2 + x*y + y = x^3 + (-36)*x + (-70)',
             'y^2 + x*y + y = x^3 + (130*e-356)*x + (2000*e-2038)',
             'y^2 + x*y + y = x^3 + 4*x + (-6)']
            sage: G2.edges(sort=True)
            [('y^2 + x*y + y = x^3 + (-130*e-356)*x + (-2000*e-2038)',
              'y^2 + x*y + y = x^3 + (-36)*x + (-70)', None),
             ('y^2 + x*y + y = x^3 + (-36)*x + (-70)',
              'y^2 + x*y + y = x^3 + (130*e-356)*x + (2000*e-2038)', None),
             ('y^2 + x*y + y = x^3 + (-36)*x + (-70)',
              'y^2 + x*y + y = x^3 + 4*x + (-6)', None)]
            sage: G3 = E.isogeny_ell_graph(3, directed=False)
            sage: G3.vertices(sort=True)
            ['y^2 + x*y + y = x^3 + (-1)*x',
             'y^2 + x*y + y = x^3 + (-171)*x + (-874)',
             'y^2 + x*y + y = x^3 + 4*x + (-6)']
            sage: G3.edges(sort=True)
            [('y^2 + x*y + y = x^3 + (-1)*x',
              'y^2 + x*y + y = x^3 + 4*x + (-6)', None),
             ('y^2 + x*y + y = x^3 + (-171)*x + (-874)',
              'y^2 + x*y + y = x^3 + 4*x + (-6)', None)]

        TESTS::

            sage: E = EllipticCurve(GF(11), j=0)
            sage: G0 = E.isogeny_ell_graph(2, directed=False)
            sage: G0.is_directed()
            False
            sage: G1 = E.isogeny_ell_graph(2, directed=True)
            sage: G1.is_directed()
            True
            sage: G2 = E.isogeny_ell_graph(2, label_by_j=False)
            sage: G2.vertices(sort=True)
            ['y^2 = x^3 + 1',
             'y^2 = x^3 + 2',
             'y^2 = x^3 + 5*x',
             'y^2 = x^3 + 7*x']
            sage: G3 = E.isogeny_ell_graph(2, label_by_j=True)
            sage: G3.vertices(sort=True)
            ['0', '0*', '1', '1*']
        """

        from warnings import warn
        from sage.matrix.constructor import Matrix

        # warn users if things are getting big
        if l == 2:
            curve_max = 1000
        if l == 3:
            curve_max = 700
        elif l < 20:
            curve_max = 200
        else:
            curve_max = 50

        Es = [self]  # list of curves in graph
        A = []  # adjacency matrix
        labels = []  # list of vertex labels
        for (i, E) in enumerate(Es):
            if 0 < curve_max and curve_max < len(Es):
                warn('Isogeny graph contains more than '
                        + str(curve_max) + ' curves.')
                curve_max = 0

            r = [0] * len(Es)  # adjacency matrix row
            for C in [I.codomain() for I in E.isogenies_prime_degree(l)]:
                j = next((k for (k, F) in enumerate(Es) if C.is_isomorphic(F)),
                        -1)  # index of curve isomorphic to codomain of isogeny
                if j >= 0:
                    r[j] += 1
                else:
                    Es.append(C)
                    r.append(1)

            # If the graph is undirected, non-symmetric values in the adjacency
            # matrix will result in Sage outputting different graphs depending
            # on the vertex ordering.  Therefore, scale down the non-loop
            # out-edges of vertices corresponding to j-invariants 0 and 1728 so
            # that the same isogeny graphs output are isomorphic as graphs
            # regardless of the starting vertex.
            if not directed and E.j_invariant() in [0, 1728]:
                m = len(E.automorphisms()) / 2  # multiplicity of out-edges
                r = [v if k == i else v / m for (k, v) in enumerate(r)]

            A.append(r)
            if label_by_j:
                s = str(E.j_invariant())
                while s in labels:
                    s += "*"
                labels.append(s)
            else:
                labels.append(E._equation_string())

        A = Matrix([r + [0] * (len(A) - len(r)) for r in A])
        if directed:
            from sage.graphs.digraph import DiGraph as GraphType
        else:
            from sage.graphs.graph import Graph as GraphType

        G = GraphType(A, format='adjacency_matrix',
                      data_structure='static_sparse')
        # inplace relabelling is necessary for static_sparse graphs
        GL = G.relabel(labels, inplace=False)
        return GL

    def endomorphism_ring_is_commutative(self):
        r"""
        Check whether the endomorphism ring of this elliptic curve
        *over its base field* is commutative.

        ALGORITHM: The endomorphism ring is always commutative in
        characteristic zero. Over finite fields, it is commutative
        if and only if the Frobenius endomorphism is not in `\ZZ`.
        All elliptic curves with non-commutative endomorphism ring
        are supersingular. (The converse holds over the algebraic
        closure, but here we consider endomorphisms *over the field
        of definition*.)

        EXAMPLES::

            sage: EllipticCurve(QQ, [1,1]).endomorphism_ring_is_commutative()
            True
            sage: EllipticCurve(QQ, [1,0]).endomorphism_ring_is_commutative()
            True
            sage: EllipticCurve(GF(19), [1,1]).endomorphism_ring_is_commutative()
            True
            sage: EllipticCurve(GF(19^2), [1,1]).endomorphism_ring_is_commutative()
            True
            sage: EllipticCurve(GF(19), [1,0]).endomorphism_ring_is_commutative()
            True
            sage: EllipticCurve(GF(19^2), [1,0]).endomorphism_ring_is_commutative()
            False
            sage: EllipticCurve(GF(19^3), [1,0]).endomorphism_ring_is_commutative()
            True
        """
        k = self.base()
        if k.characteristic() == 0 or self.is_ordinary():
            return True

        if not k.is_finite():
            raise NotImplementedError

        return self.frobenius() not in ZZ


def compute_model(E, name):
    r"""
    Return a model of an elliptic curve ``E`` of the type specified
    in the ``name`` parameter.

    Used as a helper function in
    :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`.

    INPUT:

    - ``E`` -- elliptic curve

    - ``name`` -- string; current options:

      - ``'minimal'``: Return a global minimal model of ``E`` if it
        exists, and a semi-global minimal model otherwise.
        For this choice, ``E`` must be defined over a number field.
        See :meth:`~sage.schemes.elliptic_curves.ell_number_field.EllipticCurve_number_field.global_minimal_model`.

      - ``'short_weierstrass'``: Return a short Weierstrass model of ``E``
        assuming one exists.
        See :meth:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.short_weierstrass_model`.

      - ``'montgomery'``: Return an (untwisted) Montgomery model of ``E``
        assuming one exists over this field.
        See :meth:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.montgomery_model`.

    OUTPUT: an elliptic curve of the specified type isomorphic to `E`

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_field import compute_model
        sage: E = EllipticCurve([12/7, 405/49, 0, -81/8, 135/64])
        sage: compute_model(E, 'minimal')
        Elliptic Curve defined by y^2 = x^3 - x^2 - 7*x + 10 over Rational Field
        sage: compute_model(E, 'short_weierstrass')
        Elliptic Curve defined by y^2 = x^3 - 48114*x + 4035015 over Rational Field
        sage: compute_model(E, 'montgomery')
        Elliptic Curve defined by y^2 = x^3 + 5*x^2 + x over Rational Field
    """
    if not isinstance(E, ell_generic.EllipticCurve_generic):
        raise TypeError('not an elliptic curve')

    if name == 'minimal':
        from sage.rings.number_field.number_field_base import NumberField
        if not isinstance(E.base_field(), NumberField):
            raise ValueError('can only compute minimal model for curves over number fields')
        return E.global_minimal_model(semi_global=True)

    if name == 'short_weierstrass':
        return E.short_weierstrass_model()

    if name == 'montgomery':
        return E.montgomery_model()

    raise NotImplementedError(f'cannot compute {name} model')


def point_of_order(E, n):
    r"""
    Given an elliptic curve `E` over a finite field or a number field
    and an integer `n \geq 1`, construct a point of order `n` on `E`,
    possibly defined over an extension of the base field of `E`.

    Currently only prime powers `n` are supported.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_field import point_of_order
        sage: E = EllipticCurve(GF(101), [1,2,3,4,5])
        sage: P = point_of_order(E, 5); P  # random
        (50*Y^5 + 48*Y^4 + 26*Y^3 + 37*Y^2 + 48*Y + 15 : 25*Y^5 + 31*Y^4 + 79*Y^3 + 39*Y^2 + 3*Y + 20 : 1)
        sage: P.base_ring()
        Finite Field in Y of size 101^6
        sage: P.order()
        5
        sage: P.curve().a_invariants()
        (1, 2, 3, 4, 5)

    ::

        sage: Q = point_of_order(E, 8); Q  # random
        (69*x^5 + 24*x^4 + 100*x^3 + 65*x^2 + 88*x + 97 : 65*x^5 + 28*x^4 + 5*x^3 + 45*x^2 + 42*x + 18 : 1)
        sage: 8*Q == 0 and 4*Q != 0
        True

    ::

        sage: from sage.schemes.elliptic_curves.ell_field import point_of_order
        sage: E = EllipticCurve(QQ, [7,7])
        sage: P = point_of_order(E, 3); P  # random
        (x : -Y : 1)
        sage: P.base_ring()
        Number Field in Y with defining polynomial Y^2 - x^3 - 7*x - 7 over its base field
        sage: P.base_ring().base_field()
        Number Field in x with defining polynomial x^4 + 14*x^2 + 28*x - 49/3
        sage: P.order()
        3
        sage: P.curve().a_invariants()
        (0, 0, 0, 7, 7)

    ::

        sage: Q = point_of_order(E, 4); Q  # random
        (x : Y : 1)
        sage: Q.base_ring()
        Number Field in Y with defining polynomial Y^2 - x^3 - 7*x - 7 over its base field
        sage: Q.base_ring().base_field()
        Number Field in x with defining polynomial x^6 + 35*x^4 + 140*x^3 - 245*x^2 - 196*x - 735
        sage: Q.order()
        4
    """
    # Construct the field extension defined by the given polynomial,
    # in such a way that the result is recognized by Sage as a field.
    def ffext(poly):
        rng = poly.parent()
        fld = rng.base_ring()
        if fld in FiniteFields():
            # Workaround: .extension() would return a PolynomialQuotientRing
            # rather than another FiniteField.
            return poly.splitting_field(rng.variable_name())
        return fld.extension(poly, rng.variable_name())

    n = ZZ(n)
    if n == 1:
        return E(0)

    l,m = n.is_prime_power(get_data=True)
    if not m:
        raise NotImplementedError('only prime-power orders are currently supported')

    xpoly = E.division_polynomial(n).radical()
    xpoly //= E.division_polynomial(n//l).radical()
    if xpoly.degree() < 1:  # supersingular and l == p
        raise ValueError('curve does not have any points of the specified order')

    mu = xpoly.factor()[0][0]
    FF = ffext(mu)
    xx = mu.any_root(ring=FF, assume_squarefree=True)

    Y = polygen(FF, 'Y')
    ypoly = E.defining_polynomial()(xx, Y, 1)
    if ypoly.is_irreducible():
        FF = ffext(ypoly)
        xx = FF(xx)

    EE = E.change_ring(FF)
    pt = EE.lift_x(xx)
    pt.set_order(n, check=False)
    return pt
