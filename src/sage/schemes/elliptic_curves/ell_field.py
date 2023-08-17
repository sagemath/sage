r"""
Elliptic curves over a general field

This module defines the class :class:`EllipticCurve_field`, based on
:class:`EllipticCurve_generic`, for elliptic curves over general fields.
"""
#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.rings.abc
from sage.categories.number_fields import NumberFields
from sage.categories.finite_fields import FiniteFields
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.rational_field import QQ
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_field

from .constructor import EllipticCurve
from .ell_curve_isogeny import EllipticCurveIsogeny, isogeny_codomain_from_kernel
from . import ell_generic

class EllipticCurve_field(ell_generic.EllipticCurve_generic, ProjectivePlaneCurve_field):

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

            sage: E = EllipticCurve(GF(3), [0, -1, 0, -346, 2652])                      # optional - sage.rings.finite_rings
            sage: E.genus()                                                             # optional - sage.rings.finite_rings
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

        - ``D`` (default None) the twisting parameter (see below).

        In characteristics other than 2, `D` must be nonzero, and the
        twist is isomorphic to self after adjoining `\sqrt(D)` to the
        base.

        In characteristic 2, `D` is arbitrary, and the twist is
        isomorphic to self after adjoining a root of `x^2+x+D` to the
        base.

        In characteristic 2 when `j=0`, this is not implemented.

        If the base field `F` is finite, `D` need not be specified,
        and the curve returned is the unique curve (up to isomorphism)
        defined over `F` isomorphic to the original curve over the
        quadratic extension of `F` but not over `F` itself.  Over
        infinite fields, an error is raised if `D` is not given.

        EXAMPLES::

            sage: E = EllipticCurve([GF(1103)(1), 0, 0, 107, 340]); E                   # optional - sage.rings.finite_rings
            Elliptic Curve defined by y^2 + x*y  = x^3 + 107*x + 340
             over Finite Field of size 1103
            sage: F = E.quadratic_twist(-1); F                                          # optional - sage.rings.finite_rings
            Elliptic Curve defined by y^2  = x^3 + 1102*x^2 + 609*x + 300
             over Finite Field of size 1103
            sage: E.is_isomorphic(F)                                                    # optional - sage.rings.finite_rings
            False
            sage: E.is_isomorphic(F, GF(1103^2,'a'))                                    # optional - sage.rings.finite_rings
            True

        A characteristic 2 example::

            sage: E = EllipticCurve(GF(2), [1,0,1,1,1])                                 # optional - sage.rings.finite_rings
            sage: E1 = E.quadratic_twist(1)                                             # optional - sage.rings.finite_rings
            sage: E.is_isomorphic(E1)                                                   # optional - sage.rings.finite_rings
            False
            sage: E.is_isomorphic(E1, GF(4,'a'))                                        # optional - sage.rings.finite_rings
            True

        Over finite fields, the twisting parameter may be omitted::

            sage: k.<a> = GF(2^10)                                                      # optional - sage.rings.finite_rings
            sage: E = EllipticCurve(k, [a^2,a,1,a+1,1])                                 # optional - sage.rings.finite_rings
            sage: Et = E.quadratic_twist()                                              # optional - sage.rings.finite_rings
            sage: Et  # random (only determined up to isomorphism)                      # optional - sage.rings.finite_rings
            Elliptic Curve defined
             by y^2 + x*y  = x^3 + (a^7+a^4+a^3+a^2+a+1)*x^2 + (a^8+a^6+a^4+1)
             over Finite Field in a of size 2^10
            sage: E.is_isomorphic(Et)                                                   # optional - sage.rings.finite_rings
            False
            sage: E.j_invariant() == Et.j_invariant()                                   # optional - sage.rings.finite_rings
            True

            sage: p = next_prime(10^10)                                                 # optional - sage.rings.finite_rings
            sage: k = GF(p)                                                             # optional - sage.rings.finite_rings
            sage: E = EllipticCurve(k, [1,2,3,4,5])                                     # optional - sage.rings.finite_rings
            sage: Et = E.quadratic_twist()                                              # optional - sage.rings.finite_rings
            sage: Et  # random (only determined up to isomorphism)                      # optional - sage.rings.finite_rings
            Elliptic Curve defined
             by y^2  = x^3 + 7860088097*x^2 + 9495240877*x + 3048660957
             over Finite Field of size 10000000019
            sage: E.is_isomorphic(Et)                                                   # optional - sage.rings.finite_rings
            False
            sage: k2 = GF(p^2,'a')                                                      # optional - sage.rings.finite_rings
            sage: E.change_ring(k2).is_isomorphic(Et.change_ring(k2))                   # optional - sage.rings.finite_rings
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
            sage: K.<alpha> = QQ.extension(E.division_polynomial(2).monic())            # optional - sage.rings.number_field
            sage: E.base_extend(K).two_torsion_rank()                                   # optional - sage.rings.number_field
            1
            sage: E.reduction(53).two_torsion_rank()
            2

        ::

            sage: E = EllipticCurve('14a1')
            sage: E.two_torsion_rank()
            1
            sage: f = E.division_polynomial(2).monic().factor()[1][0]
            sage: K.<alpha> = QQ.extension(f)                                           # optional - sage.rings.number_field
            sage: E.base_extend(K).two_torsion_rank()                                   # optional - sage.rings.number_field
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

        - ``D`` (must be nonzero) -- the twisting parameter

        .. NOTE::

            The characteristic must not be 2 or 3, and the `j`-invariant must be 1728.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(1728)); E                             # optional - sage.rings.finite_rings
            Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 13
            sage: E1 = E.quartic_twist(2); E1                                           # optional - sage.rings.finite_rings
            Elliptic Curve defined by y^2 = x^3 + 5*x over Finite Field of size 13
            sage: E.is_isomorphic(E1)                                                   # optional - sage.rings.finite_rings
            False
            sage: E.is_isomorphic(E1, GF(13^2,'a'))                                     # optional - sage.rings.finite_rings
            False
            sage: E.is_isomorphic(E1, GF(13^4,'a'))                                     # optional - sage.rings.finite_rings
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

        - ``D`` (must be nonzero) -- the twisting parameter

        .. NOTE::

            The characteristic must not be 2 or 3, and the `j`-invariant must be 0.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(0)); E                                # optional - sage.rings.finite_rings
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 13
            sage: E1 = E.sextic_twist(2); E1                                            # optional - sage.rings.finite_rings
            Elliptic Curve defined by y^2 = x^3 + 11 over Finite Field of size 13
            sage: E.is_isomorphic(E1)                                                   # optional - sage.rings.finite_rings
            False
            sage: E.is_isomorphic(E1, GF(13^2,'a'))                                     # optional - sage.rings.finite_rings
            False
            sage: E.is_isomorphic(E1, GF(13^4,'a'))                                     # optional - sage.rings.finite_rings
            False
            sage: E.is_isomorphic(E1, GF(13^6,'a'))                                     # optional - sage.rings.finite_rings
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

        - ``other`` -- an elliptic curves with the same base field as self.

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
            sage: type(E1.is_quadratic_twist(E1)) == type(E1.is_quadratic_twist(E2))   #trac 6574
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

            sage: F = GF(101)                                                           # optional - sage.rings.finite_rings
            sage: E1 = EllipticCurve(F, [4,7])                                          # optional - sage.rings.finite_rings
            sage: E2 = E1.quadratic_twist()                                             # optional - sage.rings.finite_rings
            sage: D = E1.is_quadratic_twist(E2); D != 0                                 # optional - sage.rings.finite_rings
            True
            sage: F = GF(101)                                                           # optional - sage.rings.finite_rings
            sage: E1 = EllipticCurve(F, [4,7])                                          # optional - sage.rings.finite_rings
            sage: E2 = E1.quadratic_twist()                                             # optional - sage.rings.finite_rings
            sage: D = E1.is_quadratic_twist(E2)                                         # optional - sage.rings.finite_rings
            sage: E1.quadratic_twist(D).is_isomorphic(E2)                               # optional - sage.rings.finite_rings
            True
            sage: E1.is_isomorphic(E2)                                                  # optional - sage.rings.finite_rings
            False
            sage: F2 = GF(101^2,'a')                                                    # optional - sage.rings.finite_rings
            sage: E1.change_ring(F2).is_isomorphic(E2.change_ring(F2))                  # optional - sage.rings.finite_rings
            True

        A characteristic 3 example::

            sage: F = GF(3^5,'a')                                                       # optional - sage.rings.finite_rings
            sage: E1 = EllipticCurve_from_j(F(1))                                       # optional - sage.rings.finite_rings
            sage: E2 = E1.quadratic_twist(-1)                                           # optional - sage.rings.finite_rings
            sage: D = E1.is_quadratic_twist(E2); D != 0                                 # optional - sage.rings.finite_rings
            True
            sage: E1.quadratic_twist(D).is_isomorphic(E2)                               # optional - sage.rings.finite_rings
            True

        ::

            sage: E1 = EllipticCurve_from_j(F(0))                                       # optional - sage.rings.finite_rings
            sage: E2 = E1.quadratic_twist()                                             # optional - sage.rings.finite_rings
            sage: D = E1.is_quadratic_twist(E2); D                                      # optional - sage.rings.finite_rings
            1
            sage: E1.is_isomorphic(E2)                                                  # optional - sage.rings.finite_rings
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
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

        - ``other`` -- an elliptic curves with the same base field as self.

        OUTPUT:

        Either 0, if the curves are not quartic twists, or `D` if
        ``other`` is ``self.quartic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        .. NOTE::

            Not fully implemented in characteristics 2 or 3.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(1728))                                # optional - sage.rings.finite_rings
            sage: E1 = E.quartic_twist(2)                                               # optional - sage.rings.finite_rings
            sage: D = E.is_quartic_twist(E1); D!=0                                      # optional - sage.rings.finite_rings
            True
            sage: E.quartic_twist(D).is_isomorphic(E1)                                  # optional - sage.rings.finite_rings
            True

        ::

            sage: E = EllipticCurve_from_j(1728)
            sage: E1 = E.quartic_twist(12345)
            sage: D = E.is_quartic_twist(E1); D
            15999120
            sage: (D/12345).is_perfect_power(4)
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
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

        - ``other`` -- an elliptic curves with the same base field as self.

        OUTPUT:

        Either 0, if the curves are not sextic twists, or `D` if
        ``other`` is ``self.sextic_twist(D)`` (up to isomorphism).
        If ``self`` and ``other`` are isomorphic, returns 1.

        .. NOTE::

            Not fully implemented in characteristics 2 or 3.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(0))                                   # optional - sage.rings.finite_rings
            sage: E1 = E.sextic_twist(2)                                                # optional - sage.rings.finite_rings
            sage: D = E.is_sextic_twist(E1); D != 0                                     # optional - sage.rings.finite_rings
            True
            sage: E.sextic_twist(D).is_isomorphic(E1)                                   # optional - sage.rings.finite_rings
            True

        ::

            sage: E = EllipticCurve_from_j(0)
            sage: E1 = E.sextic_twist(12345)
            sage: D = E.is_sextic_twist(E1); D
            575968320
            sage: (D/12345).is_perfect_power(6)
            True
        """
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
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
        Given an elliptic curve self defined over a field `L` and a
        subfield `K` of `L`, return all elliptic curves over `K` which
        are isomorphic over `L` to self.

        INPUT:

        - `K` -- a field which embeds into the base field `L` of self.

        - `f` (optional) -- an embedding of `K` into `L`.  Ignored if
          `K` is `\QQ`.

        OUTPUT:

        A list (possibly empty) of elliptic curves defined over `K`
        which are isomorphic to self over `L`, up to isomorphism over
        `K`.

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

            sage: F.<b> = QuadraticField(23)                                            # optional - sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: G.<a> = F.extension(x^3 + 5)                                          # optional - sage.rings.number_field
            sage: E = EllipticCurve(j=1728*b).change_ring(G)                            # optional - sage.rings.number_field
            sage: EF = E.descend_to(F); EF                                              # optional - sage.rings.number_field
            [Elliptic Curve defined by y^2 = x^3 + (27*b-621)*x + (-1296*b+2484)
              over Number Field in b with defining polynomial x^2 - 23
              with b = 4.795831523312720?]
            sage: all(Ei.change_ring(G).is_isomorphic(E) for Ei in EF)                  # optional - sage.rings.number_field
            True

        ::

            sage: L.<a> = NumberField(x^4 - 7)                                          # optional - sage.rings.number_field
            sage: K.<b> = NumberField(x^2 - 7, embedding=a^2)                           # optional - sage.rings.number_field
            sage: E = EllipticCurve([a^6, 0])                                           # optional - sage.rings.number_field
            sage: EK = E.descend_to(K); EK                                              # optional - sage.rings.number_field
            [Elliptic Curve defined by y^2 = x^3 + b*x over Number Field in b
              with defining polynomial x^2 - 7 with b = a^2,
             Elliptic Curve defined by y^2 = x^3 + 7*b*x over Number Field in b
              with defining polynomial x^2 - 7 with b = a^2]
            sage: all(Ei.change_ring(L).is_isomorphic(E) for Ei in EK)                  # optional - sage.rings.number_field
            True

        ::

            sage: K.<a> = QuadraticField(17)                                            # optional - sage.rings.number_field
            sage: E = EllipticCurve(j=2*a)                                              # optional - sage.rings.number_field
            sage: E.descend_to(QQ)                                                      # optional - sage.rings.number_field
            []

        TESTS:

        Check that :trac:`16456` is fixed::

            sage: K.<a> = NumberField(x^3 - 2)                                          # optional - sage.rings.number_field
            sage: E = EllipticCurve('11a1').quadratic_twist(2)                          # optional - sage.rings.number_field
            sage: EK = E.change_ring(K)                                                 # optional - sage.rings.number_field
            sage: EK2 = EK.change_weierstrass_model((a,a,a,a+1))                        # optional - sage.rings.number_field
            sage: EK2.descend_to(QQ)                                                    # optional - sage.rings.number_field
            [Elliptic Curve defined by y^2 = x^3 + x^2 - 41*x - 199 over Rational Field]

            sage: k.<i> = QuadraticField(-1)                                            # optional - sage.rings.number_field
            sage: E = EllipticCurve(k,[0,0,0,1,0])                                      # optional - sage.rings.number_field
            sage: E.descend_to(QQ)                                                      # optional - sage.rings.number_field
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

    def division_field(self, l, names='t', map=False, **kwds):
        r"""
        Given an elliptic curve over a number field or finite field `F`
        and a prime number `\ell`, construct the `\ell`-division field
        `F(E[\ell])`.

        The `\ell`-division field is the smallest extension of `F` over
        which all `\ell`-torsion points of `E` are defined.

        INPUT:

        - `\ell` -- a prime number (an element of `\ZZ`)
        - ``names`` -- (default: ``'t'``) a variable name for the division field
        - ``map`` -- (default: ``False``) also return an embedding of the
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
            field is large (e.g. when `\ell` is large or when the Galois
            representation is surjective).  The ``simplify`` flag also
            has a big influence on the running time over number fields:
            sometimes ``simplify=False`` is faster, sometimes the default
            ``simplify=True`` is faster.

        EXAMPLES:

        The 2-division field is the same as the splitting field of
        the 2-division polynomial (therefore, it has degree 1, 2, 3 or 6)::

            sage: E = EllipticCurve('15a1')
            sage: K.<b> = E.division_field(2); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial x
            sage: E = EllipticCurve('14a1')
            sage: K.<b> = E.division_field(2); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^2 + 5*x + 92
            sage: E = EllipticCurve('196b1')
            sage: K.<b> = E.division_field(2); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^3 + x^2 - 114*x - 127
            sage: E = EllipticCurve('19a1')
            sage: K.<b> = E.division_field(2); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial
             x^6 + 10*x^5 + 24*x^4 - 212*x^3 + 1364*x^2 + 24072*x + 104292

        For odd primes `\ell`, the division field is either the splitting
        field of the `\ell`-division polynomial, or a quadratic extension
        of it. ::

            sage: E = EllipticCurve('50a1')
            sage: F.<a> = E.division_polynomial(3).splitting_field(simplify_all=True); F    # optional - sage.rings.number_field
            Number Field in a
             with defining polynomial x^6 - 3*x^5 + 4*x^4 - 3*x^3 - 2*x^2 + 3*x + 3
            sage: K.<b> = E.division_field(3, simplify_all=True); K                         # optional - sage.rings.number_field
            Number Field in b
             with defining polynomial x^6 - 3*x^5 + 4*x^4 - 3*x^3 - 2*x^2 + 3*x + 3

        If we take any quadratic twist, the splitting field of the
        3-division polynomial remains the same, but the 3-division field
        becomes a quadratic extension::

            sage: E = E.quadratic_twist(5)  # 50b3                                          # optional - sage.rings.number_field
            sage: F.<a> = E.division_polynomial(3).splitting_field(simplify_all=True); F    # optional - sage.rings.number_field
            Number Field in a
             with defining polynomial x^6 - 3*x^5 + 4*x^4 - 3*x^3 - 2*x^2 + 3*x + 3
            sage: K.<b> = E.division_field(3, simplify_all=True); K                         # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^12 - 3*x^11 + 8*x^10 - 15*x^9
             + 30*x^8 - 63*x^7 + 109*x^6 - 144*x^5 + 150*x^4 - 120*x^3 + 68*x^2 - 24*x + 4

        Try another quadratic twist, this time over a subfield of `F`::

            sage: G.<c>,_,_ = F.subfields(3)[0]                                             # optional - sage.rings.number_field
            sage: E = E.base_extend(G).quadratic_twist(c); E                                # optional - sage.rings.number_field
            Elliptic Curve defined
             by y^2 = x^3 + 5*a0*x^2 + (-200*a0^2)*x + (-42000*a0^2+42000*a0+126000)
             over Number Field in a0 with defining polynomial x^3 - 3*x^2 + 3*x + 9
            sage: K.<b> = E.division_field(3, simplify_all=True); K                         # optional - sage.rings.number_field
            Number Field in b with defining polynomial
             x^12 + 5*x^10 + 40*x^8 + 315*x^6 + 750*x^4 + 675*x^2 + 2025

        Some higher-degree examples::

            sage: E = EllipticCurve('11a1')                                                 # optional - sage.rings.number_field
            sage: K.<b> = E.division_field(2); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial
             x^6 + 2*x^5 - 48*x^4 - 436*x^3 + 1668*x^2 + 28792*x + 73844
            sage: K.<b> = E.division_field(3); K  # long time (3s on sage.math, 2014)       # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^48 ...
            sage: K.<b> = E.division_field(5); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^4 - x^3 + x^2 - x + 1
            sage: E.division_field(5, 'b', simplify=False)                                  # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^4 + x^3 + 11*x^2 + 41*x + 101
            sage: E.base_extend(K).torsion_subgroup()  # long time (2s on sage.math, 2014)  # optional - sage.rings.number_field
            Torsion Subgroup isomorphic to Z/5 + Z/5 associated to the Elliptic Curve
             defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20)
             over Number Field in b with defining polynomial x^4 - x^3 + x^2 - x + 1

            sage: E = EllipticCurve('27a1')                                                 # optional - sage.rings.number_field
            sage: K.<b> = E.division_field(3); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^2 + 3*x + 9
            sage: K.<b> = E.division_field(2); K                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial
             x^6 + 6*x^5 + 24*x^4 - 52*x^3 - 228*x^2 + 744*x + 3844
            sage: K.<b> = E.division_field(2, simplify_all=True); K                         # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^6 - 3*x^5 + 5*x^3 - 3*x + 1
            sage: K.<b> = E.division_field(5); K   # long time (4s on sage.math, 2014)      # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^48 ...
            sage: K.<b> = E.division_field(7); K  # long time (8s on sage.math, 2014)       # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^72 ...

        Over a number field::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)                                              # optional - sage.rings.number_field
            sage: E = EllipticCurve([0,0,0,0,i])                                            # optional - sage.rings.number_field
            sage: L.<b> = E.division_field(2); L                                            # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^4 - x^2 + 1
            sage: L.<b>, phi = E.division_field(2, map=True); phi                           # optional - sage.rings.number_field
            Ring morphism:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Number Field in b with defining polynomial x^4 - x^2 + 1
              Defn: i |--> -b^3
            sage: L.<b>, phi = E.division_field(3, map=True)                                # optional - sage.rings.number_field
            sage: L                                                                         # optional - sage.rings.number_field
            Number Field in b with defining polynomial x^24 - 6*x^22 - 12*x^21
             - 21*x^20 + 216*x^19 + 48*x^18 + 804*x^17 + 1194*x^16 - 13488*x^15
             + 21222*x^14 + 44196*x^13 - 47977*x^12 - 102888*x^11 + 173424*x^10
             - 172308*x^9 + 302046*x^8 + 252864*x^7 - 931182*x^6 + 180300*x^5
             + 879567*x^4 - 415896*x^3 + 1941012*x^2 + 650220*x + 443089
            sage: phi                                                                       # optional - sage.rings.number_field
            Ring morphism:
              From: Number Field in i with defining polynomial x^2 + 1
              To:   Number Field in b with defining polynomial x^24 ...
              Defn: i |--> -215621657062634529/183360797284413355040732*b^23 ...

        Over a finite field::

            sage: E = EllipticCurve(GF(431^2), [1,0])                                       # optional - sage.rings.finite_rings
            sage: E.division_field(5, map=True)                                             # optional - sage.rings.finite_rings
            (Finite Field in t of size 431^4,
             Ring morphism:
               From: Finite Field in z2 of size 431^2
               To:   Finite Field in t of size 431^4
               Defn: z2 |--> 52*t^3 + 222*t^2 + 78*t + 105)

        ::

            sage: E = EllipticCurve(GF(433^2), [1,0])                                       # optional - sage.rings.finite_rings
            sage: K.<v> = E.division_field(7); K                                            # optional - sage.rings.finite_rings
            Finite Field in v of size 433^16

        .. SEEALSO::

            To compute a basis of the `\ell`-torsion once the base field
            has been extended, you may use
            :meth:`sage.schemes.elliptic_curves.ell_number_field.EllipticCurve_number_field.torsion_subgroup`
            or
            :meth:`sage.schemes.elliptic_curves.ell_finite_field.EllipticCurve_finite_field.torsion_basis`.

        TESTS:

        Some random testing::

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
            sage: q = next_prime_power(randrange(1, 10^9))                                  # optional - sage.rings.finite_rings
            sage: F.<a> = GF(q)                                                             # optional - sage.rings.finite_rings
            sage: while True:                                                               # optional - sage.rings.finite_rings
            ....:     try:
            ....:         E = EllipticCurve([F.random_element() for _ in range(5)])
            ....:     except ArithmeticError:
            ....:         continue
            ....:     break
            sage: l = random_prime(8)                                                       # optional - sage.rings.finite_rings
            sage: K = E.division_field(l)                                                   # optional - sage.rings.finite_rings
            sage: n = E.cardinality(extension_degree=K.degree()//F.degree())                # optional - sage.rings.finite_rings
            sage: (l^2 if q%l else 0 + E.is_ordinary()).divides(n)                          # optional - sage.rings.finite_rings
            True
            sage: check(E, l, K)    # long time                                             # optional - sage.rings.finite_rings

        AUTHORS:

        - Jeroen Demeyer (2014-01-06): :trac:`11905`, use
          ``splitting_field`` method, moved from ``gal_reps.py``, make
          it work over number fields.
        - Lorenz Panny (2022): extend to finite fields
        """
        from sage.misc.verbose import verbose
        l = Integer(l)
        if not l.is_prime():
            raise ValueError("l must be a prime number")

        verbose("Adjoining X-coordinates of %s-torsion points" % l)
        F = self.base_ring()
        f = self.division_polynomial(l)
        if l == 2 or f.is_constant():
            # For l = 2, the division field is the splitting field of
            # the division polynomial.
            # If f is a non-zero constant, the l-torsion is trivial:
            # This means the curve must be supersingular and l == p.
            return f.splitting_field(names, map=map, **kwds)

        # Compute splitting field of X-coordinates.
        # The Galois group of the division field is a subgroup of GL(2,l).
        # The Galois group of the X-coordinates is a subgroup of GL(2,l)/{-1,+1}.
        # We need the map to change the elliptic curve invariants to K.
        if F in NumberFields():
            deg_mult = F.degree() * l * (l+1) * (l-1)**2 // 2
            K, F_to_K = f.splitting_field(names, degree_multiple=deg_mult, map=True, **kwds)
        elif F in FiniteFields():
            K, F_to_K = f.splitting_field('u', map=True, **kwds)
        else:
            raise NotImplementedError('only number fields and finite fields are currently supported')

        verbose("Adjoining Y-coordinates of %s-torsion points" % l)

        # THEOREM (Cremona, https://github.com/sagemath/sage/issues/11905#comment:21).
        # Let K be a field, E an elliptic curve over K and p an odd
        # prime number. Assume that K contains all roots of the
        # p-division polynomial of E. Then either K contains all
        # p-torsion points on E, or it does not contain any p-torsion
        # point.
        #
        # PROOF. Let G be the absolute Galois group of K (every element
        # in it fixes all elements of K). For any p-torsion point P
        # over the algebraic closure and any sigma in G, we must have
        # either sigma(P) = P or sigma(P) = -P (since K contains the
        # X-coordinate of P). Now assume that K does not contain all
        # p-torsion points. Then there exists a point P1 and a sigma in
        # G such that sigma(P1) = -P1. Now take a different p-torsion
        # point P2. Since sigma(P2) must be P2 or -P2 and
        # sigma(P1+P2) = sigma(P1)+sigma(P2) = sigma(P1)-P2 must
        # be P1+P2 or -(P1+P2), it follows that sigma(P2) = -sigma(P2).
        # Therefore, K cannot contain any p-torsion point.
        #
        # This implies that it suffices to adjoin the Y-coordinate
        # of just one point.

        # First factor f over F and then compute a root X of f over K.
        g = f.factor()[0][0]
        X = g.map_coefficients(F_to_K).roots(multiplicities=False)[0]

        # Polynomial defining the corresponding Y-coordinate
        curve = self.defining_polynomial().map_coefficients(F_to_K)
        ypol = curve(X, polygen(K), 1)
        L = ypol.splitting_field(names, map=map, **kwds)
        if map:
            L, K_to_L = L
            return L, F_to_K.post_compose(K_to_L)
        else:
            return L

    def isogeny(self, kernel, codomain=None, degree=None, model=None, check=True, algorithm=None):
        r"""
        Return an elliptic-curve isogeny from this elliptic curve.

        The isogeny can be specified in two ways, by passing either a
        polynomial or a set of torsion points.  The methods used are:

        - Vélu's Formulas: Vélu's original formulas for computing
          isogenies.  This algorithm is selected by giving as the
          ``kernel`` parameter a single point, or a list of points,
          generating a finite subgroup.

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

        - Factored Isogenies (see
          :mod:`~sage.schemes.elliptic_curves.hom_composite`):
          Given a list of points which generate a composite-order
          subgroup, decomposes the isogeny into prime-degree steps.
          This can be used to construct isogenies of extremely large,
          smooth degree.
          This algorithm is selected using ``algorithm="factored"``.

        INPUT:

        - ``kernel`` -- a kernel: either a point on this curve, a list of
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

        - ``degree`` -- an integer (default: ``None``).

          - If ``kernel`` is ``None``, then this is the degree of the isogeny
            from this curve to ``codomain``.

          - If ``kernel`` is not ``None``, then this is used to determine
            whether or not to skip a `\gcd` of the given kernel polynomial
            with the two-torsion polynomial of this curve.

        - ``model`` -- a string (default: ``None``).  Supported values
          (cf. :func:`~sage.schemes.elliptic_curves.ell_field.compute_model`):

          - ``"minimal"``: If ``self`` is a curve over the rationals or
            over a number field, then the codomain is a global minimal
            model where this exists.

          - ``"short_weierstrass"``: The codomain is a short Weierstrass curve,
            assuming one exists.

          - ``"montgomery"``: The codomain is an (untwisted) Montgomery
            curve, assuming one exists over this field.

        - ``check`` (default: ``True``) -- check whether the input is valid.
          Setting this to ``False`` can lead to significant speedups.

        - ``algorithm`` -- string (optional). By default (when ``algorithm``
          is omitted), the "traditional" implementation
          :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`
          is used. The other choices are:

          - ``"velusqrt"``: Use
            :class:`~sage.schemes.elliptic_curves.hom_velusqrt.EllipticCurveHom_velusqrt`.

          - ``"factored"``: Use
            :class:`~sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite`
            to decompose the isogeny into prime-degree steps.

          The ``degree`` parameter is not supported when an ``algorithm``
          is specified.

        OUTPUT:

        An isogeny between elliptic curves. This is a morphism of curves.
        (In all cases, the returned object will be an instance of
        :class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom`.)

        EXAMPLES::

            sage: F = GF(2^5, 'alpha'); alpha = F.gen()                                 # optional - sage.rings.finite_rings
            sage: E = EllipticCurve(F, [1,0,1,1,1])                                     # optional - sage.rings.finite_rings
            sage: R.<x> = F[]                                                           # optional - sage.rings.finite_rings
            sage: phi = E.isogeny(x + 1)                                                # optional - sage.rings.finite_rings
            sage: phi.rational_maps()                                                   # optional - sage.rings.finite_rings
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

            sage: E = EllipticCurve(GF(19),[1,1])                                       # optional - sage.rings.finite_rings
            sage: P = E(15,3); Q = E(2,12)                                              # optional - sage.rings.finite_rings
            sage: (P.order(), Q.order())                                                # optional - sage.rings.finite_rings
            (7, 3)
            sage: phi = E.isogeny([P,Q]); phi                                           # optional - sage.rings.finite_rings
            Isogeny of degree 21
             from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19
               to Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19
            sage: phi(E.random_point())  # all points defined over GF(19) are in the kernel     # optional - sage.rings.finite_rings
            (0 : 1 : 0)

        ::

            sage: E = EllipticCurve(GF(2^32 - 5), [170246996, 2036646110])              # optional - sage.rings.finite_rings
            sage: P = E.lift_x(2)                                                       # optional - sage.rings.finite_rings
            sage: E.isogeny(P, algorithm="factored")                                    # optional - sage.rings.finite_rings
            Composite morphism of degree 1073721825 = 3^4*5^2*11*19*43*59:
              From: Elliptic Curve defined by y^2 = x^3 + 170246996*x + 2036646110
                     over Finite Field of size 4294967291
              To:   Elliptic Curve defined by y^2 = x^3 + 272790262*x + 1903695400
                     over Finite Field of size 4294967291

        Not all polynomials define a finite subgroup (:trac:`6384`)::

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])                                # optional - sage.rings.finite_rings
            sage: phi = E.isogeny([14,27,4,1])                                          # optional - sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: the polynomial x^3 + 4*x^2 + 27*x + 14 does not define a finite
            subgroup of Elliptic Curve defined by y^2 + x*y = x^3 + x + 2
            over Finite Field of size 31

        .. SEEALSO::

            - :class:`~sage.schemes.elliptic_curves.hom.EllipticCurveHom`
            - :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`
            - :class:`~sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite`

        TESTS:

        Until the checking of kernel polynomials was implemented in
        :trac:`23222`, the following raised no error but returned an
        invalid morphism.  See also :trac:`11578`::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2 - x - 1)                                      # optional - sage.rings.number_field
            sage: E = EllipticCurve(K, [-13392, -1080432])                              # optional - sage.rings.number_field
            sage: R.<x> = K[]                                                           # optional - sage.rings.number_field
            sage: phi = E.isogeny( (x-564)*(x - 396/5*a + 348/5) )                      # optional - sage.rings.number_field
            Traceback (most recent call last):
            ...
            ValueError: the polynomial x^2 + (-396/5*a - 2472/5)*x + 223344/5*a - 196272/5 does not
            define a finite subgroup of Elliptic Curve defined by y^2 = x^3 + (-13392)*x + (-1080432)
            over Number Field in a with defining polynomial x^2 - x - 1

        We check that the cached order is correctly copied over::

            sage: E = EllipticCurve(GF(2^127-1), [1,2,3,4,5])                           # optional - sage.rings.finite_rings
            sage: E.set_order(170141183460469231746191640949390434666)                  # optional - sage.rings.finite_rings
            sage: phi = E.isogeny(E.lift_x(77347718128277853096420969229987528666))     # optional - sage.rings.finite_rings
            sage: phi.codomain()._order                                                 # optional - sage.rings.finite_rings
            170141183460469231746191640949390434666
        """
        if algorithm is not None and degree is not None:
            raise TypeError('cannot pass "degree" and "algorithm" parameters simultaneously')
        if algorithm == "velusqrt":
            from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
            return EllipticCurveHom_velusqrt(self, kernel, codomain=codomain, model=model)
        if algorithm == "factored":
            from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            return EllipticCurveHom_composite(self, kernel, codomain=codomain, model=model)
        try:
            return EllipticCurveIsogeny(self, kernel, codomain, degree, model, check=check)
        except AttributeError as e:
            raise RuntimeError("Unable to construct isogeny: %s" % e)

    def isogeny_codomain(self, kernel, degree=None):
        r"""
        Return the codomain of the isogeny from ``self`` with given kernel.

        INPUT:

        - ``kernel`` -- Either a list of points in the kernel of the isogeny,
          or a kernel polynomial (specified as either a univariate polynomial or a coefficient list.)

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

            sage: E = EllipticCurve(GF(2^127 - 1), [1,2,3,4,5])                                     # optional - sage.rings.finite_rings
            sage: E.set_order(170141183460469231746191640949390434666)                              # optional - sage.rings.finite_rings
            sage: E2 = E.isogeny_codomain(E.lift_x(77347718128277853096420969229987528666))         # optional - sage.rings.finite_rings
            sage: E2._order                                                                         # optional - sage.rings.finite_rings
            170141183460469231746191640949390434666

        Test deprecation warning for obsolete argument::

            sage: E.isogeny_codomain(E.lift_x(77347718128277853096420969229987528666), degree=11)   # optional - sage.rings.finite_rings
            doctest:warning
            ...
            DeprecationWarning: The "degree" argument to .isogeny_codomain() does nothing and will be removed.
            ...
            Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 20731788786372791581385345584850817122*x + 125200507378516567345719286707201096361
             over Finite Field of size 170141183460469231731687303715884105727
        """
        if degree is not None:
            from sage.misc.superseded import deprecation
            deprecation(33619, 'The "degree" argument to .isogeny_codomain() does nothing and will be removed.')

        E = isogeny_codomain_from_kernel(self, kernel)
        if self.base_field().is_finite():
            E._fetch_cached_order(self)
        return E

    def isogenies_prime_degree(self, l=None, max_l=31):
        """
        Return a list of all separable isogenies of given prime degree(s)
        with domain equal to ``self``, which are defined over the base
        field of ``self``.

        INPUT:

        - ``l`` -- a prime or a list of primes.

        - ``max_l`` -- (default: 31) a bound on the primes to be tested.
          This is only used if ``l`` is None.

        OUTPUT:

        (list) All separable `l`-isogenies for the given `l` with domain self.

        ALGORITHM:

        Calls the generic function :func:`isogenies_prime_degree()`.
        This is generic code, valid for all fields. It requires that
        certain operations have been implemented over the base field,
        such as root-finding for univariate polynomials.

        EXAMPLES:

        Examples over finite fields::

            sage: E = EllipticCurve(GF(next_prime(1000000)), [7,8])                     # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(2)                                           # optional - sage.rings.finite_rings
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003,
             Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(3)                                           # optional - sage.rings.finite_rings
            []
            sage: E.isogenies_prime_degree(5)                                           # optional - sage.rings.finite_rings
            []
            sage: E.isogenies_prime_degree(7)                                           # optional - sage.rings.finite_rings
            []
            sage: E.isogenies_prime_degree(11)                                          # optional - sage.rings.finite_rings
            []
            sage: E.isogenies_prime_degree(13)                                          # optional - sage.rings.finite_rings
            [Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
             Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003
                to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(max_l=13)                                    # optional - sage.rings.finite_rings
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
            sage: E.isogenies_prime_degree()  # Default limit of 31                     # optional - sage.rings.finite_rings
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

            sage: E = EllipticCurve(GF(17), [2,0])                                      # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(3)                                           # optional - sage.rings.finite_rings
            []
            sage: E.isogenies_prime_degree(2)                                           # optional - sage.rings.finite_rings
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

            sage: E = EllipticCurve(GF(13), [2,8])                                      # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(max_l=3)                                     # optional - sage.rings.finite_rings
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13
                to Elliptic Curve defined by y^2 = x^3 + 7*x + 4 over Finite Field of size 13,
             Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13
                to Elliptic Curve defined by y^2 = x^3 + 9*x + 11 over Finite Field of size 13]
            sage: E = EllipticCurve(GF(13^6), [2,8])                                    # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(max_l=3)                                     # optional - sage.rings.finite_rings
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

            sage: E = EllipticCurve(GF(13), [2,8])                                      # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(13)                                          # optional - sage.rings.finite_rings
            [Isogeny of degree 13
              from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13
                to Elliptic Curve defined by y^2 = x^3 + 6*x + 5 over Finite Field of size 13]
            sage: E = EllipticCurve(GF(5), [1,1])                                       # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(5)                                           # optional - sage.rings.finite_rings
            [Isogeny of degree 5
              from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 5
                to Elliptic Curve defined by y^2 = x^3 + x + 4 over Finite Field of size 5]
            sage: k.<a> = GF(3^4)                                                       # optional - sage.rings.finite_rings
            sage: E = EllipticCurve(k, [0,1,0,0,a])                                     # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(3)                                           # optional - sage.rings.finite_rings
            [Isogeny of degree 3
              from Elliptic Curve defined by y^2 = x^3 + x^2 + a
                   over Finite Field in a of size 3^4
                to Elliptic Curve defined by y^2 = x^3 + x^2 + (2*a^3+a^2+2)*x + (a^2+2)
                   over Finite Field in a of size 3^4]

        In the supersingular case, there are no separable isogenies of
        degree equal to the characteristic::

            sage: E = EllipticCurve(GF(5), [0,1])                                       # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(5)                                           # optional - sage.rings.finite_rings
            []

        An example over a rational function field::

            sage: R.<t> = GF(5)[]                                                       # optional - sage.rings.finite_rings
            sage: K = R.fraction_field()                                                # optional - sage.rings.finite_rings
            sage: E = EllipticCurve(K, [1, t^5])                                        # optional - sage.rings.finite_rings
            sage: E.isogenies_prime_degree(5)                                           # optional - sage.rings.finite_rings
            [Isogeny of degree 5
              from Elliptic Curve defined by y^2 = x^3 + x + t^5 over Fraction Field
                   of Univariate Polynomial Ring in t over Finite Field of size 5
                to Elliptic Curve defined by y^2 = x^3 + x + 4*t over Fraction Field
                   of Univariate Polynomial Ring in t over Finite Field of size 5]

        Examples over number fields (other than QQ)::

            sage: x = polygen(ZZ, 'x')
            sage: QQroot2.<e> = NumberField(x^2 - 2)                                    # optional - sage.rings.number_field
            sage: E = EllipticCurve(QQroot2, j=8000)                                    # optional - sage.rings.number_field
            sage: E.isogenies_prime_degree()                                            # optional - sage.rings.number_field
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

            sage: E = EllipticCurve(QQroot2, [1,0,1,4, -6]); E                          # optional - sage.rings.number_field
            Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6)
             over Number Field in e with defining polynomial x^2 - 2
            sage: E.isogenies_prime_degree(2)                                           # optional - sage.rings.number_field
            [Isogeny of degree 2
              from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6)
                   over Number Field in e with defining polynomial x^2 - 2
                to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-36)*x + (-70)
                   over Number Field in e with defining polynomial x^2 - 2]
            sage: E.isogenies_prime_degree(3)                                           # optional - sage.rings.number_field
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

            sage: E = EllipticCurve(QQbar, [1,18]); E                                   # optional - sage.rings.number_field
            Elliptic Curve defined by y^2 = x^3 + x + 18 over Algebraic Field
            sage: E.isogenies_prime_degree()                                            # optional - sage.rings.number_field
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

    def is_isogenous(self, other, field=None):
        """
        Return whether or not self is isogenous to other.

        INPUT:

        - ``other`` -- another elliptic curve.

        - ``field`` (default None) -- Currently not implemented. A
          field containing the base fields of the two elliptic curves
          onto which the two curves may be extended to test if they
          are isogenous over this field. By default ``is_isogenous`` will
          not try to find this field unless one of the curves can be
          be extended into the base field of the other, in which case
          it will test over the larger base field.

        OUTPUT:

        (bool) True if there is an isogeny from curve ``self`` to
        curve ``other`` defined over ``field``.

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
        from .ell_generic import is_EllipticCurve
        if not is_EllipticCurve(other):
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

        - ``algorithm`` -- string or ``None`` (default: ``None``):
          a choice of algorithm among ``"pari"``, ``"fast"``, ``"quadratic"``;
          or ``None`` to let this function determine the best algorithm to use.

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
        a ``ValueError`` is raised.

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

            sage: EllipticCurve(GF(101), [0,0,0,0,1]).hasse_invariant()                 # optional - sage.rings.finite_rings
            0
            sage: EllipticCurve(GF(101), [0,0,0,1,1]).hasse_invariant()                 # optional - sage.rings.finite_rings
            98
            sage: EllipticCurve(GF(103), [0,0,0,0,1]).hasse_invariant()                 # optional - sage.rings.finite_rings
            20
            sage: EllipticCurve(GF(103), [0,0,0,1,1]).hasse_invariant()                 # optional - sage.rings.finite_rings
            17
            sage: F.<a> = GF(107^2)                                                     # optional - sage.rings.finite_rings
            sage: EllipticCurve(F, [0,0,0,a,1]).hasse_invariant()                       # optional - sage.rings.finite_rings
            62*a + 75
            sage: EllipticCurve(F, [0,0,0,0,a]).hasse_invariant()                       # optional - sage.rings.finite_rings
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

        OUTPUT: A :class:`Graph` or :class:`DiGraph`.

        EXAMPLES:

        Ordinary curve over finite extension field of degree 2::

            sage: x = polygen(ZZ, 'x')
            sage: E = EllipticCurve(GF(59^2, "i", x^2 + 1), j=5)                        # optional - sage.rings.finite_rings
            sage: G = E.isogeny_ell_graph(5, directed=False, label_by_j=True)           # optional - sage.graphs sage.rings.finite_rings
            sage: G                                                                     # optional - sage.graphs sage.rings.finite_rings
            Graph on 20 vertices
            sage: G.vertices(sort=True)                                                 # optional - sage.graphs sage.rings.finite_rings
            ['1',
             '12',
             ...
             'i + 55']
            sage: G.edges(sort=True)                                                    # optional - sage.graphs sage.rings.finite_rings
            [('1', '28*i + 11', None),
             ('1', '31*i + 11', None),
             ...
             ('8', 'i + 1', None)]

        Supersingular curve over prime field::

            sage: E = EllipticCurve(GF(419), j=1728)                                    # optional - sage.rings.finite_rings
            sage: G3 = E.isogeny_ell_graph(3, directed=False, label_by_j=True)          # optional - sage.graphs sage.rings.finite_rings
            sage: G3                                                                    # optional - sage.graphs sage.rings.finite_rings
            Graph on 27 vertices
            sage: G3.vertices(sort=True)                                                # optional - sage.graphs sage.rings.finite_rings
            ['0',
             '0*',
             ...
             '98*']
            sage: G3.edges(sort=True)                                                   # optional - sage.graphs sage.rings.finite_rings
            [('0', '0*', None),
             ('0', '13', None),
             ...
             ('48*', '98*', None)]
             sage: G5 = E.isogeny_ell_graph(5, directed=False, label_by_j=True)         # optional - sage.graphs sage.rings.finite_rings
             sage: G5                                                                   # optional - sage.graphs sage.rings.finite_rings
             Graph on 9 vertices
             sage: G5.vertices(sort=True)                                               # optional - sage.graphs sage.rings.finite_rings
             ['13', '13*', '407', '407*', '52', '62', '62*', '98', '98*']
             sage: G5.edges(sort=True)                                                  # optional - sage.graphs sage.rings.finite_rings
             [('13', '52', None),
              ('13', '98', None),
              ...
              ('62*', '98*', None)]

        Supersingular curve over finite extension field of degree 2::

            sage: K = GF(431^2, "i", x^2 + 1)                                           # optional - sage.rings.finite_rings
            sage: E = EllipticCurve(K, j=0)                                             # optional - sage.rings.finite_rings
            sage: E.is_supersingular()                                                  # optional - sage.rings.finite_rings
            True
            sage: G = E.isogeny_ell_graph(2, directed=True, label_by_j=True)            # optional - sage.graphs sage.rings.finite_rings
            sage: G                                                                     # optional - sage.graphs sage.rings.finite_rings
            Looped multi-digraph on 37 vertices
            sage: G.vertices(sort=True)                                                 # optional - sage.graphs sage.rings.finite_rings
            ['0',
             '102',
             ...
             '87*i + 190']
            sage: G.edges(sort=True)                                                    # optional - sage.graphs sage.rings.finite_rings
            [('0', '125', None),
             ('0', '125', None),
             ...
             '81*i + 65', None)]
            sage: H = E.isogeny_ell_graph(2, directed=False, label_by_j=True)           # optional - sage.graphs sage.rings.finite_rings
            sage: H                                                                     # optional - sage.graphs sage.rings.finite_rings
            Looped multi-graph on 37 vertices
            sage: H.vertices(sort=True)                                                 # optional - sage.graphs sage.rings.finite_rings
            ['0',
             '102',
             ...
             '87*i + 190']
            sage: H.edges(sort=True)                                                    # optional - sage.graphs sage.rings.finite_rings
            [('0', '125', None),
             ('102', '125', None),
             ...
             ('81*i + 65', '87*i + 190', None)]

        Curve over a quadratic number field::

            sage: K.<e> = NumberField(x^2 - 2)                                          # optional - sage.rings.number_field
            sage: E = EllipticCurve(K, [1, 0, 1, 4, -6])                                # optional - sage.rings.number_field
            sage: G2 = E.isogeny_ell_graph(2, directed=False)                           # optional - sage.graphs sage.rings.number_field
            sage: G2.vertices(sort=True)                                                # optional - sage.graphs sage.rings.number_field
            ['y^2 + x*y + y = x^3 + (-130*e-356)*x + (-2000*e-2038)',
             'y^2 + x*y + y = x^3 + (-36)*x + (-70)',
             'y^2 + x*y + y = x^3 + (130*e-356)*x + (2000*e-2038)',
             'y^2 + x*y + y = x^3 + 4*x + (-6)']
            sage: G2.edges(sort=True)                                                   # optional - sage.graphs sage.rings.number_field
            [('y^2 + x*y + y = x^3 + (-130*e-356)*x + (-2000*e-2038)',
              'y^2 + x*y + y = x^3 + (-36)*x + (-70)', None),
             ('y^2 + x*y + y = x^3 + (-36)*x + (-70)',
              'y^2 + x*y + y = x^3 + (130*e-356)*x + (2000*e-2038)', None),
             ('y^2 + x*y + y = x^3 + (-36)*x + (-70)',
              'y^2 + x*y + y = x^3 + 4*x + (-6)', None)]
            sage: G3 = E.isogeny_ell_graph(3, directed=False)                           # optional - sage.graphs sage.rings.number_field
            sage: G3.vertices(sort=True)                                                # optional - sage.graphs sage.rings.number_field
            ['y^2 + x*y + y = x^3 + (-1)*x',
             'y^2 + x*y + y = x^3 + (-171)*x + (-874)',
             'y^2 + x*y + y = x^3 + 4*x + (-6)']
            sage: G3.edges(sort=True)                                                   # optional - sage.graphs sage.rings.number_field
            [('y^2 + x*y + y = x^3 + (-1)*x',
              'y^2 + x*y + y = x^3 + 4*x + (-6)', None),
             ('y^2 + x*y + y = x^3 + (-171)*x + (-874)',
              'y^2 + x*y + y = x^3 + 4*x + (-6)', None)]

        TESTS::

            sage: E = EllipticCurve(GF(11), j=0)                                        # optional - sage.rings.finite_rings
            sage: G0 = E.isogeny_ell_graph(2, directed=False)                           # optional - sage.graphs sage.rings.finite_rings
            sage: G0.is_directed()                                                      # optional - sage.graphs sage.rings.finite_rings
            False
            sage: G1 = E.isogeny_ell_graph(2, directed=True)                            # optional - sage.graphs sage.rings.finite_rings
            sage: G1.is_directed()                                                      # optional - sage.graphs sage.rings.finite_rings
            True
            sage: G2 = E.isogeny_ell_graph(2, label_by_j=False)                         # optional - sage.graphs sage.rings.finite_rings
            sage: G2.vertices(sort=True)                                                # optional - sage.graphs sage.rings.finite_rings
            ['y^2 = x^3 + 1',
             'y^2 = x^3 + 2',
             'y^2 = x^3 + 5*x',
             'y^2 = x^3 + 7*x']
            sage: G3 = E.isogeny_ell_graph(2, label_by_j=True)                          # optional - sage.graphs sage.rings.finite_rings
            sage: G3.vertices(sort=True)                                                # optional - sage.graphs sage.rings.finite_rings
            ['0', '0*', '1', '1*']
        """

        from warnings import warn
        from sage.graphs.graph import DiGraph, Graph
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
        G = (DiGraph if directed else Graph)(A, format="adjacency_matrix",
                data_structure="static_sparse")
        # inplace relabelling is necessary for static_sparse graphs
        GL = G.relabel(labels, inplace=False)
        return GL


def compute_model(E, name):
    r"""
    Return a model of an elliptic curve ``E`` of the type specified
    in the ``name`` parameter.

    Used as a helper function in
    :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`.

    INPUT:

    - ``E`` (elliptic curve)

    - ``name`` (string) -- current options:

      - ``"minimal"``: Return a global minimal model of ``E`` if it
        exists, and a semi-global minimal model otherwise.
        For this choice, ``E`` must be defined over a number field.
        See :meth:`~sage.schemes.elliptic_curves.ell_number_field.EllipticCurve_number_field.global_minimal_model`.

      - ``"short_weierstrass"``: Return a short Weierstrass model of ``E``
        assuming one exists.
        See :meth:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.short_weierstrass_model`.

      - ``"montgomery"``: Return an (untwisted) Montgomery model of ``E``
        assuming one exists over this field.
        See :meth:`~sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.montgomery_model`.

    OUTPUT:

    An elliptic curve of the specified type isomorphic to `E`.

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
