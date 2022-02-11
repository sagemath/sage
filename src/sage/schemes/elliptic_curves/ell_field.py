r"""
Elliptic curves over a general field

This module defines the class ``EllipticCurve_field``, based on
``EllipticCurve_generic``, for elliptic curves over general fields.
"""
#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.rings.all as rings
import sage.rings.abc
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

            sage: E = EllipticCurve(GF(3), [0, -1, 0, -346, 2652])
            sage: E.genus()
            1

            sage: R = FractionField(QQ['z'])
            sage: E = EllipticCurve(R, [0, -1, 0, -346, 2652])
            sage: E.genus()
            1
        """
        return rings.ZZ.one()

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

            sage: E = EllipticCurve([GF(1103)(1), 0, 0, 107, 340]); E
            Elliptic Curve defined by y^2 + x*y  = x^3 + 107*x + 340 over Finite Field of size 1103
            sage: F=E.quadratic_twist(-1); F
            Elliptic Curve defined by y^2  = x^3 + 1102*x^2 + 609*x + 300 over Finite Field of size 1103
            sage: E.is_isomorphic(F)
            False
            sage: E.is_isomorphic(F,GF(1103^2,'a'))
            True

        A characteristic 2 example::

            sage: E=EllipticCurve(GF(2),[1,0,1,1,1])
            sage: E1=E.quadratic_twist(1)
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1,GF(4,'a'))
            True

        Over finite fields, the twisting parameter may be omitted::

            sage: k.<a> = GF(2^10)
            sage: E = EllipticCurve(k,[a^2,a,1,a+1,1])
            sage: Et = E.quadratic_twist()
            sage: Et # random (only determined up to isomorphism)
            Elliptic Curve defined by y^2 + x*y  = x^3 + (a^7+a^4+a^3+a^2+a+1)*x^2 + (a^8+a^6+a^4+1) over Finite Field in a of size 2^10
            sage: E.is_isomorphic(Et)
            False
            sage: E.j_invariant()==Et.j_invariant()
            True

            sage: p=next_prime(10^10)
            sage: k = GF(p)
            sage: E = EllipticCurve(k,[1,2,3,4,5])
            sage: Et = E.quadratic_twist()
            sage: Et # random (only determined up to isomorphism)
            Elliptic Curve defined by y^2  = x^3 + 7860088097*x^2 + 9495240877*x + 3048660957 over Finite Field of size 10000000019
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
                x = rings.polygen(K)
                if char == 2:
                    # We find D such that x^2+x+D is irreducible. If the
                    # degree is odd we can take D=1; otherwise it suffices to
                    # consider odd powers of a generator.
                    D = K(1)
                    if K.degree()%2==0:
                        D = K.gen()
                        a = D**2
                        while (x**2 + x + D).roots():
                            D *= a
                else:
                    # We could take a multiplicative generator but
                    # that might be expensive to compute; otherwise
                    # half the elements will do
                    D = K.random_element()
                    while (x**2 - D).roots():
                        D = K.random_element()
            else:
                raise ValueError("twisting parameter D must be specified over infinite fields.")
        else:
            try:
                D=K(D)
            except ValueError:
                raise ValueError("twisting parameter D must be in the base field.")

            if char!=2 and D.is_zero():
                raise ValueError("twisting parameter D must be nonzero when characteristic is not 2")

        if char!=2:
            b2,b4,b6,b8=self.b_invariants()
            # E is isomorphic to  [0,b2,0,8*b4,16*b6]
            return EllipticCurve(K,[0,b2*D,0,8*b4*D**2,16*b6*D**3])

        # now char==2
        if self.j_invariant() !=0: # iff a1!=0
            a1,a2,a3,a4,a6=self.ainvs()
            E0=self.change_weierstrass_model(a1,a3/a1,0,(a1**2*a4+a3**2)/a1**3)
            # which has the form = [1,A2,0,0,A6]
            assert E0.a1()==K(1)
            assert E0.a3()==K(0)
            assert E0.a4()==K(0)
            return EllipticCurve(K,[1,E0.a2()+D,0,0,E0.a6()])
        else:
            raise ValueError("Quadratic twist not implemented in char 2 when j=0")

    def two_torsion_rank(self):
        r"""
        Return the dimension of the 2-torsion subgroup of
        `E(K)`.

        This will be 0, 1 or 2.

        EXAMPLES::

            sage: E=EllipticCurve('11a1')
            sage: E.two_torsion_rank()
            0
            sage: K.<alpha>=QQ.extension(E.division_polynomial(2).monic())
            sage: E.base_extend(K).two_torsion_rank()
            1
            sage: E.reduction(53).two_torsion_rank()
            2

        ::

            sage: E = EllipticCurve('14a1')
            sage: E.two_torsion_rank()
            1
            sage: K.<alpha>=QQ.extension(E.division_polynomial(2).monic().factor()[1][0])
            sage: E.base_extend(K).two_torsion_rank()
            2

        ::

            sage: EllipticCurve('15a1').two_torsion_rank()
            2

        """
        f=self.division_polynomial(rings.Integer(2))
        n=len(f.roots())+1
        return rings.Integer(n).ord(rings.Integer(2))


    def quartic_twist(self, D):
        r"""
        Return the quartic twist of this curve by `D`.

        INPUT:

        - ``D`` (must be nonzero) -- the twisting parameter..

        .. note::

           The characteristic must not be 2 or 3, and the `j`-invariant must be 1728.

        EXAMPLES::

            sage: E=EllipticCurve_from_j(GF(13)(1728)); E
            Elliptic Curve defined by y^2  = x^3 + x over Finite Field of size 13
            sage: E1=E.quartic_twist(2); E1
            Elliptic Curve defined by y^2  = x^3 + 5*x over Finite Field of size 13
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1,GF(13^2,'a'))
            False
            sage: E.is_isomorphic(E1,GF(13^4,'a'))
            True
        """
        K=self.base_ring()
        char=K.characteristic()
        D=K(D)

        if char==2 or char==3:
            raise ValueError("Quartic twist not defined in chars 2,3")

        if self.j_invariant() !=K(1728):
            raise ValueError("Quartic twist not defined when j!=1728")

        if D.is_zero():
            raise ValueError("quartic twist requires a nonzero argument")

        c4,c6=self.c_invariants()
        # E is isomorphic to  [0,0,0,-27*c4,0]
        assert c6==0
        return EllipticCurve(K,[0,0,0,-27*c4*D,0])

    def sextic_twist(self, D):
        r"""
        Return the quartic twist of this curve by `D`.

        INPUT:

        - ``D`` (must be nonzero) -- the twisting parameter..

        .. note::

           The characteristic must not be 2 or 3, and the `j`-invariant must be 0.

        EXAMPLES::

            sage: E=EllipticCurve_from_j(GF(13)(0)); E
            Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 13
            sage: E1=E.sextic_twist(2); E1
            Elliptic Curve defined by y^2 = x^3 + 11 over Finite Field of size 13
            sage: E.is_isomorphic(E1)
            False
            sage: E.is_isomorphic(E1,GF(13^2,'a'))
            False
            sage: E.is_isomorphic(E1,GF(13^4,'a'))
            False
            sage: E.is_isomorphic(E1,GF(13^6,'a'))
            True
        """
        K=self.base_ring()
        char=K.characteristic()
        D=K(D)

        if char==2 or char==3:
            raise ValueError("Sextic twist not defined in chars 2,3")

        if self.j_invariant() !=K(0):
            raise ValueError("Sextic twist not defined when j!=0")

        if D.is_zero():
            raise ValueError("Sextic twist requires a nonzero argument")

        c4,c6=self.c_invariants()
        # E is isomorphic to  [0,0,0,0,-54*c6]
        assert c4==0
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

        .. note::

           Not fully implemented in characteristic 2, or in
           characteristic 3 when both `j`-invariants are 0.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: Et = E.quadratic_twist(-24)
            sage: E.is_quadratic_twist(Et)
            -6

            sage: E1=EllipticCurve([0,0,1,0,0])
            sage: E1.j_invariant()
            0
            sage: E2=EllipticCurve([0,0,0,0,2])
            sage: E1.is_quadratic_twist(E2)
            2
            sage: E1.is_quadratic_twist(E1)
            1
            sage: type(E1.is_quadratic_twist(E1)) == type(E1.is_quadratic_twist(E2))   #trac 6574
            True

        ::

            sage: E1=EllipticCurve([0,0,0,1,0])
            sage: E1.j_invariant()
            1728
            sage: E2=EllipticCurve([0,0,0,2,0])
            sage: E1.is_quadratic_twist(E2)
            0
            sage: E2=EllipticCurve([0,0,0,25,0])
            sage: E1.is_quadratic_twist(E2)
            5

        ::

            sage: F = GF(101)
            sage: E1 = EllipticCurve(F,[4,7])
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2); D!=0
            True
            sage: F = GF(101)
            sage: E1 = EllipticCurve(F,[4,7])
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

            sage: F = GF(3^5,'a')
            sage: E1 = EllipticCurve_from_j(F(1))
            sage: E2 = E1.quadratic_twist(-1)
            sage: D = E1.is_quadratic_twist(E2); D!=0
            True
            sage: E1.quadratic_twist(D).is_isomorphic(E2)
            True

        ::

            sage: E1 = EllipticCurve_from_j(F(0))
            sage: E2 = E1.quadratic_twist()
            sage: D = E1.is_quadratic_twist(E2); D
            1
            sage: E1.is_isomorphic(E2)
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
        j=E.j_invariant()
        if  j != F.j_invariant():
            return zero

        if E.is_isomorphic(F):
            if K is rings.QQ:
                return rings.ZZ(1)
            return K.one()

        char=K.characteristic()

        if char==2:
            raise NotImplementedError("not implemented in characteristic 2")
        elif char==3:
            if j==0:
                raise NotImplementedError("not implemented in characteristic 3 for curves of j-invariant 0")
            D = E.b2()/F.b2()

        else:
            # now char!=2,3:
            c4E,c6E = E.c_invariants()
            c4F,c6F = F.c_invariants()

            if j==0:
                um = c6E/c6F
                x=rings.polygen(K)
                ulist=(x**3-um).roots(multiplicities=False)
                if  not ulist:
                    D = zero
                else:
                    D = ulist[0]
            elif j==1728:
                um=c4E/c4F
                x=rings.polygen(K)
                ulist=(x**2-um).roots(multiplicities=False)
                if not ulist:
                    D = zero
                else:
                    D = ulist[0]
            else:
                D = (c6E*c4F)/(c6F*c4E)

        # Normalization of output:

        if D.is_zero():
            return D

        if K is rings.QQ:
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

        .. note::

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
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
            raise ValueError("arguments are not elliptic curves")
        K = E.base_ring()
        zero = K.zero()
        if not K == F.base_ring():
            return zero
        j=E.j_invariant()
        if  j != F.j_invariant() or j!=K(1728):
            return zero

        if E.is_isomorphic(F):
            return K.one()

        char=K.characteristic()

        if char==2:
            raise NotImplementedError("not implemented in characteristic 2")
        elif char==3:
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

        .. note::

           Not fully implemented in characteristics 2 or 3.

        EXAMPLES::

            sage: E = EllipticCurve_from_j(GF(13)(0))
            sage: E1 = E.sextic_twist(2)
            sage: D = E.is_sextic_twist(E1); D!=0
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
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        E = self
        F = other
        if not is_EllipticCurve(E) or not is_EllipticCurve(F):
            raise ValueError("arguments are not elliptic curves")
        K = E.base_ring()
        zero = K.zero()
        if not K == F.base_ring():
            return zero
        j=E.j_invariant()
        if  j != F.j_invariant() or not j.is_zero():
            return zero

        if E.is_isomorphic(F):
            return K.one()

        char=K.characteristic()

        if char==2:
            raise NotImplementedError("not implemented in characteristic 2")
        elif char==3:
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

            sage: F.<b> = QuadraticField(23)
            sage: G.<a> = F.extension(x^3+5)
            sage: E = EllipticCurve(j=1728*b).change_ring(G)
            sage: EF = E.descend_to(F); EF
            [Elliptic Curve defined by y^2 = x^3 + (27*b-621)*x + (-1296*b+2484) over Number Field in b with defining polynomial x^2 - 23 with b = 4.795831523312720?]
            sage: all(Ei.change_ring(G).is_isomorphic(E) for Ei in EF)
            True

        ::

            sage: L.<a> = NumberField(x^4 - 7)
            sage: K.<b> = NumberField(x^2 - 7, embedding=a^2)
            sage: E = EllipticCurve([a^6,0])
            sage: EK = E.descend_to(K); EK
            [Elliptic Curve defined by y^2 = x^3 + b*x over Number Field in b with defining polynomial x^2 - 7 with b = a^2,
             Elliptic Curve defined by y^2 = x^3 + 7*b*x over Number Field in b with defining polynomial x^2 - 7 with b = a^2]
            sage: all(Ei.change_ring(L).is_isomorphic(E) for Ei in EK)
            True

        ::

            sage: K.<a> = QuadraticField(17)
            sage: E = EllipticCurve(j = 2*a)
            sage: E.descend_to(QQ)
            []

        TESTS:

        Check that :trac:`16456` is fixed::

            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve('11a1').quadratic_twist(2)
            sage: EK = E.change_ring(K)
            sage: EK2 = EK.change_weierstrass_model((a,a,a,a+1))
            sage: EK2.descend_to(QQ)
            [Elliptic Curve defined by y^2 = x^3 + x^2 - 41*x - 199 over Rational Field]

            sage: k.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(k,[0,0,0,1,0])
            sage: E.descend_to(QQ)
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
        from sage.rings.rational_field import QQ
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

    def isogeny(self, kernel, codomain=None, degree=None, model=None, check=True, algorithm=None):
        r"""
        Return an elliptic-curve isogeny from this elliptic curve.

        The isogeny can be determined in two ways, either by a
        polynomial or a set of torsion points.  The methods used are:

        - Velu's Formulas: Velu's original formulas for computing
          isogenies.  This algorithm is selected by giving as the
          ``kernel`` parameter a point or a list of points which
          generate a finite subgroup.

        - Kohel's Formulas: Kohel's original formulas for computing
          isogenies.  This algorithm is selected by giving as the
          ``kernel`` parameter a polynomial (or a coefficient list
          (little endian)) which will define the kernel of the
          isogeny.

        - Factored Isogenies (*experimental* --- see
          :mod:`sage.schemes.elliptic_curves.hom_composite`):
          Given a list of points which generate a composite-order
          subgroup, decomposes the isogeny into prime-degree steps.
          This can be used to construct isogenies of extremely large,
          smooth degree.
          This algorithm is selected using ``algorithm="factored"``.

        INPUT:

        - ``E``         - an elliptic curve, the domain of the isogeny to
                          initialize.

        - ``kernel`` - a kernel, either a point in ``E``, a list of
                          points in ``E``, a univariate kernel
                          polynomial or ``None``.  If initiating from
                          a domain/codomain, this must be set to None.
                          Validity of input is checked (unless
                          check=False).

        - ``codomain``  - an elliptic curve (default:None).  If ``kernel`` is
                          None, then this must be the codomain of a separable
                          normalized isogeny, furthermore, ``degree`` must be
                          the degree of the isogeny from ``E`` to ``codomain``.
                          If ``kernel`` is not None, then this must be
                          isomorphic to the codomain of the normalized separable
                          isogeny defined by ``kernel``, in this case, the
                          isogeny is post composed with an isomorphism so that
                          this parameter is the codomain.

        - ``degree``    - an integer (default:None). If ``kernel`` is None,
                          then this is the degree of the isogeny from ``E`` to
                          ``codomain``. If ``kernel`` is not None, then this is
                          used to determine whether or not to skip a gcd of the
                          kernel polynomial with the two torsion polynomial of
                          ``E``.

        - ``model`` - a string (default:None).  Only supported
                          variable is "minimal", in which case if``E``
                          is a curve over the rationals or over a
                          number field, then the codomain is a global
                          minimum model where this exists.

        - ``check`` (default: True) checks that the input is valid,
                          i.e., that the polynomial provided is a
                          kernel polynomial, meaning that its roots
                          are the x-coordinates of a finite subgroup.

        - ``algorithm`` (optional): When ``algorithm="factored"`` is
          passed, decompose the isogeny into prime-degree steps.
          The ``degree`` and ``model`` parameters are not supported by
          ``algorithm="factored"``.

        OUTPUT:

        An isogeny between elliptic curves. This is a morphism of curves.

        EXAMPLES::

            sage: F = GF(2^5, 'alpha'); alpha = F.gen()
            sage: E = EllipticCurve(F, [1,0,1,1,1])
            sage: R.<x> = F[]
            sage: phi = E.isogeny(x+1)
            sage: phi.rational_maps()
            ((x^2 + x + 1)/(x + 1), (x^2*y + x)/(x^2 + 1))

        ::

            sage: E = EllipticCurve('11a1')
            sage: P = E.torsion_points()[1]
            sage: E.isogeny(P)
            Isogeny of degree 5 from Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field

        ::

            sage: E = EllipticCurve(GF(19),[1,1])
            sage: P = E(15,3); Q = E(2,12);
            sage: (P.order(), Q.order())
            (7, 3)
            sage: phi = E.isogeny([P,Q]); phi
            Isogeny of degree 21 from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19 to Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 19
            sage: phi(E.random_point()) # all points defined over GF(19) are in the kernel
            (0 : 1 : 0)

        ::

            sage: E = EllipticCurve(GF(2^32-5), [170246996, 2036646110])
            sage: P = E.lift_x(2)
            sage: E.isogeny(P, algorithm="factored")    # experimental
            doctest:warning
            ...
            Composite morphism of degree 1073721825 = 3^4*5^2*11*19*43*59:
              From: Elliptic Curve defined by y^2 = x^3 + 170246996*x + 2036646110 over Finite Field of size 4294967291
              To:   Elliptic Curve defined by y^2 = x^3 + 272790262*x + 1903695400 over Finite Field of size 4294967291

        Not all polynomials define a finite subgroup (:trac:`6384`)::

            sage: E = EllipticCurve(GF(31),[1,0,0,1,2])
            sage: phi = E.isogeny([14,27,4,1])
            Traceback (most recent call last):
            ...
            ValueError: The polynomial x^3 + 4*x^2 + 27*x + 14 does not define a finite subgroup of Elliptic Curve defined by y^2 + x*y = x^3 + x + 2 over Finite Field of size 31.

        TESTS:

        Until the checking of kernel polynomials was implemented in
        :trac:`23222`, the following raised no error but returned an
        invalid morphism.  See also :trac:`11578`::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^2-x-1)
            sage: E = EllipticCurve(K, [-13392, -1080432])
            sage: R.<x> = K[]
            sage: phi = E.isogeny( (x-564)*(x - 396/5*a + 348/5) )
            Traceback (most recent call last):
            ...
            ValueError: The polynomial x^2 + (-396/5*a - 2472/5)*x + 223344/5*a - 196272/5 does not define a finite subgroup of Elliptic Curve defined by y^2 = x^3 + (-13392)*x + (-1080432) over Number Field in a with defining polynomial x^2 - x - 1.
        """
        if algorithm == "factored":
            if degree is not None:
                raise TypeError('algorithm="factored" does not support the "degree" parameter')
            if model  is not None:
                raise TypeError('algorithm="factored" does not support the "model" parameter')
            from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            return EllipticCurveHom_composite(self, kernel, codomain=codomain)
        try:
            return EllipticCurveIsogeny(self, kernel, codomain, degree, model, check=check)
        except AttributeError as e:
            raise RuntimeError("Unable to construct isogeny: %s" % e)


    def isogeny_codomain(self, kernel, degree=None):
        r"""
        Return the codomain of the isogeny from self with given
        kernel.

        INPUT:

        - ``kernel`` - Either a list of points in the kernel of the isogeny,
                       or a kernel polynomial (specified as a either a
                       univariate polynomial or a coefficient list.)

        - ``degree`` - an integer, (default:None) optionally specified degree
                       of the kernel.

        OUTPUT:

        An elliptic curve, the codomain of the separable normalized
        isogeny from this kernel

        EXAMPLES::

            sage: E = EllipticCurve('17a1')
            sage: R.<x> = QQ[]
            sage: E2 = E.isogeny_codomain(x - 11/4); E2
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 1461/16*x - 19681/64 over Rational Field

        """
        return isogeny_codomain_from_kernel(self, kernel, degree=None)

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

            sage: E = EllipticCurve(GF(next_prime(1000000)), [7,8])
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(3)
            []
            sage: E.isogenies_prime_degree(5)
            []
            sage: E.isogenies_prime_degree(7)
            []
            sage: E.isogenies_prime_degree(11)
            []
            sage: E.isogenies_prime_degree(13)
            [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
            Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree(max_l=13)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003,
             Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
             Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003]
            sage: E.isogenies_prime_degree()  # Default limit of 31
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 970389*x + 794257 over Finite Field of size 1000003,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 29783*x + 206196 over Finite Field of size 1000003,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 999960*x + 78 over Finite Field of size 1000003,
             Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 878063*x + 845666 over Finite Field of size 1000003,
             Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 375648*x + 342776 over Finite Field of size 1000003,
             Isogeny of degree 17 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 347438*x + 594729 over Finite Field of size 1000003,
             Isogeny of degree 17 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 674846*x + 7392 over Finite Field of size 1000003,
             Isogeny of degree 23 from Elliptic Curve defined by y^2 = x^3 + 7*x + 8 over Finite Field of size 1000003 to Elliptic Curve defined by y^2 = x^3 + 390065*x + 605596 over Finite Field of size 1000003]

            sage: E = EllipticCurve(GF(17), [2,0])
            sage: E.isogenies_prime_degree(3)
            []
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17 to Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 17,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17 to Elliptic Curve defined by y^2 = x^3 + 5*x + 9 over Finite Field of size 17,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x over Finite Field of size 17 to Elliptic Curve defined by y^2 = x^3 + 5*x + 8 over Finite Field of size 17]

        The base field matters, over a field extension we find more
        isogenies::

            sage: E = EllipticCurve(GF(13), [2,8])
            sage: E.isogenies_prime_degree(max_l=3)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13 to Elliptic Curve defined by y^2 = x^3 + 7*x + 4 over Finite Field of size 13,
             Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13 to Elliptic Curve defined by y^2 = x^3 + 9*x + 11 over Finite Field of size 13]
            sage: E = EllipticCurve(GF(13^6), [2,8])
            sage: E.isogenies_prime_degree(max_l=3)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6 to Elliptic Curve defined by y^2 = x^3 + 7*x + 4 over Finite Field in z6 of size 13^6,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6 to Elliptic Curve defined by y^2 = x^3 + (2*z6^5+6*z6^4+9*z6^3+8*z6+7)*x + (3*z6^5+9*z6^4+7*z6^3+12*z6+7) over Finite Field in z6 of size 13^6,
             Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6 to Elliptic Curve defined by y^2 = x^3 + (11*z6^5+7*z6^4+4*z6^3+5*z6+9)*x + (10*z6^5+4*z6^4+6*z6^3+z6+10) over Finite Field in z6 of size 13^6,
             Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6 to Elliptic Curve defined by y^2 = x^3 + 9*x + 11 over Finite Field in z6 of size 13^6,
             Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6 to Elliptic Curve defined by y^2 = x^3 + (3*z6^5+5*z6^4+8*z6^3+11*z6^2+5*z6+12)*x + (12*z6^5+6*z6^4+8*z6^3+4*z6^2+7*z6+6) over Finite Field in z6 of size 13^6,
             Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6 to Elliptic Curve defined by y^2 = x^3 + (7*z6^4+12*z6^3+7*z6^2+4)*x + (6*z6^5+10*z6^3+12*z6^2+10*z6+8) over Finite Field in z6 of size 13^6,
             Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field in z6 of size 13^6 to Elliptic Curve defined by y^2 = x^3 + (10*z6^5+z6^4+6*z6^3+8*z6^2+8*z6)*x + (8*z6^5+7*z6^4+8*z6^3+10*z6^2+9*z6+7) over Finite Field in z6 of size 13^6]

        If the degree equals the characteristic, we find only separable
        isogenies::

            sage: E = EllipticCurve(GF(13), [2,8])
            sage: E.isogenies_prime_degree(13)
            [Isogeny of degree 13 from Elliptic Curve defined by y^2 = x^3 + 2*x + 8 over Finite Field of size 13 to Elliptic Curve defined by y^2 = x^3 + 6*x + 5 over Finite Field of size 13]
            sage: E = EllipticCurve(GF(5), [1,1])
            sage: E.isogenies_prime_degree(5)
            [Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 5 to Elliptic Curve defined by y^2 = x^3 + x + 4 over Finite Field of size 5]
            sage: k.<a> = GF(3^4)
            sage: E = EllipticCurve(k, [0,1,0,0,a])
            sage: E.isogenies_prime_degree(3)
            [Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + x^2 + a over Finite Field in a of size 3^4 to Elliptic Curve defined by y^2 = x^3 + x^2 + (2*a^3+a^2+2)*x + (a^2+2) over Finite Field in a of size 3^4]

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
            [Isogeny of degree 5 from Elliptic Curve defined by y^2 = x^3 + x + t^5 over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5 to Elliptic Curve defined by y^2 = x^3 + x + 4*t over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5]

        Examples over number fields (other than QQ)::

            sage: QQroot2.<e> = NumberField(x^2-2)
            sage: E = EllipticCurve(QQroot2, j=8000)
            sage: E.isogenies_prime_degree()
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 = x^3 + (-36750)*x + 2401000 over Number Field in e with defining polynomial x^2 - 2,
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 = x^3 + (220500*e-257250)*x + (54022500*e-88837000) over Number Field in e with defining polynomial x^2 - 2,
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + (-150528000)*x + (-629407744000) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 = x^3 + (-220500*e-257250)*x + (-54022500*e-88837000) over Number Field in e with defining polynomial x^2 - 2]

            sage: E = EllipticCurve(QQroot2, [1,0,1,4, -6]); E
            Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2
            sage: E.isogenies_prime_degree(2)
            [Isogeny of degree 2 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-36)*x + (-70) over Number Field in e with defining polynomial x^2 - 2]
            sage: E.isogenies_prime_degree(3)
            [Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-1)*x over Number Field in e with defining polynomial x^2 - 2,
            Isogeny of degree 3 from Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x + (-6) over Number Field in e with defining polynomial x^2 - 2 to Elliptic Curve defined by y^2 + x*y + y = x^3 + (-171)*x + (-874) over Number Field in e with defining polynomial x^2 - 2]

        These are not implemented yet::

            sage: E = EllipticCurve(QQbar, [1,18]); E
            Elliptic Curve defined by y^2 = x^3 + x + 18 over Algebraic Field
            sage: E.isogenies_prime_degree()
            Traceback (most recent call last):
            ...
            NotImplementedError: This code could be implemented for QQbar, but has not been yet.

            sage: E = EllipticCurve(CC, [1,18]); E
            Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 18.0000000000000 over Complex Field with 53 bits of precision
            sage: E.isogenies_prime_degree(11)
            Traceback (most recent call last):
            ...
            NotImplementedError: This code could be implemented for general complex fields, but has not been yet.

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
        if F is rings.QQbar:
            raise NotImplementedError("This code could be implemented for QQbar, but has not been yet.")

        if l is None:
            from sage.rings.all import prime_range
            L = prime_range(max_l + 1)
        else:
            try:
                l = list(l)
            except TypeError:
                L = [rings.ZZ(l)]
            else:
                L = [rings.ZZ(d) for d in l]

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
          are isogenous over this field. By default is_isogenous will
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
            Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 18.0000000000000 over Complex Field with 53 bits of precision
            sage: E2 = EllipticCurve(CC, [2,7]); E2
            Elliptic Curve defined by y^2 = x^3 + 2.00000000000000*x + 7.00000000000000 over Complex Field with 53 bits of precision
            sage: E1.is_isogenous(E2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented for isomorphic curves over general fields.

            sage: E1 = EllipticCurve(Frac(PolynomialRing(ZZ,'t')), [2,19]); E1
            Elliptic Curve defined by y^2 = x^3 + 2*x + 19 over Fraction Field of Univariate Polynomial Ring in t over Integer Ring
            sage: E2 = EllipticCurve(CC, [23,4]); E2
            Elliptic Curve defined by y^2 = x^3 + 23.0000000000000*x + 4.00000000000000 over Complex Field with 53 bits of precision
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
        Computes the Weierstrass `\wp`-function of the elliptic curve.

        INPUT:

        - ``mprec`` - precision

        - ``algorithm`` - string (default:``None``) an algorithm identifier
                      indicating using the ``pari``, ``fast`` or ``quadratic``
                      algorithm. If the algorithm is ``None``, then this
                      function determines the best algorithm to use.

        OUTPUT:

        a Laurent series in one variable `z` with coefficients in the
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
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
            sage: E.weierstrass_p(prec=20, algorithm='pari')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
            sage: E.weierstrass_p(prec=20, algorithm='quadratic')
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)
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
        a ``ValueError`` is returned.

        EXAMPLES::

            sage: E = EllipticCurve([Mod(1,2),Mod(1,2),0,0,Mod(1,2)])
            sage: E.hasse_invariant()
            1
            sage: E = EllipticCurve([0,0,Mod(1,3),Mod(1,3),Mod(1,3)])
            sage: E.hasse_invariant()
            0
            sage: E = EllipticCurve([0,0,Mod(1,5),0,Mod(2,5)])
            sage: E.hasse_invariant()
            0
            sage: E = EllipticCurve([0,0,Mod(1,5),Mod(1,5),Mod(2,5)])
            sage: E.hasse_invariant()
            2

        Some examples over larger fields::

            sage: EllipticCurve(GF(101),[0,0,0,0,1]).hasse_invariant()
            0
            sage: EllipticCurve(GF(101),[0,0,0,1,1]).hasse_invariant()
            98
            sage: EllipticCurve(GF(103),[0,0,0,0,1]).hasse_invariant()
            20
            sage: EllipticCurve(GF(103),[0,0,0,1,1]).hasse_invariant()
            17
            sage: F.<a> = GF(107^2)
            sage: EllipticCurve(F,[0,0,0,a,1]).hasse_invariant()
            62*a + 75
            sage: EllipticCurve(F,[0,0,0,0,a]).hasse_invariant()
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
            f=(x**3+E.a4()*x+E.a6())**((p-1)//2)
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

        OUTPUT:

        A ``Graph`` or ``Digraph``

        EXAMPLES:

        Ordinary curve over finite extension field of degree 2::

            sage: E = EllipticCurve(GF(59^2, "i", x^2 + 1), j=5)
            sage: G = E.isogeny_ell_graph(5, directed=False, label_by_j=True)
            sage: G
            Graph on 20 vertices
            sage: G.vertices()
            ['1',
             '12',
             ...
             'i + 55']
            sage: G.edges()
            [('1', '28*i + 11', None),
             ('1', '31*i + 11', None),
             ...
             ('8', 'i + 1', None)]

        Supersingular curve over prime field::

            sage: E = EllipticCurve(GF(419), j=1728)
            sage: G3 = E.isogeny_ell_graph(3, directed=False, label_by_j=True)
            sage: G3
            Graph on 27 vertices
            sage: G3.vertices()
            ['0',
             '0*',
             ...
             '98*']
            sage: G3.edges()
            [('0', '0*', None),
             ('0', '13', None),
             ...
             ('48*', '98*', None)]
             sage: G5 = E.isogeny_ell_graph(5, directed=False, label_by_j=True)
             sage: G5
             Graph on 9 vertices
             sage: G5.vertices()
             ['13', '13*', '407', '407*', '52', '62', '62*', '98', '98*']
             sage: G5.edges()
             [('13', '52', None),
              ('13', '98', None),
              ...
              ('62*', '98*', None)]

        Supersingular curve over finite extension field of degree 2::

            sage: K = GF(431^2, "i", x^2 + 1)
            sage: E = EllipticCurve(K, j=0)
            sage: E.is_supersingular()
            True
            sage: G = E.isogeny_ell_graph(2, directed=True, label_by_j=True)
            sage: G
            Looped multi-digraph on 37 vertices
            sage: G.vertices()
            ['0',
             '102',
             ...
             '87*i + 190']
            sage: G.edges()
            [('0', '125', None),
             ('0', '125', None),
             ...
             '81*i + 65', None)]
            sage: H = E.isogeny_ell_graph(2, directed=False, label_by_j=True)
            sage: H
            Looped multi-graph on 37 vertices
            sage: H.vertices()
            ['0',
             '102',
             ...
             '87*i + 190']
            sage: H.edges()
            [('0', '125', None),
             ('102', '125', None),
             ...
             ('81*i + 65', '87*i + 190', None)]

        Curve over a quadratic number field::

            sage: K.<e> = NumberField(x^2 - 2)
            sage: E = EllipticCurve(K, [1,0,1,4, -6])
            sage: G2 = E.isogeny_ell_graph(2, directed=False)
            sage: G2.vertices()
            ['y^2 + x*y + y = x^3 + (-130*e-356)*x + (-2000*e-2038)',
             'y^2 + x*y + y = x^3 + (-36)*x + (-70)',
             'y^2 + x*y + y = x^3 + (130*e-356)*x + (2000*e-2038)',
             'y^2 + x*y + y = x^3 + 4*x + (-6)']
            sage: G2.edges()
            [('y^2 + x*y + y = x^3 + (-130*e-356)*x + (-2000*e-2038)',
             'y^2 + x*y + y = x^3 + (-36)*x + (-70)', None),
             ('y^2 + x*y + y = x^3 + (-36)*x + (-70)',
             'y^2 + x*y + y = x^3 + (130*e-356)*x + (2000*e-2038)', None),
             ('y^2 + x*y + y = x^3 + (-36)*x + (-70)',
             'y^2 + x*y + y = x^3 + 4*x + (-6)', None)]
            sage: G3 = E.isogeny_ell_graph(3, directed=False)
            sage: G3.vertices()
            ['y^2 + x*y + y = x^3 + (-1)*x',
             'y^2 + x*y + y = x^3 + (-171)*x + (-874)',
             'y^2 + x*y + y = x^3 + 4*x + (-6)']
            sage: G3.edges()
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
            sage: G2.vertices()
            ['y^2 = x^3 + 1',
             'y^2 = x^3 + 2',
             'y^2 = x^3 + 5*x',
             'y^2 = x^3 + 7*x']
            sage: G3 = E.isogeny_ell_graph(2, label_by_j=True)
            sage: G3.vertices()
            ['0', '0*', '1', '1*']

        """

        from warnings import warn
        from sage.graphs.graph import DiGraph, Graph
        from sage.matrix.all import Matrix

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
