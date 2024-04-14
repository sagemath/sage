# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.finite_rings
"""
Specific algorithms to compute cardinality of elliptic curves over a finite field

Since point counting now uses PARI/GP, this code is only used when a
specific algorithm was specified or when the j-invariant lies in a subfield.

AUTHORS:

- John Cremona (2008-2009): Original point counting code

- Jeroen Demeyer (2017-2018): Refactored and moved to ``cardinality.py``.
"""
# ****************************************************************************
#       Copyright (C) 2008-2009 John Cremona
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from .constructor import EllipticCurve, EllipticCurve_from_j
from sage.schemes.curves.projective_curve import Hasse_bounds
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.polynomial.polynomial_ring import polygen
from sage.groups.generic import order_from_bounds


def _cardinality_with_j_invariant_1728(self):
    r"""
    Special function to compute cardinality when j=1728.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.cardinality import _cardinality_with_j_invariant_1728

    An example with q=p=1 (mod 4)::

        sage: F = GF(10009)
        sage: [_cardinality_with_j_invariant_1728(EllipticCurve(F, [0,0,0,11^i,0]))
        ....:  for i in range(4)]
        [10016, 10210, 10004, 9810]

    An example with q=p=3 (mod 4)::

        sage: F = GF(10007)
        sage: [_cardinality_with_j_invariant_1728(EllipticCurve(F, [0,0,0,5^i,0]))
        ....:  for i in range(4)]
        [10008, 10008, 10008, 10008]

    An example with `q=p^2`, p=3 (mod 4)::

        sage: F.<a> = GF(10007^2,'a')
        sage: [_cardinality_with_j_invariant_1728(EllipticCurve(F,[0,0,0,a^i,0]))
        ....:  for i in range(4)]
        [100160064, 100140050, 100120036, 100140050]

    Examples with `q=2^d`, d odd (3 isomorphism classes)::

        sage: F.<a> = GF(2**15,'a')
        sage: ais = [[0,0,1,0,0],[0,0,1,1,0],[0,0,1,1,1]]
        sage: curves = [EllipticCurve(F,ai) for ai in ais]
        sage: all((e1 == e2 or not e1.is_isomorphic(e2)) for e1 in curves for e2 in curves)
        True
        sage: [_cardinality_with_j_invariant_1728(e) for e in curves]
        [32769, 33025, 32513]

    Examples with `q=2^d`, d even (7 isomorphism classes)::

        sage: F.<a> = GF(2**16,'a')
        sage: b = a^11 # trace 1
        sage: ais = [[0,0,1,0,0], [0,0,1,0,b], [0,0,1,b,0], [0,0,a,0,0],
        ....:        [0,0,a,0,a^2*b], [0,0,a^2,0,0], [0,0,a^2,0,a^4*b]]
        sage: curves = [EllipticCurve(F, ai) for ai in ais]
        sage: all((e1 == e2 or not e1.is_isomorphic(e2)) for e1 in curves for e2 in curves)
        True
        sage: [_cardinality_with_j_invariant_1728(e) for e in curves]
        [65025, 66049, 65537, 65793, 65281, 65793, 65281]

    Examples with `q=3^d`, d odd (4 isomorphism classes)::

        sage: F.<a> = GF(3**15,'a')
        sage: b = a^7  # has trace 1
        sage: ais = [[0,0,0,1,0], [0,0,0,-1,0], [0,0,0,-1,b], [0,0,0,-1,-b]]
        sage: curves = [EllipticCurve(F, ai) for ai in ais]
        sage: all((e1 == e2 or not e1.is_isomorphic(e2)) for e1 in curves for e2 in curves)
        True
        sage: [_cardinality_with_j_invariant_1728(e) for e in curves]
        [14348908, 14348908, 14342347, 14355469]

    Examples with `q=3^d`, d even (6 isomorphism classes)::

        sage: F.<g> = GF(3^18,'g')
        sage: i = F(-1).sqrt()
        sage: a = g^8  # has trace 1
        sage: ais = [[0,0,0,1,0], [0,0,0,1,i*a], [0,0,0,g,0],
        ....:        [0,0,0,g^3,0], [0,0,0,g^2,0], [0,0,0,g^2,i*a*g^3]]
        sage: curves = [EllipticCurve(F, ai) for ai in ais]
        sage: all((e1 == e2 or not e1.is_isomorphic(e2)) for e1 in curves for e2 in curves)
        True
        sage: [_cardinality_with_j_invariant_1728(e) for e in curves]
        [387459856, 387400807, 387420490, 387420490, 387381124, 387440173]

    TESTS:

    Check that a bug noted at :issue:`15667` is fixed::

        sage: F.<a> = GF(3^6,'a')
        sage: EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1]).cardinality()
        784
        sage: EllipticCurve([a^5 + 2*a^3 + 2*a^2 + 2*a, a^4 + a^3 + 2*a + 1]).cardinality_exhaustive()
        784
    """
    try:
        return self._order
    except AttributeError:
        pass

    k = self.base_ring()
    assert self.j_invariant() == k(1728)
    q = k.cardinality()
    p = k.characteristic()
    d = k.degree()
    x = polygen(ZZ)

    # p=2, j=0=1728
    #
    # Number of isomorphism classes is 3 (in odd degree) or 7 (in even degree)
    #
    if p == 2:
        if d % 2:
            # The 3 classes are represented, independently of d,
            # by [0,0,1,0,0], [0,0,1,1,0], [0,0,1,1,1]
            E = EllipticCurve(k, [0, 0, 1, 0, 0])
            if self.is_isomorphic(E):
                t = 0
            else:
                n = (d + 1) // 2
                t = 2**n
                n = n % 4
                if n == 0 or n == 1:
                    t = -t
                E = EllipticCurve(k, [0, 0, 1, 1, 1])
                if self.is_isomorphic(E):
                    t = -t
        else:
            # The 7 classes are represented by E1=[0,0,1,0,0],
            # E2=[0,0,1,0,b], E3=[0,0,1,b,0], E4=[0,0,a,0,0],
            # E4=[0,0,a,0,a^2*b], E6=[0,0,a^2,0,0],
            # E7=[0,0,a^2,0,a^4*b], where a is a non-cube and b
            # has trace 1.  E1's Frobenius is pi=(-2)**(d//2); the
            # Frobeniuses are then pi, -pi, 0; w*pi, -w*pi;
            # w^2*pi, -w^2*pi where w is either cube root of
            # unity, so the traces are 2*pi, -2*pi, 0, -pi, +pi;
            # -pi, +pi.
            delta = self.discriminant()
            discube = (delta**((q - 1) // 3) == k(1))
            pi = (-2)**(d // 2)
            if discube:
                a = k.gen()
                b = a
                while b.trace() == 0:
                    b *= a
                if self.is_isomorphic(EllipticCurve(k, [0, 0, 1, b, 0])):
                    t = 0
                else:
                    t = 2 * pi
                    if not self.is_isomorphic(EllipticCurve(k, [0, 0, 1, 0, 0])):
                        t = -t

            else:
                t = pi
                if self.is_isomorphic(EllipticCurve(k, [0, 0, delta, 0, 0])):
                    t = -t

    # p=3, j=0=1728
    #
    # Number of isomorphism classes is 4 (odd degree) or 6 (even degree)
    #
    elif p == 3:
        if d % 2:
            # The 4 classes are represented by [0,0,0,1,0],
            # [0,0,0,-1,0], [0,0,0,-1,a], [0,0,0,-1,-a] where a
            # has trace 1
            delta = self.discriminant()
            if (-delta).is_square():
                t = 0
            else:
                u = delta.sqrt()
                if not u.is_square():
                    u = -u
                tr = ((self.a3()**2 + self.a6()) / u).trace()
                if tr == 0:
                    t = 0
                else:
                    d2 = (d + 1) // 2
                    t = 3**d2
                    if d2 % 2:
                        t = -t
                    if tr == -1:
                        t = -t
        else:
            # The 6 classes are represented by: [0,0,0,1,0],
            # [0,0,0,1,i*a]; [0,0,0,g,0], [0,0,0,g^3,0];
            # [0,0,0,g^2,0], [0,0,0,g^2,i*a*g^3]; where g
            # generates the multiplicative group modulo 4th
            # powers, and a has nonzero trace.

            # The curve is isomorphic to [0,0,0,A4,A6]

            A4 = self.a4() - self.a1() * self.a3()  # = -b4 = 2*b4
            if A4.is_square():
                u = A4.sqrt()
                t = (-3)**(d // 2)
                i = k(-1).sqrt()
                A6 = self.a3()**2 + self.a6()   # = b6
                if (A6 / (i * u * A4)).trace() == 0:
                    t *= 2
                else:
                    t *= -1
                if not u.is_square():
                    t *= -1
            else:
                t = 0

    # p>3, j=1728
    #
    # Number of isomorphism classes is 4 if q=1 (mod 4), else 2
    #
    elif p % 4 == 3:
        if d % 2:
            t = 0
        else:
            t = (-p)**(d // 2)
            w = (self.c4() / k(48))**((q - 1) // 4)
            if w == 1:
                t = 2 * t
            elif w == -1:
                t = -2 * t
            else:
                t = 0

    # p=1 (mod 4).  First find Frobenius pi=a+b*i for [0,0,0,-1,0] over GF(p):
    # N(pi)=p and N(pi-1)=0 (mod 8).
    #
    else:
        R = ZZ.extension(x**2 + 1, 'i')
        i = R.gen(1)
        pi = R.fraction_field().factor(p)[0][0].gens_reduced()[0]
        a, b = pi.list()
        if a % 2 == 0:
            a, b = -b, a
        if (a + b + 1) % 4 == 0:
            a, b = -a, -b
        pi = a + b * i        # Now pi=a+b*i with (a,b)=(1,0),(3,2) mod 4

        # Lift to Frobenius for [0,0,0,-1,0] over GF(p^d):
        if d > 1:
            pi = pi**d
            a, b = pi.list()

        # Compute appropriate quartic twist:
        w = (self.c4() / k(48))**((q - 1) // 4)
        if w == 1:
            t = 2 * a
        elif w == -1:
            t = -2 * a
        elif k(b) == w * k(a):
            t = 2 * b
        else:
            t = -2 * b

    return Integer(q + 1 - t)


def _cardinality_with_j_invariant_0(self):
    r"""
    Special function to compute cardinality when j=0.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.cardinality import _cardinality_with_j_invariant_0

    An example with q=p=1 (mod 6)::

        sage: F = GF(1009)
        sage: [_cardinality_with_j_invariant_0(EllipticCurve(F,[0,0,0,0,11^i])) for i in range(6)]
        [948, 967, 1029, 1072, 1053, 991]

    An example with q=p=5 (mod 6)::

        sage: F = GF(1013)
        sage: [_cardinality_with_j_invariant_0(EllipticCurve(F,[0,0,0,0,3^i])) for i in range(6)]
        [1014, 1014, 1014, 1014, 1014, 1014]

    An example with `q=p^2`, p=5 (mod 6)::

        sage: F.<a> = GF(1013^2,'a')
        sage: [_cardinality_with_j_invariant_0(EllipticCurve(F,[0,0,0,0,a^i])) for i in range(6)]
        [1028196, 1027183, 1025157, 1024144, 1025157, 1027183]

    For examples in characteristic 2 and 3, see the function
    :func:`_cardinality_with_j_invariant_1728`.
    """
    k = self.base_ring()
    assert self.j_invariant() == k.zero()
    p = k.characteristic()
    if p == 2 or p == 3:  # then 0==1728
        return _cardinality_with_j_invariant_1728(self)

    q = k.cardinality()
    d = k.degree()
    x = polygen(ZZ)

    # p>3, j=0
    #
    # Number of isomorphism classes is 4 if q=1 (mod 4), else 2
    #
    if p % 6 == 5:
        if d % 2:
            t = 0
        else:
            t = (-p)**(d // 2)
            w = (self.c6() / k(-864))**((q - 1) // 6)
            if w == 1:
                t = 2 * t
            elif w == -1:
                t = -2 * t
            elif w**3 == 1:
                t = -t

    # p=1 (mod 6).  First find Frobenius pi=a+b*w for [0,0,0,0,1] over GF(p):
    # N(pi)=p and N(pi-1)=0 (mod 12).
    #
    else:
        R = ZZ.extension(x**2 - x + 1, 'zeta6')
        zeta6 = R.gen(1)
        pi = R.fraction_field().factor(p)[0][0].gens_reduced()[0]
        while (pi - 1).norm() % 12:
            pi *= zeta6
        a, b = pi.list()
        z = k(-b) / k(a)  # a *specific* 6th root of unity in k

        # Now pi=a+b*zeta6 with N(pi-1)=0 (mod 12)

        # Lift to Frobenius for [0,0,0,0,1] over GF(p^d):
        if d > 1:
            pi = pi**d
            a, b = pi.list()

        # Compute appropriate sextic twist:
        w = (self.c6() / k(-864))**((q - 1) // 6)

        if w == 1:
            t = 2 * a + b  # = Trace(pi)
        elif w == -1:
            t = -2 * a - b  # = Trace(-pi)
        elif w == z:
            t = a - b    # = Trace(pi*zeta6)
        elif w == z**2:
            t = -a - 2 * b  # = Trace(pi*zeta6**2)
        elif w == z**4:
            t = b - a    # = Trace(pi*zeta6**4)
        elif w == z**5:
            t = a + 2 * b  # = Trace(pi*zeta6**5)

    return Integer(q + 1 - t)


def cardinality_exhaustive(self):
    r"""
    Return the cardinality of ``self`` over the base field.

    This simply adds up the number of points with each x-coordinate:
    only used for small field sizes!

    EXAMPLES::

        sage: p = next_prime(10^3)
        sage: E = EllipticCurve(GF(p), [3,4])
        sage: E.cardinality_exhaustive()
        1020
        sage: E = EllipticCurve(GF(3^4,'a'), [1,1])
        sage: E.cardinality_exhaustive()
        64
    """
    self._order = Integer(1 + sum(len(self.lift_x(x, all=True))
                                  for x in self.base_field()))
    return self._order


def cardinality_bsgs(self, verbose=False):
    r"""
    Return the cardinality of ``self`` over the base field.

    ALGORITHM: A variant of "Mestre's trick" extended to all finite
    fields by Cremona and Sutherland, 2008.

    .. NOTE::

        1. The Mestre-Schoof-Cremona-Sutherland algorithm may fail for
           a small finite number of curves over `F_q` for `q` at most 49, so
           for `q<50` we use an exhaustive count.

        2. Quadratic twists are not implemented in characteristic 2
           when `j=0 (=1728)`; but this case is treated separately.

    EXAMPLES::

        sage: p = next_prime(10^3)
        sage: E = EllipticCurve(GF(p), [3,4])
        sage: E.cardinality_bsgs()
        1020
        sage: E = EllipticCurve(GF(3^4,'a'), [1,1])
        sage: E.cardinality_bsgs()
        64
        sage: F.<a> = GF(101^3,'a')
        sage: E = EllipticCurve([2*a^2 + 48*a + 27, 89*a^2 + 76*a + 24])
        sage: E.cardinality_bsgs()
        1031352
    """
    E1 = self
    k = self.base_field()
    q = k.order()
    if q < 50:
        if verbose:
            print("q=", q, "< 50 so using exhaustive count")
        return cardinality_exhaustive(self)

    # Construct the quadratic twist:
    E2 = E1.quadratic_twist()
    if verbose:
        print("Quadratic twist is ", E2.ainvs())

    bounds = Hasse_bounds(q)
    lower, upper = bounds
    B = upper - q - 1  # = floor(2*sqrt(q))
    a = ZZ(0)
    N1 = N2 = M = ZZ(1)
    kmin = -B
    kmax = B
    q1 = q + 1
    # Throughout, we have #E=q+1-t where |t|<=B and t=a+k*M = a
    # (mod M) where kmin <= k <= kmax.

    # M is the lcm of the orders of all the points found on E1 and
    # E2, which will eventually exceed 2*B, at which point
    # kmin=kmax.

    if q > 2**10:
        N1 = ZZ(2)**sum([e for P, e in E1._p_primary_torsion_basis(2)])
        N2 = ZZ(2)**sum([e for P, e in E2._p_primary_torsion_basis(2)])
        if q > 2**20:
            N1 *= ZZ(3)**sum([e for P, e in E1._p_primary_torsion_basis(3)])
            N2 *= ZZ(3)**sum([e for P, e in E2._p_primary_torsion_basis(3)])
            if q > 2**40:
                N1 *= ZZ(5)**sum([e for P, e in E1._p_primary_torsion_basis(5)])
                N2 *= ZZ(5)**sum([e for P, e in E2._p_primary_torsion_basis(5)])
        # We now know that t=q+1 (mod N1) and t=-(q+1) (mod N2)
        a = q1
        M = N1
        g, u, v = M.xgcd(N2)  # g==u*M+v*N2
        if N2 > g:
            a = (a * v * N2 - q1 * u * M) // g
            M *= (N2 // g)  # = lcm(M,N2)
            a = a % M
            if verbose:
                print("(a,M)=", (a, M))
            kmin = ((-B - a) / M).ceil()
            kmax = ((B - a) / M).floor()
            if kmin == kmax:
                self._order = q1 - a - kmin * M
                if verbose:
                    print("no random points were needed")
                return self._order
        if verbose:
            print("(2,3,5)-torsion subgroup gives M=", M)

    # N1, N2 are divisors of the orders of E1, E2 separately,
    # which are used to speed up the computation of the orders of
    # points on each curve.  For large q it is worth initializing
    # these with the full order of the (2,3,5)-torsion which are
    # often non-trivial.
    while kmax != kmin:
        # Get a random point on E1 and find its order, using the
        # Hasse bounds and the fact that we know that the group
        # order is a multiple of N1:
        n = order_from_bounds(E1.random_point(), bounds, N1, operation='+')
        if verbose:
            print("New point on E has order ", n)
        # update N1 and M
        N1 = N1.lcm(n)
        g, u, v = M.xgcd(n)  # g==u*M+v*n
        if n > g:
            # update congruence a (mod M) with q+1 (mod n)
            a = (a * v * n + q1 * u * M) // g
            M *= (n // g)  # = lcm(M,n)
            a = a % M
            if verbose:
                print("(a,M)=", (a, M))
            kmin = ((-B - a) / M).ceil()
            kmax = ((B - a) / M).floor()
            if kmin == kmax:
                self._order = q1 - a - kmin * M
                return self._order
            if verbose:
                print("number of possibilities is now ", kmax - kmin + 1)

        # Get a random point on E2 and find its order, using the
        # Hasse bounds and the fact that we know that the group
        # order is a multiple of N2:
        n = order_from_bounds(E2.random_point(), bounds, N2, operation='+')
        if verbose:
            print("New point on E' has order ", n)
        # update N2 and M
        N2 = N2.lcm(n)
        g, u, v = M.xgcd(n)  # g==u*M+v*n
        if n > g:
            # update congruence a (mod M) with -(q+1) (mod n)
            a = (a * v * n - q1 * u * M) // g
            M *= (n // g)  # = lcm(M,n)
            a = a % M
            if verbose:
                print("(a,M)=", (a, M))
            kmin = ((-B - a) / M).ceil()
            kmax = ((B - a) / M).floor()
            if kmin == kmax:
                self._order = q1 - a - kmin * M
                return self._order
            if verbose:
                print("number of possibilities is now ", kmax - kmin + 1)


def _cardinality_subfield(self, jpol):
    """
    Count the number of points on the elliptic curve ``self``, assuming
    that the j-invariant lies in a subfield.

    INPUT:

    - ``jpol`` -- minimal polynomial (over the prime field) of the
      j-invariant of ``self``

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.cardinality import _cardinality_subfield
        sage: k.<a> = GF(7^5)
        sage: E = EllipticCurve(k, [1,2,3,4,5]); E
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
         over Finite Field in a of size 7^5
        sage: _cardinality_subfield(E, E.j_invariant().minimal_polynomial())
        17019

    TESTS::

        sage: k.<a> = GF(2^6)
        sage: for c in range(500):  # long time
        ....:     ainvs = [k.random_element() for i in range(5)]
        ....:     try:
        ....:         E = EllipticCurve(ainvs)
        ....:     except ArithmeticError:
        ....:         continue  # Singular curve
        ....:     jpol = E.j_invariant().minimal_polynomial()
        ....:     if jpol.degree() < 6:
        ....:         N1 = _cardinality_subfield(E, jpol)
        ....:         N2 = pari.ellcard(E)
        ....:         assert N1 == N2
    """
    k = self.base_ring()
    p = k.characteristic()
    d = k.degree()

    jdeg = jpol.degree()
    if jdeg >= d:
        raise ValueError("j-invariant does not lie in a subfield")

    # we count points on a standard curve with the same
    # j-invariant, defined over the field it generates, then lift
    # to the curve's own field, and finally allow for twists

    # Let j be the j-invariant as element of the smallest finite
    # field over which j is defined.
    GFj = GF((p, jdeg), name='j', modulus=jpol)
    j = GFj.gen()

    # Use special code for j = 0, 1728
    if j == 1728:
        return _cardinality_with_j_invariant_1728(self)
    elif j == 0:
        return _cardinality_with_j_invariant_0(self)

    # Recursive call which does all the real work:
    E0 = EllipticCurve_from_j(j)
    N = E0.cardinality(extension_degree=d // jdeg)

    # Map to the original larger field
    phi = GFj.hom([self.j_invariant()])

    # Check if we got a twist. Since j is not 0, 1728 the only twists
    # are quadratic.
    if self.is_isomorphic(E0.base_extend(phi)):
        return N
    else:
        q = k.cardinality()
        return 2 * (q + 1) - N
