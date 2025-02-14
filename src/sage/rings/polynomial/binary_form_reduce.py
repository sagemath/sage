r"""
Helper functions for reduction of binary forms.

The algorithm for reducing is from Stoll and Cremona's "On the Reduction Theory of
Binary Forms" [CS2003]_. This takes a two variable homogeneous polynomial and finds a
reduced form. This is an `SL(2,\ZZ)`-equivalent binary form whose covariant in
the upper half plane is in the fundamental domain. Further, the algorithm
from Hutz and Stoll [HS2018]_ allows the form to be further minimized so that
the coefficients have either smallest height or smallest `L_2` norm.

AUTHORS:

- Rebecca Lauren Miller -- initial version of reduction as part of GSOC 2016

- Ben Hutz (2018-7) -- improvements to reduce and implement smallest coefficient model
"""

# ****************************************************************************
#       Copyright (C) 2018 Benjamin Hutz <bn4941#gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import matrix
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.cc import CC
from sage.rings.complex_mpfr import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.integer_ring import ZZ
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField

lazy_import("sage.calculus.functions", "jacobian")
lazy_import("sage.functions.hyperbolic", ["cosh", "sinh"])
lazy_import("sage.functions.log", "exp")


def covariant_z0(F, z0_cov=False, prec=53, emb=None, error_limit=0.000001):
    r"""
    Return the covariant and Julia invariant from Cremona-Stoll [CS2003]_.

    In [CS2003]_ and [HS2018]_ the Julia invariant is denoted as `\Theta(F)`
    or `R(F, z(F))`. Note that you may get faster convergence if you first move
    `z_0(F)` to the fundamental domain before computing the true covariant

    INPUT:

    - ``F`` -- binary form of degree at least 3 with no multiple roots

    - ``z0_cov`` -- boolean, compute only the `z_0` invariant. Otherwise, solve
      the minimization problem

    - ``prec`` -- positive integer. precision to use in CC

    - ``emb`` -- embedding into CC

    - ``error_limit`` -- sets the error tolerance (default: 0.000001)

    OUTPUT: a complex number, a real number

    EXAMPLES::

        sage: from sage.rings.polynomial.binary_form_reduce import covariant_z0
        sage: R.<x,y> = QQ[]
        sage: F = 19*x^8 - 262*x^7*y + 1507*x^6*y^2 - 4784*x^5*y^3 + 9202*x^4*y^4\
        ....: - 10962*x^3*y^5 + 7844*x^2*y^6 - 3040*x*y^7 + 475*y^8
        sage: covariant_z0(F, prec=80, z0_cov=True)
        (1.3832330115323681438175 + 0.31233552177413614978744*I,
         3358.4074848663492819259)
        sage: F = -x^8 + 6*x^7*y - 7*x^6*y^2 - 12*x^5*y^3 + 27*x^4*y^4\
        ....: - 4*x^3*y^5 - 19*x^2*y^6 + 10*x*y^7 - 5*y^8
        sage: covariant_z0(F, prec=80)
        (0.64189877107807122203366 + 1.1852516565091601348355*I,
        3134.5148284344627168276)

    ::

        sage: R.<x,y> = QQ[]
        sage: covariant_z0(x^3 + 2*x^2*y - 3*x*y^2, z0_cov=True)[0]
        0.230769230769231 + 0.799408065031789*I
        sage: -1/covariant_z0(-y^3 + 2*y^2*x + 3*y*x^2, z0_cov=True)[0]
        0.230769230769231 + 0.799408065031789*I

    ::

        sage: R.<x,y> = QQ[]
        sage: covariant_z0(2*x^2*y - 3*x*y^2, z0_cov=True)[0]
        0.750000000000000 + 1.29903810567666*I
        sage: -1/covariant_z0(-x^3 - x^2*y + 2*x*y^2, z0_cov=True)[0] + 1
        0.750000000000000 + 1.29903810567666*I

    ::

        sage: R.<x,y> = QQ[]
        sage: covariant_z0(x^2*y - x*y^2, prec=100) # tol 1e-28
         (0.50000000000000000000000000003 + 0.86602540378443864676372317076*I,
         1.5396007178390020386910634147)

    TESTS::

        sage: R.<x,y>=QQ[]
        sage: covariant_z0(x^2 + 24*x*y + y^2)
        Traceback (most recent call last):
        ...
        ValueError: must be at least degree 3
        sage: covariant_z0((x+y)^3, z0_cov=True)
        Traceback (most recent call last):
        ...
        ValueError: cannot have multiple roots for z0 invariant
        sage: covariant_z0(x^3 + 3*x*y + y)
        Traceback (most recent call last):
        ...
        TypeError: must be a binary form
        sage: covariant_z0(-2*x^2*y^3 + 3*x*y^4 + 127*y^5)
        Traceback (most recent call last):
        ...
        ValueError: cannot have a root with multiplicity >= 5/2
        sage: covariant_z0((x^2+2*y^2)^2)
        Traceback (most recent call last):
        ...
        ValueError: must have at least 3 distinct roots
    """
    R = F.parent()
    d = ZZ(F.degree())
    if R.ngens() != 2 or any(sum(t) != d for t in F.exponents()):
        raise TypeError('must be a binary form')
    if d < 3:
        raise ValueError('must be at least degree 3')

    f = F.subs({R.gen(1): 1}).univariate_polynomial()
    if f.degree() < d:
        # we have a root at infinity
        if f.constant_coefficient() != 0:
            # invert so we find all roots!
            mat = matrix(ZZ, 2, 2, [0, -1, 1, 0])
        else:
            t = 0
            while f(t) == 0:
                t += 1
            mat = matrix(ZZ, 2, 2, [t, -1, 1, 0])
    else:
        mat = matrix(ZZ, 2, 2, [1, 0, 0, 1])
    f = F(list(mat * vector(R.gens()))).subs({R.gen(1): 1}).univariate_polynomial()
    # now we have a single variable polynomial with all the roots of F
    K = ComplexField(prec=prec)
    if f.base_ring() != K:
        if emb is None:
            f = f.change_ring(K)
        else:
            f = f.change_ring(emb)
    roots = f.roots()
    if max(ex for _, ex in roots) > 1 or f.degree() < d - 1:
        if z0_cov:
            raise ValueError('cannot have multiple roots for z0 invariant')
        else:
            # just need a starting point for Newton's method
            f = f.lc() * prod(p for p, ex in f.factor())  # removes multiple roots
            if f.degree() < 3:
                raise ValueError('must have at least 3 distinct roots')
            roots = f.roots()
    roots = [p for p, _ in roots]

    # finding quadratic Q_0, gives us our covariant, z_0
    dF = f.derivative()
    n = ZZ(f.degree())
    PR = PolynomialRing(K, 'x,y')
    x, y = PR.gens()
    # finds Stoll and Cremona's Q_0
    q = sum([(1/(dF(r).abs()**(2/(n-2)))) * ((x-(r*y)) * (x-(r.conjugate()*y)))
             for r in roots])
    # this is Q_0 , always positive def as long as F has distinct roots
    A = q.monomial_coefficient(x**2)
    B = q.monomial_coefficient(x * y)
    C = q.monomial_coefficient(y**2)
    # need positive root
    try:
        z = ((-B + ((B**2)-(4*A*C)).sqrt()) / (2 * A))
    except ValueError:
        raise ValueError("not enough precision")
    if z.imag() < 0:
        z = (-B - ((B**2)-(4*A*C)).sqrt()) / (2 * A)

    if z0_cov:
        FM = f  # for Julia's invariant
    else:
        # solve the minimization problem for 'true' covariant
        CF = ComplexIntervalField(prec=prec)  # keeps trac of our precision error
        z = CF(z)
        FM = F(list(mat * vector(R.gens()))).subs({R.gen(1): 1}).univariate_polynomial()
        from sage.rings.polynomial.complex_roots import complex_roots
        L1 = complex_roots(FM, min_prec=prec)
        L = []
        # making sure multiplicity isn't too large using convergence conditions in paper
        for p, e in L1:
            if e >= d / 2:
                raise ValueError('cannot have a root with multiplicity >= %s/2' % d)
            L.extend(p for _ in range(e))
        RCF = PolynomialRing(CF, 'u,t')
        a = RCF.zero()
        c = RCF.zero()
        u, t = RCF.gens()
        for l in L:
            denom = ((t - l) * (t - l.conjugate()) + u**2)
            a += u**2 / denom
            c += (t - l.real()) / denom
        # Newton's Method, to find solutions. Error bound is less than diameter of our z
        err = z.diameter()
        zz = z.diameter()
        g1 = a.numerator() - d / 2 * a.denominator()
        g2 = c.numerator()
        G = vector([g1, g2])
        J = jacobian(G, [u, t])
        v0 = vector([z.imag(), z.real()])  # z0 as starting point
        # finds our correct z
        while err <= zz:
            NJ = J.subs({u: v0[0], t: v0[1]})
            NJinv = NJ.inverse()
            # inverse for CIF matrix seems to return fractions not CIF elements, fix them
            if NJinv.base_ring() != CF:
                NJinv = matrix(CF, 2, 2, [CF(zw.numerator() / zw.denominator())
                                          for zw in NJinv.list()])
            w = z
            v0 = v0 - NJinv*G.subs({u: v0[0], t: v0[1]})
            z = v0[1].constant_coefficient() + v0[0].constant_coefficient()*CF.gen(0)
            err = z.diameter() # precision
            zz = (w - z).abs().lower() # difference in w and z
        # despite there is no break, this happens
        if err > error_limit or err.is_NaN():
            raise ValueError("accuracy of Newton's root not within tolerance(%s > %s), increase precision" % (err, error_limit))
        if z.imag().upper() <= z.diameter():
            raise ArithmeticError("Newton's method converged to z not in the upper half plane")
        z = z.center()

    # Julia's invariant
    if FM.base_ring() != ComplexField(prec=prec):
        FM = FM.change_ring(ComplexField(prec=prec))
    tF = z.real()
    uF = z.imag()
    th = FM.lc().abs()**2
    for r, ex in FM.roots():
        for _ in range(ex):
            th = th * ((((r-tF).abs())**2 + uF**2)/uF)

    # undo shift and invert (if needed)
    # since F \cdot m ~ m^(-1)\cdot z
    # we apply m to z to undo m acting on F
    l = mat * vector([z, 1])
    return l[0] / l[1], th


# // compute inverse of eps_F
# from Michael Stoll
def epsinv(F, target, prec=53, target_tol=0.001, z=None, emb=None):
    """
    Compute a bound on the hyperbolic distance.

    The true minimum will be within the computed bound.
    It is computed as the inverse of epsilon_F from [HS2018]_.

    INPUT:

    - ``F`` -- binary form of degree at least 3 with no multiple roots

    - ``target`` -- positive real number; the value we want to attain, i.e.,
      the value we are taking the inverse of

    - ``prec`` -- positive integer; precision to use in CC

    - ``target_tol`` -- positive real number; the tolerance with which we
      attain the target value

    - ``z`` -- complex number; ``z_0`` covariant for F

    - ``emb`` -- embedding into CC

    OUTPUT: a real number delta satisfying  target + target_tol > eps_F(delta) > target.

    EXAMPLES::

        sage: from sage.rings.polynomial.binary_form_reduce import epsinv
        sage: R.<x,y> = QQ[]
        sage: epsinv(-2*x^3 + 2*x^2*y + 3*x*y^2 + 127*y^3, 31.5022020249597) # tol 1e-12
        4.02520895942207
    """
    def RQ(delta):
        # this is the quotient R(F_0,z)/R(F_0,z(F)) for a generic z
        # at distance delta from j. See Lemma 4.2 in [HS2018].
        cd = cosh(delta).n(prec=prec)
        sd = sinh(delta).n(prec=prec)
        return prod([cd + (cost * phi[0] + sint * phi[1]) * sd for phi in phis])

    def epsF(delta):
        pol = RQ(delta)  # get R quotient in terms of z
        S = PolynomialRing(C, 'v')
        g = S([(i - d) * pol[i - d] for i in range(2 * d + 1)])  # take derivative
        drts = [e for e in g.roots(ring=C, multiplicities=False)
                if (e.norm() - 1).abs() < 0.1]
        # find min
        return min([pol(r / r.abs()).real() for r in drts])

    C = ComplexField(prec=prec)
    R = F.parent()
    d = F.degree()
    if z is None:
        z, th = covariant_z0(F, prec=prec, emb=emb)
    else:  # need to do our own input checking
        if R.ngens() != 2 or any(sum(t) != d for t in F.exponents()):
            raise TypeError('must be a binary form')
        if d < 3:
            raise ValueError('must be at least degree 3')

    f = F.subs({R.gen(1): 1}).univariate_polynomial()
    # now we have a single variable polynomial
    if (max(ex for p, ex in f.roots(ring=C)) >= QQ(d)/2 or
            f.degree() < QQ(d)/2):
        raise ValueError('cannot have root with multiplicity >= deg(F)/2')

    R = RealField(prec=prec)
    PR = PolynomialRing(R, 't')
    t = PR.gen(0)
    # compute phi_1, ..., phi_k
    # first find F_0 and its roots
    # this change of variables on f moves z(f) to j, i.e. produces F_0
    rts = f(z.imag()*t + z.real()).roots(ring=C)
    phis = []  # stereographic projection of roots
    for r, e in rts:
        phis.extend([[2*r.real()/(r.norm()+1), (r.norm()-1)/(r.norm()+1)]])
    if d != f.degree():  # include roots at infinity
        phis.extend([(d - f.degree()) * [0, 1]])

    # for writing RQ in terms of generic z to minimize
    LC = LaurentSeriesRing(C, 'u', default_prec=2 * d + 2)
    u = LC.gen(0)
    cost = (u + u**(-1)) / 2
    sint = (u - u**(-1)) / (2 * C.gen(0))

    # first find an interval containing the desired value
    # then use regula falsi on log eps_F
    # d -> delta value in interval [0,1]
    # v in value in interval [1,epsF(1)]
    dl = R(0.0)
    vl = R(1.0)
    du = R(1.0)
    vu = epsF(du)
    while vu < target:
        # compute the next value of epsF for delta = 2*delta
        dl = du
        vl = vu
        du *= 2
        vu = epsF(du)
    # now dl < delta <= du
    logt = target.log()
    l2 = (vu.log() - logt).n(prec=prec)
    l1 = (vl.log() - logt).n(prec=prec)
    dn = (dl*l2 - du*l1)/(l2 - l1)
    vn = epsF(dn)
    dl = du
    vl = vu
    du = dn
    vu = vn
    while (du - dl).abs() >= target_tol or max(vl, vu) < target:
        l2 = (vu.log() - logt).n(prec=prec)
        l1 = (vl.log() - logt).n(prec=prec)
        dn = (dl * l2 - du * l1) / (l2 - l1)
        vn = epsF(dn)
        dl = du
        vl = vu
        du = dn
        vu = vn
    return max(dl, du)


def get_bound_poly(F, prec=53, norm_type='norm', emb=None):
    """
    The hyperbolic distance from `j` which must contain the smallest poly.

    This defines the maximum possible distance from `j` to the `z_0` covariant
    in the hyperbolic 3-space for which the associated `F` could have smaller
    coefficients.

    INPUT:

    - ``F`` -- binary form of degree at least 3 with no multiple roots

    - ``prec`` -- positive integer. precision to use in CC

    - ``norm_type`` -- string, either norm or height

    - ``emb`` -- embedding into CC

    OUTPUT: a positive real number

    EXAMPLES::

        sage: from sage.rings.polynomial.binary_form_reduce import get_bound_poly
        sage: R.<x,y> = QQ[]
        sage: F = -2*x^3 + 2*x^2*y + 3*x*y^2 + 127*y^3
        sage: get_bound_poly(F) # tol 1e-12
        28.0049336543295
        sage: get_bound_poly(F, norm_type='height') # tol 1e-11
        111.890642019092
    """
    if F.base_ring() != ComplexField(prec=prec):
        if emb is None:
            compF = F.change_ring(ComplexField(prec=prec))
        else:
            compF = F.change_ring(emb)
    else:
        compF = F
    n = F.degree()
    assert (n > 2), "degree 2 polynomial"

    z0F, thetaF = covariant_z0(compF, prec=prec, emb=emb)
    if norm_type == 'norm':
        # euclidean norm squared
        normF = (sum([abs(i)**2 for i in compF.coefficients()]))
        target = (2**(n - 1)) * normF / thetaF
    elif norm_type == 'height':
        hF = exp(max([c.global_height(prec=prec) for c in F.coefficients()]))  # height
        target = (2**(n - 1)) * (n + 1) * (hF**2) / thetaF
    else:
        raise ValueError('type must be norm or height')
    return cosh(epsinv(F, target, prec=prec))


def smallest_poly(F, prec=53, norm_type='norm', emb=None):
    r"""
    Determine the poly with smallest coefficients in `SL(2,\Z)` orbit of ``F``.

    Smallest can be in the sense of `L_2` norm or height.
    The method is the algorithm in Hutz-Stoll [HS2018]_.

    ``F`` needs to be a binary form with no multiple roots of degree
    at least 3. It should already be reduced in the sense of
    Cremona-Stoll [CS2003]_.

    INPUT:

    - ``F`` -- binary form of degree at least 3 with no multiple roots

    - ``norm_type`` -- string; ``'norm'`` or ``'height'`` controlling what
      ``smallest`` means for the coefficients

    OUTPUT: pair [poly, matrix]

    EXAMPLES::

        sage: from sage.rings.polynomial.binary_form_reduce import smallest_poly
        sage: R.<x,y> = QQ[]
        sage: F = -x^8 + 6*x^7*y - 7*x^6*y^2 - 12*x^5*y^3 + 27*x^4*y^4\
        ....: - 4*x^3*y^5 - 19*x^2*y^6 + 10*x*y^7 - 5*y^8
        sage: smallest_poly(F, prec=100) #long time
        [
        -x^8 - 2*x^7*y + 7*x^6*y^2 + 16*x^5*y^3 + 2*x^4*y^4 - 2*x^3*y^5 + 4*x^2*y^6 - 5*y^8,
        <BLANKLINE>
        [1 1]
        [0 1]
        ]

    ::

        sage: from sage.rings.polynomial.binary_form_reduce import smallest_poly, get_bound_poly
        sage: R.<x,y> = QQ[]
        sage: F = -2*x^3 + 2*x^2*y + 3*x*y^2 + 127*y^3
        sage: smallest_poly(F)
        [
                                               [1 4]
        -2*x^3 - 22*x^2*y - 77*x*y^2 + 43*y^3, [0 1]
        ]
        sage: F0, M = smallest_poly(F, norm_type='height')
        sage: F0, M  # random
        (
                                                [5 4]
        -58*x^3 - 47*x^2*y + 52*x*y^2 + 43*y^3, [1 1]
        )
        sage: M in SL2Z, F0 == R.hom(M * vector([x, y]))(F)
        (True, True)
        sage: get_bound_poly(F0, norm_type='height')  # tol 1e-12
        23.3402702199809

    An example with a multiple root::

        sage: R.<x,y> = QQ[]
        sage: F = -16*x^7 - 114*x^6*y - 345*x^5*y^2 - 599*x^4*y^3 - 666*x^3*y^4\
        ....: - 481*x^2*y^5 - 207*x*y^6 - 40*y^7
        sage: F.reduced_form()
        (
                                                              [-1 -1]
        -x^5*y^2 - 24*x^3*y^4 - 3*x^2*y^5 - 2*x*y^6 + 16*y^7, [ 1  0]
        )
    """
    def insert_item(pts, item, index):
        # binary insertion to maintain list of points left to consider
        N = len(pts)
        if N == 0:
            return [item]
        elif N == 1:
            if item[index] > pts[0][index]:
                pts.insert(0, item)
            else:
                pts.append(item)
            return pts
        else:  # binary insertion
            left = 1
            right = N
            mid = (left + right) // 2  # these are ints so this is .floor()
            if item[index] > pts[mid][index]:  # item goes into first half
                return insert_item(pts[:mid], item, index) + pts[mid:N]
            else:  # item goes into second half
                return pts[:mid] + insert_item(pts[mid:N], item, index)

    def coshdelta(z):
        # The cosh of the hyperbolic distance from z = t+uj to j
        return (z.norm() + 1)/(2*z.imag())  # reduce in the sense of Cremona-Stoll
    G = F
    MG = matrix(ZZ, 2, 2, [1, 0, 0, 1])
    x, y = G.parent().gens()
    if norm_type == 'norm':
        current_size = sum([abs(i)**2 for i in G.coefficients()])  # euclidean norm squared
    elif norm_type == 'height':  # height
        current_size = exp(max([c.global_height(prec=prec) for c in G.coefficients()]))
    else:
        raise ValueError('type must be norm or height')
    v0, th = covariant_z0(G, prec=prec, emb=emb)
    rep = 2 * CC.gen(0)  # representative point in fundamental domain
    from math import isnan
    if isnan(v0.abs()):
        raise ValueError("invalid covariant: %s" % v0)
    R = get_bound_poly(G, prec=prec, norm_type=norm_type)

    # check orbit
    S = matrix(ZZ, 2, 2, [0, -1, 1, 0])
    T = matrix(ZZ, 2, 2, [1, 1, 0, 1])
    TI = matrix(ZZ, 2, 2, [1, -1, 0, 1])

    count = 0
    pts = [[G, v0, rep, MG, coshdelta(v0), 0]]  # label - 0:None, 1:S, 2:T, 3:T^(-1)
    current_min = [G, v0, rep, MG, coshdelta(v0)]
    while pts:
        G, v, rep, M, D, label = pts.pop()
        # apply ST and keep z, Sz
        if D > R:
            break  # all remaining pts are too far away
        # check if it is smaller. If so, we can improve the bound
        count += 1
        if norm_type == 'norm':
            new_size = sum([abs(i)**2 for i in G.coefficients()])  # euclidean norm squared
        else:  # height
            new_size = exp(max([c.global_height(prec=prec) for c in G.coefficients()]))
        if new_size < current_size:
            current_min = [G, v, rep, M, coshdelta(v)]
            current_size = new_size
            R = get_bound_poly(G, norm_type=norm_type, prec=prec, emb=emb)

        # add new points to check
        if label != 1 and min((rep+1).norm(), (rep-1).norm()) >= 1:  # don't undo S
            # the 2nd condition is equivalent to |\Re(-1/rep)| <= 1/2
            # this means that rep can have resulted from an inversion step in
            # the shift-and-invert procedure, so don't invert

            # do inversion
            z = -1 / v
            new_pt = [G.subs({x: -y, y: x}), z, -1/rep, M*S, coshdelta(z), 1]
            pts = insert_item(pts, new_pt, 4)
        if label != 3:  # don't undo TI
            # do right shift
            z = v - 1
            new_pt = [G.subs({x: x + y}), z, rep-1, M*T, coshdelta(z), 2]
            pts = insert_item(pts, new_pt, 4)
        if label != 2:  # don't undo T
            # do left shift
            z = v + 1
            new_pt = [G.subs({x: x - y}), z, rep + 1, M * TI, coshdelta(z), 3]
            pts = insert_item(pts, new_pt, 4)

    return [current_min[0], current_min[3]]
