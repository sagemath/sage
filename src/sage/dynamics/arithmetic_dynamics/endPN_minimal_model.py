r"""
Sage functions to compute minimal models of rational functions
under the conjugation action of `PGL_2(QQ)`.

AUTHORS:

- Alex Molnar (May 22, 2012)

- Brian Stout, Ben Hutz (Nov 2013): Modified code to use projective
  morphism functionality so that it can be included in Sage.

REFERENCES: [BM2012]_, [Mol2015]_
"""

# ****************************************************************************
#       Copyright (C) 2012 Alexander Molnar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.functions.hyperbolic import cosh
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.cc import CC
from sage.rings.complex_mpfr import ComplexField
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.binary_form_reduce import covariant_z0, epsinv
from sage.rings.rational_field import QQ
from sage.schemes.affine.affine_space import AffineSpace
from sage.symbolic.constants import e
from sage.arith.all import gcd
from copy import copy

def bCheck(c, v, p, b):
    r"""
    Compute a lower bound on the value of ``b``.

    This value is needed for a transformation
    `A(z) = z*p^k + b` to satisfy `ord_p(Res(\phi^A)) < ord_p(Res(\phi))` for a
    rational map `\phi`. See Theorem 3.3.5 in [Molnar]_.

    INPUT:

    - ``c`` -- a list of polynomials in `b`. See v for their use

    - ``v`` -- a list of rational numbers, where we are considering the inequalities
      `ord_p(c[i]) > v[i]`

    - ``p`` -- a prime

    - ``b`` -- local variable

    OUTPUT: ``bval`` -- Integer, lower bound in Theorem 3.3.5

    EXAMPLES::

        sage: R.<b> = PolynomialRing(QQ)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import bCheck
        sage: bCheck(11664*b^2 + 70227*b + 76059, 15/2, 3, b)
        -1
    """
    val = (v+1).floor()
    deg = c.degree()
    coeffs = c.coefficients(sparse=False)
    lcoeff = coeffs[deg]
    coeffs.remove(lcoeff)
    check1 = [(coeffs[i].valuation(p) - lcoeff.valuation(p))/(deg - i) for i in range(0,len(coeffs)) if coeffs[i] != 0]
    check2 = (val - lcoeff.valuation(p))/deg
    check1.append(check2)
    bval = min(check1)
    return (bval).ceil()


def scale(c, v, p):
    r"""
    Create scaled integer polynomial with respect to prime ``p``.

    Given an integral polynomial ``c``, we can write `c = p^i*c'`, where ``p`` does not
    divide ``c``. Return ``c'`` and `v - i` where `i` is the smallest valuation of the
    coefficients of `c`.

    INPUT:

    - ``c`` -- an integer polynomial

    - ``v`` -- an integer - the bound on the exponent from blift

    - ``p`` -- a prime

    OUTPUT:

    - boolean -- the new exponent bound is 0 or negative

    - the scaled integer polynomial

    - an integer the new exponent bound

    EXAMPLES::

        sage: R.<b> = PolynomialRing(QQ)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import scale
        sage: scale(24*b^3 + 108*b^2 + 162*b + 81, 1, 3)
        [False, 8*b^3 + 36*b^2 + 54*b + 27, 0]
    """
    scaleval = min([coeff.valuation(p) for coeff in c.coefficients()])
    if scaleval > 0:
        c = c/(p**scaleval)
        v = v - scaleval
    if v <= 0:
        flag = False
    else:
        flag = True
    return [flag, c, v]


def blift(LF, Li, p, k, S=None, all_orbits=False):
    r"""
    Search for a solution to the given list of inequalities.

    If found, lift the solution to
    an appropriate valuation. See Lemma 3.3.6 in [Molnar]_

    INPUT:

    - ``LF`` -- a list of integer polynomials in one variable (the normalized coefficients)

    - ``Li`` -- an integer, the bound on coefficients

    - ``p`` -- a prime

    - ``k`` -- the scaling factor that makes the solution a ``p``-adic integer

    - ``S`` -- polynomial ring to use

    - ``all_orbits`` -- boolean; whether or not to use ``==`` in the
      inequalities to find all orbits

    OUTPUT:

    - boolean -- whether or not the lift is successful

    - integer -- the lift

    EXAMPLES::

        sage: R.<b> = PolynomialRing(QQ)
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import blift
        sage: blift([8*b^3 + 12*b^2 + 6*b + 1, 48*b^2 + 483*b + 117, 72*b + 1341,\
        ....: -24*b^2 + 411*b + 99, -144*b + 1233, -216*b], 2, 3, 2)
        [[True, 4]]
    """

    P = LF[0].parent()
    #Determine which inequalities are trivial, and scale the rest, so that we only lift
    #as many times as needed.
    keepScaledIneqs = [scale(P(coeff),Li,p) for coeff in LF if coeff != 0]
    keptVals = [i[2] for i in keepScaledIneqs if i[0]]
    if keptVals:
        # Determine the valuation to lift until.
        # liftval = max(keptVals)
        pass
    else:
        #All inequalities are satisfied.
        if all_orbits:
            return [[True, t] for t in range(p)]
        return [[True, 1]]
    if S is None:
        S = PolynomialRing(Zmod(p), 'b')
    keptScaledIneqs = [S(i[1]) for i in keepScaledIneqs if i[0]]
    #We need a solution for each polynomial on the left hand side of the inequalities,
    #so we need only find a solution for their gcd.
    g = gcd(keptScaledIneqs)
    rts = g.roots(multiplicities=False)
    good = []
    for r in rts:
        #Recursively try to lift each root
        r_initial = QQ(r)
        newInput = P([r_initial, p])
        LG = [F(newInput) for F in LF]
        new_good = blift(LG, Li, p, k, S=S)
        for lift,lifted in new_good:
            if lift:
                #Lift successful.
                if not all_orbits:
                    return [[True, r_initial + p*lifted]]

                #only need up to SL(2,ZZ) equivalence
                #this helps control the size of the resulting coefficients
                if r_initial + p*lifted < p**k:
                    good.append([True, r_initial + p*lifted])
                else:
                    new_r = r_initial + p*lifted - p**k
                    while new_r > p**k:
                        new_r -= p**k
                    if [True, new_r] not in good:
                        good.append([True, new_r])
    if good:
        return good
    #Lift non successful.
    return [[False,0]]


def affine_minimal(vp, return_transformation=False, D=None, quick=False):
    r"""
    Determine if given map is affine minimal.

    Given vp a scheme morphism on the projective line over the rationals,
    this procedure determines if `\phi` is minimal. In particular, it determines
    if the map is affine minimal, which is enough to decide if it is minimal
    or not. See Proposition 2.10 in [BM2012]_.

    INPUT:

    - ``vp`` -- dynamical system on the projective line

    - ``D`` -- a list of primes, in case one only wants to check minimality
      at those specific primes

    - ``return_transformation`` -- (default: ``False``) boolean; this
      signals a return of the `PGL_2` transformation to conjugate
      this map to the calculated models

    - ``quick`` -- a boolean value. If true the algorithm terminates once
      algorithm determines F/G is not minimal, otherwise algorithm only
      terminates once a minimal model has been found

    OUTPUT:

    - ``newvp`` -- dynamical system on the projective line

    - ``conj`` -- linear fractional transformation which conjugates ``vp`` to ``newvp``

    EXAMPLES::

        sage: PS.<X,Y> = ProjectiveSpace(QQ, 1)
        sage: vp = DynamicalSystem_projective([X^2 + 9*Y^2, X*Y])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import affine_minimal
        sage: affine_minimal(vp, True)
        (
        Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (X : Y) to
                    (X^2 + Y^2 : X*Y)
        ,
        [3 0]
        [0 1]
        )
    """
    from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
    BR = vp.domain().base_ring()
    conj = matrix(BR, 2, 2, 1)
    flag = True
    vp.normalize_coordinates()
    d = vp.degree()

    Affvp = vp.dehomogenize(1)
    R = Affvp.coordinate_ring()
    if R.is_field():
        #want the polynomial ring not the fraction field
        R = R.ring()
    F = R(Affvp[0].numerator())
    G = R(Affvp[0].denominator())
    if R(G.degree()) == 0 or R(F.degree()) == 0:
        raise TypeError("affine minimality is only considered for maps not of the form f or 1/f for a polynomial f")

    z = F.parent().gen(0)
    minG = G
    #If the valuation of a prime in the resultant is small enough, we can say the
    #map is affine minimal at that prime without using the local minimality loop. See
    #Theorem 3.2.2 in [Molnar, M.Sc. thesis]
    if d%2 == 0:
        g = d
    else:
        g = 2*d
    Res = vp.resultant()

    #Some quantities needed for the local minimization loop, but we compute now
    #since the value is constant, so we do not wish to compute in every local loop.
    #See Theorem 3.3.3 in [Molnar, M.Sc thesis]
    H = F - z * minG
    A = AffineSpace(BR, 1, H.parent().variable_name())
    ubRes = DynamicalSystem_affine([H/minG], domain=A).homogenize(1).resultant()
    #Set the primes to check minimality at, if not already prescribed
    if D is None:
        D = ZZ(Res).prime_divisors()

    #Check minimality at all primes in D. If D is all primes dividing
    #Res(minF/minG), this is enough to show whether minF/minG is minimal or not. See
    #Propositions 3.2.1 and 3.3.7 in [Molnar, M.Sc. thesis].
    for p in D:
        while True:
            if Res.valuation(p) < g:
                #The model is minimal at p
                min = True
            else:
                #The model may not be minimal at p.
                newvp,conj = Min(vp, p, ubRes, conj, all_orbits=False)
                if newvp == vp:
                    min = True
                else:
                    vp = newvp
                    Affvp = vp.dehomogenize(1)
                    min = False
            if min:
                #The model is minimal at p
                break
            elif F == Affvp[0].numerator() and G == Affvp[0].denominator():
                #The model is minimal at p
                break
            else:
                #The model is not minimal at p
                flag = False
                if quick:
                    break
        if quick and not flag:
            break

    if quick: #only return whether the model is minimal
        return flag

    if return_transformation:
        return vp, conj
    return vp


def Min(Fun, p, ubRes, conj, all_orbits=False):
    r"""
    Local loop for Affine_minimal, where we check minimality at the prime p.

    First we bound the possible k in our transformations A = zp^k + b.
    See Theorems 3.3.2 and 3.3.3 in [Molnar]_.

    INPUT:

    - ``Fun`` -- a dynamical system on projective space

    - ``p`` - a prime

    - ``ubRes`` -- integer, the upper bound needed for Th. 3.3.3 in [Molnar]_

    - ``conj`` -- a 2x2 matrix keeping track of the conjugation

    - ``all_orbits`` -- boolean; whether or not to use ``==`` in the
      inequalities to find all orbits

    OUTPUT:

    - boolean -- ``True`` if ``Fun`` is minimal at ``p``, ``False`` otherwise

    - a dynamical system on projective space minimal at ``p``

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ, 1)
        sage: f = DynamicalSystem_projective([149*x^2 + 39*x*y + y^2, -8*x^2 + 137*x*y + 33*y^2])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import Min
        sage: Min(f, 3, -27000000, matrix(QQ,[[1, 0],[0, 1]]))
        [
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (157*x^2 + 72*x*y + 3*y^2 : -24*x^2 + 121*x*y + 54*y^2)        ,
        <BLANKLINE>
        [3 1]
        [0 1]
        ]
    """
    d = Fun.degree()
    AffFun = Fun.dehomogenize(1)
    R = AffFun.coordinate_ring()
    if R.is_field():
        #want the polynomial ring not the fraction field
        R = R.ring()
    F = R(AffFun[0].numerator())
    G = R(AffFun[0].denominator())
    dG = G.degree()
    # all_orbits scales bounds for >= and <= if searching for orbits instead of min model
    if dG > (d+1)/2:
        lowerBound = (-2*(G[dG]).valuation(p)/(2*dG - d + 1) + 1).floor() - int(all_orbits)
    else:
        lowerBound = (-2*(F[d]).valuation(p)/(d-1) + 1).floor() - int(all_orbits)
    upperBound = 2*(ubRes.valuation(p)) + int(all_orbits)

    if upperBound < lowerBound:
        # There are no possible transformations to reduce the resultant.
        if not all_orbits:
            return [Fun, conj]
        return []
    # Looping over each possible k, we search for transformations to reduce
    # the resultant of F/G
    all_found = []
    k = lowerBound
    Qb = PolynomialRing(QQ,'b')
    b = Qb.gen(0)
    Q = PolynomialRing(Qb,'z')
    z = Q.gen(0)
    while k <= upperBound:
        A = (p**k)*z + b
        Ft = Q(F(A) - b*G(A))
        Gt = Q((p**k)*G(A))
        Fcoeffs = Ft.coefficients(sparse=False)
        Gcoeffs = Gt.coefficients(sparse=False)
        coeffs = Fcoeffs + Gcoeffs
        RHS = (d + 1) * k / 2
        # If there is some b such that Res(phi^A) < Res(phi), we must have
        # ord_p(c) > RHS for each c in coeffs.
        # Make sure constant coefficients in coeffs satisfy the inequality.
        if all(QQ(c).valuation(p) > RHS - int(all_orbits)
               for c in coeffs if c.degree() == 0):
            # Constant coefficients in coeffs have large enough valuation, so
            # check the rest. We start by checking if simply picking b=0 works.
            if all(c(0).valuation(p) > RHS - int(all_orbits) for c in coeffs):
                # A = z*p^k satisfies the inequalities, and F/G is not minimal
                # "Conjugating by", p,"^", k, "*z +", 0
                newconj = matrix(QQ, 2, 2, [p**k, 0, 0, 1])
                minFun = Fun.conjugate(newconj)
                minFun.normalize_coordinates()
                if not all_orbits:
                    return [minFun, conj*newconj]
                all_found.append([p, k, 0])

            # Otherwise we search if any value of b will work. We start by
            # finding a minimum bound on the valuation of b that is necessary.
            # See Theorem 3.3.5 in [Molnar, M.Sc. thesis].
            bval = max(bCheck(coeff, RHS, p, b) for coeff in coeffs if coeff.degree() > 0)

            # We scale the coefficients in coeffs, so that we may assume
            # ord_p(b) is at least 0
            scaledCoeffs = [coeff(b*(p**bval)) for coeff in coeffs]

            # We now scale the inequalities, ord_p(coeff) > RHS, so that
            # coeff is in ZZ[b]
            scale = QQ(max(coeff.denominator() for coeff in scaledCoeffs))
            normalizedCoeffs = [coeff * scale for coeff in scaledCoeffs]
            scaleRHS = RHS + scale.valuation(p)

            # We now search for integers that satisfy the inequality
            # ord_p(coeff) > RHS. See Lemma 3.3.6 in [Molnar, M.Sc. thesis].
            bound = (scaleRHS + 1 - int(all_orbits)).floor()
            all_blift = blift(normalizedCoeffs, bound, p, k, all_orbits=all_orbits)

            # If bool is true after lifting, we have a solution b, and F/G
            # is not minimal.
            for boolval, sol in all_blift:
                if boolval:
                    #Rescale, conjugate and return new map
                    bsol = QQ(sol * (p**bval))
                    #only add 'minimal orbit element'
                    while bsol.abs() >= p**k:
                        if bsol < 0:
                            bsol += p**k
                        else:
                            bsol -= p**k
                    #"Conjugating by ", p,"^", k, "*z +", bsol
                    newconj = matrix(QQ, 2, 2, [p**k, bsol, 0, 1])
                    minFun = Fun.conjugate(newconj)

                    minFun.normalize_coordinates()
                    if not all_orbits:
                        return [minFun, conj*newconj]
                    if [p,k,bsol] not in all_found:
                        all_found.append([p, k, bsol])
        k = k + 1
    if not all_orbits:
        return [Fun, conj]
    return all_found

###################################################
# algorithms from Hutz-Stoll
###################################################

#modification of Bruin-Molnar for all representatives
def BM_all_minimal(vp, return_transformation=False, D=None):
    r"""
    Determine a representative in each `SL(2,\ZZ)` orbit with minimal
    resultant.

    This function modifies the Bruin-Molnar algorithm ([BM2012]_) to solve
    in the inequalities as ``<=`` instead of ``<``. Among the list of
    solutions is all conjugations that preserve the resultant. From that
    list the `SL(2,\ZZ)` orbits are identified and one representative from
    each orbit is returned. This function assumes that the given model is
    a minimal model.

    INPUT:

    - ``vp`` -- a minimal model of a dynamical system on the projective line

    - ``return_transformation`` -- (default: ``False``) boolean; this
      signals a return of the ``PGL_2`` transformation to conjugate ``vp``
      to the calculated minimal model

    - ``D`` -- a list of primes, in case one only wants to check minimality
      at those specific primes

    OUTPUT:

    List of pairs ``[f, m]`` where ``f`` is a dynamical system and ``m`` is a
    `2 \times 2` matrix.

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([x^3 - 13^2*y^3, x*y^2])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import BM_all_minimal
        sage: BM_all_minimal(f)
        [Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (x^3 - 169*y^3 : x*y^2),
         Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (13*x^3 - y^3 : x*y^2)]

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([x^3 - 6^2*y^3, x*y^2])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import BM_all_minimal
        sage: BM_all_minimal(f, D=[3])
        [Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (x^3 - 36*y^3 : x*y^2),
         Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (3*x^3 - 4*y^3 : x*y^2)]

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([x^3 - 4^2*y^3, x*y^2])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import BM_all_minimal
        sage: cl = BM_all_minimal(f, return_transformation=True)
        sage: all(f.conjugate(m) == g for g, m in cl)
        True
    """
    mp = copy(vp)
    mp.normalize_coordinates()
    BR = mp.domain().base_ring()
    MS = MatrixSpace(QQ, 2)
    M_Id = MS.one()
    F, G = list(mp)  #coordinate polys
    aff_map = mp.dehomogenize(1)
    f, g = aff_map[0].numerator(), aff_map[0].denominator()
    z = aff_map.domain().gen(0)
    Res = mp.resultant()

    ##### because of how the bound is compute in lemma 3.3
    from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
    h = f - z*g
    A = AffineSpace(BR, 1, h.parent().variable_name())
    res = DynamicalSystem_affine([h/g], domain=A).homogenize(1).resultant()

    if D is None:
        D = ZZ(Res).prime_divisors()

    # get the conjugations for each prime independently
    # these are returning (p,k,b) so that the matrix is [p^k,b,0,1]
    all_pM = []
    for p in D:
        # all_orbits used to scale inequalities to equalities
        all_pM.append(Min(mp, p, res, M_Id, all_orbits=True))
        # need the identity for each prime
        if [p, 0, 0] not in all_pM[-1]:
            all_pM[-1].append([p, 0, 0])

    #combine conjugations for all primes
    all_M = [M_Id]
    for prime_data in all_pM:
        #these are (p,k,b) so that the matrix is [p^k,b,0,1]
        new_M = []
        if prime_data:
            p = prime_data[0][0]
            for m in prime_data:
                mat = MS([m[0]**m[1], m[2], 0, 1])
                new_map = mp.conjugate(mat)
                new_map.normalize_coordinates()
                # make sure the resultant didn't change and that it is a different SL(2,ZZ) orbit
                if (mat == M_Id) or (new_map.resultant().valuation(p) == Res.valuation(p)
                                     and mat.det() not in [1,-1]):
                    new_M.append(m)
        if new_M:
            all_M = [m1 * MS([m[0]**m[1], m[2], 0, 1])
                     for m1 in all_M for m in new_M]

    #get all models with same resultant
    all_maps = []
    for M in all_M:
        new_map = mp.conjugate(M)
        new_map.normalize_coordinates()
        if not [new_map, M] in all_maps:
            all_maps.append([new_map, M])

    #Split into conjugacy classes
    #We just keep track of the two matrices that come from
    #the original to get the conjugation that goes between these!!
    classes = []
    for funct, mat in all_maps:
        if not classes:
            classes.append([funct, mat])
        else:
            found = False
            for Func, Mat in classes:
                #get conjugation
                M = mat.inverse() * Mat
                assert funct.conjugate(M) == Func
                if M.det() in [1,-1]:
                    #same SL(2,Z) orbit
                    found = True
                    break
            if found is False:
                classes.append([funct, mat])

    if return_transformation:
        return classes
    else:
        return [funct for funct, matr in classes]

###################################################
# enumerative algorithms from Hutz-Stoll
###################################################

#find minimal model
def HS_minimal(f, return_transformation=False, D=None):
    r"""
    Compute a minimal model for the given projective dynamical system.

    This function implements the algorithm in Hutz-Stoll [HS2018]_.
    A representative with minimal resultant in the conjugacy class
    of ``f`` returned.

    INPUT:

    - ``f`` -- dynamical system on the projective line with minimal resultant

    - ``return_transformation`` -- (default: ``False``) boolean; this
      signals a return of the `PGL_2` transformation to conjugate
      this map to the calculated models

    - ``D`` -- a list of primes, in case one only wants to check minimality
      at those specific primes

    OUTPUT:

    - a dynamical system
    - (optional) a `2 \times 2` matrix

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([x^3 - 6^2*y^3, x^2*y])
        sage: m = matrix(QQ,2,2,[5,1,0,1])
        sage: g = f.conjugate(m)
        sage: g.normalize_coordinates()
        sage: g.resultant().factor()
        2^4 * 3^4 * 5^12
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import HS_minimal
        sage: HS_minimal(g).resultant().factor()
        2^4 * 3^4
        sage: HS_minimal(g, D=[2]).resultant().factor()
        2^4 * 3^4 * 5^12
        sage: F,m = HS_minimal(g, return_transformation=True)
        sage: g.conjugate(m) == F
        True
    """
    F = copy(f)
    d = F.degree()
    F.normalize_coordinates()
    MS = MatrixSpace(ZZ, 2, 2)
    m = MS.one()
    res = ZZ(F.resultant())
    if D is None:
        D = res.prime_divisors()

    # minimize for each prime
    for p in D:
        vp = res.valuation(p)
        minimal = False
        while not minimal:
            if (d % 2 == 0 and vp < d) or (d % 2 == 1 and vp < 2 * d):
                # must be minimal
                minimal = True
                break
            minimal = True
            t = MS([1, 0, 0, p])
            F1 = F.conjugate(t)
            F1.normalize_coordinates()
            res1 = F1.resultant()
            vp1 = res1.valuation(p)
            if vp1 < vp: # check if smaller
                F = F1
                vp = vp1
                m = m * t # keep track of conjugation
                minimal = False
            else:
                # still search for smaller
                for b in range(p):
                    t = matrix(ZZ,2,2,[p, b, 0, 1])
                    F1 = F.conjugate(t)
                    F1.normalize_coordinates()
                    res1 = ZZ(F1.resultant())
                    vp1 = res1.valuation(p)
                    if vp1 < vp: # check if smaller
                        F = F1
                        m = m * t # keep track of transformation
                        minimal = False
                        vp = vp1
                        break # exit for loop
    if return_transformation:
        return F, m
    return F

#find all representatives of orbits for one prime
def HS_all_minimal_p(p, f, m=None, return_transformation=False):
    r"""
    Find a representative in each distinct `SL(2,\ZZ)` orbit with
    minimal `p`-resultant.

    This function implements the algorithm in Hutz-Stoll [HS2018]_.
    A representatives in each distinct `SL(2,\ZZ)` orbit with minimal
    valuation with respect to the prime ``p`` is returned. The input
    ``f`` must have minimal resultant in its conjugacy class.

    INPUT:

    - ``p`` -- a prime

    - ``f`` -- dynamical system on the projective line with minimal resultant

    - ``m`` -- (optional) `2 \times 2` matrix associated with ``f``

    - ``return_transformation`` -- (default: ``False``) boolean; this
      signals a return of the ``PGL_2`` transformation to conjugate ``vp``
      to the calculated minimal model

    OUTPUT:

    List of pairs ``[f, m]`` where ``f`` is a dynamical system and ``m`` is a
    `2 \times 2` matrix.

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([x^5 - 6^4*y^5, x^2*y^3])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import HS_all_minimal_p
        sage: HS_all_minimal_p(2, f)
        [Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (x^5 - 1296*y^5 : x^2*y^3),
         Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (4*x^5 - 162*y^5 : x^2*y^3)]
        sage: cl = HS_all_minimal_p(2, f, return_transformation=True)
        sage: all(f.conjugate(m) == g for g, m in cl)
        True
    """
    count = 0
    prev = 0 # no exclusions
    F = copy(f)
    res = ZZ(F.resultant())
    vp = res.valuation(p)
    MS = MatrixSpace(ZZ, 2)
    if m is None:
        m = MS.one()
    if f.degree() % 2 == 0 or vp == 0:
        # there is only one orbit for even degree
        # nothing to do if the prime doesn't divide the resultant
        if return_transformation:
            return [[f, m]]
        else:
            return [f]
    to_do = [[F, m, prev]] # repns left to check
    reps = [[F, m]] # orbit representatives for f
    while to_do:
        F, m, prev = to_do.pop()
        # there are at most two directions preserving the resultant
        if prev == 0:
            count = 0
        else:
            count = 1
        if prev != 2: # [p,a,0,1]
            t = MS([1, 0, 0, p])
            F1 = F.conjugate(t)
            F1.normalize_coordinates()
            res1 = ZZ(F1.resultant())
            vp1 = res1.valuation(p)
            if vp1 == vp:
                count += 1
                # we have a new representative
                reps.append([F1, m*t])
                # need to check if it has any neighbors
                to_do.append([F1, m*t, 1])
        for b in range(p):
            if not (b == 0 and prev == 1):
                t = MS([p, b, 0, 1])
                F1 = F.conjugate(t)
                F1.normalize_coordinates()
                res1 = ZZ(F1.resultant())
                vp1 = res1.valuation(p)
                if vp1 == vp:
                    count += 1
                    # we have a new representative
                    reps.append([F1, m*t])
                    # need to check if it has any neighbors
                    to_do.append([F1, m*t, 2])
            if count >= 2: # at most two neighbors
                break

    if return_transformation:
        return reps
    else:
        return [funct for funct, matr in reps]

#find all representatives of orbits
def HS_all_minimal(f, return_transformation=False, D=None):
    r"""
    Determine a representative in each `SL(2,\ZZ)` orbit with minimal resultant.

    This function implements the algorithm in Hutz-Stoll [HS2018]_.
    A representative in each distinct `SL(2,\ZZ)` orbit is returned.
    The input ``f`` must have minimal resultant in its conjugacy class.

    INPUT:

    - ``f`` -- dynamical system on the projective line with minimal resultant

    - ``return_transformation`` -- (default: ``False``) boolean; this
      signals a return of the ``PGL_2`` transformation to conjugate ``vp``
      to the calculated minimal model

    - ``D`` -- a list of primes, in case one only wants to check minimality
      at those specific primes

    OUTPUT:

    List of pairs ``[f, m]``, where ``f`` is a dynamical system and ``m``
    is a `2 \times 2` matrix.

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([x^3 - 6^2*y^3, x^2*y])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import HS_all_minimal
        sage: HS_all_minimal(f)
        [Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (x^3 - 36*y^3 : x^2*y),
         Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (9*x^3 - 12*y^3 : 9*x^2*y),
         Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (4*x^3 - 18*y^3 : 4*x^2*y),
         Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (36*x^3 - 6*y^3 : 36*x^2*y)]
        sage: HS_all_minimal(f, D=[3])
        [Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (x^3 - 36*y^3 : x^2*y),
         Dynamical System of Projective Space of dimension 1 over Rational Field
           Defn: Defined on coordinates by sending (x : y) to
                 (9*x^3 - 12*y^3 : 9*x^2*y)]

    ::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([x^3 - 6^2*y^3, x*y^2])
        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import HS_all_minimal
        sage: cl = HS_all_minimal(f, return_transformation=True)
        sage: all(f.conjugate(m) == g for g, m in cl)
        True
    """
    MS = MatrixSpace(ZZ, 2)
    m = MS.one()
    F = copy(f)
    F.normalize_coordinates()
    if F.degree() == 1:
        raise ValueError("function must be degree at least 2")
    if f.degree() % 2 == 0:
        #there is only one orbit for even degree
        if return_transformation:
            return [[f, m]]
        else:
            return [f]
    if D is None:
        res = ZZ(F.resultant())
        D = res.prime_divisors()
    M = [[F, m]]
    for p in D:
        # get p-orbits
        Mp = HS_all_minimal_p(p, F, m, return_transformation=True)
        # combine with previous orbits representatives
        M = [[g.conjugate(t), t*s] for g,s in M for G,t in Mp]

    if return_transformation:
        return M
    else:
        return [funct for funct, matr in M]

#######################
#functionality for smallest coefficients
#
# Ben Hutz July 2018
#####################################3

def get_bound_dynamical(F, f, m=1, dynatomic=True, prec=53, emb=None):
    """
    The hyperbolic distance from `j` which must contain the smallest map.

    This defines the maximum possible distance from `j` to the `z_0` covariant
    of the associated binary form `F` in the hyperbolic 3-space
    for which the map `f`` could have smaller coefficients.

    INPUT:

    - ``F`` -- binary form of degree at least 3 with no multiple roots associated
      to ``f``

    - ``f`` -- a dynamical system on `P^1`

    - ``m`` - positive integer. the period used to create ``F``

    - ``dynatomic`` -- boolean. whether ``F`` is the periodic points or the
      formal periodic points of period ``m`` for ``f``

    - ``prec``-- positive integer. precision to use in CC

    - ``emb`` -- embedding into CC

    OUTPUT: a positive real number

    EXAMPLES::

        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import get_bound_dynamical
        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([50*x^2 + 795*x*y + 2120*y^2, 265*x^2 + 106*y^2])
        sage: get_bound_dynamical(f.dynatomic_polynomial(1), f)
        35.5546923182219
    """
    def coshdelta(z):
        #The cosh of the hyperbolic distance from z = t+uj to j
        return (z.norm() + 1)/(2*z.imag())
    if F.base_ring() != ComplexField(prec=prec):
        if emb is None:
            compF = F.change_ring(ComplexField(prec=prec))
        else:
            compF = F.change_ring(emb)
    else:
        compF = F
    n = F.degree()

    z0F, thetaF = covariant_z0(compF, prec=prec, emb=emb)
    d = f.degree()
    hF = e**f.global_height(prec=prec)
    #get precomputed constants C,k
    if m == 1:
        C = 4*d+2
        k = 2
    else:
        Ck_values = {(False, 2, 2): (322, 6), (False, 2, 3): (385034, 14),
                     (False, 2, 4): (4088003923454, 30), (False, 3, 2): (18044, 8),
                     (False, 4, 2): (1761410, 10), (False, 5, 2): (269283820, 12),
                     (True, 2, 2): (43, 4), (True, 2, 3): (106459, 12),
                     (True, 2, 4): (39216735905, 24), (True, 3, 2): (1604, 6),
                     (True, 4, 2): (114675, 8), (True, 5, 2): (14158456, 10)}
        try:
            C, k = Ck_values[(dynatomic,d,m)]
        except KeyError:
            raise ValueError("constants not computed for this (m,d) pair")
    if n == 2 and d == 2:
        #bound with epsilonF = 1
        bound = 2*((2*C*(hF**k))/(thetaF))
    else:
        bound = cosh(epsinv(F, (2**(n-1))*C*(hF**k)/thetaF, prec=prec))
    return bound


def smallest_dynamical(f, dynatomic=True, start_n=1, prec=53, emb=None, algorithm='HS', check_minimal=True):
    r"""
    Determine the poly with smallest coefficients in `SL(2,\ZZ)` orbit of ``F``

    Smallest is in the sense of global height.
    The method is the algorithm in Hutz-Stoll [HS2018]_.
    A binary form defining the periodic points is associated to ``f``.
    From this polynomial a bound on the search space can be determined.

    ``f`` should already be a minimal model or finding the orbit
    representatives may give wrong results.

    INPUT:

    - ``f`` -- a dynamical system on `P^1`

    - ``dynatomic`` -- boolean. whether ``F`` is the periodic points or the
      formal periodic points of period ``m`` for ``f``

    - ``start_n`` - positive integer. the period used to start trying to
      create associate binary form ``F``

    - ``prec``-- positive integer. precision to use in CC

    - ``emb`` -- embedding into CC

    - ``algorithm`` -- (optional) string; either ``'BM'`` for the Bruin-Molnar
      algorithm or ``'HS'`` for the Hutz-Stoll algorithm. If not specified,
      properties of the map are utilized to choose how to compute minimal
      orbit representatives

    - ``check_minimal`` -- (default: True), boolean, whether to check
      if this map is a minimal model

    OUTPUT: pair [dynamical system, matrix]

    EXAMPLES::

        sage: from sage.dynamics.arithmetic_dynamics.endPN_minimal_model import smallest_dynamical
        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem([50*x^2 + 795*x*y + 2120*y^2, 265*x^2 + 106*y^2])
        sage: smallest_dynamical(f)  #long time
        [
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (-480*x^2 - 1125*x*y + 1578*y^2 : 265*x^2 + 1060*x*y + 1166*y^2),
        <BLANKLINE>
        [1 2]
        [0 1]
        ]
    """
    def insert_item(pts, item, index):
        # binary insertion to maintain list of points left to consider
        N = len(pts)
        if N == 0:
          return [item]
        elif N == 1:
            if item[index] > pts[0][index]:
                pts.insert(0,item)
            else:
                pts.append(item)
            return pts
        else: # binary insertion
            left = 1
            right = N
            mid = (left + right) // 2  # these are ints so this is .floor()
            if item[index] > pts[mid][index]: # item goes into first half
                return insert_item(pts[:mid], item, index) + pts[mid:N]
            else: # item goes into second half
                return pts[:mid] + insert_item(pts[mid:N], item, index)

    def coshdelta(z):
        # The cosh of the hyperbolic distance from z = t+uj to j
        return (z.norm() + 1)/(2*z.imag())

    # can't be smaller if height 0
    f.normalize_coordinates()
    if f.global_height(prec=prec) == 0:
        return [f, matrix(ZZ,2,2,[1,0,0,1])]
    all_min = f.all_minimal_models(return_transformation=True, algorithm=algorithm, check_minimal=check_minimal)

    current_min = None
    current_size = None
    # search for minimum over all orbits
    for g,M in all_min:
        PS = g.domain()
        CR = PS.coordinate_ring()
        x,y = CR.gens()
        n = start_n # sometimes you get a problem later with 0,infty as roots
        if dynatomic:
            pts_poly = g.dynatomic_polynomial(n)
        else:
            gn = g.nth_iterate_map(n)
            pts_poly = y*gn[0] - x*gn[1]
        d = ZZ(pts_poly.degree())
        max_mult = max([ex for p,ex in pts_poly.factor()])
        while ((d < 3) or (max_mult >= d/2) and (n < 5)):
            n = n+1
            if dynatomic:
                pts_poly = g.dynatomic_polynomial(n)
            else:
                gn = g.nth_iterate_map(n)
                pts_poly = y*gn[0] - x*gn[1]
            d = ZZ(pts_poly.degree())
            max_mult = max([ex for p,ex in pts_poly.factor()])
        assert(n<=4), "n > 4, failed to find usable poly"

        R = get_bound_dynamical(pts_poly, g, m=n, dynatomic=dynatomic, prec=prec, emb=emb)
        # search starts in fundamental domain
        G,MG = pts_poly.reduced_form(prec=prec, emb=emb, smallest_coeffs=False)
        red_g = f.conjugate(M*MG)
        if G != pts_poly:
            R2 = get_bound_dynamical(G, red_g, m=n, dynatomic=dynatomic, prec=prec, emb=emb)
            if R2 < R:
                # use the better bound
                R = R2
        red_g.normalize_coordinates()
        if red_g.global_height(prec=prec) == 0:
            return [red_g, M*MG]

        # height
        if current_size is None:
            current_size = e**red_g.global_height(prec=prec)
        v0, th = covariant_z0(G, prec=prec, emb=emb)
        rep = 2*CC.gen(0)
        from math import isnan
        if isnan(v0.abs()):
            raise ValueError("invalid covariant: %s"%v0)

        # get orbit
        S = matrix(ZZ,2,2,[0,-1,1,0])
        T = matrix(ZZ,2,2,[1,1,0,1])
        TI = matrix(ZZ,2,2,[1,-1,0,1])

        count = 0
        pts = [[G, red_g, v0, rep, M*MG, coshdelta(v0), 0]]  # label - 0:None, 1:S, 2:T, 3:T^(-1)
        if current_min is None:
            current_min = [G, red_g, v0, rep, M*MG, coshdelta(v0)]
        while pts != []:
            G, g, v, rep, M, D, label = pts.pop()
            # apply ST and keep z, Sz
            if D > R:
                break #all remaining pts are too far away
            # check if it is smaller. If so, we can improve the bound
            count += 1
            new_size = e**g.global_height(prec=prec)
            if new_size < current_size:
                current_min = [G ,g, v, rep, M, coshdelta(v)]
                current_size = new_size
                if new_size == 1: # early exit
                    return [current_min[1], current_min[4]]
                new_R = get_bound_dynamical(G, g, m=n, dynatomic=dynatomic, prec=prec, emb=emb)
                if new_R < R:
                    R = new_R

            # add new points to check
            if label != 1 and min((rep+1).norm(), (rep-1).norm()) >= 1: # don't undo S
                # the 2nd condition is equivalent to |\Re(-1/rep)| <= 1/2
                # this means that rep can have resulted from an inversion step in
                # the shift-and-invert procedure, so don't invert

                # do inversion
                z = -1/v
                new_pt = [G.subs({x:-y, y:x}), g.conjugate(S), z, -1/rep, M*S, coshdelta(z), 1]
                pts = insert_item(pts, new_pt, 5)
            if label != 3:  # don't undo T on g
                # do right shift
                z = v-1
                new_pt = [G.subs({x:x+y}), g.conjugate(TI), z, rep-1, M*TI, coshdelta(z), 2]
                pts = insert_item(pts, new_pt, 5)
            if label != 2:  # don't undo TI on g
                # do left shift
                z = v+1
                new_pt = [G.subs({x:x-y}), g.conjugate(T), z, rep+1, M*T, coshdelta(z), 3]
                pts = insert_item(pts, new_pt, 5)

    return [current_min[1], current_min[4]]

