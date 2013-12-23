r"""
Sage functions to compute minimal models of rational functions
under the conjugation action of `PGL_2(QQ)`.

AUTHORS:

 - Alex Molnar (May 22, 2012)

- Brian Stout, Ben Hutz (Nov 2013): Modified code to use projective morphism functionality
  so that it can be included in Sage.

REFERENCES:

.. [Bruin-Molnar] N. Bruin and A. Molnar, Minimal models for rational
   functions in a dynamical setting
   LMS Journal of Computation and Mathematics, Volume 15 (2012), pp 400-417.

.. [Molnar] A. Molnar, Fractional Linear Minimal Models of Rational Functions,
   M.Sc. Thesis.
"""

#*****************************************************************************
#       Copyright (C) 2012 Alexander Molnar
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.homset                            import End
from copy                                              import copy
from sage.matrix.constructor                           import matrix
from sage.rings.arith                                  import gcd
from sage.rings.finite_rings.integer_mod_ring          import Zmod
from sage.rings.integer_ring                           import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field                         import QQ
from sage.schemes.affine.affine_space                  import AffineSpace

def bCheck(c,v,p,b):
    r"""
    Compute a lower bound on the value of ``b`` needed, for a transformation
    `A(z) = z*p^k + b` to satisfy `ord_p(Res(\phi^A)) < ord_p(Res(\phi))` for a
    rational map `\phi`. See Theorem 3.3.5 in [Molnar, M.Sc. thesis].

    INPUT:

    - ``c`` -- a list of polynomials in `b`. See v for their use.

    - ``v`` -- a list of rational numbers, where we are considering the inequalities
           `ord_p(c[i]) > v[i]`.

    - ``p`` -- a prime.

    - ``b`` -- local variable.

    OUTPUT:

    - ``bval`` -- Integer, lower bound in Theorem 3.3.5

    EXAMPLES::

        sage: R.<b> = PolynomialRing(QQ)
        sage: from sage.schemes.projective.endPN_minimal_model import bCheck
        sage: bCheck(11664*b^2 + 70227*b + 76059, 15/2, 3, b)
        -1
    """
    val = (v+1).floor()
    deg = c.degree()
    coeffs = c.coeffs()
    lcoeff = coeffs[deg]; coeffs.remove(lcoeff)
    check1 = [(coeffs[i].valuation(p) - lcoeff.valuation(p))/(deg - i) for i in range(0,len(coeffs)) if coeffs[i] != 0]
    check2 = (val - lcoeff.valuation(p))/deg
    check1.append(check2)
    bval = min(check1)
    return (bval).ceil()


def scale(c,v,p):
    r"""
    Given an integral polynomial ``c``, we can write `c = p^i*c'`, where ``p`` does not
    divide ``c``. Returns ``c'`` and `v - i` where `i` is the smallest valuation of the
    coefficients of `c`.

    INPUT:

    - ``c`` -- an integer polynomial

    - ``v`` -- an integer - the bound on the exponent from blift

    - ``p`` -- a prime

    OUTPUT:

    - Boolean -- the new exponent bound is 0 or negative

    - the scaled integer polynomial

    - an integer the new exponent bound

    EXAMPLES::

        sage: R.<b> = PolynomialRing(QQ)
        sage: from sage.schemes.projective.endPN_minimal_model import scale
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
    return [flag,c,v]


def blift(LF,Li,p,S=None):
    r"""
    Search for a solution to the given list of inequalities. If found, lift the solution to
    an appropriate valuation. See Lemma 3.3.6 in [Molnar, M.Sc. thesis]

    INPUT:

    - ``LF`` -- a list of integer polynomials in one variable (the normalized coefficients)

    - ``Li`` -- an integer, the bound on coefficients

    - ``p`` -- a prime

    OUTPUT:

    - Boolean -- whether or not the lift is successful

    - integer -- the lift

    EXAMPLES::

        sage: R.<b> = PolynomialRing(QQ)
        sage: from sage.schemes.projective.endPN_minimal_model import blift
        sage: blift([8*b^3 + 12*b^2 + 6*b + 1, 48*b^2 + 483*b + 117, 72*b + 1341, -24*b^2 + 411*b + 99, -144*b + 1233, -216*b], 2, 3)
        (True, 4)
    """

    P = LF[0].parent()
    #Determine which inequalities are trivial, and scale the rest, so that we only lift
    #as many times as needed.
    keepScaledIneqs = [scale(P(coeff),Li,p) for coeff in LF if coeff != 0]
    keptVals = [i[2] for i in keepScaledIneqs if i[0]]
    if keptVals != []:
        #Determine the valuation to lift until.
        liftval = max(keptVals)
    else:
        #All inequalities are satisfied.
        return True,1
    if S is None:
        S = PolynomialRing(Zmod(p),'b')
    keptScaledIneqs = [S(i[1]) for i in keepScaledIneqs if i[0]]
    #We need a solution for each polynomial on the left hand side of the inequalities,
    #so we need only find a solution for their gcd.
    g = gcd(keptScaledIneqs)
    rts = g.roots(multiplicities=False)
    for r in rts:
        #Recursively try to lift each root
        r_initial=QQ(r)
        newInput = P([r_initial,p])
        LG = [F(newInput) for F in LF]
        lift,lifted = blift(LG,Li,p,S=S)
        if lift:
            #Lift successful.
            return True,r_initial+ p*lifted
    #Lift non successful.
    return False,0


def affine_minimal(vp, return_transformation = False,D=None, quick = False):
    r"""
    Given vp a scheme morphisms on the projective line over the rationals,
    this procedure determines if `\phi` is minimal. In particular, it determines
    if the map is affine minimal, which is enough to decide if it is minimal
    or not. See Proposition 2.10 in [Bruin-Molnar].

    INPUT:

    - ``vp`` -- scheme morphism on the projective line.

    - ``D`` -- a list of primes, in case one only wants to check minimality
               at those specific primes.

    - ``return_transformation`` -- a boolean value, default value True. This
                                signals a return of the ``PGL_2`` transformation
                                to conjugate ``vp`` to the calculated minimal
                                model. default: False

    - ``quick`` -- a boolean value. If true the algorithm terminates once
                   algorithm determines F/G is not minimal, otherwise algorithm
                   only terminates once a minimal model has been found.

    OUTPUT:

    - ``newvp`` -- scheme morphism on the projective line.

    - ``conj`` -- linear fractional transformation which conjugates ``vp`` to ``newvp``

    EXAMPLES::
        sage: PS.<X,Y> = ProjectiveSpace(QQ,1)
        sage: H = Hom(PS,PS)
        sage: vp = H([X^2+9*Y^2,X*Y])
        sage: from sage.schemes.projective.endPN_minimal_model import affine_minimal
        sage: affine_minimal(vp,True)
        (
        Scheme endomorphism of Projective Space of dimension 1 over Rational
        Field
          Defn: Defined on coordinates by sending (X : Y) to
                (X^2 + Y^2 : X*Y)
        ,
        [3 0]
        [0 1]
        )
    """
    BR=vp.domain().base_ring()
    conj=matrix(BR,2,2,1)
    flag = True
    d = vp.degree()

    vp.normalize_coordinates();
    Affvp=vp.dehomogenize(1)
    R=Affvp.coordinate_ring()
    if R.is_field():
        #want the polynomial ring not the fraction field
        R=R.ring()
    F=R(Affvp[0].numerator())
    G=R(Affvp[0].denominator())
    if G.degree()==0 or F.degree()==0:
        raise TypeError("Affine minimality is only considered for maps not of the form f or 1/f for a polynomial f.")
    z=F.parent().gen(0)
    minF,minG = F,G
    #If the valuation of a prime in the resultant is small enough, we can say the
    #map is affine minimal at that prime without using the local minimality loop. See
    #Theorem 3.2.2 in [Molnar, M.Sc. thesis]
    if d%2 == 0:
        g = d
    else:
        g = 2*d
    Res = vp.resultant();

    #Some quantities needed for the local minimization loop, but we compute now
    #since the value is constant, so we do not wish to compute in every local loop.
    #See Theorem 3.3.3 in [Molnar, M.Sc thesis]

    H=F-z*minG
    d1=F.degree()
    A=AffineSpace(BR,1,H.parent().variable_name())
    end_ring=End(A)

    ubRes = end_ring([H/minG]).homogenize(1).resultant()
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
                newvp,conj = Min(vp,p,ubRes,conj)
                if newvp==vp:
                    min=True
                else:
                    vp=newvp
                    Affvp=vp.dehomogenize(1)
                    min=False
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

def Min(Fun,p,ubRes,conj):
    r"""
    Local loop for Affine_minimal, where we check minimality at the prime p.
    First we bound the possible k in our transformations A = zp^k + b. See Theorems
    3.3.2 and 3.3.3 in [Molnar, M.Sc. thesis].

    INPUT:

    - ``Fun`` -- a projective space morphisms

    - ``p`` - a prime.

    - ``ubRes`` -- integer -- The upper bound needed for Theorem 3.3.3 in [Molnar, M.Sc. thesis].

    - ``conj`` -- a 2x2 matrix keeping track of the conjugation

    OUTPUT:

    - Boolean -- True if Fun is minimal at ``p``, False otherwise

    - a projective morphism minimal at ``p``

    EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: H = End(P)
        sage: f = H([149*x^2 + 39*x*y + y^2, -8*x^2 + 137*x*y + 33*y^2])
        sage: from sage.schemes.projective.endPN_minimal_model import Min
        sage: Min(f, 3, -27000000, matrix(QQ,[[1, 0],[0, 1]]))
        (
        Scheme endomorphism of Projective Space of dimension 1 over Rational
        Field
          Defn: Defined on coordinates by sending (x : y) to
                (181*x^2 + 313*x*y + 81*y^2 : -24*x^2 + 73*x*y + 151*y^2)
        ,
        [3 4]
        [0 1]
        )
    """

    d=Fun.degree()
    AffFun=Fun.dehomogenize(1)
    R=AffFun.coordinate_ring()
    if R.is_field():
        #want the polynomial ring not the fraction field
        R=R.ring()
    F=R(AffFun[0].numerator())
    G=R(AffFun[0].denominator())
    dG = G.degree()
    if dG > (d+1)/2:
        lowerBound = (-2*(G[dG]).valuation(p)/(2*dG - d + 1) + 1).floor()
    else:
        lowerBound = (-2*(F[d]).valuation(p)/(d-1) + 1).floor()
    upperBound = 2*(ubRes.valuation(p))

    if upperBound < lowerBound:
        #There are no possible transformations to reduce the resultant.
        return Fun,conj
    else:
        #Looping over each possible k, we search for transformations to reduce the
        #resultant of F/G
        k = lowerBound
        Qb = PolynomialRing(QQ,'b')
        b=Qb.gen(0)
        Q = PolynomialRing(Qb,'z')
        z=Q.gen(0)
        while k <= upperBound:
            A = (p**k)*z + b
            Ft = Q(F(A) - b*G(A))
            Gt = Q((p**k)*G(A))
            Fcoeffs = Ft.coeffs()
            Gcoeffs = Gt.coeffs()
            coeffs = Fcoeffs + Gcoeffs
            RHS = (d + 1)*k/2
            #If there is some b such that Res(phi^A) < Res(phi), we must have ord_p(c) >
            #RHS for each c in coeffs.

            #Make sure constant coefficients in coeffs satisfy the inequality.

            if all( QQ(c).valuation(p) > RHS for c in coeffs if c.degree() ==0 ):
                #Constant coefficients in coeffs have large enough valuation, so check
                #the rest.
                #We start by checking if simply picking b=0 works
                if all(c(0).valuation(p) > RHS for c in coeffs):
                    #A = z*p^k satisfies the inequalities, and F/G is not minimal

                    #"Conjugating by", p,"^", k, "*z +", 0
                    newconj=matrix(QQ,2,2,[p**k,0,0,1])
                    minFun=Fun.conjugate(newconj)
                    conj=conj*newconj

                    minFun.normalize_coordinates()
                    return minFun, conj

                #Otherwise we search if any value of b will work. We start by finding a
                #minimum bound on the valuation of b that is necessary. See Theorem 3.3.5
                #in [Molnar, M.Sc. thesis].
                bval = max([bCheck(coeff,RHS,p,b) for coeff in coeffs if coeff.degree() > 0])

                #We scale the coefficients in coeffs, so that we may assume ord_p(b) is
                #at least 0
                scaledCoeffs = [coeff(b*(p**bval)) for coeff in coeffs]

                #We now scale the inequalities, ord_p(coeff) > RHS, so that coeff is in
                #ZZ[b]
                scale = QQ(max([coeff.denominator() for coeff in scaledCoeffs]))
                normalizedCoeffs = [coeff*scale for coeff in scaledCoeffs]
                scaleRHS = RHS + scale.valuation(p)

                #We now search for integers that satisfy the inequality ord_p(coeff) >
                #RHS. See Lemma 3.3.6 in [Molnar, M.Sc. thesis].
                bound = (scaleRHS+1).floor()
                bool,sol = blift(normalizedCoeffs,bound,p)

                #If bool is true after lifting, we have a solution b, and F/G is not
                #minimal.
                if bool:
                    #Rescale, conjugate and return new map
                    bsol = QQ(sol*(p**bval))
                    #"Conjugating by ", p,"^", k, "*z +", bsol
                    newconj=matrix(QQ,2,2,[p**k,bsol,0,1])
                    minFun=Fun.conjugate(newconj)
                    conj=conj*newconj

                    minFun.normalize_coordinates()
                    return minFun, conj
            k = k + 1
        return Fun, conj
