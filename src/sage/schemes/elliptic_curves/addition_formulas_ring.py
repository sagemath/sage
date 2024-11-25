
def add(E, P, Q):
    r"""
    Addition formulas for elliptic curves over general rings
    with trivial Picard group.

    REFERENCES:

    These formulas were derived by Bosma and Lenstra [BL1995]_,
    with corrections by Best [Best2021]_ (Appendix A).

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.addition_formulas_ring import add
        sage: M = Zmod(13*17*19)
        sage: R.<U,V> = M[]
        sage: S.<u,v> = R.quotient(U*V - 17)
        sage: E = EllipticCurve(S, [1,2,3,4,5])
        sage: P = E(817, 13, 19)
        sage: Q = E(425, 123, 17)
        sage: PQ1, PQ2 = add(E, P, Q)
        sage: PQ1
        (1188, 1674, 540)
        sage: PQ2
        (582, 2347, 1028)
        sage: E(PQ1) == E(PQ2)
        True
    """
    a1, a2, a3, a4, a6 = E.a_invariants()
    b2, b4, b6, b8 = E.b_invariants()

    assert P in E
    assert Q in E
    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q

    #TODO: I've made a half-hearted attempt at simplifying the formulas
    # by caching common subexpressions. This could almost certainly be
    # sped up significantly with some more serious optimization effort.

    XYdif = X1*Y2 - X2*Y1
    XYsum = X1*Y2 + X2*Y1
    XZdif = X1*Z2 - X2*Z1
    XZsum = X1*Z2 + X2*Z1
    YZdif = Y1*Z2 - Y2*Z1
    YZsum = Y1*Z2 + Y2*Z1

    a1sq, a2sq, a3sq, a4sq = (a**2 for a in (a1, a2, a3, a4))

    X31 = XYdif*YZsum+XZdif*Y1*Y2+a1*X1*X2*YZdif+a1*XYdif*XZsum-a2*X1*X2*XZdif+a3*XYdif*Z1*Z2+a3*XZdif*YZsum-a4*XZsum*XZdif-3*a6*XZdif*Z1*Z2

    Y31 = -3*X1*X2*XYdif-Y1*Y2*YZdif-2*a1*XZdif*Y1*Y2+(a1sq+3*a2)*X1*X2*YZdif-(a1sq+a2)*XYsum*XZdif+(a1*a2-3*a3)*X1*X2*XZdif-(2*a1*a3+a4)*XYdif*Z1*Z2+a4*XZsum*YZdif+(a1*a4-a2*a3)*XZsum*XZdif+(a3sq+3*a6)*YZdif*Z1*Z2+(3*a1*a6-a3*a4)*XZdif*Z1*Z2

    Z31 = 3*X1*X2*XZdif-YZsum*YZdif+a1*XYdif*Z1*Z2-a1*XZdif*YZsum+a2*XZsum*XZdif-a3*YZdif*Z1*Z2+a4*XZdif*Z1*Z2

    yield (X31, Y31, Z31)

    X32 = Y1*Y2*XYsum+a1*(2*X1*Y2+X2*Y1)*X2*Y1+a1sq*X1*X2**2*Y1-a2*X1*X2*XYsum-a1*a2*X1**2*X2**2+a3*X2*Y1*(YZsum+Y2*Z1)+a1*a3*X1*X2*YZdif-a1*a3*XYsum*XZdif-a4*X1*X2*YZsum-a4*XYsum*XZsum-a1sq*a3*X1**2*X2*Z2-a1*a4*X1*X2*(X1*Z2+XZsum)-a2*a3*X1*X2**2*Z1-a3sq*X1*Z2*(Y2*Z1+YZsum)-3*a6*XYsum*Z1*Z2-3*a6*XZsum*YZsum-a1*a3sq*X1*Z2*(XZsum+X2*Z1)-3*a1*a6*X1*Z2*(XZsum+X2*Z1)-a3*a4*(X1*Z2+XZsum)*X2*Z1-b8*YZsum*Z1*Z2-a1*b8*X1*Z1*Z2**2-a3**3*XZsum*Z1*Z2-3*a3*a6*(XZsum+X2*Z1)*Z1*Z2-a3*b8*Z1**2*Z2**2

    Y32 = Y1**2*Y2**2+a1*X2*Y1**2*Y2+(a1*a2-3*a3)*X1*X2**2*Y1+a3*Y1**2*Y2*Z2-(a2sq-3*a4)*X1**2*X2**2+(a1*a4-a2*a3)*(2*X1*Z2+X2*Z1)*X2*Y1+(a1sq*a4-2*a1*a2*a3+3*a3sq)*X1**2*X2*Z2-(a2*a4-9*a6)*X1*X2*XZsum+(3*a1*a6-a3*a4)*(XZsum+X2*Z1)*Y1*Z2+(3*a1sq*a6-2*a1*a3*a4+a2*a3sq+3*a2*a6-a4sq)*X1*Z2*(XZsum+X2*Z1)+(3*a2*a6-a4sq)*X2*Z1*(2*X1*Z2+Z1*X2)+(a1**3*a6-a1sq*a3*a4+a1*a2*a3sq-a1*a4sq+4*a1*a2*a6-a3**3-3*a3*a6)*Y1*Z1*Z2**2+(a1**4*a6-a1**3*a3*a4+5*a1sq*a2*a6+a1sq*a2*a3sq-a1*a2*a3*a4-a1*a3**3-3*a1*a3*a6-a1sq*a4sq+a2sq*a3sq-a2*a4sq+4*a2sq*a6-a3 **2*a4-3*a4*a6)*X1*Z1*Z2**2+(a1sq*a2*a6-a1*a2*a3*a4+3*a1*a3*a6+a2sq*a3sq-a2*a4sq+4*a2sq*a6-2*a3sq*a4-3*a4*a6)*X2*Z1**2*Z2+(a1**3*a3*a6-a1sq*a3sq*a4+a1sq*a4*a6+a1*a2*a3**3+4*a1*a2*a3*a6-2*a1*a3*a4sq+a2*a3sq*a4+4*a2*a4*a6-a3**4-6*a3 **2*a6-a4**3-9*a6**2)*Z1**2*Z2**2

    Z32 = 3*X1*X2*XYsum+Y1*Y2*YZsum+3*a1*X1**2*X2**2+a1*(2*X1*Y2+Y1*X2)*Y1*Z2+a1sq*X1*Z2*(2*X2*Y1+X1*Y2)+a2*X1*X2*YZsum+a2*XYsum*XZsum+a1**3*X1**2*X2*Z2+a1*a2*X1*X2*(2*X1*Z2+X2*Z1)+3*a3*X1*X2**2*Z1+a3*Y1*Z2*(YZsum+Y2*Z1)+2*a1*a3*X1*Z2*YZsum+2*a1*a3*X2*Y1*Z1*Z2+a4*XYsum*Z1*Z2+a4*XZsum*YZsum+(a1sq*a3+a1*a4)*X1*Z2*(XZsum+X2*Z1)+a2*a3*X2*Z1*(2*X1*Z2+X2*Z1)+a3sq*Y1*Z1*Z2**2+(a3sq+3*a6)*YZsum*Z1*Z2+a1*a3sq*(2*X1*Z2+X2*Z1)*Z1*Z2+3*a1*a6*X1*Z1*Z2**2+a3*a4*(XZsum+X2*Z1)*Z1*Z2+(a3**3+3*a3*a6)*Z1**2*Z2**2

    yield (X32, Y32, Z32)

