# sage.doctest: needs sage.symbolic
r"""
Further examples from Wester's paper

These are all the problems at
http://yacas.sourceforge.net/essaysmanual.html

They come from the 1994 paper "Review of CAS mathematical
capabilities", by Michael Wester, who put forward 123 problems that
a reasonable computer algebra system should be able to solve and
tested the then current versions of various commercial CAS on this
list. Sage can do most of the problems natively now, i.e., with no
explicit calls to Maxima or other systems.

::

    sage: # (YES) factorial of 50, and factor it
    sage: factorial(50)
    30414093201713378043612608166064768844377641568960512000000000000
    sage: factor(factorial(50))
    2^47 * 3^22 * 5^12 * 7^8 * 11^4 * 13^3 * 17^2 * 19^2 * 23^2 * 29 * 31 * 37 * 41 * 43 * 47

::

    sage: # (YES) 1/2+...+1/10 = 4861/2520
    sage: sum(1/n for n in range(2,10+1)) == 4861/2520
    True

::

    sage: # (YES) Evaluate  e^(Pi*Sqrt(163)) to 50 decimal digits
    sage: a = e^(pi*sqrt(163)); a
    e^(sqrt(163)*pi)
    sage: RealField(150)(a)
    2.6253741264076874399999999999925007259719820e17

::

    sage: # (YES) Evaluate the Bessel function J[2] numerically at z=1+I.
    sage: bessel_J(2, 1+I).n()
    0.0415798869439621 + 0.247397641513306*I

::

    sage: # (YES) Obtain period of decimal fraction 1/7=0.(142857).
    sage: a = 1/7
    sage: a
    1/7
    sage: a.period()
    6

::

    sage: # (YES) Continued fraction of 3.1415926535
    sage: a = 3.1415926535
    sage: continued_fraction(a)
    [3; 7, 15, 1, 292, 1, 1, 6, 2, 13, 4]

::

    sage: # (YES) Sqrt(2*Sqrt(3)+4)=1+Sqrt(3).
    sage: # The Maxima backend equality checker does this;
    sage: # note the equality only holds for one choice of sign,
    sage: # but Maxima always chooses the "positive" one
    sage: a = sqrt(2*sqrt(3) + 4); b = 1 + sqrt(3)
    sage: float(a-b)
    0.0
    sage: bool(a == b)
    True
    sage: # We can, of course, do this in a quadratic field
    sage: k.<sqrt3> = QuadraticField(3)
    sage: asqr = 2*sqrt3 + 4
    sage: b = 1+sqrt3
    sage: asqr == b^2
    True

::

    sage: # (YES) Sqrt(14+3*Sqrt(3+2*Sqrt(5-12*Sqrt(3-2*Sqrt(2)))))=3+Sqrt(2).
    sage: a = sqrt(14+3*sqrt(3+2*sqrt(5-12*sqrt(3-2*sqrt(2)))))
    sage: b = 3+sqrt(2)
    sage: a, b
    (sqrt(3*sqrt(2*sqrt(-12*sqrt(-2*sqrt(2) + 3) + 5) + 3) + 14), sqrt(2) + 3)
    sage: bool(a==b)
    True
    sage: abs(float(a-b)) < 1e-10
    True
    sage: # 2*Infinity-3=Infinity.
    sage: 2*infinity-3 == infinity
    True

::

    sage: # (YES) Standard deviation of the sample (1, 2, 3, 4, 5).
    sage: v = vector(RDF, 5, [1,2,3,4,5])
    sage: v.standard_deviation()
    1.5811388300841898

::

    sage: # (NO) Hypothesis testing with t-distribution.
    sage: # (NO) Hypothesis testing with chi^2 distribution
    sage: # (But both are included in Scipy and R)

::

    sage: # (YES) (x^2-4)/(x^2+4*x+4)=(x-2)/(x+2).
    sage: R.<x> = QQ[]
    sage: (x^2-4)/(x^2+4*x+4) == (x-2)/(x+2)
    True
    sage: restore('x')

::

    sage: # (YES -- Maxima doesn't immediately consider them
    sage: # equal, but simplification shows that they are)
    sage: # (Exp(x)-1)/(Exp(x/2)+1)=Exp(x/2)-1.
    sage: f = (exp(x)-1)/(exp(x/2)+1)
    sage: g = exp(x/2)-1
    sage: f
    (e^x - 1)/(e^(1/2*x) + 1)
    sage: g
    e^(1/2*x) - 1
    sage: f.canonicalize_radical()
    e^(1/2*x) - 1
    sage: g
    e^(1/2*x) - 1
    sage: f(x=10.0).n(53), g(x=10.0).n(53)
    (147.413159102577, 147.413159102577)
    sage: bool(f == g)
    True

::

    sage: # (YES) Expand (1+x)^20, take derivative and factorize.
    sage: # first do it using algebraic polys
    sage: R.<x> = QQ[]
    sage: f = (1+x)^20; f
    x^20 + 20*x^19 + 190*x^18 + 1140*x^17 + 4845*x^16 + 15504*x^15 + 38760*x^14 + 77520*x^13 + 125970*x^12 + 167960*x^11 + 184756*x^10 + 167960*x^9 + 125970*x^8 + 77520*x^7 + 38760*x^6 + 15504*x^5 + 4845*x^4 + 1140*x^3 + 190*x^2 + 20*x + 1
    sage: deriv = f.derivative()
    sage: deriv
    20*x^19 + 380*x^18 + 3420*x^17 + 19380*x^16 + 77520*x^15 + 232560*x^14 + 542640*x^13 + 1007760*x^12 + 1511640*x^11 + 1847560*x^10 + 1847560*x^9 + 1511640*x^8 + 1007760*x^7 + 542640*x^6 + 232560*x^5 + 77520*x^4 + 19380*x^3 + 3420*x^2 + 380*x + 20
    sage: deriv.factor()
    (20) * (x + 1)^19
    sage: restore('x')
    sage: # next do it symbolically
    sage: var('y')
    y
    sage: f = (1+y)^20; f
    (y + 1)^20
    sage: g = f.expand(); g
    y^20 + 20*y^19 + 190*y^18 + 1140*y^17 + 4845*y^16 + 15504*y^15 + 38760*y^14 + 77520*y^13 + 125970*y^12 + 167960*y^11 + 184756*y^10 + 167960*y^9 + 125970*y^8 + 77520*y^7 + 38760*y^6 + 15504*y^5 + 4845*y^4 + 1140*y^3 + 190*y^2 + 20*y + 1
    sage: deriv = g.derivative(); deriv
    20*y^19 + 380*y^18 + 3420*y^17 + 19380*y^16 + 77520*y^15 + 232560*y^14 + 542640*y^13 + 1007760*y^12 + 1511640*y^11 + 1847560*y^10 + 1847560*y^9 + 1511640*y^8 + 1007760*y^7 + 542640*y^6 + 232560*y^5 + 77520*y^4 + 19380*y^3 + 3420*y^2 + 380*y + 20
    sage: deriv.factor()
    20*(y + 1)^19

::

    sage: # (YES) Factorize x^100-1.
    sage: factor(x^100-1)
    (x^40 - x^30 + x^20 - x^10 + 1)*(x^20 + x^15 + x^10 + x^5 + 1)*(x^20 - x^15 + x^10 - x^5 + 1)*(x^8 - x^6 + x^4 - x^2 + 1)*(x^4 + x^3 + x^2 + x + 1)*(x^4 - x^3 + x^2 - x + 1)*(x^2 + 1)*(x + 1)*(x - 1)
    sage: # Also, algebraically
    sage: x = polygen(QQ)
    sage: factor(x^100 - 1)
    (x - 1) * (x + 1) * (x^2 + 1) * (x^4 - x^3 + x^2 - x + 1) * (x^4 + x^3 + x^2 + x + 1) * (x^8 - x^6 + x^4 - x^2 + 1) * (x^20 - x^15 + x^10 - x^5 + 1) * (x^20 + x^15 + x^10 + x^5 + 1) * (x^40 - x^30 + x^20 - x^10 + 1)
    sage: restore('x')

::

    sage: # (YES) Factorize  x^4-3*x^2+1 in the field of rational numbers extended by roots of  x^2-x-1.
    sage: x = polygen(ZZ, 'x')
    sage: k.< a> = NumberField(x^2 - x -1)
    sage: R.< y> = k[]
    sage: f = y^4 - 3*y^2 + 1
    sage: f
    y^4 - 3*y^2 + 1
    sage: factor(f)
    (y - a) * (y - a + 1) * (y + a - 1) * (y + a)

::

    sage: # (YES) Factorize  x^4-3*x^2+1 mod 5.
    sage: k.< x > = GF(5) [ ]
    sage: f = x^4 - 3*x^2 + 1
    sage: f.factor()
    (x + 2)^2 * (x + 3)^2
    sage: # Alternatively, from symbol x as follows:
    sage: reset('x')
    sage: f = x^4 - 3*x^2 + 1
    sage: f.polynomial(GF(5)).factor()
    (x + 2)^2 * (x + 3)^2

::

    sage: # (YES) Partial fraction decomposition of (x^2+2*x+3)/(x^3+4*x^2+5*x+2)
    sage: f = (x^2+2*x+3)/(x^3+4*x^2+5*x+2); f
    (x^2 + 2*x + 3)/(x^3 + 4*x^2 + 5*x + 2)
    sage: f.partial_fraction()
    3/(x + 2) - 2/(x + 1) + 2/(x + 1)^2

::

    sage: # (YES) Assuming  x>=y,  y>=z,  z>=x, deduce  x=z.
    sage: forget()
    sage: var('x,y,z')
    (x, y, z)
    sage: assume(x>=y, y>=z,z>=x)
    sage: bool(x==z)
    True

::

    sage: # (YES) Assuming x>y, y>0, deduce 2*x^2>2*y^2.
    sage: forget()
    sage: assume(x>y, y>0)
    sage: sorted(assumptions())
    [x > y, y > 0]
    sage: bool(2*x^2 > 2*y^2)
    True
    sage: forget()
    sage: assumptions()
    []

::

    sage: # (NO) Solve the inequality Abs(x-1)>2.
    sage: # Maxima doesn't solve inequalities
    sage: # (but some Maxima packages do):
    sage: eqn = abs(x-1) > 2
    sage: eqn.solve(x)
    [[x < -1], [3 < x]]

::

    sage: # (NO) Solve the inequality (x-1)*...*(x-5)<0.
    sage: eqn = prod(x-i for i in range(1,5 +1)) < 0
    sage: # but don't know how to solve
    sage: eqn.solve(x)
    [[x < 1], [x > 2, x < 3], [x > 4, x < 5]]

::

    sage: # (YES) Cos(3*x)/Cos(x)=Cos(x)^2-3*Sin(x)^2 or similar equivalent combination.
    sage: f = cos(3*x)/cos(x)
    sage: g = cos(x)^2 - 3*sin(x)^2
    sage: h = f-g
    sage: h.trig_simplify()
                                           0

::

    sage: # (YES) Cos(3*x)/Cos(x)=2*Cos(2*x)-1.
    sage: f = cos(3*x)/cos(x)
    sage: g = 2*cos(2*x) - 1
    sage: h = f-g
    sage: h.trig_simplify()
                                           0

::

    sage: # (GOOD ENOUGH) Define rewrite rules to match  Cos(3*x)/Cos(x)=Cos(x)^2-3*Sin(x)^2.
    sage: # Sage has no notion of "rewrite rules", but
    sage: # it can simplify both to the same thing.
    sage: (cos(3*x)/cos(x)).simplify_full()
    4*cos(x)^2 - 3
    sage: (cos(x)^2-3*sin(x)^2).simplify_full()
    4*cos(x)^2 - 3

::

    sage: # (YES) Sqrt(997)-(997^3)^(1/6)=0
    sage: a = sqrt(997) - (997^3)^(1/6)
    sage: a.simplify()
    0
    sage: bool(a == 0)
    True

::

    sage: # (YES) Sqrt(99983)-99983^3^(1/6)=0
    sage: a = sqrt(99983) - (99983^3)^(1/6)
    sage: bool(a==0)
    True
    sage: float(a)
    1.1368683772...e-13
    sage: 13*7691
    99983

::

    sage: # (YES) (2^(1/3) + 4^(1/3))^3 - 6*(2^(1/3) + 4^(1/3))-6 = 0
    sage: a = (2^(1/3) + 4^(1/3))^3 - 6*(2^(1/3) + 4^(1/3)) - 6; a
    (4^(1/3) + 2^(1/3))^3 - 6*4^(1/3) - 6*2^(1/3) - 6
    sage: bool(a==0)
    True
    sage: abs(float(a)) < 1e-10
    True

Or we can do it using number fields. ::

    sage: reset('x')
    sage: k.<b> = NumberField(x^3-2)
    sage: a = (b  + b^2)^3 - 6*(b  + b^2) - 6
    sage: a
    0

::

    sage: # (NO, except numerically) Ln(Tan(x/2+Pi/4))-ArcSinh(Tan(x))=0
    # Sage uses the Maxima convention when comparing symbolic expressions and
    # returns True only when it can prove equality. Thus, in this case, we get
    # False even though the equality holds.
    sage: f = log(tan(x/2 + pi/4)) - arcsinh(tan(x))
    sage: bool(f == 0)
    False
    sage: [abs(float(f(x=i/10))) < 1e-15 for i in range(1,5)]
    [True, True, True, True]
    sage: # Numerically, the expression Ln(Tan(x/2+Pi/4))-ArcSinh(Tan(x))=0 and its derivative at x=0 are zero.
    sage: g = f.derivative()
    sage: abs(float(f(x=0))) < 1e-10
    True
    sage: abs(float(g(x=0))) < 1e-10
    True
    sage: g
    -sqrt(tan(x)^2 + 1) + 1/2*(tan(1/4*pi + 1/2*x)^2 + 1)/tan(1/4*pi + 1/2*x)

::

    sage: # (NO) Ln((2*Sqrt(r) + 1)/Sqrt(4*r 4*Sqrt(r) 1))=0.
    sage: var('r')
    r
    sage: f = log( (2*sqrt(r) + 1) / sqrt(4*r  + 4*sqrt(r) +  1)); f
    log((2*sqrt(r) + 1)/sqrt(4*r + 4*sqrt(r) + 1))
    sage: bool(f == 0)
    False
    sage: [abs(float(f(r=i))) < 1e-10 for i in [0.1,0.3,0.5]]
    [True, True, True]

::

    sage: # (NO)
    sage: # (4*r+4*Sqrt(r)+1)^(Sqrt(r)/(2*Sqrt(r)+1))*(2*Sqrt(r)+1)^(2*Sqrt(r)+1)^(-1)-2*Sqrt(r)-1=0, assuming r>0.
    sage: assume(r>0)
    sage: f = (4*r+4*sqrt(r)+1)^(sqrt(r)/(2*sqrt(r)+1))*(2*sqrt(r)+1)^(2*sqrt(r)+1)^(-1)-2*sqrt(r)-1
    sage: f
    (4*r + 4*sqrt(r) + 1)^(sqrt(r)/(2*sqrt(r) + 1))*(2*sqrt(r) + 1)^(1/(2*sqrt(r) + 1)) - 2*sqrt(r) - 1
    sage: bool(f == 0)
    False
    sage: [abs(float(f(r=i))) < 1e-10 for i in [0.1,0.3,0.5]]
    [True, True, True]

::

    sage: # (YES) Obtain real and imaginary parts of Ln(3+4*I).
    sage: a = log(3+4*I); a
    log(4*I + 3)
    sage: a.real()
    log(5)
    sage: a.imag()
    arctan(4/3)

::

    sage: # (YES) Obtain real and imaginary parts of Tan(x+I*y)
    sage: z = var('z')
    sage: a = tan(z); a
    tan(z)
    sage: a.real()
    sin(2*real_part(z))/(cos(2*real_part(z)) + cosh(2*imag_part(z)))
    sage: a.imag()
    sinh(2*imag_part(z))/(cos(2*real_part(z)) + cosh(2*imag_part(z)))


::

    sage: # (YES) Simplify Ln(Exp(z)) to z for -Pi<Im(z)<=Pi.
    sage: # Unfortunately (?), Maxima does this even without
    sage: # any assumptions.
    sage: # We *would* use assume(-pi < imag(z))
    sage: # and assume(imag(z) <= pi)
    sage: f = log(exp(z)); f
    log(e^z)
    sage: f.simplify()
    z
    sage: forget()

::

    sage: # (YES) Assuming Re(x)>0, Re(y)>0, deduce x^(1/n)*y^(1/n)-(x*y)^(1/n)=0.
    sage: # Maxima 5.26 has different behaviours depending on the current
    sage: # domain.
    sage: # To stick with the behaviour of previous versions, the domain is set
    sage: # to 'real' in the following.
    sage: # See Issue #10682 for further details.
    sage: n = var('n')
    sage: f = x^(1/n)*y^(1/n)-(x*y)^(1/n)
    sage: assume(real(x) > 0, real(y) > 0)
    sage: f.simplify()
    x^(1/n)*y^(1/n) - (x*y)^(1/n)
    sage: maxima = sage.calculus.calculus.maxima
    sage: maxima.set('domain', 'real') # set domain to real
    sage: f.simplify()
    0
    sage: maxima.set('domain', 'complex') # set domain back to its default value
    sage: forget()

::

    sage: # (YES) Transform equations, (x==2)/2+(1==1)=>x/2+1==2.
    sage: eq1 = x == 2
    sage: eq2 = SR(1) == SR(1)
    sage: eq1/2 + eq2
    1/2*x + 1 == 2

::

    sage: # (SOMEWHAT) Solve Exp(x)=1 and get all solutions.
    sage: # to_poly_solve in Maxima can do this.
    sage: solve(exp(x) == 1, x)
    [x == 0]

::

    sage: # (SOMEWHAT) Solve Tan(x)=1 and get all solutions.
    sage: # to_poly_solve in Maxima can do this.
    sage: solve(tan(x) == 1, x)
    [x == 1/4*pi]

::

    sage: # (YES) Solve a degenerate 3x3 linear system.
    sage: # x+y+z==6,2*x+y+2*z==10,x+3*y+z==10
    sage: # First symbolically:
    sage: solve([x+y+z==6, 2*x+y+2*z==10, x+3*y+z==10], x,y,z)
    [[x == -r1 + 4, y == 2, z == r1]]

::

    sage: # (YES) Invert a 2x2 symbolic matrix.
    sage: # [[a,b],[1,a*b]]
    sage: # Using multivariate poly ring -- much nicer
    sage: R.<a,b> = QQ[]
    sage: m = matrix(2,2,[a,b,  1, a*b])
    sage: zz = m^(-1)
    sage: zz
    [     a/(a^2 - 1)   (-1)/(a^2 - 1)]
    [(-1)/(a^2*b - b)    a/(a^2*b - b)]

::

    sage: # (YES) Compute and factor the determinant of the 4x4 Vandermonde matrix in a, b, c, d.
    sage: var('a,b,c,d')
    (a, b, c, d)
    sage: m = matrix(SR, 4, 4, [[z^i for i in range(4)] for z in [a,b,c,d]])
    sage: m
    [  1   a a^2 a^3]
    [  1   b b^2 b^3]
    [  1   c c^2 c^3]
    [  1   d d^2 d^3]
    sage: d = m.determinant()
    sage: d.factor()
    (a - b)*(a - c)*(a - d)*(b - c)*(b - d)*(c - d)

::

    sage: # (YES) Compute and factor the determinant of the 4x4 Vandermonde matrix in a, b, c, d.
    sage: # Do it instead in a multivariate ring
    sage: R.<a,b,c,d> = QQ[]
    sage: m = matrix(R, 4, 4, [[z^i for i in range(4)] for z in [a,b,c,d]])
    sage: m
    [  1   a a^2 a^3]
    [  1   b b^2 b^3]
    [  1   c c^2 c^3]
    [  1   d d^2 d^3]
    sage: d = m.determinant()
    sage: d
    a^3*b^2*c - a^2*b^3*c - a^3*b*c^2 + a*b^3*c^2 + a^2*b*c^3 - a*b^2*c^3 - a^3*b^2*d + a^2*b^3*d + a^3*c^2*d - b^3*c^2*d - a^2*c^3*d + b^2*c^3*d + a^3*b*d^2 - a*b^3*d^2 - a^3*c*d^2 + b^3*c*d^2 + a*c^3*d^2 - b*c^3*d^2 - a^2*b*d^3 + a*b^2*d^3 + a^2*c*d^3 - b^2*c*d^3 - a*c^2*d^3 + b*c^2*d^3
    sage: d.factor()
    (-1) * (c - d) * (-b + c) * (b - d) * (-a + c) * (-a + b) * (a - d)

::

    sage: # (YES) Find the eigenvalues of a 3x3 integer matrix.
    sage: m = matrix(QQ, 3, [5,-3,-7, -2,1,2, 2,-3,-4])
    sage: m.eigenspaces_left()
    [(3,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [ 1  0 -1]),
     (1,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [ 1  1 -1]),
     (-2,
      Vector space of degree 3 and dimension 1 over Rational Field
      User basis matrix:
      [0 1 1])]

::

    sage: # (YES) Verify some standard limits found by L'Hopital's rule:
    sage: #   Verify(Limit(x,Infinity) (1+1/x)^x, Exp(1));
    sage: #   Verify(Limit(x,0) (1-Cos(x))/x^2, 1/2);
    sage: limit( (1+1/x)^x, x = oo)
    e
    sage: limit( (1-cos(x))/(x^2), x = 1/2)
    -4*cos(1/2) + 4

::

    sage: # (OK-ish) D(x)Abs(x)
    sage: #    Verify(D(x) Abs(x), Sign(x));
    sage: diff(abs(x))
    1/2*(x + conjugate(x))/abs(x)
    sage: _.simplify_full()
    x/abs(x)
    sage: _ = var('x', domain='real')
    sage: diff(abs(x))
    x/abs(x)
    sage: forget()

::

    sage: # (YES) (Integrate(x)Abs(x))=Abs(x)*x/2
    sage: integral(abs(x), x)
    1/2*x*abs(x)

::

    sage: #  (YES) Compute derivative of Abs(x), piecewise defined.
    sage: #     Verify(D(x)if(x<0) (-x) else x,
    sage: #        Simplify(if(x<0) -1 else 1))
    Piecewise defined function with 2 parts, [[(-10, 0), -1], [(0, 10), 1]]
    sage: #  (NOT really) Integrate Abs(x), piecewise defined.
    sage: #      Verify(Simplify(Integrate(x)
    sage: #        if(x<0) (-x) else x),
    sage: #        Simplify(if(x<0) (-x^2/2) else x^2/2));
    sage: f = piecewise([ ((-10,0), -x), ((0,10), x)])
    sage: f.integral(definite=True)
    100

::

    sage: # (YES) Taylor series of 1/Sqrt(1-v^2/c^2) at v=0.
    sage: var('v,c')
    (v, c)
    sage: taylor(1/sqrt(1-v^2/c^2), v, 0, 7)
    1/2*v^2/c^2 + 3/8*v^4/c^4 + 5/16*v^6/c^6 + 1

::

    sage: # (OK-ish) (Taylor expansion of Sin(x))/(Taylor expansion of Cos(x)) = (Taylor expansion of Tan(x)).
    sage: #      TestYacas(Taylor(x,0,5)(Taylor(x,0,5)Sin(x))/
    sage: #        (Taylor(x,0,5)Cos(x)), Taylor(x,0,5)Tan(x));
    sage: f = taylor(sin(x), x, 0, 8)
    sage: g = taylor(cos(x), x, 0, 8)
    sage: h = taylor(tan(x), x, 0, 8)
    sage: f = f.power_series(QQ)
    sage: g = g.power_series(QQ)
    sage: h = h.power_series(QQ)
    sage: f - g*h
    O(x^8)

::

    sage: # (YES) Taylor expansion of Ln(x)^a*Exp(-b*x) at x=1.
    sage: a,b = var('a,b')
    sage: taylor(log(x)^a*exp(-b*x), x, 1, 3)
    -1/48*(a^3*(x - 1)^a + a^2*(6*b + 5)*(x - 1)^a + 8*b^3*(x - 1)^a + 2*(6*b^2 + 5*b + 3)*a*(x - 1)^a)*(x - 1)^3*e^(-b) + 1/24*(3*a^2*(x - 1)^a + a*(12*b + 5)*(x - 1)^a + 12*b^2*(x - 1)^a)*(x - 1)^2*e^(-b) - 1/2*(a*(x - 1)^a + 2*b*(x - 1)^a)*(x - 1)*e^(-b) + (x - 1)^a*e^(-b)

::

    sage: # (YES) Taylor expansion of Ln(Sin(x)/x) at x=0.
    sage: taylor(log(sin(x)/x), x, 0, 10)
    -1/467775*x^10 - 1/37800*x^8 - 1/2835*x^6 - 1/180*x^4 - 1/6*x^2

::

    sage: # (NO) Compute n-th term of the Taylor series of Ln(Sin(x)/x) at x=0.
    sage: # need formal functions

::

    sage: # (NO) Compute n-th term of the Taylor series of Exp(-x)*Sin(x) at x=0.
    sage: # (Sort of, with some work)
    sage: # Solve x=Sin(y)+Cos(y) for y as Taylor series in x at x=1.
    sage: #      TestYacas(InverseTaylor(y,0,4) Sin(y)+Cos(y),
    sage: #        (y-1)+(y-1)^2/2+2*(y-1)^3/3+(y-1)^4);
    sage: #       Note that InverseTaylor does not give the series in terms of x but in terms of y which is semantically
    sage: # wrong. But other CAS do the same.
    sage: f = sin(y) + cos(y)
    sage: g = f.taylor(y, 0, 10)
    sage: h = g.power_series(QQ)
    sage: k = (h - 1).reverse()
    sage: k
    y + 1/2*y^2 + 2/3*y^3 + y^4 + 17/10*y^5 + 37/12*y^6 + 41/7*y^7 + 23/2*y^8 + 1667/72*y^9 + 3803/80*y^10 + O(y^11)

::

    sage: # (OK) Compute Legendre polynomials directly from Rodrigues's formula, P[n]=1/(2^n*n!) *(Deriv(x,n)(x^2-1)^n).
    sage: #      P(n,x) := Simplify( 1/(2*n)!! *
    sage: #        Deriv(x,n) (x^2-1)^n );
    sage: #      TestYacas(P(4,x), (35*x^4)/8+(-15*x^2)/4+3/8);
    sage: P = lambda n, x: simplify(diff((x^2-1)^n,x,n) / (2^n * factorial(n)))
    sage: P(4,x).expand()
    35/8*x^4 - 15/4*x^2 + 3/8

::

    sage: # (YES) Define the polynomial p=Sum(i,1,5,a[i]*x^i).
    sage: # symbolically
    sage: ps = sum(var('a%s'%i)*x^i for i in range(1,6)); ps
    a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x
    sage: ps.parent()
    Symbolic Ring
    sage: # algebraically
    sage: R = PolynomialRing(QQ,5,names='a')
    sage: S.<x> = PolynomialRing(R)
    sage: p = S(list(R.gens()))*x; p
    a4*x^5 + a3*x^4 + a2*x^3 + a1*x^2 + a0*x
    sage: p.parent()
    Univariate Polynomial Ring in x over Multivariate Polynomial Ring in a0, a1, a2, a3, a4 over Rational Field

::

    sage: # (YES) Convert the above to Horner's form.
    sage: #      Verify(Horner(p, x), ((((a[5]*x+a[4])*x
    sage: #        +a[3])*x+a[2])*x+a[1])*x);
    sage: restore('x')
    sage: SR(p).horner(x)
    ((((a4*x + a3)*x + a2)*x + a1)*x + a0)*x

::

    sage: # (NO) Convert the result of problem 127 to Fortran syntax.
    sage: #      CForm(Horner(p, x));

::

    sage: # (YES) Verify that True And False=False.
    sage: (True and False) is False
    True

::

    sage: # (YES) Prove x Or Not x.
    sage: for x in [True, False]:
    ....:    print(x or (not x))
    True
    True

::

    sage: # (YES) Prove x Or y Or x And y=>x Or y.
    sage: for x in [True, False]:
    ....:   for y in [True, False]:
    ....:       if x or y or x and y:
    ....:           if not (x or y):
    ....:              print("failed!")
"""
