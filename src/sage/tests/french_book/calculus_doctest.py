## -*- encoding: utf-8 -*-
"""
Doctests from French Sage book
Test file for chapter "Analyse et algèbre avec Sage" ("Calculus and
algebra with Sage")

Tests extracted from ./calculus.tex.

Sage example in ./calculus.tex, line 37::

    sage: bool(x^2 + 3*x + 1 == (x+1)*(x+2))
    False

Sage example in ./calculus.tex, line 74::

    sage: a, x = var('a, x'); y = cos(x+a) * (x+1); y
    (x + 1)*cos(a + x)
    sage: y.subs(a=-x); y.subs(x=pi/2, a=pi/3); y.subs(x=0.5, a=2.3)
    x + 1
    -1/4*sqrt(3)*(pi + 2)
    -1.41333351100299
    sage: y(a=-x); y(x=pi/2, a=pi/3); y(x=0.5, a=2.3)
    x + 1
    -1/4*sqrt(3)*(pi + 2)
    -1.41333351100299

Sage example in ./calculus.tex, line 91::

    sage: x, y, z = var('x, y, z') ; q = x*y + y*z + z*x
    sage: bool(q(x=y, y=z, z=x) == q), bool(q(z=y)(y=x) == 3*x^2)
    (True, True)

Sage example in ./calculus.tex, line 99::

    sage: y, z = var('y, z'); f = x^3 + y^2 + z
    sage: f.subs(x^3 == y^2, z==1)
    2*y^2 + 1

Sage example in ./calculus.tex, line 110::

    sage: f(x)=(2*x+1)^3 ; f(-3)
    -125
    sage: f.expand()
    x |--> 8*x^3 + 12*x^2 + 6*x + 1

Sage example in ./calculus.tex, line 122::

    sage: y = var('y'); u = sin(x) + x*cos(y)
    sage: v = u.function(x,y); v
    (x, y) |--> x*cos(y) + sin(x)
    sage: w(x, y) = u; w
    (x, y) |--> x*cos(y) + sin(x)

Sage example in ./calculus.tex, line 153::

    sage: x, y = SR.var('x,y')
    sage: p = (x+y)*(x+1)^2
    sage: p2 = p.expand(); p2
    x^3 + x^2*y + 2*x^2 + 2*x*y + x + y

Sage example in ./calculus.tex, line 160::

    sage: p2.collect(x)
    x^3 + x^2*(y + 2) + x*(2*y + 1) + y

Sage example in ./calculus.tex, line 165::

    sage: ((x+y+sin(x))^2).expand().collect(sin(x))
    x^2 + 2*x*y + y^2 + 2*(x + y)*sin(x) + sin(x)^2

Sage example in ./calculus.tex, line 254::

    sage: (x^x/x).simplify()
    x^(x - 1)

Sage example in ./calculus.tex, line 260::

    sage: f = (e^x-1) / (1+e^(x/2)); f.canonicalize_radical()
    e^(1/2*x) - 1

Sage example in ./calculus.tex, line 266::

    sage: f = cos(x)^6 + sin(x)^6 + 3 * sin(x)^2 * cos(x)^2
    sage: f.simplify_trig()
    1

Sage example in ./calculus.tex, line 273::

    sage: f = cos(x)^6; f.reduce_trig()
    1/32*cos(6*x) + 3/16*cos(4*x) + 15/32*cos(2*x) + 5/16
    sage: f = sin(5 * x); f.expand_trig()
    5*cos(x)^4*sin(x) - 10*cos(x)^2*sin(x)^3 + sin(x)^5

Sage example in ./calculus.tex, line 306::

    sage: n = var('n'); f = factorial(n+1)/factorial(n)
    sage: f.simplify_factorial()
    n + 1

Sage example in ./calculus.tex, line 318::

    sage: f = sqrt(abs(x)^2); f.canonicalize_radical()
    abs(x)
    sage: f = log(x*y); f.canonicalize_radical()
    log(x) + log(y)

Sage example in ./calculus.tex, line 371::

    sage: assume(x > 0); bool(sqrt(x^2) == x)
    True
    sage: forget(x > 0); bool(sqrt(x^2) == x)
    False
    sage: n = var('n'); assume(n, 'integer'); sin(n*pi).simplify()
    0

Sage example in ./calculus.tex, line 420::

    sage: a = var('a')
    sage: c = (a+1)^2 - (a^2+2*a+1)

Sage example in ./calculus.tex, line 425::

    sage: eq =  c * x == 0

Sage example in ./calculus.tex, line 430::

    sage: eq2 = eq / c; eq2
    x == 0
    sage: solve(eq2, x)
    [x == 0]

Sage example in ./calculus.tex, line 437::

    sage: solve(eq, x)
    [x == x]

Sage example in ./calculus.tex, line 444::

    sage: expand(c)
    0

Sage example in ./calculus.tex, line 452::

    sage: c = cos(a)^2 + sin(a)^2 - 1
    sage: eq = c*x == 0
    sage: solve(eq, x)
    [x == 0]

Sage example in ./calculus.tex, line 460::

    sage: c.simplify_trig()
    0
    sage: c.is_zero()
    True

Sage example in ./calculus.tex, line 516::

    sage: z, phi = var('z, phi')
    sage: eq =  z**2 - 2/cos(phi)*z + 5/cos(phi)**2 - 4 == 0; eq
    z^2 - 2*z/cos(phi) + 5/cos(phi)^2 - 4 == 0

Sage example in ./calculus.tex, line 523::

    sage: eq.lhs()
    z^2 - 2*z/cos(phi) + 5/cos(phi)^2 - 4
    sage: eq.rhs()
    0

Sage example in ./calculus.tex, line 531::

    sage: solve(eq, z)
    [z == -(2*sqrt(cos(phi)^2 - 1) - 1)/cos(phi),
    z ==  (2*sqrt(cos(phi)^2 - 1) + 1)/cos(phi)]

Sage example in ./calculus.tex, line 537::

    sage: y = var('y'); solve(y^6==y, y)
    [y == 1/4*sqrt(5) + 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4, y == -1/4*sqrt(5) + 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4, y == -1/4*sqrt(5) - 1/4*I*sqrt(-2*sqrt(5) + 10) - 1/4, y == 1/4*sqrt(5) - 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4, y == 1, y == 0]

Sage example in ./calculus.tex, line 544::

    sage: solve(x^2-1, x, solution_dict=True)
    [{x: -1}, {x: 1}]

Sage example in ./calculus.tex, line 550::

    sage: solve([x+y == 3, 2*x+2*y == 6], x, y)
    [[x == -r1 + 3, y == r1]]

Sage example in ./calculus.tex, line 560::

    sage: solve([cos(x)*sin(x) == 1/2, x+y == 0], x, y)
    [[x == 1/4*pi + pi*z..., y == -1/4*pi - pi*z...]]

Sage example in ./calculus.tex, line 565::

    sage: solve(x^2+x-1 > 0, x)
    [[x < -1/2*sqrt(5) - 1/2], [x > 1/2*sqrt(5) - 1/2]]

Sage example in ./calculus.tex, line 583::

    sage: x, y, z = var('x, y, z')
    sage: solve([x^2 * y * z == 18, x * y^3 * z == 24,\
    ....:        x * y * z^4 == 3], x, y, z)
    [[x == (-2.767364733... - 1.713479699...*I), y == (-0.5701035039... + 2.003705978...*I), z == (-0.8016843376... - 0.1498607749...*I)], ...]

Sage example in ./calculus.tex, line 597::

    sage: expr = sin(x) + sin(2 * x) + sin(3 * x)
    sage: solve(expr, x)
    [sin(3*x) == -sin(2*x) - sin(x)]

Sage example in ./calculus.tex, line 605::

    sage: find_root(expr, 0.1, pi)
    2.094395102393195...

Sage example in ./calculus.tex, line 610::

    sage: f = expr.simplify_trig(); f
    2*(2*cos(x)^2 + cos(x))*sin(x)
    sage: solve(f, x)
    [x == 0, x == 2/3*pi, x == 1/2*pi]

Sage example in ./calculus.tex, line 629::

    sage: (x^3+2*x+1).roots(x)
    [(-1/2*(1/18*sqrt(59)*sqrt(3) - 1/2)^(1/3)*(I*sqrt(3) + 1)
    - 1/3*(I*sqrt(3) - 1)/(1/18*sqrt(59)*sqrt(3) - 1/2)^(1/3), 1),
      (-1/2*(1/18*sqrt(59)*sqrt(3) - 1/2)^(1/3)*(-I*sqrt(3) + 1)
    - 1/3*(-I*sqrt(3) - 1)/(1/18*sqrt(59)*sqrt(3) - 1/2)^(1/3), 1),
      ((1/18*sqrt(59)*sqrt(3) - 1/2)^(1/3) - 2/3/(1/18*sqrt(59)*sqrt(3)
    - 1/2)^(1/3), 1)]

Sage example in ./calculus.tex, line 658::

    sage: (x^3+2*x+1).roots(x, ring=RR)
    [(-0.453397651516404, 1)]

Sage example in ./calculus.tex, line 662::

    sage: (x^3+2*x+1).roots(x, ring=CC)
    [(-0.453397651516404, 1),
    (0.226698825758202 - 1.46771150871022*I, 1),
    (0.226698825758202 + 1.46771150871022*I, 1)]

Sage example in ./calculus.tex, line 680::

    sage: solve(x^(1/x)==(1/x)^x, x)
    [(1/x)^x == x^(1/x)]

Sage example in ./calculus.tex, line 706::

    sage: y = function('y')(x)
    sage: desolve(diff(y,x,x) + x*diff(y,x) + y == 0, y, [0,0,1])
    -1/2*I*sqrt(2)*sqrt(pi)*erf(1/2*I*sqrt(2)*x)*e^(-1/2*x^2)

Sage example in ./calculus.tex, line 733::

    sage: k, n = var('k, n')
    sage: sum(k, k, 1, n).factor()
    1/2*(n + 1)*n

Sage example in ./calculus.tex, line 739::

    sage: n, k, y = var('n, k, y')
    sage: sum(binomial(n,k) * x^k * y^(n-k), k, 0, n)
    (x + y)^n

Sage example in ./calculus.tex, line 745::

    sage: k, n = var('k, n')
    sage: sum(binomial(n,k), k, 0, n),\
    ....: sum(k * binomial(n, k), k, 0, n),\
    ....: sum((-1)^k*binomial(n,k), k, 0, n)
    (2^n, 2^(n - 1)*n, 0)

Sage example in ./calculus.tex, line 753::

    sage: a, q, k, n = var('a, q, k, n')
    sage: sum(a*q^k, k, 0, n)
    (a*q^(n + 1) - a)/(q - 1)

Sage example in ./calculus.tex, line 760::

    sage: assume(abs(q) < 1)
    sage: sum(a*q^k, k, 0, infinity)
    -a/(q - 1)

Sage example in ./calculus.tex, line 766::

    sage: forget(); assume(q > 1); sum(a*q^k, k, 0, infinity)
    Traceback (most recent call last):
    ...
    ValueError: Sum is divergent.

Sage example in ./calculus.tex, line 842::

    sage: limit((x**(1/3) - 2) / ((x + 19)**(1/3) - 3), x = 8)
    9/4
    sage: f(x) = (cos(pi/4-x)-tan(x))/(1-sin(pi/4 + x))
    sage: limit(f(x), x = pi/4)
    Infinity

Sage example in ./calculus.tex, line 855::

    sage: limit(f(x), x = pi/4, dir='minus')
    +Infinity
    sage: limit(f(x), x = pi/4, dir='plus')
    -Infinity

Sage example in ./calculus.tex, line 898::

    sage: u(n) = n^100 / 100^n
    sage: u(2.);u(3.);u(4.);u(5.);u(6.);u(7.);u(8.);u(9.);u(10.)
    1.26765060022823e26
    5.15377520732011e41
    1.60693804425899e52
    7.88860905221012e59
    6.53318623500071e65
    3.23447650962476e70
    2.03703597633449e74
    2.65613988875875e77
    1.00000000000000e80

Sage example in ./calculus.tex, line 914::

    sage: plot(u(x), x, 1, 40)
    Graphics object consisting of 1 graphics primitive

Sage example in ./calculus.tex, line 929::

    sage: v(x) = diff(u(x), x); sol = solve(v(x) == 0, x); sol
    [x == 100/log(100), x == 0]
    sage: floor(sol[0].rhs())
    21

Sage example in ./calculus.tex, line 938::

    sage: limit(u(n), n=infinity)
    0
    sage: n0 = find_root(u(n) - 1e-8 == 0, 22, 1000); n0
    105.07496210187252

Sage example in ./calculus.tex, line 988::

    sage: taylor((1+arctan(x))**(1/x), x, 0, 3)
    1/16*x^3*e + 1/8*x^2*e - 1/2*x*e + e

Sage example in ./calculus.tex, line 993::

    sage: (ln(2*sin(x))).series(x==pi/6, 3)
    (sqrt(3))*(-1/6*pi + x) + (-2)*(-1/6*pi + x)^2
    + Order(-1/216*(pi - 6*x)^3)

Sage example in ./calculus.tex, line 1002::

    sage: (ln(2*sin(x))).series(x==pi/6, 3).truncate()
    -1/18*(pi - 6*x)^2 - 1/6*sqrt(3)*(pi - 6*x)

Sage example in ./calculus.tex, line 1017::

    sage: taylor((x**3+x)**(1/3) - (x**3-x)**(1/3), x, infinity, 2)
    2/3/x

Sage example in ./calculus.tex, line 1041::

    sage: tan(4*arctan(1/5)).simplify_trig()
    120/119
    sage: tan(pi/4+arctan(1/239)).simplify_trig()
    120/119

Sage example in ./calculus.tex, line 1052::

    sage: f = arctan(x).series(x, 10); f
    1*x + (-1/3)*x^3 + 1/5*x^5 + (-1/7)*x^7 + 1/9*x^9 + Order(x^10)
    sage: (16*f.subs(x==1/5) - 4*f.subs(x==1/239)).n(); pi.n()
    3.14159268240440
    3.14159265358979

Sage example in ./calculus.tex, line 1093::

    sage: k = var('k')
    sage: sum(1/k^2, k, 1, infinity),\
    ....: sum(1/k^4, k, 1, infinity),\
    ....: sum(1/k^5, k, 1, infinity)
    (1/6*pi^2, 1/90*pi^4, zeta(5))

Sage example in ./calculus.tex, line 1111::

    sage: s = 2*sqrt(2)/9801*(sum((factorial(4*k)) * (1103+26390*k) /
    ....:     ((factorial(k)) ^ 4 * 396 ^ (4 * k)) for k in (0..11)))
    sage: (1/s).n(digits=100)
    3.141592653589793238462643383279502884197169399375105820974...
    sage: (pi-1/s).n(digits=100).n()
    -4.36415445739398e-96

Sage example in ./calculus.tex, line 1139::

    sage: n = var('n'); u = sin(pi*(sqrt(4*n^2+1)-2*n))
    sage: taylor(u, n, infinity, 3)
    1/4*pi/n - 1/384*(6*pi + pi^3)/n^3

Sage example in ./calculus.tex, line 1163::

    sage: diff(sin(x^2), x)
    2*x*cos(x^2)
    sage: function('f')(x); function('g')(x); diff(f(g(x)), x)
    f(x)
    g(x)
    D[0](f)(g(x))*diff(g(x), x)
    sage: diff(ln(f(x)), x)
    diff(f(x), x)/f(x)

Sage example in ./calculus.tex, line 1180::

    sage: f(x,y) = x*y + sin(x^2) + e^(-x); derivative(f, x)
    (x, y) |--> 2*x*cos(x^2) + y - e^(-x)
    sage: derivative(f, y)
    (x, y) |--> x

Sage example in ./calculus.tex, line 1195::

    sage: x, y = var('x, y'); f = ln(x**2+y**2) / 2
    sage: delta = diff(f,x,2) + diff(f,y,2)
    sage: delta.simplify_full()
    0

Sage example in ./calculus.tex, line 1231::

    sage: sin(x).integral(x, 0, pi/2)
    1
    sage: integrate(1/(1+x^2), x)
    arctan(x)
    sage: integrate(1/(1+x^2), x, -infinity, infinity)
    pi
    sage: integrate(exp(-x**2), x, 0, infinity)
    1/2*sqrt(pi)

Sage example in ./calculus.tex, line 1241::

    sage: integrate(exp(-x), x, -infinity, infinity)
    Traceback (most recent call last):
    ...
    ValueError: Integral is divergent.

Sage example in ./calculus.tex, line 1254::

    sage: u = var('u'); f = x * cos(u) / (u^2 + x^2)
    sage: assume(x>0); f.integrate(u, 0, infinity)
    1/2*pi*e^(-x)
    sage: forget(); assume(x<0); f.integrate(u, 0, infinity)
    -1/2*pi*e^x

Sage example in ./calculus.tex, line 1270::

    sage: integral_numerical(sin(x)/x, 0, 1)   # abs tol 1e-12
    (0.94608307036718287, 1.0503632079297086e-14)
    sage: g = integrate(exp(-x**2), x, 0, infinity)
    sage: g, g.n()                             # abs tol 1e-12
    (1/2*sqrt(pi), 0.886226925452758)
    sage: approx = integral_numerical(exp(-x**2), 0, infinity)
    sage: approx                               # abs tol 1e-12
    (0.88622692545275705, 1.7147744320162414e-08)
    sage: approx[0]-g.n()                      # abs tol 1e-12
    -8.88178419700125e-16

Sage example in ./calculus.tex, line 1482::

    sage: A = matrix(QQ, [[1,2],[3,4]]); A
    [1 2]
    [3 4]

Sage example in ./calculus.tex, line 1629::

    sage: A = matrix(QQ, [[2,4,3],[-4,-6,-3],[3,3,1]])
    sage: A.characteristic_polynomial()
    x^3 + 3*x^2 - 4
    sage: A.eigenvalues()
    [1, -2, -2]
    sage: A.minimal_polynomial().factor()
    (x - 1) * (x + 2)^2

Sage example in ./calculus.tex, line 1641::

    sage: A.eigenvectors_right()
    [(1, [
    (1, -1, 1)
    ], 1), (-2, [
    (1, -1, 0)
    ], 2)]

Sage example in ./calculus.tex, line 1652::

    sage: A.jordan_form(transformation=True)
    (
    [ 1| 0  0]
    [--+-----]  [ 1  1  1]
    [ 0|-2  1]  [-1 -1  0]
    [ 0| 0 -2], [ 1  0 -1]
    )

Sage example in ./calculus.tex, line 1686::

    sage: A = matrix(QQ, [[1,-1/2],[-1/2,-1]])
    sage: A.jordan_form()
    Traceback (most recent call last):
    ...
    RuntimeError: Some eigenvalue does not exist in Rational Field.

Sage example in ./calculus.tex, line 1695::

    sage: A = matrix(QQ, [[1,-1/2],[-1/2,-1]])
    sage: A.minimal_polynomial()
    x^2 - 5/4

Sage example in ./calculus.tex, line 1701::

    sage: R = QQ[sqrt(5)]
    sage: A = A.change_ring(R)
    sage: A.jordan_form(transformation=True, subdivide=False)
    (
    [ 1/2*sqrt5          0]  [         1          1]
    [         0 -1/2*sqrt5], [-sqrt5 + 2  sqrt5 + 2]
    )

Sage example in ./calculus.tex, line 1734::

    sage: K.<sqrt2> = NumberField(x^2 - 2)
    sage: L.<sqrt3> = K.extension(x^2 - 3)
    sage: A = matrix(L, [[2, sqrt2*sqrt3, sqrt2],  \
    ....:                [sqrt2*sqrt3, 3, sqrt3],  \
    ....:                [sqrt2, sqrt3, 1]])
    sage: A.jordan_form(transformation=True)
    (
    [6|0|0]
    [-+-+-]
    [0|0|0]  [              1               1               0]
    [-+-+-]  [1/2*sqrt2*sqrt3               0               1]
    [0|0|0], [      1/2*sqrt2          -sqrt2          -sqrt3]
    )

"""

"""
Tests extracted from sol/calculus.tex.

Sage example in ./sol/calculus.tex, line 3::

    sage: reset()

Sage example in ./sol/calculus.tex, line 9::

    sage: n, k = var('n, k'); p = 4; s = [n + 1]
    sage: for k in (1..p):
    ....:  s += [factor((((n+1)^(k+1) \
    ....:       - sum(binomial(k+1, j)\
    ....:       * s[j] for j in (0..k-1))) / (k+1)))]
    ...
    sage: s
    [n + 1, 1/2*(n + 1)*n, 1/6*(2*n + 1)*(n + 1)*n,
    1/4*(n + 1)^2*n^2, 1/30*(3*n^2 + 3*n - 1)*(2*n + 1)*(n + 1)*n]

Sage example in ./sol/calculus.tex, line 34::

    sage: x, h, a = var('x, h, a'); f = function('f')
    sage: g(x) = taylor(f(x), x, a, 3)
    sage: phi(h) = (g(a+3*h) - 3*g(a+2*h) \
    ....:           + 3*g(a+h) - g(a)) / h^3
    sage: phi(h).expand()
    diff(f(a), a, a, a)

Sage example in ./sol/calculus.tex, line 57::

    sage: n = 7; x, h, a = var('x h a')
    sage: f = function('f')
    sage: g(x) = taylor(f(x), x, a, n)
    sage: phi(h) = sum(binomial(n,k)*(-1)^(n-k) \
    ....:          * g(a+k*h) for k in (0..n)) / h^n
    sage: phi(h).expand()
    diff(f(a), a, a, a, a, a, a, a)

Sage example in ./sol/calculus.tex, line 82::

    sage: theta = 12*arctan(1/38) + 20*arctan(1/57) \
    ....:       + 7*arctan(1/239) + 24*arctan(1/268)
    sage: x = tan(theta)
    sage: y = x.trig_expand()
    sage: y.trig_simplify()
    1

Sage example in ./sol/calculus.tex, line 94::

    sage: M = 12*(1/38)+20*(1/57)+ 7*(1/239)+24*(1/268)
    sage: M
    37735/48039

Sage example in ./sol/calculus.tex, line 113::

    sage: x = var('x')
    sage: f(x) = taylor(arctan(x), x, 0, 21)
    sage: approx = 4 * (12 * f(1/38) + 20 * f(1/57)
    ....:           + 7 * f(1/239) + 24 * f(1/268))
    sage: approx.n(digits = 50); pi.n(digits = 50)
    3.1415926535897932384626433832795028851616168852864
    3.1415926535897932384626433832795028841971693993751
    sage: approx.n(digits = 50) - pi.n(digits = 50)
    9.6444748591132486785420917537404705292978817080880e-37

Sage example in ./sol/calculus.tex, line 143::

    sage: n = var('n')
    sage: phi = lambda x: n*pi+pi/2-arctan(1/x)
    sage: x = pi*n
    sage: for i in range(4):
    ....:    x = taylor(phi(x), n, oo, 2*i); x
    ...
    1/2*pi + pi*n
    1/2*pi + pi*n - 1/(pi*n) + 1/2/(pi*n^2)
    1/2*pi + pi*n - 1/(pi*n) + 1/2/(pi*n^2)
    - 1/12*(3*pi^2 + 8)/(pi^3*n^3) + 1/8*(pi^2 + 8)/(pi^3*n^4)
      1/2*pi + pi*n - 1/(pi*n) + 1/2/(pi*n^2)
    - 1/12*(3*pi^2 + 8)/(pi^3*n^3) + 1/8*(pi^2 + 8)/(pi^3*n^4)
    - 1/240*(15*pi^4 + 240*pi^2 + 208)/(pi^5*n^5)
    + 1/96*(3*pi^4 + 80*pi^2 + 208)/(pi^5*n^6)

Sage example in ./sol/calculus.tex, line 192::

    sage: h = var('h')
    sage: f(x, y) = x * y * (x**2 - y**2) / (x**2 + y**2)
    sage: D1f(x, y) = diff(f(x,y), x)
    sage: limit((D1f(0,h) - 0) / h, h=0)
    -1
    sage: D2f(x, y) = diff(f(x,y), y)
    sage: limit((D2f(h,0) - 0) / h, h=0)
    1
    sage: g = plot3d(f(x, y), (x, -3, 3), (y, -3, 3))

Sage example in ./sol/calculus.tex, line 230::

    sage: n, t = var('n, t')
    sage: v(n)=(4/(8*n+1)-2/(8*n+4)-1/(8*n+5)-1/(8*n+6))*1/16^n
    sage: assume(8*n+1>0)
    sage: u(n) = integrate((4*sqrt(2)-8*t^3-4*sqrt(2)*t^4\
    ....:                  -8*t^5) * t^(8*n), t, 0, 1/sqrt(2))
    sage: (u(n)-v(n)).canonicalize_radical()
    0

Sage example in ./sol/calculus.tex, line 258::

    sage: t = var('t')
    sage: J = integrate((4*sqrt(2)-8*t^3 \
    ....:      - 4*sqrt(2)*t^4-8*t^5)\
    ....:      / (1-t^8), t, 0, 1/sqrt(2))
    sage: J.canonicalize_radical()
    pi + 2*log(sqrt(2) + 1) + 2*log(sqrt(2) - 1)

Sage example in ./sol/calculus.tex, line 272::

    sage: ln(exp(J).simplify_log())
    pi

Sage example in ./sol/calculus.tex, line 281::

    sage: l = sum(v(n) for n in (0..40)); l.n(digits=60)
    3.14159265358979323846264338327950288419716939937510581474759
    sage: pi.n(digits=60)
    3.14159265358979323846264338327950288419716939937510582097494
    sage: print("%e" % (l-pi).n(digits=60))
    -6.227358e-54

Sage example in ./sol/calculus.tex, line 302::

    sage: X = var('X')
    sage: ps = lambda f,g : integral(f * g, X, -pi, pi)
    sage: n = 5; Q = sin(X)
    sage: a, a0, a1, a2, a3, a4, a5 = var('a a0 a1 a2 a3 a4 a5')
    sage: a= [a0, a1, a2, a3, a4, a5]
    sage: P = sum(a[k] * X^k for k in (0..n))
    sage: equ = [ps(P - Q, X^k) for k in (0..n)]
    sage: sol = solve(equ, a)
    sage: P = sum(sol[0][k].rhs() * X^k for k in (0..n))
    sage: g = plot(P,X,-6,6,color='red') + plot(Q,X,-6,6,color='blue')

Sage example in ./sol/calculus.tex, line 353::

    sage: p, e = var('p e')
    sage: theta1, theta2, theta3 = var('theta1 theta2 theta3')
    sage: r(theta) = p / (1-e * cos(theta))
    sage: r1 = r(theta1); r2 = r(theta2); r3 = r(theta3)
    sage: R1 = vector([r1 * cos(theta1), r1 * sin(theta1), 0])
    sage: R2 = vector([r2 * cos(theta2), r2 * sin(theta2), 0])
    sage: R3 = vector([r3 * cos(theta3), r3 * sin(theta3), 0])

Sage example in ./sol/calculus.tex, line 365::

    sage: D = R1.cross_product(R2) + R2.cross_product(R3) \
    ....:   + R3.cross_product(R1)
    sage: i = vector([1, 0, 0])
    sage: S = (r1 - r3) * R2 + (r3 - r2) * R1 +   (r2 - r1) * R3
    sage: V =  S + e * i.cross_product(D)
    sage: [x.simplify_full() for x in V]
    [0, 0, 0]

Sage example in ./sol/calculus.tex, line 390::

    sage: N = r3 * R1.cross_product(R2) + r1 * R2.cross_product(R3)\
    ....:   + r2 * R3.cross_product(R1)
    sage: W =  p * S + e * i.cross_product(N)
    sage: [x.simplify_full() for x in W]
    [0, 0, 0]

Sage example in ./sol/calculus.tex, line 409::

    sage: R1=vector([0,1.,0]);R2=vector([2.,2.,0]);R3=vector([3.5,0,0])
    sage: r1 = R1.norm(); r2 = R2.norm(); r3 = R3.norm()
    sage: D = R1.cross_product(R2) + R2.cross_product(R3) \
    ....:   + R3.cross_product(R1)
    sage: S = (r1 - r3) * R2 + (r3 - r2) * R1 + (r2 - r1) * R3
    sage: V =  S + e * i.cross_product(D)
    sage: N = r3 * R1.cross_product(R2) + r1 * R2.cross_product(R3) \
    ....:   + r2 * R3.cross_product(R1)
    sage: i = vector([1, 0, 0]); W =  p * S + e * i.cross_product(N)
    sage: e = S.norm() / D.norm(); p = N.norm() / D.norm()
    sage: a = p/(1-e^2); c = a * e; b = sqrt(a^2 - c^2)
    sage: X = S.cross_product(D); i = X / X.norm()
    sage: phi = atan2(i[1],i[0]) * 180 / pi.n()
    sage: print("%.3f %.3f %.3f %.3f %.3f %.3f" % (a, b, c, e, p, phi))
    2.360 1.326 1.952 0.827 0.746 17.917

Sage example in ./sol/calculus.tex, line 445::

    sage: A = matrix(QQ, [[2, -3, 2, -12, 33],
    ....:                 [ 6, 1, 26, -16, 69],
    ....:                 [10, -29, -18, -53, 32],
    ....:                 [2, 0, 8, -18, 84]])
    sage: A.right_kernel()
    Vector space of degree 5 and dimension 2 over Rational Field
    Basis matrix:
    [     1      0  -7/34   5/17   1/17]
    [     0      1  -3/34 -10/17  -2/17]

Sage example in ./sol/calculus.tex, line 463::

    sage: H = A.echelon_form()

Sage example in ./sol/calculus.tex, line 484::

    sage: A.column_space()
    Vector space of degree 4 and dimension 3 over Rational Field
    Basis matrix:
    [       1        0        0 1139/350]
    [       0        1        0    -9/50]
    [       0        0        1   -12/35]

Sage example in ./sol/calculus.tex, line 496::

    sage: S.<x, y, z, t>=QQ[]
    sage: C = matrix(S, 4, 1, [x, y, z, t])
    sage: B = block_matrix([A, C], ncols=2)
    sage: C = B.echelon_form()
    sage: C[3,5]*350
    -1139*x + 63*y + 120*z + 350*t

Sage example in ./sol/calculus.tex, line 511::

    sage: K = A.kernel(); K
    Vector space of degree 4 and dimension 1 over Rational Field
    Basis matrix:
    [        1  -63/1139 -120/1139 -350/1139]

Sage example in ./sol/calculus.tex, line 519::

    sage: matrix(K.0).right_kernel()
    Vector space of degree 4 and dimension 3 over Rational Field
    Basis matrix:
    [       1        0        0 1139/350]
    [       0        1        0    -9/50]
    [       0        0        1   -12/35]

Sage example in ./sol/calculus.tex, line 533::

    sage: A = matrix(QQ, [[-2, 1, 1], [8, 1, -5], [4, 3, -3]])
    sage: C = matrix(QQ, [[1, 2, -1], [2, -1, -1], [-5, 0, 3]])

Sage example in ./sol/calculus.tex, line 540::

    sage: B = C.solve_left(A); B
    [ 0 -1  0]
    [ 2  3  0]
    [ 2  1  0]

Sage example in ./sol/calculus.tex, line 548::

    sage: C.left_kernel()
    Vector space of degree 3 and dimension 1 over Rational Field
    Basis matrix:
    [1 2 1]

Sage example in ./sol/calculus.tex, line 560::

    sage: x, y, z = var('x, y, z'); v = matrix([[1, 2, 1]])
    sage: B = B+(x*v).stack(y*v).stack(z*v); B
    [      x 2*x - 1       x]
    [  y + 2 2*y + 3       y]
    [  z + 2 2*z + 1       z]

Sage example in ./sol/calculus.tex, line 568::

    sage: A == B*C
    True

"""
# This file was *autogenerated* from the file calculus_doctest.sage.
