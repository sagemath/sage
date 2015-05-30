## -*- encoding: utf-8 -*-
"""
Doctests from French Sage book
Test file for chapter "Polynômes" ("Polynomials")

Sage example in ./polynomes.tex, line 55 (in svn rev 1261)::

  sage: x = var('x'); p = (2*x+1)*(x+2)*(x^4-1)
  sage: print p, "est de degré", p.degree(x)
  (x^4 - 1)*(2*x + 1)*(x + 2) est de degré 6

Sage example in ./polynomes.tex, line 69::

  sage: x = polygen(QQ, 'x'); p = (2*x+1)*(x+2)*(x^4-1)
  sage: print p, "est de degré", p.degree()
  2*x^6 + 5*x^5 + 2*x^4 - 2*x^2 - 5*x - 2 est de degré 6

Sage example in ./polynomes.tex, line 97::

  sage: R = PolynomialRing(QQ, 'x')
  sage: x = R.gen()

Sage example in ./polynomes.tex, line 104::

  sage: x.parent()
  Univariate Polynomial Ring in x over Rational Field

Sage example in ./polynomes.tex, line 125::

  sage: x = polygen(QQ, 'y'); y = polygen(QQ, 'x')

Sage example in ./polynomes.tex, line 128::

  sage: x^2 + 1
  y^2 + 1
  sage: (y^2 + 1).parent()
  Univariate Polynomial Ring in x over Rational Field

Sage example in ./polynomes.tex, line 135::

  sage: Q.<x> = QQ[]; p = x + 1; x = 2; p = p + x

Sage example in ./polynomes.tex, line 154::

  sage: R.<x,y,z,t> = QQ[]; p = (x+y+z*t)^2
  sage: p.polynomial(t).reverse()
  (x^2 + 2*x*y + y^2)*t^2 + (2*x*z + 2*y*z)*t + z^2

Sage example in ./polynomes.tex, line 162::

  sage: x = polygen(QQ); y = polygen(QQ[x], 'y')
  sage: p = x^3 + x*y + y + y^2; p
  y^2 + (x + 1)*y + x^3
  sage: q = QQ['x,y'](p); q
  x^3 + x*y + y^2 + y
  sage: QQ['x']['y'](q)
  y^2 + (x + 1)*y + x^3

Sage example in ./polynomes.tex, line 217::

  sage: def rook_polynomial(n, var='x'):
  ...       return ZZ[var]([binomial(n, k)^2 * factorial(k)
  ...                                      for k in (0..n) ])
  ...

Sage example in ./polynomes.tex, line 259::

  sage: x = polygen(QQ)
  sage: p = x^2 - 16*x + 3
  sage: p.factor()
  x^2 - 16*x + 3
  sage: p.change_ring(RDF).factor()
  (x - 15.8102496...) * (x - 0.189750324...)

Sage example in ./polynomes.tex, line 270::

  sage: p.change_ring(GF(3))
  x^2 + 2*x

Sage example in ./polynomes.tex, line 317::

  sage: list(GF(2)['x'].polynomials(of_degree=2))
  [x^2, x^2 + 1, x^2 + x, x^2 + x + 1]

Sage example in ./polynomes.tex, line 326::

  sage: A = QQ['x']
  sage: A.is_ring() and A.is_noetherian()
  True

Sage example in ./polynomes.tex, line 332::

  sage: ZZ.is_subring(A)
  True
  sage: [n for n in range(20)
  ...       if Integers(n)['x'].is_integral_domain()]
  [0, 2, 3, 5, 7, 11, 13, 17, 19]

Sage example in ./polynomes.tex, line 395::

  sage: R.<t> = Integers(42)[]; (t^20-1) % (t^5+8*t+7)
  22*t^4 + 14*t^3 + 14*t + 6

Sage example in ./polynomes.tex, line 407::

  sage: ((t^2+t)//t).parent()
  Univariate Polynomial Ring in t over Ring of integers modulo 42
  sage: (t^2+t)/t
  Traceback (most recent call last):
  ...
  TypeError: self must be an integral domain.

Sage example in ./polynomes.tex, line 420::

  sage: x = polygen(QQ); [chebyshev_T(n, x) for n in (0..4)]
  [1, x, 2*x^2 - 1, 4*x^3 - 3*x, 8*x^4 - 8*x^2 + 1]

Sage example in ./polynomes.tex, line 436::

  sage: S.<x> = ZZ[]; p = 2*(x^10-1)*(x^8-1)
  sage: p.gcd(p.derivative())
  2*x^2 - 2

Sage example in ./polynomes.tex, line 447::

  sage: R.<x> = QQ[]; p = (x^5-1); q = (x^3-1)
  sage: print "le pgcd est %s = (%s)*p + (%s)*q" % p.xgcd(q)
  le pgcd est x - 1 = (-x)*p + (x^3 + 1)*q

Sage example in ./polynomes.tex, line 484::

  sage: R.<x> = QQ[]; p = 3*x^2 - 6
  sage: p.is_irreducible(), p.change_ring(ZZ).is_irreducible()
  (True, False)

Sage example in ./polynomes.tex, line 527::

  sage: x = polygen(ZZ); p = 54*x^4+36*x^3-102*x^2-72*x-12
  sage: p.factor()
  2 * 3 * (3*x + 1)^2 * (x^2 - 2)

Sage example in ./polynomes.tex, line 533::

  sage: for A in [QQ, ComplexField(16), GF(5), QQ[sqrt(2)]]:
  ...       print A, ":"; print A['x'](p).factor()
  Rational Field :
  (54) * (x + 1/3)^2 * (x^2 - 2)
  Complex Field with 16 bits of precision :
  (54.00) * (x - 1.414) * (x + 0.3333)^2 * (x + 1.414)
  Finite Field of size 5 :
  (4) * (x + 2)^2 * (x^2 + 3)
  Number Field in sqrt2 with defining polynomial x^2 - 2 :
  (54) * (x - sqrt2) * (x + sqrt2) * (x + 1/3)^2

Sage example in ./polynomes.tex, line 601::

  sage: R.<x> = ZZ[]; p = (2*x^2-5*x+2)^2 * (x^4-7); p.roots()
  [(2, 2)]

Sage example in ./polynomes.tex, line 608::

  sage: p.roots(QQ)
  [(2, 2), (1/2, 2)]
  sage: p.roots(Zp(19, print_max_terms=3))
  [(7 + 16*19 + 17*19^2 + ... + O(19^20), 1),
  (12 + 2*19 + 19^2 + ... + O(19^20), 1),
  (10 + 9*19 + 9*19^2 + ... + O(19^20), 2),
  (2 + O(19^20), 2)]

Sage example in ./polynomes.tex, line 623::

  sage: racines = p.roots(AA); racines
  [(-1.626576561697786?, 1), (0.500000000000000?, 2),
  (1.626576561697786?, 1), (2.000000000000000?, 2)]

Sage example in ./polynomes.tex, line 629::

  sage: a = racines[0][0]^4; a.simplify(); a
  7

Sage example in ./polynomes.tex, line 646::

  sage: x = polygen(ZZ); (x-12).resultant(x-20)
  -8

Sage example in ./polynomes.tex, line 701::

  sage: R.<a,b,c,d> = QQ[]; x = polygen(R); p = a*x^2+b*x+c
  sage: p.resultant(p.derivative())
  -a*b^2 + 4*a^2*c
  sage: p.discriminant()
  b^2 - 4*a*c
  sage: (a*x^3 + b*x^2 + c*x + d).discriminant()
  b^2*c^2 - 4*a*c^3 - 4*b^3*d + 18*a*b*c*d - 27*a^2*d^2

Sage example in ./polynomes.tex, line 770::

  sage: R.<x> = QQ[]
  sage: J1 = (x^2 - 2*x + 1, 2*x^2 + x - 3)*R; J1
  Principal ideal (x - 1) of Univariate Polynomial Ring in x
  over Rational Field

Sage example in ./polynomes.tex, line 778::

  sage: J2 = R.ideal(x^5 + 2)
  sage: ((3*x+5)*J1*J2).reduce(x^10)
  421/81*x^6 - 502/81*x^5 + 842/81*x - 680/81

Sage example in ./polynomes.tex, line 785::

  sage: B = R.quo((3*x+5)*J1*J2) # quo nomme automatiquement 'xbar'
  sage: B(x^10)                  #   le générateur de B image de x
  421/81*xbar^6 - 502/81*xbar^5 + 842/81*xbar - 680/81
  sage: B(x^10).lift()
  421/81*x^6 - 502/81*x^5 + 842/81*x - 680/81

Sage example in ./polynomes.tex, line 856::

  sage: x = polygen(RR); r = (1 + x)/(1 - x^2); r.parent()
  Fraction Field of Univariate Polynomial Ring in x over Real
  Field with 53 bits of precision
  sage: r
  (x + 1.00000000000000)/(-x^2 + 1.00000000000000)

Sage example in ./polynomes.tex, line 864::

  sage: r.reduce(); r
  1.00000000000000/(-x + 1.00000000000000)

Sage example in ./polynomes.tex, line 907::

  sage: R.<x> = QQ[]; r = x^10/((x^2-1)^2*(x^2+3))
  sage: poly, parts = r.partial_fraction_decomposition()
  sage: print poly
  x^4 - x^2 + 6
  sage: for part in parts: print part.factor()
  (17/32) * (x - 1)^-1
  (1/16) * (x - 1)^-2
  (-17/32) * (x + 1)^-1
  (1/16) * (x + 1)^-2
  (-243/16) * (x^2 + 3)^-1

Sage example in ./polynomes.tex, line 931::

  sage: C = ComplexField(15)
  sage: Frac(C['x'])(r).partial_fraction_decomposition()
  (x^4 - x^2 + 6.000,
  [0.5312/(x - 1.000),
   0.06250/(x^2 - 2.000*x + 1.000),
    4.385*I/(x - 1.732*I),
    (-4.385*I)/(x + 1.732*I),
    (-0.5312)/(x + 1.000),
    0.06250/(x^2 + 2.000*x + 1.000)])

Sage example in ./polynomes.tex, line 966::

  sage: A = Integers(101); R.<x> = A[]
  sage: f6 = sum( (i+1)^2 * x^i for i in (0..5) ); f6
  36*x^5 + 25*x^4 + 16*x^3 + 9*x^2 + 4*x + 1
  sage: num, den = f6.rational_reconstruct(x^6, 1, 3); num/den
  (100*x + 100)/(x^3 + 98*x^2 + 3*x + 100)

Sage example in ./polynomes.tex, line 974::

  sage: S = PowerSeriesRing(A, 'x', 7); S(num)/S(den)
  1 + 4*x + 9*x^2 + 16*x^3 + 25*x^4 + 36*x^5 + 49*x^6 + O(x^7)

Sage example in ./polynomes.tex, line 1015::

  sage: x = var('x'); s = tan(x).taylor(x, 0, 20)
  sage: p = previous_prime(2^30); ZpZx = Integers(p)['x']
  sage: Qx = QQ['x']

Sage example in ./polynomes.tex, line 1020::

  sage: num, den = ZpZx(s).rational_reconstruct(ZpZx(x)^10,4,5)
  sage: num/den
  (1073741779*x^3 + 105*x)/(x^4 + 1073741744*x^2 + 105)

Sage example in ./polynomes.tex, line 1026::

  sage: def lift_sym(a):
  ...       m = a.parent().defining_ideal().gen()
  ...       n = a.lift()
  ...       if n <= m // 2: return n
  ...       else: return n - m

Sage example in ./polynomes.tex, line 1034::

  sage: Qx(map(lift_sym, num))/Qx(map(lift_sym, den))
  (-10*x^3 + 105*x)/(x^4 - 45*x^2 + 105)

Sage example in ./polynomes.tex, line 1042::

  sage: def mypade(pol, n, k):
  ...       x = ZpZx.gen();
  ...       n,d = ZpZx(pol).rational_reconstruct(x^n, k-1, n-k)
  ...       return Qx(map(lift_sym, n))/Qx(map(lift_sym, d))

Sage example in ./polynomes.tex, line 1109::

  sage: R.<x> = PowerSeriesRing(QQ)

Sage example in ./polynomes.tex, line 1126::

  sage: R.<x> = QQ[[]]
  sage: f = 1 + x + O(x^2); g = x + 2*x^2 + O(x^4)
  sage: f + g
  1 + 2*x + O(x^2)
  sage: f * g
  x + 3*x^2 + O(x^3)

Sage example in ./polynomes.tex, line 1136::

  sage: (1 + x^3).prec()
  +Infinity

Sage example in ./polynomes.tex, line 1141::

  sage: R.<x> = PowerSeriesRing(Reals(24), default_prec=4)
  sage: 1/(1 + RR.pi() * x)^2
  1.00000 - 6.28319*x + 29.6088*x^2 - 124.025*x^3 + O(x^4)

Sage example in ./polynomes.tex, line 1148::

  sage: R.<x> = QQ[[]]
  sage: 1 + x + O(x^2) == 1 + x + x^2 + O(x^3)
  True

Sage example in ./polynomes.tex, line 1158::

  sage: (1/(1+x)).sqrt().integral().exp()/x^2 + O(x^4)
  x^-2 + x^-1 + 1/4 + 1/24*x - 1/192*x^2 + 11/1920*x^3 + O(x^4)

Sage example in ./polynomes.tex, line 1178::

  sage: (1+x^2).sqrt().solve_linear_de(prec=6, b=x.exp())
  1 + 2*x + 3/2*x^2 + 5/6*x^3 + 1/2*x^4 + 7/30*x^5 + O(x^6)

Sage example in ./polynomes.tex, line 1186::

  sage: S.<x> = PowerSeriesRing(QQ, default_prec=5)
  sage: f = S(1)
  sage: for i in range(5):
  ...       f = (x*f).exp()
  ...       print f
  1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + O(x^5)
  1 + x + 3/2*x^2 + 5/3*x^3 + 41/24*x^4 + O(x^5)
  1 + x + 3/2*x^2 + 8/3*x^3 + 101/24*x^4 + O(x^5)
  1 + x + 3/2*x^2 + 8/3*x^3 + 125/24*x^4 + O(x^5)
  1 + x + 3/2*x^2 + 8/3*x^3 + 125/24*x^4 + O(x^5)

Sage example in ./polynomes.tex, line 1233::

  sage: L.<x> = LazyPowerSeriesRing(QQ)
  sage: lazy_exp = x.exponential(); lazy_exp
  O(1)

Sage example in ./polynomes.tex, line 1239::

  sage: lazy_exp[5]
  1/120
  sage: lazy_exp
  1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + O(x^6)

Sage example in ./polynomes.tex, line 1247::

  sage: f = L(1)  # la série paresseuse constante 1
  sage: for i in range(5):
  ...       f = (x*f).exponential()
  ...       f.compute_coefficients(5)  # force le calcul des
  ...       print f                    #   premiers coefficients
  1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + O(x^6)
  1 + x + 3/2*x^2 + 5/3*x^3 + 41/24*x^4 + 49/30*x^5 + O(x^6)
  1 + x + 3/2*x^2 + 8/3*x^3 + 101/24*x^4 + 63/10*x^5 + O(x^6)
  1 + x + 3/2*x^2 + 8/3*x^3 + 125/24*x^4 + 49/5*x^5 + O(x^6)
  1 + x + 3/2*x^2 + 8/3*x^3 + 125/24*x^4 + 54/5*x^5 + O(x^6)

Sage example in ./polynomes.tex, line 1264::

  sage: f[7]
  28673/630

Sage example in ./polynomes.tex, line 1271::

  sage: from sage.combinat.species.series import LazyPowerSeries
  sage: f = LazyPowerSeries(L, name='f')
  sage: f.define((x*f).exponential())
  sage: f.coefficients(8)
  [1, 1, 3/2, 8/3, 125/24, 54/5, 16807/720, 16384/315]

Sage example in ./polynomes.tex, line 1309::

  sage: R = PolynomialRing(ZZ, 'x', sparse=True)
  sage: p = R.cyclotomic_polynomial(2^50); p, p.derivative()
  (x^562949953421312 + 1, 562949953421312*x^562949953421311)

"""
from sage.all_cmdline import *   # import sage library
