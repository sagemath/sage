r"""
Tests for the Sage <-> PARI interface

The default precision is 64 bits, see :issue:`21425`::

    sage: pari("bitprecision(Pi)")
    64

Sage-specific API checks:

Creating PARI objects::

    sage: pari(Matrix(2,2,range(4)))                                                    # needs sage.modules
    [0, 1; 2, 3]
    sage: pari(x^2-3)                                                                   # needs sage.symbolic
    x^2 - 3

The following example caused Sage to crash before
:issue:`20630`::

    sage: R.<theta> = QQ[]
    sage: K.<a> = NumberField(theta^2 + 1)                                              # needs sage.rings.number_field
    sage: K.absolute_polynomial().galois_group(pari_group=True)                         # needs sage.groups sage.rings.number_field
    PARI group [2, -1, 1, "S2"] of degree 2

Before :issue:`15654`, this used to take a very long time.
Now it takes much less than a second::

    sage: pari.allocatemem(200000)
    PARI stack size set to 200000 bytes, maximum size set to ...
    sage: x = polygen(ZpFM(3,10))                                                       # needs sage.rings.padics
    sage: pol = ((x-1)^50 + x)                                                          # needs sage.rings.padics
    sage: pari(pol).poldisc()                                                           # needs sage.rings.padics
    2*3 + 3^4 + 2*3^6 + 3^7 + 2*3^8 + 2*3^9 + O(3^10)

This used to give the wrong answer before :issue:`23259`::

    sage: R.<x> = QQ[]
    sage: f = pari(x^12 + x^7 - 1/5*x^6 - 3*x^5 + 13/5*x^4 + 11/5*x^3 + 2/5*x^2 + 2/5*x + 1/5)
    sage: g,h = f.polredabs(1)
    sage: f.subst(x,h)
    Mod(0, x^12 - 2*x^11 + 2*x^10 - 11*x^9 + 13*x^8 + 15*x^7 - x^6 - 5*x^5 + 5)

Getting the coefficients of a Laurent series behaves differently
in Sage and PARI. In PARI we get all coefficients starting
from the lowest degree term.  This includes trailing zeros::

    sage: R.<x> = LaurentSeriesRing(QQ)
    sage: s = x^2 + O(x^8)
    sage: s.list()
    [1]
    sage: pari(s).list()
    [1, 0, 0, 0, 0, 0]
    sage: s = x^-2 + O(x^0)
    sage: s.list()
    [1]
    sage: pari(s).list()
    [1, 0]

Number fields::

    sage: # needs sage.rings.number_field
    sage: x = polygen(QQ)
    sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
    sage: pari(K).nf_get_pol()
    y^4 - 4*y^2 + 1
    sage: L.<b> = K.extension(x^2 - 5)
    sage: pari(L).nf_get_pol()              # Absolute
    y^8 - 28*y^6 + 208*y^4 - 408*y^2 + 36
    sage: L.pari_rnf().nf_get_pol()         # Relative
    x^2 - 5
    sage: K.pari_nf().nf_get_pol()
    y^4 - 4*y^2 + 1
    sage: K.pari_bnf().nf_get_pol()
    y^4 - 4*y^2 + 1

    sage: # needs sage.rings.number_field
    sage: K.<a> = QuadraticField(-65)
    sage: G = K.pari_bnf().bnf_get_gen(); G
    [[3, 2; 0, 1], [2, 1; 0, 1]]
    sage: [K.ideal(J) for J in G]
    [Fractional ideal (3, a + 2), Fractional ideal (2, a + 1)]

Conversions::

    sage: # needs sage.rings.number_field
    sage: K.<i> = QuadraticField(-1)
    sage: F = pari(K).idealfactor(K.ideal(5)); F
    [[5, [-2, 1]~, 1, 1, [2, -1; 1, 2]], 1; [5, [2, 1]~, 1, 1, [-2, -1; 1, -2]], 1]
    sage: F[0,0].pr_get_p()
    5
    sage: J = pari(K).idealstar(K.ideal(4*i + 2))
    sage: J.bid_get_cyc()
    [4, 2]

    sage: int(pari(RealField(63)(2^63 - 1)))                                            # needs sage.rings.real_mpfr
    9223372036854775807
    sage: int(pari(RealField(63)(2^63 + 2)))                                            # needs sage.rings.real_mpfr
    9223372036854775810

    sage: # needs sage.rings.padics
    sage: K = Qp(11,5)
    sage: x = K(11^-10 + 5*11^-7 + 11^-6)
    sage: y = pari(x)
    sage: y.padicprime()
    11
    sage: y.padicprime().type()
    't_INT'

    sage: # needs sage.rings.finite_rings
    sage: x = polygen(GF(3))
    sage: k.<a> = GF(9, modulus=x^2 + 1)
    sage: b = pari(a).ffprimroot()
    sage: b  # random
    a + 1
    sage: b.fforder()
    8

    sage: pari(4).Zn_issquare(30.factor())
    True
    sage: pari(4).Zn_sqrt(30.factor())
    22

    sage: a = pari(1/2); a, a.type()
    (1/2, 't_FRAC')

Conversion from matrices and vectors is supported::

    sage: a = pari(matrix(2,3,[1,2,3,4,5,6])); a, a.type()                              # needs sage.modules
    ([1, 2, 3; 4, 5, 6], 't_MAT')
    sage: v = vector([1.2, 3.4, 5.6])                                                   # needs sage.modules
    sage: pari(v)                                                                       # needs sage.modules
    [1.20000000000000, 3.40000000000000, 5.60000000000000]

Some more exotic examples::

    sage: # needs sage.rings.number_field
    sage: K.<a> = NumberField(polygen(QQ)^3 - 2)
    sage: pari(K)
    [y^3 - 2, [1, 1], -108, 1, [[1, 1.25992104989487, 1.58740105196820; 1, -0.629960524947437 + 1.09112363597172*I, -0.793700525984100 - 1.37472963699860*I], [1, 1.25992104989487, 1.58740105196820; 1, 0.461163111024285, -2.16843016298270; 1, -1.72108416091916, 0.581029111014503], [16, 20, 25; 16, 7, -35; 16, -28, 9], [3, 0, 0; 0, 0, 6; 0, 6, 0], [6, 0, 0; 0, 6, 0; 0, 0, 3], [2, 0, 0; 0, 0, 1; 0, 1, 0], [2, [0, 0, 2; 1, 0, 0; 0, 1, 0]], [2, 3]], [1.25992104989487, -0.629960524947437 + 1.09112363597172*I], [1, y, y^2], [1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 0, 0, 0, 0, 2, 0, 2, 0; 0, 1, 0, 1, 0, 0, 0, 0, 2; 0, 0, 1, 0, 1, 0, 1, 0, 0]]

    sage: E = EllipticCurve('37a1')                                                     # needs sage.schemes
    sage: pari(E)                                                                       # needs sage.schemes
    [0, 0, 1, -1, 0, 0, -2, 1, -1, 48, -216, 37, 110592/37, Vecsmall([1]), [Vecsmall([64, 1])], [0, 0, 0, 0, 0, 0, 0, 0]]

Deprecation checks::

    sage: pari.poltchebi(10)
    doctest:...: DeprecationWarning: the PARI/GP function poltchebi is obsolete (2013-04-03)
    512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
    sage: pari("x^3 + 1").polsturm(-1, 1)
    doctest:...: DeprecationWarning: argument 2 of the PARI/GP function polsturm is undocumented and deprecated
    1
    sage: x = pari('10^100')
    sage: x.Str().length()
    101
    sage: x.sizedigit()
    doctest:...: DeprecationWarning: the PARI/GP function sizedigit is obsolete (2015-01-13)
    101
    sage: x = pari('1.234')
    sage: x
    1.23400000000000
    sage: x.sizedigit()
    1
    sage: pari('7234.1').sizedigit()
    4
    sage: pari('9234.1').sizedigit()
    5
    sage: [pari(2*n).bernfrac() for n in range(9)]
    [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
    sage: e = pari([0,1,1,-2,0]).ellinit()
    sage: e.elllseries(2.1)
    doctest:...: DeprecationWarning: the PARI/GP function elllseries is obsolete (2016-08-08)
    0.402838047956645
    sage: e.elllseries(1, precision=128)
    -6.17606670058278 E-39
    sage: e.elllseries(1, precision=256)
    -2.05598131842639 E-77
    sage: e.elllseries(-2)
    0
    sage: e.elllseries(2.1, A=1.1)
    0.402838047956645

A long list of doctests which used to be part of manually written code
which is now automatically generated:

Reading a gp file::

    sage: import tempfile
    sage: gpfile = tempfile.NamedTemporaryFile(mode='w')
    sage: __ = gpfile.file.write("mysquare(n) = {\n")
    sage: __ = gpfile.file.write("    n^2;\n")
    sage: __ = gpfile.file.write("}\n")
    sage: __ = gpfile.file.write("polcyclo(5)\n")
    sage: gpfile.file.flush()
    sage: pari.read(gpfile.name)
    x^4 + x^3 + x^2 + x + 1
    sage: pari('mysquare(12)')
    144

Constants::

    sage: pari.euler()
    0.577215664901533
    sage: pari.euler(precision=100).sage()
    0.577215664901532860606512090082...
    sage: pari.pi()
    3.14159265358979
    sage: pari.pi(precision=100).sage()
    3.1415926535897932384626433832...

Polynomial functions::

    sage: R.<x> = PolynomialRing(ZZ)
    sage: pari(2*x^2 + 2).content()
    2
    sage: pari("4*x^3 - 2*x/3 + 2/5").content()
    2/15

    sage: x = pari('y^8+6*y^6-27*y^5+1/9*y^2-y+1')
    sage: x.newtonpoly(3)
    [1, 1, -1/3, -1/3, -1/3, -1/3, -1/3, -1/3]

    sage: f = pari("x^2 + y^3 + x*y")
    sage: f
    x^2 + y*x + y^3
    sage: f.polcoef(1)
    y
    sage: f.polcoef(3)
    0
    sage: f.polcoef(3, "y")
    1
    sage: f.polcoef(1, "y")
    x

    sage: pari("x^2 + 1").poldisc()
    -4

    sage: pari.pollegendre(7)
    429/16*x^7 - 693/16*x^5 + 315/16*x^3 - 35/16*x
    sage: pari.pollegendre(7, 'z')
    429/16*z^7 - 693/16*z^5 + 315/16*z^3 - 35/16*z
    sage: pari.pollegendre(0)
    1

    sage: pari.polcyclo(8)
    x^4 + 1
    sage: pari.polcyclo(7, 'z')
    z^6 + z^5 + z^4 + z^3 + z^2 + z + 1
    sage: pari.polcyclo(1)
    x - 1

Power series::

    sage: f = pari('x+x^2+x^3+O(x^4)'); f
    x + x^2 + x^3 + O(x^4)
    sage: g = f.serreverse(); g
    x - x^2 + x^3 + O(x^4)
    sage: f.subst('x',g)
    x + O(x^4)
    sage: g.subst('x',f)
    x + O(x^4)

Random seed::

    sage: a = pari.getrand()
    sage: a.type()
    't_INT'

Constructors::

    sage: v = pari([1,2,3])
    sage: v
    [1, 2, 3]
    sage: v.type()
    't_VEC'
    sage: w = v.List()
    sage: w
    List([1, 2, 3])
    sage: w.type()
    't_LIST'

    sage: x = pari(5)
    sage: x.type()
    't_INT'
    sage: y = x.Mat()
    sage: y
    Mat(5)
    sage: y.type()
    't_MAT'
    sage: x = pari('[1,2;3,4]')
    sage: x.type()
    't_MAT'
    sage: x = pari('[1,2,3,4]')
    sage: x.type()
    't_VEC'
    sage: y = x.Mat()
    sage: y
    Mat([1, 2, 3, 4])
    sage: y.type()
    't_MAT'

    sage: v = pari('[1,2;3,4]').Vec(); v
    [[1, 3]~, [2, 4]~]
    sage: v.Mat()
    [1, 2; 3, 4]
    sage: v = pari('[1,2;3,4]').Col(); v
    [[1, 2], [3, 4]]~
    sage: v.Mat()
    [1, 2; 3, 4]

    sage: z = pari(3)
    sage: x = z.Mod(pari(7))
    sage: x
    Mod(3, 7)
    sage: x^2
    Mod(2, 7)
    sage: x^100
    Mod(4, 7)
    sage: x.type()
    't_INTMOD'
    sage: f = pari("x^2 + x + 1")
    sage: g = pari("x")
    sage: a = g.Mod(f)
    sage: a
    Mod(x, x^2 + x + 1)
    sage: a*a
    Mod(-x - 1, x^2 + x + 1)
    sage: a.type()
    't_POLMOD'

    sage: v = pari("[1,2,3,4]")
    sage: f = v.Pol()
    sage: f
    x^3 + 2*x^2 + 3*x + 4
    sage: f*f
    x^6 + 4*x^5 + 10*x^4 + 20*x^3 + 25*x^2 + 24*x + 16

    sage: v = pari("[1,2;3,4]")
    sage: v.Pol()
    [1, 3]~*x + [2, 4]~

    sage: v = pari("[1,2,3,4]")
    sage: f = v.Polrev()
    sage: f
    4*x^3 + 3*x^2 + 2*x + 1
    sage: v.Pol()
    x^3 + 2*x^2 + 3*x + 4
    sage: v.Polrev('y')
    4*y^3 + 3*y^2 + 2*y + 1

    sage: f
    4*x^3 + 3*x^2 + 2*x + 1
    sage: f.Polrev()
    4*x^3 + 3*x^2 + 2*x + 1
    sage: v = pari("[1,2;3,4]")
    sage: v.Polrev()
    [2, 4]~*x + [1, 3]~

    sage: pari(3).Qfb(7, 1)
    Qfb(3, 7, 1)
    sage: pari(3).Qfb(7, 2)
    Traceback (most recent call last):
    ...
    PariError: domain error in Qfb: issquare(disc) = 1

    sage: pari([1,5,2]).Set()
    [1, 2, 5]
    sage: pari([]).Set()
    []
    sage: pari([1,1,-1,-1,3,3]).Set()
    [-1, 1, 3]
    sage: pari(1).Set()
    [1]
    sage: pari('1/(x*y)').Set()
    [1/(y*x)]
    sage: pari('["bc","ab","bc"]').Set()
    ["ab", "bc"]

    sage: pari([65,66,123]).strchr()
    "AB{"
    sage: pari('"Sage"').Vecsmall()
    Vecsmall([83, 97, 103, 101])
    sage: _.strchr()
    "Sage"
    sage: pari([83, 97, 103, 101]).strchr()
    "Sage"

Basic functions::

    sage: pari(0).binary()
    []
    sage: pari(-5).binary()
    [1, 0, 1]
    sage: pari(5).binary()
    [1, 0, 1]
    sage: pari(2005).binary()
    [1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]
    sage: pari('"2"').binary()
    Traceback (most recent call last):
    ...
    PariError: incorrect type in binary (t_STR)

    sage: pari(1.4).ceil()
    2
    sage: pari(-1.4).ceil()
    -1
    sage: pari(3/4).ceil()
    1

    sage: # needs sage.symbolic
    sage: x = SR.symbol('x')
    sage: pari(x).ceil()
    x
    sage: pari((x^2+x+1)/x).ceil()
    x + 1
    sage: pari(x^2+5*x+2.5).ceil()
    x^2 + 5*x + 2.50000000000000

    sage: x = pari(-2).Mod(5)
    sage: x.centerlift()
    -2
    sage: x.lift()
    3
    sage: f = pari('x-1').Mod('x^2 + 1')
    sage: f.centerlift()
    x - 1
    sage: f.lift()
    x - 1
    sage: f = pari('x-y').Mod('x^2+1')
    sage: f
    Mod(x - y, x^2 + 1)
    sage: f.centerlift('x')
    x - y
    sage: f.centerlift('y')
    Mod(x - y, x^2 + 1)
    sage: pari("Mod(3,5)").lift_centered()
    -2

    sage: pari([0,1,2,3,4]).component(1)
    0
    sage: pari([0,1,2,3,4]).component(2)
    1
    sage: pari([0,1,2,3,4]).component(4)
    3
    sage: pari('x^3 + 2').component(1)
    2
    sage: pari('x^3 + 2').component(2)
    0
    sage: pari('x^3 + 2').component(4)
    1
    sage: pari('x').component(0)
    Traceback (most recent call last):
    ...
    PariError: nonexistent component: index < 1

    sage: pari('x+1').conj()
    x + 1
    sage: pari('x+I').conj()
    x - I
    sage: pari('1/(2*x+3*I)').conj()
    1/(2*x - 3*I)
    sage: pari([1,2,'2-I','Mod(x,x^2+1)', 'Mod(x,x^2-2)']).conj()
    [1, 2, 2 + I, Mod(-x, x^2 + 1), Mod(-x, x^2 - 2)]
    sage: pari('Mod(x,x^2-2)').conj()
    Mod(-x, x^2 - 2)
    sage: pari('Mod(x,x^3-3)').conj()
    Traceback (most recent call last):
    ...
    PariError: incorrect type in gconj (t_POLMOD)

    sage: pari('Mod(1+x,x^2-2)').conjvec()
    [-0.414213562373095, 2.41421356237310]~
    sage: pari('Mod(x,x^3-3)').conjvec()
    [1.44224957030741, -0.721124785153704 - 1.24902476648341*I, -0.721124785153704 + 1.24902476648341*I]~
    sage: pari('Mod(1+x,x^2-2)').conjvec(precision=192)[0].sage()
    -0.414213562373095048801688724209698078569671875376948073177

    sage: pari('5/9').denominator()
    9
    sage: pari('(x+1)/(x-2)').denominator()
    x - 2
    sage: pari('2/3 + 5/8*x + 7/3*x^2 + 1/5*y').denominator()
    1
    sage: pari('2/3*x').denominator()
    1
    sage: pari('[2/3, 5/8, 7/3, 1/5]').denominator()
    120

    sage: pari(5/9).floor()
    0
    sage: pari(11/9).floor()
    1
    sage: pari(1.17).floor()
    1
    sage: pari([1.5,2.3,4.99]).floor()
    [1, 2, 4]
    sage: pari([[1.1,2.2],[3.3,4.4]]).floor()
    [[1, 2], [3, 4]]

    sage: # needs sage.symbolic
    sage: x = SR.symbol('x')
    sage: pari(x).floor()
    x
    sage: pari((x^2+x+1)/x).floor()
    x + 1
    sage: pari(x^2+5*x+2.5).floor()
    x^2 + 5*x + 2.50000000000000

    sage: pari('"hello world"').floor()
    Traceback (most recent call last):
    ...
    PariError: incorrect type in gfloor (t_STR)

    sage: pari(1.75).frac()
    0.750000000000000
    sage: pari(sqrt(2)).frac()                                                          # needs sage.symbolic
    0.414213562373095
    sage: pari('sqrt(-2)').frac()
    Traceback (most recent call last):
    ...
    PariError: incorrect type in gfrac (t_COMPLEX)

    sage: pari('1+2*I').imag()
    2
    sage: pari(sqrt(-2)).imag()                                                         # needs sage.symbolic
    1.41421356237310
    sage: pari('x+I').imag()
    1
    sage: pari('x+2*I').imag()
    2
    sage: pari('(1+I)*x^2+2*I').imag()
    x^2 + 2
    sage: pari('[1,2,3] + [4*I,5,6]').imag()
    [4, 0, 0]

    sage: x = pari("x")
    sage: a = x.Mod('x^3 + 17*x + 3')
    sage: a
    Mod(x, x^3 + 17*x + 3)
    sage: b = a^4; b
    Mod(-17*x^2 - 3*x, x^3 + 17*x + 3)
    sage: b.lift()
    -17*x^2 - 3*x

    sage: pari(pi).sign()                                                               # needs sage.symbolic
    1
    sage: pari(0).sign()
    0
    sage: pari(-1/2).sign()
    -1
    sage: pari(SR(I)).sign()                                                            # needs sage.symbolic
    Traceback (most recent call last):
    ...
    PariError: incorrect type in gsigne (t_COMPLEX)
    sage: pari(I).sign()
    Traceback (most recent call last):
    ...
    PariError: incorrect type in gsigne (t_POLMOD)

    sage: y = pari('y')
    sage: x = pari('9') + y - y
    sage: x
    9
    sage: x.type()
    't_POL'
    sage: x.factor()
    matrix(0,2)
    sage: pari('9').factor()
    Mat([3, 2])
    sage: x.simplify()
    9
    sage: x.simplify().factor()
    Mat([3, 2])
    sage: x = pari('1.5 + 0*I')
    sage: x.type()
    't_REAL'
    sage: x.simplify()
    1.50000000000000
    sage: y = x.simplify()
    sage: y.type()
    't_REAL'

    sage: pari(2).sqr()
    4
    sage: pari("1+O(2^5)").sqr()
    1 + O(2^6)
    sage: pari("1+O(2^5)")*pari("1+O(2^5)")
    1 + O(2^5)
    sage: x = pari("1+O(2^5)"); x*x
    1 + O(2^6)

    sage: x = pari("x"); y = pari("y")
    sage: f = pari('x^3 + 17*x + 3')
    sage: f.subst(x, y)
    y^3 + 17*y + 3
    sage: f.subst(x, "z")
    z^3 + 17*z + 3
    sage: f.subst(x, "z")^2
    z^6 + 34*z^4 + 6*z^3 + 289*z^2 + 102*z + 9
    sage: f.subst(x, "x+1")
    x^3 + 3*x^2 + 20*x + 21
    sage: f.subst(x, "xyz")
    xyz^3 + 17*xyz + 3
    sage: f.subst(x, "xyz")^2
    xyz^6 + 34*xyz^4 + 6*xyz^3 + 289*xyz^2 + 102*xyz + 9

    sage: pari(9).valuation(3)
    2
    sage: pari(9).valuation(9)
    1
    sage: x = pari(9).Mod(27); x.valuation(3)
    2
    sage: pari('5/3').valuation(3)
    -1
    sage: pari('9 + 3*x + 15*x^2').valuation(3)
    1
    sage: pari([9,3,15]).valuation(3)
    1
    sage: pari('9 + 3*x + 15*x^2 + O(x^5)').valuation(3)
    1
    sage: pari('x^2*(x+1)^3').valuation(pari('x+1'))
    3
    sage: pari('x + O(x^5)').valuation('x')
    1
    sage: pari('2*x^2 + O(x^5)').valuation('x')
    2
    sage: pari(0).valuation(3)
    +oo

    sage: pari('x^2 + x -2').variable()
    x
    sage: pari('1+2^3 + O(2^5)').variable()
    2
    sage: pari('x+y0').variable()
    x
    sage: pari('y0+z0').variable()
    y0

Bitwise functions::

    sage: pari(8).bitand(4)
    0
    sage: pari(8).bitand(8)
    8
    sage: pari(6).binary()
    [1, 1, 0]
    sage: pari(7).binary()
    [1, 1, 1]
    sage: pari(6).bitand(7)
    6
    sage: pari(19).bitand(-1)
    19
    sage: pari(-1).bitand(-1)
    -1

    sage: pari(10).bitneg()
    -11
    sage: pari(1).bitneg()
    -2
    sage: pari(-2).bitneg()
    1
    sage: pari(-1).bitneg()
    0
    sage: pari(569).bitneg()
    -570
    sage: pari(569).bitneg(10)
    454
    sage: 454 % 2^10
    454
    sage: -570 % 2^10
    454

    sage: pari(14).bitnegimply(0)
    14
    sage: pari(8).bitnegimply(8)
    0
    sage: pari(8+4).bitnegimply(8)
    4

    sage: pari(14).bitor(0)
    14
    sage: pari(8).bitor(4)
    12
    sage: pari(12).bitor(1)
    13
    sage: pari(13).bitor(1)
    13

    sage: pari(6).bitxor(4)
    2
    sage: pari(0).bitxor(4)
    4
    sage: pari(6).bitxor(0)
    6

Transcendental functions::

    sage: x = pari("-27.1")
    sage: x.abs()
    27.1000000000000
    sage: pari('1 + I').abs(precision=128).sage()
    1.4142135623730950488016887242096980786
    sage: pari('x-1.2*x^2').abs()
    1.20000000000000*x^2 - x
    sage: pari('-2 + t + O(t^2)').abs()
    2 - t + O(t^2)

    sage: pari(0.5).acos()
    1.04719755119660
    sage: pari(1/2).acos()
    1.04719755119660
    sage: pari(1.1).acos()
    0.443568254385115*I
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1.1+i).acos()                                                            # needs sage.rings.real_mpfr
    0.849343054245252 - 1.09770986682533*I

    sage: pari(2).acosh()
    1.31695789692482
    sage: pari(0).acosh()
    1.57079632679490*I
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(i).acosh()                                                               # needs sage.rings.real_mpfr
    0.881373587019543 + 1.57079632679490*I

    sage: pari(2).agm(2)
    2.00000000000000
    sage: pari(0).agm(1)
    0
    sage: pari(1).agm(2)
    1.45679103104691
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).agm(-3)                                                             # needs sage.rings.real_mpfr
    -0.964731722290876 + 1.15700282952632*I

    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(2+i).arg()                                                               # needs sage.rings.real_mpfr
    0.463647609000806

    sage: pari(pari(0.5).sin()).asin()
    0.500000000000000
    sage: pari(2).asin()
    1.57079632679490 - 1.31695789692482*I

    sage: pari(2).asinh()
    1.44363547517881
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(2+i).asinh()                                                             # needs sage.rings.real_mpfr
    1.52857091948100 + 0.427078586392476*I

    sage: pari(1).atan()
    0.785398163397448
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1.5+i).atan()                                                            # needs sage.rings.real_mpfr
    1.10714871779409 + 0.255412811882995*I

    sage: pari(0).atanh()
    0.E-19
    sage: pari(2).atanh()
    0.549306144334055 - 1.57079632679490*I

    sage: pari(2).besselh1(3)
    0.486091260585891 - 0.160400393484924*I
    sage: pari(2).besselh2(3)
    0.486091260585891 + 0.160400393484924*I
    sage: pari(2).besselj(3)
    0.486091260585891
    sage: pari(2).besseljh(3)
    0.412710032209716
    sage: pari(2).besseli(3)
    2.24521244092995

    sage: # needs sage.rings.real_mpfr
    sage: C.<i> = ComplexField()
    sage: pari(2).besseli(3+i)
    1.12539407613913 + 2.08313822670661*I
    sage: C.<i> = ComplexField()
    sage: pari(2+i).bessely(3)
    -0.280775566958244 - 0.486708533223726*I

    sage: pari(1.5).cos()
    0.0707372016677029
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).cos()                                                               # needs sage.rings.real_mpfr
    0.833730025131149 - 0.988897705762865*I
    sage: pari('x+O(x^8)').cos()
    1 - 1/2*x^2 + 1/24*x^4 - 1/720*x^6 + 1/40320*x^8 + O(x^9)

    sage: pari(1.5).cosh()
    2.35240961524325
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).cosh()                                                              # needs sage.rings.real_mpfr
    0.833730025131149 + 0.988897705762865*I
    sage: pari('x+O(x^8)').cosh()
    1 + 1/2*x^2 + 1/24*x^4 + 1/720*x^6 + ...O(...)

    sage: pari(5).cotan()
    -0.295812915532746
    sage: x = RR(pi)                                                                    # needs sage.symbolic
    sage: pari(x).cotan()  # random
    -8.17674825 E15

    sage: pari(1).dilog()
    1.64493406684823
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).dilog()                                                             # needs sage.rings.real_mpfr
    0.616850275068085 + 1.46036211675312*I

    sage: pari(1).erfc()
    0.157299207050285

    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(i).eta()                                                                 # needs sage.rings.real_mpfr
    0.998129069925959

    sage: pari(0).exp()
    1.00000000000000
    sage: pari(1).exp()
    2.71828182845905
    sage: pari('x+O(x^8)').exp()
    1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + 1/720*x^6 + 1/5040*x^7 + O(x^8)

    sage: pari(2).gamma()
    1.00000000000000
    sage: pari(5).gamma()
    24.0000000000000
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).gamma()                                                             # needs sage.rings.real_mpfr
    0.498015668118356 - 0.154949828301811*I
    sage: pari(-1).gamma()
    Traceback (most recent call last):
    ...
    PariError: domain error in gamma: argument = nonpositive integer

    sage: pari(2).gammah()
    1.32934038817914
    sage: pari(5).gammah()
    52.3427777845535
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).gammah()                                                            # needs sage.rings.real_mpfr
    0.575315188063452 + 0.0882106775440939*I

    sage: pari(1).hyperu(2,3)
    0.333333333333333

    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).incgam(3-i)                                                         # needs sage.rings.real_mpfr
    -0.0458297859919946 + 0.0433696818726677*I
    sage: pari(1).incgamc(2)
    0.864664716763387

    sage: pari(5).log()
    1.60943791243410
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(i).log()                                                                 # needs sage.rings.real_mpfr
    0.E-19 + 1.57079632679490*I

    sage: pari(100).lngamma()
    359.134205369575
    sage: pari(100).log_gamma()
    359.134205369575

    sage: pari(1).psi()
    -0.577215664901533

    sage: pari(1).sin()
    0.841470984807897
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).sin()                                                               # needs sage.rings.real_mpfr
    1.29845758141598 + 0.634963914784736*I

    sage: pari(0).sinh()
    0.E-19
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(1+i).sinh()                                                              # needs sage.rings.real_mpfr
    0.634963914784736 + 1.29845758141598*I

    sage: pari(2).sqrt()
    1.41421356237310

    sage: pari(8).sqrtint()
    2
    sage: pari(10^100).sqrtint()
    100000000000000000000000000000000000000000000000000

    sage: pari(2).tan()
    -2.18503986326152
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: pari(i).tan()                                                                 # needs sage.rings.real_mpfr
    0.761594155955765*I

    sage: pari(1).tanh()
    0.761594155955765
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: z = pari(i); z                                                                # needs sage.rings.real_mpfr
    1.00000000000000*I
    sage: result = z.tanh()
    sage: result.real() <= 1e-18
    True
    sage: result.imag()
    1.55740772465490

    sage: pari('2+O(7^5)').teichmuller()
    2 + 4*7 + 6*7^2 + 3*7^3 + O(7^5)

    sage: pari(0.5).theta(2)
    1.63202590295260

    sage: pari(0.5).thetanullk(1)
    0.548978532560341

    sage: # needs sage.rings.real_mpfr
    sage: C.<i> = ComplexField()
    sage: pari(i).weber()
    1.18920711500272
    sage: pari(i).weber(1)
    1.09050773266526
    sage: pari(i).weber(2)
    1.09050773266526

    sage: pari(2).zeta()
    1.64493406684823
    sage: x = RR(pi)^2/6                                                                # needs sage.symbolic
    sage: pari(x)                                                                       # needs sage.symbolic
    1.64493406684823
    sage: pari(3).zeta()
    1.20205690315959
    sage: pari('1+5*7+2*7^2+O(7^3)').zeta()
    4*7^-2 + 5*7^-1 + O(7^0)

Linear algebra::

    sage: pari('[1,2,3; 4,5,6;  7,8,9]').matadjoint()
    [-3, 6, -3; 6, -12, 6; -3, 6, -3]
    sage: pari('[a,b,c; d,e,f; g,h,i]').matadjoint()
    [(i*e - h*f), (-i*b + h*c), (f*b - e*c); (-i*d + g*f), i*a - g*c, -f*a + d*c; (h*d - g*e), -h*a + g*b, e*a - d*b]

    sage: pari('[1,1;1,-1]').matsolve(pari('[1;0]'))
    [1/2; 1/2]

    sage: D = pari('[3,4]~')
    sage: B = pari('[1,2]~')
    sage: M = pari('[1,2;3,4]')
    sage: M.matsolvemod(D, B)
    [10, 0]~
    sage: M.matsolvemod(3, 1)
    [2, 1]~
    sage: M.matsolvemod(pari('[3,0]~'), pari('[1,2]~'))
    [6, -4]~
    sage: M2 = pari('[1,10;9,18]')
    sage: M2.matsolvemod(3, pari('[2,3]~'), 1)
    [[2, 0]~, [3, 2; 0, 1]]
    sage: M2.matsolvemod(9, pari('[2,3]~'))
    0
    sage: M2.matsolvemod(9, pari('[2,45]~'), 1)
    [[2, 0]~, [9, 8; 0, 1]]

    sage: pari('[1,2,3;4,5,6;7,8,9]').matker()
    [1; -2; 1]
    sage: pari('[1,2,3;4,5,6;7,8,9]').matker(1)
    [1; -2; 1]
    sage: pari('matrix(3,3,i,j,i)').matker()
    [-1, -1; 1, 0; 0, 1]
    sage: pari('[1,2,3;4,5,6;7,8,9]*Mod(1,2)').matker()
    [Mod(1, 2); Mod(0, 2); Mod(1, 2)]

    sage: pari('[1,2; 3,4]').matdet(0)
    -2
    sage: pari('[1,2; 3,4]').matdet(1)
    -2

    sage: pari('[1,2; 3,4]').trace()
    5

    sage: pari('[1,2,3; 4,5,6;  7,8,9]').mathnf()
    [6, 1; 3, 1; 0, 1]

    sage: # needs sage.modules
    sage: M = matrix([[1,2,3],[4,5,6],[7,8,11]])
    sage: d = M.det()
    sage: pari(M).mathnfmod(d)
    [6, 4, 3; 0, 1, 0; 0, 0, 1]
    sage: M = matrix([[1,0,0],[0,2,0],[0,0,6]])
    sage: pari(M).mathnfmod(6)
    [1, 0, 0; 0, 1, 0; 0, 0, 6]
    sage: pari(M).mathnfmod(12)
    [1, 0, 0; 0, 2, 0; 0, 0, 6]

    sage: # needs sage.modules
    sage: M = matrix([[1,0,0],[0,2,0],[0,0,6]])
    sage: pari(M).mathnfmodid(6)
    [1, 0, 0; 0, 2, 0; 0, 0, 6]
    sage: pari(M).mathnfmod(6)
    [1, 0, 0; 0, 1, 0; 0, 0, 6]

    sage: pari('[1,2,3; 4,5,6;  7,8,9]').matsnf()
    [0, 3, 1]

    sage: a = pari('[1,2;3,4]')
    sage: a.matfrobenius()
    [0, 2; 1, 5]
    sage: a.matfrobenius(flag=1)
    [x^2 - 5*x - 2]
    sage: a.matfrobenius(2)
    [[0, 2; 1, 5], [1, -1/3; 0, 1/3]]
    sage: v = a.matfrobenius(2)
    sage: v[0]
    [0, 2; 1, 5]
    sage: v[1]^(-1)*v[0]*v[1]
    [1, 2; 3, 4]
    sage: t = pari('[3, -2, 0, 0; 0, -2, 0, 1; 0, -1, -2, 2; 0, -2, 0, 2]')
    sage: t.matfrobenius()
    [0, 0, 0, -12; 1, 0, 0, -2; 0, 1, 0, 8; 0, 0, 1, 1]
    sage: t.charpoly('x')
    x^4 - x^3 - 8*x^2 + 2*x + 12
    sage: t.matfrobenius(1)
    [x^4 - x^3 - 8*x^2 + 2*x + 12]

Quadratic forms::

    sage: # needs sage.modules
    sage: A = Matrix(3,3,[1,2,3,2,5,5,3,5,11])
    sage: A.is_positive_definite()
    True
    sage: pari(A).qfminim(10, 5).sage()
    [
             [17 14 15 16 13]
             [-4 -3 -3 -3 -2]
    146, 10, [-3 -3 -3 -3 -3]
    ]
    sage: pari(A).qfminim().sage()
    [
          [ 5  2  1]
          [-1 -1  0]
    6, 1, [-1  0  0]
    ]
    sage: pari(A.change_ring(RR)).qfminim(5, m=5, flag=2).sage()
    [
                             [ -5 -10  -2  -7   3]
                             [  1   2   1   2   0]
    10, 5.00000000000000000, [  1   2   0   1  -1]
    ]

    sage: # needs sage.modules
    sage: M = diagonal_matrix([1,1,-1])
    sage: P = M.__pari__().qfparam([0,1,-1]); P
    [0, -2, 0; 1, 0, -1; -1, 0, -1]
    sage: R.<x,y> = QQ[]
    sage: v = P.sage() * vector([x^2, x*y, y^2]); v
    (-2*x*y, x^2 - y^2, -x^2 - y^2)
    sage: v(x=2, y=1)
    (-4, 3, -5)
    sage: v(x=3,y=8)
    (-48, -55, -73)
    sage: 48^2 + 55^2 == 73^2
    True

    sage: # needs sage.modules
    sage: M = diagonal_matrix([1,2,3,4,-5])
    sage: M.__pari__().qfsolve()
    [0, 1, -1, 0, -1]~
    sage: M = diagonal_matrix([4,-9])
    sage: M.__pari__().qfsolve()
    [6, 4]~
    sage: M = diagonal_matrix([1,1,1,1,1])
    sage: M.__pari__().qfsolve()
    -1
    sage: M = diagonal_matrix([1,1,-3])
    sage: M.__pari__().qfsolve()
    3
    sage: M = diagonal_matrix([1,-42])
    sage: M.__pari__().qfsolve()
    -2
    sage: M = diagonal_matrix([1,-1,0,0])
    sage: M.__pari__().qfsolve().sage()
    [0 0]
    [0 0]
    [1 0]
    [0 1]

Number-theoretical functions::

    sage: n = pari.set_real_precision(210)
    sage: w1 = pari('z1=2-sqrt(26); (z1+I)/(z1-I)')
    sage: f = w1.algdep(12); f
    545*x^11 - 297*x^10 - 281*x^9 + 48*x^8 - 168*x^7 + 690*x^6 - 168*x^5 + 48*x^4 - 281*x^3 - 297*x^2 + 545*x
    sage: f(w1).abs() < 1.0e-200
    True
    sage: f.factor()
    [x, 1; x + 1, 2; x^2 + 1, 1; x^2 + x + 1, 1; 545*x^4 - 1932*x^3 + 2790*x^2 - 1932*x + 545, 1]
    sage: pari.set_real_precision(n)
    210

    sage: pari(6).binomial(2)
    15
    sage: pari('x+1').binomial(3)
    1/6*x^3 - 1/6*x
    sage: pari('2+x+O(x^2)').binomial(3)
    1/3*x + O(x^2)

    sage: pari(10).eulerphi()
    4

    sage: # needs sage.symbolic
    sage: x = SR.symbol('x')
    sage: pari(10).gcd(15)
    5
    sage: pari([5, 'y']).gcd()
    1
    sage: pari([x, x^2]).gcd()
    x
    sage: pari(10).lcm(15)
    30
    sage: pari([5, 'y']).lcm()
    5*y
    sage: pari([10, x, x^2]).lcm()
    10*x^2

    sage: pari(20).numbpart()
    627
    sage: pari(100).numbpart()
    190569292

    sage: pari(10).numdiv()
    4

    sage: pari(7).primepi()
    4
    sage: pari(100).primepi()
    25
    sage: pari(1000).primepi()
    168
    sage: pari(100000).primepi()
    9592
    sage: pari(0).primepi()
    0
    sage: pari(-15).primepi()
    0
    sage: pari(500509).primepi()
    41581
    sage: pari(10^7).primepi()
    664579

    sage: pari(4).znprimroot()
    Mod(3, 4)
    sage: pari(10007^3).znprimroot()
    Mod(5, 1002101470343)
    sage: pari(2*109^10).znprimroot()
    Mod(236736367459211723407, 473472734918423446802)

    sage: pari(0).znstar()
    [2, [2], [-1]]
    sage: pari(96).znstar()
    [32, [8, 2, 2], [Mod(37, 96), Mod(31, 96), Mod(65, 96)]]
    sage: pari(-5).znstar()
    [4, [4], [Mod(2, 5)]]

Finite fields::

    sage: x = GF(2)['x'].gen()
    sage: pari(x^2+x+2).ffgen()
    x
    sage: pari(x^2+x+1).ffgen('a')
    a

    sage: pari(7).ffinit(11)
    Mod(1, 7)*x^11 + Mod(1, 7)*x^10 + Mod(4, 7)*x^9 + Mod(5, 7)*x^8 + Mod(1, 7)*x^7 + Mod(1, 7)*x^2 + Mod(1, 7)*x + Mod(6, 7)
    sage: pari(2003).ffinit(3)
    Mod(1, 2003)*x^3 + Mod(1, 2003)*x^2 + Mod(1993, 2003)*x + Mod(1995, 2003)

    sage: # needs sage.rings.finite_rings
    sage: k.<a> = GF(2^12)
    sage: g = pari(a).ffprimroot()
    sage: (g^1234).fflog(g)
    1234
    sage: pari(k(1)).fflog(g)
    0
    sage: b = g^5
    sage: ord = b.fforder(); ord
    819
    sage: (b^555).fflog(b, ord)
    555
    sage: (b^555).fflog(b, (ord, ord.factor()) )
    555

    sage: # needs sage.rings.finite_rings
    sage: k.<a> = GF(5^80)
    sage: g = pari(a).ffprimroot()
    sage: g.fforder()
    82718061255302767487140869206996285356581211090087890624
    sage: g.fforder( (5^80-1, factor(5^80-1)) )
    82718061255302767487140869206996285356581211090087890624
    sage: k(2).__pari__().fforder(o=4)
    4

`p`-adic functions::

    sage: # needs sage.rings.padics
    sage: K = Qp(11,5)
    sage: x = K(11^-10 + 5*11^-7 + 11^-6)
    sage: y = pari(x)
    sage: y.padicprec(11)
    -5
    sage: y.padicprec(17)
    Traceback (most recent call last):
    ...
    PariError: inconsistent moduli in padicprec: 11 != 17
    sage: R.<t> = PolynomialRing(Zp(3))
    sage: pol = R([O(3^4), O(3^6), O(3^5)])
    sage: pari(pol).padicprec(3)
    4

Elliptic curves::

    sage: e = pari([0,1,0,1,0]).ellinit(); e
    [0, 1, 0, 1, 0, 4, 2, 0, -1, -32, 224, -48, 2048/3, Vecsmall([1]), [Vecsmall([64, -1])], [0, 0, 0, 0, 0, 0, 0, 0]]

    sage: pari([0,1/2,0,-3/4,0]).ellinit()
    [0, 1/2, 0, -3/4, 0, 2, -3/2, 0, -9/16, 40, -116, 117/4, 256000/117, Vecsmall([1]), [Vecsmall([64, 1])], [0, 0, 0, 0, 0, 0, 0, 0]]
    sage: pari([0,0.5,0,-0.75,0]).ellinit()
    [0, 0.500000000000000, 0, -0.750000000000000, 0, 2.00000000000000, -1.50000000000000, 0, -0.562500000000000, 40.0000000000000, -116.000000000000, 29.2500000000000, 2188.03418803419, Vecsmall([0]), [Vecsmall([64, 1])], [0, 0, 0, 0]]
    sage: pari([0,SR(I),0,1,0]).ellinit()                                               # needs sage.symbolic
    [0, I, 0, 1, 0, 4*I, 2, 0, -1, -64, 352*I, -80, 16384/5, Vecsmall([0]), [Vecsmall([64, 0])], [0, 0, 0, 0]]
    sage: x = SR.symbol('x')                                                            # needs sage.symbolic
    sage: pari([0,x,0,2*x,1]).ellinit()                                                 # needs sage.symbolic
    [0, x, 0, 2*x, 1, 4*x, 4*x, 4, -4*x^2 + 4*x, 16*x^2 - 96*x, -64*x^3 + 576*x^2 - 864, 64*x^4 - 576*x^3 + 576*x^2 - 432, (256*x^6 - 4608*x^5 + 27648*x^4 - 55296*x^3)/(4*x^4 - 36*x^3 + 36*x^2 - 27), Vecsmall([0]), [Vecsmall([64, 0])], [0, 0, 0, 0]]

    sage: e = pari([0,1,1,-2,0]).ellinit()
    sage: e.ellheight([1,0])
    0.476711659343740
    sage: e.ellheight([1,0], precision=128).sage()
    0.476711659343739537379486058884653059459022942211150879336
    sage: e.ellheight([1, 0], [-1, 1])
    0.418188984498861

    sage: e = pari([0,1,1,-2,0]).ellinit()
    sage: x = pari([1,0])
    sage: e.ellisoncurve([1,4])
    False
    sage: e.ellisoncurve(x)
    True
    sage: f = e.ellchangecurve([1,2,3,-1])
    sage: f[:5]   # show only first five entries
    [6, -2, -1, 17, 8]
    sage: x.ellchangepoint([1,2,3,-1])
    [-1, 4]
    sage: f.ellisoncurve([-1,4])
    True

    sage: e = pari([0, 5, 2, -1, 1]).ellinit()
    sage: e.ellglobalred()
    [20144, [1, -2, 0, -1], 1, [2, 4; 1259, 1], [[4, 2, 0, 1], [1, 5, 0, 1]]]

    sage: # needs sage.schemes
    sage: e = pari(EllipticCurve('17a').a_invariants()).ellinit()
    sage: e.ellglobalred()
    [17, [1, 0, 0, 0], 4, Mat([17, 1]), [[1, 8, 0, 4]]]

    sage: e = pari([0, 1, 1, -2, 0]).ellinit()
    sage: e.elladd([1,0], [-1,1])
    [-3/4, -15/8]

    sage: e = pari([0, -1, 1, -10, -20]).ellinit()
    sage: e.ellak(6)
    2
    sage: e.ellak(2005)
    2
    sage: e.ellak(-1)
    0
    sage: e.ellak(0)
    0

    sage: # needs sage.schemes
    sage: E = EllipticCurve('389a1')
    sage: pari(E).ellanalyticrank()                                                     # needs sage.rings.number_field
    [2, 1.51863300057685]

    sage: e = pari([0, -1, 1, -10, -20]).ellinit()
    sage: e.ellap(2)
    -2
    sage: e.ellap(2003)
    4

    sage: e = pari([1,2,3,4,5]).ellinit()
    sage: e.ellglobalred()
    [10351, [1, -1, 0, -1], 1, [11, 1; 941, 1], [[1, 5, 0, 1], [1, 5, 0, 1]]]
    sage: f = e.ellchangecurve([1,-1,0,-1])
    sage: f[:5]
    [1, -1, 0, 4, 3]

    sage: e = pari([0,0,0,-82,0]).ellinit()
    sage: e.elleta()
    [3.60546360143265, 3.60546360143265*I]
    sage: w1, w2 = e.omega()
    sage: eta1, eta2 = e.elleta()
    sage: w1*eta2 - w2*eta1
    6.28318530717959*I

    sage: e = pari([0,1,1,-2,0]).ellinit().ellminimalmodel()[0]
    sage: e.ellheightmatrix([[1,0], [-1,1]])
    [0.476711659343740, 0.418188984498861; 0.418188984498861, 0.686667083305587]

    sage: e = pari([0,1,1,-2,0]).ellinit()
    sage: om = e.omega()
    sage: om
    [2.49021256085506, -1.97173770155165*I]
    sage: om.elleisnum(2)
    10.0672605281120
    sage: om.elleisnum(4)
    112.000000000000
    sage: om.elleisnum(100)
    2.15314248576078 E50

    sage: # needs sage.schemes
    sage: e = pari([0,0,0,0,1]).ellinit()
    sage: e.elllocalred(7)
    [0, 1, [1, 0, 0, 0], 1]
    sage: e = pari(EllipticCurve('27a3').a_invariants()).ellinit()
    sage: e.elllocalred(3)
    [3, 2, [1, 0, 0, 0], 1]
    sage: e = pari(EllipticCurve('24a4').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [3, 3, [1, 0, 0, 0], 2]
    sage: e = pari(EllipticCurve('20a2').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [2, 4, [1, 0, 0, 0], 3]
    sage: e = pari(EllipticCurve('11a2').a_invariants()).ellinit()
    sage: e.elllocalred(11)
    [1, 5, [1, 0, 0, 0], 1]
    sage: e = pari(EllipticCurve('14a4').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [1, 6, [1, 0, 0, 0], 2]
    sage: e = pari(EllipticCurve('14a1').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [1, 10, [1, 0, 0, 0], 2]
    sage: e = pari(EllipticCurve('32a3').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [5, -1, [1, 0, 0, 0], 1]
    sage: e = pari(EllipticCurve('24a5').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [3, -2, [1, 0, 0, 0], 1]
    sage: e = pari(EllipticCurve('24a2').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [3, -3, [1, 0, 0, 0], 2]
    sage: e = pari(EllipticCurve('20a1').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [2, -4, [1, 0, 0, 0], 3]
    sage: e = pari(EllipticCurve('24a1').a_invariants()).ellinit()
    sage: e.elllocalred(2)
    [3, -5, [1, 0, 0, 0], 4]
    sage: e = pari(EllipticCurve('90c2').a_invariants()).ellinit()
    sage: e.elllocalred(3)
    [2, -10, [1, 0, 0, 0], 4]

    sage: # needs sage.schemes
    sage: e = pari(EllipticCurve('65a1').a_invariants()).ellinit()
    sage: e.ellorder([0,0])
    2
    sage: e.ellorder([1,0])
    0

    sage: e = pari([0,1,1,-2,0]).ellinit()
    sage: e.ellordinate(0)
    [0, -1]
    sage: e.ellordinate(SR(I))                                                          # needs sage.symbolic
    [0.582203589721741 - 1.38606082464177*I, -1.58220358972174 + 1.38606082464177*I]
    sage: e.ellordinate(SR(I), precision=128)[0].sage()                                 # needs sage.symbolic
    0.58220358972174117723338947874993600727 - 1.3860608246417697185311834209833653345*I
    sage: e.ellordinate(1+3*5^1+O(5^3))                                                 # needs sage.rings.padics
    [4*5 + 5^2 + O(5^3), 4 + 3*5^2 + O(5^3)]
    sage: e.ellordinate('z+2*z^2+O(z^4)')
    [-2*z - 7*z^2 - 23*z^3 + O(z^4), -1 + 2*z + 7*z^2 + 23*z^3 + O(z^4)]
    sage: e.ellordinate(5)
    []
    sage: e.ellordinate(5.0)
    [11.3427192823270, -12.3427192823270]

    sage: e = pari([0,0,0,1,0]).ellinit()
    sage: e.ellpointtoz([0,0])
    1.85407467730137
    sage: e.ellpointtoz([0])
    0

    sage: e = pari([0,0,0,3,0]).ellinit()
    sage: p = [1,2]  # Point of infinite order
    sage: e.ellmul([0,0], 2)
    [0]
    sage: e.ellmul(p, 2)
    [1/4, -7/8]
    sage: q = e.ellmul(p, SR(1+I)); q                                                   # needs sage.symbolic
    [-2*I, 1 + I]
    sage: e.ellmul(q, SR(1-I))                                                          # needs sage.symbolic
    [1/4, -7/8]
    sage: for D in [-7, -8, -11, -12, -16, -19, -27, -28]:  # long time (1s)
    ....:     hcpol = hilbert_class_polynomial(D)
    ....:     j = hcpol.roots(multiplicities=False)[0]
    ....:     t = (1728-j)/(27*j)
    ....:     E = EllipticCurve([4*t,16*t^2])
    ....:     P = E.point([0, 4*t])
    ....:     assert(E.j_invariant() == j)
    ....:     #
    ....:     # Compute some CM number and its minimal polynomial
    ....:     #
    ....:     cm = pari('cm = (3*quadgen(%s)+2)'%D)
    ....:     cm_minpoly = pari('minpoly(cm)')
    ....:     #
    ....:     # Evaluate cm_minpoly(cm)(P), which should be zero
    ....:     #
    ....:     e = pari(E)  # Convert E to PARI
    ....:     P2 = e.ellmul(P, cm_minpoly[2]*cm + cm_minpoly[1])
    ....:     P0 = e.elladd(e.ellmul(P, cm_minpoly[0]), e.ellmul(P2, cm))
    ....:     assert(P0 == E(0))

    sage: e = pari([0,0,0,-82,0]).ellinit()
    sage: e.ellrootno()
    -1
    sage: e.ellrootno(2)
    1
    sage: e.ellrootno(1009)
    1

    sage: # needs sage.rings.real_mpfr
    sage: e = pari([0,0,0,1,0]).ellinit()
    sage: C.<i> = ComplexField()
    sage: e.ellsigma(2+i)
    1.43490215804166 + 1.80307856719256*I

    sage: e = pari([0, 1, 1, -2, 0]).ellinit()
    sage: e.ellsub([1,0], [-1,1])
    [0, 0]

    sage: e = pari([0,0,0,1,0]).ellinit()
    sage: e.ellzeta(1)
    1.06479841295883
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: e.ellzeta(i-1)                                                                # needs sage.rings.real_mpfr
    -0.350122658523049 - 0.350122658523049*I

    sage: e = pari([0,0,0,1,0]).ellinit()
    sage: C.<i> = ComplexField()                                                        # needs sage.rings.real_mpfr
    sage: e.ellztopoint(1+i)                                                            # needs sage.rings.real_mpfr
    [0.E-... - 1.02152286795670*I, -0.149072813701096 - 0.149072813701096*I]
    sage: e.ellztopoint(0)
    [0]

    sage: pari(SR(I)).ellj()                                                            # needs sage.symbolic
    1728.00000000000
    sage: pari(SR(3*I)).ellj()                                                          # needs sage.symbolic
    153553679.396729
    sage: pari('quadgen(-3)').ellj()
    0.E-54
    sage: pari('quadgen(-7)').ellj(precision=256).sage()
    -3375.000000000000000000000000000000000000000000000000000000000000000000000000
    sage: pari(SR(-I)).ellj()                                                           # needs sage.symbolic
    Traceback (most recent call last):
    ...
    PariError: domain error in modular function: Im(argument) <= 0

Quadratic class numbers::

    sage: pari(10009).qfbhclassno()
    0
    sage: pari(2).qfbhclassno()
    0
    sage: pari(0).qfbhclassno()
    -1/12
    sage: pari(4).qfbhclassno()
    1/2
    sage: pari(3).qfbhclassno()
    1/3
    sage: pari(23).qfbhclassno()
    3

    sage: pari(-4).qfbclassno()
    1
    sage: pari(-23).qfbclassno()
    3
    sage: pari(-104).qfbclassno()
    6
    sage: pari(109).qfbclassno()
    1
    sage: pari(10001).qfbclassno()
    16
    sage: pari(10001).qfbclassno(flag=1)
    16
    sage: pari(3).qfbclassno()
    Traceback (most recent call last):
    ...
    PariError: domain error in classno2: disc % 4 > 1
    sage: pari(4).qfbclassno()
    Traceback (most recent call last):
    ...
    PariError: domain error in classno2: issquare(disc) = 1

    sage: pari(-4).quadclassunit()
    [1, [], [], 1]
    sage: pari(-23).quadclassunit()
    [3, [3], [Qfb(2, 1, 3)], 1]
    sage: pari(-104).quadclassunit()
    [6, [6], [Qfb(5, -4, 6)], 1]
    sage: pari(109).quadclassunit()
    [1, [], [], 5.56453508676047, -1]
    sage: pari(10001).quadclassunit() # random generators
    [16, [16], [Qfb(10, 99, -5, 0.E-38)], 5.29834236561059]
    sage: pari(10001).quadclassunit()[0]
    16
    sage: pari(10001).quadclassunit()[1]
    [16]
    sage: pari(10001).quadclassunit()[3]
    5.29834236561059
    sage: pari(3).quadclassunit()
    Traceback (most recent call last):
    ...
    PariError: domain error in Buchquad: disc % 4 > 1
    sage: pari(4).quadclassunit()
    Traceback (most recent call last):
    ...
    PariError: domain error in Buchquad: issquare(disc) = 1

General number fields::

    sage: x = polygen(QQ)
    sage: K.<a> = NumberField(x^2 - 1/8)                                                # needs sage.rings.number_field
    sage: pari(x^2 - 2).factornf(K.pari_polynomial("a"))                                # needs sage.rings.number_field
    doctest:...: DeprecationWarning: the PARI/GP function factornf is obsolete (2016-08-08)
    [x + Mod(-a, a^2 - 2), 1; x + Mod(a, a^2 - 2), 1]

    sage: K.<z> = QuadraticField(-23)                                                   # needs sage.rings.number_field
    sage: p = K.primes_above(3)[0]                                                      # needs sage.rings.number_field
    sage: K.pari_bnf().bnrclassno(p._pari_bid_())                                       # needs sage.rings.number_field
    3

    sage: # needs sage.symbolic
    sage: x = SR.symbol('x')
    sage: P = pari(x^6 + 108)
    sage: G = P.galoisinit()
    sage: G[0] == P
    True
    sage: len(G[5]) == prod(G[7])
    True

    sage: # needs sage.symbolic
    sage: G = pari(x^6 + 108).galoisinit()
    sage: G.galoispermtopol(G[5])
    [x, 1/12*x^4 - 1/2*x, -1/12*x^4 - 1/2*x, 1/12*x^4 + 1/2*x, -1/12*x^4 + 1/2*x, -x]
    sage: G.galoispermtopol(G[5][1])
    1/12*x^4 - 1/2*x
    sage: G.galoispermtopol(G[5][1:4])
    [1/12*x^4 - 1/2*x, -1/12*x^4 - 1/2*x, 1/12*x^4 + 1/2*x]

    sage: # needs sage.symbolic
    sage: G = pari(x^4 + 1).galoisinit()
    sage: G.galoisfixedfield(G[5][1], flag=2)
    [y^2 - 2, Mod(-x^3 + x, x^4 + 1), [x^2 - y*x + 1, x^2 + y*x + 1]]
    sage: G.galoisfixedfield(G[5][5:7])
    [x^4 + 1, Mod(x, x^4 + 1)]
    sage: L = G.galoissubgroups()
    sage: G.galoisfixedfield(L[3], flag=2, v='z')
    [z^2 + 2, Mod(x^3 + x, x^4 + 1), [x^2 - z*x - 1, x^2 + z*x - 1]]

    sage: G = pari(x^6 + 108).galoisinit()
    sage: L = G.galoissubgroups()
    sage: list(L[0][1])
    [3, 2]

    sage: # needs sage.symbolic
    sage: G = pari(x^6 + 108).galoisinit()
    sage: G.galoisisabelian()
    0
    sage: H = G.galoissubgroups()[2]
    sage: H.galoisisabelian()
    Mat(2)
    sage: H.galoisisabelian(flag=1)
    1

    sage: # needs sage.symbolic
    sage: G = pari(x^6 + 108).galoisinit()
    sage: L = G.galoissubgroups()
    sage: G.galoisisnormal(L[0])
    1
    sage: G.galoisisnormal(L[2])
    0

    sage: # needs sage.rings.number_field
    sage: F = QuadraticField(5, 'alpha')
    sage: nf = F.__pari__()
    sage: P = F.ideal(F.gen())
    sage: Q = F.ideal(2)
    sage: moduli = pari.matrix(2,2,[P.pari_prime(),4,Q.pari_prime(),4])                 # needs sage.symbolic
    sage: residues = pari.vector(2,[0,1])
    sage: b = F(nf.idealchinese(moduli,residues))                                       # needs sage.symbolic
    sage: b.valuation(P) >= 4                                                           # needs sage.symbolic
    True
    sage: (b-1).valuation(Q) >= 2
    True

    sage: # needs sage.symbolic
    sage: F = NumberField(x^3-2, 'alpha')
    sage: nf = F.__pari__()
    sage: x = pari('[1, -1, 2]~')
    sage: y = pari('[1, -1, 3]~')
    sage: nf.idealcoprime(x, y)
    1

    sage: y = pari('[2, -2, 4]~')
    sage: nf.idealcoprime(x, y)                                                         # needs sage.rings.number_field
    [5/43, 9/43, -1/43]~

    sage: # needs sage.rings.number_field
    sage: R.<x> = PolynomialRing(QQ)
    sage: K.<a> = NumberField(x^2 + 1)
    sage: L = K.pari_nf().ideallist(100)
    sage: L[0]   # One ideal of norm 1.                                                 # needs sage.symbolic
    [[1, 0; 0, 1]]
    sage: L[64]  # 4 ideals of norm 65.                                                 # needs sage.symbolic
    [[65, 8; 0, 1], [65, 47; 0, 1], [65, 18; 0, 1], [65, 57; 0, 1]]

    sage: # needs sage.rings.number_field
    sage: F = NumberField(x^3-2, 'alpha')
    sage: nf = F.__pari__()
    sage: I = pari('[1, -1, 2]~')
    sage: bid = nf.idealstar(I)
    sage: nf.ideallog(5, bid)
    [25]~

    sage: # needs sage.rings.number_field
    sage: K.<i> = QuadraticField(-1)
    sage: F = pari(K).idealprimedec(5); F
    [[5, [-2, 1]~, 1, 1, [2, -1; 1, 2]], [5, [2, 1]~, 1, 1, [-2, -1; 1, -2]]]
    sage: F[0].pr_get_p()
    5

    sage: # needs sage.rings.number_field
    sage: x = polygen(ZZ)
    sage: F = NumberField(x^3 - 2, 'alpha')
    sage: nf = F.__pari__()
    sage: I = pari('[1, -1, 2]~')
    sage: nf.idealstar(I)
    [[[43, 9, 5; 0, 1, 0; 0, 0, 1], [0]], [42, [42]], [Mat([[43, [9, 1, 0]~, 1, 1, [-5, 2, -18; -9, -5, 2; 1, -9, -5]], 1]), Mat([[43, [9, 1, 0]~, 1, 1, [-5, 2, -18; -9, -5, 2; 1, -9, -5]], 1])], [[[[42], [3], [43, 9, 5; 0, 1, 0; 0, 0, 1], [[[-14, -8, 20]~, [1, 34, 38], [43, [9, 1, 0]~, 1, 1, [-5, 2, -18; -9, -5, 2; 1, -9, -5]]]~, 3, [42, [2, 1; 3, 1; 7, 1]]]]], [[], Vecsmall([])]], [Mat(1)]]

    sage: # needs sage.rings.number_field
    sage: x = polygen(QQ)
    sage: K.<a> = NumberField(x^3 - 17)
    sage: Kpari = K.pari_nf()
    sage: Kpari.getattr('zk')
    [1, 1/3*y^2 - 1/3*y + 1/3, y]
    sage: Kpari.nfbasistoalg(42)
    Mod(42, y^3 - 17)
    sage: Kpari.nfbasistoalg("[3/2, -5, 0]~")
    Mod(-5/3*y^2 + 5/3*y - 1/6, y^3 - 17)
    sage: Kpari.getattr('zk') * pari("[3/2, -5, 0]~")
    -5/3*y^2 + 5/3*y - 1/6

    sage: # needs sage.rings.number_field
    sage: k.<a> = NumberField(x^2 + 5)
    sage: x = 10
    sage: y = a + 1
    sage: pari(k).nfeltdiveuc(x, y)
    [2, -2]~

    sage: # needs sage.rings.number_field
    sage: x = polygen(ZZ)
    sage: k.<a> = NumberField(x^2 + 5)
    sage: I = k.ideal(a)
    sage: kp = pari(k)
    sage: kp.nfeltreduce(12, I.pari_hnf())
    [2, 0]~
    sage: 12 - k(kp.nfeltreduce(12, I.pari_hnf())) in I
    True

    sage: x = QQ['x'].0; nf = pari(x^2 + 2).nfinit()
    sage: nf.nfgaloisconj()
    [-x, x]~
    sage: nf = pari(x^3 + 2).nfinit()
    sage: nf.nfgaloisconj()
    [x]~
    sage: nf = pari(x^4 + 2).nfinit()
    sage: nf.nfgaloisconj()
    [-x, x]~

    sage: # needs sage.rings.number_field
    sage: x = polygen(QQ)
    sage: K.<t> = NumberField(x^3 - x + 1)
    sage: pari(K).nfhilbert(t, t + 2)
    -1
    sage: P = K.ideal(t^2 + t - 2)   # Prime ideal above 5
    sage: pari(K).nfhilbert(t, t + 2, P.pari_prime())
    -1
    sage: P = K.ideal(t^2 + 3*t - 1) # Prime ideal above 23, ramified
    sage: pari(K).nfhilbert(t, t + 2, P.pari_prime())
    1

    sage: # needs sage.rings.number_field
    sage: F.<a> = NumberField(x^2-x-1)
    sage: Fp = pari(F)
    sage: A = matrix(F,[[1,2,a,3],[3,0,a+2,0],[0,0,a,2],[3+a,a,0,1]])
    sage: I = [F.ideal(-2*a+1),F.ideal(7), F.ideal(3),F.ideal(1)]
    sage: Fp.nfhnf([pari(A),[pari(P) for P in I]])
    [[1, [-969/5, -1/15]~, [15, -2]~, [-1938, -3]~; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1], [[3997, 1911; 0, 7], [15, 6; 0, 3], 1, 1]]
    sage: K.<b> = NumberField(x^3-2)
    sage: Kp = pari(K)
    sage: A = matrix(K,[[1,0,0,5*b],[1,2*b^2,b,57],[0,2,1,b^2-3],[2,0,0,b]])
    sage: I = [K.ideal(2),K.ideal(3+b^2),K.ideal(1),K.ideal(1)]
    sage: Kp.nfhnf([pari(A),[pari(P) for P in I]])
    [[1, -225, 72, -31; 0, 1, [0, -1, 0]~, [0, 0, -1/2]~; 0, 0, 1, [0, 0, -1/2]~; 0, 0, 0, 1], [[1116, 756, 612; 0, 18, 0; 0, 0, 18], 2, 1, [2, 0, 0; 0, 1, 0; 0, 0, 1]]]
    sage: K.<b> = NumberField(x^2+5)
    sage: Kp = pari(K)
    sage: A = matrix(K,[[1,0,0,5*b],[1,2*b^2,b,57],[0,2,1,b^2-3],[2,0,0,b]])
    sage: I = [K.ideal(2),K.ideal(3+b^2),K.ideal(1),K.ideal(1)]
    sage: Kp.nfhnf([pari(A),[pari(P) for P in I]])
    [[1, [15, 6]~, [0, -54]~, [113, 72]~; 0, 1, [-4, -1]~, [0, -1]~; 0, 0, 1, 0; 0, 0, 0, 1], [[360, 180; 0, 180], [6, 4; 0, 2], 1, 1]]
    sage: A = matrix(K,[[1,0,0,5*b],[1,2*b,b,57],[0,2,1,b-3],[2,0,b,b]])
    sage: I = [K.ideal(2).factor()[0][0],K.ideal(3+b),K.ideal(1),K.ideal(1)]
    sage: Kp.nfhnf([pari(A),[pari(P) for P in I]])
    [[1, [7605, 4]~, [5610, 5]~, [7913, -6]~; 0, 1, 0, -1; 0, 0, 1, 0; 0, 0, 0, 1], [[19320, 13720; 0, 56], [2, 1; 0, 1], 1, 1]]

    sage: pari('x^3 - 17').nfinit()
    [x^3 - 17, [1, 1], -867, 3, [[1, 1.68006914259990, 2.57128159065824; 1, -0.340034571299952 - 2.65083754153991*I, -1.28564079532912 + 2.22679517779329*I], [1, 1.68006914259990, 2.57128159065824; 1, -2.99087211283986, 0.941154382464174; 1, 2.31080297023995, -3.51243597312241], [16, 27, 41; 16, -48, 15; 16, 37, -56], [3, 1, 0; 1, -11, 17; 0, 17, 0], [51, 0, 16; 0, 17, 3; 0, 0, 1], [17, 0, -1; 0, 0, 3; -1, 3, 2], [51, [-17, 6, -1; 0, -18, 3; 1, 0, -16]], [3, 17]], [2.57128159065824, -1.28564079532912 + 2.22679517779329*I], [3, x^2 - x + 1, 3*x], [1, 0, -1; 0, 0, 3; 0, 1, 1], [1, 0, 0, 0, -4, 6, 0, 6, -1; 0, 1, 0, 1, 1, -1, 0, -1, 3; 0, 0, 1, 0, 2, 0, 1, 0, 1]]
    sage: pari('x^2 + 10^100 + 1').nfinit()
    [...]
    sage: pari('1.0').nfinit()
    Traceback (most recent call last):
    ...
    PariError: incorrect type in checknf [please apply nfinit()] (t_REAL)

    sage: # needs sage.rings.number_field
    sage: F = NumberField(x^3-2,'alpha')
    sage: G = NumberField(x^3-2,'beta')
    sage: F.__pari__().nfisisom(G.__pari__())
    [y]
    sage: GG = NumberField(x^3-4,'gamma')
    sage: F.__pari__().nfisisom(GG.__pari__())
    [1/2*y^2]
    sage: F.__pari__().nfisisom(GG.pari_nf())
    [1/2*y^2]
    sage: F.pari_nf().nfisisom(GG.__pari__()[0])
    [1/2*y^2]
    sage: H = NumberField(x^2-2,'alpha')
    sage: F.__pari__().nfisisom(H.__pari__())
    0
    sage: K.<a> = NumberField(x^2 + x + 1)
    sage: L.<b> = NumberField(x^2 + 3)
    sage: pari(K).nfisisom(L)
    [-1/2*y - 1/2, 1/2*y - 1/2]

    sage: y = QQ['yy'].0; _ = pari(y) # pari has variable ordering rules
    sage: x = QQ['zz'].0; nf = pari(x^2 + 2).nfinit()
    sage: nf.nfroots(y^2 + 2)
    [Mod(-zz, zz^2 + 2), Mod(zz, zz^2 + 2)]~
    sage: nf = pari(x^3 + 2).nfinit()
    sage: nf.nfroots(y^3 + 2)
    [Mod(zz, zz^3 + 2)]~
    sage: nf = pari(x^4 + 2).nfinit()
    sage: nf.nfroots(y^4 + 2)
    [Mod(-zz, zz^4 + 2), Mod(zz, zz^4 + 2)]~

    sage: nf = pari('x^2 + 1').nfinit()
    sage: nf.nfrootsof1()
    [4, [0, 1]~]

    sage: x = ZZ['xx1'].0; pari(x)
    xx1
    sage: y = ZZ['yy1'].0; pari(y)
    yy1
    sage: nf = pari(y^2 - 6*y + 24).nfinit()
    sage: rnf = nf.rnfinit(x^2 - pari(y))
    sage: P = pari('[[[1, 0]~, [0, 0]~; [0, 0]~, [1, 0]~], [[2, 0; 0, 2], [2, 0; 0, 1/2]]]')
    sage: rnf.rnfidealdown(P)
    2

    sage: f = pari('y^3+y+1')
    sage: K = f.nfinit()
    sage: x = pari('x'); y = pari('y')
    sage: g = x^5 - x^2 + y
    sage: L = K.rnfinit(g)

    sage: pari(-23).quadhilbert()
    x^3 - x^2 + 1
    sage: pari(145).quadhilbert()
    x^4 - x^3 - 5*x^2 - x + 1
    sage: pari(-12).quadhilbert()   # Not fundamental
    Traceback (most recent call last):
    ...
    PariError: domain error in quadray: isfundamental(D) = 0

    sage: # needs sage.rings.number_field sage.symbolic
    sage: x = SR.symbol('x')
    sage: F = NumberField(x^3-2,'alpha')
    sage: F.__pari__()[0].nfdisc()
    -108
    sage: G = NumberField(x^5-11,'beta')
    sage: G.__pari__()[0].nfdisc()
    45753125
    sage: f = x^3-2
    sage: f.__pari__()
    x^3 - 2
    sage: f.__pari__().nfdisc()
    -108

These are some doctests that used to be part of Sage and were removed from the cypari2
library::

    sage: e = pari([0,0,0,-82,0]).ellinit()
    sage: eta1 = e.elleta(precision=50)[0]
    sage: eta1.sage()
    3.6054636014326520859158205642077267748
    sage: eta1 = e.elleta(precision=150)[0]
    sage: eta1.sage()
    3.605463601432652085915820564207726774810268996598024745444380641429820491740 # 64-bit
    3.605463601432652085915820564207726774810268996598024745444380641430          # 32-bit
    sage: from cypari2 import Pari
    sage: pari = Pari()

    sage: # needs sage.symbolic
    sage: f = pari('(2/3)*x^3 + x - 5/7 + y'); f
    2/3*x^3 + x + (y - 5/7)
    sage: var('x,y')
    (x, y)
    sage: f.sage({'x':x, 'y':y})
    2/3*x^3 + x + y - 5/7

    sage: pari.default("debug")
    0
    sage: pari.default("debug", 3)
    sage: pari(2**67+1).factor()
    IFAC: cracking composite
            49191317529892137643
    IFAC: factor 6713103182899
            is prime
    IFAC: factor 7327657
            is prime
    IFAC: prime 7327657
            appears with exponent = 1
    IFAC: prime 6713103182899
            appears with exponent = 1
    IFAC: found 2 large prime (power) factors.
    [3, 1; 7327657, 1; 6713103182899, 1]
    sage: pari.default("debug", 0)
    sage: pari(2**67+1).factor()
    [3, 1; 7327657, 1; 6713103182899, 1]

    sage: pari(18).bernreal(precision=192).sage()
    54.9711779448621553884711779448621553884711779448621553885
"""
