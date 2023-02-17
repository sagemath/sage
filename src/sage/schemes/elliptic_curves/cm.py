"""
Complex multiplication for elliptic curves

This module implements the functions

- ``hilbert_class_polynomial``
- ``is_HCP``
- ``cm_j_invariants``
- ``cm_orders``
- ``discriminants_with_bounded_class_number``
- ``cm_j_invariants_and_orders``
- ``largest_fundamental_disc_with_class_number``
- ``is_cm_j_invariant``

AUTHORS:

- Robert Bradshaw
- John Cremona
- William Stein

"""

# ****************************************************************************
#       Copyright (C) 2005-2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.interfaces.magma import magma
from sage.rings.all import (Integer,
                            QQ,
                            ZZ,
                            IntegerRing,
                            GF,
                            RR)

from sage.schemes.elliptic_curves.all import EllipticCurve
from sage.misc.cachefunc import cached_function

@cached_function
def hilbert_class_polynomial(D, algorithm=None):
    r"""
    Return the Hilbert class polynomial for discriminant `D`.

    INPUT:

    - ``D`` (int) -- a negative integer congruent to 0 or 1 modulo 4.

    - ``algorithm`` (string, default None).

    OUTPUT:

    (integer polynomial) The Hilbert class polynomial for the
    discriminant `D`.

    ALGORITHM:

    - If ``algorithm`` = "arb" (default): Use Arb's implementation which uses complex interval arithmetic.

    - If ``algorithm`` = "sage": Use complex approximations to the roots.

    - If ``algorithm`` = "magma": Call the appropriate Magma function (if available).

    AUTHORS:

    - Sage implementation originally by Eduardo Ocampo Alvarez and
      AndreyTimofeev

    - Sage implementation corrected by John Cremona (using corrected precision bounds from Andreas Enge)

    - Magma implementation by David Kohel

    EXAMPLES::

        sage: hilbert_class_polynomial(-4)
        x - 1728
        sage: hilbert_class_polynomial(-7)
        x + 3375
        sage: hilbert_class_polynomial(-23)
        x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
        sage: hilbert_class_polynomial(-37*4)
        x^2 - 39660183801072000*x - 7898242515936467904000000
        sage: hilbert_class_polynomial(-37*4, algorithm="magma") # optional - magma
        x^2 - 39660183801072000*x - 7898242515936467904000000
        sage: hilbert_class_polynomial(-163)
        x + 262537412640768000
        sage: hilbert_class_polynomial(-163, algorithm="sage")
        x + 262537412640768000
        sage: hilbert_class_polynomial(-163, algorithm="magma") # optional - magma
        x + 262537412640768000

    TESTS::

        sage: all(hilbert_class_polynomial(d, algorithm="arb") ==
        ....:      hilbert_class_polynomial(d, algorithm="sage")
        ....:        for d in range(-1,-100,-1) if d % 4 in [0, 1])
        True
    """
    if algorithm is None:
        algorithm = "arb"

    D = Integer(D)
    if D >= 0:
        raise ValueError("D (=%s) must be negative" % D)
    if not (D % 4 in [0, 1]):
        raise ValueError("D (=%s) must be a discriminant" % D)

    if algorithm == "arb":
        import sage.libs.arb.arith
        return sage.libs.arb.arith.hilbert_class_polynomial(D)

    if algorithm == "magma":
        magma.eval("R<x> := PolynomialRing(IntegerRing())")
        f = str(magma.eval("HilbertClassPolynomial(%s)" % D))
        return IntegerRing()['x'](f)

    if algorithm != "sage":
        raise ValueError("%s is not a valid algorithm" % algorithm)

    from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
    from sage.rings.all import RR, ComplexField
    from sage.functions.all import elliptic_j

    # get all primitive reduced quadratic forms, (necessary to exclude
    # imprimitive forms when D is not a fundamental discriminant):

    rqf = BinaryQF_reduced_representatives(D, primitive_only=True)

    # compute needed precision
    #
    # NB: [https://arxiv.org/abs/0802.0979v1], quoting Enge (2006), is
    # incorrect.  Enge writes (2009-04-20 email to John Cremona) "The
    # source is my paper on class polynomials
    # [https://hal.inria.fr/inria-00001040] It was pointed out to me by
    # the referee after ANTS that the constant given there was
    # wrong. The final version contains a corrected constant on p.7
    # which is consistent with your example. It says:

    # "The logarithm of the absolute value of the coefficient in front
    # of X^j is bounded above by
    #
    # log (2*k_2) * h + pi * sqrt(|D|) * sum (1/A_i)
    #
    # independently of j", where k_2 \approx 10.163.

    h = len(rqf) # class number
    c1 = 3.05682737291380 # log(2*10.63)
    c2 = sum([1/RR(qf[0]) for qf in rqf], RR(0))
    prec =  c2*RR(3.142)*RR(D).abs().sqrt() + h*c1  # bound on log
    prec = prec * 1.45   # bound on log_2 (1/log(2) = 1.44..)
    prec = 10 + prec.ceil()  # allow for rounding error

    # set appropriate precision for further computing

    Dsqrt = D.sqrt(prec=prec)
    R = ComplexField(prec)['t']
    t = R.gen()
    pol = R(1)
    for qf in rqf:
        a, b, c = list(qf)
        tau = (b + Dsqrt) / (a << 1)
        pol *= (t - elliptic_j(tau))

    coeffs = [cof.real().round() for cof in pol.coefficients(sparse=False)]
    return IntegerRing()['x'](coeffs)


def is_HCP(f, check_monic_irreducible=True):
    r"""
    Return ``True``, `D` if ``f`` is the Hilbert Class Polynomial `H_D`, else ``False``, 0.

    INPUT:

    - ``f`` -- a polynomial in `\ZZ[X]`.
    - ``check_monic_irreducible`` (boolean, default ``True``) -- if
      ``True``, check that ``f`` is a monic, irreducible, integer
      polynomial.

    OUTPUT:

    (integer) -- `D` if ``f`` is the Hilbert Class Polynomial (HCP) `H_D`, esle 0.

    ALGORITHM:

    Cremona and Sutherland: Algorithm 2 of [CreSuth2023]_.

    EXAMPLES::

    Even for large degrees this is fast.  We test the largest
    discriminant of class number 100, for which the HCP has coefficients
    with thousands of digits::

        sage: from sage.schemes.elliptic_curves.cm import is_HCP
        sage: D = -1856563
        sage: D.class_number()
        100
        sage: H = hilbert_class_polynomial(D)
        sage: H.degree()
        100
        sage: max(c for c in H).ndigits()
        2774
        sage: is_HCP(H)
        (True, -1856563)

    Testing polynomials which are not HCPs is faster::

        sage: is_HCP(H+1)
        (False, 0)


    TESTS::

        sage: from sage.schemes.elliptic_curves.cm import is_HCP
        sage: all(is_HCP(hilbert_class_polynomial(D))==(True,D) for D in srange(-4,-100,-1) if D.is_discriminant())
        True
        sage: all(is_HCP(hilbert_class_polynomial(D)+1)==(False,0) for D in srange(-4,-100,-1) if D.is_discriminant())
        True

    """
    zero = ZZ(0)
    # optional check that input is monic and irreducible
    if check_monic_irreducible:
        try:
            if not (all(c in ZZ for c in f) and f.is_monic()):
                return (False, zero)
            f = f.change_ring(ZZ)
        except AttributeError:
            return (False, zero)

    h = f.degree()
    h2list = [d for d in h.divisors() if (d-h)%2==0 and d.prime_to_m_part(2)==1]
    pmin = 33 * (h**2 * (RR(h+2).log().log()+2)**2).ceil()
    # Guarantees 4*p > |D| for fundamental D under GRH
    p = pmin-1
    n = 0
    from sage.arith.all import next_prime
    while True:
        p = next_prime(p)
        n += 1
        fp = f.change_ring(GF(p))
        # Compute X^p-X mod fp manually, avoiding quotient ring which is slower
        r = zpow = z = fp.parent().gen()
        m = p>>1
        while m:
            zpow = (zpow**2) % fp
            if m & 1:
                r = (zpow * r) % fp
            m >>= 1
        # now r = X^p mod fp
        d = (r-z).gcd(fp).degree()  # number of roots mod p
        if d==0:
            continue
        if not fp.is_squarefree():
            continue
        if d<h and d not in h2list:
            return (False, zero)
        jp = fp.any_root(degree=-1, assume_squarefree=True)
        E = EllipticCurve(j=jp)
        if E.is_supersingular():
            continue
        D = E.endomorphism_discriminant_from_class_number(h)
        if not D:
            return (False, zero)
        return (True, D) if f == hilbert_class_polynomial(D) else (False, zero)

def OrderClassNumber(D0,h0,f):
    r"""
    Return the class number h(f**2 * D0), given h(D0)=h0.

    INPUT:

    - ``D0`` (integer) -- a negative fundamental discriminant
    - ``h0`` (integer) -- the class number of the (maximal) imaginary quadratic order of discriminant ``D0``
    - ``f`` (integer) -- a positive integer

    OUTPUT:

    (integer) the class number of the imaginary quadratic order of discriminant ``D0*f**2``

    ALGORITHM:

    We use the formula for the class number of the order `\mathcal{O}_{D}` in terms of the class number of the
     maximal order  `\mathcal{O}_{D_0}`; see [Cox1989]_ Theorem 7.24:

    .. MATH::

    h(D) = \frac{h(D_0)f}{[\mathcal{O}_{D_0}^\times:\mathcal{O}_{D}^\times]}\prod_{p\,|\,f}\left(1-\left(\frac{D_0}{p}\right)\frac{1}{p}\right)

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.cm import OrderClassNumber
        sage: D0 = -4
        sage: h = D0.class_number()
        sage: [OrderClassNumber(D0,h,f) for f in srange(1,20)]
        [1, 1, 2, 2, 2, 4, 4, 4, 6, 4, 6, 8, 6, 8, 8, 8, 8, 12, 10]
        sage: all([OrderClassNumber(D0,h,f) == (D0*f**2).class_number() for f in srange(1,20)])
        True

    """
    if not D0.is_fundamental_discriminant():
        raise ValueError("{} is not a fundamental discriminant".format(D0))
    if not D0.is_fundamental_discriminant() or f <= 0:
        raise ValueError("{} is not a positive integer".format(f))
    if f == 1:
        return h0
    ps = f.prime_divisors()
    from sage.misc.misc_c import prod
    from sage.arith.all import kronecker_symbol
    n = (f // prod(ps)) * prod(p-kronecker_symbol(D0,p) for p in ps)
    if D0 == -3:
        #assert h0 == 1 and n%3==0
        return n//3
    if D0 == -4:
        #assert h0 == 1 and n%2==0
        return n//2
    return n*h0

@cached_function
def cm_j_invariants(K, proof=None):
    r"""
    Return a list of all CM `j`-invariants in the field `K`.

    INPUT:

    - ``K`` -- a number field
    - ``proof`` -- (default: proof.number_field())

    OUTPUT:

    (list) -- A list of CM `j`-invariants in the field `K`.

    EXAMPLES::

        sage: cm_j_invariants(QQ)
        [-262537412640768000, -147197952000, -884736000, -12288000, -884736, -32768, -3375, 0, 1728, 8000, 54000, 287496, 16581375]

    Over imaginary quadratic fields there are no more than over `QQ`::

        sage: cm_j_invariants(QuadraticField(-1, 'i'))
        [-262537412640768000, -147197952000, -884736000, -12288000, -884736, -32768, -3375, 0, 1728, 8000, 54000, 287496, 16581375]

    Over real quadratic fields there may be more, for example::

        sage: len(cm_j_invariants(QuadraticField(5, 'a')))
        31

    Over number fields K of many higher degrees this also works::

        sage: K.<a> = NumberField(x^3 - 2)
        sage: cm_j_invariants(K)
        [-262537412640768000, -147197952000, -884736000,
         -884736, -32768, 8000, -3375, 16581375, 1728, 287496, 0,
         54000, -12288000,
         31710790944000*a^2 + 39953093016000*a + 50337742902000]
        sage: K.<a> = NumberField(x^4 - 2)
        sage: len(cm_j_invariants(K))
        23
    """
    return sorted(j for D, f, j in cm_j_invariants_and_orders(K, proof=proof))


@cached_function
def cm_j_invariants_and_orders(K, proof=None):
    r"""
    Return a list of all CM `j`-invariants in the field `K`, together with the associated orders.

    INPUT:

    - ``K`` -- a number field
    - ``proof`` -- (default: proof.number_field())

    OUTPUT:

    (list) A list of 3-tuples `(D,f,j)` where `j` is a CM
    `j`-invariant in `K` with quadratic fundamental discriminant `D`
    and conductor `f`.

    EXAMPLES::

        sage: cm_j_invariants_and_orders(QQ)
        [(-3, 3, -12288000), (-3, 2, 54000), (-3, 1, 0), (-4, 2, 287496), (-4, 1, 1728), (-7, 2, 16581375), (-7, 1, -3375), (-8, 1, 8000), (-11, 1, -32768), (-19, 1, -884736), (-43, 1, -884736000), (-67, 1, -147197952000), (-163, 1, -262537412640768000)]

    Over an imaginary quadratic field there are no more than over `QQ`::

        sage: cm_j_invariants_and_orders(QuadraticField(-1, 'i'))
        [(-163, 1, -262537412640768000), (-67, 1, -147197952000),
         (-43, 1, -884736000), (-19, 1, -884736), (-11, 1, -32768),
         (-8, 1, 8000), (-7, 1, -3375), (-7, 2, 16581375), (-4, 1, 1728),
         (-4, 2, 287496), (-3, 1, 0), (-3, 2, 54000), (-3, 3, -12288000)]

    Over real quadratic fields there may be more::

        sage: v = cm_j_invariants_and_orders(QuadraticField(5,'a')); len(v)
        31
        sage: [(D, f) for D, f, j in v if j not in QQ]
        [(-235, 1), (-235, 1), (-115, 1), (-115, 1), (-40, 1), (-40, 1),
         (-35, 1), (-35, 1), (-20, 1), (-20, 1), (-15, 1), (-15, 1), (-15, 2),
         (-15, 2), (-4, 5), (-4, 5), (-3, 5), (-3, 5)]

    Over number fields K of many higher degrees this also works::

        sage: K.<a> = NumberField(x^3 - 2)
        sage: cm_j_invariants_and_orders(K)
        [(-163, 1, -262537412640768000), (-67, 1, -147197952000),
         (-43, 1, -884736000), (-19, 1, -884736), (-11, 1, -32768),
         (-8, 1, 8000), (-7, 1, -3375), (-7, 2, 16581375), (-4, 1, 1728),
         (-4, 2, 287496), (-3, 1, 0), (-3, 2, 54000), (-3, 3, -12288000),
         (-3, 6, 31710790944000*a^2 + 39953093016000*a + 50337742902000)]
    """
    if K == QQ:
        return [(ZZ(d), ZZ(f), ZZ(j)) for d, f, j in [
            (-3, 3, -12288000),
            (-3, 2, 54000),
            (-3, 1, 0),
            (-4, 2, 287496),
            (-4, 1, 1728),
            (-7, 2, 16581375),
            (-7, 1, -3375),
            (-8, 1, 8000),
            (-11, 1, -32768),
            (-19, 1, -884736),
            (-43, 1, -884736000),
            (-67, 1, -147197952000),
            (-163, 1, -262537412640768000)]]

    # Get the list of CM orders that could possibly have Hilbert class
    # polynomial F(x) with a root in K.  If F(x) has a root alpha in K,
    # then F is the minimal polynomial of alpha in K, so the degree of
    # F(x) is at most [K:QQ].
    dlist = sorted(Df for v in discriminants_with_bounded_class_number(K.degree(), proof=proof).values() for Df in v)

    return [(D, f, j) for D, f in dlist
            for j in hilbert_class_polynomial(D*f*f).roots(K, multiplicities=False)]


@cached_function
def cm_orders(h, proof=None):
    """
    Return a list of all pairs `(D,f)` where there is a CM order of
    discriminant `D f^2` with class number h, with `D` a fundamental
    discriminant.

    INPUT:

    - `h` -- positive integer
    - ``proof`` -- (default: proof.number_field())

    OUTPUT:

    - list of 2-tuples `(D,f)`

    EXAMPLES::

        sage: cm_orders(0)
        []
        sage: v = cm_orders(1); v
        [(-3, 3), (-3, 2), (-3, 1), (-4, 2), (-4, 1), (-7, 2), (-7, 1), (-8, 1), (-11, 1), (-19, 1), (-43, 1), (-67, 1), (-163, 1)]
        sage: type(v[0][0]), type(v[0][1])
        (<... 'sage.rings.integer.Integer'>, <... 'sage.rings.integer.Integer'>)
        sage: v = cm_orders(2); v
         [(-3, 7), (-3, 5), (-3, 4), (-4, 5), (-4, 4), (-4, 3), (-7, 4), (-8, 3), (-8, 2), (-11, 3), (-15, 2), (-15, 1), (-20, 1), (-24, 1), (-35, 1), (-40, 1), (-51, 1), (-52, 1), (-88, 1), (-91, 1), (-115, 1), (-123, 1), (-148, 1), (-187, 1), (-232, 1), (-235, 1), (-267, 1), (-403, 1), (-427, 1)]
        sage: len(v)
        29
        sage: set([hilbert_class_polynomial(D*f^2).degree() for D,f in v])
        {2}

    Any degree up to 100 is implemented, but may be prohibitively slow::

        sage: cm_orders(3)
        [(-3, 9), (-3, 6), (-11, 2), (-19, 2), (-23, 2), (-23, 1), (-31, 2), (-31, 1), (-43, 2), (-59, 1), (-67, 2), (-83, 1), (-107, 1), (-139, 1), (-163, 2), (-211, 1), (-283, 1), (-307, 1), (-331, 1), (-379, 1), (-499, 1), (-547, 1), (-643, 1), (-883, 1), (-907, 1)]
        sage: len(cm_orders(4))
        84
    """
    h = Integer(h)
    if h <= 0:
        # trivial case
        return []
    # Get information for all discriminants then throw away everything
    # but for h.  If this is replaced by a table it will be faster,
    # but not now.   (David Kohel is rumored to have a large table.)
    return discriminants_with_bounded_class_number(h, proof=proof)[h]

# Table from Mark Watkins paper "Class numbers of imaginary quadratic fields".
# I extracted this by cutting/pasting from the pdf, and running this program:
# z = {}
# for X in open('/Users/wstein/tmp/a.txt').readlines():
#    if len(X.strip()):
#        v = [int(a) for a in X.split()]
#        for i in range(5):
#            z[v[3*i]]=(v[3*i+2], v[3*i+1])


watkins_table = {1: (163, 9), 2: (427, 18), 3: (907, 16), 4: (1555, 54), 5: (2683, 25),
                 6: (3763, 51), 7: (5923, 31), 8: (6307, 131), 9: (10627, 34), 10:
                 (13843, 87), 11: (15667, 41), 12: (17803, 206), 13: (20563, 37), 14:
                 (30067, 95), 15: (34483, 68), 16: (31243, 322), 17: (37123, 45), 18:
                 (48427, 150), 19: (38707, 47), 20: (58507, 350), 21: (61483, 85), 22:
                 (85507, 139), 23: (90787, 68), 24: (111763, 511), 25: (93307, 95), 26:
                 (103027, 190), 27: (103387, 93), 28: (126043, 457), 29: (166147, 83),
                 30: (134467, 255), 31: (133387, 73), 32: (164803, 708), 33: (222643, 101),
                 34: (189883, 219), 35: (210907, 103), 36: (217627, 668), 37:
                 (158923, 85), 38: (289963, 237), 39: (253507, 115), 40: (260947, 912),
                 41: (296587, 109), 42: (280267, 339), 43: (300787, 106), 44: (319867, 691),
                 45: (308323, 154), 46: (462883, 268), 47: (375523, 107), 48:
                 (335203, 1365), 49: (393187, 132), 50: (389467, 345), 51: (546067, 159),
                 52: (439147, 770), 53: (425107, 114), 54: (532123, 427), 55: (452083,163),
                 56: (494323, 1205), 57: (615883, 179), 58: (586987, 291),
                 59:(474307, 128), 60: (662803, 1302), 61: (606643, 132), 62: (647707, 323),
                 63: (991027, 216), 64: (693067, 1672), 65: (703123, 164), 66: (958483, 530),
                 67: (652723, 120), 68: (819163, 976), 69: (888427, 209), 70:(811507, 560),
                 71: (909547, 150), 72: (947923, 1930), 73: (886867, 119),
                 74: (951043, 407), 75: (916507, 237), 76: (1086187, 1075), 77: (1242763, 216),
                 78: (1004347, 561), 79: (1333963, 175), 80: (1165483, 2277), 81: (1030723, 228),
                 82: (1446547, 402), 83: (1074907, 150), 84: (1225387,1715),
                 85: (1285747, 221), 86: (1534723, 472), 87: (1261747, 222),
                 88:(1265587, 1905), 89: (1429387, 192), 90: (1548523, 801),
                 91: (1391083,214), 92: (1452067, 1248), 93: (1475203, 262), 94: (1587763, 509),
                 95:(1659067, 241), 96: (1684027, 3283), 97: (1842523, 185), 98: (2383747,580),
                 99: (1480627, 289), 100: (1856563, 1736)}


def largest_fundamental_disc_with_class_number(h):
    """
    Return largest absolute value of any fundamental discriminant with
    class number `h`, and the number of fundamental discriminants with
    that class number.  This is known for `h` up to 100, by work of Mark
    Watkins.

    INPUT:

    - `h` -- integer

    EXAMPLES::

        sage: sage.schemes.elliptic_curves.cm.largest_fundamental_disc_with_class_number(0)
        (0, 0)
        sage: sage.schemes.elliptic_curves.cm.largest_fundamental_disc_with_class_number(1)
        (163, 9)
        sage: sage.schemes.elliptic_curves.cm.largest_fundamental_disc_with_class_number(2)
        (427, 18)
        sage: sage.schemes.elliptic_curves.cm.largest_fundamental_disc_with_class_number(10)
        (13843, 87)
        sage: sage.schemes.elliptic_curves.cm.largest_fundamental_disc_with_class_number(100)
        (1856563, 1736)
        sage: sage.schemes.elliptic_curves.cm.largest_fundamental_disc_with_class_number(101)
        Traceback (most recent call last):
        ...
        NotImplementedError: largest discriminant not known for class number 101
    """
    h = Integer(h)
    if h <= 0:
        # very easy special case
        return Integer(0), Integer(0)
    try:
        # simply look up the answer in Watkins's table.
        B, c = watkins_table[h]
        return (Integer(B), Integer(c))
    except KeyError:
        # nobody knows, since I guess Watkins's is state of the art.
        raise NotImplementedError("largest discriminant not known for class number %s" % h)


@cached_function
def discriminants_with_bounded_class_number(hmax, B=None, proof=None):
    r"""
    Return dictionary with keys class numbers `h\le hmax` and values the
    list of all pairs `(D, f)`, with `D<0` a fundamental discriminant such
    that `Df^2` has class number `h`.  If the optional bound `B` is given,
    return only those pairs with fundamental `|D| \le B`, though `f` can
    still be arbitrarily large.

    INPUT:

    - ``hmax`` -- integer
    - `B` -- integer or None; if None returns all pairs
    - ``proof`` -- this code calls the PARI function :pari:`qfbclassno`, so it
      could give wrong answers when ``proof``==``False``.  The default is
      whatever ``proof.number_field()`` is.  If ``proof==False`` and `B` is
      ``None``, at least the number of discriminants is correct, since it
      is double checked with Watkins's table.

    OUTPUT:

    - dictionary

    In case `B` is not given, we use Mark Watkins's: "Class numbers of
    imaginary quadratic fields" to compute a `B` that captures all `h`
    up to `hmax` (only available for `hmax\le100`).

    EXAMPLES::

        sage: v = sage.schemes.elliptic_curves.cm.discriminants_with_bounded_class_number(3)
        sage: sorted(v)
        [1, 2, 3]
        sage: v[1]
        [(-3, 3), (-3, 2), (-3, 1), (-4, 2), (-4, 1), (-7, 2), (-7, 1), (-8, 1), (-11, 1), (-19, 1), (-43, 1), (-67, 1), (-163, 1)]
        sage: v[2]
        [(-3, 7), (-3, 5), (-3, 4), (-4, 5), (-4, 4), (-4, 3), (-7, 4), (-8, 3), (-8, 2), (-11, 3), (-15, 2), (-15, 1), (-20, 1), (-24, 1), (-35, 1), (-40, 1), (-51, 1), (-52, 1), (-88, 1), (-91, 1), (-115, 1), (-123, 1), (-148, 1), (-187, 1), (-232, 1), (-235, 1), (-267, 1), (-403, 1), (-427, 1)]
        sage: v[3]
        [(-3, 9), (-3, 6), (-11, 2), (-19, 2), (-23, 2), (-23, 1), (-31, 2), (-31, 1), (-43, 2), (-59, 1), (-67, 2), (-83, 1), (-107, 1), (-139, 1), (-163, 2), (-211, 1), (-283, 1), (-307, 1), (-331, 1), (-379, 1), (-499, 1), (-547, 1), (-643, 1), (-883, 1), (-907, 1)]
        sage: v = sage.schemes.elliptic_curves.cm.discriminants_with_bounded_class_number(8, proof=False)
        sage: sorted(len(v[h]) for h in v)
        [13, 25, 29, 29, 38, 84, 101, 208]

    Find all class numbers for discriminant up to 50::

        sage: sage.schemes.elliptic_curves.cm.discriminants_with_bounded_class_number(hmax=5, B=50)
        {1: [(-3, 3), (-3, 2), (-3, 1), (-4, 2), (-4, 1), (-7, 2), (-7, 1), (-8, 1), (-11, 1), (-19, 1), (-43, 1)], 2: [(-3, 7), (-3, 5), (-3, 4), (-4, 5), (-4, 4), (-4, 3), (-7, 4), (-8, 3), (-8, 2), (-11, 3), (-15, 2), (-15, 1), (-20, 1), (-24, 1), (-35, 1), (-40, 1)], 3: [(-3, 9), (-3, 6), (-11, 2), (-19, 2), (-23, 2), (-23, 1), (-31, 2), (-31, 1), (-43, 2)], 4: [(-3, 13), (-3, 11), (-3, 8), (-4, 10), (-4, 8), (-4, 7), (-4, 6), (-7, 8), (-7, 6), (-7, 3), (-8, 6), (-8, 4), (-11, 5), (-15, 4), (-19, 5), (-19, 3), (-20, 3), (-20, 2), (-24, 2), (-35, 3), (-39, 2), (-39, 1), (-40, 2), (-43, 3)], 5: [(-47, 2), (-47, 1)]}
    """
    # imports that are needed only for this function
    from sage.structure.proof.proof import get_flag

    # deal with input defaults and type checking
    proof = get_flag(proof, 'number_field')
    hmax = Integer(hmax)

    # T stores the output
    T = {}

    # Easy case -- instead of giving error, give meaningful output
    if hmax < 1:
        return T

    if B is None:
        # Determine how far we have to go by applying Watkins's theorem.
        v = [largest_fundamental_disc_with_class_number(h) for h in range(1, hmax+1)]
        B = max([b for b,_ in v])
        fund_count = [0] + [cnt for _,cnt in v]
    else:
        # Nothing to do -- set to None so we can use this later to know not
        # to do a double check about how many we find.
        fund_count = None
        B = Integer(B)

    if B <= 2:
        # This is an easy special case, since there are no fundamental discriminants
        # this small.
        return T

    # This lower bound gets used in an inner loop below.
    from math import log

    def lb(f):
        """Lower bound on euler_phi."""
        # 1.79 > e^gamma = 1.7810724...
        if f <= 1:
            return 0  # don't do log(log(1)) = log(0)
        llf = log(log(f))
        return f/(1.79*llf + 3.0/llf)

    for D in range(-B, -2):
        D = Integer(D)
        if D.is_fundamental_discriminant():
            h_D = D.class_number(proof)
            # For each fundamental discriminant D, loop through the f's such
            # that h(D*f^2) could possibly be <= hmax.  As explained to me by Cremona,
            # we have h(D*f^2) >= (1/c)*h(D)*phi_D(f) >= (1/c)*h(D)*euler_phi(f), where
            # phi_D(f) is like euler_phi(f) but the factor (1-1/p) is replaced
            # by a factor of (1-kr(D,p)*p), where kr(D/p) is the Kronecker symbol.
            # The factor c is 1 unless D=-4 and f>1 (when c=2) or D=-3 and f>1 (when c=3).
            # Since (1-1/p) <= 1 and (1-1/p) <= (1+1/p), we see that
            #     euler_phi(f) <= phi_D(f).
            #
            # We have the following analytic lower bound on euler_phi:
            #
            #     euler_phi(f) >= lb(f) = f / (exp(euler_gamma)*log(log(n)) + 3/log(log(n))).
            #
            # See Theorem 8 of Peter Clark's
            #   http://math.uga.edu/~pete/4400arithmeticorders.pdf
            # which is a consequence of Theorem 15 of
            # [Rosser and Schoenfeld, 1962].
            #
            # By Calculus, we see that the lb(f) is an increasing function of f >= 2.
            #
            # NOTE: You can visibly "see" that it is a lower bound in Sage with
            #   lb(n) = n/(exp(euler_gamma)*log(log(n)) + 3/log(log(n)))
            #   plot(lb, (n, 1, 10^4), color='red') + plot(lambda x: euler_phi(int(x)), 1, 10^4).show()
            #
            # So we consider f=1,2,..., until the first f with lb(f)*h_D > c*h_max.
            # (Note that lb(f) is <= 0 for f=1,2, so nothing special is needed there.)
            #
            # TODO: Maybe we could do better using a bound for phi_D(f).
            #
            f = Integer(1)
            chmax=hmax
            if D==-3:
                chmax*=3
            else:
                if D==-4:
                    chmax*=2
            while lb(f)*h_D <= chmax:
                h = OrderClassNumber(D,h_D,f)
                # If the class number of this order is within the range, then use it.
                if h <= hmax:
                    z = (D, f)
                    if h in T:
                        T[h].append(z)
                    else:
                        T[h] = [z]
                f += 1

    for h in T:
        T[h] = list(reversed(T[h]))

    if fund_count is not None:
        # Double check that we found the right number of fundamental
        # discriminants; we might as well, since Watkins provides this
        # data.
        for h in T:
            if len([DD for DD, ff in T[h] if ff == 1]) != fund_count[h]:
                raise RuntimeError("number of discriminants inconsistent with Watkins's table")

    return T


@cached_function
def is_cm_j_invariant(j, method='CremonaSutherland'):
    r"""Return whether or not this is a CM `j`-invariant, and the CM discriminant if it is.

    INPUT:

    - ``j`` -- an element of a number field `K`

    - ``method`` (string, default 'CremonaSutherland') -- the method
      used, either 'CremonaSutherland' (the default, very much faster
      for all but very small degrees), 'exhaustive' or 'reduction'

    OUTPUT:

    A pair (bool, (d,f)) which is either (False, None) if `j` is not a
    CM j-invariant or (True, (d,f)) if `j` is the `j`-invariant of the
    imaginary quadratic order of discriminant `D=df^2` where `d` is
    the associated fundamental discriminant and `f` the index.

    ALGORITHM:

    The default algorithm used is to test whether the minimal
    polynomial of ``j`` is a Hilbert CLass Polynomail, using
    :func:`is_HCP` which implements Algorithm 2 of [CreSuth2023]_ by
    Cremona and Sutherland.

    Two older methods are available, both of which are much slower
    except for very small degrees.

    Method 'exhaustive' makes use of the complete and unconditionsl classification of
    all orders of class number up to 100, and hence will raise an
    error if `j` is an algebraic integer of degree greater than
    this.

    Method 'reduction' constructs an elliptic curve over the number
    field `\QQ(j)` and computes its traces of Frobenius at several
    primes of degree 1.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.cm import is_cm_j_invariant
        sage: is_cm_j_invariant(0)
        (True, (-3, 1))
        sage: is_cm_j_invariant(8000)
        (True, (-8, 1))

        sage: K.<a> = QuadraticField(5)
        sage: is_cm_j_invariant(282880*a + 632000)
        (True, (-20, 1))
        sage: K.<a> = NumberField(x^3 - 2)
        sage: is_cm_j_invariant(31710790944000*a^2 + 39953093016000*a + 50337742902000)
        (True, (-3, 6))

    An example of large degree.  This is only possible using the default algorithm::

        sage: from sage.schemes.elliptic_curves.cm import is_cm_j_invariant
        sage: D = -1856563
        sage: H = hilbert_class_polynomial(D)
        sage: H.degree()
        100
        sage: K.<j> = NumberField(H)
        sage: is_cm_j_invariant(j)
        (True, (-1856563, 1))

    TESTS::

        sage: from sage.schemes.elliptic_curves.cm import is_cm_j_invariant
        sage: all(is_cm_j_invariant(j) == (True, (d,f)) for d,f,j in cm_j_invariants_and_orders(QQ))
        True

    """
    # First we check that j is an algebraic number:
    from sage.rings.all import NumberFieldElement, NumberField
    if not isinstance(j, NumberFieldElement) and j not in QQ:
        raise NotImplementedError("is_cm_j_invariant() is only implemented for number field elements")

    # for j in ZZ we have a lookup-table:

    if j in ZZ:
        j = ZZ(j)
        table = dict([(jj,(d,f)) for d,f,jj in cm_j_invariants_and_orders(QQ)])
        if j in table:
            return True, table[j]
        return False, None

    # Otherwise if j is in Q then it is not integral so is not CM:

    if j in QQ:
        return False, None

    # Next we find its minimal polynomial of j:

    if j.parent().absolute_degree()==2:
        jpol = j.absolute_minpoly() # no algorithm parameter
    else:
        jpol = j.absolute_minpoly(algorithm='pari')

    # If it does not have integer coefficients then j is not integral, hence not CM:

    if not all(c in ZZ for c in jpol):
        return False, None

    # Otherwise test whether it is a Hilbert Class Polynomial
    # (using the fact that we know that it is monic and irreducible):

    if method == 'CremonaSutherland':
        flag, D = is_HCP(jpol, check_monic_irreducible=False)
        if flag:
            D0 = D.squarefree_part()
            if D0%4 !=1:
                D0 *= 4
            f = ZZ(D//D0).isqrt()
            return (True, (D0,f))
        else:
            return (False, None)

    h = jpol.degree()
    if method in ['exhaustive', 'old']:
        if h>100:
            raise NotImplementedError("CM data only available for class numbers up to 100")
        for d,f in cm_orders(h):
            if jpol == hilbert_class_polynomial(d*f**2):
                return (True, (d,f))
        return (False, None)

    if method not in ['reduction', 'new']:
        raise ValueError("Invalid method {} in is_cm_j_invariant".format(method))

    # Now we use the reduction method

    # If the degree h is less than the degree of j.parent() we recreate j as an element
    # of Q(j, and replace j by a clone whose parent is Q(j), if necessary:

    K = j.parent()
    if h < K.absolute_degree():
        K = NumberField(jpol, 'j')
        j = K.gen()

    # Construct an elliptic curve with j-invariant j, with
    # integral model:

    E = EllipticCurve(j=j).integral_model()
    D = E.discriminant()
    prime_bound = 1000 # test primes of degree 1 up to this norm
    max_primes =    20 # test at most this many primes
    num_prime = 0
    cmd = 0
    cmf = 0

    # Test primes of good reduction.  If E has CM then for half the
    # primes P we will have a_P=0, and for all other prime P the CM
    # field is Q(sqrt(a_P^2-4N(P))).  Hence if these fields are
    # different for two primes then E does not have CM.  If they are
    # all equal for the primes tested, then we have a candidate CM
    # field.  Moreover the discriminant of the endomorphism ring
    # divides all the values a_P^2-4N(P), since that is the
    # discriminant of the order containing the Frobenius at P.  So we
    # end up with a finite number (usually one) of candidate
    # discriminants to test.  Each is tested by checking that its class
    # number is h, and if so then that j is a root of its Hilbert
    # class polynomial.  In practice non CM curves will be eliminated
    # by the local test at a small number of primes (probably just 2).

    for P in K.primes_of_degree_one_iter(prime_bound):
        if num_prime > max_primes:
            if cmd: # we have a candidate CM field already
                break
            else:   # we need to try more primes
                max_primes *=2
        if D.valuation(P)>0: # skip bad primes
            continue
        aP = E.reduction(P).trace_of_frobenius()
        if aP == 0: # skip supersingular primes
            continue
        num_prime += 1
        DP = aP**2 - 4*P.norm()
        dP = DP.squarefree_part()
        fP = ZZ(DP//dP).isqrt()
        if cmd==0:      # first one, so store d and f
            cmd = dP
            cmf = fP
        elif cmd != dP: # inconsistent with previous
            return (False, None)
        else:           # consistent d, so update f
            cmf = cmf.gcd(fP)

    if cmd==0: # no conclusion, we found no degree 1 primes, revert to default method
        return is_cm_j_invariant(j, method='CremonaSutherland')

    # it looks like cm by disc cmd * f**2 where f divides cmf

    if cmd % 4 != 1:
        cmd = cmd * 4
        cmf = cmf // 2

    # Now we must check if h(cmd*f**2)==h for f|cmf; if so we check
    # whether j is a root of the associated Hilbert class polynomial.
    h0 = cmd.class_number()
    for f in cmf.divisors(): # only positive divisors
        if h != OrderClassNumber(cmd,h0,f):
            continue
        if jpol == hilbert_class_polynomial(cmd*f**2):
            return (True, (cmd, f))
    return (False, None)
