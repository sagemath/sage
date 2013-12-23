"""
Complex multiplication for elliptic curves

This module implements the functions

- ``hilbert_class_polynomial``
- ``cm_j_invariants``
- ``cm_orders``
- ``discriminants_with_bounded_class_number``
- ``cm_j_invariants_and_orders``
- ``largest_fundamental_disc_with_class_number``

AUTHORS:

- Robert Bradshaw
- John Cremona
- William Stein

"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.interfaces.all import magma
from sage.rings.all import (Integer,
                            QQ,
                            IntegerRing,
                            is_fundamental_discriminant,
                            PolynomialRing)

def hilbert_class_polynomial(D, algorithm=None):
    r"""
    Returns the Hilbert class polynomial for discriminant `D`.

    INPUT:

    - ``D`` (int) -- a negative integer congruent to 0 or 1 modulo 4.

    - ``algorithm`` (string, default None) -- if "sage" then use the Sage implementation; if "magma" then call Magma (if available).

    OUTPUT:

    (integer polynomial) The Hilbert class polynomial for the
    discriminant `D`.

    ALGORITHM:

    - If ``algorithm`` = "sage": Use complex approximations to the roots.

    - If ``algorithm`` = "magma": Call the appropriate Magma function.

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
        sage: hilbert_class_polynomial(-163, algorithm="magma") # optional - magma
        x + 262537412640768000

    """
    if algorithm is None:
        algorithm = "sage"

    D = Integer(D)
    if D >= 0:
        raise ValueError, "D (=%s) must be negative"%D
    if not (D%4 in [0,1]):
         raise ValueError, "D (=%s) must be a discriminant"%D

    if algorithm == "magma":
        magma.eval("R<x> := PolynomialRing(IntegerRing())")
        f = str(magma.eval("HilbertClassPolynomial(%s)"%D))
        return IntegerRing()['x'](f)

    if algorithm != "sage":
        raise ValueError, "%s is not a valid algorithm"%algorithm

    from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
    from sage.rings.all import RR, ZZ, ComplexField
    from sage.functions.all import elliptic_j

    # get all primitive reduced quadratic forms, (necessary to exclude
    # imprimitive forms when D is not a fundamental discriminant):

    rqf = BinaryQF_reduced_representatives(D, primitive_only=True)

    # compute needed precision
    #
    # NB: [http://arxiv.org/abs/0802.0979v1], quoting Enge (2006), is
    # incorrect.  Enge writes (2009-04-20 email to John Cremona) "The
    # source is my paper on class polynomials
    # [http://hal.inria.fr/inria-00001040] It was pointed out to me by
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
        tau = (b+Dsqrt)/(a<<1)
        pol *=  (t - elliptic_j(tau))

    coeffs = [cof.real().round() for cof in pol.coeffs()]
    return IntegerRing()['x'](coeffs)


def cm_j_invariants(K, proof=None):
    r"""
    Return a list of all CM `j`-invariants in the field `K`.

    INPUT:

    - ``K`` -- a number field
    - ``proof`` -- (default: proof.number_field())

    OUTPUT:

    (list) -- A list of CM `j`-invariants in the field `K`.

    EXAMPLE::

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
        [-12288000, 54000, 0, 287496, 1728, 16581375, -3375, 8000, -32768, -884736, -884736000, -147197952000, -262537412640768000, 31710790944000*a^2 + 39953093016000*a + 50337742902000]
        sage: K.<a> = NumberField(x^4 - 2)
        sage: len(cm_j_invariants(K))
        23
    """
    return list(sorted([j for D,f,j in cm_j_invariants_and_orders(K, proof=proof)]))

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

    EXAMPLE::

        sage: cm_j_invariants_and_orders(QQ)
        [(-3, 3, -12288000), (-3, 2, 54000), (-3, 1, 0), (-4, 2, 287496), (-4, 1, 1728), (-7, 2, 16581375), (-7, 1, -3375), (-8, 1, 8000), (-11, 1, -32768), (-19, 1, -884736), (-43, 1, -884736000), (-67, 1, -147197952000), (-163, 1, -262537412640768000)]

    Over an imaginary quadratic field there are no more than over `QQ`::

        sage: cm_j_invariants_and_orders(QuadraticField(-1, 'i'))
        [(-3, 3, -12288000), (-3, 2, 54000), (-3, 1, 0), (-4, 2, 287496), (-4, 1, 1728), (-7, 2, 16581375), (-7, 1, -3375), (-8, 1, 8000), (-11, 1, -32768), (-19, 1, -884736), (-43, 1, -884736000), (-67, 1, -147197952000), (-163, 1, -262537412640768000)]

    Over real quadratic fields there may be more::

        sage: v = cm_j_invariants_and_orders(QuadraticField(5,'a')); len(v)
        31
        sage: [(D,f) for D,f,j in v if j not in QQ]
        [(-3, 5), (-3, 5), (-4, 5), (-4, 5), (-15, 2), (-15, 2), (-15, 1), (-15, 1), (-20, 1), (-20, 1), (-35, 1), (-35, 1), (-40, 1), (-40, 1), (-115, 1), (-115, 1), (-235, 1), (-235, 1)]

    Over number fields K of many higher degrees this also works::

        sage: K.<a> = NumberField(x^3 - 2)
        sage: cm_j_invariants_and_orders(K)
        [(-3, 3, -12288000), (-3, 2, 54000), (-3, 1, 0), (-4, 2, 287496), (-4, 1, 1728), (-7, 2, 16581375), (-7, 1, -3375), (-8, 1, 8000), (-11, 1, -32768), (-19, 1, -884736), (-43, 1, -884736000), (-67, 1, -147197952000), (-163, 1, -262537412640768000), (-3, 6, 31710790944000*a^2 + 39953093016000*a + 50337742902000)]
    """
    # Get the list of CM orders that could possibly have Hilbert class
    # polynomial F(x) with a root in K.  If F(x) has a root alpha in K,
    # then F is the minimal polynomial of alpha in K, so the degree of
    # F(x) is at most [K:QQ].
    dlist = sum([v for h,v in discriminants_with_bounded_class_number(K.degree(), proof=proof).iteritems()], [])

    return [(D,f,j) for D, f in dlist
             for j in hilbert_class_polynomial(D*f*f).roots(K, multiplicities=False)]


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
        (<type 'sage.rings.integer.Integer'>, <type 'sage.rings.integer.Integer'>)
        sage: v = cm_orders(2); v
         [(-3, 7), (-3, 5), (-3, 4), (-4, 5), (-4, 4), (-4, 3), (-7, 4), (-8, 3), (-8, 2), (-11, 3), (-15, 2), (-15, 1), (-20, 1), (-24, 1), (-35, 1), (-40, 1), (-51, 1), (-52, 1), (-88, 1), (-91, 1), (-115, 1), (-123, 1), (-148, 1), (-187, 1), (-232, 1), (-235, 1), (-267, 1), (-403, 1), (-427, 1)]
        sage: len(v)
        29
        sage: set([hilbert_class_polynomial(D*f^2).degree() for D,f in v])
        set([2])

    Any degree up to 100 is implemented, but may be prohibitively slow::

        sage: cm_orders(3)
        [(-3, 9), (-3, 6), (-11, 2), (-19, 2), (-23, 2), (-23, 1), (-31, 2), (-31, 1), (-43, 2), (-59, 1), (-67, 2), (-83, 1), (-107, 1), (-139, 1), (-163, 2), (-211, 1), (-283, 1), (-307, 1), (-331, 1), (-379, 1), (-499, 1), (-547, 1), (-643, 1), (-883, 1), (-907, 1)]
        sage: len(cm_orders(4))
        84
    """
    h = Integer(h)
    T = None
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
        raise NotImplementedError, "largest discriminant not known for class number %s"%h

def discriminants_with_bounded_class_number(hmax, B=None, proof=None):
    """
    Return dictionary with keys class numbers `h\le hmax` and values the
    list of all pairs `(D, f)`, with `D<0` a fundamental discriminant such
    that `Df^2` has class number `h`.  If the optional bound `B` is given,
    return only those pairs with fundamental `|D| \le B`, though `f` can
    still be arbitrarily large.

    INPUT:

    - ``hmax`` -- integer
    - `B` -- integer or None; if None returns all pairs
    - ``proof`` -- this code calls the PARI function ``qfbclassno``, so it
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
        sage: v.keys()
        [1, 2, 3]
        sage: v[1]
        [(-3, 3), (-3, 2), (-3, 1), (-4, 2), (-4, 1), (-7, 2), (-7, 1), (-8, 1), (-11, 1), (-19, 1), (-43, 1), (-67, 1), (-163, 1)]
        sage: v[2]
        [(-3, 7), (-3, 5), (-3, 4), (-4, 5), (-4, 4), (-4, 3), (-7, 4), (-8, 3), (-8, 2), (-11, 3), (-15, 2), (-15, 1), (-20, 1), (-24, 1), (-35, 1), (-40, 1), (-51, 1), (-52, 1), (-88, 1), (-91, 1), (-115, 1), (-123, 1), (-148, 1), (-187, 1), (-232, 1), (-235, 1), (-267, 1), (-403, 1), (-427, 1)]
        sage: v[3]
        [(-3, 9), (-3, 6), (-11, 2), (-19, 2), (-23, 2), (-23, 1), (-31, 2), (-31, 1), (-43, 2), (-59, 1), (-67, 2), (-83, 1), (-107, 1), (-139, 1), (-163, 2), (-211, 1), (-283, 1), (-307, 1), (-331, 1), (-379, 1), (-499, 1), (-547, 1), (-643, 1), (-883, 1), (-907, 1)]
        sage: v = sage.schemes.elliptic_curves.cm.discriminants_with_bounded_class_number(8, proof=False)
        sage: [len(v[h]) for h in v.keys()]
        [13, 29, 25, 84, 29, 101, 38, 208]

    Find all class numbers for discriminant up to 50::

        sage: sage.schemes.elliptic_curves.cm.discriminants_with_bounded_class_number(hmax=5, B=50)
        {1: [(-3, 3), (-3, 2), (-3, 1), (-4, 2), (-4, 1), (-7, 2), (-7, 1), (-8, 1), (-11, 1), (-19, 1), (-43, 1)], 2: [(-3, 7), (-3, 5), (-3, 4), (-4, 5), (-4, 4), (-4, 3), (-7, 4), (-8, 3), (-8, 2), (-11, 3), (-15, 2), (-15, 1), (-20, 1), (-24, 1), (-35, 1), (-40, 1)], 3: [(-3, 9), (-3, 6), (-11, 2), (-19, 2), (-23, 2), (-23, 1), (-31, 2), (-31, 1), (-43, 2)], 4: [(-3, 13), (-3, 11), (-3, 8), (-4, 10), (-4, 8), (-4, 7), (-4, 6), (-7, 8), (-7, 6), (-7, 3), (-8, 6), (-8, 4), (-11, 5), (-15, 4), (-19, 5), (-19, 3), (-20, 3), (-20, 2), (-24, 2), (-35, 3), (-39, 2), (-39, 1), (-40, 2), (-43, 3)], 5: [(-47, 2), (-47, 1)]}
    """
    # imports that are needed only for this function
    from sage.structure.proof.proof import get_flag
    from sage.libs.pari.all import pari
    import math
    from sage.misc.functional import round

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
        if f <= 1: return 0  # don't do log(log(1)) = log(0)
        return f/(1.79*log(log(f)) + 3.0/log(log(f)))

    # We define a little function to compute the class number of
    # discriminant d quickly, using pari, with proof inherited from
    # the containing scope.  Fast classno functionality should
    # probably be moved elsewhere in Sage, but I'm not sure where.
    # The same line also occurs in the number fields code, but to use
    # it requires making a number field, which is very slow, and this
    # function must be fast, since it is the main bottleneck.
    def classno(d):
        """Return the class number of the order of discriminant d."""
        # There is currently no qfbclassno method in gen.pyx, hence the string.
        return Integer(pari('qfbclassno(%s,%s)'%(d, 1 if proof else 0)))

    for D in range(-B, -2):
        if is_fundamental_discriminant(D):
            h_D = classno(D)
            # For each fundamental discrimant D, loop through the f's such
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
            # TODO: Maybe we could do better using a bound for for phi_D(f).
            #
            f = 1
            chmax=hmax
            if D==-3:
                chmax*=3
            else:
                if D==-4:
                    chmax*=2
            while lb(f)*h_D <= chmax:
                if f == 1:
                    h = h_D
                else:
                    h = classno(D*f*f)
                # If the class number of this order is within the range, then
                # use it.  (NOTE: In some cases there is a simple relation between
                # the class number for D and D*f^2, and this could be used to
                # optimize this inner loop a little.)
                if h <= hmax:
                    z = (Integer(D), Integer(f))
                    if T.has_key(h):
                        T[h].append(z)
                    else:
                        T[h] = [z]
                f += 1

    for h in T.keys():
        T[h] = list(reversed(T[h]))

    if fund_count is not None:
        # Double check that we found the right number of fundamental
        # discriminants; we might as well, since Watkins provides this
        # data.
        for h in T.keys():
            if len([D for D,f in T[h] if f==1]) != fund_count[h]:
                raise RuntimeError, "number of discriminants inconsistent with Watkins's table"

    return T




