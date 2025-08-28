r"""
Global and semi-global minimal models for elliptic curves over number fields

When `E` is an elliptic curve defined over a number field `K` of class
number 1, then it has a global minimal model, and we have a method to
compute it, namely ``E.global_minimal_model()``.  Until Sage-6.7 this was
done using Tate's algorithm to minimise one prime at a time without
affecting the other primes.  When the class number is not 1, a
different approach is used.

In the general case global minimal models may or may not exist. This
module includes functions to determine this, and to find a global
minimal model when it does exist.  The obstruction to the existence of
a global minimal model is encoded in an ideal class, which is trivial
if and only if a global minimal model exists: we provide a function
which returns this class. When the obstruction is not trivial, there
exist models which are minimal at all primes except at a single prime
in the obstruction class, where the discriminant valuation is 12 more
than the minimal valuation at that prime; we provide a function to
return such a model.

The implementation of this functionality is based on work of Kraus
[Kra1989]_ which gives a local condition for when a pair of number field
elements `c_4`, `c_6` belong to a Weierstrass model which is
integral at a prime `P`, together with a global version. Only primes
dividing 2 or 3 are hard to deal with. In order to compute the
corresponding integral model one then needs to combine together the
local transformations implicit in [Kra1989]_ into a single global one.

Various utility functions relating to the testing and use of Kraus's
conditions are included here.

AUTHORS:

- John Cremona (2015)
"""

##############################################################################
#       Copyright (C) 2012-2014 John Cremona <john.cremona@gmail.com>
#                          William Stein <wstein@gmail.com>
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
##############################################################################

from sage.schemes.elliptic_curves.constructor import EllipticCurve


def c4c6_nonsingular(c4, c6):
    r"""
    Check if `c_4`, `c_6` are integral with valid associated discriminant.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    OUTPUT:

    Boolean, ``True`` if `c_4`, `c_6` are both integral and `c_4^3-c_6^2` is a
    nonzero multiple of 1728.

    EXAMPLES:

    Over `\QQ`::

        sage: from sage.schemes.elliptic_curves.kraus import c4c6_nonsingular
        sage: c4c6_nonsingular(0,0)
        False
        sage: c4c6_nonsingular(0,1/2)
        False
        sage: c4c6_nonsingular(2,3)
        False
        sage: c4c6_nonsingular(4,8)
        False
        sage: all(c4c6_nonsingular(*E.c_invariants()) for E in cremona_curves([    11..100]))
        True

    Over number fields::

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: c4c6_nonsingular(-217728*a - 679104, 141460992*a + 409826304)
        True
        sage: K.<a> = NumberField(x^3 - 10)
        sage: c4c6_nonsingular(-217728*a - 679104, 141460992*a + 409826304)
        True
    """
    if not (c4.is_integral() and c6.is_integral()):
        return False
    D = (c4**3-c6**2)/1728
    return not D.is_zero() and D.is_integral()


def c4c6_model(c4, c6, assume_nonsingular=False):
    r"""
    Return the elliptic curve [0,0,0,-c4/48,-c6/864] with given c-invariants.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``assume_nonsingular`` -- boolean (default: ``False``); if ``True``,
      check for integrality and nosingularity

    OUTPUT:

    The elliptic curve with a-invariants `[0,0,0,-c_4/48,-c_6/864]`, whose
    c-invariants are the given `c_4`, `c_6`.  If the supplied invariants are
    singular, returns ``None`` when ``assume_nonsingular`` is ``False`` and
    raises an :exc:`ArithmeticError` otherwise.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.kraus import c4c6_model
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^3 - 10)                                             # needs sage.rings.number_field
        sage: c4c6_model(-217728*a - 679104, 141460992*a + 409826304)                   # needs sage.rings.number_field
        Elliptic Curve defined by y^2 = x^3 + (4536*a+14148)*x + (-163728*a-474336)
         over Number Field in a with defining polynomial x^3 - 10

        sage: c4, c6 = EllipticCurve('389a1').c_invariants()
        sage: c4c6_model(c4,c6)
        Elliptic Curve defined by y^2 = x^3 - 7/3*x + 107/108 over Rational Field
    """
    if not assume_nonsingular:
        if not c4c6_nonsingular(c4,c6):
            return None
    return EllipticCurve([0,0,0,-c4/48,-c6/864])

# Arithmetic utility functions


def make_integral(a, P, e):
    r"""
    Return `b` in `O_K` with `P^e|(a-b)`, given `a` in `O_{K,P}`.

    INPUT:

    - ``a`` -- a number field element integral at `P`

    - ``P`` -- a prime ideal of the number field

    - ``e`` -- positive integer

    OUTPUT:

    A globally integral number field element `b` which is congruent to
    `a` modulo `P^e`.

    ALGORITHM:

    Totally naive, we simply test residues modulo `P^e` until one
    works.  We will only use this when `P` is a prime dividing 2 and `e`
    is the ramification degree, so the number of residues to check is
    at worst `2^d` where `d` is the degree of the field.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.kraus import make_integral

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: P = K.primes_above(2)[0]
        sage: e = P.ramification_index(); e
        2
        sage: x = 1/5
        sage: b = make_integral(x, P, e); b
        1
        sage: (b-x).valuation(P) >= e
        True
        sage: make_integral(1/a, P, e)
        Traceback (most recent call last):
        ...
        ArithmeticError: Cannot lift 1/10*a to O_K mod (Fractional ideal (2, a))^2
    """
    for b in (P**e).residues():
        if (a-b).valuation(P) >= e:
            return b
    raise ArithmeticError("Cannot lift %s to O_K mod (%s)^%s" % (a, P, e))


def sqrt_mod_4(x, P):
    r"""
    Return a local square root mod 4, if it exists.

    INPUT:

    - ``x`` -- an integral number field element

    - ``P`` -- a prime ideal of the number field dividing 2

    OUTPUT:

    A pair ``(True, r)`` where that `r^2-x` has valuation at least `2e`,
    or ``(False, 0)`` if there is no such `r`.  Note that
    `r^2\mod{P^{2e}}` only depends on `r\mod{P^e}`.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import sqrt_mod_4
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: P = K.primes_above(2)[0]
        sage: sqrt_mod_4(1 + 2*a, P)
        (False, 0)
        sage: sqrt_mod_4(-1 + 2*a, P)
        (True, a + 1)
        sage: (1+a)^2 - (-1 + 2*a)
        12
        sage: e = P.ramification_index()
        sage: ((1+a)^2 - (-1+2*a)).mod(P**e)
        0
    """
    e = P.ramification_index()
    P2 = P**e
    for r in P2.residues():
        if (r*r-x).valuation(P) >= 2*e:
            return True, r
    return False, 0

# Kraus test and check for primes dividing 3:


def check_b2_local(c4, c6, P, b2, debug=False):
    r"""
    Test if `b_2` gives a valid model at a prime dividing 3.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``P`` -- a prime ideal of the number field which divides 3

    - ``b2`` -- an element of the number field

    OUTPUT:

    The elliptic curve which is the `(b_2/12,0,0)`-transform of
    `[0,0,0,-c_4/48,-c_6/864]` if this is integral at `P`, else ``False``.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: c4 = -60544*a + 385796
        sage: c6 = -55799680*a + 262126328
        sage: P3a, P3b = K.primes_above(3)
        sage: from sage.schemes.elliptic_curves.kraus import check_b2_local

    `b_2=0` works at the first prime but not the second::

        sage: b2 = 0
        sage: check_b2_local(c4,c6,P3a,b2)                                               # needs sage.rings.number_field
        Elliptic Curve defined by
         y^2 = x^3 + (3784/3*a-96449/12)*x + (1743740/27*a-32765791/108)
         over Number Field in a with defining polynomial x^2 - 10
        sage: check_b2_local(c4,c6,P3b,b2)                                               # needs sage.rings.number_field
        False

    `b_2=-a` works at the second prime but not the first::

        sage: b2 = -a                                                                   # needs sage.rings.number_field
        sage: check_b2_local(c4,c6,P3a,b2,debug=True)                                    # needs sage.rings.number_field
        check_b2_local: not integral at Fractional ideal (3, a + 1)
        False
        sage: check_b2_local(c4,c6,P3b,b2)                                               # needs sage.rings.number_field
        Elliptic Curve defined by
         y^2 = x^3 + (-1/4*a)*x^2 + (3784/3*a-192893/24)*x + (56378369/864*a-32879311/108)
         over Number Field in a with defining polynomial x^2 - 10

    Using CRT we can do both with the same `b_2`::

        sage: b2 = K.solve_CRT([0,-a],[P3a,P3b]); b2                                    # needs sage.rings.number_field
        a + 1
        sage: check_b2_local(c4,c6,P3a,b2)                                               # needs sage.rings.number_field
        Elliptic Curve defined by
         y^2 = x^3 + (1/4*a+1/4)*x^2 + (10091/8*a-128595/16)*x + (4097171/64*a-19392359/64)
         over Number Field in a with defining polynomial x^2 - 10
        sage: check_b2_local(c4,c6,P3b,b2)                                               # needs sage.rings.number_field
        Elliptic Curve defined
         by y^2 = x^3 + (1/4*a+1/4)*x^2 + (10091/8*a-128595/16)*x + (4097171/64*a-19392359/64)
         over Number Field in a with defining polynomial x^2 - 10
    """
    E = c4c6_model(c4,c6).rst_transform(b2/12,0,0)
    if not (c4,c6) == E.c_invariants():
        if debug:
            print("check_b2_local: wrong c-invariants at P=%s" % P)
        return False
    if not E.is_local_integral_model(P):
        if debug:
            print("check_b2_local: not integral at %s" % P)
        return False
    return E


def check_b2_global(c4, c6, b2, debug=False):
    r"""
    Test if `b_2` gives a valid model at all primes dividing 3.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``b2`` -- an element of the number field

    OUTPUT:

    The elliptic curve which is the `(b_2/12,0,0)`-transform of
    `[0,0,0,-c_4/48,-c_6/864]` if this is integral at all primes `P`
    dividing 3, else ``False``.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: c4 = -60544*a + 385796
        sage: c6 = -55799680*a + 262126328
        sage: b2 = a+1
        sage: from sage.schemes.elliptic_curves.kraus import check_b2_global
        sage: check_b2_global(c4,c6,b2)
        Elliptic Curve defined by
         y^2 = x^3 + (1/4*a+1/4)*x^2 + (10091/8*a-128595/16)*x + (4097171/64*a-19392359/64)
         over Number Field in a with defining polynomial x^2 - 10
        sage: check_b2_global(c4,c6,0,debug=True)
        check_b2_global: not integral at all primes dividing 3
        False
        sage: check_b2_global(c4,c6,-a,debug=True)
        check_b2_global: not integral at all primes dividing 3
        False
    """
    E = c4c6_model(c4,c6).rst_transform(b2/12,0,0)
    if not (c4,c6) == E.c_invariants():
        if debug:
            print("check_b2_global: wrong c-invariants")
        return False
    if not all(E.is_local_integral_model(P)
               for P in c4.parent().primes_above(3)):
        if debug:
            print("check_b2_global: not integral at all primes dividing 3")
        return False
    return E


def check_Kraus_local_3(c4, c6, P, assume_nonsingular=False, debug=False):
    r"""
    Test if `c_4`, `c_6` satisfy Kraus's conditions at a prime `P` dividing 3.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``P`` -- a prime ideal of the number field which divides 3

    - ``assume_nonsingular`` -- boolean (default: ``False``); if ``True``,
      check for integrality and nosingularity

    OUTPUT:

    Either ``(False, 0)`` if Kraus's conditions fail, or ``(True, b2)`` if
    they pass, in which case the elliptic curve which is the
    `(b_2/12,0,0)`-transform of `[0,0,0,-c_4/48,-c_6/864]` is integral at `P`.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import check_Kraus_local_3
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: c4 = -60544*a + 385796
        sage: c6 = -55799680*a + 262126328
        sage: P3a, P3b = K.primes_above(3)
        sage: check_Kraus_local_3(c4,c6,P3a)
        (True, 0)
        sage: check_Kraus_local_3(c4,c6,P3b)
        (True, -a)

    An example in a field where 3 is ramified::

        sage: # needs sage.rings.number_field
        sage: K.<a> = NumberField(x^2 - 15)
        sage: c4 = -60504*a + 386001
        sage: c6 = -55346820*a + 261045153
        sage: P3 = K.primes_above(3)[0]
        sage: check_Kraus_local_3(c4,c6,P3)
        (True, a)
    """
    if not assume_nonsingular:
        if not c4c6_nonsingular(c4,c6):
            return False, 0
    e = P.ramification_index()
    P3 = P**e
    if c4.valuation(P) == 0:
        b2 = (-c6*c4.inverse_mod(P3)).mod(P3)
        if debug:
            assert check_b2_local(c4,c6,P,b2)
        return True, b2
    if c6.valuation(P) >= 3*e:
        b2 = c6.parent().zero()
        if debug:
            assert check_b2_local(c4,c6,P,b2)
        return True, b2
    # check for a solution x to x^3-3*x*c4-26=0 (27), such an x must
    # also satisfy x*c4+c6=0 (3) and x^2=c4 (3) and x^3=-c6 (9), and
    # if x is a solution then so is any x'=x (3) so it is enough to
    # check residues mod 3.
    for x in P3.residues():
        if (x*c4+c6).valuation(P) >= e:
            if (x*(x*x-3*c4)-2*c6).valuation(P) >= 3*e:
                if debug:
                    assert check_b2_local(c4,c6,P,x)
                return True, x
    return False, 0

# Kraus test and check for primes dividing 2:


def check_a1a3_local(c4, c6, P, a1, a3, debug=False):
    r"""
    Test if `a_1`, `a_3` are valid at a prime `P` dividing `2`.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``P`` -- a prime ideal of the number field which divides 2

    - ``a1``, ``a3`` -- elements of the number field

    OUTPUT:

    The elliptic curve which is the `(a_1^2/12,a_1/2,a_3/2)`-transform of
    `[0,0,0,-c_4/48,-c_6/864]` if this is integral at `P`, else ``False``.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import check_a1a3_local
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: c4 = -60544*a + 385796
        sage: c6 = -55799680*a + 262126328
        sage: P = K.primes_above(2)[0]
        sage: check_a1a3_local(c4,c6,P,a,0)
        Elliptic Curve defined by
         y^2 + a*x*y = x^3 + (3784/3*a-24106/3)*x + (1772120/27*a-2790758/9)
         over Number Field in a with defining polynomial x^2 - 10
        sage: check_a1a3_local(c4,c6,P,a,a,debug=True)
        check_a1a3_local: not integral at Fractional ideal (2, a)
        False
    """
    E = c4c6_model(c4,c6).rst_transform(a1**2/12,a1/2,a3/2)
    if not (c4,c6) == E.c_invariants():
        if debug:
            print("check_a1a3_local: wrong c-invariants at P=%s" % P)
        return False
    if not E.is_local_integral_model(P):
        if debug:
            print("check_a1a3_local: not integral at %s" % P)
        return False
    return E


def check_a1a3_global(c4, c6, a1, a3, debug=False):
    r"""
    Test if `a_1`, `a_3` are valid at all primes `P` dividing 2.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``a1``, ``a3`` -- elements of the number field

    OUTPUT:

    The elliptic curve which is the `(a_1^2/12,a_1/2,a_3/2)`-transform of
    `[0,0,0,-c_4/48,-c_6/864]` if this is integral at all primes `P`
    dividing 2, else ``False``.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import check_a1a3_global
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: c4 = -60544*a + 385796
        sage: c6 = -55799680*a + 262126328
        sage: check_a1a3_global(c4,c6,a,a,debug=False)
        False
        sage: check_a1a3_global(c4,c6,a,0)
        Elliptic Curve defined by
         y^2 + a*x*y = x^3 + (3784/3*a-24106/3)*x + (1772120/27*a-2790758/9)
         over Number Field in a with defining polynomial x^2 - 10
    """
    E = c4c6_model(c4,c6).rst_transform(a1**2/12,a1/2,a3/2)
    if not (c4, c6) == E.c_invariants():
        if debug:
            print("wrong c-invariants")
        return False
    if not all(E.is_local_integral_model(P)
               for P in c4.parent().primes_above(2)):
        if debug:
            print("not integral at all primes above 2")
        return False
    return E


def check_rst_global(c4, c6, r, s, t, debug=False):
    r"""
    Test if the `(r,s,t)`-transform of the standard `c_4,c_6`-model is integral.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``r``, ``s``, ``t`` -- elements of the number field

    OUTPUT:

    The elliptic curve which is the `(r,s,t)`-transform of
    `[0,0,0,-c_4/48,-c_6/864]` if this is integral at all primes `P`, else
    ``False``.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import check_rst_global
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2-10)
        sage: c4 = -60544*a + 385796
        sage: c6 = -55799680*a + 262126328
        sage: check_rst_global(c4,c6,1/3*a - 133/6, 3/2*a, -89/2*a + 5)
        Elliptic Curve defined by
         y^2 + 3*a*x*y + (-89*a+10)*y = x^3 + (a-89)*x^2 + (1202*a-5225)*x + (34881*a-151813)
         over Number Field in a with defining polynomial x^2 - 10
        sage: check_rst_global(c4,c6,a, 3, -89*a, debug=False)
        False
    """
    E = c4c6_model(c4,c6).rst_transform(r,s,t)
    if not (c4,c6) == E.c_invariants():
        if debug:
            print("test_rst_global: wrong c-invariants")
        return False
    if not E.is_global_integral_model():
        if debug:
            print("test_rst_global: not integral at some prime")
            print(E.ainvs())
            K = E.base_field()
            for P in K.primes_above(2)+K.primes_above(3):
                if not E.is_local_integral_model(P):
                    print(" -- not integral at P=%s" % P)
        return False
    return E

# When a1 is None this function finds a pair a1, a3 such that there is
# a model with these invariants and a2=0 with the given c4, c6,
# integral at P.  The value of a1 is unique modulo 2 (i.e. mod P^e
# where e is the ramification degree of 2); the value of a3 is unique
# mod 2 once a1 is fixed, but not otherwise: when a1 is replaced by
# a1+2s we must replace a3 by a3 + a1*s*(a1+s).

# Because of the latter point, and to fix the bug at #19665, we allow
# the user to specify a1, in which case a3 is computed from it.  This
# is important in the global application where we need to put together
# the transforms for all the primes above 2: we must first get the
# local a1's, then CRT these to get a global a1, then go back to get
# the local a3's and finally CRT these.


def check_Kraus_local_2(c4, c6, P, a1=None, assume_nonsingular=False):
    r"""
    Test if `c_4`, `c_6` satisfy Kraus's conditions at a prime `P` dividing 2.

    INPUT:

    - ``c4``, ``c6`` -- integral elements of a number field

    - ``P`` -- a prime ideal of the number field which divides 2

    - ``a1`` -- an integral elements of a number field, or ``None`` (default)

    - ``assume_nonsingular`` -- boolean (default: ``False``); if ``True``,
      check for integrality and nonsingularity

    OUTPUT:

    Either ``(False, 0, 0)`` if Kraus's conditions fail, or ``(True, a1,
    a3)`` if they pass, in which case the elliptic curve which is the
    `(a_1^2/12,a_1/2,a_3/2)`-transform of `[0,0,0,-c_4/48,-c_6/864]` is
    integral at `P`.  If `a_1` is provided and valid then the output will
    be ``(True, a1, a3)`` for suitable `a_3`.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import check_Kraus_local_2
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: c4 = -60544*a + 385796  # EllipticCurve([a,a,0,1263*a-8032,62956*a-305877])
        sage: c6 = -55799680*a + 262126328
        sage: P = K.primes_above(2)[0]
        sage: check_Kraus_local_2(c4,c6,P)
        (True, a, 0)
    """
    if not assume_nonsingular:
        if not c4c6_nonsingular(c4,c6):
            return False,0,0
    e = P.ramification_index()
    P2 = P**e
    c4val = c4.valuation(P)

    if c4val == 0:
        if a1 is None:
            flag, t = sqrt_mod_4(-c6,P)
            if not flag:
                return False,0,0
            # In the assignment to a1, a3 we divide by units at P,
            # (note that c6+a1**6 = 0 mod P**e so dividing by 4 is OK)
            # but the results, which are well-defined modulo P^e, may
            # not be globally integral
            a1 = make_integral(c4/t,P,e)
        a13 = a1**3
        a3 = make_integral((c6+a13**2)/(4*a13),P,2*e)
        if check_a1a3_local(c4,c6,P,a1,a3):
            return True, a1,a3
        else:
            raise RuntimeError("check_Kraus_local_2 fails")

    if c4val >= 4*e:
        if a1 is None:
            a1 = c4.parent().zero() # 0
        flag, a3 = sqrt_mod_4(c6/8,P)
        if flag:
            if check_a1a3_local(c4,c6,P,a1,a3):
                return True, a1,a3
            else:
                raise RuntimeError("check_Kraus_local_2 fails")
        else:
            return False,0,0

    # val(c4) strictly between 0 and 4e; a1 unique mod 2, with 3 conditions to be satisfied:

    P2res = [a1] if a1 else P2.residues()
    for a1 in P2res:
        Px = -a1**6+3*a1**2*c4+2*c6
        if Px.valuation(P) >= 4*e:                                   # (i)
            flag, a3 = sqrt_mod_4(Px/16,P)                           # (ii)
            if flag:
                a1sq = a1*a1
                if (4*a1sq*Px-(a1sq**2-c4)**2).valuation(P) >= 8*e:  # (iii)
                    if check_a1a3_local(c4,c6,P,a1,a3):
                        return True, a1, a3
                    else:
                        raise RuntimeError("check_Kraus_local_2 fails")
    # end of loop, but no a1 found
    return False, 0, 0

# Wrapper function for local Kraus check, outsources the real work to
# other functions for primes dividing 2 or 3:


def check_Kraus_local(c4, c6, P, assume_nonsingular=False):
    r"""
    Check Kraus's conditions locally at a prime `P`.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``P`` -- a prime ideal of the number field

    - ``assume_nonsingular`` -- boolean (default: ``False``); if ``True``,
      check for integrality and nonsingularity

    OUTPUT:

    Tuple: either ``(True, E)`` if there is a Weierstrass model `E` integral
    at `P` and with invariants `c_4`, `c_6`, or ``(False, None)`` if there is
    none.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import check_Kraus_local
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 15)
        sage: P2 = K.primes_above(2)[0]
        sage: P3 = K.primes_above(3)[0]
        sage: P5 = K.primes_above(5)[0]
        sage: E = EllipticCurve([a,a,0,1263*a-8032,62956*a-305877])
        sage: c4, c6 = E.c_invariants()
        sage: flag, E = check_Kraus_local(c4,c6,P2); flag
        True
        sage: E.is_local_integral_model(P2) and (c4,c6)==E.c_invariants()
        True
        sage: flag, E = check_Kraus_local(c4,c6,P3); flag
        True
        sage: E.is_local_integral_model(P3) and (c4,c6)==E.c_invariants()
        True
        sage: flag, E = check_Kraus_local(c4,c6,P5); flag
        True
        sage: E.is_local_integral_model(P5) and (c4,c6)==E.c_invariants()
        True

        sage: # needs sage.rings.number_field
        sage: c4 = 123+456*a
        sage: c6 = 789+101112*a
        sage: check_Kraus_local(c4,c6,P2)
        (False, None)
        sage: check_Kraus_local(c4,c6,P3)
        (False, None)
        sage: check_Kraus_local(c4,c6,P5)
        (False, None)
    """
    if not assume_nonsingular:
        if not c4c6_nonsingular(c4,c6):
            return False, None
    K = c4.parent()
    if K(2).valuation(P) > 0:
        flag, a1, a3 = check_Kraus_local_2(c4,c6,P,None,True)
        if flag:
            E = check_a1a3_local(c4,c6,P,a1,a3)
            if E:
                return (True, E)
        return (False, None)
    if K(3).valuation(P) > 0:
        flag, b2 = check_Kraus_local_3(c4,c6,P,True)
        if flag:
            E = check_b2_local(c4,c6,P,b2)
            if E:
                return (True, E)
        return (False, None)
    return (True, c4c6_model(c4,c6))


def check_Kraus_global(c4, c6, assume_nonsingular=False, debug=False):
    r"""
    Test if `c_4`, `c_6` satisfy Kraus's conditions at all primes.

    INPUT:

    - ``c4``, ``c6`` -- elements of a number field

    - ``assume_nonsingular`` -- boolean (default: ``False``); if ``True``,
      check for integrality and nonsingularity

    OUTPUT:

    Either ``False`` if Kraus's conditions fail, or, if they pass, an
    elliptic curve `E` which is integral and has c-invariants `c_4`, `c_6`.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.schemes.elliptic_curves.kraus import check_Kraus_global
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: E = EllipticCurve([a,a,0,1263*a-8032,62956*a-305877])
        sage: c4, c6 = E.c_invariants()
        sage: check_Kraus_global(c4,c6,debug=True)
        Local Kraus conditions for (c4,c6)=(-60544*a + 385796,-55799680*a + 262126328) pass at all primes dividing 3
        Using b2=a + 1 gives a model integral at 3:
        (0, 1/4*a + 1/4, 0, 10091/8*a - 128595/16, 4097171/64*a - 19392359/64)
        Local Kraus conditions for (c4,c6)=(-60544*a + 385796,-55799680*a + 262126328) pass at all primes dividing 2
        Using (a1,a3)=(3*a,0) gives a model integral at 2:
        (3*a, 0, 0, 3784/3*a - 23606/3, 1999160/27*a - 9807634/27)
        (a1, b2, a3) = (3*a, a + 1, 0)
        Using (r, s, t)=(1/3*a - 133/6, 3/2*a, -89/2*a + 5) should give a global integral model...
        ...and it does!
        Elliptic Curve defined by
         y^2 + 3*a*x*y + (-89*a+10)*y = x^3 + (a-89)*x^2 + (1202*a-5225)*x + (34881*a-151813)
         over Number Field in a with defining polynomial x^2 - 10

        sage: # needs sage.rings.number_field
        sage: K.<a> = NumberField(x^2 - 15)
        sage: E = EllipticCurve([0, 0, 0, 4536*a + 14148, -163728*a - 474336])
        sage: c4, c6 = E.c_invariants()
        sage: check_Kraus_global(c4,c6)
        Elliptic Curve defined by y^2 = x^3 + (4536*a+14148)*x + (-163728*a-474336)
         over Number Field in a with defining polynomial x^2 - 15

    TESTS (see :issue:`17295`)::

        sage: # needs sage.rings.number_field
        sage: K.<a> = NumberField(x^3 - 7*x - 5)
        sage: E = EllipticCurve([a, 0, 1, 2*a^2 + 5*a + 3, -a^2 - 3*a - 2])
        sage: assert E.conductor().norm() == 8
        sage: G = K.galois_group(names='b')
        sage: def conj_curve(E, sigma): return EllipticCurve([sigma(a) for a in E.ainvs()])
        sage: EL = conj_curve(E,G[0])
        sage: L = EL.base_field()
        sage: assert L.class_number() == 2
        sage: EL.isogeny_class()  # long time (~10s)
        Isogeny class of Elliptic Curve defined by
         y^2 + (-1/90*b^4+7/18*b^2-1/2*b-98/45)*x*y + y = x^3 + (1/45*b^5-1/18*b^4-7/9*b^3+41/18*b^2+167/90*b-29/9)*x + (-1/90*b^5+1/30*b^4+7/18*b^3-4/3*b^2-61/90*b+11/5)
         over Number Field in b with defining polynomial x^6 - 42*x^4 + 441*x^2 - 697
    """
    if not assume_nonsingular:
        if not c4c6_nonsingular(c4,c6):
            return False

    # Check all primes dividing 3; for each get the value of b2
    K = c4.parent()
    three = K.ideal(3)
    Plist3 = K.primes_above(3)
    dat = [check_Kraus_local_3(c4,c6,P,True) for P in Plist3]
    if not all(d[0] for d in dat):
        if debug:
            print("Local Kraus condition for (c4,c6)=(%s,%s) fails at some prime dividing 3" % (c4,c6))
        return False
    if debug:
        print("Local Kraus conditions for (c4,c6)=(%s,%s) pass at all primes dividing 3" % (c4,c6))

    # OK at all primes dividing 3; now use CRT to combine the b2
    # values to get a single residue class for b2 mod 3:

    b2list = [d[1] for d in dat]
    P3list = [P**three.valuation(P) for P in Plist3]
    b2 = K.solve_CRT(b2list,P3list, check=True).mod(three)

    # test that this b2 value works at all P|3:
    if debug:
        E = check_b2_global(c4,c6,b2)
        if E:
            print("Using b2=%s gives a model integral at 3:\n%s" % (b2,E.ainvs()))
        else:
            raise RuntimeError("Error in check_Kraus_global at some prime dividing 3")

    # Check all primes dividing 2; for each get the value of a1, then
    # CRT these to get a single a1 (mod 2) and use these to obtain
    # local a3; finally CRT these
    two = K.ideal(2)
    Plist2 = K.primes_above(2)
    dat = [check_Kraus_local_2(c4, c6, P, None, True) for P in Plist2]
    if not all(d[0] for d in dat):
        if debug:
            print("Local Kraus condition for (c4,c6)=(%s,%s) fails at some prime dividing 2" % (c4, c6))
        return False
    if debug:
        print("Local Kraus conditions for (c4,c6)=(%s,%s) pass at all primes dividing 2" % (c4,c6))

    # OK at all primes dividing 2; now use CRT to combine the a1
    # values to get the residue classes of a1 mod 2:
    P2list = [P**(two.valuation(P)) for P in Plist2]
    a1list = [d[1] for d in dat]
    a1 = K.solve_CRT(a1list,P2list, check=True)
    # See comment below: this is needed for when we combine with the primes above 3.
    if a1 not in three:  # three.divides(a1) causes a segfault
        a1 = 3*a1

    # Using this a1, recompute the local a3's:
    dat = [check_Kraus_local_2(c4,c6,P,a1,True) for P in Plist2]
    # Use CRT to combine these:
    a3list = [d[2] for d in dat]
    a3 = K.solve_CRT(a3list,P2list, check=True)

    # test that these a1,a3 values work at all P|2:
    if debug:
        E = check_a1a3_global(c4,c6,a1,a3,debug)
        if E:
            print("Using (a1,a3)=(%s,%s) gives a model integral at 2:\n%s" % (a1,a3,E.ainvs()))
        else:
            raise RuntimeError("Error in check_Kraus_global at some prime dividing 2")

    # Now we put together the 2-adic and 3-adic transforms to get a
    # global (r,s,t)-transform from [0,0,0,-c4/48,-c6/864] to a global
    # integral model.

    # We need the combined transform (r,s,t) to have both the forms
    # (r,s,t) = (a1^2/12,a1/2,a3/2)*(r2,0,0) with 2-integral r2, and
    # (r,s,t) = (b2/12,0,0,0)*(r3,s3,t3) with 3-integral r3,s3,t3.

    # A solution (assuming that a1,a3,b2 are globally integral) is
    # r2=(b2-a1^2)/3, r3=(b2-a1^2)/4, s3=a1/2, t3=(a1*r2+a3)/2,
    # provided that a1 =0 (mod 3), to make t3 3-integral.  Since a1
    # was only determined mod 2 this can be fixed first, simply by
    # multiplying a1 by 3 if necessary.  We did this above.

    if debug:
        print("(a1, b2, a3) = (%s, %s, %s)" % (a1,b2,a3))
    assert a1.is_integral()
    assert a3.is_integral()
    assert b2.is_integral()
    s = a1/2
    r = b2/3 - s**2
    t = s*(b2-a1**2)/3 + a3/2
    if debug:
        print("Using (r, s, t)=(%s, %s, %s) should give a global integral model..." % (r,s,t))

    # Final computation of the curve E:
    E = check_rst_global(c4,c6,r,s,t,debug)
    if not E:
        if debug:
            print("Error in check_Kraus_global with combining mod-2 and mod-3 transforms")
            E = c4c6_model(c4,c6).rst_transform(r,s,t)
            print("Transformed model is %a" % (E.ainvs(),))
            for P in Plist2+Plist3:
                if not E.is_local_integral_model(P):
                    print("Not integral at P=%s" % P)
        raise RuntimeError("Error in check_Kraus_global combining transforms at 2 and 3")

    # Success!
    if debug:
        print("...and it does!")
    return E


def semi_global_minimal_model(E, debug=False):
    r"""
    Return a global minimal model for this elliptic curve if it
    exists, or a model minimal at all but one prime otherwise.

    INPUT:

    - ``E`` -- an elliptic curve over a number field

    - ``debug`` -- boolean (default: ``False``); if ``True``, prints some
      messages about the progress of the computation

    OUTPUT:

    A tuple ``(Emin, I)`` where ``Emin`` is an elliptic curve which is either a
    global minimal model of `E` if one exists (i.e., an integral model
    which is minimal at every prime), or a semi-global minimal model
    (i.e., an integral model which is minimal at every prime except
    one).  `I` is the unit ideal of ``Emin`` is a global minimal model, else
    is the unique prime at which ``Emin`` is not minimal.  Thus in all
    cases,
    ``Emin.minimal_discriminant_ideal() * I**12 == (E.discriminant())``.

    .. NOTE::

        This function is normally not called directly by users, who
        will use the elliptic curve method :meth:`global_minimal_model`
        instead; that method also applied various reductions after
        minimising the model.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 - 10)
        sage: K.class_number()
        2
        sage: E = EllipticCurve([0,0,0,-186408*a - 589491, 78055704*a + 246833838])
        sage: from sage.schemes.elliptic_curves.kraus import semi_global_minimal_model
        sage: Emin, P = semi_global_minimal_model(E)
        sage: Emin
        Elliptic Curve defined by
         y^2 + 3*x*y + (2*a-11)*y = x^3 + (a-10)*x^2 + (-152*a-415)*x + (1911*a+5920)
         over Number Field in a with defining polynomial x^2 - 10
        sage: E.minimal_discriminant_ideal()*P**12 == K.ideal(Emin.discriminant())
        True

    TESTS:

    Consider (see :issue:`20737`) a curve with no global minimal model
    whose non-minimality class has order 3 in the class group, which
    has order 3315. The smallest prime in that ideal class has norm
    23567::

        sage: # long time, needs sage.rings.number_field
        sage: K.<a> = NumberField(x^2 - x + 31821453)
        sage: ainvs = (0, 0, 0, -382586771000351226384*a - 2498023791133552294513515,
        ....:          358777608829102441023422458989744*a + 1110881475104109582383304709231832166)
        sage: E = EllipticCurve(ainvs)
        sage: from sage.schemes.elliptic_curves.kraus import semi_global_minimal_model
        sage: Emin, p = semi_global_minimal_model(E)  # 25s
        sage: p
        Fractional ideal (23567, a + 2270)
        sage: p.norm()
        23567
        sage: Emin.discriminant().norm().factor()
        23567^12
    """
    c = E.global_minimality_class()
    I = c.ideal()
    c4, c6 = E.c_invariants()
    if c.is_one():
        P = E.base_field().ideal(1)
    else:
        if debug:
            print("No global minimal model, obstruction class = %s of order %s" % (c, c.order()))
        bound = E.base_field().minkowski_bound().round() * 5
        have_prime = False
        while not have_prime:
            try:
                P = c.representative_prime(norm_bound=bound)
                have_prime = True
            except RuntimeError:
                bound *= 2
        if debug:
            print("Using prime {} of norm {} in that class".format(P, P.norm()))
        I = I / P
    u = I.gens_reduced()[0]
    rc4 = c4 / u**4
    rc6 = c6 / u**6
    Emin = check_Kraus_global(rc4, rc6, assume_nonsingular=True, debug=debug)
    if Emin:
        return Emin, P
    raise RuntimeError("failed to compute global minimal model")
