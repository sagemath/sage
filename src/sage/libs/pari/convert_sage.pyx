# sage_setup: distribution = sagemath-pari
r"""
Convert PARI objects to Sage types
"""

# ****************************************************************************
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#       Copyright (C) 2016 Luca De Feo <luca.defeo@polytechnique.edu>
#       Copyright (C) 2016 Vincent Delecroix <vincent.delecroix@u-bordeaux.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from cypari2.types cimport (GEN, typ, t_INT, t_FRAC, t_REAL, t_COMPLEX,
                            t_INTMOD, t_PADIC, t_INFINITY, t_VEC, t_COL,
                            t_VECSMALL, t_MAT, t_STR,
                            lg, precp)
from cypari2.paridecl cimport *
from cypari2.stack cimport new_gen
from sage.libs.pari.convert_gmp cimport INT_to_mpz, new_gen_from_mpz_t, new_gen_from_mpq_t, INTFRAC_to_mpq

from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_sgn, mpz_get_ui, mpz_set, mpz_set_si, mpz_set_ui
from sage.libs.gmp.mpq cimport mpq_denref, mpq_numref
from sage.misc.lazy_import import LazyImport
from sage.rings.integer cimport smallInteger
from sage.rings.infinity import Infinity

try:
    from sage.rings.real_mpfr import RealField
    from sage.rings.complex_mpfr import ComplexField
except ImportError:
    pass

try:
    from sage.rings.number_field.number_field import QuadraticField
except ImportError:
    pass
else:
    QQi = QuadraticField(-1, 'i')

try:
    from sage.rings.padics.factory import Qp
except ImportError:
    pass

pari_typ_to_entries_type = LazyImport('sage.libs.pari.convert_sage_matrix', 'pari_typ_to_entries_type')


cpdef gen_to_sage(Gen z, locals=None):
    """
    Convert a PARI gen to a Sage/Python object.

    INPUT:

    - ``z`` -- PARI ``gen``

    - ``locals`` -- (optional) dictionary used in fallback cases that
      involve :func:`sage_eval`

    OUTPUT: one of the following depending on the PARI type of ``z``

    - a :class:`~sage.rings.integer.Integer` if ``z`` is an integer (type ``t_INT``)

    - a :class:`~sage.rings.rational.Rational` if ``z`` is a rational (type ``t_FRAC``)

    - a :class:`~sage.rings.real_mpfr.RealNumber` if ``z`` is a real
      number (type ``t_REAL``). The precision will be equivalent.

    - a :class:`~sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic`
      or a :class:`~sage.rings.complex_mpfr.ComplexNumber` if ``z`` is a complex
      number (type ``t_COMPLEX``). The former is used when the real and imaginary parts are
      integers or rationals and the latter when they are floating point numbers. In that
      case The precision will be the maximal precision of the real and imaginary parts.

    - a Python list if ``z`` is a vector or a list (type ``t_VEC``, ``t_COL``)

    - a Python string if ``z`` is a string (type ``t_STR``)

    - a Python list of Python integers if ``z`` is a small vector (type ``t_VECSMALL``)

    - a matrix if ``z`` is a matrix (type ``t_MAT``)

    - a `p`-adic element (type ``t_PADIC``)

    - a :class:`~sage.rings.infinity.Infinity` if ``z`` is an infinity
      (type ``t_INF``)

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage import gen_to_sage

    Converting an integer::

        sage: z = pari('12'); z
        12
        sage: z.type()
        't_INT'
        sage: a = gen_to_sage(z); a
        12
        sage: a.parent()
        Integer Ring

        sage: gen_to_sage(pari('7^42'))
        311973482284542371301330321821976049

    Converting a rational number::

        sage: z = pari('389/17'); z
        389/17
        sage: z.type()
        't_FRAC'
        sage: a = gen_to_sage(z); a
        389/17
        sage: a.parent()
        Rational Field

        sage: gen_to_sage(pari('5^30 / 3^50'))
        931322574615478515625/717897987691852588770249

    Converting a real number::

        sage: pari.set_real_precision(70)
        15
        sage: z = pari('1.234'); z
        1.234000000000000000000000000000000000000000000000000000000000000000000
        sage: a = gen_to_sage(z); a                                                     # needs sage.rings.real_mpfr
        1.234000000000000000000000000000000000000000000000000000000000000000000000000
        sage: a.parent()                                                                # needs sage.rings.real_mpfr
        Real Field with 256 bits of precision
        sage: pari.set_real_precision(15)
        70
        sage: a = gen_to_sage(pari('1.234')); a                                         # needs sage.rings.real_mpfr
        1.23400000000000000
        sage: a.parent()                                                                # needs sage.rings.real_mpfr
        Real Field with 64 bits of precision

    For complex numbers, the parent depends on the PARI type::

        sage: z = pari('(3+I)'); z
        3 + I
        sage: z.type()
        't_COMPLEX'
        sage: a = gen_to_sage(z); a                                                     # needs sage.rings.number_field
        i + 3
        sage: a.parent()                                                                # needs sage.rings.number_field
        Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        sage: z = pari('(3+I)/2'); z
        3/2 + 1/2*I
        sage: a = gen_to_sage(z); a                                                     # needs sage.rings.number_field
        1/2*i + 3/2
        sage: a.parent()                                                                # needs sage.rings.number_field
        Number Field in i with defining polynomial x^2 + 1 with i = 1*I

        sage: z = pari('1.0 + 2.0*I'); z
        1.00000000000000 + 2.00000000000000*I
        sage: a = gen_to_sage(z); a                                                     # needs sage.rings.real_mpfr
        1.00000000000000000 + 2.00000000000000000*I
        sage: a.parent()                                                                # needs sage.rings.real_mpfr
        Complex Field with 64 bits of precision

        sage: z = pari('1 + 1.0*I'); z
        1 + 1.00000000000000*I
        sage: a = gen_to_sage(z); a                                                     # needs sage.rings.real_mpfr
        1.00000000000000000 + 1.00000000000000000*I
        sage: a.parent()                                                                # needs sage.rings.real_mpfr
        Complex Field with 64 bits of precision

        sage: z = pari('1.0 + 1*I'); z
        1.00000000000000 + I
        sage: a = gen_to_sage(z); a                                                     # needs sage.rings.real_mpfr
        1.00000000000000000 + 1.00000000000000000*I
        sage: a.parent()                                                                # needs sage.rings.real_mpfr
        Complex Field with 64 bits of precision

    Converting polynomials::

        sage: f = pari('(2/3)*x^3 + x - 5/7 + y')
        sage: f.type()
        't_POL'

        sage: R.<x,y> = QQ[]
        sage: gen_to_sage(f, {'x': x, 'y': y})
        2/3*x^3 + x + y - 5/7
        sage: parent(gen_to_sage(f, {'x': x, 'y': y}))
        Multivariate Polynomial Ring in x, y over Rational Field

        sage: # needs sage.symbolic
        sage: x,y = SR.var('x,y')
        sage: gen_to_sage(f, {'x': x, 'y': y})
        2/3*x^3 + x + y - 5/7
        sage: parent(gen_to_sage(f, {'x': x, 'y': y}))
        Symbolic Ring

        sage: gen_to_sage(f)
        Traceback (most recent call last):
        ...
        NameError: name 'x' is not defined

    Converting vectors::

        sage: # needs sage.rings.number_field sage.rings.real_mpfr
        sage: z1 = pari('[-3, 2.1, 1+I]'); z1
        [-3, 2.10000000000000, 1 + I]
        sage: z2 = pari('[1.0*I, [1,2]]~'); z2
        [1.00000000000000*I, [1, 2]]~
        sage: z1.type(), z2.type()
        ('t_VEC', 't_COL')
        sage: a1 = gen_to_sage(z1)
        sage: a2 = gen_to_sage(z2)
        sage: type(a1), type(a2)
        (<... 'list'>, <... 'list'>)
        sage: [parent(b) for b in a1]
        [Integer Ring,
         Real Field with 64 bits of precision,
         Number Field in i with defining polynomial x^2 + 1 with i = 1*I]
        sage: [parent(b) for b in a2]
        [Complex Field with 64 bits of precision, <... 'list'>]

        sage: z = pari('Vecsmall([1,2,3,4])')
        sage: z.type()
        't_VECSMALL'
        sage: a = gen_to_sage(z); a
        [1, 2, 3, 4]
        sage: type(a)
        <... 'list'>
        sage: [parent(b) for b in a]
        [<... 'int'>, <... 'int'>, <... 'int'>, <... 'int'>]

    Matrices::

        sage: z = pari('[1,2;3,4]')
        sage: z.type()
        't_MAT'

        sage: # needs sage.modules
        sage: a = gen_to_sage(z); a
        [1 2]
        [3 4]
        sage: a.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

    Conversion of `p`-adics::

        sage: # needs sage.rings.padics
        sage: z = pari('569 + O(7^8)'); z
        2 + 4*7 + 4*7^2 + 7^3 + O(7^8)
        sage: a = gen_to_sage(z); a
        2 + 4*7 + 4*7^2 + 7^3 + O(7^8)
        sage: a.parent()
        7-adic Field with capped relative precision 8

    Conversion of infinities::

        sage: gen_to_sage(pari('oo'))
        +Infinity
        sage: gen_to_sage(pari('-oo'))
        -Infinity

    Conversion of strings::

        sage: s = pari('"foo"').sage(); s
        'foo'
        sage: type(s)
        <class 'str'>
    """
    cdef GEN g = z.g
    cdef long t = typ(g)
    cdef long tx, ty
    cdef Gen real, imag, prec, xprec, yprec
    cdef Py_ssize_t i, j, nr, nc

    if t == t_INT:
        return Integer(z)
    elif t == t_FRAC:
        return Rational(z)
    elif t == t_REAL:
        prec = z.bitprecision()
        if typ(prec.g) == t_INFINITY:
            sage_prec = 53
        else:
            sage_prec = prec
        return RealField(sage_prec)(z)
    elif t == t_COMPLEX:
        real = z.real()
        imag = z.imag()
        tx = typ(real.g)
        ty = typ(imag.g)
        if tx in [t_INTMOD, t_PADIC] or ty in [t_INTMOD, t_PADIC]:
            raise NotImplementedError("No conversion to python available for t_COMPLEX with t_INTMOD or t_PADIC components")
        if tx == t_REAL or ty == t_REAL:
            xprec = real.bitprecision()  # will be infinite if exact
            yprec = imag.bitprecision()  # will be infinite if exact
            if typ(xprec.g) == t_INFINITY:
                if typ(yprec.g) == t_INFINITY:
                    sage_prec = 53
                else:
                    sage_prec = yprec
            elif typ(yprec.g) == t_INFINITY:
                sage_prec = xprec
            else:
                sage_prec = max(xprec, yprec)

            R = RealField(sage_prec)
            C = ComplexField(sage_prec)
            return C(R(real), R(imag))
        else:
            return QQi([gen_to_sage(real), gen_to_sage(imag)])
    elif t == t_VEC or t == t_COL:
        return [gen_to_sage(x, locals) for x in z.python_list()]
    elif t == t_VECSMALL:
        return z.python_list_small()
    elif t == t_MAT:
        from sage.libs.pari.convert_sage_matrix import gen_to_sage_matrix
        return gen_to_sage_matrix(z, locals)
    elif t == t_PADIC:
        p = z.padicprime()
        K = Qp(Integer(p), precp(g))
        return K(z.lift())
    elif t == t_STR:
        return str(z)
    elif t == t_INFINITY:
        if inf_get_sign(g) >= 0:
            return Infinity
        else:
            return -Infinity

    # Fallback (e.g. polynomials): use string representation
    from sage.misc.sage_eval import sage_eval
    locals = {} if locals is None else locals
    return sage_eval(str(z), locals=locals)


cpdef set_integer_from_gen(Integer self, Gen x):
    r"""
    EXAMPLES::

        sage: [Integer(pari(x)) for x in [1, 2^60, 2., GF(3)(1), GF(9,'a')(2)]]         # needs sage.rings.finite_rings
        [1, 1152921504606846976, 2, 1, 2]
        sage: Integer(pari(2.1)) # indirect doctest
        Traceback (most recent call last):
        ...
        TypeError: Attempt to coerce non-integral real number to an Integer
    """
    # Simplify and lift until we get an integer
    while typ((<Gen>x).g) != t_INT:
        x = x.simplify()
        paritype = typ((<Gen>x).g)
        if paritype == t_INT:
            break
        elif paritype == t_REAL:
            # Check that the fractional part is zero
            if not x.frac().gequal0():
                raise TypeError("Attempt to coerce non-integral real number to an Integer")
            # floor yields an integer
            x = x.floor()
            break
        elif paritype == t_PADIC:
            if x._valp() < 0:
                raise TypeError("cannot convert p-adic with negative valuation to an integer")
            # Lifting a PADIC yields an integer
            x = x.lift()
            break
        elif paritype == t_INTMOD:
            # Lifting an INTMOD yields an integer
            x = x.lift()
            break
        elif paritype == t_POLMOD:
            x = x.lift()
        elif paritype == t_FFELT:
            # x = (f modulo defining polynomial of finite field);
            # we extract f.
            sig_on()
            x = new_gen(FF_to_FpXQ_i((<Gen>x).g))
        else:
            raise TypeError("Unable to coerce PARI %s to an Integer"%x)

    # Now we have a true PARI integer, convert it to Sage
    INT_to_mpz(self.value, (<Gen>x).g)


cpdef Gen new_gen_from_integer(Integer self):
    """
    TESTS::

        sage: Rational(pari(2))  # indirect doctest
        2
        sage: Rational(pari(-1))
        -1
    """
    return new_gen_from_mpz_t(self.value)


cpdef set_rational_from_gen(Rational self, Gen x):
    r"""
    EXAMPLES::

        sage: [Rational(pari(x)) for x in [1, 1/2, 2^60, 2., GF(3)(1), GF(9,'a')(2)]]   # needs sage.rings.finite_rings
        [1, 1/2, 1152921504606846976, 2, 1, 2]
        sage: Rational(pari(2.1)) # indirect doctest
        Traceback (most recent call last):
        ...
        TypeError: Attempt to coerce non-integral real number to an Integer
    """
    x = x.simplify()
    if is_rational_t(typ((<Gen>x).g)):
        INTFRAC_to_mpq(self.value, (<Gen>x).g)
    else:
        a = Integer(x)
        mpz_set(mpq_numref(self.value), a.value)
        mpz_set_si(mpq_denref(self.value), 1)


cpdef Gen new_gen_from_rational(Rational self):
    """
    TESTS::

        sage: Integer(pari(2/2))  # indirect doctest
        1
        sage: Rational(pari(-1/2))
        -1/2
    """
    return new_gen_from_mpq_t(self.value)


cpdef list pari_divisors_small(Integer self):
    r"""
    Return the list of divisors of this number using PARI ``divisorsu``.

    .. SEEALSO::

        This method is better used through :meth:`sage.rings.integer.Integer.divisors`.

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage import pari_divisors_small
        sage: pari_divisors_small(4)
        [1, 2, 4]

    The integer must fit into an unsigned long::

        sage: pari_divisors_small(-4)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: pari_divisors_small(2**65)
        Traceback (most recent call last):
        ...
        AssertionError
    """
    # we need n to fit into a long and not a unsigned long in order to use
    # smallInteger
    assert mpz_fits_slong_p(self.value) and mpz_sgn(self.value) > 0

    cdef unsigned long n = mpz_get_ui(self.value)

    global avma
    cdef pari_sp ltop = avma
    cdef GEN d
    cdef list output

    try:
        sig_on()
        d = divisorsu(n)
        sig_off()
        output = [smallInteger(d[i]) for i in range(1,lg(d))]
        return output
    finally:
        avma = ltop


cpdef pari_is_prime(Integer p):
    r"""
    Return whether ``p`` is a prime.

    The caller must ensure that ``p.value`` fits in a long.

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage import pari_is_prime
        sage: pari_is_prime(2)
        True
        sage: pari_is_prime(3)
        True
        sage: pari_is_prime(1)
        False
        sage: pari_is_prime(4)
        False

    Its recommended to use :meth:`sage.rings.integer.Integer.is_prime`, which checks overflow.
    The following is incorrect, because the number does not fit into a long::

        sage: pari_is_prime(2**64 + 2)
        True
    """
    return bool(uisprime(mpz_get_ui(p.value)))


cpdef pari_is_prime_power(Integer q, bint get_data):
    r"""
    Return whether ``q`` is a prime power.

    The caller must ensure that ``q.value`` fits in a long.

    OUTPUT:

    If ``get_data`` return a tuple of the prime and the exponent.
    Otherwise return a boolean.

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage import pari_is_prime_power
        sage: pari_is_prime_power(2, False)
        True
        sage: pari_is_prime_power(2, True)
        (2, 1)
        sage: pari_is_prime_power(4, False)
        True
        sage: pari_is_prime_power(4, True)
        (2, 2)
        sage: pari_is_prime_power(6, False)
        False
        sage: pari_is_prime_power(6, True)
        (6, 0)

    Its recommended to use :meth:`sage.rings.integer.Integer.is_prime_power`, which checks overflow.
    The following is incorrect, because the number does not fit into a long::

        sage: pari_is_prime_power(2**64 + 2, False)
        True
    """
    cdef long p, n
    n = uisprimepower(mpz_get_ui(q.value), <ulong*>(&p))
    if n:
        return (smallInteger(p), smallInteger(n)) if get_data else True
    else:
        return (q, smallInteger(0)) if get_data else False


cpdef unsigned long pari_maxprime() noexcept:
    """
    Return to which limit PARI has computed the primes.

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage import pari_maxprime
        sage: a = pari_maxprime()
        sage: res = prime_range(2, 2*a)
        sage: b = pari_maxprime()
        sage: b >= 2*a
        True
    """
    return maxprime()


cpdef list pari_prime_range(long c_start, long c_stop, bint py_ints=False):
    """
    Return a list of all primes between ``start`` and ``stop - 1``, inclusive.

    .. SEEALSO::

        :func:`sage.rings.fast_arith.prime_range`

    TESTS::

        sage: from sage.libs.pari.convert_sage import pari_prime_range
        sage: pari_prime_range(2, 19)
        [2, 3, 5, 7, 11, 13, 17]
    """
    cdef long p = 0
    cdef byteptr pari_prime_ptr = diffptr
    res = []
    while p < c_start:
        NEXT_PRIME_VIADIFF(p, pari_prime_ptr)
    while p < c_stop:
        if py_ints:
            res.append(p)
        else:
            z = <Integer>PY_NEW(Integer)
            mpz_set_ui(z.value, p)
            res.append(z)
        NEXT_PRIME_VIADIFF(p, pari_prime_ptr)
    return res
