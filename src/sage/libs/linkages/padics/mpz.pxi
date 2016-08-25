"""
This linkage file implements the padics API using MPIR mpz_t
multiprecision integers.

AUTHORS:

- David Roe (2012-3-1) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2007-2012 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
include "cysignals/signals.pxi"
from cpython.list cimport *

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.pylong cimport mpz_pythonhash
from sage.libs.gmp.rational_reconstruction cimport mpq_rational_reconstruction
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.common_conversion cimport cconv_mpz_t_out_shared, cconv_mpz_t_shared, cconv_mpq_t_out_shared, cconv_mpq_t_shared, cconv_shared
import sage.rings.finite_rings.integer_mod

cdef Integer holder = PY_NEW(Integer)
cdef Integer holder2 = PY_NEW(Integer)

cdef inline int cconstruct(mpz_t value, PowComputer_ prime_pow) except -1:
    """
    Construct a new element.

    INPUT:

    - ``unit`` -- an ``mpz_t`` to be initialized.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_init(value)

cdef inline int cdestruct(mpz_t value, PowComputer_ prime_pow) except -1:
    """
    Deallocate an element.

    INPUT:

    - ``unit`` -- an ``mpz_t`` to be cleared.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_clear(value)

cdef inline int ccmp(mpz_t a, mpz_t b, long prec, bint reduce_a, bint reduce_b, PowComputer_ prime_pow) except -2:
    """
    Comparison of two elements.

    INPUT:

    - ``a`` -- an ``mpz_t``.
    - ``b`` -- an ``mpz_t``.
    - ``prec`` -- a long, the precision of the comparison.
    - ``reduce_a`` -- a bint, whether a needs to be reduced.
    - ``reduce_b`` -- a bint, whether b needs to be reduced.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - If neither a nor be needs to be reduced, returns
      -1 (`a < b`), 0 (`a = b`) or 1 (`a > b`)
    - If at least one needs to be reduced, returns
      0 (``a == b mod p^prec``) or 1 (otherwise)
    """
    cdef int ans
    if reduce_a or reduce_b:
        mpz_sub(holder.value, a, b)
        mpz_mod(holder.value, holder.value, prime_pow.pow_mpz_t_tmp(prec))
        return mpz_sgn(holder.value)
    else:
        ans = mpz_cmp(a,b)
        if ans > 0:
            return 1
        elif ans < 0:
            return -1
        return 0

cdef inline int cneg(mpz_t out, mpz_t a, long prec, PowComputer_ prime_pow) except -1:
    """
    Negation.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the negation.
    - ``a`` -- an ``mpz_t`` to be negated.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_neg(out, a)

cdef inline int cadd(mpz_t out, mpz_t a, mpz_t b, long prec, PowComputer_ prime_pow) except -1:
    """
    Addition.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the sum.
    - ``a`` -- an ``mpz_t``, the first summand.
    - ``b`` -- an ``mpz_t``, the second summand.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_add(out, a, b)

cdef inline bint creduce(mpz_t out, mpz_t a, long prec, PowComputer_ prime_pow) except -1:
    """
    Reduce modulo a power of the maximal ideal.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the reduction.
    - ``a`` -- the element to be reduced.
    - ``prec`` -- a long, the precision to reduce modulo.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if the reduction is zero; False otherwise.
    """
    mpz_mod(out, a, prime_pow.pow_mpz_t_tmp(prec))
    return mpz_sgn(out) == 0

cdef inline bint creduce_small(mpz_t out, mpz_t a, long prec, PowComputer_ prime_pow) except -1:
    """
    Reduce modulo a power of the maximal ideal.

    This function assumes that the input satisfies `-p <= a < 2p`, so
    that it doesn't need any divisions.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the reduction.
    - ``a`` -- the element to be reduced.
    - ``prec`` -- a long, the precision to reduce modulo.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if the reduction is zero; False otherwise.
    """
    if mpz_sgn(a) < 0:
        mpz_add(out, a, prime_pow.pow_mpz_t_tmp(prec))
    elif mpz_cmp(a, prime_pow.pow_mpz_t_tmp(prec)) >= 0:
        mpz_sub(out, a, prime_pow.pow_mpz_t_tmp(prec))
    else:
        mpz_set(out, a)
    return mpz_sgn(out) == 0

cdef inline long cremove(mpz_t out, mpz_t a, long prec, PowComputer_ prime_pow) except -1:
    """
    Extract the maximum power of the uniformizer dividing this
    element.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the unit.
    - ``a`` -- the element whose valuation and unit are desired.
    - ``prec`` -- a long, used if `a = 0`.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - if `a = 0`, returns prec.  Otherwise, returns the number of
      times p divides a.
    """
    if mpz_sgn(a) == 0:
        mpz_set_ui(out, 0)
        return prec
    return mpz_remove(out, a, prime_pow.prime.value)

cdef inline long cvaluation(mpz_t a, long prec, PowComputer_ prime_pow) except -1:
    """
    Returns the maximum power of the uniformizer dividing this
    element.

    This function differs from :meth:`cremove` in that the unit is
    discarded.

    INPUT:

    - ``a`` -- the element whose valuation is desired.
    - ``prec`` -- a long, used if `a = 0`.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - if `a = 0`, returns prec.  Otherwise, returns the number of
      times p divides a.
    """
    if mpz_sgn(a) == 0:
        return prec
    return mpz_remove(holder.value, a, prime_pow.prime.value)

cdef inline bint cisunit(mpz_t a, PowComputer_ prime_pow) except -1:
    """
    Returns whether this element has valuation zero.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a` has valuation 0, and False otherwise.
    """
    return mpz_divisible_p(a, prime_pow.prime.value) == 0

cdef inline int cshift(mpz_t out, mpz_t a, long n, long prec, PowComputer_ prime_pow, bint reduce_afterward) except -1:
    """
    Multiplies by a power of the uniformizer.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the result.  If `n >= 0` then
                 out will be set to `a * p^n`.  If `n < 0`, out will
                 be set to `a // p^n`.
    - ``a`` -- the element to shift.
    - ``n`` -- long, the amount to shift by.
    - ``prec`` -- long, a precision modulo which to reduce.
    - ``prime_pow`` -- the PowComputer for the ring.
    - ``reduce_afterward`` -- whether to reduce afterward.
    """
    if n > 0:
        mpz_mul(out, a, prime_pow.pow_mpz_t_tmp(n))
    elif n < 0:
        sig_on()
        mpz_fdiv_q(out, a, prime_pow.pow_mpz_t_tmp(-n))
        sig_off()
    else: # elif a != out:
        mpz_set(out, a)
    if reduce_afterward:
        creduce(out, out, prec, prime_pow)

cdef inline int cshift_notrunc(mpz_t out, mpz_t a, long n, long prec, PowComputer_ prime_pow) except -1:
    """
    Multiplies by a power of the uniformizer, assuming that the
    valuation of a is at least -n.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the result.  If `n >= 0` then
                 out will be set to `a * p^n`.  If `n < 0`, out will
                 be set to `a // p^n`.
    - ``a`` -- the element to shift.  Assumes that the valuation of a
      is at least -n.
    - ``n`` -- long, the amount to shift by.
    - ``prec`` -- long, a precision modulo which to reduce.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    if n > 0:
        mpz_mul(out, a, prime_pow.pow_mpz_t_tmp(n))
    elif n < 0:
        sig_on()
        mpz_divexact(out, a, prime_pow.pow_mpz_t_tmp(-n))
        sig_off()
    else:
        mpz_set(out, a)

cdef inline int csub(mpz_t out, mpz_t a, mpz_t b, long prec, PowComputer_ prime_pow) except -1:
    """
    Subtraction.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the difference.
    - ``a`` -- an ``mpz_t``, the first input.
    - ``b`` -- an ``mpz_t``, the second input.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_sub(out, a, b)

cdef inline int cinvert(mpz_t out, mpz_t a, long prec, PowComputer_ prime_pow) except -1:
    """
    Inversion.

    The result will be reduced modulo p^prec.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the inverse.
    - ``a`` -- an ``mpz_t``, the element to be inverted.
    - ``prec`` -- a long, the precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    cdef bint success
    success = mpz_invert(out, a, prime_pow.pow_mpz_t_tmp(prec))
    if not success:
        raise ZeroDivisionError
    
cdef inline int cmul(mpz_t out, mpz_t a, mpz_t b, long prec, PowComputer_ prime_pow) except -1:
    """
    Multiplication.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the product.
    - ``a`` -- an ``mpz_t``, the first input.
    - ``b`` -- an ``mpz_t``, the second input.
    - ``prec`` -- a long, the precision: ignored.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_mul(out, a, b)

cdef inline int cdivunit(mpz_t out, mpz_t a, mpz_t b, long prec, PowComputer_ prime_pow) except -1:
    """
    Division.

    The inversion is perfomed modulo p^prec.  Note that no reduction
    is performed after the product.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the quotient.
    - ``a`` -- an ``mpz_t``, the first input.
    - ``b`` -- an ``mpz_t``, the second input.
    - ``prec`` -- a long, the precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    cdef bint success
    success = mpz_invert(out, b, prime_pow.pow_mpz_t_tmp(prec))
    if not success:
        raise ZeroDivisionError
    mpz_mul(out, a, out)

cdef inline int csetone(mpz_t out, PowComputer_ prime_pow) except -1:
    """
    Sets to 1.

    INPUT:

    - ``out`` -- the ``mpz_t`` in which to store 1.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_set_ui(out, 1)
    
cdef inline int csetzero(mpz_t out, PowComputer_ prime_pow) except -1:
    """
    Sets to 0.

    INPUT:

    - ``out`` -- the ``mpz_t`` in which to store 0.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_set_ui(out, 0)
    
cdef inline bint cisone(mpz_t out, PowComputer_ prime_pow) except -1:
    """
    Returns whether this element is equal to 1.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a = 1`, and False otherwise.
    """
    return mpz_cmp_ui(out, 1) == 0
    
cdef inline bint ciszero(mpz_t out, PowComputer_ prime_pow) except -1:
    """
    Returns whether this element is equal to 0.

    INPUT:

    - ``a`` -- the element to test.
    - ``prime_pow`` -- the PowComputer for the ring.

    OUTPUT:

    - returns True if `a = 0`, and False otherwise.
    """
    return mpz_cmp_ui(out, 0) == 0
    
cdef inline int cpow(mpz_t out, mpz_t a, mpz_t n, long prec, PowComputer_ prime_pow) except -1:
    """
    Exponentiation.

    INPUT:

    - ``out`` -- the ``mpz_t`` in which to store the result.
    - ``a`` -- the base.
    - ``n`` -- an ``mpz_t``, the exponent.
    - ``prec`` -- a long, the working absolute precision.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_powm(out, a, n, prime_pow.pow_mpz_t_tmp(prec))

cdef inline int ccopy(mpz_t out, mpz_t a, PowComputer_ prime_pow) except -1:
    """
    Copying.

    INPUT:

    - ``out`` -- the ``mpz_t`` to store the result.
    - ``a`` -- the element to copy.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_set(out, a)

cdef inline cpickle(mpz_t a, PowComputer_ prime_pow):
    """
    Serialization into objects that Sage knows how to pickle.

    INPUT:

    - ``a`` the element to pickle.
    - ``prime_pow`` the PowComputer for the ring.

    OUTPUT:

    - an Integer storing ``a``.
    """
    cdef Integer pic = PY_NEW(Integer)
    mpz_set(pic.value, a)
    return pic

cdef inline int cunpickle(mpz_t out, x, PowComputer_ prime_pow) except -1:
    """
    Reconstruction from the output of meth:`cpickle`.

    INPUT:

    - ``out`` -- the ``mpz_t`` in which to store the result.
    - ``x`` -- the result of `meth`:cpickle.
    - ``prime_pow`` -- the PowComputer for the ring.
    """
    mpz_set(out, (<Integer?>x).value)

cdef inline long chash(mpz_t a, long ordp, long prec, PowComputer_ prime_pow) except -1:
    """
    Hashing.

    INPUT:

    - ``a`` -- an ``mpz_t`` storing the underlying element to hash.
    - ``ordp`` -- a long storing the valuation.
    - ``prec`` -- a long storing the precision.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    # This implementation is for backward compatibility and may be changed in the future
    cdef long n, d
    if ordp == 0:
        return mpz_pythonhash(a)
    elif ordp > 0:
        mpz_mul(holder.value, a, prime_pow.pow_mpz_t_tmp(ordp))
        return mpz_pythonhash(holder.value)
    else:
        n = mpz_pythonhash(a)
        d = mpz_pythonhash(prime_pow.pow_mpz_t_tmp(-ordp))
        if d == 1:
            return n
        n = n ^ d
        if n == -1:
            return -2
        return n

cdef clist(mpz_t a, long prec, bint pos, PowComputer_ prime_pow):
    """
    Returns a list of digits in the series expansion.

    This function is used in printing, and expresses ``a`` as a series
    in the standard uniformizer ``p``.

    INPUT:

    - ``a`` -- an ``mpz_t`` giving the underlying `p`-adic element.
    - ``prec`` -- a precision giving the number of digits desired.
    - ``pos`` -- if True then representatives in 0..(p-1) are used;
                 otherwise the range (-p/2..p/2) is used.
    - ``prime_pow`` -- a PowComputer for the ring.

    OUTPUT:

    - A list of p-adic digits `[a_0, a_1, \ldots]` so that
      `a = a_0 + a_1*p + \cdots` modulo `p^{prec}`.
    """
    cdef mpz_t tmp, halfp
    cdef bint neg
    cdef long curpower
    cdef Integer list_elt
    ans = PyList_New(0)
    mpz_set(holder.value, a)
    if pos:
        curpower = prec
        while mpz_sgn(holder.value) != 0 and curpower >= 0:
            list_elt = PY_NEW(Integer)
            mpz_mod(list_elt.value, holder.value, prime_pow.prime.value)
            mpz_sub(holder.value, holder.value, list_elt.value)
            mpz_divexact(holder.value, holder.value, prime_pow.prime.value)
            PyList_Append(ans, list_elt)
            curpower -= 1
    else:
        neg = False
        curpower = prec
        mpz_fdiv_q_2exp(holder2.value, prime_pow.prime.value, 1)
        while mpz_sgn(holder.value) != 0 and curpower > 0:
            curpower -= 1
            list_elt = PY_NEW(Integer)
            mpz_mod(list_elt.value, holder.value, prime_pow.prime.value)
            if mpz_cmp(list_elt.value, holder2.value) > 0:
                mpz_sub(list_elt.value, list_elt.value, prime_pow.prime.value)
                neg = True
            else:
                neg = False
            mpz_sub(holder.value, holder.value, list_elt.value)
            mpz_divexact(holder.value, holder.value, prime_pow.prime.value)
            if neg:
                if mpz_cmp(holder.value, prime_pow.pow_mpz_t_tmp(curpower)) >= 0:
                    mpz_sub(holder.value, holder.value, prime_pow.pow_mpz_t_tmp(curpower))
            PyList_Append(ans, list_elt)
    return ans

# The element is filled in for zero in the output of clist if necessary.
# It could be [] for some other linkages.
_list_zero = Integer(0)

cdef int cteichmuller(mpz_t out, mpz_t value, long prec, PowComputer_ prime_pow) except -1:
    """
    Teichmuller lifting.

    INPUT:

    - ``out`` -- an ``mpz_t`` which is set to a `p-1` root of unity
                 congruent to `value` mod `p`; or 0 if `a \equiv 0
                 \pmod{p}`.
    - ``value`` -- an ``mpz_t``, the element mod `p` to lift.
    - ``prec`` -- a long, the precision to which to lift.
    - ``prime_pow`` -- the Powcomputer of the ring.
    """
    if mpz_divisible_p(value, prime_pow.prime.value) != 0:
        mpz_set_ui(out, 0)
        return 0
    if prec <= 0:
        raise ValueError
    if mpz_sgn(value) < 0 or mpz_cmp(value, prime_pow.pow_mpz_t_tmp(prec)) >= 0:
        mpz_mod(out, value, prime_pow.pow_mpz_t_tmp(prec))
    else:
        mpz_set(out, value)
    # holder.value = 1 / Mod(1 - p, prime_pow.pow_mpz_t_tmp(prec))
    mpz_sub(holder.value, prime_pow.pow_mpz_t_tmp(prec), prime_pow.prime.value)
    mpz_add_ui(holder.value, holder.value, 1)
    mpz_invert(holder.value, holder.value, prime_pow.pow_mpz_t_tmp(prec))
    # Consider x as Mod(value, prime_pow.pow_mpz_t_tmp(prec))
    # holder2.value = x + holder.value*(x^p - x)
    mpz_powm(holder2.value, out, prime_pow.prime.value, prime_pow.pow_mpz_t_tmp(prec))
    mpz_sub(holder2.value, holder2.value, out)
    mpz_mul(holder2.value, holder2.value, holder.value)
    mpz_add(holder2.value, holder2.value, out)
    mpz_mod(holder2.value, holder2.value, prime_pow.pow_mpz_t_tmp(prec))
    # while x != holder2.value:
    #     x = holder2.value
    #     holder2.value = x + holder.value*(x^p - x)
    while mpz_cmp(out, holder2.value) != 0:
        mpz_set(out, holder2.value)
        mpz_powm(holder2.value, out, prime_pow.prime.value, prime_pow.pow_mpz_t_tmp(prec))
        mpz_sub(holder2.value, holder2.value, out)
        mpz_mul(holder2.value, holder2.value, holder.value)
        mpz_add(holder2.value, holder2.value, out)
        mpz_mod(holder2.value, holder2.value, prime_pow.pow_mpz_t_tmp(prec))

cdef int cconv(mpz_t out, x, long prec, long valshift, PowComputer_ prime_pow) except -2:
    """
    Conversion from other Sage types.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the output.

    - ``x`` -- a Sage element that can be converted to a `p`-adic element.

    - ``prec`` -- a long, giving the precision desired: absolute if
                  `valshift = 0`, relative if `valshift != 0`.

    - ``valshift`` -- the power of the uniformizer to divide by before
      storing the result in ``out``.

    - ``prime_pow`` -- a PowComputer for the ring.
    """
    return cconv_shared(out, x, prec, valshift, prime_pow)

cdef inline long cconv_mpq_t(mpz_t out, mpq_t x, long prec, bint absolute, PowComputer_ prime_pow) except? -10000:
    """
    A fast pathway for conversion of rationals that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the output.
    - ``x`` -- an ``mpq_t`` giving the integer to be converted.
    - ``prec`` -- a long, giving the precision desired: absolute or
      relative depending on the ``absolute`` input.
    - ``absolute`` -- if False then extracts the valuation and returns
                      it, storing the unit in ``out``; if True then
                      just reduces ``x`` modulo the precision.
    - ``prime_pow`` -- a PowComputer for the ring.

    OUTPUT:

    - If ``absolute`` is False then returns the valuation that was
      extracted (``maxordp`` when `x = 0`).
    """
    return cconv_mpq_t_shared(out, x, prec, absolute, prime_pow)

cdef inline int cconv_mpq_t_out(mpq_t out, mpz_t x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    """
    Converts the underlying `p`-adic element into a rational

    - ``out`` -- gives a rational approximating the input.  Currently uses rational reconstruction but
                 may change in the future to use a more naive method
    - ``x`` -- an ``mpz_t`` giving the underlying `p`-adic element
    - ``valshift`` -- a long giving the power of `p` to shift `x` by
    -` ``prec`` -- a long, the precision of ``x``, used in rational reconstruction
    - ``prime_pow`` -- a PowComputer for the ring
    """
    return cconv_mpq_t_out_shared(out, x, valshift, prec, prime_pow)

cdef inline long cconv_mpz_t(mpz_t out, mpz_t x, long prec, bint absolute, PowComputer_ prime_pow) except -2:
    """
    A fast pathway for conversion of integers that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``mpz_t`` to store the output.
    - ``x`` -- an ``mpz_t`` giving the integer to be converted.
    - ``prec`` -- a long, giving the precision desired: absolute or
                  relative depending on the ``absolute`` input.
    - ``absolute`` -- if False then extracts the valuation and returns
                      it, storing the unit in ``out``; if True then
                      just reduces ``x`` modulo the precision.
    - ``prime_pow`` -- a PowComputer for the ring.

    OUTPUT:

    - If ``absolute`` is False then returns the valuation that was
      extracted (``maxordp`` when `x = 0`).
    """
    return cconv_mpz_t_shared(out, x, prec, absolute, prime_pow)

cdef inline int cconv_mpz_t_out(mpz_t out, mpz_t x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    """
    Converts the underlying `p`-adic element into an integer if
    possible.

    - ``out`` -- stores the resulting integer as an integer between 0
      and `p^{prec + valshift}`.
    - ``x`` -- an ``mpz_t`` giving the underlying `p`-adic element.
    - ``valshift`` -- a long giving the power of `p` to shift `x` by.
    -` ``prec`` -- a long, the precision of ``x``: currently not used.
    - ``prime_pow`` -- a PowComputer for the ring.
    """
    return cconv_mpz_t_out_shared(out, x, valshift, prec, prime_pow)
