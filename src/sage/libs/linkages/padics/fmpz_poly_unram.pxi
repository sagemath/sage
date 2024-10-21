# sage_setup: distribution = sagemath-pari
"""
This linkage file implements the API for unramified extensions of the padics
using FLINT's fmpz_poly_t.

AUTHORS:

- David Roe, Julian Rueth (2013-03-21) -- initial version
"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from cpython.list cimport PyList_Check, PyList_New, PyList_Append

from sage.ext.stdsage cimport PY_NEW

from sage.rings.padics.common_conversion cimport cconv_mpz_t_out_shared, cconv_mpz_t_shared, cconv_mpq_t_out_shared, cconv_mpq_t_shared, cconv_shared

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.finite_rings.integer_mod cimport IntegerMod_abstract
from sage.rings.finite_rings.integer_mod_ring import Zmod

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpz_poly_sage cimport *


DEF CELEMENT_IS_PY_OBJECT = False

cdef inline int cconstruct(celement value, PowComputer_ prime_pow) except -1:
    """
    Construct a new element.

    INPUT:

    - ``unit`` -- an ``celement`` to be initialized
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_init(value)

cdef inline int cdestruct(celement value, PowComputer_ prime_pow) except -1:
    """
    Deallocate an element.

    INPUT:

    - ``unit`` -- an ``celement`` to be cleared
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_clear(value)

cdef inline int ccmp(celement a, celement b, long prec, bint reduce_a, bint reduce_b, PowComputer_ prime_pow) except -2:
    """
    Comparison of two elements.

    INPUT:

    - ``a`` -- an ``celement``
    - ``b`` -- an ``celement``
    - ``prec`` -- a long, the precision of the comparison
    - ``reduce_a`` -- a bint, whether ``a`` needs to be reduced
    - ``reduce_b`` -- a bint, whether ``b`` needs to be reduced
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT:

    - If neither ``a`` nor ``b`` needs to be reduced, returns
      -1 (if `a < b`), 0 (if `a == b`) or 1 (if `a > b`)

    - If at least one needs to be reduced, returns
      0 (if ``a == b mod p^prec``) or 1 (otherwise)
    """
    csub(prime_pow.poly_ccmp, a, b, prec, prime_pow)
    creduce(prime_pow.poly_ccmp, prime_pow.poly_ccmp, prec, prime_pow)

    if reduce_a or reduce_b:
        return not ciszero(prime_pow.poly_ccmp, prime_pow)

    if prec == 0:
        return 0

    if ciszero(prime_pow.poly_ccmp, prime_pow):
        return 0

    cdef long da = fmpz_poly_degree(a)
    cdef long db = fmpz_poly_degree(b)
    if da < db:
        return -1
    elif da > db:
        return 1

    cdef long cmp
    cdef long i
    for i in range(da+1):
        fmpz_poly_get_coeff_fmpz(prime_pow.fmpz_ccmp, prime_pow.poly_ccmp, i)
        cmp = fmpz_cmp_si(prime_pow.fmpz_ccmp, 0)
        if cmp < 0:
            return -1
        elif cmp > 0:
            return 1
    assert False

cdef inline int cneg(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    """
    Negation

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the negation
    - ``a`` -- an ``celement`` to be negated
    - ``prec`` -- a long, the precision: ignored
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_neg(out, a)

cdef inline int cadd(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    """
    Addition

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the sum
    - ``a`` -- an ``celement``, the first summand
    - ``b`` -- an ``celement``, the second summand
    - ``prec`` -- a long, the precision: ignored
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_add(out, a, b)

cdef inline bint creduce(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    """
    Reduce modulo a power of the maximal ideal.

    INPUT:

    - ``out`` -- an ``celement`` to store the reduction
    - ``a`` -- the element to be reduced
    - ``prec`` -- a long, the precision to reduce modulo
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT: ``True`` if the reduction is zero; ``False`` otherwise
    """
    if prec == 0:
        csetzero(out, prime_pow)
        return True
    sig_on()
    fmpz_poly_rem(out, a, prime_pow.get_modulus(prec)[0])
    fmpz_poly_scalar_mod_fmpz(out, out, prime_pow.pow_fmpz_t_tmp(prec)[0])
    sig_off()
    return ciszero(out, prime_pow)

cdef inline bint creduce_small(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    """
    Reduce modulo a power of the maximal ideal.

    This function assumes that at most one addition/subtraction has
    happened on reduced inputs.  For integral inputs this translates
    to the assumption that `-p^prec < a < 2p^prec`.

    INPUT:

    - ``out`` -- an ``celement`` to store the reduction
    - ``a`` -- the element to be reduced
    - ``prec`` -- a long, the precision to reduce modulo
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT: ``True`` if the reduction is zero; ``False`` otherwise
    """
    return creduce(out, a, prec, prime_pow)

cdef inline long cremove(celement out, celement a, long prec, PowComputer_ prime_pow, bint reduce_relative=False) except -1:
    """
    Extract the maximum power of the uniformizer dividing this element.

    INPUT:

    - ``out`` -- an ``celement`` to store the unit
    - ``a`` -- the element whose valuation and unit are desired
    - ``prec`` -- a long, used if `a = 0`
    - ``prime_pow`` -- the PowComputer for the ring
    - ``reduce_relative`` -- a bint: whether the final result
      should be reduced at precision ``prec`` (case ``False``)
      or ``prec - valuation`` (case ``True``)

    OUTPUT:

    - if `a = 0`, returns prec (the value of ``out`` is undefined).
      Otherwise, returns the number of times `p` divides `a`.
    """
    if ciszero(a, prime_pow):
        return prec
    cdef long ret = cvaluation(a, prec, prime_pow)
    if ret:
        sig_on()
        fmpz_poly_scalar_divexact_fmpz(out, a, (<PowComputer_flint_unram>prime_pow).pow_fmpz_t_tmp(ret)[0])
        sig_off()
    else:
        fmpz_poly_set(out, a)
    return ret

cdef inline long cvaluation(celement a, long prec, PowComputer_ prime_pow) except -1:
    """
    Return the maximum power of the uniformizer dividing this
    element.

    This function differs from :meth:`cremove` in that the unit is
    discarded.

    INPUT:

    - ``a`` -- the element whose valuation is desired
    - ``prec`` -- a long, used if `a = 0`
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT:

    - if `a = 0`, returns prec.  Otherwise, returns the number of
      times p divides a.
    """
    if ciszero(a, prime_pow):
        return prec
    cdef long ret = maxordp
    cdef long val
    cdef long i
    for i in range(fmpz_poly_length(a)):
        fmpz_poly_get_coeff_fmpz(prime_pow.fmpz_cval, a, i)
        if fmpz_is_zero(prime_pow.fmpz_cval):
            continue
        val = fmpz_remove(prime_pow.fmpz_cval, prime_pow.fmpz_cval, prime_pow.fprime)
        if val < ret:
            ret = val
    return ret

cdef inline bint cisunit(celement a, PowComputer_ prime_pow) except -1:
    """
    Return whether this element has valuation zero.

    INPUT:

    - ``a`` -- the element to test
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT: ``True`` if `a` has valuation 0, and ``False`` otherwise
    """
    fmpz_poly_scalar_mod_fmpz(prime_pow.poly_cisunit, a, prime_pow.fprime)
    return not ciszero(prime_pow.poly_cisunit, prime_pow)

cdef inline int cshift(celement out, celement rem, celement a, long n, long prec, PowComputer_ prime_pow, bint reduce_afterward) except -1:
    """
    Multiplies by a power of the uniformizer.

    INPUT:

    - ``out`` -- a ``celement`` to store the result.  If `n >= 0`
      then out will be set to `a * p^n`
      If `n < 0`, out will be set to `a // p^-n`.
    - ``rem`` -- a ``celement`` to store the remainder of the division
      Should not be aliased with `a`
    - ``a`` -- the element to shift
    - ``n`` -- long, the amount to shift by
    - ``prec`` -- long, a precision modulo which to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    - ``reduce_afterward`` -- whether to reduce afterward
    """
    if n > 0:
        fmpz_poly_zero(rem)
        fmpz_poly_scalar_mul_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(n)[0])
    elif n < 0:
        sig_on()
        fmpz_poly_scalar_mod_fmpz(rem, a, prime_pow.pow_fmpz_t_tmp(-n)[0])
        fmpz_poly_scalar_fdiv_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(-n)[0])
        sig_off()
    else:
        fmpz_poly_zero(rem)
        fmpz_poly_set(out, a)
    if reduce_afterward:
        creduce(out, out, prec, prime_pow)

cdef inline int cshift_notrunc(celement out, celement a, long n, long prec, PowComputer_ prime_pow, bint reduce_afterward) except -1:
    """
    Multiplies by a power of the uniformizer, assuming that the
    valuation of a is at least -n.

    INPUT:

    - ``out`` -- an ``celement`` to store the result.  If `n >= 0`
      then out will be set to `a * p^n`
      If `n < 0`, out will be set to `a // p^-n`.
    - ``a`` -- the element to shift.  Assumes that the valuation of a
      is at least -n.
    - ``n`` -- long, the amount to shift by
    - ``prec`` -- long, a precision modulo which to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    - ``reduce_afterward`` -- whether to reduce afterward
    """
    if n > 0:
        fmpz_poly_scalar_mul_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(n)[0])
    elif n < 0:
        sig_on()
        fmpz_poly_scalar_divexact_fmpz(out, a, prime_pow.pow_fmpz_t_tmp(-n)[0])
        sig_off()
    else:
        fmpz_poly_set(out, a)
    if reduce_afterward:
        creduce(out, out, prec, prime_pow)

cdef inline int csub(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    """
    Subtraction.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the difference
    - ``a`` -- an ``celement``, the first input
    - ``b`` -- an ``celement``, the second input
    - ``prec`` -- a long, the precision: ignored
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_sub(out, a, b)

cdef inline int cinvert(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    """
    Inversion

    The result will be reduced modulo p^prec.

    INPUT:

    - ``out`` -- an ``celement`` to store the inverse
    - ``a`` -- an ``celement``, the element to be inverted
    - ``prec`` -- a long, the precision
    - ``prime_pow`` -- the PowComputer for the ring
    """
    sig_on()
    try:
        fmpz_poly_set(prime_pow.poly_cinv, prime_pow.get_modulus(prec)[0])
        fmpz_poly_primitive_part(prime_pow.poly_cinv, prime_pow.poly_cinv)

        fmpz_poly_content(prime_pow.fmpz_cinv, a)
        fmpz_poly_scalar_divexact_fmpz(out, a, prime_pow.fmpz_cinv)

        fmpz_poly_xgcd(prime_pow.fmpz_cinv2, out, prime_pow.poly_cinv2, out, prime_pow.poly_cinv)
        if fmpz_is_zero(prime_pow.fmpz_cinv2):
            raise ValueError("polynomials are not coprime")

        fmpz_mul(prime_pow.fmpz_cinv2, prime_pow.fmpz_cinv, prime_pow.fmpz_cinv2)
        if not fmpz_invmod(prime_pow.fmpz_cinv2, prime_pow.fmpz_cinv2, prime_pow.pow_fmpz_t_tmp(prec)[0]):
            raise ValueError("content or xgcd is not a unit")
        fmpz_poly_scalar_mul_fmpz(out, out, prime_pow.fmpz_cinv2)

        creduce(out, out, prec, prime_pow)
    finally:
        sig_off()

cdef inline int cmul(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    """
    Multiplication.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the product
    - ``a`` -- an ``celement``, the first input
    - ``b`` -- an ``celement``, the second input
    - ``prec`` -- a long, the precision: ignored
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_mul(out, a, b)

cdef inline int cdivunit(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    """
    Division.

    The inversion is performed modulo p^prec.  Note that no reduction
    is performed after the product.

    INPUT:

    - ``out`` -- an ``celement`` to store the quotient
    - ``a`` -- an ``celement``, the first input
    - ``b`` -- an ``celement``, the second input
    - ``prec`` -- a long, the precision
    - ``prime_pow`` -- the PowComputer for the ring
    """
    cinvert(prime_pow.aliasing, b, prec, prime_pow)
    cmul(out, a, prime_pow.aliasing, prec, prime_pow)

cdef inline int csetone(celement out, PowComputer_ prime_pow) except -1:
    """
    Set to 1.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 1
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_set_ui(out, 1)

cdef inline int csetzero(celement out, PowComputer_ prime_pow) except -1:
    """
    Set to 0.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 0
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_set_ui(out, 0)

cdef inline bint cisone(celement a, PowComputer_ prime_pow) except -1:
    """
    Return whether this element is equal to 1.

    INPUT:

    - ``a`` -- the element to test
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT: ``True`` if `a = 1`, and ``False`` otherwise
    """
    return fmpz_poly_is_one(a)

cdef inline bint ciszero(celement a, PowComputer_ prime_pow) except -1:
    """
    Return whether this element is equal to 0.

    INPUT:

    - ``a`` -- the element to test
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT: ``True`` if `a = 0`, and ``False`` otherwise
    """
    return fmpz_poly_is_zero(a)

cdef inline int cpow(celement out, celement a, mpz_t n, long prec, PowComputer_ prime_pow) except -1:
    """
    Exponentiation.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result
    - ``a`` -- the base
    - ``n`` -- an ``mpz_t``, the exponent
    - ``prec`` -- a long, the working absolute precision
    - ``prime_pow`` -- the PowComputer for the ring
    """
    if mpz_sgn(n) < 0:
        raise NotImplementedError("negative exponent")
    elif mpz_sgn(n) == 0:
        csetone(out, prime_pow)
    elif mpz_even_p(n):
        mpz_divexact_ui(prime_pow.mpz_cpow, n, 2)
        cpow(out, a, prime_pow.mpz_cpow, prec, prime_pow)
        fmpz_poly_sqr(out, out)
    else:
        mpz_sub_ui(prime_pow.mpz_cpow, n, 1)
        cpow(out, a, prime_pow.mpz_cpow, prec, prime_pow)
        fmpz_poly_mul(out, out, a)

    creduce(out, out, prec, prime_pow)

cdef inline int ccopy(celement out, celement a, PowComputer_ prime_pow) except -1:
    """
    Copying.

    INPUT:

    - ``out`` -- the ``celement`` to store the result
    - ``a`` -- the element to copy
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_poly_set(out, a)

cdef inline cpickle(celement a, PowComputer_ prime_pow):
    """
    Serialization into objects that Sage knows how to pickle.

    INPUT:

    - ``a`` -- the element to pickle
    - ``prime_pow`` the PowComputer for the ring

    OUTPUT:

    - a serializable object storing ``a``.
    """
    return fmpz_poly_get_str(a).decode("UTF-8")

cdef inline int cunpickle(celement out, x, PowComputer_ prime_pow) except -1:
    """
    Reconstruction from the output of :meth:`cpickle`.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result
    - ``x`` -- the result of :meth:`cpickle`
    - ``prime_pow`` -- the PowComputer for the ring
    """
    byte_string = x.encode("UTF-8")
    cdef char* c_str = byte_string
    fmpz_poly_set_str(out, c_str)

cdef inline long chash(celement a, long ordp, long prec, PowComputer_ prime_pow) except -1:
    """
    Hashing.

    INPUT:

    - ``a`` -- a ``celement`` storing the underlying element to hash
    - ``ordp`` -- a long storing the valuation
    - ``prec`` -- a long storing the precision
    - ``prime_pow`` -- a PowComputer for the ring
    """
    if ciszero(a, prime_pow):
        return 0

    cdef Integer h = PY_NEW(Integer)
    fmpz_poly_get_coeff_mpz(h.value, a, 0)
    return hash(h)

cdef inline cmodp_rep(fmpz_poly_t rep, fmpz_poly_t value, expansion_mode mode, bint return_list, PowComputer_ prime_pow):
    """
    Compute a polynomial that is reduced modulo p and equivalent to the given value.

    INPUT:

    - ``rep`` -- the reduction mod p
    - ``value`` -- the element to be reduced
    - ``mode`` -- if ``smallest_mode``, the coefficients of the reduction
`     will be between -p/2 and p/2 instead of between 0 and p
    - ``return_list`` -- boolean, whether to return a list of integers giving the coefficients of the expansion
    - ``prime_pow`` -- a PowComputer for the ring
    """
    cdef long i
    cdef fmpz* c
    cdef Integer digit
    sig_on()
    fmpz_poly_scalar_mod_fmpz(rep, value, prime_pow.fprime)
    sig_off()
    if return_list or mode == smallest_mode:
        L = []
        for i in range(fmpz_poly_length(rep)):
            c = fmpz_poly_get_coeff_ptr(rep, i)
            if mode == smallest_mode and fmpz_cmp(c, prime_pow.half_prime) > 0:
                fmpz_sub(c, c, prime_pow.fprime)
            if return_list:
                digit = PY_NEW(Integer)
                fmpz_get_mpz(digit.value, c)
                L.append(digit)
        if return_list:
            return L

# the expansion_mode enum is defined in padic_template_element_header.pxi
cdef inline cexpansion_next(fmpz_poly_t value, expansion_mode mode, long curpower, PowComputer_ prime_pow):
    """
    Return the next digit in a `p`-adic expansion of ``value``.

    INPUT:

    - ``value`` -- the `p`-adic element whose expansion is desired
    - ``mode`` -- either ``simple_mode`` or ``smallest_mode``
    - ``curpower`` -- the current power of `p` for which the coefficient
      is being found.  Only used in ``smallest_mode``
    - ``prime_pow`` -- a ``PowComputer`` holding `p`-adic data
    """
    if mode == teichmuller_mode:
        raise NotImplementedError
    ans = []
    cdef fmpz* c
    cdef long i
    cdef Integer digit
    for i in range(fmpz_poly_length(value)):
        c = fmpz_poly_get_coeff_ptr(value, i)
        fmpz_fdiv_qr(c, prime_pow.fmpz_cexp, c, prime_pow.fprime)
        if mode == smallest_mode and fmpz_cmp(prime_pow.fmpz_cexp, prime_pow.half_prime) > 0:
            fmpz_sub(prime_pow.fmpz_cexp, prime_pow.fmpz_cexp, prime_pow.fprime)
            fmpz_add_ui(c, c, 1)
        digit = PY_NEW(Integer)
        fmpz_get_mpz(digit.value, prime_pow.fmpz_cexp)
        ans.append(digit)
    _fmpz_poly_normalise(value)
    return trim_zeros(ans) # defined in sage.rings.padics.misc and imported in padic_template_element

cdef inline cexpansion_getitem(fmpz_poly_t value, long m, PowComputer_ prime_pow):
    """
    Return the `m`-th `p`-adic digit in the ``simple_mode`` expansion.

    INPUT:

    - ``value`` -- the `p`-adic element whose expansion is desired
    - ``m`` -- nonnegative integer: which entry in the `p`-adic expansion to return
    - ``prime_pow`` -- a ``PowComputer`` holding `p`-adic data
    """
    ans = []
    cdef fmpz* c
    cdef long i
    cdef Integer digit
    for i in range(fmpz_poly_length(value)):
        c = fmpz_poly_get_coeff_ptr(value, i)
        if m > 0:
            fmpz_fdiv_q(prime_pow.fmpz_cexp, c, prime_pow.pow_fmpz_t_tmp(m)[0])
            fmpz_mod(prime_pow.fmpz_cexp, prime_pow.fmpz_cexp, prime_pow.fprime)
        else:
            fmpz_mod(prime_pow.fmpz_cexp, c, prime_pow.fprime)
        digit = PY_NEW(Integer)
        fmpz_get_mpz(digit.value, prime_pow.fmpz_cexp)
        ans.append(digit)
    _fmpz_poly_normalise(value)
    return trim_zeros(ans) # defined in sage.rings.padics.misc and imported in padic_template_element

# The element is filled in for zero in the p-adic expansion if necessary.
_expansion_zero = []

cdef list ccoefficients(celement x, long valshift, long prec, PowComputer_ prime_pow):
    """
    Return a list of coefficients, as elements that can be converted into the base ring.

    INPUT:

    - ``x`` -- a ``celement`` giving the underlying `p`-adic element, or possibly its unit part
    - ``valshift`` -- a long giving the power of the uniformizer to shift `x` by
    - ``prec`` -- a long, the (relative) precision desired, used in rational reconstruction
    - ``prime_pow`` -- the ``PowComputer`` of the ring
    """
    cdef Integer ansz
    cdef Rational ansq
    cdef long i
    ans = []
    for i in range(fmpz_poly_length(x)):
        if valshift >= 0:
            ansz = PY_NEW(Integer)
            fmpz_poly_get_coeff_mpz(ansz.value, x, i)
            if valshift > 0:
                mpz_mul(ansz.value, ansz.value, prime_pow.pow_mpz_t_tmp(valshift))
            ans.append(ansz)
        else:
            ansq = Rational.__new__(Rational)
            fmpz_poly_get_coeff_mpz(mpq_numref(ansq.value), x, i)
            mpz_set(mpq_denref(ansq.value), prime_pow.pow_mpz_t_tmp(-valshift))
            mpq_canonicalize(ansq.value)
            ans.append(ansq)
    return ans

cdef int cteichmuller(celement out, celement value, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Teichmuller lifting.

    INPUT:

    - ``out`` -- an ``celement`` which is set to a `q-1` root of unity
                 congruent to `value` mod `\pi`; or 0 if `a \equiv 0
                 \pmod{\pi}`.
    - ``value`` -- an ``celement``, the element mod `\pi` to lift
    - ``prec`` -- a long, the precision to which to lift
    - ``prime_pow`` -- the ``PowComputer`` of the ring

    ALGORITHM:

    We use Hensel lifting to solve the equation `f(T)=T^q-T`. Instead of
    dividing by the derivative of `f`, we divide by `( q - 1 )` whose first
    digits coincide with `f'`. This does probably not yield quadratic
    convergence but taking inverses would be much more expensive than what is
    done here.
    """
    fmpz_poly_set(out, value)

    if prec == 0:
        return 0

    # fmpz_ctm = 1 / (1 - q) (mod p^prec)
    fmpz_set_ui(prime_pow.fmpz_ctm, 1)
    fmpz_sub(prime_pow.fmpz_ctm, prime_pow.fmpz_ctm, prime_pow.q)
    fmpz_invmod(prime_pow.fmpz_ctm, prime_pow.fmpz_ctm, prime_pow.pow_fmpz_t_tmp(prec)[0])
    while True:
        # poly_ctm = out + fmpz_ctm*(out^q - out)
        fmpz_get_mpz(prime_pow.mpz_ctm, prime_pow.q)
        cpow(prime_pow.poly_ctm, out, prime_pow.mpz_ctm, prec, prime_pow)
        csub(prime_pow.poly_ctm, prime_pow.poly_ctm, out, prec, prime_pow)
        fmpz_poly_scalar_mul_fmpz(prime_pow.poly_ctm, prime_pow.poly_ctm, prime_pow.fmpz_ctm)
        cadd(prime_pow.poly_ctm, prime_pow.poly_ctm, out, prec, prime_pow)
        creduce(prime_pow.poly_ctm, prime_pow.poly_ctm, prec, prime_pow)
        # break if out == poly_ctm
        if ccmp(prime_pow.poly_ctm, out, prec, False, False, prime_pow) == 0:
            return 0
        # out = poly_ctm
        fmpz_poly_set(out, prime_pow.poly_ctm)

cdef int cconv(celement out, x, long prec, long valshift, PowComputer_ prime_pow) except -2:
    """
    Conversion from other Sage types.

    INPUT:

    - ``out`` -- an ``celement`` to store the output

    - ``x`` -- a Sage element that can be converted to a `p`-adic element

    - ``prec`` -- a long, giving the precision desired; absolute if
      `valshift = 0`, relative if `valshift != 0`

    - ``valshift`` -- the power of the uniformizer to divide by before
      storing the result in ``out``

    - ``prime_pow`` -- a PowComputer for the ring
    """
    cdef long i
    cdef long degree

    if PyList_Check(x):
        for i in range(len(x)):
            cconv(prime_pow.poly_cconv, x[i], prec, valshift, prime_pow)
            degree = fmpz_poly_degree(prime_pow.poly_cconv)
            if degree == -1:
                continue
            elif degree == 0:
                fmpz_poly_get_coeff_fmpz(prime_pow.fmpz_cconv, prime_pow.poly_cconv, 0)
                fmpz_poly_set_coeff_fmpz(out, i, prime_pow.fmpz_cconv)
            else:
                raise ValueError
        creduce(out, out, prec, prime_pow)
    else:
        cconv_shared(prime_pow.mpz_cconv, x, prec, valshift, prime_pow)
        fmpz_poly_set_mpz(out, prime_pow.mpz_cconv)

cdef inline long cconv_mpq_t(celement out, mpq_t x, long prec, bint absolute, PowComputer_ prime_pow) except? -10000:
    """
    A fast pathway for conversion of rationals that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``celement`` to store the output
    - ``x`` -- an ``mpq_t`` giving the integer to be converted
    - ``prec`` -- a long, giving the precision desired: absolute or
      relative depending on the ``absolute`` input
    - ``absolute`` -- if ``False`` then extracts the valuation and returns
      it, storing the unit in ``out``; if ``True`` then just reduces ``x``
      modulo the precision.
    - ``prime_pow`` -- a PowComputer for the ring

    OUTPUT:

    If ``absolute`` is ``False`` then returns the valuation that was
    extracted (``maxordp`` when `x = 0`).
    """
    cdef long val
    val = cconv_mpq_t_shared(prime_pow.mpz_cconv, x, prec, absolute, prime_pow)
    fmpz_poly_set_mpz(out, prime_pow.mpz_cconv)
    return val

cdef inline int cconv_mpq_t_out(mpq_t out, celement x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    """
    Convert the underlying `p`-adic element into a rational.

    - ``out`` -- gives a rational approximating the input.  Currently uses rational reconstruction but
                 may change in the future to use a more naive method
    - ``x`` -- an ``celement`` giving the underlying `p`-adic element
    - ``valshift`` -- a long giving the power of `p` to shift `x` by
    -` ``prec`` -- a long, the precision of ``x``, used in rational reconstruction
    - ``prime_pow`` -- a PowComputer for the ring
    """
    cdef long degree = fmpz_poly_degree(x)
    if degree > 0:
        raise ValueError
    elif degree == -1:
        mpz_set_ui(prime_pow.mpz_cconv, 0)
    else:
        fmpz_poly_get_coeff_mpz(prime_pow.mpz_cconv, x, 0)

    cconv_mpq_t_out_shared(out, prime_pow.mpz_cconv, valshift, prec, prime_pow)

cdef inline long cconv_mpz_t(celement out, mpz_t x, long prec, bint absolute, PowComputer_ prime_pow) except -2:
    """
    A fast pathway for conversion of integers that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``celement`` to store the output
    - ``x`` -- an ``mpz_t`` giving the integer to be converted
    - ``prec`` -- a long, giving the precision desired: absolute or
      relative depending on the ``absolute`` input
    - ``absolute`` -- if False then extracts the valuation and returns
      it, storing the unit in ``out``; if ``True`` then just reduces ``x``
      modulo the precision
    - ``prime_pow`` -- a PowComputer for the ring

    OUTPUT:

    - If ``absolute`` is False then returns the valuation that was
      extracted (``maxordp`` when `x = 0`).
    """
    cdef long val
    val = cconv_mpz_t_shared(prime_pow.mpz_cconv, x, prec, absolute, prime_pow)
    fmpz_poly_set_mpz(out, prime_pow.mpz_cconv)
    return val

cdef inline int cconv_mpz_t_out(mpz_t out, celement x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    """
    Convert the underlying `p`-adic element into an integer if
    possible.

    - ``out`` -- stores the resulting integer as an integer between 0
      and `p^{prec + valshift}`
    - ``x`` -- an ``celement`` giving the underlying `p`-adic element
    - ``valshift`` -- a long giving the power of `p` to shift `x` by
    -` ``prec`` -- a long, the precision of ``x``: currently not used
    - ``prime_pow`` -- a PowComputer for the ring
    """
    cdef long degree = fmpz_poly_degree(x)
    if degree > 0:
        raise ValueError
    elif degree == -1:
        mpz_set_ui(prime_pow.mpz_cconv, 0)
    else:
        fmpz_poly_get_coeff_mpz(prime_pow.mpz_cconv, x, 0)

    cconv_mpz_t_out_shared(out, prime_pow.mpz_cconv, valshift, prec, prime_pow)


## Extra functions ##

cdef cmatrix_mod_pn(celement a, long aprec, long valshift, PowComputer_ prime_pow):
    r"""
    Return the matrix of right multiplication by the element on
    the power basis `1, x, x^2, \ldots, x^{d-1}` for this
    extension field.  Thus the *rows* of this matrix give the
    images of each of the `x^i`.  The entries of the matrices are
    IntegerMod elements, defined modulo ``p^(self.absprec() / e)``.
    """
    cdef Py_ssize_t i, j, d, deg = prime_pow.deg
    cdef int fail
    R = Zmod(prime_pow.pow_Integer(aprec))
    cdef IntegerMod_abstract zero = R(0)
    cdef IntegerMod_abstract item
    L = []
    cshift_notrunc(prime_pow.poly_matmod, a, valshift, aprec, prime_pow, True)
    for i in range(deg):
        L.append([])
        d = fmpz_poly_degree(prime_pow.poly_matmod)
        for j in range(d+1):
            item = zero._new_c_from_long(0)
            fmpz_poly_get_coeff_mpz(prime_pow.mpz_matmod, prime_pow.poly_matmod, j)
            item.set_from_mpz(prime_pow.mpz_matmod)
            L[-1].append(item)
        for j in range(d+1,deg):
            L[-1].append(zero)
        fmpz_poly_shift_left(prime_pow.poly_matmod, prime_pow.poly_matmod, 1)
        creduce(prime_pow.poly_matmod, prime_pow.poly_matmod, aprec, prime_pow)
    from sage.matrix.constructor import matrix
    return matrix(R, deg, deg, L)
