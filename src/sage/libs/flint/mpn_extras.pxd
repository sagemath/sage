# distutils: libraries = flint
# distutils: depends = flint/mpn_extras.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void flint_mpn_debug(mp_srcptr x, mp_size_t xsize)
    # Prints debug information about ``(x, xsize)`` to ``stdout``.
    # In particular, this will print binary representations of all the limbs.

    int flint_mpn_zero_p(mp_srcptr x, mp_size_t xsize)
    # Returns `1` if all limbs of ``(x, xsize)`` are zero, otherwise `0`.

    mp_limb_t flint_mpn_mul(mp_ptr z, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
    # Sets ``(z, xn+yn)`` to the product of ``(x, xn)`` and ``(y, yn)``
    # and returns the top limb of the result.
    # We require `xn \ge yn \ge 1`
    # and that ``z`` is not aliased with either input operand.
    # This function uses FFT multiplication if the operands are large enough
    # and otherwise calls ``mpn_mul``.

    void flint_mpn_mul_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n)
    # Sets ``z`` to the product of ``(x, n)`` and ``(y, n)``.
    # We require `n \ge 1`
    # and that ``z`` is not aliased with either input operand.
    # This function uses FFT multiplication if the operands are large enough
    # and otherwise calls ``mpn_mul_n``.

    void flint_mpn_sqr(mp_ptr z, mp_srcptr x, mp_size_t n)
    # Sets ``z`` to the square of ``(x, n)``.
    # We require `n \ge 1`
    # and that ``z`` is not aliased with either input operand.
    # This function uses FFT multiplication if the operands are large enough
    # and otherwise calls ``mpn_sqr``.

    mp_size_t flint_mpn_fmms1(mp_ptr y, mp_limb_t a1, mp_srcptr x1, mp_limb_t a2, mp_srcptr x2, mp_size_t n)
    # Given not-necessarily-normalized `x_1` and `x_2` of length `n > 0` and output `y` of length `n`, try to compute `y = a_1\cdot x_1 - a_2\cdot x_2`.
    # Return the normalized length of `y` if `y \ge 0` and `y` fits into `n` limbs. Otherwise, return `-1`.
    # `y` may alias `x_1` but is not allowed to alias `x_2`.

    int flint_mpn_divisible_1_odd(mp_srcptr x, mp_size_t xsize, mp_limb_t d)
    # Expression determining whether ``(x, xsize)`` is divisible by the
    # ``mp_limb_t d`` which is assumed to be odd-valued and at least `3`.
    # This function is implemented as a macro.

    mp_size_t flint_mpn_divexact_1(mp_ptr x, mp_size_t xsize, mp_limb_t d)
    # Divides `x` once by a known single-limb divisor, returns the new size.

    mp_size_t flint_mpn_remove_2exp(mp_ptr x, mp_size_t xsize, flint_bitcnt_t *bits)
    # Divides ``(x, xsize)`` by `2^n` where `n` is the number of trailing
    # zero bits in `x`. The new size of `x` is returned, and `n` is stored in
    # the bits argument. `x` may not be zero.

    mp_size_t flint_mpn_remove_power_ascending(mp_ptr x, mp_size_t xsize, mp_ptr p, mp_size_t psize, ulong *exp)
    # Divides ``(x, xsize)`` by the largest power `n` of ``(p, psize)``
    # that is an exact divisor of `x`. The new size of `x` is returned, and
    # `n` is stored in the ``exp`` argument. `x` may not be zero, and `p`
    # must be greater than `2`.
    # This function works by testing divisibility by ascending squares
    # `p, p^2, p^4, p^8, \dotsc`, making it efficient for removing potentially
    # large powers. Because of its high overhead, it should not be used as
    # the first stage of trial division.

    int flint_mpn_factor_trial(mp_srcptr x, mp_size_t xsize, slong start, slong stop)
    # Searches for a factor of ``(x, xsize)`` among the primes in positions
    # ``start, ..., stop-1`` of ``flint_primes``. Returns `i` if
    # ``flint_primes[i]`` is a factor, otherwise returns `0` if no factor
    # is found. It is assumed that ``start >= 1``.

    int flint_mpn_factor_trial_tree(slong * factors, mp_srcptr x, mp_size_t xsize, slong num_primes)
    # Searches for a factor of ``(x, xsize)`` among the primes in positions
    # approximately in the range ``0, ..., num_primes - 1`` of ``flint_primes``.
    # Returns the number of prime factors found and fills ``factors`` with their
    # indices in ``flint_primes``. It is assumed that ``num_primes`` is in the
    # range ``0, ..., 3512``.
    # If the input fits in a small ``fmpz`` the number is fully factored instead.
    # The algorithm used is a tree based gcd with a product of primes, the tree
    # for which is cached globally (it is threadsafe).

    int flint_mpn_divides(mp_ptr q, mp_srcptr array1, mp_size_t limbs1, mp_srcptr arrayg, mp_size_t limbsg, mp_ptr temp)
    # If ``(arrayg, limbsg)`` divides ``(array1, limbs1)`` then
    # ``(q, limbs1 - limbsg + 1)`` is set to the quotient and 1 is
    # returned, otherwise 0 is returned. The temporary space ``temp``
    # must have space for ``limbsg`` limbs.
    # Assumes ``limbs1 >= limbsg > 0``.

    mp_limb_t flint_mpn_preinv1(mp_limb_t d, mp_limb_t d2)
    # Computes a precomputed inverse from the leading two limbs of the
    # divisor ``b, n`` to be used with the ``preinv1`` functions.
    # We require the most significant bit of ``b, n`` to be 1.

    mp_limb_t flint_mpn_divrem_preinv1(mp_ptr q, mp_ptr a, mp_size_t m, mp_srcptr b, mp_size_t n, mp_limb_t dinv)
    # Divide ``a, m`` by ``b, n``, returning the high limb of the
    # quotient (which will either be 0 or 1), storing the remainder in-place
    # in ``a, n`` and the rest of the quotient in ``q, m - n``.
    # We require the most significant bit of ``b, n`` to be 1.
    # ``dinv`` must be computed from ``b[n - 1]``, ``b[n - 2]`` by
    # ``flint_mpn_preinv1``. We also require ``m >= n >= 2``.

    void flint_mpn_mulmod_preinv1(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_limb_t dinv, ulong norm)
    # Given a normalised integer `d` with precomputed inverse ``dinv``
    # provided by ``flint_mpn_preinv1``, computes `ab \pmod{d}` and
    # stores the result in `r`. Each of `a`, `b` and `r` is expected to
    # have `n` limbs of space, with zero padding if necessary.
    # The value ``norm`` is provided for convenience. If `a`, `b` and
    # `d` have been shifted left by ``norm`` bits so that `d` is
    # normalised, then `r` will be shifted right by ``norm`` bits
    # so that it has the same shift as all the inputs.
    # We require `a` and `b` to be reduced modulo `n` before calling the
    # function.

    void flint_mpn_preinvn(mp_ptr dinv, mp_srcptr d, mp_size_t n)
    # Compute an `n` limb precomputed inverse ``dinv`` of the `n` limb
    # integer `d`.
    # We require that `d` is normalised, i.e. with the most significant
    # bit of the most significant limb set.

    void flint_mpn_mod_preinvn(mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv)
    # Given a normalised integer `d` of `n` limbs, with precomputed inverse
    # ``dinv`` provided by ``flint_mpn_preinvn`` and integer `a` of `m`
    # limbs, computes `a \pmod{d}` and stores the result in-place in the lower
    # `n` limbs of `a`. The remaining limbs of `a` are destroyed.
    # We require `m \geq n`. No aliasing of `a` with any of the other operands
    # is permitted.
    # Note that this function is not always as fast as ordinary division.

    mp_limb_t flint_mpn_divrem_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv)
    # Given a normalised integer `d` with precomputed inverse ``dinv``
    # provided by ``flint_mpn_preinvn``, computes the quotient of `a` by `d`
    # and stores the result in `q` and the remainder in the lower `n` limbs of
    # `a`. The remaining limbs of `a` are destroyed.
    # The value `q` is expected to have space for `m - n` limbs and we require
    # `m \ge n`. No aliasing is permitted between `q` and `a` or between these
    # and any of the other operands.
    # Note that this function is not always as fast as ordinary division.

    void flint_mpn_mulmod_preinvn(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)
    # Given a normalised integer `d` with precomputed inverse ``dinv``
    # provided by ``flint_mpn_preinvn``, computes `ab \pmod{d}` and
    # stores the result in `r`. Each of `a`, `b` and `r` is expected to
    # have `n` limbs of space, with zero padding if necessary.
    # The value ``norm`` is provided for convenience. If `a`, `b` and
    # `d` have been shifted left by ``norm`` bits so that `d` is
    # normalised, then `r` will be shifted right by ``norm`` bits
    # so that it has the same shift as all the inputs.
    # We require `a` and `b` to be reduced modulo `n` before calling the
    # function.
    # Note that this function is not always as fast as ordinary division.

    mp_size_t flint_mpn_gcd_full2(mp_ptr arrayg, mp_srcptr array1, mp_size_t limbs1, mp_srcptr array2, mp_size_t limbs2, mp_ptr temp)
    # Sets ``(arrayg, retvalue)`` to the gcd of ``(array1, limbs1)`` and
    # ``(array2, limbs2)``.
    # The only assumption is that neither ``limbs1`` nor ``limbs2`` is
    # zero.
    # The function must be supplied with ``limbs1 + limbs2`` limbs of temporary
    # space, or ``NULL`` must be passed to ``temp`` if the function should
    # allocate its own space.

    mp_size_t flint_mpn_gcd_full(mp_ptr arrayg, mp_srcptr array1, mp_size_t limbs1, mp_srcptr array2, mp_size_t limbs2)
    # Sets ``(arrayg, retvalue)`` to the gcd of ``(array1, limbs1)`` and
    # ``(array2, limbs2)``.
    # The only assumption is that neither ``limbs1`` nor ``limbs2`` is
    # zero.

    void flint_mpn_rrandom(mp_limb_t *rp, gmp_randstate_t state, mp_size_t n)
    # Generates a random number with ``n`` limbs and stores
    # it on ``rp``. The number it generates will tend to have
    # long strings of zeros and ones in the binary representation.
    # Useful for testing functions and algorithms, since this kind of random
    # numbers have proven to be more likely to trigger corner-case bugs.

    void flint_mpn_urandomb(mp_limb_t *rp, gmp_randstate_t state, flint_bitcnt_t n)
    # Generates a uniform random number of ``n`` bits and stores
    # it on ``rp``.
