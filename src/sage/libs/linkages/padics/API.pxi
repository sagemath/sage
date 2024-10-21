# sage_setup: distribution = sagemath-pari
"""
This file defines the common API for p-adic elements.
Elements using different precision models (e.g. capped relative,
capped absolute and fixed modulus) can plug in different linkage
files.  By separating the linkage file from the precision templates we
reduce code duplication.

Each linkage file should implement the following functions.  These
functions will be implemented for many backends, usually in C (thus
the c at the beginning of the function names).  The remainder of the
file gives the function signatures.

- :meth:`cconstruct` -- construct a new element
- :meth:`cdestruct` -- deallocate an element
- :meth:`ccmp` -- comparison of two elements
- :meth:`cneg` -- negation
- :meth:`cadd` -- addition
- :meth:`creduce` -- reduce modulo a power of the maximal ideal
- :meth:`creduce_small` -- reduce modulo a power of the maximal ideal,
  assuming only a single addition/subtraction
- :meth:`cremove` -- extract the maximum power of the uniformizer
  dividing this element
- :meth:`cvaluation` -- return the maximum power of the uniformizer
  dividing this element
- :meth:`cisunit` -- returns whether this element has valuation zero
- :meth:`cshift` -- multiplies by a power of the uniformizer
- :meth:`cshift_notrunc` -- multiplies by a power of the uniformizer,
  assuming no truncation
- :meth:`csub` -- subtraction
- :meth:`cinvert` -- inversion
- :meth:`cmul` -- multiplication
- :meth:`cdivunit` -- division
- :meth:`csetone` -- sets to 1
- :meth:`csetzero` -- sets to 0
- :meth:`cisone` -- tests for equality with 1
- :meth:`ciszero` -- tests for equality with 0
- :meth:`cpow` -- exponentiation
- :meth:`ccopy` -- copying
- :meth:`cpickle` -- serialization into objects that Sage knows how to
  pickle
- :meth:`cunpickle` -- reconstruction from the output of :meth:`cpickle`
- :meth:`chash` -- hashing
- :meth:`cexpansion_next` -- gets the next digit in the series expansion
- :meth:`cexpansion_getitem` -- gets a specified digit in the series expansion
- :meth:`ccoefficients` -- list of coefficients as elements of the base ring
- :meth:`cteichmuller` -- Teichmuller lifting
- :meth:`cconv` -- conversion from other types in Sage
- :meth:`cconv_mpz_t` -- conversion from mpz_t, separated for speed and
  convenience
- :meth:`cconv_mpq_t` -- conversion from mpq_t, separated for speed and
  convenience
- :meth:`cconv_mpz_t_out` -- conversion into an mpz_t
- :meth:`cconv_mpq_t_out` -- conversion into an mpq_t
- _expansion_zero -- the entry that should be used for zero in p-adic expansions

The gluing file should ctypedef celement as appropriate.

.. NOTE::

    This particular file is not included anywhere.  It just defines
    the interface implemented by other pxi files in this directory.

AUTHORS:

- David Roe (2012-3-1) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef inline int cconstruct(celement value, PowComputer_class prime_pow) except -1:
    """
    Construct a new element.

    INPUT:

    - ``unit`` -- an ``celement`` to be initialized
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline int cdestruct(celement value, PowComputer_class prime_pow) except -1:
    """
    Deallocate an element.

    INPUT:

    - ``unit`` -- an ``celement`` to be cleared
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline int ccmp(celement a, celement b, long prec, bint reduce_a, bint reduce_b, PowComputer_class prime_pow) except -2:
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
    pass

cdef inline int cneg(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Negation

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- an ``celement`` to store the negation
    - ``a`` -- an ``celement`` to be negated
    - ``prec`` -- a long, the precision: ignored
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline int cadd(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
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
    pass

cdef inline bint creduce(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Reduce modulo a power of the maximal ideal.

    INPUT:

    - ``out`` -- an ``celement`` to store the reduction
    - ``a`` -- the element to be reduced
    - ``prec`` -- a long, the precision to reduce modulo
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT: returns ``True`` if the reduction is zero; ``False`` otherwise
    """
    pass

cdef inline bint creduce_small(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
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

    OUTPUT:

    - returns ``True`` if the reduction is zero; ``False`` otherwise.
    """
    pass

cdef inline long cremove(celement out, celement a, long prec, PowComputer_class prime_pow, bint reduce_relative=False) except -1:
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
    pass

cdef inline long cvaluation(celement a, long prec, PowComputer_class prime_pow) except -1:
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
    pass

cdef inline bint cisunit(celement a, PowComputer_class prime_pow) except -1:
    """
    Return whether this element has valuation zero.

    INPUT:

    - ``a`` -- the element to test
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT: returns ``True`` if `a` has valuation 0, and ``False`` otherwise
    """
    pass

cdef inline int cshift(celement out, celement rem, celement a, long n, long prec, PowComputer_class prime_pow, bint reduce_afterward) except -1:
    """
    Multiply by a power of the uniformizer.

    INPUT:

    - ``out`` -- an ``celement`` to store the result.  If `n >= 0`
      then out will be set to `a * p^n`
      If `n < 0`, out will be set to `a // p^-n`.
    - ``rem`` -- a ``celement`` to store the remainder, when `n < 0`
    - ``a`` -- the element to shift
    - ``n`` -- long, the amount to shift by
    - ``prec`` -- long, a precision modulo which to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    - ``reduce_afterward`` -- whether to reduce afterward
    """
    pass

cdef inline int cshift_notrunc(celement out, celement a, long n, long prec, PowComputer_class prime_pow, bint reduce_afterward) except -1:
    """
    Multiply by a power of the uniformizer, assuming that the
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
    pass

cdef inline int csub(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
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
    pass

cdef inline int cinvert(celement out, celement a, long prec, PowComputer_class prime_pow) except -1:
    """
    Inversion

    The result will be reduced modulo p^prec.

    INPUT:

    - ``out`` -- an ``celement`` to store the inverse
    - ``a`` -- an ``celement``, the element to be inverted
    - ``prec`` -- a long, the precision
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline int cmul(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
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
    pass

cdef inline int cdivunit(celement out, celement a, celement b, long prec, PowComputer_class prime_pow) except -1:
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
    pass

cdef inline int csetone(celement out, PowComputer_class prime_pow) except -1:
    """
    Set to 1.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 1
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline int csetzero(celement out, PowComputer_class prime_pow) except -1:
    """
    Set to 0.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 0
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline bint cisone(celement a, PowComputer_class prime_pow) except -1:
    """
    Return whether this element is equal to 1.

    INPUT:

    - ``a`` -- the element to test
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT:

    - returns ``True`` if `a = 1`, and ``False`` otherwise.
    """
    pass

cdef inline bint ciszero(celement a, PowComputer_class prime_pow) except -1:
    """
    Return whether this element is equal to 0.

    INPUT:

    - ``a`` -- the element to test
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT:

    - returns ``True`` if `a = 0`, and ``False`` otherwise.
    """
    pass

cdef inline int cpow(celement out, celement a, mpz_t n, long prec, PowComputer_class prime_pow) except -1:
    """
    Exponentiation.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result
    - ``a`` -- the base
    - ``n`` -- an ``mpz_t``, the exponent
    - ``prec`` -- a long, the working absolute precision
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline int ccopy(celement out, celement a, PowComputer_class prime_pow) except -1:
    """
    Copying.

    INPUT:

    - ``out`` -- the ``celement`` to store the result
    - ``a`` -- the element to copy
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline cpickle(celement a, PowComputer_class prime_pow):
    """
    Serialization into objects that Sage knows how to pickle.

    INPUT:

    - ``a`` -- the element to pickle
    - ``prime_pow`` the PowComputer for the ring

    OUTPUT:

    - a serializable object storing ``a``.
    """
    pass

cdef inline int cunpickle(celement out, x, PowComputer_class prime_pow) except -1:
    """
    Reconstruction from the output of :meth:`cpickle`.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result
    - ``x`` -- the result of :meth:`cpickle`
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline long chash(celement a, long ordp, long prec, PowComputer_class prime_pow) except -1:
    """
    Hashing.

    INPUT:

    - ``a`` -- an ``celement`` storing the underlying element to hash
    - ``ordp`` -- a long storing the valuation
    - ``prec`` -- a long storing the precision
    - ``prime_pow`` -- a PowComputer for the ring
    """
    pass

# the expansion_mode enum is defined in padic_template_element_header.pxi
cdef inline cexpansion_next(celement value, expansion_mode mode, long curpower, PowComputer_ prime_pow):
    r"""
    Return the next digit in a `\pi`-adic expansion of ``value``.

    INPUT:

    - ``value`` -- the `\pi`-adic element whose expansion is desired
    - ``mode`` -- either ``simple_mode`` or ``smallest_mode``
    - ``curpower`` -- the current power of `\pi` for which the coefficient
      is being found.  Only used in ``smallest_mode``
    - ``prime_pow`` -- a ``PowComputer`` holding `\pi`-adic data
    """
    pass

cdef inline cexpansion_getitem(celement value, long m, PowComputer_ prime_pow):
    r"""
    Return the `m`-th `\pi`-adic digit in the ``simple_mode`` expansion.

    INPUT:

    - ``value`` -- the `\pi`-adic element whose expansion is desired
    - ``m`` -- nonnegative integer: which entry in the `\pi`-adic expansion to return
    - ``prime_pow`` -- a ``PowComputer`` holding `\pi`-adic data
    """
    pass

# The element is filled in for zero in the p-adic expansion if necessary.
# It could be [] for some other linkages.
_expansion_zero = Integer(0)

cdef list ccoefficients(celement x, long valshift, PowComputer_class prime_pow):
    """
    Return a list of coefficients, as elements that can be converted into the base ring.

    INPUT:

    - ``x`` -- a ``celement`` giving the underlying `p`-adic element, or possibly its unit part
    - ``valshift`` -- a long giving the power of the uniformizer to shift `x` by
    - ``prec`` -- a long, the (relative) precision desired, used in rational reconstruction
    - ``prime_pow`` -- the ``PowComputer`` of the ring
    """
    pass

cdef int cteichmuller(celement out, celement value, long prec, PowComputer_class prime_pow) except -1:
    r"""
    Teichmuller lifting.

    INPUT:

    - ``out`` -- an ``celement`` which is set to a `q-1` root of unity
      congruent to `value` mod `\pi`; or 0 if `a \equiv 0 \pmod{\pi}`
    - ``value`` -- an ``celement``, the element mod `\pi` to lift
    - ``prec`` -- a long, the precision to which to lift
    - ``prime_pow`` -- the ``PowComputer`` of the ring
    """
    pass

cdef int cconv(celement out, x, long prec, long valshift, PowComputer_class prime_pow) except -2:
    """
    Conversion from other Sage types.

    INPUT:

    - ``out`` -- an ``celement`` to store the output
    - ``x`` -- a Sage element that can be converted to a `p`-adic element
    - ``prec`` -- a long, giving the precision desired: absolute if
      `valshift = 0`, relative if `valshift > 0`
    - ``valshift`` -- the power of the uniformizer to divide by before
      storing the result in ``out``
    - ``prime_pow`` -- a PowComputer for the ring
    """
    pass

cdef inline long cconv_mpz_t(celement out, mpz_t x, long prec, bint absolute, PowComputer_class prime_pow) except -2:
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
    pass

cdef inline int cconv_mpz_t_out(mpz_t out, celement x, long valshift, long prec, PowComputer_class prime_pow) except -1:
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
    pass

cdef inline long cconv_mpq_t(celement out, mpq_t x, long prec, bint absolute, PowComputer_class prime_pow) except? -10000:
    """
    A fast pathway for conversion of rationals that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- an ``celement`` to store the output
    - ``x`` -- an ``mpq_t`` giving the rational to be converted
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
    pass

cdef inline int cconv_mpq_t_out(mpq_t out, celement x, long valshift, long prec, PowComputer_class prime_pow) except -1:
    """
    Convert the underlying `p`-adic element into a rational.

    - ``out`` -- gives a rational approximating the input.  Currently
      uses rational reconstruction but may change in the future to use a more
      naive method
    - ``x`` -- an ``celement`` giving the underlying `p`-adic element
    - ``valshift`` -- a long giving the power of the uniformizer to shift `x` by
    - ``prec`` -- a long, the precision of ``x``, used in rational
      reconstruction
    - ``prime_pow`` -- a PowComputer for the ring
    """
    pass

# In order to ensure that this file is not actually included, we make this file cause a compile error.
cdef UndefinedType raise_compile_error
