r"""
This linkage file implements the padics API using Sage Polynomials.

It is included into Polynomial_ram.pxi and Polynomial_unram.pxi,
where functions that depend on the ramification of the defining polynomial are placed.

.. NOTE::

    There are no doctests in this file since the functions here cannot be
    called directly from Python. Testing of this function is necessarily
    indirect and mostly done through arithmetic black-box tests that are part
    of the test suites of the `p`-adic parents.

AUTHORS:

- David Roe, Julian Rüth (2017-06-11): initial version
"""
#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#                     2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from cpython.list cimport *
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.ext.stdsage cimport PY_NEW
from copy import copy
from sage.rings.padics.common_conversion cimport cconv_mpz_t_out_shared, cconv_mpz_t_shared, cconv_mpq_t_out_shared, cconv_mpq_t_shared, cconv_shared

DEF CELEMENT_IS_PY_OBJECT = True

cdef inline int cconstruct(celement value, PowComputer_ prime_pow) except -1:
    r"""
    Construct a new element.

    INPUT:

    - ``value`` -- a ``celement`` to be initialized

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    .. NOTE::

        For Polynomial elements, this function calls the ``__init__`` method of
        ``celement``. When creating ``celement`` instances, we use ``__new__`` which
        does not call ``__init__``, because of the metaclass that we are using.
    """
    value.__init__(prime_pow.poly_ring)

cdef inline int cdestruct(celement value, PowComputer_ prime_pow) except -1:
    r"""
    Deallocate ``value``.

    INPUT:

    - ``value`` -- a ``celement`` to be cleared

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    .. NOTE::

        This function has no effect for polynomials. There is no manual garbage
        collection necessary, as ``celement`` is a managed Python type.
    """
    pass

cdef inline int ccmp(celement a, celement b, long prec, bint reduce_a, bint reduce_b, PowComputer_ prime_pow) except -2:
    r"""
    Return the comparison of ``a`` and ``b``.

    This function returns -1, 0, or 1 for use with ``cmp`` of Python 2. In
    particular, this function return 0 if ``a`` and ``b`` differ by an element
    of valuation at least ``prec``.

    .. NOTE::

        You should not rely on whether this function returns -1 or 1 as this
        might not be consistent.

    INPUT:

    - ``a`` -- a ``celement``

    - ``b`` -- a ``celement``

    - ``prec`` -- a ``long``, the precision to which the comparison should be
      performed

    - ``reduce_a`` -- whether ``a`` is already reduced with respect to ``prec``

    - ``reduce_b`` -- whether ``b`` is already reduced with respect to ``prec``

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    if not (reduce_a or reduce_b):
        return 0 if a == b else 1
    csub(prime_pow.tmp_ccmp_a, a, b, prec, prime_pow)
    coeffs = prime_pow.tmp_ccmp_a._coeffs
    cdef long i, coeff_prec, break_pt
    if prime_pow.e == 1:
        for i in range(prime_pow.tmp_ccmp_a.degree()+1):
            if coeffs[i] and coeffs[i].valuation() < prec:
                return 1
    else:
        coeff_prec = prec / prime_pow.e + 1
        break_pt = prec % prime_pow.e
        for i in range(len(coeffs)):
            if coeffs[i] and (i < break_pt and coeffs[i].valuation() < coeff_prec or
                              i >= break_pt and coeffs[i].valuation() < coeff_prec - 1):
                return 1
    return 0

cdef inline long cremove(celement out, celement a, long prec, PowComputer_ prime_pow, bint reduce_relative=False) except -1:
    """
    Extract the maximum power of the uniformizer dividing this element.

    INPUT:

    - ``out`` -- a ``celement`` to store the unit part

    - ``a`` -- the ``celement`` whose valuation and unit are desired

    - ``prec`` -- a ``long``, the return value if ``a`` is zero

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    - ``reduce_relative`` -- a bint: whether the final result
      should be reduced at precision ``prec`` (case ``False``)
      or ``prec - valuation`` (case ``True``)

    OUTPUT:

    The number of times the uniformizer divides ``a``, or ``prec`` if ``a`` is
    zero.
    """
    if a == 0:
        return prec
    cdef long v = cvaluation(a, prec, prime_pow)
    if reduce_relative:
        cshift_notrunc(out, a, -v, prec-v, prime_pow, True)
    else:
        cshift_notrunc(out, a, -v, prec, prime_pow, True)
    return v

cdef inline bint cisunit(celement a, PowComputer_ prime_pow) except -1:
    r"""
    Return whether ``a`` has valuation zero.

    INPUT:

    - ``a`` -- the ``celement`` to test

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    return cvaluation(a, 1, prime_pow) == 0

cdef inline int cneg(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to the negative of ``a``.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- a ``celement`` to store the negation

    - ``a`` -- a ``celement`` to be negated

    - ``prec`` -- a ``long``, the precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    cdef celement ma = -a
    if ma is a:
        out._coeffs = ma._coeffs[:]
    else:
        out._coeffs = ma._coeffs

cdef inline int cadd(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to the sum of ``a + b``.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- a ``celement`` to store the sum

    - ``a`` -- a ``celement``, the first summand

    - ``b`` -- a ``celement``, the second summand

    - ``prec`` -- a ``long``, the precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    cdef celement sm = a + b
    if sm is a or sm is b:
        out._coeffs = sm._coeffs[:]
    else:
        out._coeffs = sm._coeffs

cdef inline int csub(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to the difference of ``a - b``.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- a ``celement`` to store the difference

    - ``a`` -- a ``celement``, the first input

    - ``b`` -- a ``celement``, the second input

    - ``prec`` -- a ``long``, the precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    cdef celement df = a - b
    if df is a or df is b:
        out._coeffs = df._coeffs[:]
    else:
        out._coeffs = df._coeffs

cdef inline int cmul(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to the product of ``a*b``.

    Note that no reduction is performed.

    INPUT:

    - ``out`` -- a ``celement`` to store the product

    - ``a`` -- a ``celement``, the first input

    - ``b`` -- a ``celement``, the second input

    - ``prec`` -- a ``long``, the precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    cdef celement pd = a*b
    if pd is a or pd is b:
        out._coeffs = pd._coeffs[:]
    else:
        out._coeffs = pd._coeffs

cdef inline int csetone(celement out, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to 1.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 1

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    out._coeffs = [prime_pow.base_ring(1)]

cdef inline int csetzero(celement out, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to 0.

    INPUT:

    - ``out`` -- the ``celement`` in which to store 0

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    out._coeffs = []

cdef inline bint cisone(celement a, PowComputer_ prime_pow) except -1:
    r"""
    Return whether ``a`` is equal to 1.

    INPUT:

    - ``a`` -- the element to test

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    return a.is_one()

cdef inline bint ciszero(celement a, PowComputer_ prime_pow) except -1:
    r"""
    Return whether ``a`` is equal to 0.

    INPUT:

    - ``a`` -- the element to test

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    return a.is_zero()

cdef inline int ccopy(celement out, celement a, PowComputer_ prime_pow) except -1:
    r"""
    Copy ``a`` into ``out``.

    INPUT:

    - ``out`` -- the ``celement`` to store the result

    - ``a`` -- the element from which to copy

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    out._coeffs = a._coeffs[:]

cdef inline cpickle(celement a, PowComputer_ prime_pow):
    r"""
    Return a representation of ``a`` for pickling.

    INPUT:

    - ``a`` -- the element to pickle

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    return a._coeffs

cdef inline int cunpickle(celement out, x, PowComputer_ prime_pow) except -1:
    r"""
    Reconstruct from the output of :meth:`cpickle`.

    INPUT:

    - ``out`` -- the ``celement`` in which to store the result

    - ``x`` -- the result of `meth`:cpickle

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    out._coeffs = x

cdef inline long chash(celement a, long ordp, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Return a hash value for ``a``.

    INPUT:

    - ``a`` -- the element to hash

    - ``ordp`` -- the valuation as a ``long``

    - ``prec`` -- the precision as a ``long``

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    if ciszero(a, prime_pow):
        return 0

    return hash((a._cache_key(), ordp, prec))

cdef int cconv(celement out, x, long prec, long valshift, PowComputer_ prime_pow) except -2:
    r"""
    Convert ``x`` to a `p`-adic element and store it in ``out``.

    INPUT:

    - ``out`` -- a ``celement`` to store the output

    - ``x`` -- a Sage element that can be converted to a `p`-adic element

    - ``prec`` -- a ``long``, giving the precision desired: absolute if `valshift =
      0`, relative if `valshift > 0`

    - ``valshift`` -- the power of the uniformizer to divide by before storing
      the result in ``out``

    - ``prime_pow`` -- a ``PowComputer`` for the ring
    """
    cdef celement xx, shift
    if valshift < 0:
        if not isinstance(x, list):
            raise NotImplementedError
        # Since our polynomials are stored with ring coefficients, we need to shift the entries of x
        baseshift = -valshift // prime_pow.e
        shifted_x = [c << baseshift for c in x]
        R = prime_pow.poly_ring
        xx = R(shifted_x)
        shift = R([R.base_ring().uniformizer_pow(baseshift)])
        cshift_notrunc(shift, shift, valshift, prec, prime_pow, False)
    else:
        xx = prime_pow.poly_ring(x)
    if xx is x:
        out._coeffs = xx._coeffs[:]
    else:
        out._coeffs = xx._coeffs
    if valshift > 0:
        cshift_notrunc(out, out, -valshift, prec, prime_pow, True)
    elif valshift == 0:
        creduce(out, out, prec, prime_pow)
    elif valshift < 0:
        cdivunit(out, out, shift, prec, prime_pow)
        creduce(out, out, prec, prime_pow)

cdef inline long cconv_mpz_t(celement out, mpz_t x, long prec, bint absolute, PowComputer_ prime_pow) except -2:
    r"""
    Set ``out`` to the integer ``x``.

    This function provides a fast pathway for conversion of integers that
    doesn't require precomputation of the valuation.

    INPUT:

    - ``out`` -- a ``celement`` to store the output

    - ``x`` -- n ``mpz_t`` giving the integer to be converted

    - ``prec`` -- a ``long``, giving the precision desired: absolute or relative
      depending on the ``absolute`` input

    - ``absolute`` -- if ``False`` then extracts the valuation and returns it,
      storing the unit in ``out``; if ``True`` then just reduces ``x`` modulo
      the precision.

    - ``prime_pow`` -- a ``PowComputer`` for the ring

    OUTPUT:

    If ``absolute`` is ``False``, the valuation that was extracted (``maxordp``
    when `x = 0`).
    """
    cdef long valuation = maxordp

    cdef Integer n = PY_NEW(Integer)
    mpz_set(n.value, x)

    if n:
        out._coeffs = [prime_pow.base_ring(n)]
        if not absolute:
            valuation = cremove(out, out, prec, prime_pow)
        creduce(out, out, prec, prime_pow)
    else:
        out._coeffs = []

    return valuation

cdef inline int cconv_mpz_t_out(mpz_t out, celement x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to the integer approximated by ``x``.

    INPUT:

    - ``out`` -- stores the resulting integer as an integer between 0 and
      `p^{\mathrm{prec} + \mathrm{valshift}}`

    - ``x`` -- a ``celement`` giving the underlying `p`-adic element

    - ``valshift`` -- a ``long`` giving the power of the uniformizer to shift ``x`` by

    -` ``prec`` -- a ``long``, the precision of ``x`` and ``out``

    - ``prime_pow`` -- a ``PowComputer`` for the ring
    """
    cdef Integer n

    if valshift:
        cshift_notrunc(prime_pow.powhelper_cconv_out, x, -valshift, prec, prime_pow, True)
    else:
        prime_pow.powhelper_cconv_out = x

    if len(prime_pow.powhelper_cconv_out._coeffs) == 0:
        mpz_set_ui(out, 0)
    elif len(prime_pow.powhelper_cconv_out._coeffs) == 1:
        # recursively let the underlying polynomial convert the constant
        # coefficient to an integer (if possible)
        n = ZZ(prime_pow.powhelper_cconv_out._coeffs[0])
        mpz_set(out, n.value)
    else:
        raise ValueError("cannot convert to integer")

cdef inline long cconv_mpq_t(celement out, mpq_t x, long prec, bint absolute, PowComputer_ prime_pow) except? -10000:
    r"""
    Set ``out`` to the rational ``x``.

    A fast pathway for conversion of rationals that doesn't require
    precomputation of the valuation.

    INPUT:

    - ``out`` -- a ``celement`` to store the output

    - ``x`` -- an ``mpq_t`` giving the rational to be converted

    - ``prec`` -- a ``long``, giving the precision desired: absolute or relative
      depending on ``absolute``

    - ``absolute`` -- if ``False`` then extracts the valuation and returns it,
      storing the unit in ``out``; if ``True`` then just reduces ``x`` modulo the
      precision.

    - ``prime_pow`` -- a ``PowComputer`` for the ring

    OUTPUT:

    If ``absolute`` is ``False``, the valuation that was extracted (``maxordp``
    when `x = 0`)
    """
    cdef Rational r = PY_NEW(Rational)
    mpq_set(r.value, x)
    out._coeffs = [prime_pow.base_ring(r)]

    if not absolute:
        return cremove(out, out, prec, prime_pow)
    creduce(out, out, prec, prime_pow)

cdef inline int cconv_mpq_t_out(mpq_t out, celement x, long valshift, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Set ``out`` to the rational approximated by ``x``.

    INPUT:

    - ``out`` -- set to a rational approximating the input.  Currently uses
      rational reconstruction but may change in the future to use a more naive
      method.

    - ``x`` -- a ``celement``, the element to convert

    - ``valshift`` -- a ``long`` giving the power of the uniformizer to shift `x`
      by before converting

    - ``prec`` -- a long, the precision of ``x``; ignored

    - ``prime_pow`` -- a ``PowComputer`` for the ring
    """
    cdef Rational c

    if valshift:
        cshift_notrunc(prime_pow.powhelper_cconv_out, x, -valshift, prec, prime_pow, True)
    else:
        prime_pow.powhelper_cconv_out = x

    if len(prime_pow.powhelper_cconv_out._coeffs) == 0:
        mpq_set_ui(out, 0, 1)
    elif len(prime_pow.powhelper_cconv_out._coeffs) == 1:
        # recursively let the underlying polynomial convert the constant
        # coefficient to a rational (if possible)
        c = QQ(prime_pow.powhelper_cconv_out._coeffs[0])
        mpq_set(out, c.value)
    else:
        raise ValueError("cannot convert to rational")
