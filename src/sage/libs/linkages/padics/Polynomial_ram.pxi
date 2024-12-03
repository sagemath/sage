r"""
This linkage file implements the padics API for ramified extensions using Sage
Polynomials.

It contains the bits that are specific for ramified extensions. Everything that
is independent of ramification is in Polynomial_shared.pxi.

.. NOTE::

    There are no doctests in this file since the functions here cannot be
    called directly from Python. Testing of this function is necessarily
    indirect and mostly done through arithmetic black-box tests that are part
    of the test suites of the `p`-adic parents.

AUTHORS:

- David Roe, Julian Rüth (2017-06-11): initial version
"""
# ****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#                     2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.integer cimport Integer
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport *

include "sage/libs/linkages/padics/Polynomial_shared.pxi"

cdef inline bint creduce(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Reduce ``a`` modulo a power of the maximal ideal.

    INPUT:

    - ``out`` -- a ``celement`` to store the reduction

    - ``a`` -- the ``celement`` to be reduced

    - ``prec`` -- a ``long``, the precision to reduce modulo

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    OUTPUT: ``True`` if the reduction is zero, ``False`` otherwise
    """
    cdef celement ared = a % prime_pow.modulus
    if ared is a and out is not a:
        out._coeffs = ared._coeffs[:]
    else:
        out._coeffs = ared._coeffs
    cdef long coeff_prec = prec / prime_pow.e + 1
    cdef long break_pt = prec % prime_pow.e
    for i in range(len(out._coeffs)):
        if i == break_pt:
            coeff_prec -= 1
        out._coeffs[i] = out._coeffs[i].add_bigoh(coeff_prec)
    out._normalize()
    return out == 0

cdef inline bint creduce_small(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Reduce ``a`` modulo a power of the maximal ideal.

    Similar to ``creduce`` but this function assumes that at most one
    addition/subtraction has happened on reduced inputs.  For integral inputs
    this translates to the assumption that `-p^\mathrm{prec} < a < 2p^\mathrm{prec}`.

    INPUT:

    - ``out`` -- a ``celement`` to store the reduction

    - ``a`` -- the ``celement`` to be reduced

    - ``prec`` -- a ``long``, the precision to reduce modulo

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    OUTPUT: ``True`` if the reduction is zero, ``False`` otherwise
    """
    return creduce(out, a, prec, prime_pow)

cdef inline long cvaluation(celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Return the maximum power of the uniformizer dividing ``a``.

    This function differs from :meth:`cremove` in that the unit is discarded.

    INPUT:

    - ``a`` -- the element whose valuation is desired

    - ``prec`` -- a ``long``; if the valuation of ``a`` exceeds ``prec``, this
      function returns ``prec``. In particular, ``prec`` is returned if ``a``
      is zero.

    - ``prec`` -- a ``long``, the return value if ``a`` is zero

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    OUTPUT:

    The number of times the uniformizer divides ``a``, or ``prec`` if that is
    higher.
    """
    C = a._coeffs
    if not C:
        return prec
    cdef long ret = maxordp

    for i,c in enumerate(C):
        ret = min(ret, c.valuation()*prime_pow.e + i)

    return ret

cdef inline int cshift(celement shifted, celement rem, celement a, long n, long prec, PowComputer_ prime_pow, bint reduce_afterward) except -1:
    r"""
    Multiply ``a`` by an ``n``-th power of the uniformizer.

    This function shifts the π-adic expansion of ``a`` by ``n``, i.e., it
    multiplies ``a`` by the `n`-th power of the uniformizer and drops any terms
    with negative powers of the uniformizer in the π-adic expansion.

    INPUT:

    - ``shifted`` -- a ``celement`` to store the product with π^n

    - ``rem`` -- a ``celement`` to store the remainder, when `n < 0`

    - ``a`` -- the ``celement`` to shift

    - ``n`` -- a ``long``, the amount to shift by

    - ``prec`` -- a ``long``, a precision modulo which to reduce if
      ``reduce_afterward`` is set

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    - ``reduce_afterward`` -- whether to :meth:`creduce` ``shifted`` before
      returning
    """
    cdef long v
    notrunc = False
    if n >= 0:
        notrunc = True
    else:
        v = cvaluation(a, prec, prime_pow)
        notrunc = (v >= -n)
    if notrunc:
        rem._coeffs = []
        return cshift_notrunc(shifted, a, n, prec, prime_pow, reduce_afterward)
    if v > 0:
        b = prime_pow.poly_ring(0)
        cshift_notrunc(b, a, -v, prec, prime_pow, False)
        a = b
        n += v
    # Because the π-adic and p-adic expansions are different,
    # there isn't an obvious way to translate shift-with-truncation
    # in the π-adic expansion in terms of shifting by p in the base.
    # We thus have to compute the first few terms of the π-adic expansion
    # and subtract them.
    _rem = prime_pow.poly_ring(0)
    for i in range(-n):
        rep = a[0]._modp_rep(return_list=False)
        _rem += rep * prime_pow.uniformizer_pow(i)
        a -= rep
        a *= prime_pow.px_pow(1)
        a %= prime_pow.modulus()
        a = a.map_coefficients(lambda c: c>>1)
    cshift_notrunc(rem, _rem, v, prime_pow.ram_prec_cap, prime_pow, True)
    if reduce_afterward:
        creduce(shifted, a, prec, prime_pow)
    else:
        shifted._coeffs = a._coeffs[:]

cdef inline int cshift_notrunc(celement out, celement a, long n, long prec, PowComputer_ prime_pow, bint reduce_afterward) except -1:
    r"""
    Multiply ``a`` with an ``n``-th power of the uniformizer.

    This method is identical to :meth:`cshift` but assumes that the valuation
    of ``a`` is at least ``-n``.

    INPUT:

    - ``out`` -- a ``celement`` to store the result

    - ``a`` -- the ``celement`` to shift

    - ``n`` -- a ``long``, the amount to shift by

    - ``prec`` -- a ``long``, a precision modulo which to reduce if
      ``reduce_afterward`` is set

    - ``prime_pow`` -- the ``PowComputer`` for the ring

    - ``reduce_afterward`` -- whether to :meth:`creduce` ``out`` before
      returning
    """
    cdef long q, r

    if n > 0:
        a *= prime_pow.uniformizer_pow(n)
    elif n < 0:
        q = -n / prime_pow.e # ≥ 0
        r = -n % prime_pow.e # ≥ 0
        # As 0 > n = -q*e - r, π^n = p^-q * (p/π^e)^q * π^-r
        if q:
            # Multiply by p^-q.  Important to do first for fixed-mod.
            a = a.map_coefficients(lambda c: c>>q)
            # and adjust to the π-adic expansion by multiplying by (p/π^e)^q.
            a *= prime_pow.pxe_pow(q)
        if r:
            # Note that in fixed mod this operation can cause errors in the last
            # few digits, since we must divide by p afterward.
            # Multiply by (p/π^r).
            a *= prime_pow.px_pow(r)
            a %= prime_pow.modulus()
            # Divide by p.  This must happen afterward since not all
            # coefficients were multiples of p before.
            a = a.map_coefficients(lambda c: c>>1)

    if reduce_afterward:
        creduce(out, a, prec, prime_pow)
    else:
        out._coeffs = a._coeffs[:]

cdef inline int cinvert(celement out, celement a, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Compute the inverse of ``a``.

    INPUT:

    - ``out`` -- a ``celement`` to store the inverse

    - ``a`` -- a ``celement``, the element to be inverted

    - ``prec`` -- a ``long``, ``out`` is reduced to this precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    out._coeffs = prime_pow.invert(a, prec)._coeffs
    creduce(out, out, prec, prime_pow)

cdef inline int cdivunit(celement out, celement a, celement b, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Divide ``a`` by ``b``.

    This function computes ``a*(1/b)`` where the inverse of ``b`` is determined
    to precision ``prec``. No reduction is performed after the product.

    INPUT:

    - ``out`` -- a ``celement`` to store the quotient

    - ``a`` -- a ``celement``, the dividend

    - ``b`` -- a ``celement``, the divisor, an element of valuation zero

    - ``prec`` -- a ``long``, the precision to which the inverse of ``b`` is
      determined

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    binv = prime_pow.invert(b, prec)
    cmul(out, a, binv, prec, prime_pow)

cdef inline int cpow(celement out, celement a, mpz_t n, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Raise ``a`` to the ``n``-th power.

    INPUT:

    - ``out`` -- a ``celement`` in which to store the result

    - ``a`` -- a ``celement``, the base

    - ``n`` -- an ``mpz_t``, the exponent

    - ``prec`` -- a ``long``, the working absolute precision

    - ``prime_pow`` -- the ``PowComputer`` for the ring
    """
    cdef Integer zn = PY_NEW(Integer)
    mpz_set(zn.value, n)

    csetone(out, prime_pow)
    if zn == 0:
        return 0

    cmul(out, out, a, prec, prime_pow)

    for digit in zn.binary()[1:]:
        cmul(out, out, out, prec, prime_pow)
        if digit == '1':
            cmul(out, out, a, prec, prime_pow)
        # we should probably not creduce that frequently to increase the performance
        creduce(out, out, prec, prime_pow)

# The element is filled in for zero in the p-adic expansion if necessary.
# This WON'T work if the absolute inertia degree is 1.
_expansion_zero = []

# the expansion_mode enum is defined in padic_template_element_header.pxi
cdef inline cexpansion_next(celement value, expansion_mode mode, long curpower, PowComputer_ prime_pow):
    if mode == teichmuller_mode:
        raise NotImplementedError
    # This is not very efficient, but there's no clear better way.
    # We assume this is only called on two-step extensions (for more general
    # extensions, convert to the absolute field).
    R = value.base_ring()
    p = R.prime()
    if R.absolute_degree() == 1:
        raise NotImplementedError("Absolute extensions using Sage polynomials not completely supported")
    if R.base_ring().absolute_degree() != 1:
        raise TypeError("cexpansion only allowed on towers of height 2")
    ans = []
    p2 = (p-1)//2
    # the following is specific to the ramified over unramified case.
    modp_rep, term = value[0]._modp_rep(mode == smallest_mode)
    if term:
        value._coeffs[0] -= modp_rep
    cshift_notrunc(value, value, -1, curpower, prime_pow, False)
    return term

cdef inline cexpansion_getitem(celement value, long m, PowComputer_ prime_pow):
    """
    Return the `m`-th `p`-adic digit in the ``simple_mode`` expansion.

    INPUT:

    - ``value`` -- the `p`-adic element whose expansion is desired
    - ``m`` -- nonnegative integer; which entry in the `p`-adic expansion to return
    - ``prime_pow`` -- a ``PowComputer`` holding `p`-adic data
    """
    R = value.base_ring()
    p = R.prime()
    while m >= 0:
        modp_rep, term = value[0]._modp_rep()
        if m:
            if len(value._coeffs):
                value._coeffs[0] -= modp_rep
            else:
                value._coeffs.append(-modp_rep)
            cshift_notrunc(value, value, -1, 1, prime_pow, False)
        m -= 1
    return term

cdef int cteichmuller(celement out, celement value, long prec, PowComputer_ prime_pow) except -1:
    r"""
    Compute a Teichmüller representative congruent to ``value``.

    INPUT:

    - ``out`` -- a ``celement`` which is set to a `q-1`-th root of unity
      congruent to ``value`` modulo π; or 0 if `a \equiv 0 \pmod{π}`

    - ``value`` -- n ``celement``, the element mod π to lift

    - ``prec`` -- a ``long``, the precision to which to lift

    - ``prime_pow`` -- the ``PowComputer`` of the ring
    """
    if value[0].valuation() > 0:
        out._coeffs = []
    else:
        out._coeffs = [value[0].parent().teichmuller(value[0])]

cdef list ccoefficients(celement x, long valshift, long prec, PowComputer_ prime_pow):
    """
    Return a list of coefficients, as elements that can be converted into the base ring.

    INPUT:

    - ``x`` -- a ``celement`` giving the underlying `p`-adic element, or possibly its unit part
    - ``valshift`` -- a long giving the power of the uniformizer to shift `x` by
    - ``prec`` -- a long, the (relative) precision desired, used in rational reconstruction
    - ``prime_pow`` -- the Powcomputer of the ring
    """
    if valshift == 0:
        return x.list()
    else:
        cshift_notrunc(prime_pow.tmp_ccoeffs, x, valshift, valshift+prec, prime_pow, True)
        return prime_pow.tmp_ccoeffs.list()
