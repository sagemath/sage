r"""
The functions in this file are used in creating new p-adic elements.

When creating a p-adic element, the user can specify that the absolute
precision be bounded and/or that the relative precision be bounded.
Moreover, different p-adic parents impose their own bounds on the
relative or absolute precision of their elements.  The precision
determines to what power of `p` the defining data will be reduced, but
the valuation of the resulting element needs to be determined before
the element is created.  Moreover, some defining data can impose their
own precision bounds on the result.

AUTHORS:

- David Roe (2012-03-01)
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.int cimport *
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
import sage.rings.finite_rings.integer_mod
from sage.libs.pari.types cimport *
from sage.libs.pari.gen cimport gen as pari_gen
from sage.rings.infinity import infinity

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1
# The following Integer is used so that the functions here don't need to initialize an mpz_t.
cdef Integer tmp = PY_NEW(Integer)


cdef long get_ordp(x, PowComputer_class prime_pow) except? -10000:
    """
    This function determines the valuation of the `p`-adic element
    that will result from the given data ``x``.

    Note that the valuation can differ depending on the ring: if the
    new `p`-adic element is being created in a ring with ramification
    then the valuation will be larger.

    Some preprocessing is done in the initialization methods of the element, so
    the type of x is restricted.

    Also note that for some kinds of elements conversion to and from
    Integers and Rationals is done in a custom morphism rather than
    through this function.

    INPUT:

    - ``x`` -- data defining a new p-adic element: a Python int, an
      Integer, Rational, an element of Zp or Qp with the same prime, a
      PARI p-adic element, or an IntegerMod.

    - a PowComputer associated to a `p`-adic ring, which determines
      the prime and the ramification degree.

    OUTPUT:

    - a long, giving the valuation of the resulting `p`-adic element.
      If the input is zero, returns ``maxordp``
    """
    cdef long k, n, p
    cdef Integer value
    cdef GEN pari_tmp
    if PyInt_Check(x):
        if x == 0:
            return maxordp
        try:
            n = PyInt_AsLong(x)
        except OverflowError:
            x = Integer(x)
        else:
            if mpz_fits_slong_p(prime_pow.prime.value) == 0:
                # x is not divisible by p
                return 0
            p = mpz_get_si(prime_pow.prime.value)
            k = 0
            while n % p == 0:
                k += 1
                n = n / p
    if isinstance(x, Integer):
        if mpz_sgn((<Integer>x).value) == 0:
            return maxordp
        k = mpz_remove(tmp.value, (<Integer>x).value, prime_pow.prime.value)
    elif isinstance(x, Rational):
        if mpq_sgn((<Rational>x).value) == 0:
            return maxordp
        k = mpz_remove(tmp.value, mpq_numref((<Rational>x).value), prime_pow.prime.value)
        if k == 0:
            k = -mpz_remove(tmp.value, mpq_denref((<Rational>x).value), prime_pow.prime.value)
    elif isinstance(x, pAdicGenericElement) and (<pAdicGenericElement>x)._is_base_elt(prime_pow.prime):
        k = (<pAdicGenericElement>x).valuation_c()
    elif isinstance(x, pari_gen):
        pari_tmp = (<pari_gen>x).g
        if typ(pari_tmp) == t_PADIC:
            k = valp(pari_tmp)
        else: # t_INT and t_FRAC were converted before this function
            raise TypeError("unsupported coercion from pari: only p-adics, integers and rationals allowed")
    elif sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
        value = <Integer>x.lift()
        if mpz_sgn(value.value) == 0:
            return maxordp
        k = mpz_remove(tmp.value, value.value, prime_pow.prime.value)
    else:
        raise RuntimeError
    # Should check for overflow
    return k * prime_pow.e

cdef long get_preccap(x, PowComputer_class prime_pow) except? -10000:
    """
    This function determines the maximum absolute precision possible
    for an element created from the given data ``x``.

    Note that the valuation can differ depending on the ring: if the
    new `p`-adic element is being created in a ring with ramification
    then the precision bound will be larger.

    Some preprocessing is done in the initialization methods of the element, so
    the type of x is restricted.

    Also note that for some kinds of elements conversion to and from
    Integers and Rationals is done in a custom morphism rather than
    through this function.

    INPUT:

    - ``x`` -- data defining a new p-adic element: an Integer,
      Rational, an element of Zp or Qp with the same prime, a PARI
      p-adic element, or an IntegerMod.
    - ``prime_pow`` -- the PowComputer for the ring into which ``x``
      is being converted.  This is used to determine the prime and the
      ramification degree.

    OUTPUT:

    - a long, giving the absolute precision modulo which the input is
      defined.  If the input is exact, returns ``maxordp``
    """
    cdef long k
    cdef Integer prec
    cdef GEN pari_tmp
    if PyInt_Check(x) or isinstance(x, Integer) or isinstance(x, Rational):
        return maxordp
    elif isinstance(x, pAdicGenericElement) and (<pAdicGenericElement>x)._is_base_elt(prime_pow.prime):
        if (<pAdicGenericElement>x)._is_exact_zero():
            return maxordp
        prec = <Integer>x.precision_absolute()
        k = mpz_get_si(prec.value)
    elif isinstance(x, pari_gen):
        pari_tmp = (<pari_gen>x).g
        # since get_ordp has been called typ(x.g) == t_PADIC
        k = valp(pari_tmp) + precp(pari_tmp)
    elif sage.rings.finite_rings.integer_mod.is_IntegerMod(x):
        k = mpz_remove(tmp.value, (<Integer>x.modulus()).value, prime_pow.prime.value)
        if mpz_cmp_ui(tmp.value, 1) != 0:
            raise TypeError("cannot coerce from the given integer mod ring (not a power of the same prime)")
    else:
        raise RuntimeError
    return k * prime_pow.e

cdef long comb_prec(iprec, long prec) except? -10000:
    """
    This function returns the minimum of iprec and prec.

    INPUT:

    - ``iprec`` -- either infinity, an Integer, a Python int or
      something that can be converted to an Integer.

    - ``prec`` -- a long.
    """
    if iprec is infinity: return prec
    cdef Integer intprec
    if isinstance(iprec, Integer):
        intprec = <Integer>iprec
        if mpz_cmp_si(intprec.value, prec) >= 0:
            return prec
        if mpz_fits_slong_p(intprec.value) == 0:
            raise OverflowError("precision overflow")
        return mpz_get_si(intprec.value)
    if isinstance(iprec, int):
        return min(PyInt_AS_LONG(iprec), prec)
    return comb_prec(Integer(iprec), prec)

cdef int _process_args_and_kwds(long *aprec, long *rprec, args, kwds, bint absolute, PowComputer_class prime_pow) except -1:
    """
    This function obtains values for absprec and relprec from a
    combination of positional and keyword arguments.

    When creating a p-adic element, the user can pass in two arguments: ``absprec`` and ``relprec``.
    In implementing morphisms to speed up conversion from Integers and Rationals,
    we need to determine ``absprec`` and ``relprec`` from the ``args`` and ``kwds`` arguments of
    ``_call_with_args``.  This function collects the code to do so.

    INPUT:

    - ``args`` -- a tuple of positional arguments (at most two)

    - ``kwds`` -- a dictionary of keyword arguments (only
      ``'relprec'`` and ``'absprec'`` are used)

    - ``absolute`` -- (boolean) True if the precision cap of the ring
      is a cap on absolute precision, False if a cap on relative
      precision.

    - ``prime_pow`` -- a
      :class:`sage.rings.padics.pow_computer.PowComputer_class`
      instance

    OUTPUT:

    - ``aprec`` -- (first argument) the maximum absolute precision of
      the resulting element.

    - ``rprec`` -- (second argument) the maximum relative precision of
      the resulting element.

    - error status
    """
    if "empty" in kwds:
        # For backward compatibility
        aprec[0] = 0
        rprec[0] = 0
        return 0
    if len(args) > 2:
        raise TypeError("too many positional arguments")
    if len(args) == 2:
        if "relprec" in kwds:
            raise TypeError("_call_with_args() got multiple values for keyword argument 'relprec'")
        relprec = args[1]
    else:
        relprec = kwds.get("relprec",infinity)
    if len(args) >= 1:
        if "absprec" in kwds:
            raise TypeError("_call_with_args() got multiple values for keyword argument 'absprec'")
        absprec = args[0]
    else:
        absprec = kwds.get("absprec",infinity)
    if absolute:
        aprec[0] = comb_prec(absprec, prime_pow.prec_cap)
        rprec[0] = comb_prec(relprec, maxordp)
    else:
        rprec[0] = comb_prec(relprec, prime_pow.prec_cap)
        aprec[0] = comb_prec(absprec, maxordp)
