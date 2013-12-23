"""
`p`-Adic ``ZZ_pX`` CA Element

This file implements elements of eisenstein and unramified extensions
of ``Zp`` with capped absolute precision.

For the parent class see padic_extension_leaves.pyx.

The underlying implementation is through NTL's ``ZZ_pX`` class.  Each
element contains the following data:

- ``absprec`` (long) -- An integer giving the precision to which this
  element is defined.  This is the power of the uniformizer modulo
  which the element is well defined.

- ``value`` (``ZZ_pX_c``) -- An ntl ``ZZ_pX`` storing the value.  The
  variable `x` is the uniformizer in the case of eisenstein extensions.
  This ZZ_pX is created with global ntl modulus determined by absprec.
  Let `a` be absprec and `e` be the ramification index over
  `\mathbb{Q}_p` or `\mathbb{Z}_p`.  Then the modulus is given by
  `p^{ceil(a/e)}`.  Note that all kinds of problems arise if you try
  to mix moduli.  ``ZZ_pX_conv_modulus`` gives a semi-safe way to
  convert between different moduli without having to pass through ZZX
  (see ``sage/libs/ntl/decl.pxi`` and ``c_lib/src/ntl_wrap.cpp``)

- ``prime_pow`` (some subclass of ``PowComputer_ZZ_pX``) -- a class,
  identical among all elements with the same parent, holding common
  data.

  + ``prime_pow.deg`` -- The degree of the extension

  + ``prime_pow.e``   -- The ramification index

  + ``prime_pow.f``   -- The inertia degree

  + ``prime_pow.prec_cap`` -- the unramified precision cap.  For
    eisenstein extensions this is the smallest power of p that is
    zero.

  + ``prime_pow.ram_prec_cap`` -- the ramified precision cap.  For
    eisenstein extensions this will be the smallest power of `x` that
    is indistinguishable from zero.

  + ``prime_pow.pow_ZZ_tmp``, prime_pow.pow_mpz_t_tmp``,
    ``prime_pow.pow_Integer`` -- functions for accessing powers of
    `p`.  The first two return pointers.  See
    ``sage/rings/padics/pow_computer_ext`` for examples and important
    warnings.

  + ``prime_pow.get_context``, ``prime_pow.get_context_capdiv``,
    ``prime_pow.get_top_context`` -- obtain an
    ``ntl_ZZ_pContext_class`` corresponding to `p^n`.  The capdiv
    version divides by ``prime_pow.e`` as appropriate.
    ``top_context`` corresponds to `p^{prec_cap}`.

  + ``prime_pow.restore_context``,
    ``prime_pow.restore_context_capdiv``,
    ``prime_pow.restore_top_context`` -- restores the given context.

  + ``prime_pow.get_modulus``, ``get_modulus_capdiv``,
    ``get_top_modulus`` -- Returns a ``ZZ_pX_Modulus_c*`` pointing to
    a polynomial modulus defined modulo `p^n` (appropriately divided
    by ``prime_pow.e`` in the capdiv case).

EXAMPLES:

An eisenstein extension::

    sage: R = ZpCA(5,5)
    sage: S.<x> = ZZ[]
    sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
    sage: W.<w> = R.ext(f); W
    Eisenstein Extension of 5-adic Ring with capped absolute precision 5 in w defined by (1 + O(5^5))*x^5 + (O(5^5))*x^4 + (3*5^2 + O(5^5))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5))*x^2 + (5^3 + O(5^5))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5))
    sage: z = (1+w)^5; z
    1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
    sage: y = z >> 1; y
    w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
    sage: y.valuation()
    4
    sage: y.precision_relative()
    20
    sage: y.precision_absolute()
    24
    sage: z - (y << 1)
    1 + O(w^25)
    sage: (1/w)^12+w
    w^-12 + w + O(w^12)
    sage: (1/w).parent()
    Eisenstein Extension of 5-adic Field with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (O(5^6))*x^4 + (3*5^2 + O(5^6))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^6))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))

An unramified extension::

    sage: g = x^3 + 3*x + 3
    sage: A.<a> = R.ext(g)
    sage: z = (1+a)^5; z
    (2*a^2 + 4*a) + (3*a^2 + 3*a + 1)*5 + (4*a^2 + 3*a + 4)*5^2 + (4*a^2 + 4*a + 4)*5^3 + (4*a^2 + 4*a + 4)*5^4 + O(5^5)
    sage: z - 1 - 5*a - 10*a^2 - 10*a^3 - 5*a^4 - a^5
    O(5^5)
    sage: y = z >> 1; y
    (3*a^2 + 3*a + 1) + (4*a^2 + 3*a + 4)*5 + (4*a^2 + 4*a + 4)*5^2 + (4*a^2 + 4*a + 4)*5^3 + O(5^4)
    sage: 1/a
    (3*a^2 + 4) + (a^2 + 4)*5 + (3*a^2 + 4)*5^2 + (a^2 + 4)*5^3 + (3*a^2 + 4)*5^4 + O(5^5)
    sage: FFA = A.residue_field()
    sage: a0 = FFA.gen(); A(a0^3)
    (2*a + 2) + O(5)

Different printing modes::

    sage: R = ZpCA(5, print_mode='digits'); S.<x> = ZZ[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; repr(z)
    '...4110403113210310442221311242000111011201102002023303214332011214403232013144001400444441030421100001'
    sage: R = ZpCA(5, print_mode='bars'); S.<x> = ZZ[]; g = x^3 + 3*x + 3; A.<a> = R.ext(g)
    sage: z = (1+a)^5; repr(z)
    '...[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 3, 4]|[1, 3, 3]|[0, 4, 2]'
    sage: R = ZpCA(5, print_mode='terse'); S.<x> = ZZ[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; z
    6 + 95367431640505*w + 25*w^2 + 95367431640560*w^3 + 5*w^4 + O(w^100)
    sage: R = ZpCA(5, print_mode='val-unit'); S.<x> = ZZ[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: y = (1+w)^5 - 1; y
    w^5 * (2090041 + 19073486126901*w + 1258902*w^2 + 674*w^3 + 16785*w^4) + O(w^100)

You can get at the underlying ntl representation::

    sage: z._ntl_rep()
    [6 95367431640505 25 95367431640560 5]
    sage: y._ntl_rep()
    [5 95367431640505 25 95367431640560 5]
    sage: y._ntl_rep_abs()
    ([5 95367431640505 25 95367431640560 5], 0)

NOTES:

If you get an error 'internal error: can't grow this
_ntl_gbigint,' it indicates that moduli are being mixed
inappropriately somewhere.

For example, when calling a function with a ZZ_pX_c as an
argument, it copies.  If the modulus is not set to the modulus of
the ZZ_pX_c, you can get errors.

AUTHORS:

- David Roe (2008-01-01): initial version

- Robert Harron (2011-09): fixes/enhancements

- Julian Rueth (2012-10-15): fixed an initialization bug

"""

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#                     2012 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.libs.pari.all import pari_gen
from sage.interfaces.gp import GpElement
from sage.rings.finite_rings.integer_mod import is_IntegerMod
from sage.rings.all import IntegerModRing
from sage.rings.padics.padic_ext_element cimport pAdicExtElement
from sage.rings.padics.precision_error import PrecisionError

from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX_small_Eis
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX_big_Eis

cdef object infinity
from sage.rings.infinity import infinity

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) -1

cdef class pAdicZZpXCAElement(pAdicZZpXElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, empty = False):
        """
        Creates an element of a capped absolute precision, unramified or eisenstein extension of Zp or Qp.

        INPUT:

        - ``parent`` -- either an ``EisensteinRingCappedAbsolute`` or
          ``UnramifiedRingCappedAbsolute``

        - `x` -- an integer, rational, `p`-adic element, polynomial,
          list, integer_mod, pari int/frac/poly_t/pol_mod, an
          ``ntl_ZZ_pX``, an ``ntl_ZZ``, an ``ntl_ZZ_p``, an
          ``ntl_ZZX``, or something convertible into parent.residue_field()

        - ``absprec`` -- an upper bound on the absolute precision of
          the element created

        - ``relprec`` -- an upper bound on the relative precision of
          the element created

        - ``empty`` -- whether to return after initializing to zero.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1+w)^5; z # indirect doctest
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: W(R(3,3))
            3 + O(w^15)
            sage: W(pari('3 + O(5^3)'))
            3 + O(w^15)
            sage: W(w, 14)
            w + O(w^14)

        TESTS:

        Check that :trac:`13600` is fixed::

            sage: K = W.fraction_field()
            sage: W(K.zero())
            O(w^25)
            sage: W(K.one())
            1 + O(w^25)
            sage: W(K.zero().add_bigoh(3))
            O(w^3)

        Check that :trac:`3865` is fixed:

            sage: W(gp('5 + O(5^2)'))
            w^5 + 2*w^7 + 4*w^9 + O(w^10)
        """
        pAdicZZpXElement.__init__(self, parent)
        cdef long aprec, rprec, ctx_prec
        if empty:
            self.absprec = 0
            return
        self.absprec = -1 # to signal that self is uninitialized
        if absprec is infinity:
            aprec = self.prime_pow.ram_prec_cap
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_sgn((<Integer>absprec).value) < 0:
                aprec = 0
            elif mpz_fits_slong_p((<Integer>absprec).value) == 0:
                aprec = self.prime_pow.ram_prec_cap
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
                if aprec > self.prime_pow.ram_prec_cap:
                    aprec = self.prime_pow.ram_prec_cap
        if relprec is infinity:
            # This might not be the right default
            rprec = self.prime_pow.ram_prec_cap
        else:
            if not PY_TYPE_CHECK(relprec, Integer):
                rprec = Integer(relprec)
            if mpz_cmp_ui((<Integer>relprec).value, aprec) >= 0:
                rprec = self.prime_pow.ram_prec_cap
            elif relprec < 0:
                rprec = 0
            else:
                rprec = mpz_get_si((<Integer>relprec).value)
        cdef mpz_t tmp
        cdef ZZ_c tmp_z
        cdef Py_ssize_t i
        cdef Integer tmp_Int
        if PY_TYPE_CHECK(x, pAdicGenericElement):
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes."
            if x.valuation() < 0:
                raise ValueError, "element has negative valuation"
        if PY_TYPE_CHECK(x, pAdicBaseGenericElement):
            mpz_init(tmp)
            (<pAdicBaseGenericElement>x)._set_mpz_into(tmp)
            if mpz_sgn(tmp) == 0:
                if (<pAdicBaseGenericElement>x)._is_exact_zero():
                    self._set_inexact_zero(aprec)
                    mpz_clear(tmp)
                    return
            ltmp = mpz_get_si((<Integer>x.precision_absolute()).value) * self.prime_pow.e
            if ltmp < aprec:
                aprec = ltmp
            if relprec is infinity:
                self._set_from_mpz_abs(tmp, aprec)
            else:
                self._set_from_mpz_both(tmp, aprec, rprec)
            mpz_clear(tmp)
            return
        if isinstance(x, pari_gen) or isinstance(x, GpElement):
            if isinstance(x, GpElement):
                x = x._pari_()
            if x.type() == "t_PADIC":
                if x.variable() != self.prime_pow.prime:
                    raise TypeError, "Cannot coerce a pari p-adic with the wrong prime."
                ltmp = x.padicprec(self.prime_pow.prime) * self.prime_pow.e
                if ltmp < aprec:
                    aprec = ltmp
                x = x.lift()
            if x.type() == 't_INT':
                x = Integer(x)
            elif x.type() == 't_FRAC':
                x = Rational(x)
            elif x.type() == 't_POLMOD' or x.type == 't_POL':
                # This code doesn't check to see if the primes are the same.
                L = []
                x = x.lift().lift()
                for i from 0 <= i <= x.poldegree():
                    L.append(Integer(x.polcoeff(i)))
                x = L
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers, rationals, polynomials and pol_mods allowed"
        elif is_IntegerMod(x):
            mpz_init(tmp)
            ctx_prec = mpz_remove(tmp, (<Integer>x.modulus()).value, self.prime_pow.prime.value)
            if mpz_cmp_ui(tmp, 1) == 0:
                mpz_clear(tmp)
                x = x.lift()
                if ctx_prec < aprec:
                    aprec = ctx_prec
            else:
                mpz_clear(tmp)
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"
        elif PY_TYPE_CHECK(x, ntl_ZZ_p):
            ctx_prec = ZZ_remove(tmp_z, (<ntl_ZZ>x.modulus()).x, self.prime_pow.pow_ZZ_tmp(1)[0])
            if ZZ_IsOne(tmp_z):
                x = x.lift()
                tmp_Int = PY_NEW(Integer)
                ZZ_to_mpz(&tmp_Int.value, &(<ntl_ZZ>x).x)
                x = tmp_Int
                if ctx_prec < aprec:
                    aprec = ctx_prec
            else:
                raise TypeError, "cannot coerce the given ntl_ZZ_p (modulus not a power of the same prime)"
        elif PY_TYPE_CHECK(x, ntl_ZZ):
            tmp_Int = PY_NEW(Integer)
            ZZ_to_mpz(&tmp_Int.value, &(<ntl_ZZ>x).x)
            x = tmp_Int
        elif isinstance(x, (int, long)):
            x = Integer(x)
        elif x in parent.residue_field():
            # Should only reach here if x is not in F_p
            z = parent.gen()
            poly = x.polynomial().list()
            x = sum([poly[i].lift() * (z ** i) for i in range(len(poly))])
            if 1 < aprec:
                aprec = 1
        cdef pAdicZZpXCAElement _x
        cdef pAdicZZpXCRElement __x
        if PY_TYPE_CHECK(x, Integer):
            if relprec is infinity:
                self._set_from_mpz_abs((<Integer>x).value, aprec)
            else:
                self._set_from_mpz_both((<Integer>x).value, aprec, rprec)
        elif PY_TYPE_CHECK(x, Rational):
            if relprec is infinity:
                self._set_from_mpq_abs((<Rational>x).value, aprec)
            else:
                self._set_from_mpq_both((<Rational>x).value, aprec, rprec)
        elif PY_TYPE_CHECK(x, ntl_ZZ_pX):
            if relprec is infinity:
                self._set_from_ZZ_pX_abs(&(<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c, aprec)
            else:
                self._set_from_ZZ_pX_both(&(<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c, aprec, rprec)
        elif PY_TYPE_CHECK(x, ntl_ZZX):
            if relprec is infinity:
                self._set_from_ZZX_abs((<ntl_ZZX>x).x, aprec)
            else:
                self._set_from_ZZX_both((<ntl_ZZX>x).x, aprec, rprec)
        elif PY_TYPE_CHECK(x, pAdicExtElement):
            if x.parent() is parent:
                _x = <pAdicZZpXCAElement>x
                if _x.absprec < aprec:
                    aprec = _x.absprec
                if rprec < self.prime_pow.ram_prec_cap:
                    self._set_from_ZZ_pX_both(&_x.value, None, aprec, rprec)
                else:
                    self._set(&_x.value, aprec)
            elif x.parent() is parent.fraction_field():
                __x = <pAdicZZpXCRElement>x
                if __x.relprec < 0:
                    __x._normalize()
                if __x._is_exact_zero():
                    self._set_inexact_zero(self.prime_pow.ram_prec_cap)
                elif __x.ordp < 0:
                    raise ValueError, "x has negative valuation"
                elif __x._is_inexact_zero():
                    if __x.ordp <= self.prime_pow.ram_prec_cap:
                        self._set_inexact_zero(__x.ordp)
                    else:
                        self._set_inexact_zero(self.prime_pow.ram_prec_cap)
                else:
                    poly = __x._ntl_rep_abs()[0]
                    if __x.relprec < rprec:
                        rprec = __x.relprec
                    if rprec + __x.ordp < aprec:
                        aprec = rprec + __x.ordp
                    self._set(&(<ntl_ZZ_pX>poly).x, aprec)
            else:
                raise NotImplementedError, "Conversion from different p-adic extensions not yet supported"
        else:
            try:
                x = list(x)
            except TypeError:
                try:
                    x = x.list()
                except AttributeError:
                    raise TypeError, "cannot convert x to a p-adic element"
            self._set_from_list_both(x, aprec, rprec)

    cdef int _set_inexact_zero(self, long absprec) except -1:
        """
        Sets ``self`` to be zero with valuation ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(0,6); z # indirect doctest
            O(w^6)
            sage: z.valuation()
            6
            sage: z.precision_absolute()
            6
            sage: z.precision_relative()
            0
        """
        if absprec == self.absprec:
            ZZ_pX_clear(self.value)
        else:
            self._set_prec_abs(absprec)

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Tests if ``self`` is an inexact zero.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(0)
            sage: z._is_inexact_zero() #indirect doctest
            True
            sage: z = W(0,6)
            sage: z._is_inexact_zero()
            True
        """
        return self.absprec == 0 or ZZ_pX_IsZero(self.value) or self.valuation_c() == self.absprec

    cdef int _set(self, ZZ_pX_c* value, long absprec) except -1:
        """
        Sets ``value`` and ``absprec`` directly.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: F = W.fraction_field()
            sage: z = F(1+w); z # indirect doctest
            1 + w + O(w^25)
            sage: W.precision_cap()
            25
            sage: F.precision_cap()
            25
        """
        self._set_prec_abs(absprec) # restores context
        if self.absprec != 0:
            ZZ_pX_conv_modulus(self.value, value[0], self.prime_pow.get_context_capdiv(absprec).x)

    cdef int _set_from_mpz_abs(self, mpz_t x, long absprec) except -1:
        """
        Sets ``self`` from an ``mpz_t`` with absolute precision
        bounded by ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, absprec = 13) # indirect doctest
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
            sage: W(70, absprec = 4)
            O(w^4)
            sage: W(70, absprec = 0)
            O(w^0)
        """
        self._set_prec_abs(absprec) # restores context
        cdef mpz_t tmp_m
        cdef ZZ_c tmp_z
        if self.absprec != 0:
            mpz_init_set(tmp_m, x)
            mpz_to_ZZ(&tmp_z, &tmp_m)
            mpz_clear(tmp_m)
            ZZ_pX_SetCoeff(self.value, 0, ZZ_to_ZZ_p(tmp_z))

    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, long relprec) except -1:
        """
        Sets ``self`` from an ``mpz_t`` with relative precision
        bounded by ``relprec`` and absolute precision bounded by
        ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, relprec = 3) # indirect doctest
            4*w^5 + 3*w^7 + O(w^8)
            sage: W(70, absprec = 4, relprec = 2)
            O(w^4)
            sage: W(70, absprec = 0, relprec = 3)
            O(w^0)
        """
        if mpz_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return 0
        cdef mpz_t tmp_m
        cdef ZZ_c tmp_z
        cdef long shift
        mpz_init(tmp_m)
        sig_on()
        shift = mpz_remove(tmp_m, x, self.prime_pow.prime.value)
        mpz_set(tmp_m, x)
        sig_off()
        self._set_prec_both_with_ordp(shift * self.prime_pow.e, absprec, relprec)
        mpz_to_ZZ(&tmp_z, &tmp_m)
        mpz_clear(tmp_m)
        if self.absprec != 0:
            ZZ_pX_SetCoeff(self.value, 0, ZZ_to_ZZ_p(tmp_z))

    cdef int _set_from_mpq_abs(self, mpq_t x, long absprec) except -1:
        """
        Sets ``self`` from an ``mpq_t`` with absolute precision bounded by
        ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(70/3, 14); z # indirect doctest
            3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + O(w^14)
            sage: z * 3
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + O(w^14)
            sage: W(70)
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23 + O(w^25)
            sage: W(70/3, absprec = 4)
            O(w^4)
            sage: W(70/3, absprec = 0)
            O(w^0)
        """
        if mpq_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return 0
        if mpz_divisible_p(mpq_denref(x), self.prime_pow.prime.value):
            raise ValueError, "p divides the denominator"
        self._set_prec_abs(absprec) # restores context
        self._set_from_mpq_part2(x)

    cdef int _set_from_mpq_both(self, mpq_t x, long absprec, long relprec) except -1:
        """
        Sets ``self`` from an ``mpq_t`` with relative precision bounded by
        ``relprec`` and absolute precision bounded by ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(70/3, 14); z # indirect doctest
            3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + O(w^14)
            sage: z * 3
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + O(w^14)
            sage: W(70)
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23 + O(w^25)
            sage: W(70/3, absprec = 0, relprec = 1)
            O(w^0)
        """
        if mpq_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return 0
        cdef long num_ordp
        cdef mpz_t num_unit
        if mpz_divisible_p(mpq_denref(x), self.prime_pow.prime.value):
            raise ValueError, "p divides the denominator"
        sig_on()
        mpz_init(num_unit)
        num_ordp = mpz_remove(num_unit, mpq_numref(x), self.prime_pow.prime.value)
        mpz_clear(num_unit)
        sig_off()
        self._set_prec_both_with_ordp(num_ordp * self.prime_pow.e, absprec, relprec) # restores context
        self._set_from_mpq_part2(x)

    cdef int _set_from_mpq_part2(self, mpq_t x) except -1:
        """
        Given that the appropriate context has been restored, sets
        ``self`` from ``x``.

        TESTS::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(QQ(4), 23) # indirect doctest
            4 + O(w^23)
        """
        cdef mpz_t tmp_m
        cdef ZZ_c den_z, num_z
        cdef ZZ_p_c value
        if self.absprec != 0:
            mpz_init_set(tmp_m, mpq_numref(x))
            mpz_to_ZZ(&num_z, &tmp_m)
            mpz_set(tmp_m, mpq_denref(x))
            mpz_to_ZZ(&den_z, &tmp_m)
            mpz_clear(tmp_m)
            ZZ_p_div(value, ZZ_to_ZZ_p(num_z), ZZ_to_ZZ_p(den_z))
            ZZ_pX_SetCoeff(self.value, 0, value)

    cdef int _set_from_ZZX_abs(self, ZZX_c poly, long absprec) except -1:
        """
        Sets ``self`` from a ``ZZX`` with absolute precision bounded by
        ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZX([4,1,16]), absprec = 14); z # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + 4*w^13 + O(w^14)
            sage: z._ntl_rep()
            [4 1 16]
            sage: W(ntl.ZZX([4,1,16]), absprec = 0)
            O(w^0)
        """
        self._set_prec_abs(absprec)
        cdef ZZ_pX_c poly_p
        if self.absprec != 0:
            ZZX_to_ZZ_pX(poly_p, poly)
            self._set_from_ZZ_pX_abs(&poly_p, None, absprec)

    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, long relprec) except -1:
        """
        Sets ``self`` from a ``ZZX`` with relative precision bounded by
        ``relprec`` and absolute precision bounded by ``absprec``.

        EXAMPLES::
            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZX([4,1,16]), relprec = 12); z # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + O(w^12)
            sage: z._ntl_rep()
            [4 1 16]
            sage: W(ntl.ZZX([4,1,16]), absprec = 0, relprec = 4)
            O(w^0)
        """
        self._set_prec_abs(absprec)
        cdef ZZ_pX_c poly_p
        if self.absprec != 0:
            ZZX_to_ZZ_pX(poly_p, poly)
            self._set_from_ZZ_pX_both(&poly_p, None, absprec, relprec)

    cdef int _set_from_ZZ_pX_abs(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec) except -1:
        """
        Sets ``self`` from a ``ZZ_pX`` with absolute precision bounded by ``absprec`` (and by ``ctx``).

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2)); z # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + O(w^10)
            sage: z._ntl_rep()
            [4 1 16]
            sage: W(ntl.ZZ_pX([4,1,16],5^2), absprec = 0)
            O(w^0)
        """
        cdef long ctx_prec = -1
        if ctx is not None:
            ctx_prec = self._check_ZZ_pContext(ctx) * self.prime_pow.e
            if ctx_prec < absprec:
                absprec = ctx_prec
        if ZZ_pX_IsZero(poly[0]):
            self._set_inexact_zero(absprec)
            return 0
        self._set_prec_abs(absprec) # restores context
        if self.absprec != 0:
            ZZ_pX_conv_modulus(self.value, poly[0], self.prime_pow.get_context_capdiv(absprec).x)

    cdef int _set_from_ZZ_pX_both(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec, long relprec) except -1:
        """
        Sets ``self`` from a ``ZZ_pX`` with relative precision bounded by
        ``relprec`` and absolute precision bounded by ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2), absprec = 8, relprec = 12); z # indirect doctest
            4 + w + w^2 + 3*w^7 + O(w^8)
            sage: z._ntl_rep()
            [4 1 16]
            sage: W(ntl.ZZ_pX([4,1,16],5^2), absprec = 0, relprec = 5)
            O(w^0)
        """
        cdef long ctx_prec
        if ctx is not None:
            ctx_prec = self._check_ZZ_pContext(ctx)
            if ctx_prec * self.prime_pow.e < absprec:
                absprec = ctx_prec * self.prime_pow.e
        if ZZ_pX_IsZero(poly[0]):
            self._set_inexact_zero(absprec)
            return 0
        cdef long val, index
        ZZ_pX_min_val_coeff(val, index, poly[0], self.prime_pow.pow_ZZ_tmp(1)[0])
        if self.prime_pow.e == 1:
            self._set_prec_both_with_ordp(val, absprec, relprec) #restores context
        else:
            self._set_prec_both_with_ordp(val * self.prime_pow.e + index, absprec, relprec) # restores context
        if self.absprec != 0:
            ZZ_pX_conv_modulus(self.value, poly[0], self.prime_pow.get_context_capdiv(self.absprec).x)

    cdef bint _set_prec_abs(self, long absprec) except -1:
        """
        Safely sets the absolute precision of self to ``absprec``.

        Returns ``True`` iff ``self.absprec`` was reset.

        Note that this will wipe out anything in ``self.value``.  Be
        careful resetting ``self.value`` directly: if you set it to a
        different modulus, NTL may have problems.  The safest way to
        reset ``self.value`` to a different modulus is::

            self.prime_pow.restore_context_capdiv(self.absprec)
            cdef ZZ_pX_c tmp = self.value
            self._set_prec_abs(new_abs_prec)
            ZZ_pX_conv_modulus(self.value, tmp, self.prime_pow.get_context_capdiv(self.absprec).x)

        You may be able to just set ``self.absprec`` and
        ``ZZ_pX_conv_modulus`` if you're decreasing precision.  I'm not
        sure.

        TESTS::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, 13) # indirect doctest
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
        """
        if absprec < 0:
            raise ValueError, "absprec must be non-negative"
        if self.absprec == absprec:
            return False
        if self.absprec > 0:
            ZZ_pX_destruct(&self.value)
        if absprec > 0:
            self.prime_pow.restore_context_capdiv(absprec)
            ZZ_pX_construct(&self.value)
        self.absprec = absprec
        return True

    cdef bint _set_prec_both(self, long absprec, long relprec) except -1:
        raise TypeError, "use _set_prec_both_with_ord"

    cdef bint _set_prec_both_with_ordp(self, long ordp, long absprec, long relprec) except -1:
        """
        Sets the absolute precision of ``self`` to the minimum of ``absprec``
        and ``ordp + relprec``.

        Note that this will wipe out anything in ``self.value``.  Be
        careful resetting ``self.value`` directly: if you set it to a
        different modulus, NTL may have problems.  The safest way to
        reset ``self.value`` to a different modulus is::

            self.prime_pow.restore_context_capdiv(self.absprec)
            cdef ZZ_pX_c tmp = self.value
            self._set_prec_abs(new_abs_prec)
            ZZ_pX_conv_modulus(self.value, tmp, self.prime_pow.get_context_capdiv(self.relprec).x)

        You may be able to just set ``self.absprec`` and
        ``ZZ_pX_conv_modulus`` if you're decreasing precision.  I'm not
        sure.

        TESTS::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, relprec = 3) # indirect doctest
            4*w^5 + 3*w^7 + O(w^8)
        """
        if absprec <= ordp + relprec:
            self._set_prec_abs(absprec)
        else:
            self._set_prec_abs(ordp + relprec)

    def __dealloc__(self):
        """
        Deallocates ``self.value``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(25/3, relprec = 6)
            sage: del z #indirect doctest
        """
        if self.absprec != 0:
            ZZ_pX_destruct(&self.value)

    cdef pAdicZZpXCAElement _new_c(self, long absprec):
        """
        Returns a new element with the same parent as ``self`` and
        absolute precision ``absprec``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w^5 + 1 # indirect doctest
            1 + w^5 + O(w^25)
        """
        cdef pAdicZZpXCAElement ans = PY_NEW(pAdicZZpXCAElement)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        ans.absprec = absprec
        if absprec > 0:
            self.prime_pow.restore_context_capdiv(absprec)
            ZZ_pX_construct(&ans.value)
        elif absprec < 0:
            raise ValueError, "absprec must be positive"
        return ans

    def __reduce__(self):
        """
        Pickles ``self``.

        EXAMPLES::

            sage: R = Qp(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5 - 1
            sage: loads(dumps(z)) == z
            True
        """
        cdef Integer absprec
        absprec = PY_NEW(Integer)
        mpz_set_si(absprec.value, self.absprec)
        if self.absprec == 0:
            return make_ZZpXCAElement, (self.parent(), None, absprec, 0)
        self.prime_pow.restore_context_capdiv(self.absprec)
        cdef ntl_ZZ_pX holder = PY_NEW(ntl_ZZ_pX)
        holder.c = self.prime_pow.get_context_capdiv(self.absprec)
        holder.x = self.value
        return make_ZZpXCAElement, (self.parent(), holder, absprec, 0)

    cdef int _cmp_units(left, pAdicGenericElement right) except -2:
        """
        For units ``left`` and ``right``, returns 0 if they are equal up to
        the lesser of the two precisions, or 1 if they are not.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w == 1 # indirect doctest
            False
            sage: y = 1 + w + O(w^7)
            sage: z = 1 + w + w^10 + O(w^13)
            sage: y == z
            True
        """
        # This function needs improvement.  In particular, there are a lot of
        # speed improvements to be had, and it should be changed so that it
        # returns 1 only half the time (and -1 the other half) when left and
        # right are not equal.
        cdef pAdicZZpXCAElement diff = <pAdicZZpXCAElement> (left - right)
        if diff._is_inexact_zero():
            return 0
        # for now, just return 1
        return 1

    def __invert__(self):
        """
        Returns the inverse of ``self``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: y = ~z; y # indirect doctest
            1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14 + 3*w^16 + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24 + O(w^25)
            sage: y.parent()
            Eisenstein Extension of 5-adic Field with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (O(5^6))*x^4 + (3*5^2 + O(5^6))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^6))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
            sage: z = z - 1
            sage: ~z
            w^-5 + 4*w^-4 + 4*w^-3 + 4*w^-2 + 2*w^-1 + 1 + w + 4*w^2 + 4*w^3 + 4*w^4 + w^5 + w^6 + w^7 + 4*w^8 + 4*w^9 + 2*w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^14 + O(w^15)
        """
        return self.to_fraction_field()._invert_c_impl()

    cpdef RingElement _invert_c_impl(self):
        """
        Returns the inverse of ``self``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: y = ~z; y # indirect doctest
            1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14 + 3*w^16 + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24 + O(w^25)
            sage: y.parent()
            Eisenstein Extension of 5-adic Field with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (O(5^6))*x^4 + (3*5^2 + O(5^6))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^6))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
            sage: z = z - 1
            sage: ~z
            w^-5 + 4*w^-4 + 4*w^-3 + 4*w^-2 + 2*w^-1 + 1 + w + 4*w^2 + 4*w^3 + 4*w^4 + w^5 + w^6 + w^7 + 4*w^8 + 4*w^9 + 2*w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^14 + O(w^15)
        """
        return self.to_fraction_field()._invert_c_impl()

    cpdef pAdicZZpXCRElement to_fraction_field(self):
        """
        Returns ``self`` cast into the fraction field of ``self.parent()``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5; z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: y = z.to_fraction_field(); y #indirect doctest
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: y.parent()
            Eisenstein Extension of 5-adic Field with capped relative precision 5 in w defined by (1 + O(5^5))*x^5 + (O(5^6))*x^4 + (3*5^2 + O(5^6))*x^3 + (2*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))*x^2 + (5^3 + O(5^6))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + O(5^6))
        """
        cdef pAdicZZpXCRElement ans = PY_NEW(pAdicZZpXCRElement)
        ans._parent = self._parent.fraction_field()
        ans.prime_pow = ans._parent.prime_pow
        ans.ordp = 0
        ans.relprec = -self.absprec
        if self.absprec != 0:
            self.prime_pow.restore_context_capdiv(self.absprec)
            ZZ_pX_construct(&ans.unit)
            ans.unit = self.value
        return ans

    cdef pAdicZZpXCAElement _lshift_c(self, long n):
        """
        Multiplies ``self`` by the uniformizer raised to the power ``n``.  If
        ``n`` is negative, right shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: z << 17 # indirect doctest
            w^17 + w^22 + w^23 + 2*w^24 + O(w^25)
            sage: z << (-1)
            w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
        """
        self._rshift_c(-n)

    def __lshift__(pAdicZZpXCAElement self, shift):
        """
        Multiplies ``self`` by the uniformizer raised to the power ``n``.  If
        ``n`` is negative, right shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: z << 17 # indirect doctest
            w^17 + w^22 + w^23 + 2*w^24 + O(w^25)
            sage: z << (-1)
            w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
        """
        cdef pAdicZZpXCAElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            if mpz_sgn((<Integer>shift).value) > 0:
                ans = self._new_c(self.prime_pow.ram_prec_cap)
            else:
                ans = self._new_c(0)
            return ans
        return self._rshift_c(-mpz_get_si((<Integer>shift).value))

    cdef pAdicZZpXCAElement _rshift_c(self, long n):
        """
        Divides ``self`` by the uniformizer raised to the power ``n``.  If
        parent is not a field, throws away the non-positive part of
        the series expansion.  If ``n`` is negative, left shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpCA(5,5,print_mode='digits')
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: for m in range(26): repr(z >> m) # indirect doctest
            '...4001400444441030421100001'
            '...400140044444103042110000'
            '...40014004444410304211000'
            '...4001400444441030421100'
            '...400140044444103042110'
            '...40014004444410304211'
            '...4001400444441030421'
            '...400140044444103042'
            '...40014004444410304'
            '...4001400444441030'
            '...400140044444103'
            '...40014004444410'
            '...4001400444441'
            '...400140044444'
            '...40014004444'
            '...4001400444'
            '...400140044'
            '...40014004'
            '...4001400'
            '...400140'
            '...40014'
            '...4001'
            '...400'
            '...40'
            '...4'
            '...'
            sage: repr(z >> (-4))
            '...4004444410304211000010000'
        """
        cdef long absprec
        if n == 0:
            return self
        elif n > self.absprec: # we do these checks first in case adding would cause an overflow
            absprec = 0
        elif n < self.absprec - self.prime_pow.ram_prec_cap:
            absprec = self.prime_pow.ram_prec_cap
        else:
            absprec = self.absprec - n
        cdef pAdicZZpXCAElement ans
        if absprec > 0:
            ans = self._new_c(absprec)
            if n > -self.prime_pow.ram_prec_cap: # the result might not be zero.
                if self.prime_pow.e == 1:
                    if n > 0:
                        ZZ_pX_right_pshift(ans.value, self.value, self.prime_pow.pow_ZZ_tmp(n)[0], self.prime_pow.get_context(ans.absprec).x)
                    else:
                        ZZ_pX_left_pshift(ans.value, self.value, self.prime_pow.pow_ZZ_tmp(-n)[0], self.prime_pow.get_context(ans.absprec).x)
                else:
                    self.prime_pow.eis_shift_capdiv(&ans.value, &self.value, n, ans.absprec)
        else:
            ans = self._new_c(0)
        return ans

    def __rshift__(pAdicZZpXCAElement self, shift):
        """
        Divides ``self`` by the uniformizer raised to the power ``n``.  If
        parent is not a field, throws away the non-positive part of
        the series expansion.  If ``n`` is negative, left shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: z >> (6) # indirect doctest
            1 + 2*w + 4*w^2 + 3*w^4 + w^6 + 4*w^7 + 4*w^8 + 4*w^9 + 4*w^10 + 4*w^11 + 4*w^14 + w^15 + 4*w^18 + O(w^19)
            sage: z >> (-4)
            w^4 + w^9 + w^10 + 2*w^11 + 4*w^12 + 3*w^14 + w^16 + 4*w^17 + 4*w^18 + 4*w^19 + 4*w^20 + 4*w^21 + 4*w^24 + O(w^25)
        """
        cdef pAdicZZpXCAElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            if mpz_sgn((<Integer>shift).value) < 0:
                ans = self._new_c(self.prime_pow.ram_prec_cap)
            else:
                ans = self._new_c(0)
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cpdef ModuleElement _neg_(self):
        """
        Returns ``-self``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5; z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: -z # indirect doctest
            4 + 3*w^5 + 4*w^6 + w^7 + w^8 + w^9 + w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^15 + 3*w^16 + w^17 + 2*w^18 + 3*w^19 + 2*w^21 + 4*w^23 + 4*w^24 + O(w^25)
            sage: y = z + (-z); y
            O(w^25)
            sage: -y
            O(w^25)
        """
        cdef pAdicZZpXCAElement ans = self._new_c(self.absprec)
        if self.absprec != 0:
            self.prime_pow.restore_context_capdiv(self.absprec)
            ZZ_pX_negate(ans.value, self.value)
        return ans

#                                             / 1 + \alpha^p \pi_K^{p \lambda}                      mod \mathfrak{p}_K^{p \lambda + 1}   if 1 \le \lambda < \frac{e_K}{p-1}
#        (1 + \alpha \pi^{\lambda})^p \equiv {  1 + (\alpha^p - \epsilon \alpha) \pi_K^{p \lambda}  mod \mathfrak{p}_K^{p \lambda + 1}   if \lambda = \frac{e_K}{p-1}
#                                             \ 1 - \epsilon \alpha \pi_K^{\lambda + e}             mod \mathfrak{p}_K^{\lambda + e + 1} if \lambda > \frac{e_K}{p-1}

    def __pow__(pAdicZZpXCAElement self, _right, m): # m ignored
        r"""
        Computes ``self^right``.

        Note: when right is divisible by `p` then one can get more
        precision than expected.

        Lemma 2.1 (Constructing Class Fields over Local Fields,
        Sebastian Pauli): Let `\alpha` be in `\mathcal{O}_K`.  Let

        ..math ::

            p = -\pi_K^{e_K} \epsilon

        be the factorization of `p` where `\epsilon` is a unit.  Then
        the `p`-th power of `1 + \alpha \pi_K^{\lambda}` satisfies

        ..math ::

            (1 + \alpha \pi^{\lambda})^p \equiv \left{ \begin{array}{lll}
            1 + \alpha^p \pi_K^{p \lambda} & \mod \mathfrak{p}_K^{p \lambda + 1} & \mbox{if $1 \le \lambda < \frac{e_K}{p-1}$} \\
            1 + (\alpha^p - \epsilon \alpha) \pi_K^{p \lambda} & \mod \mathfrak{p}_K^{p \lambda + 1} & \mbox{if $\lambda = \frac{e_K}{p-1}$} \\
            1 - \epsilon \alpha \pi_K^{\lambda + e} & \mod \mathfrak{p}_K^{\lambda + e + 1} & \mbox{if $\lambda > \frac{e_K}{p-1}$}
            \end{array} \right.


        So if right is divisible by `p^k` we can multiply the relative
        precision by `p` until we exceed `e/(p-1)`, then add `e` until
        we have done a total of `k` things: the precision of the
        result can therefore be greater than the precision of self.

        There is also the issue of `p`-adic exponents, and determining
        how the precision of the exponent affects the precision of the
        result.

        In computing `(a + O(\pi^k))^{b + O(p^m)}`, one needs that the
        reduction of `a` mod `\pi` is in the prime field `\mathbb{F}_p` (so
        that the `p^m` power of the Teichmuller part is constant as
        `m` increases).  Given this restriction, we can factor out the
        Teichmuller part and use the above lemma to find the first
        spot where

        ..math ::

            (1 + \alpha \pi^{\lambda})^{p^m}

        differs from 1.  We compare this with the precision bound
        given by computing `(a + O(\pi^k))^b` and take the lesser of
        the two.

        In order to do this we need to compute the valuation of ``(self
        / self.parent().teichmuller(self)) - 1``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: (1 + w)^5 # indirect doctest
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: (1 + w + O(w^19))^5
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + O(w^24)
            sage: (1 + O(w))^5
            1 + O(w^5)
            sage: (1 + w + O(w^3))^25
            1 + w^10 + w^11 + 4*w^12 + O(w^13)
            sage: (3 + 2*w + w^2 + O(w^6))^(15 + O(125))
            2 + 4*w^6 + w^7 + 3*w^8 + 3*w^9 + 4*w^10 + O(w^11)
            sage: (3 + 2*w + w^2 + O(w^6))^(15 + O(25))
            2 + 4*w^6 + w^7 + 3*w^8 + 3*w^9 + O(w^10)
            sage: (3 + w^2 + O(w^6))^(15+O(25))
            2 + w^5 + 4*w^7 + w^9 + 3*w^10 + O(w^11)
            sage: R = ZpCA(2, 10)
            sage: S.<x> = ZZ[]
            sage: f = x^34 + 18*x^5 - 72*x^3 + 2
            sage: W.<w> = R.ext(f)
            sage: (1+w+O(w^2))^8
            1 + w^8 + O(w^16)
            sage: (1+w+O(w^2))^16
            1 + w^16 + O(w^32)
            sage: (1+w+O(w^2))^32
            1 + w^32 + w^50 + w^55 + w^60 + O(w^64)
            sage: (1+w+O(w^2))^64
            1 + w^64 + w^66 + w^71 + w^76 + w^81 + w^84 + w^86 + w^91 + w^94 + w^96 + O(w^98)
            sage: U.<a> = Zq(17^4, 6, print_mode='val-unit'); b = (a^3-a+14)^-6; b
            12003242 + 4839703*a + 2697351*a^2 + 11717046*a^3 + O(17^6)
            sage: b*(a^3-a+14)^6
            1 + O(17^6)
        """
        cdef Integer right
        cdef bint padic_exp
        cdef long exp_prec
        cdef long exp_val
        cdef long ans_relprec, ans_ordp
        cdef long self_ordp = self.valuation_c()
        cdef long self_relprec = self.absprec - self_ordp
        cdef long threshold # e / (p-1)
        cdef long prime_long
        cdef mpz_t tmp, tmp2
        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            threshold = 0
        else:
            threshold = self.prime_pow.e / (mpz_get_si(self.prime_pow.prime.value) - 1)
        cdef Integer base_level
        cdef pAdicZZpXCAElement ans
        cdef long i
        if self._is_inexact_zero():
            # If an integer exponent, return an inexact zero of valuation right * self_ordp.  Otherwise raise an error.
            if isinstance(_right, (int, long)):
                _right = Integer(_right)
            if PY_TYPE_CHECK(_right, Integer):
                mpz_init_set_si(tmp, self_ordp)
                mpz_mul(tmp, tmp, (<Integer>_right).value)
                if mpz_cmp_si(tmp, self.prime_pow.ram_prec_cap) >= 0:
                    ans = self._new_c(self.prime_pow.ram_prec_cap)
                elif mpz_sgn(tmp) <= 0:
                    ans = self._new_c(0)
                else:
                    ans = self._new_c(mpz_get_si(tmp))
                mpz_clear(tmp)
                return ans
            elif PY_TYPE_CHECK(_right, Rational) or (PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime):
                raise ValueError, "Need more precision"
            else:
                raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        if isinstance(_right, (int, long)):
            _right = Integer(_right)
        cdef pAdicZZpXCAElement unit
        if PY_TYPE_CHECK(_right, Integer):
            right = <Integer> _right
            if right < 0 and self_ordp > 0:
                return self.to_fraction_field()**right
            if right == 0:
                # return 1 to maximum precision
                ans = self._new_c(self.prime_pow.ram_prec_cap)
                ZZ_pX_SetCoeff_long(ans.value, 0, 1)
                return ans
            padic_exp = False
            exp_val = _right.valuation(self.prime_pow.prime) ##
        elif PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime:
            if self_ordp != 0:
                raise ValueError, "in order to raise to a p-adic exponent, base must be a unit"
            right = Integer(_right)
            padic_exp = True
            exp_prec = _right.precision_absolute() ##
            exp_val = _right.valuation() ##
            if exp_val < 0:
                raise NotImplementedError, "negative valuation exponents not yet supported"
            # checks to see if the residue of self's unit is in the prime field.
            if self.prime_pow.e == 1:
                unit = self.unit_part()
                for i from 1 <= i <= ZZ_pX_deg(unit.value):
                    if not ZZ_divide_test(ZZ_p_rep(ZZ_pX_coeff(unit.value, i)), self.prime_pow.pow_ZZ_tmp(1)[0]):
                        raise ValueError, "in order to raise to a p-adic exponent, base must reduce to an element of F_p mod the uniformizer"
            # compute the "level"
            teich_part = self.parent().teichmuller(self)
            base_level = (self / teich_part - 1).valuation() ##
        elif PY_TYPE_CHECK(_right, Rational):
            raise NotImplementedError
        else:
            raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        # Now we compute the increased relprec due to the exponent having positive p-adic valuation
        if exp_val > 0:
            mpz_init_set_si(tmp, self_relprec)
            while mpz_cmp_si(tmp, threshold) <= 0 and exp_val > 0:
                mpz_mul(tmp, tmp, self.prime_pow.prime.value)
                exp_val -= 1
            if exp_val > 0:
                mpz_init_set_si(tmp2, self.prime_pow.e)
                mpz_addmul_ui(tmp, tmp2, exp_val)
                mpz_clear(tmp2)
            if mpz_cmp_si(tmp, self.prime_pow.ram_prec_cap) > 0:
                ans_relprec = self.prime_pow.ram_prec_cap
            else:
                ans_relprec = mpz_get_si(tmp)
            mpz_clear(tmp)
        else:
            ans_relprec = self_relprec
        # Now we compute the limit on relprec due to a non-infinite precision on the exponent.
        if padic_exp:
            if exp_prec > 0:
                # I can freely change base_level, so I use it in place of tmp above.
                while mpz_cmp_si(base_level.value, threshold) <= 0 and exp_prec > 0:
                    mpz_mul(base_level.value, base_level.value, self.prime_pow.prime.value)
                    exp_prec -= 1
                if exp_prec > 0:
                    mpz_init_set_si(tmp2, self.prime_pow.e)
                    mpz_addmul_ui(base_level.value, tmp2, exp_prec)
                    mpz_clear(tmp2)
                if mpz_cmp_si(base_level.value, ans_relprec) < 0:
                    ans_relprec = mpz_get_si(base_level.value)
            else:
                return self._new_c(0)
        if self_ordp == 0:
            ans_ordp = 0
        else:
            mpz_init_set(tmp, right.value)
            mpz_mul_si(tmp, tmp, self_ordp)
            if mpz_cmp_si(tmp, self.prime_pow.ram_prec_cap) >= 0:
                return self._new_c(self.prime_pow.ram_prec_cap)
            # we already checked for negative tmp above
            ans_ordp = mpz_get_si(tmp)
            mpz_clear(tmp)
            if ans_ordp >= self.prime_pow.ram_prec_cap:
                return self._new_c(self.prime_pow.ram_prec_cap)
        cdef ntl_ZZ rZZ = PY_NEW(ntl_ZZ)
        mpz_to_ZZ(&rZZ.x, &right.value)
        if ans_ordp + ans_relprec <= self.prime_pow.ram_prec_cap:
            ans = self._new_c(ans_ordp + ans_relprec) # restores context
        else:
            ans = self._new_c(self.prime_pow.ram_prec_cap) # restores context
        cdef ZZ_pX_c self_value
        sig_on()
        if ans.absprec != self.absprec:
            ZZ_pX_conv_modulus(self_value, self.value, self.prime_pow.get_context_capdiv(ans.absprec).x)
            if mpz_sgn(right.value) < 0: # only happens when self.ordp == 0
                if self.prime_pow.e == 1:
                    ZZ_pX_InvMod_newton_unram(ans.value, self_value, self.prime_pow.get_modulus(ans.absprec)[0], self.prime_pow.get_context(ans.absprec).x, self.prime_pow.get_context(1).x)
                else:
                    ZZ_pX_InvMod_newton_ram(ans.value, self_value, self.prime_pow.get_modulus_capdiv(ans.absprec)[0], self.prime_pow.get_context_capdiv(ans.absprec).x)
                ZZ_negate(rZZ.x, rZZ.x)
                ZZ_pX_PowerMod_pre(ans.value, ans.value, rZZ.x, self.prime_pow.get_modulus_capdiv(ans.absprec)[0])
            else:
                ZZ_pX_PowerMod_pre(ans.value, self_value, rZZ.x, self.prime_pow.get_modulus_capdiv(ans.absprec)[0])
        else:
            if mpz_sgn(right.value) < 0: # only happens when self.ordp == 0
                if self.prime_pow.e == 1:
                    ZZ_pX_InvMod_newton_unram(ans.value, self.value, self.prime_pow.get_modulus(ans.absprec)[0], self.prime_pow.get_context(ans.absprec).x, self.prime_pow.get_context(1).x)
                else:
                    ZZ_pX_InvMod_newton_ram(ans.value, self.value, self.prime_pow.get_modulus_capdiv(ans.absprec)[0], self.prime_pow.get_context_capdiv(ans.absprec).x)
                ZZ_negate(rZZ.x, rZZ.x)
                ZZ_pX_PowerMod_pre(ans.value, ans.value, rZZ.x, self.prime_pow.get_modulus_capdiv(ans.absprec)[0])
            else:
                ZZ_pX_PowerMod_pre(ans.value, self.value, rZZ.x, self.prime_pow.get_modulus_capdiv(ans.absprec)[0])
        sig_off()
        return ans

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        Computes the sum of ``self`` and ``right``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)) - 69 # indirect doctest
            1 + O(w^13)
            sage: -69 + (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13))
            1 + O(w^13)
            sage: y = (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13))
            sage: y - 70
            O(w^13)
            sage: y + 0
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
        """
        cdef pAdicZZpXCAElement right = <pAdicZZpXCAElement>_right
        cdef pAdicZZpXCAElement ans
        cdef long tmpL
        cdef ZZ_pX_c tmpP
        if self.absprec == 0 or right.absprec == 0:
            return self._new_c(0)
        elif self.absprec == right.absprec:
            ans = self._new_c(self.absprec)
            ZZ_pX_add(ans.value, self.value, right.value)
        elif self.absprec < right.absprec:
            ans = self._new_c(self.absprec)
            ZZ_pX_conv_modulus(tmpP, right.value, self.prime_pow.get_context_capdiv(ans.absprec).x)
            ZZ_pX_add(ans.value, self.value, tmpP)
        else:
            ans = self._new_c(right.absprec)
            ZZ_pX_conv_modulus(tmpP, self.value, self.prime_pow.get_context_capdiv(ans.absprec).x)
            ZZ_pX_add(ans.value, tmpP, right.value)
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        Returns the difference of ``self`` and ``right``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(329)
            sage: b = W(111)
            sage: a - b #indirect doctest
            3 + 3*w^5 + w^7 + 2*w^9 + 3*w^10 + 4*w^11 + 2*w^13 + 2*w^14 + w^15 + 4*w^16 + 2*w^18 + 3*w^19 + 2*w^20 + 3*w^21 + w^22 + w^24 + O(w^25)
            sage: W(218)
            3 + 3*w^5 + w^7 + 2*w^9 + 3*w^10 + 4*w^11 + 2*w^13 + 2*w^14 + w^15 + 4*w^16 + 2*w^18 + 3*w^19 + 2*w^20 + 3*w^21 + w^22 + w^24 + O(w^25)
            sage: a - O(w^14)
            4 + 3*w^10 + 2*w^12 + O(w^14)
            sage: a - 0
            4 + 3*w^10 + 2*w^12 + w^14 + 2*w^15 + w^16 + 3*w^17 + 3*w^18 + w^19 + 2*w^21 + 4*w^22 + w^23 + 4*w^24 + O(w^25)
            sage: O(w^14) - a
            1 + 4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + O(w^14)
        """
        cdef pAdicZZpXCAElement right = <pAdicZZpXCAElement>_right
        cdef pAdicZZpXCAElement ans
        cdef long tmpL
        cdef ZZ_pX_c tmpP
        if self.absprec == 0 or right.absprec == 0:
            return self._new_c(0)
        elif self.absprec == right.absprec:
            ans = self._new_c(self.absprec)
            ZZ_pX_sub(ans.value, self.value, right.value)
        elif self.absprec < right.absprec:
            ans = self._new_c(self.absprec)
            ZZ_pX_conv_modulus(tmpP, right.value, self.prime_pow.get_context_capdiv(ans.absprec).x)
            ZZ_pX_sub(ans.value, self.value, tmpP)
        else:
            ans = self._new_c(right.absprec)
            ZZ_pX_conv_modulus(tmpP, self.value, self.prime_pow.get_context_capdiv(ans.absprec).x)
            ZZ_pX_sub(ans.value, tmpP, right.value)
        return ans

    cpdef RingElement _mul_(self, RingElement _right):
        """
        Returns the product of ``self`` and ``right``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(329)
            sage: b = W(111)
            sage: a*b #indirect doctest
            4 + 3*w^5 + w^7 + 2*w^9 + 4*w^11 + 3*w^12 + 2*w^13 + w^14 + 2*w^15 + 3*w^16 + 4*w^17 + 4*w^18 + 2*w^19 + 2*w^21 + 4*w^22 + 2*w^23 + w^24 + O(w^25)
            sage: a * 0
            O(w^25)
            sage: a * O(w^14)
            O(w^14)
        """
        cdef pAdicZZpXCAElement right = <pAdicZZpXCAElement>_right
        cdef pAdicZZpXCAElement ans
        cdef ZZ_pX_c self_adapted, right_adapted
        cdef long self_ordp = self.valuation_c()
        cdef long right_ordp = right.valuation_c()
        cdef long ans_ordp = self_ordp + right_ordp
        if ans_ordp >= self.prime_pow.ram_prec_cap:
            return self._new_c(self.prime_pow.ram_prec_cap)
        if self._is_inexact_zero() or right._is_inexact_zero():
            return self._new_c(ans_ordp)
        cdef long self_relprec = self.absprec - self_ordp
        cdef long right_relprec = right.absprec - right_ordp
        cdef long ans_absprec
        if self_relprec <= right_relprec:
            ans_absprec = ans_ordp + self_relprec
        else:
            ans_absprec = ans_ordp + right_relprec
        if ans_absprec > self.prime_pow.ram_prec_cap:
            ans_absprec = self.prime_pow.ram_prec_cap
        ans = self._new_c(ans_absprec) # restores the context
        if self.absprec == ans_absprec and right.absprec == ans_absprec:
            ZZ_pX_MulMod_pre(ans.value, self.value, right.value, self.prime_pow.get_modulus_capdiv(ans_absprec)[0])
        elif self.absprec == ans_absprec:
            ZZ_pX_conv_modulus(right_adapted, right.value, self.prime_pow.get_context_capdiv(ans_absprec).x)
            ZZ_pX_MulMod_pre(ans.value, self.value, right_adapted, self.prime_pow.get_modulus_capdiv(ans_absprec)[0])
        elif right.absprec == ans_absprec:
            ZZ_pX_conv_modulus(self_adapted, self.value, self.prime_pow.get_context_capdiv(ans_absprec).x)
            ZZ_pX_MulMod_pre(ans.value, self_adapted, right.value, self.prime_pow.get_modulus_capdiv(ans_absprec)[0])
        else:
            ZZ_pX_conv_modulus(self_adapted, self.value, self.prime_pow.get_context_capdiv(ans_absprec).x)
            ZZ_pX_conv_modulus(right_adapted, right.value, self.prime_pow.get_context_capdiv(ans_absprec).x)
            ZZ_pX_MulMod_pre(ans.value, self_adapted, right_adapted, self.prime_pow.get_modulus_capdiv(ans_absprec)[0])
        return ans

    cpdef RingElement _div_(self, RingElement right):
        """
        Returns the quotient of ``self`` by ``right``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(14) / W(125) #indirect doctest
            4*w^-15 + w^-13 + 3*w^-11 + 2*w^-10 + 3*w^-9 + 4*w^-8 + 4*w^-7 + 3*w^-6 + O(w^-5)
            sage: 1 / w
            w^-1 + O(w^23)
            sage: W.<w> = R.ext(x^20 - 165*x + 5)
            sage: a = (1 + w)^25 - 1
            sage: b = (1 + w)^5 - 1
            sage: c = (1 + w)^20 + (1 + w)^15 + (1 + w)^10 + (1 + w)^5 + 1
            sage: d = a / b; d == c
            True
            sage: d.precision_absolute()
            95
            sage: c.precision_absolute()
            100
            sage: 1 / a == ~a
            True
        """
        return self.to_fraction_field() * (~right)

    def _integer_(self, Z=None):
        """
        Returns an integer congruent to this element modulo
        `\pi`^``self.absolute_precision()``, if possible.

        EXAMPLES::

            sage: ZZ(ZqCA(125,names='a')(-1)) #indirect doctest
            95367431640624
            sage: R = ZpCA(5); S.<x> = ZZ[]; f = x^5 + 25*x^3 - 5; W.<w> = R.ext(f)
            sage: ZZ(W(-1))
            95367431640624
            sage: ZZ(W(0))
            0
            sage: ZZ(W(0,7))
            0
            sage: ZZ(w)
            Traceback (most recent call last):
            ...
            ValueError: This element not well approximated by an integer.
            sage: ZZ(W(5))
            5
        """
        cdef Integer ans
        cdef ZZ_c tmp_z
        if ZZ_pX_deg(self.value) > 0:
            raise ValueError, "This element not well approximated by an integer."
        ans = PY_NEW(Integer)
        tmp_z = ZZ_p_rep(ZZ_pX_ConstTerm(self.value))
        ZZ_to_mpz(&ans.value, &tmp_z)
        return ans

    def __copy__(self):
        """
        Returns a copy of ``self``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: b = W(45, 17); b
            4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16 + O(w^17)
            sage: c = copy(b); c
            4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16 + O(w^17)
            sage: c is b
            False
        """
        cdef pAdicZZpXCAElement ans = self._new_c(self.absprec) # restores context
        ans.value = self.value
        return ans

    def is_zero(self, absprec = None):
        """
        Returns whether the valuation of ``self`` is at least ``absprec``.  If
        ``absprec`` is ``None``, returns if ``self`` is indistinguishable from
        zero.

        If ``self`` is an inexact zero of valuation less than ``absprec``,
        raises a PrecisionError.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: O(w^189).is_zero()
            True
            sage: W(0).is_zero()
            True
            sage: a = W(675)
            sage: a.is_zero()
            False
            sage: a.is_zero(7)
            True
            sage: a.is_zero(21)
            False
        """
        cdef bint ans
        cdef long aprec
        if absprec is None:
            ans = ZZ_pX_IsZero(self.value)
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0:
                    ans = True
                elif ZZ_pX_IsZero(self.value):
                    raise PrecisionError, "Not enough precision to determine if element is zero"
                else:
                    ans = False
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
                if ZZ_pX_IsZero(self.value) and aprec > self.absprec:
                    raise PrecisionError, "Not enough precision to determine if element is zero"
                else:
                    ans = (self.valuation_c() >= aprec)
        return ans

    cpdef ntl_ZZ_pX _ntl_rep(self):
        """
        Returns an ``ntl_ZZ_pX`` that holds the value of ``self``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566); b = W(209)
            sage: c = a + b; c._ntl_rep() # indirect doctest
            [775]
        """
        if self.absprec == 0:
            raise ValueError, "self has 0 absolute precision"
        self.prime_pow.restore_context_capdiv(self.absprec)
        cdef ntl_ZZ_pX ans = PY_NEW(ntl_ZZ_pX)
        ans.c = self.prime_pow.get_context_capdiv(self.absprec)
        ans.x = self.value
        return ans

    cpdef _ntl_rep_abs(self):
        """
        Returns a pair ``(f, 0)`` where ``f = self._ntl_rep()``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566); b = W(209)
            sage: c = a + b; c._ntl_rep_abs()
            ([775], 0)
            sage: c
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24 + O(w^25)
            sage: c._ntl_rep_abs()
            ([775], 0)
        """
        return self._ntl_rep(), Integer(0)

    cdef ZZ_p_c _const_term(self):
        """
        Returns the constant term of ``self.value``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566)
            sage: a._const_term_test() #indirect doctest
            566
        """
        return ZZ_pX_ConstTerm(self.value)

    def is_equal_to(self, right, absprec = None):
        """
        Returns whether ``self`` is equal to ``right`` modulo
        ``self.uniformizer()^absprec``.

        If ``absprec`` is ``None``, returns if ``self`` is equal to ``right`` modulo
        the lower of their two precisions.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(47); b = W(47 + 25)
            sage: a.is_equal_to(b)
            False
            sage: a.is_equal_to(b, 7)
            True
        """
        # Should be sped up later
        return (self - right).is_zero(absprec)

#     def lift(self):
#         """
#         Returns an element of a number field defined by the same polynomial as self's parent that is congruent to self modulo an appropriate ideal.

#         Not currently implemented.
#         """

#         raise NotImplementedError

    cpdef pAdicZZpXCAElement lift_to_precision(self, absprec=None):
        """
        Returns a ``pAdicZZpXCAElement`` congruent to ``self`` but with
        absolute precision at least ``absprec``.

        INPUT:

        - ``absprec`` -- (default ``None``) the absolute precision of
          the result.  If ``None``, lifts to the maximum precision
          allowed.

        .. NOTE::

            If setting ``absprec`` that high would violate the
            precision cap, raises a precision error.

            Note that the new digits will not necessarily be zero.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(345, 17); a
            4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15 + O(w^17)
            sage: b = a.lift_to_precision(19); b # indirect doctest
            4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15 + w^17 + 2*w^18 + O(w^19)
            sage: c = a.lift_to_precision(24); c
            4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15 + w^17 + 2*w^18 + 4*w^19 + 4*w^20 + 2*w^21 + 4*w^23 + O(w^24)
            sage: a._ntl_rep()
            [345]
            sage: b._ntl_rep()
            [345]
            sage: c._ntl_rep()
            [345]
            sage: a.lift_to_precision().precision_absolute() == W.precision_cap()
            True
        """
        cdef pAdicZZpXCAElement ans
        cdef long aprec, rprec
        if absprec is not None and not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if absprec is None:
            aprec = self.prime_pow.ram_prec_cap
        elif mpz_fits_slong_p((<Integer>absprec).value) == 0:
            if mpz_sgn((<Integer>absprec).value) < 0:
                return self
            else:
                raise PrecisionError("Precision higher than allowed by the precision cap.")
        else:
            aprec = mpz_get_si((<Integer>absprec).value)
            if aprec > self.prime_pow.ram_prec_cap:
                raise PrecisionError("Precision higher than allowed by the precision cap.")
        if aprec <= self.absprec:
            return self
        ans = self._new_c(aprec) # restores context
        ZZ_pX_conv_modulus(ans.value, self.value, self.prime_pow.get_context_capdiv(aprec).x)
        return ans

    def list(self, lift_mode = 'simple'):
        """
        Returns a list giving a series representation of ``self``.

        - If ``lift_mode == 'simple'`` or ``'smallest'``, the returned
          list will consist of integers (in the eisenstein case) or a
          list of lists of integers (in the unramified case).
          ``self`` can be reconstructed as a sum of elements of the
          list times powers of the uniformiser (in the eisenstein
          case), or as a sum of powers of `p` times polynomials in the
          generator (in the unramified case).

          + If ``lift_mode == 'simple'``, all integers will be in the
            interval `[0,p-1]`

          + If ``lift_mod == 'smallest'`` they will be in the
            interval `[(1-p)/2, p/2]`.

        - If ``lift_mode == 'teichmuller'``, returns a list of
          ``pAdicZZpXCAElements``, all of which are Teichmuller
          representatives and such that ``self`` is the sum of that
          list times powers of the uniformizer.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W(775, 19); y
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: (y>>9).list()
            [0, 1, 0, 4, 0, 2, 1, 2, 4, 1]
            sage: (y>>9).list('smallest')
            [0, 1, 0, -1, 0, 2, 1, 2, 0, 1]
            sage: w^10 - w^12 + 2*w^14 + w^15 + 2*w^16 + w^18 + O(w^19)
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: y = 75 + 45*a + 1200*a^2; y
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^5)
            sage: y.list()
            [[], [0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1]]
            sage: y.list('smallest')
            [[], [0, -1], [-2, 2, -2], [1], [0, 0, 2]]
            sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^5)
            sage: W(0).list()
            [0]
            sage: A(0,4).list()
            [[]]
        """
        if lift_mode == 'simple':
            ulist = self.ext_p_list(1)
        elif lift_mode == 'smallest':
            ulist = self.ext_p_list(0)
        elif lift_mode == 'teichmuller':
            ulist = self.teichmuller_list()
        else:
            raise ValueError, "lift mode must be one of 'simple', 'smallest' or 'teichmuller'"
        ordp = self.valuation()
        if self.is_zero():
            ordp = 1
        if lift_mode == 'teichmuller':
            zero = self._new_c(0)
            return [zero]*ordp + ulist
        elif self.prime_pow.e == 1:
            return [[]] * ordp + ulist
        else:
            return [Integer(0)] * ordp + ulist

    def matrix_mod_pn(self):
        """
        Returns the matrix of right multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for this
        extension field.  Thus the *rows* of this matrix give the
        images of each of the `x^i`.  The entries of the matrices are
        ``IntegerMod`` elements, defined modulo ``p^(self.absprec() / e)``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = (3+w)^7
            sage: a.matrix_mod_pn()
            [2757  333 1068  725 2510]
            [  50 1507  483  318  725]
            [ 500   50 3007 2358  318]
            [1590 1375 1695 1032 2358]
            [2415  590 2370 2970 1032]
        """
        from sage.matrix.all import matrix
        # this may be the wrong precision when ram_prec_cap is not divisible by e.
        R = IntegerModRing(self.prime_pow.pow_Integer(self.prime_pow.capdiv(self.absprec)))
        n = self.prime_pow.deg
        L = []
        cdef ntl_ZZ_pX cur = <ntl_ZZ_pX>self._ntl_rep()
        cur.c.restore_c()
        cdef ZZ_pX_Modulus_c* m = self.prime_pow.get_modulus_capdiv(self.absprec)
        cdef ZZ_pX_c x
        ZZ_pX_SetX(x)
        cdef Py_ssize_t i, j
        zero = int(0)
        for i from 0 <= i < n:
            curlist = cur.list()
            L.extend(curlist + [zero]*(n - len(curlist)))
            ZZ_pX_MulMod_pre(cur.x, cur.x, x, m[0])
        return matrix(R, n, n,  L)

#     def matrix(self, base = None):
#         """
#         If base is None, return the matrix of right multiplication by
#         the element on the power basis `1, x, x^2, \ldots, x^{d-1}`
#         for this extension field.  Thus the \emph{rows} of this matrix
#         give the images of each of the `x^i`.

#         If base is not None, then base must be either a field that
#         embeds in the parent of self or a morphism to the parent of
#         self, in which case this function returns the matrix of
#         multiplication by self on the power basis, where we view the
#         parent field as a field over base.

#         INPUT::

#             - base -- field or morphism
#         """
#         raise NotImplementedError

#     def multiplicative_order(self, prec=None):
#         """
#         Returns the multiplicative order of self, ie the smallest
#         positive n so that there is an exact p-adic element congruent
#         to self modulo self's precision that is an nth root of unity.

#         Note: unlike the case for Qp and Zp, it is possible to have
#         non-teichmuller elements with finite orders.  This can happen
#         only if (p-1) divides the ramification index (see the
#         documentation on __pow__).

#         INPUT::

#             - self -- a p-adic element
#             - prec -- an integer

#         OUTPUT::

#             - integer -- the multiplicative order of self
#         """
#         raise NotImplementedError

    def teichmuller_list(self):
        r"""
        Returns a list [`a_0`, `a_1`,..., `a_n`] such that

        - `a_i^q = a_i`
        - ``self.unit_part()`` = `\sum_{i = 0}^n a_i \pi^i`, where `\pi` is a
          uniformizer of self.parent()
        - if `a_i \ne 0`, the absolute precision of `a_i` is
          ``self.precision_relative() - i``

        EXAMPLES::

            sage: R.<a> = Zq(5^4,4)
            sage: L = a.teichmuller_list(); L
            [a + (2*a^3 + 2*a^2 + 3*a + 4)*5 + (4*a^3 + 3*a^2 + 3*a + 2)*5^2 + (4*a^2 + 2*a + 2)*5^3 + O(5^4), (3*a^3 + 3*a^2 + 2*a + 1) + (a^3 + 4*a^2 + 1)*5 + (a^2 + 4*a + 4)*5^2 + O(5^3), (4*a^3 + 2*a^2 + a + 1) + (2*a^3 + 2*a^2 + 2*a + 4)*5 + O(5^2), (a^3 + a^2 + a + 4) + O(5)]
            sage: sum([5^i*L[i] for i in range(4)])
            a + O(5^4)
            sage: all([L[i]^625 == L[i] for i in range(4)])
            True

            sage: S.<x> = ZZ[]
            sage: f = x^3 - 98*x + 7
            sage: W.<w> = ZpCA(7,3).ext(f)
            sage: b = (1+w)^5; L = b.teichmuller_list(); L
            [1 + O(w^9), 5 + 5*w^3 + w^6 + 4*w^7 + O(w^8), 3 + 3*w^3 + O(w^7), 3 + 3*w^3 + O(w^6), O(w^5), 4 + 5*w^3 + O(w^4), 3 + O(w^3), 6 + O(w^2), 6 + O(w)]
            sage: sum([w^i*L[i] for i in range(9)]) == b
            True
            sage: all([L[i]^(7^3) == L[i] for i in range(9)])
            True

            sage: L = W(3).teichmuller_list(); L
            [3 + 3*w^3 + w^7 + O(w^9), O(w^8), O(w^7), 4 + 5*w^3 + O(w^6), O(w^5), O(w^4), 3 + O(w^3), 6 + O(w^2)]
            sage: sum([w^i*L[i] for i in range(len(L))])
            3 + O(w^9)
        """
        L = []
        cdef long rp = self.absprec - self.valuation_c()
        if rp == 0:
            return L
        cdef pAdicZZpXCAElement u = self.unit_part()
        if u is self: u = self.__copy__()
        cdef pAdicZZpXCAElement v
        while not ZZ_pX_IsZero(u.value):
            v = self._new_c(rp)
            self.prime_pow.teichmuller_set_c(&v.value, &u.value, rp)
            L.append(v)
            if rp == 1: break
            ZZ_pX_sub(u.value, u.value, v.value)
            rp -= 1
            if self.prime_pow.e == 1:
                ZZ_pX_right_pshift(u.value, u.value, self.prime_pow.pow_ZZ_tmp(1)[0], self.prime_pow.get_context(rp).x)
            else:
                self.prime_pow.eis_shift_capdiv(&u.value, &u.value, 1, rp)
        return L


    def _teichmuller_set(self):
        """
        Sets self to the teichmuller representative congruent to self
        modulo `\pi`, with the same relative precision as ``self``.

        This function should not be used externally: elements are
        supposed to be immutable.

        EXAMPLES::

            sage: R = ZpCA(11,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 33*x^3 - 121*x^2 - 77
            sage: W.<w> = R.ext(f)
            sage: y = W.teichmuller(3, 19); y #indirect doctest
            3 + 9*w^10 + 3*w^13 + 3*w^15 + 9*w^16 + 3*w^17 + w^18 + O(w^19)

            sage: y^11 == y
            True
            sage: g = x^3 + 9*x^2 + 7
            sage: A.<a> = R.ext(g)
            sage: b = A.teichmuller(1 + 2*a - a^2); b
            (10*a^2 + 2*a + 1) + (4*a^2 + 7)*11 + (5*a^2 + a + 3)*11^2 + (a^2 + 9*a + 6)*11^3 + (7*a^2 + 2*a + 3)*11^4 + O(11^5)
            sage: b^1331 == b
            True
        """
        if self.absprec == 0:
            raise ValueError, "not enough precision known"
        elif self.valuation_c() > 0:
            return self._new_c(self.prime_pow.ram_prec_cap)
        else:
            self.prime_pow.teichmuller_set_c(&self.value, &self.value, self.absprec)

#     def padded_list(self, n, lift_mode = 'simple'):
#         """
#         Returns a list of coefficients of pi starting with `pi^0` up to
#         `pi^n` exclusive (padded with zeros if needed)

#         """
#         raise NotImplementedError

    def precision_absolute(self):
        """
        Returns the absolute precision of ``self``, ie the power of the
        uniformizer modulo which this element is defined.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75, 19); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
            sage: a.valuation()
            10
            sage: a.precision_absolute()
            19
            sage: a.precision_relative()
            9
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.absprec)
        return ans

    def precision_relative(self):
        """
        Returns the relative precision of ``self``, ie the power of
        the uniformizer modulo which the unit part of ``self`` is
        defined.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75, 19); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
            sage: a.valuation()
            10
            sage: a.precision_absolute()
            19
            sage: a.precision_relative()
            9
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.absprec - self.valuation_c())
        return ans

#    def residue(self, n):
#        """
#        Reduces this element modulo pi^n.
#        """
#        raise NotImplementedError

    cdef long valuation_c(self):
        """
        Returns the valuation of ``self``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75, 19); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
            sage: a.valuation() # indirect doctest
            10
            sage: a.precision_absolute()
            19
            sage: a.precision_relative()
            9
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        if ZZ_pX_IsZero(self.value):
            return self.absprec
        cdef long minval, mini, val
        ZZ_pX_min_val_coeff(minval, mini, self.value, self.prime_pow.pow_ZZ_tmp(1)[0])
        if self.prime_pow.e == 1:
            if minval <= self.absprec:
                return minval
            else:
                return self.absprec
        else:
            val = minval * self.prime_pow.e + mini
            if val <= self.absprec:
                return val
            else:
                return self.absprec

    cpdef pAdicZZpXCAElement unit_part(self):
        """
        Returns the unit part of ``self``, ie ``self / uniformizer^(self.valuation())``

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75, 19); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
            sage: a.valuation()
            10
            sage: a.precision_absolute()
            19
            sage: a.precision_relative()
            9
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        return self._rshift_c(self.valuation_c())

    cdef ext_p_list(self, bint pos):
        """
        Returns a list of integers (in the eisenstein case) or a list
        of lists of integers (in the unramified case).  ``self`` can
        be reconstructed as a sum of elements of the list times powers
        of the uniformizer (in the eisenstein case), or as a sum of
        powers of `p` times polynomials in the generator (in the
        unramified case).

        If ``pos`` is ``True``, all integers will be in the interval
        `[0,p-1]`, otherwise they will be in the range
        `[(1-p)/2,p/2]`.

        Note that zeros are truncated from the returned list, so you
        must use the ``valuation()`` function to completely recover
        ``self``.

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W(775, 19); y
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: y._ext_p_list(True)
            [1, 0, 4, 0, 2, 1, 2, 4, 1]
            sage: y._ext_p_list(False)
            [1, 0, -1, 0, 2, 1, 2, 0, 1]
            sage: w^10 - w^12 + 2*w^14 + w^15 + 2*w^16 + w^18 + O(w^19)
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: y = 75 + 45*a + 1200*a^2; y
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^5)
            sage: y._ext_p_list(True)
            [[0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1]]
            sage: y._ext_p_list(False)
            [[0, -1], [-2, 2, -2], [1], [0, 0, 2]]
            sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^5)
        """
        return self.ext_p_list_precs(pos, self.absprec)

def make_ZZpXCAElement(parent, value, absprec, version):
    """
    For pickling.  Makes a ``pAdicZZpXCAElement`` with given ``parent``, ``value``, ``absprec``.

    EXAMPLES::

        sage: from sage.rings.padics.padic_ZZ_pX_CA_element import make_ZZpXCAElement
        sage: R = ZpCA(5,5)
        sage: S.<x> = ZZ[]
        sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
        sage: W.<w> = R.ext(f)
        sage: make_ZZpXCAElement(W, ntl.ZZ_pX([3,2,4],5^3),13,0)
        3 + 2*w + 4*w^2 + O(w^13)
    """
    cdef pAdicZZpXCAElement ans
    cdef ZZ_pX_c poly
    if version == 0:
        ans = pAdicZZpXCAElement(parent, [], empty = True)
        if mpz_sgn((<Integer>absprec).value) == 0:
            ans._set_inexact_zero(0)
        else:
            ans.prime_pow.restore_context_capdiv(mpz_get_si((<Integer>absprec).value))
            poly = (<ntl_ZZ_pX>value).x
            ans._set(&poly, mpz_get_si((<Integer>absprec).value))
        return ans
    else:
        raise ValueError, "unknown unpickling version"
