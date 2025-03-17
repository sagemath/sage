# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
# sage.doctest: needs sage.rings.padics
r"""
`p`-adic ``ZZ_pX`` CR Element

This file implements elements of Eisenstein and unramified extensions
of `\ZZ_p` and `\QQ_p` with capped relative precision.

For the parent class see :mod:`sage.rings.padics.padic_extension_leaves`.

The underlying implementation is through NTL's ``ZZ_pX`` class.  Each
element contains the following data:

- ``ordp`` -- ``long``; a power of the uniformizer to scale the unit
  by.  For unramified extensions this uniformizer is `p`, for Eisenstein
  extensions it is not.  A value equal to the maximum value of a ``long``
  indicates that the element is an exact zero.

- ``relprec`` -- ``long``; a signed integer giving the precision to
  which this element is defined.  For nonzero ``relprec``, the
  absolute value gives the power of the uniformizer modulo which the
  unit is defined.  A positive value indicates that the element is
  normalized (ie ``unit`` is actually a unit: in the case of
  Eisenstein extensions the constant term is not divisible by `p`, in
  the case of unramified extensions that there is at least one
  coefficient that is not divisible by `p`).  A negative value
  indicates that the element may or may not be normalized.  A zero
  value indicates that the element is zero to some precision.  If so,
  ``ordp`` gives the absolute precision of the element.  If ``ordp``
  is greater than ``maxordp``, then the element is an exact zero.

- ``unit`` -- ``ZZ_pX_c``; an ntl ``ZZ_pX`` storing the unit part
  The variable `x` is the uniformizer in the case of Eisenstein
  extensions. If the element is not normalized, the ``unit`` may or
  may not actually be a unit.  This ``ZZ_pX`` is created with global
  ntl modulus determined by the absolute value of ``relprec``.  If
  ``relprec`` is 0, ``unit`` **is not initialized**, or destructed if
  normalized and found to be zero.  Otherwise, let `r` be ``relprec`` and
  `e` be the ramification index over `\QQ_p` or `\ZZ_p`.
  Then the modulus of unit is given by `p^{ceil(r/e)}`.  Note that all
  kinds of problems arise if you try to mix moduli.
  ``ZZ_pX_conv_modulus`` gives a semi-safe way to convert between
  different moduli without having to pass through ``ZZX``.

- ``prime_pow`` (some subclass of :class:`PowComputer_ZZ_pX`) -- a class,
  identical among all elements with the same parent, holding common
  data.

  + ``prime_pow.deg`` -- the degree of the extension

  + ``prime_pow.e`` -- the ramification index

  + ``prime_pow.f`` -- the inertia degree

  + ``prime_pow.prec_cap`` -- the unramified precision cap.  For
    Eisenstein extensions this is the smallest power of `p` that is
    zero.

  + ``prime_pow.ram_prec_cap`` -- the ramified precision cap.  For
    Eisenstein extensions this will be the smallest power of `x` that
    is indistinguishable from zero.

  + ``prime_pow.pow_ZZ_tmp``, prime_pow.pow_mpz_t_tmp``,
    ``prime_pow.pow_Integer`` -- functions for accessing powers of
    `p`.  The first two return pointers.  See
    :mod:`sage.rings.padics.pow_computer_ext` for examples and important
    warnings.

  + ``prime_pow.get_context``, ``prime_pow.get_context_capdiv``,
    ``prime_pow.get_top_context`` -- obtain an
    ``ntl_ZZ_pContext_class`` corresponding to `p^n`.  The capdiv
    version divides by ``prime_pow.e`` as appropriate.
    ``top_context`` corresponds to `p^{\texttt{prec\_cap}}`.

  + ``prime_pow.restore_context``,
    ``prime_pow.restore_context_capdiv``,
    ``prime_pow.restore_top_context`` -- restores the given context.

  + ``prime_pow.get_modulus``, ``get_modulus_capdiv``,
    ``get_top_modulus`` -- returns a ``ZZ_pX_Modulus_c*`` pointing to
    a polynomial modulus defined modulo `p^n` (appropriately divided
    by ``prime_pow.e`` in the capdiv case).

EXAMPLES:

An Eisenstein extension::

    sage: R = Zp(5,5)
    sage: S.<x> = R[]
    sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
    sage: W.<w> = R.ext(f); W
    5-adic Eisenstein Extension Ring in w defined by x^5 + 75*x^3 - 15*x^2 + 125*x - 5
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
    w^-12 + w + O(w^13)
    sage: (1/w).parent()
    5-adic Eisenstein Extension Field in w defined by x^5 + 75*x^3 - 15*x^2 + 125*x - 5

Unramified extensions::

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
    sage: FFp = R.residue_field()
    sage: R(FFp(3))
    3 + O(5)
    sage: QQq.<zz> = Qq(25,4)
    sage: QQq(FFp(3))
    3 + O(5)
    sage: FFq = QQq.residue_field(); QQq(FFq(3))
    3 + O(5)
    sage: zz0 = FFq.gen(); QQq(zz0^2)
    (zz + 3) + O(5)

Different printing modes::

    sage: R = Zp(5, print_mode='digits'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; repr(z)
    '...4110403113210310442221311242000111011201102002023303214332011214403232013144001400444441030421100001'
    sage: R = Zp(5, print_mode='bars'); S.<x> = R[]; g = x^3 + 3*x + 3; A.<a> = R.ext(g)
    sage: z = (1+a)^5; repr(z)
    '...[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 3, 4]|[1, 3, 3]|[0, 4, 2]'
    sage: R = Zp(5, print_mode='terse'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; z
    6 + 95367431640505*w + 25*w^2 + 95367431640560*w^3 + 5*w^4 + O(w^100)
    sage: R = Zp(5, print_mode='val-unit'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: y = (1+w)^5 - 1; y
    w^5 * (2090041 + 19073486126901*w + 1258902*w^2 + 674*w^3 + 16785*w^4) + O(w^100)

You can get at the underlying ntl unit::

    sage: z._ntl_rep()
    [6 95367431640505 25 95367431640560 5]
    sage: y._ntl_rep()
    [2090041 19073486126901 1258902 674 16785]
    sage: y._ntl_rep_abs()
    ([5 95367431640505 25 95367431640560 5], 0)

.. NOTE::

    If you get an error ``internal error: can't grow this _ntl_gbigint,``
    it indicates that moduli are being mixed inappropriately somewhere.

    For example, when calling a function with a ``ZZ_pX_c`` as an
    argument, it copies.  If the modulus is not
    set to the modulus of the ``ZZ_pX_c``, you can get errors.

AUTHORS:

- David Roe (2008-01-01): initial version

- Robert Harron (2011-09): fixes/enhancements

- Julian Rueth (2014-05-09): enable caching through ``_cache_key``
"""
# ****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from sage.ext.stdsage cimport PY_NEW
include "sage/libs/ntl/decl.pxi"

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from cypari2.gen cimport Gen as pari_gen
from sage.interfaces.abc import GpElement
from sage.rings.finite_rings.integer_mod import IntegerMod_abstract
from sage.rings.padics.padic_ext_element cimport pAdicExtElement
from sage.rings.padics.precision_error import PrecisionError

from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX_small_Eis
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX_big_Eis
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

cdef object infinity
from sage.rings.infinity import infinity

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) -1
cdef long minusmaxordp = -maxordp

cdef inline int check_ordp(long a) except -1:
    if a > maxordp or a < minusmaxordp:
        raise ValueError("valuation overflow")


cdef class pAdicZZpXCRElement(pAdicZZpXElement):
    def __init__(self, parent, x, absprec=infinity, relprec=infinity, empty=False):
        r"""
        Create an element of a capped relative precision, unramified
        or Eisenstein extension of `\ZZ_p` or `\QQ_p`.

        INPUT:

        - ``parent`` -- either an ``EisensteinRingCappedRelative`` or
          ``UnramifiedRingCappedRelative``

        - ``x`` -- integer; rational, `p`-adic element, polynomial,
          list, integer_mod, pari int/frac/poly_t/pol_mod, an
          ``ntl_ZZ_pX``, an ``ntl_ZZ``, an ``ntl_ZZ_p``, an
          ``ntl_ZZX``, or something convertible into parent.residue_field()

        - ``absprec`` -- an upper bound on the absolute precision of the
          element created

        - ``relprec`` -- an upper bound on the relative precision of
          the element created

        - ``empty`` -- whether to return after initializing to zero
          (without setting the valuation)

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1+w)^5; z  # indirect doctest
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: W(pari('3 + O(5^3)'))
            3 + O(w^15)
            sage: W(R(3,3))
            3 + O(w^15)
            sage: W.<w> = R.ext(x^625 + 915*x^17 - 95)
            sage: W(3)
            3 + O(w^3125)
            sage: W(w, 14)
            w + O(w^14)

        TESTS:

        Check that :issue:`3865` is fixed::

            sage: W(gp('3 + O(5^10)'))
            3 + O(w^3125)


        Check that :issue:`13612` has been fixed::

            sage: R = Zp(3)
            sage: S.<a> = R[]
            sage: W.<a> = R.extension(a^2 + 1)
            sage: W(W.residue_field().zero())
            O(3)

            sage: K = Qp(3)
            sage: S.<a> = K[]
            sage: L.<a> = K.extension(a^2 + 1)
            sage: L(L.residue_field().zero())
            O(3)
        """
        pAdicZZpXElement.__init__(self, parent)
        self.relprec = 0
        if empty:
            return
        cdef long aprec, rprec, ctx_prec, ltmp
        if relprec is not infinity and not isinstance(relprec, Integer):
            relprec = Integer(relprec)
        if (relprec is infinity) or (relprec > parent.precision_cap()):
            rprec = self.prime_pow.ram_prec_cap
        else:
            rprec = mpz_get_si((<Integer>relprec).value)
        if rprec < 0:
            rprec = 0
        if absprec is not infinity:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                absprec = infinity
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
        cdef mpz_t tmp
        cdef ZZ_c tmp_z
        cdef Py_ssize_t i
        cdef Integer tmp_Int
        cdef Rational xlift
        if isinstance(x, pAdicGenericElement):
            if self.prime_pow.in_field == 0 and x.valuation() < 0:
                raise ValueError("element has negative valuation")
            if x._is_base_elt(self.prime_pow.prime):
                xlift = Rational(x.lift())
                if mpq_sgn(xlift.value) == 0:
                    if (<pAdicGenericElement>x)._is_exact_zero():
                        if absprec is infinity:
                            self._set_exact_zero()
                        else:
                            self._set_inexact_zero(aprec)
                        return
                ltmp = mpz_get_si((<Integer>x.precision_absolute()).value) * self.prime_pow.e
                if absprec is infinity or ltmp < aprec:
                    aprec = ltmp
                self._set_from_mpq_both(xlift.value, aprec, rprec)
                return
        if isinstance(x, GpElement):
            x = x.__pari__()
        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                if x.variable() != self.prime_pow.prime:
                    raise TypeError("cannot coerce a pari p-adic with the wrong prime")
                ltmp = x.padicprec(self.prime_pow.prime) * self.prime_pow.e
                if absprec is infinity or ltmp < aprec:
                    aprec = ltmp
                    absprec = 0 # absprec just has to be non-infinite: everything else uses aprec
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
                    L.append(Integer(x.polcoef(i)))
                x = L
            else:
                raise TypeError("unsupported coercion from pari: only p-adics, integers, rationals, polynomials and pol_mods allowed")
        elif isinstance(x, IntegerMod_abstract):
            mpz_init(tmp)
            ctx_prec = mpz_remove(tmp, (<Integer>x.modulus()).value, self.prime_pow.prime.value)
            if mpz_cmp_ui(tmp, 1) == 0:
                mpz_clear(tmp)
                x = x.lift()
                if absprec is infinity or ctx_prec < aprec:
                    aprec = ctx_prec
                    absprec = 0 # absprec just has to be non-infinite: everything else uses aprec
            else:
                mpz_clear(tmp)
                raise TypeError("cannot coerce from the given integer mod ring (not a power of the same prime)")
        elif isinstance(x, ntl_ZZ_p):
            ctx_prec = ZZ_remove(tmp_z, (<ntl_ZZ>x.modulus()).x, self.prime_pow.pow_ZZ_tmp(1)[0])
            if ZZ_IsOne(tmp_z):
                x = x.lift()
                tmp_Int = PY_NEW(Integer)
                ZZ_to_mpz(tmp_Int.value, &(<ntl_ZZ>x).x)
                x = tmp_Int
                if absprec is infinity or ctx_prec < aprec:
                    aprec = ctx_prec
                    absprec = 0 # absprec just has to be non-infinite: everything else uses aprec
            else:
                raise TypeError("cannot coerce the given ntl_ZZ_p (modulus not a power of the same prime)")
        elif isinstance(x, ntl_ZZ):
            tmp_Int = PY_NEW(Integer)
            ZZ_to_mpz(tmp_Int.value, &(<ntl_ZZ>x).x)
            x = tmp_Int
        elif isinstance(x, int):
            x = Integer(x)
        elif x in parent.residue_field() and x.parent().is_finite():
            # Should only reach here if x is not in F_p
            z = parent.gen()
            poly = x.polynomial().list()
            x = sum([poly[i].lift() * (z ** i) for i in range(len(poly))], parent.zero())
            if absprec is infinity or 1 < aprec:
                aprec = 1
                absprec = 0 # absprec just has to be non-infinite: everything else uses aprec
        cdef pAdicZZpXCRElement _x
        if isinstance(x, Integer):
            if absprec is infinity:
                self._set_from_mpz_rel((<Integer>x).value, rprec)
            else:
                self._set_from_mpz_both((<Integer>x).value, aprec, rprec)
        elif isinstance(x, Rational):
            if absprec is infinity:
                self._set_from_mpq_rel((<Rational>x).value, rprec)
            else:
                self._set_from_mpq_both((<Rational>x).value, aprec, rprec)
        elif isinstance(x, ntl_ZZ_pX):
            if absprec is infinity:
                self._set_from_ZZ_pX_rel(&(<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c, rprec)
            else:
                self._set_from_ZZ_pX_both(&(<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c, aprec, rprec)
        elif isinstance(x, ntl_ZZX):
            if absprec is infinity:
                self._set_from_ZZX_rel((<ntl_ZZX>x).x, rprec)
            else:
                self._set_from_ZZX_both((<ntl_ZZX>x).x, aprec, rprec)
        elif isinstance(x, pAdicExtElement):
            if x.parent() is parent:
                _x = <pAdicZZpXCRElement>x
                if _x.relprec == 0:
                    if absprec is infinity or aprec > _x.ordp:
                        self._set_inexact_zero(_x.ordp) # this works for exact zeros too.
                    else:
                        self._set_inexact_zero(aprec)
                elif _x.relprec < 0:
                    if -_x.relprec < rprec:
                        rprec = _x.relprec
                    else:
                        rprec = -rprec
                    if absprec is infinity or aprec > _x.ordp - rprec:
                        self._set(&_x.unit, _x.ordp, rprec)
                    elif aprec > _x.ordp:
                        self._set(&_x.unit, _x.ordp, _x.ordp - aprec) #negating relprec to indicate non-normalized.
                    else:
                        self._set_inexact_zero(aprec)
                else:
                    if _x.relprec < rprec:
                        rprec = _x.relprec
                    if absprec is infinity or aprec > _x.ordp + rprec:
                        self._set(&_x.unit, _x.ordp, rprec)
                    elif aprec > _x.ordp:
                        self._set(&_x.unit, _x.ordp, aprec - _x.ordp)
                    else:
                        self._set_inexact_zero(aprec)
            elif x.parent().fraction_field() is parent:
                if isinstance(x, pAdicZZpXCRElement):
                    _x = <pAdicZZpXCRElement>x
                    if _x.relprec < 0:
                        _x._normalize()
                    if _x._is_exact_zero():
                        self._set_exact_zero()
                    elif _x._is_inexact_zero():
                        self._set_inexact_zero(_x.ordp)
                    else:
                        if _x.relprec < rprec:
                            rprec = _x.relprec
                        self._set(&_x.unit, _x.ordp, rprec)
                else:
                    # x is a pAdicZZpXCAElement
                    xordp = x.valuation()
                    xprec = x.precision_absolute()
                    if xordp == xprec:
                        self._set_inexact_zero(mpz_get_si((<Integer>xordp).value))
                    else:
                        poly = x._ntl_rep_abs()[0]
                        if absprec is infinity:
                            self._set_from_ZZ_pX_rel(&(<ntl_ZZ_pX>poly).x,(<ntl_ZZ_pX>poly).c, rprec)
                        else:
                            self._set_from_ZZ_pX_both(&(<ntl_ZZ_pX>poly).x,(<ntl_ZZ_pX>poly).c, aprec, rprec)
            elif x.parent() is parent.fraction_field():
                _x = <pAdicZZpXCRElement>x
                if _x.relprec < 0:
                    _x._normalize()
                if _x._is_exact_zero():
                    self._set_exact_zero()
                elif _x._is_inexact_zero():
                    self._set_inexact_zero(_x.ordp)
                else:
                    if _x.relprec < rprec:
                        rprec = _x.relprec
                    self._set(&_x.unit, _x.ordp, rprec)
            else:
                raise NotImplementedError("conversion from different p-adic extensions not yet supported")
        else:
            try:
                x = list(x)
            except TypeError:
                try:
                    x = x.list()
                except AttributeError:
                    raise TypeError("cannot convert x to a p-adic element")
            if absprec is infinity:
                self._set_from_list_rel(x, rprec)
            else:
                self._set_from_list_both(x, aprec, rprec)

    def _cache_key(self):
        r"""
        Return a hashable key which identifies this element.

        This makes it possible to use this element in caches such as
        functions or methods decorated with ``@cached_function`` or
        ``@cached_method`` respectively.

        EXAMPLES:

        In the following example, ``a`` and ``b`` compare equal. They cannot
        have a meaningful hash value since then their hash value would have to
        be the same::

            sage: K.<a> = Qq(9)
            sage: b = a + O(3)
            sage: a == b
            True
            sage: hash(a)
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CR.qAdicCappedRelativeElement'

        However, we want to cache computations which depend on them. Therefore
        they define a ``_cache_key`` which is hashable and uniquely identifies
        them::

            sage: a._cache_key()
            (..., ((0, 1),), 0, 20)
            sage: b._cache_key()
            (..., ((0, 1),), 0, 1)

        TESTS:

        Check that zero values are handled correctly::

            sage: K.zero()._cache_key()
            (..., 0)
            sage: K(0,1)._cache_key()
            (..., 1, 0)
        """
        if self._is_exact_zero():
            return (self.parent(), 0)
        elif self._is_inexact_zero():
            return (self.parent(), 0, self.valuation())
        else:
            return (self.parent(),
                    tuple(tuple(c) if isinstance(c, list) else c
                          for c in self.unit_part().expansion()),
                    self.valuation(), self.precision_relative())

    cdef int _set_inexact_zero(self, long absprec) except -1:
        """
        Set ``self`` to be zero with valuation absprec.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(0,6); z  # indirect doctest
            O(w^6)
            sage: z.valuation()
            6
            sage: z.precision_absolute()
            6
            sage: z.precision_relative()
            0

        TESTS::

            sage: R = Zp(17, 3)
            sage: S.<x> = R[]
            sage: W.<w> = R.ext(x^34 - 289*x^5 + 17)
            sage: z = W(0, 6); z
            O(w^6)
            sage: z.valuation()
            6
            sage: z.precision_absolute()
            6
            sage: z.precision_relative()
            0
        """
        self.ordp = absprec
        self.relprec = 0

    cdef int _set_exact_zero(self) except -1:
        """
        Set ``self`` to be an exact zero.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = R(0); z  # indirect doctest
            0
            sage: z.valuation()
            +Infinity
            sage: z.precision_absolute()
            +Infinity
            sage: z.precision_relative()
            0

        TESTS::

            sage: R = Zp(89, 3)
            sage: S.<x> = R[]
            sage: W.<w> = R.ext(x^34 - 2*89*x^5 + 89)
            sage: z = R(0); z  # indirect doctest
            0
            sage: z.valuation()
            +Infinity
            sage: z.precision_absolute()
            +Infinity
            sage: z.precision_relative()
            0
        """
        self.ordp = maxordp
        self.relprec = 0

    cpdef bint _is_exact_zero(self) except -1:
        """
        Test if ``self`` is an exact zero.

        EXAMPLES::

            sage: R = Qp(3,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(0)
            sage: z._is_exact_zero()
            True
            sage: z = W(0,6)
            sage: z._is_exact_zero()
            False

        TESTS::

            sage: R = Qp(53, 3)
            sage: S.<x> = R[]
            sage: W.<w> = R.ext(x^34 - 2*53^5*x^9 + 53)
            sage: z = W(0)
            sage: z._is_exact_zero()
            True
            sage: z = W(0,6)
            sage: z._is_exact_zero()
            False
        """
        if self.ordp == maxordp:
            return 1
        else:
            return 0

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Test if ``self`` is an inexact zero.

        EXAMPLES::

            sage: R = Zp(7,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(0)
            sage: z._is_inexact_zero()
            False
            sage: z = W(0,6)
            sage: z._is_inexact_zero()
            True

        TESTS::

            sage: R = Qp(29, 3)
            sage: S.<x> = R[]
            sage: W.<w> = R.ext(x^29 - 2*29^5*x - 29)
            sage: z = W(0)
            sage: z._is_inexact_zero()
            False
            sage: z = W(0,6)
            sage: z._is_inexact_zero()
            True
        """
        self._normalize()
        if self.relprec == 0:
            return not self._is_exact_zero()
        else:
            return False

    cdef int _set(self, ZZ_pX_c* unit, long ordp, long relprec) except -1:
        """
        Set ``unit``, ``ordp`` and ``relprec`` directly.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: F = W.fraction_field()
            sage: z = F(1 + w); z  # indirect doctest
            1 + w + O(w^25)

        TESTS::

            sage: R = Zp(17,30)
            sage: S.<x> = R[]
            sage: f = x^51 - 34
            sage: W.<w> = R.ext(f)
            sage: F = W.fraction_field()
            sage: z = F(1 + w); z  # indirect doctest
            1 + w + O(w^1530)
            sage: z = F(w + w^2, relprec=0); z
            O(w)
        """
        self.ordp = ordp
        self._set_prec_rel(relprec)
        if self.relprec != 0:
            ZZ_pX_conv_modulus(self.unit, unit[0], self.prime_pow.get_context_capdiv(relprec).x)

    cdef int _set_from_mpz_rel(self, mpz_t x, long relprec) except -1:
        """
        Set ``self`` from an ``mpz_t`` with relative precision bounded by ``relprec``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, relprec=8)  # indirect doctest
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
            sage: W(70, relprec=0)
            O(w^5)

        TESTS::

            sage: R = Qp(13,50)
            sage: S.<x> = R[]
            sage: f = x^169 - 13
            sage: W.<w> = R.ext(f)
            sage: a = W(65, relprec=8); a.valuation()  # indirect doctest
            169
            sage: W(65, relprec=0)
            O(w^169)
        """
        if mpz_sgn(x) == 0:
            self._set_exact_zero()
            return 0
        cdef mpz_t tmp_m
        cdef ZZ_c tmp_z
        cdef long shift
        mpz_init(tmp_m)
        sig_on()
        shift = mpz_remove(tmp_m, x, self.prime_pow.prime.value)
        sig_off()
        self._set_prec_rel(relprec)
        mpz_to_ZZ(&tmp_z, tmp_m)
        mpz_clear(tmp_m)
        if self.relprec != 0:
            ZZ_pX_SetCoeff(self.unit, 0, ZZ_to_ZZ_p(tmp_z))
            self.ordp = 0
            self._pshift_self(shift)
        else:
            self.ordp = shift * self.prime_pow.e

    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, long relprec) except -1:
        """
        Set ``self`` from an ``mpz_t`` with relative precision bounded by ``relprec``
        and absolute precision bounded by ``absprec``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, 8)  # indirect doctest
            4*w^5 + 3*w^7 + O(w^8)
            sage: W(70, absprec=4)
            O(w^4)

        TESTS::

            sage: R = Zp(7,3)
            sage: S.<x> = R[]
            sage: f = x^49 + 7*x^21 - 14
            sage: W.<w> = R.ext(f)
            sage: W(70, 100)  # indirect doctest
            5*w^49 + 6*w^70 + 3*w^91 + O(w^100)
            sage: W(70, absprec=4)
            O(w^4)
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
        sig_off()
        self.ordp = shift * self.prime_pow.e
        if self._set_prec_both(absprec, relprec) == 1:
            # This indicates that self._set_inexact_zero was called
            mpz_clear(tmp_m)
            return 0
        mpz_to_ZZ(&tmp_z, tmp_m)
        mpz_clear(tmp_m)
        if self.relprec != 0:
            ZZ_pX_SetCoeff(self.unit, 0, ZZ_to_ZZ_p(tmp_z))
            self.ordp = 0
            self._pshift_self(shift)

    cdef int _set_from_mpq_rel(self, mpq_t x, long relprec) except -1:
        """
        Set ``self`` from an ``mpq_t`` with relative precision
        bounded by ``relprec``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(70/3, relprec=9); z  # indirect doctest
            3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + O(w^14)
            sage: z * 3
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + O(w^14)
            sage: W(70)
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23 + 3*w^25 + w^27 + O(w^30)
            sage: F = W.fraction_field()
            sage: y = F(3/700); y
            w^-10 + w^-8 + 4*w^-6 + w^-3 + 4*w^-2 + 3*w^-1 + 3 + 4*w + w^3 + 4*w^4 + w^5 + 4*w^6 + 2*w^7 + 3*w^8 + 4*w^9 + 3*w^10 + 4*w^11 + w^12 + O(w^15)
            sage: y * 700
            3 + O(w^25)
            sage: W(70/3, relprec=0)
            O(w^5)
            sage: c = F(5^-1 + O(5^2)); c
            w^-5 + 3*w^-3 + 2*w^3 + 4*w^5 + 4*w^6 + 3*w^7 + w^9 + O(w^10)
            sage: c * 5
            1 + O(w^15)

        TESTS::

            sage: R = Zp(11, 8, print_mode='digits')
            sage: S.<x> = R[]
            sage: f = x^3 + 1331 * x^2 - 11 * x + 11
            sage: W.<w> = R.ext(f)
            sage: z = W(77/3, relprec=11); repr(z)[3:]
            '304107A2555000'
            sage: repr(z*3)[3:]
            '56698765444000'
            sage: repr(W(77))[3:]
            '5800A6604678856698765444000'
            sage: F = W.fraction_field()
            sage: y = F(3/847); repr(y)[3:]
            '5563A4105291255628.148272'
            sage: repr(y*847)[3:]
            '000000000000000000000003'
            sage: repr(W(77/3, relprec=0))[3:]
            '000'
            sage: c = F(11^-1 + O(11^2)); repr(c)[3:]
            '011111.01A'
            sage: repr(c * 11)[3:]
            '000000001'
        """
        if mpq_sgn(x) == 0:
            self._set_exact_zero()
            return 0
        cdef mpz_t num_unit, den_unit
        self._set_from_mpq_part1(num_unit, den_unit, x)
        self._set_prec_rel(relprec)
        self._set_from_mpq_part2(num_unit, den_unit)

    cdef int _set_from_mpq_both(self, mpq_t x, long absprec, long relprec) except -1:
        """
        Set ``self`` from an ``mpq_t`` with relative precision
        bounded by ``relprec`` and absolute precision bounded by
        ``absprec``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(70/3, 14); z  # indirect doctest
            3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + O(w^14)
            sage: z * 3
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + O(w^14)
            sage: W(70)
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23 + 3*w^25 + w^27 + O(w^30)
            sage: F = W.fraction_field()
            sage: y = F(3/700,-2); y
            w^-10 + w^-8 + 4*w^-6 + w^-3 + O(w^-2)
            sage: y * 700
            3 + O(w^8)
            sage: W(70/3, absprec=4)
            O(w^4)
        """
        if mpq_sgn(x) == 0:
            self._set_inexact_zero(absprec)
            return 0
        cdef mpz_t num_unit, den_unit
        self._set_from_mpq_part1(num_unit, den_unit, x)
        if self._set_prec_both(absprec, relprec) == 1:
            # indicates an inexact zero
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            return 0
        self._set_from_mpq_part2(num_unit, den_unit)

    cdef int _set_from_mpq_part1(self, mpz_t num_unit, mpz_t den_unit, mpq_t x) except -1:
        """
        Set ``num_unit`` to be the unit of the numerator, ``den_unit`` to be
        the unit of the denominator and sets ``self.ordp`` correctly.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(7000/3, 23); z  # indirect doctest
            2*w^15 + 2*w^17 + 3*w^19 + w^22 + O(w^23)
        """
        cdef long num_ordp, den_ordp
        sig_on()
        mpz_init(num_unit)
        mpz_init(den_unit)
        num_ordp = mpz_remove(num_unit, mpq_numref(x), self.prime_pow.prime.value)
        den_ordp = mpz_remove(den_unit, mpq_denref(x), self.prime_pow.prime.value)
        sig_off()
        self.ordp = (num_ordp - den_ordp) * self.prime_pow.e
        if self.ordp < 0 and self.prime_pow.in_field == 0:
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            raise ValueError("p divides the denominator")

    cdef int _set_from_mpq_part2(self, mpz_t num_unit, mpz_t den_unit) except -1:
        """
        Given that ``self.ordp`` and ``self.relprec`` have been set, takes
        ``num_unit`` and ``den_unit`` and sets ``self.unit``.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(QQ(0), 23)  # indirect doctest
            O(w^23)
            sage: W(QQ(0))
            0
        """
        cdef ZZ_c num_zz, den_zz
        cdef ZZ_p_c tmp_zp
        cdef long val = self.ordp / self.prime_pow.e
        cdef mpz_t tmp_m
        if self.relprec != 0:
            mpz_init(tmp_m)
            mpz_set(tmp_m, num_unit)
            mpz_to_ZZ(&num_zz, tmp_m)
            mpz_set(tmp_m, den_unit)
            mpz_to_ZZ(&den_zz, tmp_m)
            mpz_clear(tmp_m)
            #The context has been restored in setting self.relprec
            ZZ_p_div(tmp_zp, ZZ_to_ZZ_p(num_zz), ZZ_to_ZZ_p(den_zz))
            ZZ_pX_SetCoeff(self.unit, 0, tmp_zp)
            self.ordp = 0
            self._pshift_self(val)

    cdef int _set_from_ZZX_rel(self, ZZX_c poly, long relprec) except -1:
        """
        Set ``self`` from a ``ZZX`` with relative precision bounded by
        ``relprec``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZX([4,1,16]), relprec=14); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + 4*w^13 + O(w^14)
            sage: z._ntl_rep()
            [4 1 16]
            sage: z = W(ntl.ZZX([5^40,5^42,3*5^41]), relprec=14); z
            w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + O(w^214)
            sage: W(5)^40 + w*W(5)^42 + w^2 * W(3) * W(5)^41
            w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + 2*w^215 + w^217 + 2*w^218 + w^220 + w^221 + w^222 + 3*w^224 + O(w^225)
            sage: z = W(ntl.ZZX([5^40,5^42,3*5^41]), relprec=0); z
            O(w^200)
        """
        if ZZX_IsZero(poly):
            self._set_exact_zero()
            return 0
        if ZZX_deg(poly) >= self.prime_pow.deg:
            raise NotImplementedError
        # the -1 in the next line signals that there is no absprec specified
        self._set_from_ZZX_part1(poly, -1, relprec)
        # context was restored in _set_from_ZZX_part1
        if relprec == 0:
            self._set_prec_rel(relprec)
            return 0
        if self.relprec + self.ordp != 0:
            self.prime_pow.restore_context_capdiv(self.relprec + self.ordp)
            ZZX_to_ZZ_pX(self.unit, poly)
            self._internal_lshift(-self.ordp)

    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, long relprec) except -1:
        """
        Set ``self`` from a ``ZZX`` with relative precision bounded by
        ``relprec`` and absolute precision bounded by ``absprec``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZX([4,1,16]), 12); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + O(w^12)
            sage: z._ntl_rep()
            [4 1 16]
            sage: z = W(ntl.ZZX([5^40,5^42,3*5^41]), 212); z
            w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + O(w^212)
            sage: z = W(ntl.ZZX([5^40,5^42,3*5^41]), 197); z
            O(w^197)
        """
        if ZZX_IsZero(poly) or absprec <= 0:
            self._set_inexact_zero(absprec)
            return 0
        if ZZX_deg(poly) >= self.prime_pow.deg:
            raise NotImplementedError
        if self._set_from_ZZX_part1(poly, absprec, relprec) == -2:
            # indicates _set_inexact_zero was called
            return 0
        # context was restored in _set_from_ZZX_part1
        if self.relprec + self.ordp != 0:
            self.prime_pow.restore_context_capdiv(self.relprec + self.ordp)
            ZZX_to_ZZ_pX(self.unit, poly)
            self._internal_lshift(-self.ordp)

    cdef int _set_from_ZZX_part1(self, ZZX_c poly, long absprec, long relprec) except -1:
        """
        Set ``self.ordp`` from ``poly`` and restores the context. ``poly`` must
        have degree less than ``self.prime_pow.deg``.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZX([4,1,16]), 12); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + O(w^12)
        """
        cdef long i = 0
        cdef long deg = ZZX_deg(poly)
        cdef long mini = -1
        cdef long minval
        cdef long curval
        cdef ZZ_c tmp_z
        while mini == -1:
            if not ZZ_IsZero(ZZX_coeff(poly,i)):
                minval = ZZ_remove(tmp_z, ZZX_coeff(poly, i), self.prime_pow.pow_ZZ_tmp(1)[0])
                mini = i
            i += 1
        while i <= deg:
            if not ZZ_IsZero(ZZX_coeff(poly,i)):
                curval = ZZ_remove(tmp_z, ZZX_coeff(poly, i), self.prime_pow.pow_ZZ_tmp(1)[0])
                if curval < minval:
                    minval = curval
                    mini = i
            i += 1
        if self.prime_pow.e == 1:
            self.ordp = minval
        else:
            self.ordp = minval * self.prime_pow.e + mini
        if absprec == -1: # indicates that _set_from_ZZX_rel is calling
            self._set_prec_rel(relprec)
        elif self._set_prec_both(absprec, relprec):
            # indicates self._set_inexact_zero was called
            return -2
        # _set_prec_rel or both has restored the context so that part2 works.

    cdef int _set_from_ZZ_pX_rel(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long relprec) except -1:
        """
        Set ``self`` from a ``ZZ_pX`` with relative precision bounded by
        ``relprec``.

        If ``ctx`` is ``None`` and ``poly`` is 0 this function will raise an
        error (a ``ZZ_pX`` cannot represent something with infinite absolute
        precision).

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2)); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + O(w^10)
            sage: z._ntl_rep()
            [4 1 16]
            sage: z = W(ntl.ZZ_pX([5^40,5^42,3*5^41], 5^44)); z
            w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + 2*w^215 + w^217 + 2*w^218 + O(w^220)
            sage: z = W(ntl.ZZ_pX([5^40,5^42,3*5^41], 5^44), relprec=0); z
            O(w^200)
        """
        cdef long ctx_prec = -1
        if ctx is not None:
            ctx_prec = self._check_ZZ_pContext(ctx) * self.prime_pow.e
        if ZZ_pX_IsZero(poly[0]):
            if ctx_prec == -1:
                raise ValueError("must specify either a context or an absolute precision bound")
            else:
                self._set_inexact_zero(ctx_prec)
            return 0
        self._set_from_ZZ_pX_part1(poly)
        if relprec == 0:
            self._set_prec_rel(relprec)
            return 0
        if ctx_prec == -1:
            self._set_prec_rel(self.ordp + relprec)
        else:
            self._set_prec_rel(min(ctx_prec, self.ordp + relprec))
        self._set_from_ZZ_pX_part2(poly)

    cdef int _set_from_ZZ_pX_both(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec, long relprec) except -1:
        """
        Set ``self`` from a ``ZZ_pX`` with relative precision bounded by
        ``relprec`` and absolute precision bounded by ``absprec``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2), absprec=8, relprec=12); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + O(w^8)
            sage: z._ntl_rep()
            [4 1 16]
            sage: z = W(ntl.ZZ_pX([5^40,5^42,3*5^41], 5^50), 220); z
            w^200 + 4*w^207 + 4*w^209 + w^210 + 2*w^211 + 2*w^213 + 2*w^215 + w^217 + 2*w^218 + O(w^220)
            sage: z = W(ntl.ZZ_pX([5^40,5^42,3*5^41], 5^44), absprec=77); z
            O(w^77)
        """
        cdef long ctx_prec
        if ctx is not None:
            ctx_prec = self._check_ZZ_pContext(ctx)
            if ctx_prec * self.prime_pow.e < absprec:
                absprec = ctx_prec * self.prime_pow.e
        if ZZ_pX_IsZero(poly[0]):
            self._set_inexact_zero(absprec)
            return 0
        self._set_from_ZZ_pX_part1(poly)
        if absprec <= self.ordp:
            self._set_inexact_zero(absprec)
        else:
            self._set_prec_rel(min(absprec, self.ordp + relprec))
            self._set_from_ZZ_pX_part2(poly)

    cdef int _set_from_ZZ_pX_part1(self, ZZ_pX_c* poly) except -1:
        """
        Set ``self.ordp`` based on ``poly``.  ``poly`` must not be 0.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2), absprec=8, relprec=12); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + O(w^8)
        """
        cdef long val = 0, index = 0
        ZZ_pX_min_val_coeff(val, index, poly[0], self.prime_pow.pow_ZZ_tmp(1)[0])
        if self.prime_pow.e == 1:
            self.ordp = val
        else:
            self.ordp = val * self.prime_pow.e + index

    cdef int _set_from_ZZ_pX_part2(self, ZZ_pX_c* poly) except -1:
        """
        Assuming that ``self.ordp`` and ``self.relprec`` have been set, sets
        ``self.unit`` to ``poly`` and then normalizes.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2), absprec=8, relprec=12); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + O(w^8)
        """
        # We've set self.relprec to what is actually the absolute precision.
        if self.relprec != 0:
            ZZ_pX_conv_modulus(self.unit, poly[0], self.prime_pow.get_context_capdiv(self.relprec).x)
            self.relprec -= self.ordp
            self._internal_lshift(-self.ordp)

    cdef bint _set_prec_rel(self, long relprec) except -1:
        """
        Safely set the relative precision of ``self`` to be the absolute
        value of ``relprec``.

        Returns ``True`` iff ``self.relprec`` was reset.

        Note that this will wipe out anything in ``self.unit``.  Be
        careful resetting ``self.unit`` directly: if you set it to a
        different modulus, NTL may have problems.  The safest way to
        reset ``self.unit`` to a different modulus is::

            self.prime_pow.restore_context_capdiv(self.relprec)
            cdef ZZ_pX_c tmp = self.unit
            self._set_prec_rel(new_rel_prec)
            ZZ_pX_conv_modulus(self.unit, tmp, self.prime_pow.get_context_capdiv(self.relprec).x)

        You may be able to just set ``self.relprec`` and
        ``ZZ_pX_conv_modulus`` if you're decreasing precision.  I'm
        not sure.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, relprec=8)  # indirect doctest
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
        """
        if self.relprec == relprec:
            return False
        if relprec != 0:
            self.prime_pow.restore_context_capdiv(relprec)
            self.unit = ZZ_pX_c()
        self.relprec = relprec
        return True

    cdef bint _set_prec_both(self, long absprec, long relprec) except -1:
        """
        Assuming ``self.ordp`` is set, sets the relative precision of ``self``
        to the minimum of ``abs(relprec)`` and ``absprec-self.ordp``.

        If ``relprec`` is negative, will set ``self.relprec`` to be negative
        (indicating unnormalized unit)

        Returns ``True`` iff ``self.relprec = 0``, ie ``self`` was set to an
        inexact zero.

        Note that this will wipe out anything in ``self.unit``.  Be
        careful resetting ``self.unit`` directly: if you set it to a
        different modulus, NTL may have problems.  The safest way to
        reset ``self.unit`` to a different modulus is:

            self.prime_pow.restore_context_capdiv(self.relprec)
            cdef ZZ_pX_c tmp = self.unit
            self._set_prec_rel(new_rel_prec)
            ZZ_pX_conv_modulus(self.unit, tmp, self.prime_pow.get_context_capdiv(self.relprec).x)

        You may be able to just set ``self.relprec`` and
        ``ZZ_pX_conv_modulus`` if you're decreasing precision.  I'm
        not sure.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70, 8)  # indirect doctest
            4*w^5 + 3*w^7 + O(w^8)
        """
        self.relprec = absprec - self.ordp
        cdef long arelprec
        if relprec < 0:
            arelprec = -relprec
        else:
            arelprec = relprec
        if self.relprec <= 0:
            self._set_inexact_zero(absprec)
        else:
            if arelprec < self.relprec:
                self.relprec = arelprec
            if self.relprec != 0:
                self.prime_pow.restore_context_capdiv(self.relprec)
                self.unit = ZZ_pX_c()
                if relprec < 0:
                    self.relprec = -self.relprec
        return self.relprec == 0

    cdef int _normalize(self) except -1:
        """
        Normalize ``self``, adjusting ``self.ordp``, ``self.relprec``, and
        ``self.unit`` so that ``self.unit`` actually represents a unit.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1+w)^5
            sage: y = z - 1
            sage: y._ntl_rep_unnormalized()
            [5 3005 25 3060 5]
            sage: y  # indirect doctest
            w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: y._ntl_rep_unnormalized()
            [41 26 152 49 535]
        """
        cdef long minval = 0, mini = 0, shift
        if self.relprec < 0:
            if ZZ_pX_IsZero(self.unit):
                self.ordp -= self.relprec # note that self.relprec < 0
                self.relprec = 0
            else:
                ZZ_pX_min_val_coeff(minval, mini, self.unit, self.prime_pow.pow_ZZ_tmp(1)[0])
                if self.prime_pow.e == 1:
                    shift = minval
                else:
                    shift = minval * self.prime_pow.e + mini
                if shift >= -self.relprec:
                    self.ordp -= self.relprec # note that self.relprec < 0
                    self.relprec = 0
                elif shift > 0:
                    self.relprec = -self.relprec - shift
                    self.ordp += shift
                    self._internal_lshift(-shift)
                else:
                    self.relprec = -self.relprec

    def _is_normalized(self):
        """
        Whether this element is currently normalized.

        EXAMPLES::

            sage: R.<a> = ZqCR(125,implementation="NTL"); b = 5*a + 4; c = 10*a^2 + 6; d = b + c
            sage: d._is_normalized()
            False
            sage: d.valuation()
            1
            sage: d._is_normalized()
            True
        """
        return self.relprec >= 0

    cdef int _internal_lshift(self, long shift) except -1:
        """
        Multiply ``self.unit`` by ``x^shift``.

        Note that ``self.relprec`` must be set before calling this
        function and should not be 0, and self.unit must be defined to
        precision ``self.relprec - shift``

        This function does not alter ``self.ordp`` even though it WILL
        change the valuation of ``self.unit``

        Also note that if you call this function you should usually
        manually set ``self.relprec = -self.relprec`` since this function
        will usually unnormalize ``self``.

        TESTS::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1+w)^5
            sage: y = z - 1
            sage: y  # indirect doctest
            w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
        """
        if self.relprec == 0:
            raise ValueError("p-adic Internal l-shift called with relative precision 0")
        cdef ZZ_pX_c tmpP
        cdef ZZ_pX_Modulus_c* mod
        if self.prime_pow.e == 1:
            if shift > 0:
                ZZ_pX_left_pshift(self.unit, self.unit, self.prime_pow.pow_ZZ_tmp(shift)[0], self.prime_pow.get_context(self.relprec).x)
            else:
                ZZ_pX_right_pshift(self.unit, self.unit, self.prime_pow.pow_ZZ_tmp(-shift)[0], self.prime_pow.get_context(self.relprec).x)
        else:
            if shift > 0:
                self.prime_pow.restore_context_capdiv(self.relprec)
                mod = self.prime_pow.get_modulus_capdiv(self.relprec)
                ZZ_pX_PowerXMod_long_pre(tmpP, shift, mod[0])
                ZZ_pX_MulMod_pre(self.unit, self.unit, tmpP, mod[0])
            elif shift < 0:
                self.prime_pow.eis_shift_capdiv(&self.unit, &self.unit, -shift, self.relprec)

    cdef int _pshift_self(self, long shift) except -1:
        """
        Multiply ``self`` by ``p^shift``.

        This function assumes that ``self.relprec``, ``self.ordp`` and
        ``self.unit`` are already set (in the case ``self.prime_pow.e
        != 1``), and is more reasonable to call externally than
        ``_internal_lshift``

        EXAMPLES::

            sage: R = Qp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(3/25, relprec=6); z
            3*w^-10 + 3*w^-8 + 2*w^-6 + O(w^-4)
            sage: z * 25
            3 + O(w^6)
        """
        cdef ZZ_pX_c high_shifter
        cdef ZZ_pX_Modulus_c *modulus
        cdef ntl_ZZ_pContext_class c
        cdef ZZ_pX_c* high_array
        cdef long i, high_length
        if self.prime_pow.e == 1:
            self.ordp += shift
        else:
            self.ordp += shift * self.prime_pow.e
            if shift < 0:
                shift = -shift
                c = self.prime_pow.get_context_capdiv(self.relprec)
                c.restore_c()
                modulus = self.prime_pow.get_modulus_capdiv(self.relprec)
                if isinstance(self.prime_pow, PowComputer_ZZ_pX_big_Eis):
                    high_array = (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_shifter
                elif isinstance(self.prime_pow, PowComputer_ZZ_pX_small_Eis):
                    high_array = (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_shifter
                else:
                    raise TypeError("unrecognized PowComputer type")
                ZZ_pX_conv_modulus(high_shifter, high_array[0], c.x)
                ZZ_pX_InvMod_newton_ram(high_shifter, high_shifter, modulus[0], c.x)
                ZZ_pX_PowerMod_long_pre(high_shifter, high_shifter, shift, modulus[0])
                ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus[0])
            elif shift > 0:
                i = 0
                c = self.prime_pow.get_context_capdiv(self.relprec)
                c.restore_c()
                modulus = self.prime_pow.get_modulus_capdiv(self.relprec)
                if isinstance(self.prime_pow, PowComputer_ZZ_pX_big_Eis):
                    high_array = (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_shifter
                    high_length = (<PowComputer_ZZ_pX_big_Eis>self.prime_pow).high_length
                elif isinstance(self.prime_pow, PowComputer_ZZ_pX_small_Eis):
                    high_array = (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_shifter
                    high_length = (<PowComputer_ZZ_pX_small_Eis>self.prime_pow).high_length
                else:
                    raise TypeError("unrecognized PowComputer type")
                if shift >= self.prime_pow.prec_cap:
                    # high_shifter = p^(2^(high_length - 1))/x^(e*2^(high_length - 1))
                    ZZ_pX_conv_modulus(high_shifter, high_array[high_length-1], c.x)
                    # if shift = r + s * 2^(high_length - 1)
                    # then high_shifter = p^(s*2^(high_length - 1))/x^(e*s*2^(high_length - 1))
                    ZZ_pX_PowerMod_long_pre(high_shifter, high_shifter, (shift / (1L << (high_length - 1))), modulus[0])
                    ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus[0])
                    # Now we only need to multiply self.unit by p^r/x^(e*r) where r < 2^(high_length - 1), which is tractable.
                    shift = shift % (1L << (high_length - 1))
                while shift > 0:
                    if shift & 1:
                        ZZ_pX_conv_modulus(high_shifter, high_array[i], c.x)
                        ZZ_pX_MulMod_pre(self.unit, self.unit, high_shifter, modulus[0])
                    shift = shift >> 1
                    i += 1

    cdef pAdicZZpXCRElement _new_c(self, long relprec):
        """
        Return a new element with the same parent as ``self`` and
        relative precision ``relprec``

        Note that if ``relprec`` is nonpositive, the convention is that
        ``relprec = 0`` indicates an exact or inexact zero, ``relprec < 0``
        indicates an unnormalized element.

        EXAMPLES::

            sage: R = Qp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w^5 + 1  # indirect doctest
            1 + w^5 + O(w^25)
        """
        cdef pAdicZZpXCRElement ans = pAdicZZpXCRElement.__new__(pAdicZZpXCRElement)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        if relprec > 0:
            self.prime_pow.restore_context_capdiv(relprec)
            ans.relprec = relprec
        elif relprec == 0:
            ans._set_exact_zero()
        else:
            self.prime_pow.restore_context_capdiv(-relprec)
            ans.relprec = relprec
        return ans

    def __reduce__(self):
        """
        Pickles ``self``.

        EXAMPLES::

            sage: R = Qp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5 - 1
            sage: loads(dumps(z)) == z
            True
        """
        cdef Integer relprec, ordp
        relprec = PY_NEW(Integer)
        ordp = PY_NEW(Integer)
        mpz_set_si(relprec.value, self.relprec)
        mpz_set_si(ordp.value, self.ordp)
        if self.relprec == 0:
            return make_ZZpXCRElement, (self.parent(), None, ordp, relprec, 0)
        self.prime_pow.restore_context_capdiv(self.relprec)
        cdef ntl_ZZ_pX holder = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        holder.c = self.prime_pow.get_context_capdiv(self.relprec)
        holder.x = self.unit
        return make_ZZpXCRElement, (self.parent(), holder, ordp, relprec, 0)

    cdef int _cmp_units(left, pAdicGenericElement right) except -2:
        """
        For units ``left`` and ``right``, return 0 if they are equal up to
        the lesser of the two precisions, or 1 if they are not.

        EXAMPLES::

            sage: R = Qp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w == 1  # indirect doctest
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
        cdef pAdicZZpXCRElement diff = <pAdicZZpXCRElement> (left - right)
        diff._normalize()
        if diff.relprec == 0:
            return 0
        # for now, just return 1
        return 1

    def __invert__(self):
        """
        Return the inverse of this element.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: y = ~z; y  # indirect doctest
            1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14 + 3*w^16 + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24 + O(w^25)
            sage: y.parent()
            5-adic Eisenstein Extension Field in w defined by x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: z = z - 1
            sage: ~z
            w^-5 + 4*w^-4 + 4*w^-3 + 4*w^-2 + 2*w^-1 + 1 + w + 4*w^2 + 4*w^3 + 4*w^4 + w^5 + w^6 + w^7 + 4*w^8 + 4*w^9 + 2*w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^14 + O(w^15)
            sage: ~z * z
            1 + O(w^20)
        """
        if self._is_exact_zero():
            raise ZeroDivisionError("cannot divide by zero")
        if self._is_inexact_zero(): # this calls _normalize
            raise PrecisionError("cannot divide by something indistinguishable from zero")
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        if not ans.prime_pow.in_field:
            ans._parent = self._parent.fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        ans.ordp = -self.ordp
        sig_on()
        if self.prime_pow.e == 1:
            ZZ_pX_InvMod_newton_unram(ans.unit, self.unit, self.prime_pow.get_modulus(ans.relprec)[0], self.prime_pow.get_context(ans.relprec).x, self.prime_pow.get_context(1).x)
        else:
            ZZ_pX_InvMod_newton_ram(ans.unit, self.unit, self.prime_pow.get_modulus_capdiv(ans.relprec)[0], self.prime_pow.get_context_capdiv(ans.relprec).x)
        sig_off()
        return ans

    cdef pAdicZZpXCRElement _lshift_c(self, long n):
        """
        Multiply ``self`` by the uniformizer raised to the power ``n``.  If
        ``n`` is negative, right shifts by ``-n``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: z << 17  # indirect doctest
            w^17 + w^22 + w^23 + 2*w^24 + 4*w^25 + 3*w^27 + w^29 + 4*w^30 + 4*w^31 + 4*w^32 + 4*w^33 + 4*w^34 + 4*w^37 + w^38 + 4*w^41 + O(w^42)
            sage: z << (-1)
            w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
        """
        if not self.prime_pow.in_field and n < -self.ordp:
            return self._rshift_c(-n)
        check_ordp(n)
        cdef pAdicZZpXCRElement ans
        if self._is_exact_zero() or n == 0:
            return self
        elif self._is_inexact_zero():
            ans = self._new_c(0)
        else:
            ans = self._new_c(self.relprec)
            ans.unit = self.unit
        ans.ordp = self.ordp + n
        check_ordp(ans.ordp)
        return ans

    def __lshift__(pAdicZZpXCRElement self, shift):
        """
        Multiply ``self`` by the uniformizer raised to the power ``n``.  If
        ``n`` is negative, right shifts by ``-n``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: z << 17  # indirect doctest
            w^17 + w^22 + w^23 + 2*w^24 + 4*w^25 + 3*w^27 + w^29 + 4*w^30 + 4*w^31 + 4*w^32 + 4*w^33 + 4*w^34 + 4*w^37 + w^38 + 4*w^41 + O(w^42)
            sage: z << (-1)
            w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + O(w^24)
        """
        cdef pAdicZZpXCRElement ans
        if not isinstance(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            if self._is_exact_zero():
                return self
            if self.prime_pow.in_field or mpz_sgn((<Integer>shift).value) > 0:
                raise ValueError("shift does not fit in long")
            else:
                ans = self._new_c(0)
                ans.ordp = 0
                return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicZZpXCRElement _rshift_c(self, long n):
        """
        Divide ``self`` by the uniformizer raised to the power ``n``.  If
        parent is not a field, throws away the nonpositive part of
        the series expansion.  If ``n`` is negative, left shifts by ``-n``.

        EXAMPLES::

            sage: R = Zp(5,5,print_mode='digits')
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: for m in range(26): repr(z >> m)  # indirect doctest
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
            '...40014004444410304211000010000'
        """
        if self.prime_pow.in_field or n <= self.ordp:
            return self._lshift_c(-n)
        if self._is_exact_zero() or n == 0:
            return self
        cdef long arelprec
        if self.relprec < 0:
            arelprec = -self.relprec
        else:
            arelprec = self.relprec
        cdef pAdicZZpXCRElement ans
        if arelprec > n - self.ordp:
            ans = self._new_c(arelprec - (n - self.ordp))
            if self.prime_pow.e == 1:
                ZZ_pX_right_pshift(ans.unit, self.unit, self.prime_pow.pow_ZZ_tmp(n - self.ordp)[0], self.prime_pow.get_context(ans.relprec).x)
            else:
                self.prime_pow.eis_shift_capdiv(&ans.unit, &self.unit, n - self.ordp, ans.relprec)
        else:
            ans = self._new_c(0)
        ans.ordp = 0
        ans.relprec = -ans.relprec
        return ans

    def __rshift__(pAdicZZpXCRElement self, shift):
        """
        Divide ``self`` by the uniformizer raised to the power ``n``.  If
        parent is not a field, throws away the nonpositive part of
        the series expansion.  If ``n`` is negative, left shifts by ``-n``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: z >> (6)  # indirect doctest
            1 + 2*w + 4*w^2 + 3*w^4 + w^6 + 4*w^7 + 4*w^8 + 4*w^9 + 4*w^10 + 4*w^11 + 4*w^14 + w^15 + 4*w^18 + O(w^19)
            sage: z >> (-4)
            w^4 + w^9 + w^10 + 2*w^11 + 4*w^12 + 3*w^14 + w^16 + 4*w^17 + 4*w^18 + 4*w^19 + 4*w^20 + 4*w^21 + 4*w^24 + w^25 + 4*w^28 + O(w^29)
            sage: F = W.fraction_field()
            sage: z = F(z)
            sage: z >> 7
            w^-7 + w^-2 + w^-1 + 2 + 4*w + 3*w^3 + w^5 + 4*w^6 + 4*w^7 + 4*w^8 + 4*w^9 + 4*w^10 + 4*w^13 + w^14 + 4*w^17 + O(w^18)
        """
        cdef pAdicZZpXCRElement ans
        if not isinstance(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            if self._is_exact_zero():
                return self
            if self.prime_pow.in_field or mpz_sgn((<Integer>shift).value) < 0:
                raise ValueError("valuation overflow")
            else:
                ans = self._new_c(0)
                ans.ordp = 0
                return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cpdef _neg_(self):
        """
        Negation.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5; z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: -z  # indirect doctest
            4 + 3*w^5 + 4*w^6 + w^7 + w^8 + w^9 + w^10 + w^11 + 2*w^12 + 4*w^13 + 4*w^15 + 3*w^16 + w^17 + 2*w^18 + 3*w^19 + 2*w^21 + 4*w^23 + 4*w^24 + O(w^25)
            sage: y = z + (-z); y
            O(w^25)
            sage: -y
            O(w^25)
            sage: -W(0)
            0
        """
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        ans.ordp = self.ordp
        if self.relprec != 0:
            self.prime_pow.restore_context_capdiv(self.relprec)
            ZZ_pX_negate(ans.unit, self.unit)
        return ans

#                                             / 1 + \alpha^p \pi_K^{p \lambda}                      mod \mathfrak{p}_K^{p \lambda + 1}   if 1 \le \lambda < \frac{e_K}{p-1}
#        (1 + \alpha \pi^{\lambda})^p \equiv {  1 + (\alpha^p - \epsilon \alpha) \pi_K^{p \lambda}  mod \mathfrak{p}_K^{p \lambda + 1}   if \lambda = \frac{e_K}{p-1}
#                                             \ 1 - \epsilon \alpha \pi_K^{\lambda + e}             mod \mathfrak{p}_K^{\lambda + e + 1} if \lambda > \frac{e_K}{p-1}

    def __pow__(pAdicZZpXCRElement self, _right, m): # m ignored
        r"""
        Compute ``self^right``.

        Note: when ``right`` is divisible by `p` then one can get more
        precision than expected.

        Lemma 2.1 (Constructing Class Fields over Local Fields, Sebastian Pauli):

        Let `\alpha` be in `\mathcal{O}_K`.  Let

        .. MATH::

            p = -\pi_K^{e_K} \epsilon

        be the factorization of `p` where `\epsilon` is a unit.  Then
        the `p`-th power of `1 + \alpha \pi_K^{\lambda}` satisfies

        .. MATH::

            (1 + \alpha \pi^{\lambda})^p \equiv \left{ \begin{array}{lll}
            1 + \alpha^p \pi_K^{p \lambda} & \mod \mathfrak{p}_K^{p \lambda + 1} & \mbox{if $1 \le \lambda < \frac{e_K}{p-1}$} \\
            1 + (\alpha^p - \epsilon \alpha) \pi_K^{p \lambda} & \mod \mathfrak{p}_K^{p \lambda + 1} & \mbox{if $\lambda = \frac{e_K}{p-1}$} \\
            1 - \epsilon \alpha \pi_K^{\lambda + e} & \mod \mathfrak{p}_K^{\lambda + e + 1} & \mbox{if $\lambda > \frac{e_K}{p-1}$}
            \end{array} \right.


        So if ``right`` is divisible by `p^k` we can multiply the
        relative precision by `p` until we exceed `e/(p-1)`, then add
        `e` until we have done a total of `k` things: the precision of
        the result can therefore be greater than the precision of
        ``self``.

        There is also the issue of `p`-adic exponents, and determining
        how the precision of the exponent affects the precision of the
        result.

        In computing `(a + O(\pi^k))^{b + O(p^m)}`, one needs that the
        reduction of `a` mod `\pi` is in the prime field
        `\GF{p}` (so that the `p^m` power of the Teichmuller
        part is constant as `m` increases).  Given this restriction,
        we can factor out the Teichmuller part and use the above lemma
        to find the first spot where

        .. MATH::

            (1 + \alpha \pi^{\lambda})^{p^m}

        differs from 1.  We compare this with the precision bound
        given by computing `(a + O(\pi^k))^b` and take the lesser of
        the two.

        In order to do this we need to compute the valuation of ``(self
        / self.parent().teichmuller(self)) - 1``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: (1 + w)^5  # indirect doctest
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24 + O(w^25)
            sage: (1 + w)^-5
            1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14 + 3*w^16 + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24 + O(w^25)
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
            sage: R = Zp(2, 10)
            sage: S.<x> = R[]
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

        TESTS:

            We define ``0^0`` to be unity, :issue:`13786`::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: type(W(0))
            <class 'sage.rings.padics.padic_ZZ_pX_CR_element.pAdicZZpXCRElement'>
            sage: W(0)^0
            1 + O(w^25)
            sage: W(0)^0 == W(1)
            True

        The value returned from ``0^0`` should belong to our ring::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: type(W(0)^0) == type(W(0))
            True
        """
        self._normalize()
        cdef Integer right
        cdef bint padic_exp
        cdef long exp_prec
        cdef long exp_val
        cdef long relprec
        cdef long threshold # e / (p-1)
        cdef mpz_t tmp, tmp2
        if mpz_fits_slong_p(self.prime_pow.prime.value) == 0:
            threshold = 0
        else:
            threshold = self.prime_pow.e / (mpz_get_si(self.prime_pow.prime.value) - 1)
        cdef Integer base_level
        cdef pAdicZZpXCRElement ans
        cdef long i
        if self._is_exact_zero():
            # Return 0 except for 0^0 error or type error on the exponent.
            if isinstance(_right, (Integer, Rational, int)) or (isinstance(_right, pAdicGenericElement) and _right._is_base_elt(self.prime_pow.prime)):
                if _right == 0:
                    return self.parent(1)
                return self
            else:
                raise TypeError("exponent must be an integer, rational or base p-adic with the same prime")
        elif self._is_inexact_zero():
            # If an integer exponent, return an inexact zero of valuation right * self.ordp.  Otherwise raise an error.
            if isinstance(_right, int):
                _right = Integer(_right)
            if isinstance(_right, Integer):
                ans = self._new_c(0)
                mpz_init_set_si(tmp, self.ordp)
                mpz_mul(tmp, tmp, (<Integer>_right).value)
                if mpz_cmp_si(tmp, maxordp) >= 0 or mpz_cmp_si(tmp, -maxordp) <= 0:
                    raise ValueError("valuation overflow")
                ans.ordp = mpz_get_si(tmp)
                mpz_clear(tmp)
                return ans
            elif isinstance(_right, Rational) or (isinstance(_right, pAdicGenericElement) and _right._is_base_elt(self.prime_pow.prime)):
                raise ValueError("Need more precision")
            else:
                raise TypeError("exponent must be an integer, rational or base p-adic with the same prime")
        if isinstance(_right, int):
            _right = Integer(_right)
        if isinstance(_right, Integer):
            right = <Integer> _right
            if right == 0:
                # return 1 to maximum precision
                ans = self._new_c(self.prime_pow.ram_prec_cap)
                ans.ordp = 0
                ZZ_pX_SetCoeff_long(ans.unit, 0, 1)
                return ans
            padic_exp = False
            exp_val = _right.valuation(self.prime_pow.prime) ##
        elif isinstance(_right, pAdicGenericElement) and _right._is_base_elt(self.prime_pow.prime):
            if self.ordp != 0:
                raise ValueError("in order to raise to a p-adic exponent, base must be a unit")
            right = Integer(_right)
            padic_exp = True
            exp_prec = _right.precision_absolute() ##
            exp_val = _right.valuation() ##
            if exp_val < 0:
                raise NotImplementedError("negative valuation exponents not yet supported")
            # checks to see if the residue of self.unit is in the prime field.
            if self.prime_pow.e == 1:
                for i from 1 <= i <= ZZ_pX_deg(self.unit):
                    if not ZZ_divide_test(ZZ_p_rep(ZZ_pX_coeff(self.unit, i)), self.prime_pow.pow_ZZ_tmp(1)[0]):
                        raise ValueError("in order to raise to a p-adic exponent, base must reduce to an element of F_p mod the uniformizer")
            # compute the "level"
            teich_part = self.parent().teichmuller(self)
            base_level = (self / teich_part - 1).valuation() ##
        elif isinstance(_right, Rational):
            raise NotImplementedError
        else:
            raise TypeError("exponent must be an integer, rational or base p-adic with the same prime")
        # Now we compute the increased relprec due to the exponent having positive p-adic valuation
        if exp_val > 0:
            mpz_init_set_si(tmp, self.relprec)
            while mpz_cmp_si(tmp, threshold) <= 0 and exp_val > 0:
                mpz_mul(tmp, tmp, self.prime_pow.prime.value)
                exp_val -= 1
            if exp_val > 0:
                mpz_init_set_si(tmp2, self.prime_pow.e)
                mpz_addmul_ui(tmp, tmp2, exp_val)
                mpz_clear(tmp2)
            if mpz_cmp_si(tmp, self.prime_pow.ram_prec_cap) > 0:
                relprec = self.prime_pow.ram_prec_cap
            else:
                relprec = mpz_get_si(tmp)
            mpz_clear(tmp)
        else:
            relprec = self.relprec
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
                if mpz_cmp_si(base_level.value, relprec) < 0:
                    relprec = mpz_get_si(base_level.value)
            else:
                ans = self._new_c(0)
                ans.ordp = 0
                return ans
        ans = self._new_c(relprec)
        if self.ordp == 0:
            ans.ordp = 0
        else:
            mpz_init_set(tmp, right.value)
            mpz_mul_si(tmp, tmp, self.ordp)
            if mpz_cmp_si(tmp, maxordp) >= 0 or mpz_cmp_si(tmp, -maxordp) <= 0:
                raise ValueError("valuation overflow")
            ans.ordp = mpz_get_si(tmp)
            mpz_clear(tmp)
        cdef ntl_ZZ rZZ = ntl_ZZ.__new__(ntl_ZZ)
        mpz_to_ZZ(&rZZ.x, right.value)
        sig_on()
        if mpz_sgn(right.value) < 0:
            if self.prime_pow.e == 1:
                ZZ_pX_InvMod_newton_unram(ans.unit, self.unit, self.prime_pow.get_modulus(ans.relprec)[0], self.prime_pow.get_context(ans.relprec).x, self.prime_pow.get_context(1).x)
            else:
                ZZ_pX_InvMod_newton_ram(ans.unit, self.unit, self.prime_pow.get_modulus_capdiv(ans.relprec)[0], self.prime_pow.get_context_capdiv(ans.relprec).x)
            ZZ_negate(rZZ.x, rZZ.x)
            ZZ_pX_PowerMod_pre(ans.unit, ans.unit, rZZ.x, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
        else:
            ZZ_pX_PowerMod_pre(ans.unit, self.unit, rZZ.x, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
        sig_off()
        return ans

    cpdef _add_(self, _right):
        """
        Compute the sum of ``self`` and ``right``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)) - 69  # indirect doctest
            1 + O(w^13)
            sage: -69 + (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13))
            1 + O(w^13)
            sage: y = (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13))
            sage: y - 70
            O(w^13)
            sage: y + 0
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + O(w^13)
        """
        cdef pAdicZZpXCRElement right = <pAdicZZpXCRElement>_right
        cdef pAdicZZpXCRElement ans
        cdef long tmpL
        cdef ZZ_pX_c tmpP
        if self.relprec == 0:
            if self.ordp >= right.ordp + right.relprec: # or self._is_exact_zero()
                return right
            elif self.ordp <= right.ordp:
                ans = self._new_c(0)
                ans.ordp = self.ordp
            else:
                ans = self._new_c(self.ordp - right.ordp)
                ZZ_pX_conv_modulus(ans.unit, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                if right.relprec < 0:
                    ans.relprec = -ans.relprec
                ans.ordp = right.ordp
            return ans
        if right.relprec == 0:
            if right.ordp >= self.ordp + self.relprec: # or right._is_exact_zero()
                return self
            elif right.ordp <= self.ordp:
                ans = self._new_c(0)
                ans.ordp = right.ordp
            else:
                ans = self._new_c(right.ordp - self.ordp)
                ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                if self.relprec < 0:
                    ans.relprec = -ans.relprec
                ans.ordp = self.ordp
            return ans
        cdef long srprec = self.relprec
        if srprec < 0:
            srprec = -srprec
        cdef long rrprec = right.relprec
        if rrprec < 0:
            rrprec = -rrprec
        if self.ordp == right.ordp:
            # The relative precision of the sum is the minimum of the relative precisions in this case, possibly decreasing if we got cancellation
            # Since the valuations are the same, we could just add the units, if they had the same modulus.
            # But they don't necessarily, so we may have to conv_modulus
            if srprec == rrprec:
                ans = self._new_c(-srprec) # -srprec indicates that ans is not normalized
                self.prime_pow.restore_context_capdiv(srprec)
                ZZ_pX_add(ans.unit, self.unit, right.unit)
            elif srprec < rrprec:
                ans = self._new_c(-srprec)
                ZZ_pX_conv_modulus(ans.unit, right.unit, self.prime_pow.get_context_capdiv(srprec).x)
                self.prime_pow.restore_context_capdiv(srprec)
                # conv_modulus should have restored the context, so we don't need to again.
                ZZ_pX_add(ans.unit, ans.unit, self.unit)
            else:
                ans = self._new_c(-rrprec)
                ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(rrprec).x)
                self.prime_pow.restore_context_capdiv(rrprec)
                # conv_modulus should have restored the context, so we don't need to again.
                ZZ_pX_add(ans.unit, ans.unit, right.unit)
            ans.ordp = self.ordp
        elif self.ordp < right.ordp:
            tmpL = right.ordp - self.ordp
            if tmpL >= srprec:
                return self
            if srprec <= tmpL + rrprec:
                ans = self._new_c(-srprec)
            else:
                ans = self._new_c(-tmpL - rrprec)
            ans.ordp = self.ordp
            ZZ_pX_conv_modulus(ans.unit, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ans._internal_lshift(tmpL)
            if srprec <= tmpL + rrprec:
                ZZ_pX_add(ans.unit, ans.unit, self.unit)
            else:
                ZZ_pX_conv_modulus(tmpP, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                ZZ_pX_add(ans.unit, ans.unit, tmpP)
            # if self is normalized, then the valuations are actually different so the sum will be normalized.
            if self.relprec > 0:
                ans.relprec = -ans.relprec
        else:
            tmpL = self.ordp - right.ordp
            if tmpL >= rrprec:
                return right
            if rrprec <= tmpL + srprec:
                ans = self._new_c(-rrprec)
            else:
                ans = self._new_c(-tmpL - srprec)
            ans.ordp = right.ordp
            ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ans._internal_lshift(tmpL)
            if rrprec <= tmpL + srprec:
                ZZ_pX_add(ans.unit, ans.unit, right.unit)
            else:
                ZZ_pX_conv_modulus(tmpP, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
                ZZ_pX_add(ans.unit, ans.unit, tmpP)
            # if right is normalized, then the valuations are actually different so the sum will be normalized.
            if right.relprec > 0:
                ans.relprec = -ans.relprec
        return ans

    cpdef _sub_(self, right):
        """
        Return the difference of two elements.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(329)
            sage: b = W(111)
            sage: a - b  # indirect doctest
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
        # For now, a simple implementation
        return self + (-right)

    cpdef _mul_(self, _right):
        """
        Return the product of two elements.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(329)
            sage: b = W(111)
            sage: a*b  # indirect doctest
            4 + 3*w^5 + w^7 + 2*w^9 + 4*w^11 + 3*w^12 + 2*w^13 + w^14 + 2*w^15 + 3*w^16 + 4*w^17 + 4*w^18 + 2*w^19 + 2*w^21 + 4*w^22 + 2*w^23 + w^24 + O(w^25)
            sage: a * 0
            0
            sage: a * O(w^14)
            O(w^14)
        """
        cdef pAdicZZpXCRElement right = <pAdicZZpXCRElement>_right
        cdef ZZ_pX_c modulus_corrected
        cdef pAdicZZpXCRElement ans
        if self._is_exact_zero():
            return self
        if right._is_exact_zero():
            return right
        self._normalize()
        right._normalize()
        if self.relprec <= right.relprec:
            ans = self._new_c(self.relprec)
        else:
            ans = self._new_c(right.relprec)
        ans.ordp = self.ordp + right.ordp
        check_ordp(ans.ordp)
        if ans.relprec == 0:
            return ans
        if self.relprec == right.relprec:
            self.prime_pow.restore_context_capdiv(ans.relprec)
            sig_on()
            ZZ_pX_MulMod_pre(ans.unit, self.unit, right.unit, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
            sig_off()
        elif self.relprec < right.relprec:
            sig_on()
            ZZ_pX_conv_modulus(modulus_corrected, right.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ZZ_pX_MulMod_pre(ans.unit, self.unit, modulus_corrected, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
            sig_off()
        else:
            sig_on()
            ZZ_pX_conv_modulus(modulus_corrected, self.unit, self.prime_pow.get_context_capdiv(ans.relprec).x)
            ZZ_pX_MulMod_pre(ans.unit, right.unit, modulus_corrected, self.prime_pow.get_modulus_capdiv(ans.relprec)[0])
            sig_off()
        return ans

    cpdef _div_(self, right):
        """
        Return the quotient of two elements.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(14) / W(125)  # indirect doctest
            4*w^-15 + w^-13 + 3*w^-11 + 2*w^-10 + 3*w^-9 + 4*w^-8 + 4*w^-7 + 3*w^-6 + 2*w^-5 + 4*w^-4 + 3*w^-3 + 2*w^-2 + 4*w^-1 + 2 + w^2 + w^4 + 4*w^5 + w^6 + w^7 + 3*w^9 + O(w^10)
            sage: 1 / w
            w^-1 + O(w^24)
            sage: W.<w> = R.ext(x^25 - 165*x + 5)
            sage: a = (1 + w)^25 - 1
            sage: b = (1 + w)^5 - 1
            sage: c = (1 + w)^20 + (1 + w)^15 + (1 + w)^10 + (1 + w)^5 + 1
            sage: d = a / b; d == c
            True
            sage: d.precision_absolute()
            120
            sage: c.precision_absolute()
            125
            sage: 1 / a == ~a
            True
        """
        # for now, a simple implementation
        return self * (~right)

    def __copy__(self):
        """
        Return a copy of this element.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: b = W(45, 17); b
            4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16 + O(w^17)
            sage: c = copy(b); c
            4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16 + O(w^17)
            sage: c is b
            False
        """
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        ans.ordp = self.ordp
        ans.unit = self.unit
        return ans

    def _integer_(self, Z=None):
        r"""
        Return an integer congruent to this element modulo
        `\pi^{\mathrm{self.absolute_precision()}`, if possible.

        EXAMPLES::

            sage: ZZ(ZqCR(125,names='a')(-1))  # indirect doctest
            95367431640624
            sage: R = Zp(5); S.<x> = ZZ[]; f = x^5 + 25*x^3 - 5; W.<w> = R.ext(f)
            sage: ZZ(W(-1))
            95367431640624
            sage: ZZ(W(0))
            0
            sage: ZZ(W(0,7))
            0
            sage: ZZ(w)
            Traceback (most recent call last):
            ...
            ValueError: this element not well approximated by an integer
            sage: ZZ(W(5))  # todo: this should be different...
            381469726562505
        """
        cdef Integer ans
        cdef ZZ_c tmp_z
        if self._is_exact_zero() or self.relprec == 0:
            ans = PY_NEW(Integer)
            return ans
        if self.ordp < 0:
            self._normalize()
            if self.ordp < 0:
                raise ValueError("this element has negative valuation")
        cdef ntl_ZZ_pX f = <ntl_ZZ_pX>self._ntl_rep_abs()[0]
        if f.degree() > 0:
            raise ValueError("this element not well approximated by an integer")
        ans = PY_NEW(Integer)
        tmp_z = ZZ_p_rep(ZZ_pX_ConstTerm(f.x))
        ZZ_to_mpz(ans.value, &tmp_z)
        return ans

    def is_zero(self, absprec=None):
        r"""
        Return whether the valuation of this element is at least
        ``absprec``.  If ``absprec`` is ``None``, checks if this element
        is indistinguishable from zero.

        If this element is an inexact zero of valuation less than ``absprec``,
        raises a :exc:`PrecisionError`.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
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
        self._normalize()
        if self._is_exact_zero():
            ans = True
        elif absprec is None:
            ans = (self.relprec == 0)
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0:
                    ans = True
                elif self.relprec == 0:
                    raise PrecisionError("not enough precision to determine if element is zero")
                else:
                    ans = False
            else:
                aprec = mpz_get_si((<Integer>absprec).value)
                if self.relprec == 0 and aprec > self.ordp:
                    raise PrecisionError("not enough precision to determine if element is zero")
                else:
                    ans = (self.ordp >= aprec)
        return ans

    cpdef ntl_ZZ_pX _ntl_rep_unnormalized(self):
        """
        Return an ``ntl_ZZ_pX`` holding the current unit part of this element.

        The element is not normalized before this operation, so the polynomial
        returned may not actually be a unit.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566); b = W(209)
            sage: c = a + b; c._ntl_rep_unnormalized()
            [775]
            sage: c
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24 + O(w^25)
            sage: c._ntl_rep_unnormalized()
            [106 60 114 35 112]
        """
        if self.relprec == 0:
            raise ValueError("self == 0")
        self.prime_pow.restore_context_capdiv(self.relprec)
        cdef ntl_ZZ_pX ans = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        ans.c = self.prime_pow.get_context_capdiv(self.relprec)
        ans.x = self.unit
        return ans

    cpdef ntl_ZZ_pX _ntl_rep(self):
        """
        Return an ``ntl_ZZ_pX`` that holds the unit part of this element.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566); b = W(209)
            sage: c = a + b; c._ntl_rep()
            [106 60 114 35 112]
            sage: c
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24 + O(w^25)
            sage: c._ntl_rep()
            [106 60 114 35 112]
        """
        self._normalize()
        return self._ntl_rep_unnormalized()

    cpdef _ntl_rep_abs(self):
        """
        Return a pair ``(f, k)`` where ``f`` is an ``ntl_ZZ_pX`` and ``k`` is a
        nonpositive integer such that ``self = f(self.parent.gen())*p^k``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566); b = W(209)
            sage: c = a + b; c._ntl_rep_abs()
            ([775], 0)
            sage: c
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24 + O(w^25)
            sage: c._ntl_rep_abs()
            ([775], 0)
            sage: (~c)._ntl_rep_abs()
            ([121], -2)
            sage: ~c
            w^-10 + w^-8 + 4*w^-6 + 4*w^-5 + 3*w^-3 + 4*w^-2 + 3*w^-1 + 4 + 4*w + 2*w^2 + 4*w^3 + 3*w^4 + O(w^5)
            sage: ~c * 25
            1 + 4*w^5 + 3*w^7 + w^9 + 4*w^10 + 2*w^11 + 3*w^12 + w^13 + 4*w^14 + O(w^15)
            sage: W(121)
            1 + 4*w^5 + 3*w^7 + w^9 + 4*w^10 + 2*w^11 + 3*w^12 + w^13 + 4*w^14 + 2*w^16 + 3*w^17 + 3*w^18 + 4*w^19 + 4*w^20 + 3*w^21 + w^22 + w^23 + 4*w^24 + O(w^25)
        """
        self._normalize()
        if self.ordp == 0:
            return self._ntl_rep(), Integer(0)
        cdef ntl_ZZ_pContext_class ctx
        cdef long little_shift, ppow
        if self.ordp > 0:
            ctx = self.prime_pow.get_context_capdiv(self.ordp + self.relprec)
        else:
            little_shift = ((-self.ordp) % self.prime_pow.e)
            if little_shift != 0:
                little_shift = self.prime_pow.e - little_shift
            ctx = self.prime_pow.get_context_capdiv(self.relprec + little_shift)
        ctx.restore_c()
        cdef pAdicZZpXCRElement dummy = pAdicZZpXCRElement.__new__(pAdicZZpXCRElement)
        cdef ntl_ZZ_pX ans = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        cdef Integer ans_k = PY_NEW(Integer)
        dummy.unit = self.unit
        dummy.prime_pow = self.prime_pow
        if self.ordp > 0:
            dummy.relprec = self.ordp + self.relprec
            dummy._internal_lshift(self.ordp)
            ans.x = dummy.unit
        else:
            ppow = (self.ordp - little_shift) / self.prime_pow.e
            mpz_set_si(ans_k.value, ppow)
            dummy.ordp = 0 # _pshift_self wants ordp set
            dummy.relprec = self.relprec + little_shift
            # self = x^(self.prime_pow.e * ppow) * x^(little_shift) * self.unit
            # so we want to _internal_lshift dummy.unit by little_shift
            dummy._internal_lshift(little_shift)
            # and then write
            # self = p^(ppow) * (x^e/p)^(ppow) * dummy.unit
            # so we need to multiply dummy.unit by (p/x^e)^(-ppow) in the Eisenstein case
            # which we can do by _pshift_self
            dummy._pshift_self(-ppow)
            ans.x = dummy.unit
        ans.c = ctx
        return ans, ans_k

    def _polynomial_list(self, pad=False):
        """
        Return the coefficient list for a polynomial over the base ring
        yielding this element.

        INPUT:

        - ``pad`` -- whether to pad the result with zeros of the appropriate precision

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: W.<w> = Qp(5).extension(x^3 - 5)
            sage: (1 + w + O(w^11))._polynomial_list()
            [1 + O(5^4), 1 + O(5^4)]
            sage: (1 + w + O(w^11))._polynomial_list(pad=True)
            [1 + O(5^4), 1 + O(5^4), O(5^3)]
            sage: W(0)._polynomial_list()
            []
            sage: W(0)._polynomial_list(pad=True)
            [0, 0, 0]
            sage: W(O(w^7))._polynomial_list()
            []
            sage: W(O(w^7))._polynomial_list(pad=True)
            [O(5^3), O(5^2), O(5^2)]
            sage: T.<a> = Qp(5).extension(x^2 - 5)
            sage: T(1/5)._polynomial_list()
            [5^-1 + O(5^19)]
            sage: T(a^-800)._polynomial_list()
            [5^-400 + O(5^-380)]
            sage: T(O(a^-2))._polynomial_list()
            []
        """
        R = self.base_ring()
        if self.is_zero():
            L = []
            k = 0
        else:
            f, k = self._ntl_rep_abs()
            L = [Integer(c) for c in f.list()]
        if pad:
            n = self.parent().degree()
            L.extend([R.zero()] * (n - len(L)))
        if self._is_exact_zero():
            return L
        prec = self.relprec + self.ordp
        e = self.parent().e()
        if e == 1:
            return [R(c, prec-k) << k for c in L]
        else:
            return [R(c, (((prec - i - 1) // e) + 1) - k) << k for i, c in enumerate(L)]

    def polynomial(self, var='x'):
        """
        Return a polynomial over the base ring that yields this element
        when evaluated at the generator of the parent.

        INPUT:

        - ``var`` -- string, the variable name for the polynomial

        EXAMPLES::

            sage: S.<x> = ZZ[]
            sage: W.<w> = Zp(5).extension(x^2 - 5)
            sage: (w + W(5, 7)).polynomial()
            (1 + O(5^3))*x + 5 + O(5^4)
        """
        R = self.base_ring()
        S = R[var]
        return S(self._polynomial_list())

    cdef ZZ_p_c _const_term(self) noexcept:
        """
        Return the constant term of ``self.unit``.

        Note: this may be divisible by `p` if ``self`` is not normalized.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566)
            sage: a._const_term_test()  # indirect doctest
            566
        """
        return ZZ_pX_ConstTerm((<pAdicZZpXCRElement>self).unit)

    def is_equal_to(self, right, absprec=None):
        """
        Return whether this element is equal to ``right`` modulo ``self.uniformizer()^absprec``.

        If ``absprec`` is ``None``, checks whether this element is equal to ``right``
        modulo the lower of their two precisions.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(47); b = W(47 + 25)
            sage: a.is_equal_to(b)
            False
            sage: a.is_equal_to(b, 7)
            True
        """
        # Should be sped up later
        return (self - right).is_zero(absprec)

#    def lift(self):
#        """
#        Return an element of a number field defined by the same polynomial as
#        ``self``'s parent that is congruent to self modulo an appropriate ideal.

#        Not currently implemented.
#        """
#        raise NotImplementedError

    cpdef pAdicZZpXCRElement lift_to_precision(self, absprec=None):
        """
        Return a ``pAdicZZpXCRElement`` congruent to this element but with
        absolute precision at least ``absprec``.

        INPUT:

        - ``absprec`` -- (default: ``None``) the absolute precision of
          the result.  If ``None``, lifts to the maximum precision
          allowed.

        .. NOTE::

            If setting ``absprec`` that high would violate the
            precision cap, raises a precision error.  If ``self`` is an
            inexact zero and ``absprec`` is greater than the maximum
            allowed valuation, raises an error.

            Note that the new digits will not necessarily be zero.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(345, 17); a
            4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15 + O(w^17)
            sage: b = a.lift_to_precision(19); b
            4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15
             + w^17 + 2*w^18 + O(w^19)
            sage: c = a.lift_to_precision(24); c
            4*w^5 + 3*w^7 + w^9 + 3*w^10 + 2*w^11 + 4*w^12 + w^13 + 2*w^14 + 2*w^15
             + w^17 + 2*w^18 + 4*w^19 + 4*w^20 + 2*w^21 + 4*w^23 + O(w^24)
            sage: a._ntl_rep()
            [19 35 118 60 121]
            sage: b._ntl_rep()
            [19 35 118 60 121]
            sage: c._ntl_rep()
            [19 35 118 60 121]
            sage: a.lift_to_precision().precision_relative() == W.precision_cap()
            True
        """
        cdef pAdicZZpXCRElement ans
        cdef long aprec, rprec
        self._normalize()
        if self._is_exact_zero():
            return self
        if absprec is not None and not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if absprec is None:
            if self.relprec == 0:
                # return an exact zero
                ans = self._new_c(0)
                ans._set_exact_zero()
                return ans
            aprec = self.prime_pow.ram_prec_cap + self.ordp
        elif mpz_fits_slong_p((<Integer>absprec).value) == 0:
            if mpz_sgn((<Integer>absprec).value) < 0 or self.relprec == self.prime_pow.ram_prec_cap:
                return self
            else:
                if self.relprec == 0:
                    raise ValueError("absprec larger than maximum allowable valuation")
                else:
                    raise PrecisionError("precision higher than allowed by the precision cap")
        else:
            aprec = mpz_get_si((<Integer>absprec).value)
        if aprec <= self.ordp + self.relprec:
            return self
        if self.relprec == 0:
            if self.ordp >= aprec:
                return self
            elif aprec >= maxordp:
                raise ValueError("absprec larger than maximum allowable valuation")
            else:
                ans = self._new_c(0)
                ans._set_inexact_zero(aprec)
                return ans
        # Now we're done handling all the special cases.
        rprec = aprec - self.ordp
        if rprec > self.prime_pow.ram_prec_cap:
            raise PrecisionError("precision higher than allowed by the precision cap")
        ans = self._new_c(rprec)
        ans.ordp = self.ordp
        ZZ_pX_conv_modulus(ans.unit, self.unit, self.prime_pow.get_context_capdiv(rprec).x)
        return ans

    def expansion(self, n=None, lift_mode='simple'):
        """
        Return a list giving a series representation of ``self``.

        - If ``lift_mode == 'simple'`` or ``'smallest'``, the returned
          list will consist of integers (in the Eisenstein case) or a
          list of lists of integers (in the unramified case).  ``self``
          can be reconstructed as a sum of elements of the list times
          powers of the uniformiser (in the Eisenstein case), or as a
          sum of powers of the `p` times polynomials in the generator
          (in the unramified case).

          + If ``lift_mode == 'simple'``, all integers will be in the interval
            `[0,p-1]`.

          + If ``lift_mode == 'smallest'`` they will be in the
            interval `[(1-p)/2, p/2]`.

        - If ``lift_mode == 'teichmuller'``, returns a list of
          ``pAdicZZpXCRElements``, all of which are Teichmuller
          representatives and such that ``self`` is the sum of that list
          times powers of the uniformizer.

        Note that zeros are truncated from the returned list if
        ``self.parent()`` is a field, so you must use the
        ``valuation`` function to fully reconstruct ``self``.

        INPUT:

        - ``n`` -- integer (default: ``None``); if given, returns the
          corresponding entry in the expansion

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W(775, 19); y
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: (y>>9).expansion()
            [0, 1, 0, 4, 0, 2, 1, 2, 4, 1]
            sage: (y>>9).expansion(lift_mode='smallest')
            [0, 1, 0, -1, 0, 2, 1, 2, 0, 1]
            sage: w^10 - w^12 + 2*w^14 + w^15 + 2*w^16 + w^18 + O(w^19)
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: y = 75 + 45*a + 1200*a^2; y
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
            sage: E = y.expansion(); E
            5-adic expansion of 4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
            sage: list(E)
            [[], [0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1], []]
            sage: list(y.expansion(lift_mode='smallest'))
            [[], [0, -1], [-2, 2, -2], [1], [0, 0, 2], []]
            sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
            sage: list(W(0).expansion())
            []
            sage: list(W(0,4).expansion())
            []
            sage: list(A(0,4).expansion())
            []

        TESTS:

        We check that :issue:`24949` is fixed::

            sage: R = Zp(2)
            sage: S.<x> = R[]
            sage: A.<a> = R.extension(x^10 + 2)
            sage: u = a^4 + a^5
            sage: v = a^2 + a^3
            sage: w = u - v^2
            sage: w.expansion(4)
            0
        """
        self._normalize()
        if lift_mode == 'teichmuller':
            zero = self.parent()(0)
        elif self.prime_pow.e == 1:
            zero = []
        else:
            zero = Integer(0)
        if isinstance(n, slice):
            return self.slice(n.start, n.stop, n.step)
        elif n is not None:
            if self._is_exact_zero() or n < self.ordp:
                return zero
            elif n >= self.ordp + self.relprec:
                raise PrecisionError
        if self.relprec == 0: # cannot have n = None
            return []
        if lift_mode == 'simple':
            ulist = self.ext_p_list(pos=True)
        elif lift_mode == 'smallest':
            ulist = self.ext_p_list(pos=False)
        elif lift_mode == 'teichmuller':
            if n is None:
                ulist = self.teichmuller_expansion()
            else:
                return self.teichmuller_expansion(n)
        else:
            raise ValueError("lift mode must be one of 'simple', 'smallest' or 'teichmuller'")
        if n is not None:
            try:
                return ulist[n - self.ordp]
            except IndexError:
                return zero
        if self.prime_pow.in_field == 0 and self.ordp > 0:
            ulist = [zero] * self.ordp + ulist
        return ulist

    def matrix_mod_pn(self):
        r"""
        Return the matrix of right multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for this
        extension field.  Thus the *rows* of this matrix give the
        images of each of the `x^i`.  The entries of the matrices are
        :class:`IntegerMod` elements, defined modulo `p^{N / e}` where `N` is
        the absolute precision of this element (unless this element is
        zero to arbitrary precision; in that case the entries are
        integer zeros.)

        Raises an error if this element has negative valuation.

        EXAMPLES::

            sage: R = ZpCR(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = (3+w)^7
            sage: a.matrix_mod_pn()
            [2757  333 1068  725 2510]
            [  50 1507  483  318  725]
            [ 500   50 3007 2358  318]
            [1590 1375 1695 1032 2358]
            [2415  590 2370 2970 1032]

        TESTS:

        Check that :issue:`13617` has been fixed::

            sage: W.zero().matrix_mod_pn()
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
            [0 0 0 0 0]
        """
        if self.valuation_c() < 0:
            raise ValueError("self must be integral")
        n = self.prime_pow.deg
        from sage.matrix.constructor import matrix
        if self._is_exact_zero():
            from sage.rings.integer_ring import IntegerRing
            return matrix(IntegerRing(), n, n)
        R = IntegerModRing(self.prime_pow.pow_Integer(self.prime_pow.capdiv(self.ordp + self.relprec)))
        L = []
        cdef ntl_ZZ_pX cur = <ntl_ZZ_pX>self._ntl_rep_abs()[0]
        cur.c.restore_c()
        cdef ZZ_pX_Modulus_c* m = self.prime_pow.get_modulus_capdiv(self.ordp + self.relprec)
        cdef ZZ_pX_c x
        ZZ_pX_SetX(x)
        cdef Py_ssize_t i
        zero = int(0)
        for i in range(n):
            curlist = cur.list()
            L.extend(curlist + [zero] * (n - len(curlist)))
            ZZ_pX_MulMod_pre(cur.x, cur.x, x, m[0])
        return matrix(R, n, n,  L)

#     def matrix(self, base=None):
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

#         INPUT:
#             base -- field or morphism
#         """
#         raise NotImplementedError

#     def multiplicative_order(self, prec=None):
#         """
#         Return the multiplicative order of ``self``, ie the smallest
#         positive `n` so that there is an exact `p`-adic element congruent
#         to ``self`` modulo ``self``'s precision that is an `n`-th root of unity.

#         Note: unlike the case for Qp and Zp, it is possible to have
#         non-teichmuller elements with finite orders.  This can happen
#         only if (p-1) divides the ramification index (see the
#         documentation on __pow__).

#         INPUT:

#         - self -- a `p`-adic element
#         - ``prec`` -- integer

#         OUTPUT: integer; the multiplicative order of ``self``
#         """
#         raise NotImplementedError

    def teichmuller_expansion(self, n=None):
        r"""
        Return a list [`a_0`, `a_1`,..., `a_n`] such that:

        - `a_i^q = a_i`
        - ``self.unit_part()`` = `\sum_{i = 0}^n a_i \pi^i`, where `\pi` is a
          uniformizer of ``self.parent()``
        - if `a_i \ne 0`, the absolute precision of `a_i` is
          ``self.precision_relative() - i``

        INPUT:

        - ``n`` -- integer (default: ``None``); if given, returns the
          corresponding entry in the expansion

        EXAMPLES::

            sage: R.<a> = ZqCR(5^4,4)
            sage: E = a.teichmuller_expansion(); E
            5-adic expansion of a + O(5^4) (teichmuller)
            sage: list(E)
            [a + (2*a^3 + 2*a^2 + 3*a + 4)*5 + (4*a^3 + 3*a^2 + 3*a + 2)*5^2
               + (4*a^2 + 2*a + 2)*5^3 + O(5^4),
             (3*a^3 + 3*a^2 + 2*a + 1) + (a^3 + 4*a^2 + 1)*5 + (a^2 + 4*a + 4)*5^2 + O(5^3),
             (4*a^3 + 2*a^2 + a + 1) + (2*a^3 + 2*a^2 + 2*a + 4)*5 + O(5^2),
             (a^3 + a^2 + a + 4) + O(5)]
            sage: sum([c * 5^i for i, c in enumerate(E)])
            a + O(5^4)
            sage: all(c^625 == c for c in E)
            True

            sage: S.<x> = ZZ[]
            sage: f = x^3 - 98*x + 7
            sage: W.<w> = ZpCR(7,3).ext(f)
            sage: b = (1+w)^5; L = b.teichmuller_expansion(); L
            [1 + O(w^9), 5 + 5*w^3 + w^6 + 4*w^7 + O(w^8), 3 + 3*w^3 + O(w^7),
             3 + 3*w^3 + O(w^6), O(w^5), 4 + 5*w^3 + O(w^4), 3 + O(w^3),
             6 + O(w^2), 6 + O(w)]
            sage: sum([w^i*L[i] for i in range(9)]) == b
            True
            sage: all(L[i]^(7^3) == L[i] for i in range(9))
            True

            sage: L = W(3).teichmuller_expansion(); L
            [3 + 3*w^3 + w^7 + O(w^9), O(w^8), O(w^7), 4 + 5*w^3 + O(w^6),
             O(w^5), O(w^4), 3 + O(w^3), 6 + O(w^2)]
            sage: sum([w^i*L[i] for i in range(len(L))])
            3 + O(w^9)
        """
        cdef long rp = self.relprec
        cdef pAdicZZpXCRElement v
        if n is None:
            L = []
            if rp == 0:
                return L
        elif self._is_exact_zero() or n < self.ordp:
            return self.parent()(0)
        elif n >= self.ordp + rp:
            raise PrecisionError
        else:
            v = self._new_c(rp)
        cdef pAdicZZpXCRElement u = self.unit_part()
        cdef long goal
        if n is not None: goal = rp - n + self.ordp
        while u.relprec > 0:
            v = self._new_c(rp)
            self.prime_pow.teichmuller_set_c(&v.unit, &u.unit, rp)
            v.ordp = 0
            if n is None:
                L.append(v)
            elif rp == goal:
                return v
            if rp == 1: break
            ZZ_pX_sub(u.unit, u.unit, v.unit)
            u.relprec = -u.relprec
            u._normalize()
            if u.relprec == 0: break
            rp -= 1
            u.ordp -= 1
            while u.ordp > 0:
                if n is None:
                    v = self._new_c(0)
                    v._set_inexact_zero(rp)
                    L.append(v)
                elif rp == goal:
                    v = self._new_c(0)
                    v._set_inexact_zero(rp)
                    return v
                rp -= 1
                u.ordp -= 1
        if n is None:
            return L
        else:
            v = self._new_c(0)
            v._set_inexact_zero(rp)
            return v

    def _teichmuller_set_unsafe(self):
        """
        Set this element to the Teichmuller representative with the
        same residue.

        .. WARNING::

            This function modifies the element, which is not safe.
            Elements are supposed to be immutable.

        EXAMPLES::

            sage: R = Zp(7,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 77*x^3 - 98*x^2 - 7
            sage: W.<w> = R.ext(f)
            sage: y = W.teichmuller(3, 15); y  # indirect doctest
            3 + 4*w^5 + 2*w^8 + 6*w^10 + w^11 + 6*w^12 + 5*w^13 + 4*w^14 + O(w^15)

            sage: y^7 == y
            True
            sage: g = x^3 + 3*x^2 + 4
            sage: A.<a> = R.ext(g)
            sage: b = A.teichmuller(1 + 2*a - a^2); b
            (6*a^2 + 2*a + 1) + (5*a + 3)*7 + (5*a + 5)*7^2 + (4*a^2 + 4*a + 2)*7^3 + (2*a + 1)*7^4 + O(7^5)
            sage: b^343 == b
            True

        TESTS:

        We check that :issue:`8239` is resolved::

            sage: K.<a> = Qq(25)
            sage: K.teichmuller(K(2/5))
            Traceback (most recent call last):
            ...
            ValueError: cannot set negative valuation element to Teichmuller representative
        """
        self._normalize()
        if self.ordp > 0:
            self._set_exact_zero()
        elif self.ordp < 0:
            raise ValueError("cannot set negative valuation element to Teichmuller representative")
        elif self.relprec == 0:
            raise ValueError("not enough precision known")
        else:
            self.prime_pow.teichmuller_set_c(&self.unit, &self.unit, self.relprec)

#    def padded_list(self, n, lift_mode='simple'):
#        """
#        Return a list of coefficients of pi starting with `pi^0` up to
#        `pi^n` exclusive (padded with zeros if needed)

#        """
#        raise NotImplementedError

    def precision_absolute(self):
        """
        Return the absolute precision of this element, i.e., the power of the
        uniformizer modulo which this element is defined.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
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
            sage: (a.unit_part() - 3).precision_absolute()
            9
        """
        cdef Integer ans
        if self.ordp == maxordp:
            return infinity
        else:
            ans = PY_NEW(Integer)
            if self.relprec > 0:
                mpz_set_si(ans.value, self.relprec + self.ordp)
            else:
                mpz_set_si(ans.value, -self.relprec + self.ordp)
            return ans

    def precision_relative(self):
        """
        Return the relative precision of this element, i.e., the power of the
        uniformizer modulo which the unit part of ``self`` is defined.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
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
        self._normalize()
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.relprec)
        return ans

#    def residue(self, n=1):
#        """
#        Reduces this element modulo pi^n.
#        """
#        raise NotImplementedError

    cdef long valuation_c(self) noexcept:
        """
        Return the valuation of this element.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75, 19); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + O(w^19)
            sage: a.valuation()  # indirect doctest
            10
            sage: a.precision_absolute()
            19
            sage: a.precision_relative()
            9
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + O(w^9)
        """
        self._normalize()
        return self.ordp

    cpdef pAdicZZpXCRElement unit_part(self):
        """
        Return the unit part of this element, ie ``self / uniformizer^(self.valuation())``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
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

        TESTS:

        We check that :issue:`13616` is resolved::

            sage: z = (1+w)^5
            sage: y = z - 1
            sage: t = y - y
            sage: t.unit_part()
            O(w^0)
        """
        self._normalize()
        cdef pAdicZZpXCRElement ans = self._new_c(self.relprec)
        ans.ordp = 0
        if self.relprec != 0:
            ans.unit = self.unit
        return ans

    cdef ext_p_list(self, bint pos):
        """
        Return a list of integers (in the Eisenstein case) or a list
        of lists of integers (in the unramified case).  ``self`` can be
        reconstructed as a sum of elements of the list times powers of
        the uniformiser (in the Eisenstein case), or as a sum of
        powers of `p` times polynomials in the generator (in the
        unramified case).

        If ``pos`` is ``True``, all integers will be in the interval `[0,p-1]`,
        otherwise they will be in the interval `[(1-p)/2, p/2]`.

        Note that zeros are truncated from the returned list, so you
        must use the method :meth:`valuation` to completely recover ``self``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
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
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
            sage: y._ext_p_list(True)
            [[0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1]]
            sage: y._ext_p_list(False)
            [[0, -1], [-2, 2, -2], [1], [0, 0, 2]]
            sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4 + O(5^6)
        """
        self._normalize()
        return self.ext_p_list_precs(pos, self.relprec)


def make_ZZpXCRElement(parent, unit, ordp, relprec, version):
    """
    Unpickling.

    EXAMPLES::

        sage: R = Zp(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
        sage: W.<w> = R.ext(f)
        sage: y = W(775, 19); y
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
        sage: loads(dumps(y))  # indirect doctest
        w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)

        sage: from sage.rings.padics.padic_ZZ_pX_CR_element import make_ZZpXCRElement
        sage: make_ZZpXCRElement(W, y._ntl_rep(), 3, 9, 0)
        w^3 + 4*w^5 + 2*w^7 + w^8 + 2*w^9 + 4*w^10 + w^11 + O(w^12)
    """
    cdef pAdicZZpXCRElement ans
    cdef ZZ_pX_c poly
    if version == 0:
        ans = pAdicZZpXCRElement(parent, [], empty = True)
        if relprec == 0:
            ans._set_inexact_zero(mpz_get_si((<Integer>ordp).value))
        else:
            ans.prime_pow.restore_context_capdiv(mpz_get_si((<Integer>relprec).value))
            poly = (<ntl_ZZ_pX>unit).x
            ans._set(&poly, mpz_get_si((<Integer>ordp).value), mpz_get_si((<Integer>relprec).value))
        return ans
    else:
        raise ValueError("unknown unpickling version")
