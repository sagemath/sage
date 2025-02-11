# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
# sage.doctest: needs sage.rings.padics
r"""
`p`-adic ``ZZ_pX`` FM Element

This file implements elements of Eisenstein and unramified extensions
of `\ZZ_p` with fixed modulus precision.

For the parent class see ``padic_extension_leaves.pyx``.

The underlying implementation is through NTL's ``ZZ_pX`` class.  Each
element contains the following data:

- ``value`` (``ZZ_pX_c``) -- an ntl ``ZZ_pX`` storing the value.  The
  variable `x` is the uniformizer in the case of Eisenstein extensions.
  This ``ZZ_pX`` is created with global ntl modulus determined by the
  parent's precision cap and shared among all elements.

- ``prime_pow`` (some subclass of :class:`PowComputer_ZZ_pX`) -- a class,
  identical among all elements with the same parent, holding common
  data.

  * ``prime_pow.deg`` -- the degree of the extension

  * ``prime_pow.e`` -- the ramification index

  * ``prime_pow.f`` -- the inertia degree

  * ``prime_pow.prec_cap`` -- the unramified precision cap: for
    Eisenstein extensions this is the smallest power of `p` that is
    zero

  * ``prime_pow.ram_prec_cap`` -- the ramified precision cap: for
    Eisenstein extensions this will be the smallest power of `x` that
    is indistinguishable from zero

  * ``prime_pow.pow_ZZ_tmp``, prime_pow.pow_mpz_t_tmp``,
    ``prime_pow.pow_Integer`` -- functions for accessing powers of
    `p`.  The first two return pointers.  See
    ``sage/rings/padics/pow_computer_ext`` for examples and important
    warnings.

  * ``prime_pow.get_context``, ``prime_pow.get_context_capdiv``,
    ``prime_pow.get_top_context`` -- obtain an
    ``ntl_ZZ_pContext_class`` corresponding to `p^n`.  The capdiv
    version divides by ``prime_pow.e`` as appropriate.
    ``top_context`` corresponds to `p^{\texttt{prec\_cap}}`.

  * ``prime_pow.restore_context``,
    ``prime_pow.restore_context_capdiv``,
    ``prime_pow.restore_top_context`` -- restores the given context

  * ``prime_pow.get_modulus``, ``get_modulus_capdiv``,
    ``get_top_modulus`` -- returns a ``ZZ_pX_Modulus_c*`` pointing to
    a polynomial modulus defined modulo `p^n` (appropriately divided
    by ``prime_pow.e`` in the capdiv case).

EXAMPLES:

An Eisenstein extension::

    sage: R = ZpFM(5,5)
    sage: S.<x> = R[]
    sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
    sage: W.<w> = R.ext(f); W
    5-adic Eisenstein Extension Ring in w defined by x^5 + 75*x^3 - 15*x^2 + 125*x - 5
    sage: z = (1+w)^5; z
    1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24
    sage: y = z >> 1; y
    w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + 4*w^24
    sage: y.valuation()
    4
    sage: y.precision_relative()
    21
    sage: y.precision_absolute()
    25
    sage: z - (y << 1)
    1

An unramified extension::

    sage: # needs sage.libs.flint
    sage: g = x^3 + 3*x + 3
    sage: A.<a> = R.ext(g)
    sage: z = (1+a)^5; z
    (2*a^2 + 4*a) + (3*a^2 + 3*a + 1)*5 + (4*a^2 + 3*a + 4)*5^2 + (4*a^2 + 4*a + 4)*5^3 + (4*a^2 + 4*a + 4)*5^4
    sage: z - 1 - 5*a - 10*a^2 - 10*a^3 - 5*a^4 - a^5
    0
    sage: y = z >> 1; y
    (3*a^2 + 3*a + 1) + (4*a^2 + 3*a + 4)*5 + (4*a^2 + 4*a + 4)*5^2 + (4*a^2 + 4*a + 4)*5^3
    sage: 1/a
    (3*a^2 + 4) + (a^2 + 4)*5 + (3*a^2 + 4)*5^2 + (a^2 + 4)*5^3 + (3*a^2 + 4)*5^4

Different printing modes::

    sage: # needs sage.libs.flint
    sage: R = ZpFM(5, print_mode='digits'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; repr(z)
    '...4110403113210310442221311242000111011201102002023303214332011214403232013144001400444441030421100001'
    sage: R = ZpFM(5, print_mode='bars'); S.<x> = R[]; g = x^3 + 3*x + 3; A.<a> = R.ext(g)
    sage: z = (1+a)^5; repr(z)
    '...[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 4, 4]|[4, 3, 4]|[1, 3, 3]|[0, 4, 2]'
    sage: R = ZpFM(5, print_mode='terse'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: z = (1+w)^5; z
    6 + 95367431640505*w + 25*w^2 + 95367431640560*w^3 + 5*w^4
    sage: R = ZpFM(5, print_mode='val-unit'); S.<x> = R[]; f = x^5 + 75*x^3 - 15*x^2 + 125*x -5; W.<w> = R.ext(f)
    sage: y = (1+w)^5 - 1; y
    w^5 * (2090041 + 19073486126901*w + 1258902*w^2 + 57220458985049*w^3 + 16785*w^4)

AUTHORS:

- David Roe  (2008-01-01) initial version
"""
# ****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

include "sage/libs/ntl/decl.pxi"

from sage.structure.richcmp cimport rich_to_bool
from sage.rings.integer cimport Integer
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.padic_ext_element cimport pAdicExtElement
from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.rings.rational cimport Rational
from cypari2.gen cimport Gen as pari_gen
from sage.interfaces.abc import GpElement
from sage.rings.finite_rings.integer_mod import IntegerMod_abstract
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing


cdef class pAdicZZpXFMElement(pAdicZZpXElement):
    def __init__(self, parent, x, absprec=None, relprec=None, empty=False):
        r"""
        Create an element of a fixed modulus, unramified or
        eisenstein extension of `\ZZ_p` or `\QQ_p`.

        INPUT:

        - ``parent`` -- either an ``EisensteinRingFixedMod`` or
          ``UnramifiedRingFixedMod``

        - ``x`` -- integer, rational, `p`-adic element, polynomial,
          list, integer_mod, pari int/frac/poly_t/pol_mod, an
          ``ntl_ZZ_pX``, an ``ntl_ZZX``, an ``ntl_ZZ``, or an
          ``ntl_ZZ_p``

        - ``absprec`` -- not used

        - ``relprec`` -- not used

        - ``empty`` -- whether to return after initializing to zero
          (without setting anything)

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1+w)^5; z  # indirect doctest
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24

        TESTS:

        Check that :issue:`3865` is fixed::

            sage: W(gp('2 + O(5^2)'))
            2

        Check that :issue:`13612` has been fixed::

            sage: # needs sage.libs.flint
            sage: R = ZpFM(3)
            sage: S.<a> = R[]
            sage: W.<a> = R.extension(a^2 + 1)
            sage: W(W.residue_field().zero())
            0
        """
        pAdicZZpXElement.__init__(self, parent)
        if empty:
            return
        cdef ZZ_c tmp_z
        cdef Integer tmp_Int
        cdef Py_ssize_t i
        if isinstance(x, pAdicGenericElement):
            if x.valuation() < 0:
                raise ValueError("element has negative valuation")
            if x._is_base_elt(self.prime_pow.prime):
                xlift = <Integer>x.lift()
                self._set_from_mpz(xlift.value)
                return
            if parent.prime() != x.parent().prime():
                raise TypeError("Cannot coerce between p-adic parents with different primes.")
        if isinstance(x, GpElement):
            x = x.__pari__()
        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                if x.variable() != self.prime_pow.prime:
                    raise TypeError("Cannot coerce a pari p-adic with the wrong prime.")
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
            if (<Integer>x.modulus())._is_power_of(<Integer>parent.prime()):
                x = x.lift()
            else:
                raise TypeError("cannot coerce from the given integer mod ring (not a power of the same prime)")
        elif x in parent.residue_field() and x.parent().is_finite():
            # Should only reach here if x is not in F_p
            z = parent.gen()
            poly = x.polynomial().list()
            x = sum([poly[i].lift() * (z ** i) for i in range(len(poly))], parent.zero())
        elif isinstance(x, ntl_ZZ_p):
            ZZ_remove(tmp_z, (<ntl_ZZ>x.modulus()).x, self.prime_pow.pow_ZZ_tmp(1)[0])
            if ZZ_IsOne(tmp_z):
                x = x.lift()
                tmp_Int = Integer.__new__(Integer)
                ZZ_to_mpz(tmp_Int.value, &(<ntl_ZZ>x).x)
                x = tmp_Int
            else:
                raise TypeError("cannot coerce the given ntl_ZZ_p (modulus not a power of the same prime)")
        elif isinstance(x, ntl_ZZ):
            tmp_Int = Integer.__new__(Integer)
            ZZ_to_mpz(tmp_Int.value, &(<ntl_ZZ>x).x)
            x = tmp_Int
        elif isinstance(x, int):
            x = Integer(x)
        if isinstance(x, Integer):
            self._set_from_mpz((<Integer>x).value)
        elif isinstance(x, Rational):
            self._set_from_mpq((<Rational>x).value)
        elif isinstance(x, ntl_ZZ_pX):
            self._set_from_ZZ_pX(&(<ntl_ZZ_pX>x).x, (<ntl_ZZ_pX>x).c)
        elif isinstance(x, ntl_ZZX):
            self._set_from_ZZX((<ntl_ZZX>x).x)
        elif isinstance(x, pAdicExtElement):
            if x.parent() is parent:
                self._set_from_ZZ_pX(&(<pAdicZZpXFMElement>x).value, self.prime_pow.get_top_context())
            else:
                raise NotImplementedError("Conversion from different p-adic extensions not yet supported")
        else:
            try:
                x = list(x)
            except TypeError:
                try:
                    x = x.list()
                except AttributeError:
                    raise TypeError("cannot convert x to a p-adic element")
            self._set_from_list(x)

    cdef int _set_from_mpz(self, mpz_t x) except -1:
        """
        Set ``self`` from an ``mpz_t``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(70)  # indirect doctest
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23
        """
        self.prime_pow.restore_top_context()
        cdef ZZ_c tmp
        cdef mpz_t tmp_m
        sig_on()
        mpz_init(tmp_m)
        mpz_set(tmp_m, x)
        mpz_to_ZZ(&tmp, tmp_m)
        mpz_clear(tmp_m)
        ZZ_pX_SetCoeff(self.value, 0, ZZ_to_ZZ_p(tmp))
        sig_off()

    cdef int _set_from_mpq(self, mpq_t x) except -1:
        """
        Set ``self`` from an ``mpq_t``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(70/3); z  # indirect doctest
            3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + 3*w^15 + 2*w^16
             + 3*w^17 + w^18 + 3*w^19 + 3*w^20 + 2*w^21 + 2*w^22 + 3*w^23 + 4*w^24
            sage: z * 3
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23
            sage: W(70)
            4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11 + w^13 + 3*w^16 + w^17 + w^18 + 4*w^20 + 4*w^21 + w^22 + 2*w^23
        """
        self.prime_pow.restore_top_context()
        if mpz_divisible_p(mpq_denref(x), self.prime_pow.prime.value):
            raise ValueError("p divides denominator")
        cdef mpz_t tmp_m
        cdef ZZ_c tmp_z
        sig_on()
        mpz_init(tmp_m)
        mpz_invert(tmp_m, mpq_denref(x), self.prime_pow.pow_mpz_t_top())
        mpz_mul(tmp_m, tmp_m, mpq_numref(x))
        mpz_mod(tmp_m, tmp_m, self.prime_pow.pow_mpz_t_top())
        mpz_to_ZZ(&tmp_z, tmp_m)
        ZZ_pX_SetCoeff(self.value, 0, ZZ_to_ZZ_p(tmp_z))
        mpz_clear(tmp_m)
        sig_off()
        return 0

    cdef int _set_from_ZZ_pX(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx) except -1:
        """
        Set ``self`` from a ``ZZ_pX_c``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2)); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + 4*w^13 + 3*w^14 + 2*w^15 + w^16
             + 3*w^18 + 2*w^19 + 4*w^20 + 4*w^21 + 2*w^22 + 2*w^23 + 4*w^24
            sage: z._ntl_rep()
            [4 1 16]
        """
        self.prime_pow.restore_top_context()
        self._check_ZZ_pContext(ctx)
        ZZ_pX_conv_modulus(self.value, poly[0], self.prime_pow.get_top_context().x)

    cdef int _set_from_ZZX(self, ZZX_c poly) except -1:
        """
        Set ``self`` from a ``ZZX`` with relative precision bounded
        by ``relprec`` and absolute precision bounded by ``absprec``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZX([4,1,16])); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + 2*w^11 + 4*w^13 + 3*w^14 + 2*w^15 + w^16
             + 3*w^18 + 2*w^19 + 4*w^20 + 4*w^21 + 2*w^22 + 2*w^23 + 4*w^24
            sage: z._ntl_rep()
            [4 1 16]
        """
        self.prime_pow.restore_top_context()
        ZZX_to_ZZ_pX(self.value, poly)

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Test if ``self`` is an inexact zero.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(0)
            sage: z._is_inexact_zero()
            True
            sage: z = W(5^6)
            sage: z._is_inexact_zero()
            True
        """
        return ZZ_pX_IsZero(self.value) or (self.prime_pow.e * self.prime_pow.prec_cap != self.prime_pow.ram_prec_cap
                                            and self.valuation_c() >= self.prime_pow.ram_prec_cap)

    def __reduce__(self):
        """
        Pickle ``self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5 - 1
            sage: loads(dumps(z)) == z  # indirect doctest
            True
        """
        self.prime_pow.restore_top_context()
        cdef ntl_ZZ_pX holder = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        holder.c = self.prime_pow.get_top_context()
        holder.x = self.value
        return make_ZZpXFMElement, (self.parent(), holder)

    cdef pAdicZZpXFMElement _new_c(self):
        """
        Return a new element with the same parent as ``self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w^5 + 1  # indirect doctest
            1 + w^5
        """
        self.prime_pow.restore_top_context()
        cdef pAdicZZpXFMElement ans = pAdicZZpXFMElement.__new__(pAdicZZpXFMElement)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        return ans

    cpdef _richcmp_(left, right, int op):
        """
        First compare valuations, then compare the values.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w == 1  # indirect doctest
            False
            sage: y = 1 + w
            sage: z = 1 + w + w^27
            sage: y == z
            True
        """
        cdef pAdicZZpXFMElement _left = <pAdicZZpXFMElement>left
        cdef pAdicZZpXFMElement _right = <pAdicZZpXFMElement>right
        cdef long x_ordp = _left.valuation_c()
        cdef long y_ordp = _right.valuation_c()
        if x_ordp < y_ordp:
            return rich_to_bool(op, -1)
        elif x_ordp > y_ordp:
            return rich_to_bool(op, 1)
        else:  # equal ordp
            _left.prime_pow.restore_top_context()
            if x_ordp == left.prime_pow.ram_prec_cap:
                return rich_to_bool(op, 0)  # since both are zero
            elif _left.value == _right.value:
                return rich_to_bool(op, 0)
            else:
                # for now just return 1
                return rich_to_bool(op, 1)

    def __invert__(self):
        """
        Return the inverse of ``self``, as long as ``self`` is a unit.

        If ``self`` is not a unit, raises a :exc:`ValueError`.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: y = ~z; y  # indirect doctest
            1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14 + 3*w^16
             + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24
            sage: y.parent()
            5-adic Eisenstein Extension Ring in w defined by x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: z = z - 1
            sage: ~z
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        if self.valuation_c() > 0:
            raise ValueError("cannot invert non-unit")
        cdef pAdicZZpXFMElement ans = self._new_c()
        sig_on()
        if self.prime_pow.e == 1:
            ZZ_pX_InvMod_newton_unram(ans.value, self.value, self.prime_pow.get_top_modulus()[0], self.prime_pow.get_top_context().x, self.prime_pow.get_context(1).x)
        else:
            ZZ_pX_InvMod_newton_ram(ans.value, self.value, self.prime_pow.get_top_modulus()[0], self.prime_pow.get_top_context().x)
        sig_off()
        return ans

    cdef pAdicZZpXFMElement _lshift_c(self, long n):
        """
        Multiply ``self`` by the uniformizer raised to the power ``n``.

        If ``n`` is negative, right shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24
            sage: z << 17  # indirect doctest
            w^17 + w^22 + w^23 + 2*w^24
            sage: z << (-1)
            w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + 4*w^24
        """
        if n < 0:
            return self._rshift_c(-n)
        elif n == 0:
            return self
        cdef pAdicZZpXFMElement ans = self._new_c()
        if n < self.prime_pow.ram_prec_cap:
            if self.prime_pow.e == 1:
                ZZ_pX_left_pshift(ans.value, self.value, self.prime_pow.pow_ZZ_tmp(n)[0], self.prime_pow.get_top_context().x)
            else:
                self.prime_pow.eis_shift(&ans.value, &self.value, -n, self.prime_pow.prec_cap)
        return ans

    def __lshift__(pAdicZZpXFMElement self, shift):
        """
        Multiply ``self`` by the uniformizer raised to the power ``n``.

        If ``n`` is negative, right shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24
            sage: z << 17  # indirect doctest
            w^17 + w^22 + w^23 + 2*w^24
            sage: z << (-1)
            w^4 + w^5 + 2*w^6 + 4*w^7 + 3*w^9 + w^11 + 4*w^12 + 4*w^13 + 4*w^14 + 4*w^15 + 4*w^16 + 4*w^19 + w^20 + 4*w^23 + 4*w^24
        """
        cdef pAdicZZpXFMElement ans
        if not isinstance(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            #Assuming that _new_c() initializes to zero.
            return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicZZpXFMElement _rshift_c(self, long n):
        """
        Divide ``self`` by the uniformizer raised to the power ``n``.

        Throws away the nonpositive part of the series expansion.
        The top digits will be garbage.  If ``n`` is negative, left
        shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpFM(5,5,print_mode='digits')
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: for m in range(26): '...' + repr(z >> m)[len(repr(z >> m)) - 25 + m:]  # indirect doctest
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
        if n < 0:
            return self._lshift_c(-n)
        if n == 0:
            return self
        cdef pAdicZZpXFMElement ans = self._new_c()
        if n < self.prime_pow.ram_prec_cap:
            if self.prime_pow.e == 1:
                ZZ_pX_right_pshift(ans.value, self.value, self.prime_pow.pow_ZZ_tmp(n)[0], self.prime_pow.get_top_context().x)
            else:
                # Why is this not eis_shift_capdiv?!!
                self.prime_pow.eis_shift(&ans.value, &self.value, n, self.prime_pow.prec_cap)
        return ans

    def __rshift__(pAdicZZpXFMElement self, shift):
        """
        Divide ``self`` by the uniformizer raised to the power ``n``.

        Throws away the nonpositive part of the series expansion.
        The top digits will be garbage.  If ``n`` is negative, left
        shifts by ``-n``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5
            sage: z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14
             + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24
            sage: z >> (6)  # indirect doctest
            1 + 2*w + 4*w^2 + 3*w^4 + w^6 + 4*w^7 + 4*w^8 + 4*w^9 + 4*w^10 + 4*w^11
             + 4*w^14 + w^15 + 4*w^18 + 4*w^19 + 2*w^20 + 3*w^21 + 2*w^22 + 3*w^24
            sage: z >> (-4)
            w^4 + w^9 + w^10 + 2*w^11 + 4*w^12 + 3*w^14 + w^16 + 4*w^17
             + 4*w^18 + 4*w^19 + 4*w^20 + 4*w^21 + 4*w^24
        """
        cdef pAdicZZpXFMElement ans
        if not isinstance(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cpdef _neg_(self):
        """
        Return ``-self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = (1 + w)^5; z
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14
             + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24
            sage: -z  # indirect doctest
            4 + 3*w^5 + 4*w^6 + w^7 + w^8 + w^9 + w^10 + w^11 + 2*w^12 + 4*w^13
             + 4*w^15 + 3*w^16 + w^17 + 2*w^18 + 3*w^19 + 2*w^21 + 4*w^23 + 4*w^24
            sage: y = z + (-z); y
            0
            sage: -y
            0
        """
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_negate(ans.value, self.value)
        return ans

    def __pow__(pAdicZZpXFMElement self, right, m): # m ignored
        """
        Compute ``self`` ^ ``right``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: (1 + w)^5  # indirect doctest
            1 + w^5 + w^6 + 2*w^7 + 4*w^8 + 3*w^10 + w^12 + 4*w^13 + 4*w^14
             + 4*w^15 + 4*w^16 + 4*w^17 + 4*w^20 + w^21 + 4*w^24
            sage: (1 + w)^-5
            1 + 4*w^5 + 4*w^6 + 3*w^7 + w^8 + 2*w^10 + w^11 + w^12 + 2*w^14
             + 3*w^16 + 3*w^17 + 4*w^18 + 4*w^19 + 2*w^20 + 2*w^21 + 4*w^22 + 3*w^23 + 3*w^24

        TESTS:

        We define ``0^0`` to be unity, :issue:`13786`::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: type(W(0))
            <class 'sage.rings.padics.padic_ZZ_pX_FM_element.pAdicZZpXFMElement'>
            sage: W(0)^0
            1
            sage: W(0)^0 == W(1)
            True

        The value returned from ``0^0`` should belong to our ring::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: type(W(0)^0) == type(W(0))
            True
        """
        if not isinstance(right, Integer):
            right = Integer(right)
        if right == 0 and self == 0:
            return self.parent(1)
        cdef pAdicZZpXFMElement ans = self._new_c()
        cdef ntl_ZZ rZZ = ntl_ZZ.__new__(ntl_ZZ)
        mpz_to_ZZ(&rZZ.x, (<Integer>right).value)
        if mpz_sgn((<Integer>right).value) < 0:
            if self.valuation_c() > 0:
                raise ValueError("cannot invert non-unit")
            sig_on()
            if self.prime_pow.e == 1:
                ZZ_pX_InvMod_newton_unram(ans.value, self.value, self.prime_pow.get_top_modulus()[0], self.prime_pow.get_top_context().x, self.prime_pow.get_context(1).x)
            else:
                ZZ_pX_InvMod_newton_ram(ans.value, self.value, self.prime_pow.get_top_modulus()[0], self.prime_pow.get_top_context().x)
            ZZ_negate(rZZ.x, rZZ.x)
            ZZ_pX_PowerMod_pre(ans.value, ans.value, rZZ.x, self.prime_pow.get_top_modulus()[0])
            sig_off()
        else:
            sig_on()
            ZZ_pX_PowerMod_pre(ans.value, self.value, rZZ.x, self.prime_pow.get_top_modulus()[0])
            sig_off()
        return ans

    cpdef _add_(self, right):
        """
        Return ``self`` + ``right``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11) - 69  # indirect doctest
            1 + 4*w^13 + 2*w^16 + 4*w^17 + 3*w^18 + 4*w^20 + 4*w^22
            sage: -69 + (4*w^5 + 3*w^7 + w^9 + 2*w^10 + 2*w^11)
            1 + 4*w^13 + 2*w^16 + 4*w^17 + 3*w^18 + 4*w^20 + 4*w^22
        """
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_add(ans.value, self.value, (<pAdicZZpXFMElement>right).value)
        return ans

    cpdef _mul_(self, right):
        """
        Return the product of ``self`` and ``right``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(329)
            sage: b = W(111)
            sage: a*b  # indirect doctest
            4 + 3*w^5 + w^7 + 2*w^9 + 4*w^11 + 3*w^12 + 2*w^13 + w^14 + 2*w^15
             + 3*w^16 + 4*w^17 + 4*w^18 + 2*w^19 + 2*w^21 + 4*w^22 + 2*w^23 + w^24
            sage: a * 0
            0
            sage: W(125) * W(375)
            0
        """
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_MulMod_pre(ans.value, self.value, (<pAdicZZpXFMElement>right).value, self.prime_pow.get_top_modulus()[0])
        return ans

    cpdef _sub_(self, right):
        """
        Return the difference of ``self`` and ``right``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(329)
            sage: b = W(111)
            sage: a - b  # indirect doctest
            3 + 3*w^5 + w^7 + 2*w^9 + 3*w^10 + 4*w^11 + 2*w^13 + 2*w^14 + w^15
             + 4*w^16 + 2*w^18 + 3*w^19 + 2*w^20 + 3*w^21 + w^22 + w^24
            sage: W(218)
            3 + 3*w^5 + w^7 + 2*w^9 + 3*w^10 + 4*w^11 + 2*w^13 + 2*w^14 + w^15
             + 4*w^16 + 2*w^18 + 3*w^19 + 2*w^20 + 3*w^21 + w^22 + w^24
        """
        cdef pAdicZZpXFMElement ans = self._new_c()
        ZZ_pX_sub(ans.value, self.value, (<pAdicZZpXFMElement>right).value)
        return ans

    cpdef _div_(self, _right):
        """
        Return the quotient of ``self`` by ``right``.

        If ``right`` is not a unit, raises a :exc:`ValueError`.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: W(125) / W(14)  # indirect doctest
            4*w^15 + 4*w^17 + w^19 + w^20 + w^23 + 2*w^24
            sage: 1 / W(14) == ~W(14)
            True
            sage: 1 / w
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit

        We check that :issue:`11403` has been resolved::

            sage: R.<t> = Zq(8,2,'fixed-mod')
            sage: 1/(t+t^2)
            (t + 1) + t^2*2
        """
        cdef pAdicZZpXFMElement right = <pAdicZZpXFMElement>_right
        if right.valuation_c() > 0:
            raise ValueError("cannot invert non-unit")
        cdef pAdicZZpXFMElement ans = self._new_c()
        sig_on()
        if self.prime_pow.e == 1:
            ZZ_pX_InvMod_newton_unram(ans.value, right.value, self.prime_pow.get_top_modulus()[0],
                                      self.prime_pow.get_top_context().x, self.prime_pow.get_context(1).x)
        else:
            ZZ_pX_InvMod_newton_ram(ans.value, right.value, self.prime_pow.get_top_modulus()[0], self.prime_pow.get_top_context().x)
        ZZ_pX_MulMod_pre(ans.value, self.value, ans.value, self.prime_pow.get_top_modulus()[0])
        sig_off()
        return ans

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: b = W(45); b
            4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16
             + 2*w^17 + w^19 + 4*w^20 + w^21 + 3*w^22 + 3*w^23 + 4*w^24
            sage: c = copy(b); c
            4*w^5 + 3*w^7 + w^9 + w^10 + 2*w^11 + w^12 + w^13 + 3*w^14 + w^16
             + 2*w^17 + w^19 + 4*w^20 + w^21 + 3*w^22 + 3*w^23 + 4*w^24
            sage: c is b
            False
        """
        cdef pAdicZZpXFMElement ans = self._new_c()
        ans.value = self.value # does this actually copy correctly
        return ans

    def is_zero(self, absprec=None):
        """
        Return whether the valuation of ``self`` is at least
        ``absprec``; if ``absprec`` is ``None``, return whether
        ``self`` is indistinguishable from zero.

        EXAMPLES::

            sage: R = ZpFM(5,5)
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
        if absprec is None:
            ans = ZZ_pX_IsZero(self.value)
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0:
                    return True
                else:
                    ans = ZZ_pX_IsZero(self.value)
            else:
                ans = (self.valuation_c() >= mpz_get_si((<Integer>absprec).value))
        return ans

    def add_bigoh(self, absprec):
        r"""
        Return a new element truncated modulo `\pi^{\text{absprec}}`.

        This is only implemented for unramified extension at
        this point.

        INPUT:

        - ``absprec`` -- integer

        OUTPUT: a new element truncated modulo `\pi^{\mbox{absprec}}`

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod')
            sage: a = R(1 + 7 + 7^2)
            sage: a.add_bigoh(1)
            1
        """
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_cmp_ui((<Integer>absprec).value, self.prime_pow.prec_cap) >= 0:
            return self

        cdef pAdicZZpXFMElement ans = self._new_c()
        cdef ZZ_pX_c tmp
        cdef ntl_ZZ_pContext_class c
        cdef unsigned long aprec
        if self.prime_pow.e == 1:
            if mpz_fits_ulong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0:
                    return ans # assumes _new_c() initializes to 0
                return self # absprec > prec_cap
            aprec = mpz_get_ui((<Integer>absprec).value)
            if aprec >= self.prime_pow.prec_cap:
                return self
            c = self.prime_pow.get_context(aprec)
            c.restore_c()
            ZZ_pX_conv_modulus(tmp, self.value, c.x)
            ZZ_pX_conv_modulus(ans.value, tmp, (<ntl_ZZ_pContext_class>self.prime_pow.get_top_context()).x)
        else:
            raise NotImplementedError
        return ans

    def _integer_(self, Z=None):
        r"""
        Return an integer congruent to this element modulo
        `\pi^a`, where `a` is ``self.absolute_precision()``, if possible.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: ZZ(ZqFM(125,names='a')(-1))  # indirect doctest
            95367431640624
            sage: R = ZpFM(5); S.<x> = ZZ[]; f = x^5 + 25*x^3 - 5; W.<w> = R.ext(f)
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
        cdef ZZ_c tmp_z
        if ZZ_pX_deg(self.value) > 0:
            raise ValueError("This element not well approximated by an integer.")
        cdef Integer ans = Integer.__new__(Integer)
        tmp_z = ZZ_p_rep(ZZ_pX_ConstTerm(self.value))
        ZZ_to_mpz(ans.value, &tmp_z)
        return ans

    def matrix_mod_pn(self):
        r"""
        Return the matrix of right multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for this
        extension field.

        The **rows** of this matrix give the images of each of the `x^i`.
        The entries of the matrices are ``IntegerMod`` elements,
        defined modulo ``p^(self.absprec() / e)``.

        Raises an error if ``self`` has negative valuation.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = (3+w)^7
            sage: a.matrix_mod_pn()                                                     # needs sage.geometry.polyhedron
            [2757  333 1068  725 2510]
            [  50 1507  483  318  725]
            [ 500   50 3007 2358  318]
            [1590 1375 1695 1032 2358]
            [2415  590 2370 2970 1032]
        """
        from sage.matrix.constructor import matrix
        R = IntegerModRing(self.prime_pow.pow_Integer(self.prime_pow.prec_cap))
        n = self.prime_pow.deg
        L = []
        cdef ntl_ZZ_pX cur = <ntl_ZZ_pX>self._ntl_rep()
        cur.c.restore_c()
        cdef ZZ_pX_Modulus_c* m = self.prime_pow.get_top_modulus()
        cdef ZZ_pX_c x
        ZZ_pX_SetX(x)
        cdef Py_ssize_t i
        zero = int(0)
        for i in range(n):
            curlist = cur.list()
            L.extend(curlist + [zero]*(n - len(curlist)))
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

    def norm(self, base=None):
        r"""
        Return the absolute or relative norm of this element.

        .. NOTE::

            This is not the `p`-adic absolute value.  This is a
            field theoretic norm down to a ground ring.

        If you want the `p`-adic absolute value, use the :func:`abs()`
        function instead.

        If `K` is given then `K` must be a subfield of the parent `L` of
        ``self``, in which case the norm is the relative norm from `L` to `K`.
        In all other cases, the norm is the absolute norm down to `\QQ_p`
        or `\ZZ_p`.

        EXAMPLES::

            sage: R = ZpCR(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: ((1+2*w)^5).norm()
            1 + 5^2 + O(5^5)
            sage: ((1+2*w)).norm()^5
            1 + 5^2 + O(5^5)
        """
        if base is not None:
            if base is self.parent():
                return self
            else:
                raise NotImplementedError
        if self._is_exact_zero():
            return self.parent().ground_ring()(0)
        elif self._is_inexact_zero():
            return self.ground_ring(0, self.valuation())
        norm_of_uniformizer = (-1)**self.parent().degree() * self.parent().defining_polynomial()[0]
        return self.parent().ground_ring()(self.unit_part().matrix_mod_pn().det()) * norm_of_uniformizer**self.valuation()

    def trace(self, base=None):
        r"""
        Return the absolute or relative trace of this element.

        If `K` is given then `K` must be a subfield of the parent `L` of
        ``self``, in which case the norm is the relative norm from `L` to `K`.
        In all other cases, the norm is the absolute norm down to `\QQ_p`
        or `\ZZ_p`.

        EXAMPLES::

            sage: R = ZpCR(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = (2+3*w)^7
            sage: b = (6+w^3)^5
            sage: a.trace()
            3*5 + 2*5^2 + 3*5^3 + 2*5^4 + O(5^5)
            sage: a.trace() + b.trace()
            4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
            sage: (a+b).trace()
            4*5 + 5^2 + 5^3 + 2*5^4 + O(5^5)
        """
        if base is not None:
            if base is self.parent():
                return self
            else:
                raise NotImplementedError
        if self._is_exact_zero():
            return self.parent().ground_ring()(0)
        elif self._is_inexact_zero():
            return self.ground_ring(0, (self.valuation() - 1) // self.parent().e() + 1)
        if self.valuation() >= 0:
            return self.parent().ground_ring()(self.matrix_mod_pn().trace())
        else:
            shift = -(self.valuation() // self.parent().e())
            return self.parent().ground_ring()((self * self.parent().prime() ** shift).matrix_mod_pn().trace()) / self.parent().prime()**shift

    def _ntl_rep(self):
        """
        Return an ``ntl_ZZ_pX`` holding ``self.value``.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = 72 + 4*w^2; b = 17 + 9*w + w^3; c = a + b
            sage: c._ntl_rep()
            [89 9 4 1]
        """
        self.prime_pow.restore_top_context()
        cdef ntl_ZZ_pX ans = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        ans.c = self.prime_pow.get_top_context()
        ans.x = self.value
        return ans

    def _polynomial_list(self, pad=False):
        """
        Return the coefficient list for a polynomial over the base ring
        yielding this element.

        INPUT:

        - ``pad`` -- whether to pad the result with zeros of the
          appropriate precision

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: W.<w> = ZpFM(5).extension(x^3 - 5)
            sage: (1 + w)._polynomial_list()
            [1, 1]
            sage: (1 + w)._polynomial_list(pad=True)
            [1, 1, 0]
        """
        R = self.base_ring()
        if self.is_zero():
            L = []
        else:
            L = [Integer(c) for c in self._ntl_rep().list()]
        if pad:
            n = self.parent().degree()
            L.extend([R.zero()] * (n - len(L)))
        return L

    def polynomial(self, var='x'):
        """
        Return a polynomial over the base ring that yields this element
        when evaluated at the generator of the parent.

        INPUT:

        - ``var`` -- string, the variable name for the polynomial

        EXAMPLES::

            sage: S.<x> = ZZ[]
            sage: W.<w> = ZpFM(5).extension(x^2 - 5)
            sage: (w + 5).polynomial()
            x + 5
        """
        R = self.base_ring()
        S = R[var]
        return S(self._polynomial_list())

    cdef ZZ_p_c _const_term(self) noexcept:
        """
        Return the constant term of ``self.unit``.

        .. NOTE::

            This may be divisible by `p` if ``self`` is not normalized.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566)
            sage: a._const_term_test()  # indirect doctest
            566
        """
        return ZZ_pX_ConstTerm(self.value)

    def is_equal_to(self, right, absprec=None):
        """
        Return whether ``self`` is equal to ``right`` modulo
        ``self.uniformizer()^absprec``.

        If ``absprec`` is ``None``, returns if ``self`` is equal to
        ``right`` modulo the precision cap.

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
#        Returns an element of a number field defined by the same polynomial as
#        ``self``'s parent that is congruent to self modulo an appropriate ideal.

#        Not currently implemented.
#        """
#        raise NotImplementedError

    def lift_to_precision(self, absprec=None):
        """
        Return ``self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w.lift_to_precision(10000)
            w
        """
        return self

    def expansion(self, n=None, lift_mode='simple'):
        r"""
        Return a list giving a series representation of this element.

        - If ``lift_mode == 'simple' or 'smallest'``, the returned list will
          consist of

          * integers (in the eisenstein case) or

          * lists of integers (in the unramified case).

        - this element can be reconstructed as

          * a sum of elements of the list times powers of the
            uniformiser (in the eisenstein case), or

          * as a sum of powers of the `p` times polynomials in the
            generator (in the unramified case).

        - If ``lift_mode == 'simple'``, all integers will be in the range
          `[0,p-1]`,

        - If ``lift_mode == 'smallest'`` they will be in the range `[(1-p)/2,
          p/2]`.

        - If ``lift_mode == 'teichmuller'``, returns a list of
          ``pAdicZZpXCRElements``, all of which are Teichmuller representatives
          and such that this element is the sum of that list times powers of the
          uniformizer.

        INPUT:

        - ``n`` -- integer (default: ``None``); if given, returns the
          corresponding entry in the expansion

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W(775); y
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24
            sage: (y>>9).expansion()
            [0, 1, 0, 4, 0, 2, 1, 2, 4, 1, 0, 1, 2, 3, 1, 1, 4, 1, 2, 4, 1, 0, 0, 3]
            sage: (y>>9).expansion(lift_mode='smallest')
            [0, 1, 0, -1, 0, 2, 1, 2, 0, 1, 2, 1, 1, -1, -1, 2, -2, 0, -2, -2, -2, 0, -2, -2, 2]
            sage: w^10 - w^12 + 2*w^14 + w^15 + 2*w^16 + w^18 + 2*w^19 + w^20 + w^21 - w^22 - w^23 + 2*w^24
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24
            sage: g = x^3 + 3*x + 3

            sage: # needs sage.libs.flint
            sage: A.<a> = R.ext(g)
            sage: y = 75 + 45*a + 1200*a^2; y
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4
            sage: E = y.expansion(); E
            5-adic expansion of 4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4
            sage: list(E)
            [[], [0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1]]
            sage: list(y.expansion(lift_mode='smallest'))
            [[], [0, -1], [-2, 2, -2], [1], [0, 0, 2]]
            sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4
            sage: W(0).expansion()
            []
            sage: list(A(0,4).expansion())
            []

        Check that :issue:`25879` has been resolved::

            sage: K = ZpCA(3,5)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^2 - 3)
            sage: a.residue()
            0
        """
        if lift_mode == 'teichmuller':
            zero = self.parent()(0)
        elif self.prime_pow.e == 1:
            zero = []
        else:
            zero = Integer(0)
        ordp = self.valuation()
        if isinstance(n, slice):
            return self.slice(n.start, n.stop, n.step)
        elif n is not None:
            if self.is_zero() or n >= self.prime_pow.ram_prec_cap or n < ordp:
                return zero
        if self.is_zero():
            return []
        if lift_mode == 'simple':
            ulist = self.ext_p_list(pos=True)
        elif lift_mode == 'smallest':
            ulist = self.ext_p_list(pos=False)
        elif lift_mode == 'teichmuller':
            if n is None:
                ulist = self.teichmuller_list()
            else:
                return self.teichmuller_expansion(n)
        else:
            raise ValueError("lift mode must be one of 'simple', 'smallest' or 'teichmuller'")
        if n is not None:
            try:
                return ulist[n - ordp]
            except IndexError:
                return zero
        return [zero] * ordp + ulist

    def teichmuller_expansion(self, n=None):
        r"""
        Return a list `[a_0, a_1, \ldots, a_n]` such that.

        - `a_i^q = a_i`
        - ``self.unit_part()`` = `\sum_{i = 0}^n a_i \pi^i`, where `\pi` is a
          uniformizer of ``self.parent()``

        INPUT:

        - ``n`` -- integer (default: ``None``); if given, returns the corresponding
          entry in the expansion

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<a> = ZqFM(5^4,4)
            sage: E = a.teichmuller_expansion(); E
            5-adic expansion of a (teichmuller)
            sage: list(E)
            [a + (2*a^3 + 2*a^2 + 3*a + 4)*5 + (4*a^3 + 3*a^2 + 3*a + 2)*5^2 + (4*a^2 + 2*a + 2)*5^3,
             (3*a^3 + 3*a^2 + 2*a + 1) + (a^3 + 4*a^2 + 1)*5 + (a^2 + 4*a + 4)*5^2 + (4*a^2 + a + 3)*5^3,
             (4*a^3 + 2*a^2 + a + 1) + (2*a^3 + 2*a^2 + 2*a + 4)*5 + (3*a^3 + 2*a^2 + a + 1)*5^2 + (a^3 + a^2 + 2)*5^3,
             (a^3 + a^2 + a + 4) + (3*a^3 + 1)*5 + (3*a^3 + a + 2)*5^2 + (3*a^3 + 3*a^2 + 3*a + 1)*5^3]
            sage: sum([c * 5^i for i, c in enumerate(E)])
            a
            sage: all(c^625 == c for c in E)
            True

            sage: S.<x> = ZZ[]
            sage: f = x^3 - 98*x + 7
            sage: W.<w> = ZpFM(7,3).ext(f)
            sage: b = (1+w)^5; L = b.teichmuller_expansion(); L
            [1,
             5 + 5*w^3 + w^6 + 4*w^7,
             3 + 3*w^3 + w^7,
             3 + 3*w^3 + w^7,
             0,
             4 + 5*w^3 + w^6 + 4*w^7,
             3 + 3*w^3 + w^7,
             6 + w^3 + 5*w^7,
             6 + w^3 + 5*w^7]
            sage: sum([w^i*L[i] for i in range(len(L))]) == b
            True
            sage: all(L[i]^(7^3) == L[i] for i in range(9))
            True

            sage: L = W(3).teichmuller_expansion(); L
            [3 + 3*w^3 + w^7,
             0,
             0,
             4 + 5*w^3 + w^6 + 4*w^7,
             0,
             0,
             3 + 3*w^3 + w^7,
             6 + w^3 + 5*w^7]
            sage: sum([w^i*L[i] for i in range(len(L))])
            3
        """
        cdef long ordp = self.valuation_c()
        cdef long rp = self.prime_pow.ram_prec_cap - ordp
        cdef long goal
        if n is not None:
            goal = self.ram_prec_cap - n
        cdef pAdicZZpXFMElement v
        if n is None:
            L = []
            if ZZ_pX_IsZero(self.value):
                return L
        elif n < ordp or n >= self.prime_pow.ram_prec_cap:
            return self.parent()(0)
        else:
            v = self._new_c()
        cdef pAdicZZpXFMElement u = self.unit_part()
        if u is self:
            u = self.__copy__()
        while u.valuation_c() < rp:
            if n is None:
                v = self._new_c()
            self.prime_pow.teichmuller_set_c(&v.value, &u.value, self.prime_pow.ram_prec_cap)
            if n is None:
                L.append(v)
            elif rp == goal:
                return v
            if rp == 1:
                break
            ZZ_pX_sub(u.value, u.value, v.value)
            rp -= 1
            if self.prime_pow.e == 1:
                ZZ_pX_right_pshift(u.value, u.value, self.prime_pow.pow_ZZ_tmp(1)[0], self.prime_pow.get_top_context().x)
            else:
                self.prime_pow.eis_shift(&u.value, &u.value, 1, self.prime_pow.ram_prec_cap)
        if n is None:
            return L
        else:
            return self.parent()(0)

    def _teichmuller_set_unsafe(self):
        """
        Set this element to the Teichmuller representative with the
        same residue.

        .. WARNING::

            This function modifies the element, which is not safe.
            Elements are supposed to be immutable.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W.teichmuller(3); y  # indirect doctest
            3 + 3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + 3*w^15
             + 2*w^16 + 3*w^17 + w^18 + 3*w^19 + 3*w^20 + 2*w^21 + 2*w^22 + 3*w^23 + 4*w^24
            sage: y^5 == y
            True

            sage: # needs sage.libs.flint
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = A.teichmuller(1 + 2*a - a^2); b
            (4*a^2 + 2*a + 1) + 2*a*5 + (3*a^2 + 1)*5^2 + (a + 4)*5^3 + (a^2 + a + 1)*5^4
            sage: b^125 == b
            True
        """
        if self.valuation_c() > 0:
            ZZ_pX_clear(self.value)
        else:
            self.prime_pow.teichmuller_set_c(&self.value, &self.value, self.prime_pow.ram_prec_cap)

#     def multiplicative_order(self):
#         """
#         Returns the multiplicative order of self, ie the smallest
#         positive `n` so that there is an exact `p`-adic element congruent
#         to self modulo ``self``'s precision that is an `n`-th root of unity.

#         Note: unlike the case for Qp and Zp, it is possible to have
#         non-teichmuller elements with finite orders.  This can happen
#         only if (p-1) divides the ramification index (see the
#         documentation on __pow__).

#         INPUT:

#         - ``self`` -- a `p`-adic element
#         - ``prec`` -- integer

#         OUTPUT: integer; the multiplicative order of self
#         """
#         raise NotImplementedError

#     def padded_list(self, n, lift_mode='simple'):
#         """
#         Returns a list of coefficients of pi starting with `pi^0` up to
#         `pi^n` exclusive (padded with zeros if needed)

#         """
#         raise NotImplementedError

    def precision_absolute(self):
        """
        Return the absolute precision of ``self``, ie the precision cap
        of ``self.parent()``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + 3*w^19 + 2*w^21 + 3*w^22 + 3*w^23
            sage: a.valuation()
            10
            sage: a.precision_absolute()
            25
            sage: a.precision_relative()
            15
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + 3*w^9 + 2*w^11 + 3*w^12
             + 3*w^13 + w^15 + 4*w^16 + 2*w^17 + w^18 + 3*w^21 + w^22 + 3*w^24
        """
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set_ui(ans.value, self.prime_pow.ram_prec_cap)
        return ans

    def precision_relative(self):
        """
        Return the relative precision of ``self``, ie the precision cap
        of ``self.parent()`` minus the ``valuation of self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + 3*w^19 + 2*w^21 + 3*w^22 + 3*w^23
            sage: a.valuation()
            10
            sage: a.precision_absolute()
            25
            sage: a.precision_relative()
            15
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + 3*w^9 + 2*w^11 + 3*w^12
             + 3*w^13 + w^15 + 4*w^16 + 2*w^17 + w^18 + 3*w^21 + w^22 + 3*w^24
        """
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set_ui(ans.value, self.prime_pow.ram_prec_cap - self.valuation_c())
        return ans

    cpdef pAdicZZpXFMElement unit_part(self):
        """
        Return the unit part of ``self``, ie
        ``self / uniformizer^(self.valuation())``

        .. WARNING::

            If this element has positive valuation then the unit part
            is not defined to the full precision of the ring.  Asking
            for the unit part of ``ZpFM(5)(0)`` will not raise an error,
            but rather return itself.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + 3*w^19 + 2*w^21 + 3*w^22 + 3*w^23
            sage: a.valuation()
            10
            sage: a.precision_absolute()
            25
            sage: a.precision_relative()
            15
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + 3*w^9 + 2*w^11 + 3*w^12
             + 3*w^13 + w^15 + 4*w^16 + 2*w^17 + w^18 + 3*w^21 + w^22 + 3*w^24

        The unit part inserts nonsense digits if this element has
        positive valuation::

            sage: (a-a).unit_part()
            0
        """
        return self._rshift_c(self.valuation_c())

    cdef long valuation_c(self) noexcept:
        """
        Return the valuation of ``self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(75); a
            3*w^10 + 2*w^12 + w^14 + w^16 + w^17 + 3*w^18 + 3*w^19 + 2*w^21 + 3*w^22 + 3*w^23
            sage: a.valuation()
            10
            sage: a.precision_absolute()
            25
            sage: a.precision_relative()
            15
            sage: a.unit_part()
            3 + 2*w^2 + w^4 + w^6 + w^7 + 3*w^8 + 3*w^9 + 2*w^11 + 3*w^12
             + 3*w^13 + w^15 + 4*w^16 + 2*w^17 + w^18 + 3*w^21 + w^22 + 3*w^24
        """
        cdef long valuation = 0, index = 0
        ZZ_pX_min_val_coeff(valuation, index, self.value, self.prime_pow.pow_ZZ_tmp(1)[0])
        if index == -1: # self == 0
            return self.prime_pow.ram_prec_cap
        if self.prime_pow.e == 1:
            return valuation
        else:
            if index + valuation * self.prime_pow.e >= self.prime_pow.ram_prec_cap:
                return self.prime_pow.ram_prec_cap
            else:
                return index + valuation * self.prime_pow.e

    cdef ext_p_list(self, bint pos):
        r"""
        Return a list giving a series representation of ``self``.

        - The returned list will consist of

          * integers (in the eisenstein case) or

          * a lists of integers (in the unramified case).

        - ``self`` can be reconstructed

          * as a sum of elements of the list times powers of the
            uniformiser (in the eisenstein case), or

          * as a sum of powers of `p` times polynomials in the
            generator (in the unramified case).

        - If ``pos`` is ``True``, all integers will be in the range `[0,p-1]`,
        otherwise they will be in the range `[(1-p)/2, p/2]`.

        Note that zeros are truncated from the returned list, so you
        must use the valuation function to fully reconstruct ``self``.

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W(775); y
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24
            sage: (y>>9).expansion()  # indirect doctest
            [0, 1, 0, 4, 0, 2, 1, 2, 4, 1, 0, 1, 2, 3, 1, 1, 4, 1, 2, 4, 1, 0, 0, 3]
            sage: (y>>9).expansion(lift_mode='smallest')  # indirect doctest
            [0, 1, 0, -1, 0, 2, 1, 2, 0, 1, 2, 1, 1, -1, -1, 2, -2, 0, -2, -2, -2, 0, -2, -2, 2]
            sage: w^10 - w^12 + 2*w^14 + w^15 + 2*w^16 + w^18 + 2*w^19 + w^20 + w^21 - w^22 - w^23 + 2*w^24
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + w^20 + 2*w^21 + 3*w^22 + w^23 + w^24
            sage: g = x^3 + 3*x + 3

            sage: # needs sage.libs.flint
            sage: A.<a> = R.ext(g)
            sage: y = 75 + 45*a + 1200*a^2; y
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4
            sage: list(y.expansion())  # indirect doctest
            [[], [0, 4], [3, 1, 3], [0, 0, 4], [0, 0, 1]]
            sage: list(y.expansion(lift_mode='smallest'))  # indirect doctest
            [[], [0, -1], [-2, 2, -2], [1], [0, 0, 2]]
            sage: 5*((-2*5 + 25) + (-1 + 2*5)*a + (-2*5 + 2*125)*a^2)
            4*a*5 + (3*a^2 + a + 3)*5^2 + 4*a^2*5^3 + a^2*5^4
        """
        return self.ext_p_list_precs(pos, self.prime_pow.ram_prec_cap)


def make_ZZpXFMElement(parent, f):
    """
    Create a new ``pAdicZZpXFMElement`` out of an ``ntl_ZZ_pX`` ``f``, with
    parent ``parent``.  For use with pickling.

    EXAMPLES::

        sage: R = ZpFM(5,5)
        sage: S.<x> = R[]
        sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
        sage: W.<w> = R.ext(f)
        sage: z = (1 + w)^5 - 1
        sage: loads(dumps(z)) == z  # indirect doctest
        True
    """
    return pAdicZZpXFMElement(parent, f)
