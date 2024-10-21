# sage_setup: distribution = sagemath-ntl
# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
# sage.doctest: needs sage.rings.padics
r"""
`p`-adic ``ZZ_pX Element``

A common superclass implementing features shared by all elements that
use NTL's ``ZZ_pX`` as the fundamental data type.

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.ext.stdsage cimport PY_NEW
from cpython.list cimport *
include "sage/libs/ntl/decl.pxi"

from sage.libs.gmp.mpz cimport *
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.finite_rings.integer_mod import IntegerMod_abstract
from sage.rings.padics.padic_printing cimport pAdicPrinter_class
from sage.rings.padics.pow_computer_ext cimport PowComputer_ext
from sage.rings.rational_field import QQ

zero = Integer(0)
one = Integer(1)
two = Integer(2)
big = two**128 + one
#this should not fit in a long, since it's supposed to be bigger than any valid absolute precision.

cdef class pAdicZZpXElement(pAdicExtElement):
    def __init__(self, parent):
        """
        Initialization.

        EXAMPLES::

            sage: A = Zp(next_prime(50000),10)
            sage: S.<x> = A[]
            sage: B.<t> = A.ext(x^2 + next_prime(50000))  # indirect doctest
        """
        self.prime_pow = <PowComputer_ZZ_pX>parent.prime_pow
        pAdicExtElement.__init__(self, parent)

    cdef int _set_from_list(self, L) except -1:
        """
        Set ``self`` from a list.

        The list can contain integers, ``IntegerMods``, rationals, or
        `p`-adic base elements

        INPUT:

        - ``L`` -- list

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: S.<x> = ZZ[]
            sage: W.<w> = R.ext(x^5 + 25*x^3 - 15*x - 5)
            sage: W([1,2,3,4])  # indirect doctest
            1 + 2*w + 3*w^2 + 4*w^3
            sage: W([5,10,15,20])
            w^5 + 4*w^6 + w^7 + w^8 + 2*w^9 + 4*w^10 + 2*w^11 + 3*w^13 + 2*w^15 + w^16 + 2*w^17 + 2*w^18 + w^19 + 4*w^20 + w^21 + 4*w^22 + 4*w^23 + 2*w^24
        """
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX((<ntl_ZZX>ntl_ZZX(L)).x)
        else:
            self._set_from_ZZ_pX(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx)
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef int _set_from_list_rel(self, L, long relprec) except -1:
        """
        Set ``self`` from a list.

        The list can contain integers, ``IntegerMods``, rationals, or
        `p`-adic base elements

        INPUT:

        - ``L`` -- list

        - ``relprec`` -- integer; capping the relative precision of ``self``

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = ZZ[]
            sage: W.<w> = R.ext(x^5 + 25*x^3 - 15*x - 5)
            sage: W([1,2,3,4])  # indirect doctest
            1 + 2*w + 3*w^2 + 4*w^3 + O(w^25)
            sage: W([5,10,15,20], relprec=16)
            w^5 + 4*w^6 + w^7 + w^8 + 2*w^9 + 4*w^10 + 2*w^11 + 3*w^13 + 2*w^15 + w^16 + 2*w^17 + 2*w^18 + w^19 + 4*w^20 + O(w^21)
        """
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX_rel((<ntl_ZZX>ntl_ZZX(L)).x, relprec)
        else:
            self._set_from_ZZ_pX_rel(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx, relprec)
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef int _set_from_list_abs(self, L, long absprec) except -1:
        """
        Set this element from a list.

        INPUT:

        - ``L`` -- list of integers, ``IntegerMod``s, rationals, or `p`-adic
          base elements

        - ``absprec`` -- integer at which the absolute precision of the
          result will be capped

        EXAMPLES::

            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
            sage: W.<w> = R.ext(x^5 + 25*x^3 - 15*x - 5)
            sage: W([1,2,3,4])
            1 + 2*w + 3*w^2 + 4*w^3 + O(w^25)
            sage: W([5,10,15,20], absprec=16)  # indirect doctest
            w^5 + 4*w^6 + w^7 + w^8 + 2*w^9 + 4*w^10 + 2*w^11 + 3*w^13 + 2*w^15 + O(w^16)
        """
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX_abs((<ntl_ZZX>ntl_ZZX(L)).x, absprec)
        else:
            self._set_from_ZZ_pX_abs(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx, absprec - (min_val * self.parent().e()))
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef int _set_from_list_both(self, L, long absprec, long relprec) except -1:
        """
        Set this element from a list.

        The list can contain integers, ``IntegerMods``, rationals, or
        `p`-adic base elements

        INPUT:

        - ``L`` -- list of integers, ``IntegerMod``s, rationals, or `p`-adic
          base elements

        - ``absprec`` -- integer at which the absolute precision of the result will be capped

        - ``relprec`` -- integer at which the relative precision of the result will be capped

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = ZZ[]
            sage: W.<w> = R.ext(x^5 + 25*x^3 - 15*x - 5)
            sage: W([1,2,3,4])
            1 + 2*w + 3*w^2 + 4*w^3 + O(w^25)
            sage: W([5,10,15,20], absprec=16)  # indirect doctest
            w^5 + 4*w^6 + w^7 + w^8 + 2*w^9 + 4*w^10 + 2*w^11 + 3*w^13 + 2*w^15 + O(w^16)
            sage: T.<a> = Qp(5).extension(x^2 - 5)
            sage: T([5^-2], absprec=-1)
            a^-4 + O(a^-1)
            sage: G.<g> = Qp(5).extension(x^2 - 5)
            sage: G(a^-41)
            g^-41 + O(g^-2)
        """
        cdef ntl_ZZ_pContext_class ctx
        L, min_val, ctx = preprocess_list(self, L)
        if ctx is None:
            self._set_from_ZZX_both((<ntl_ZZX>ntl_ZZX(L)).x, absprec, relprec)
        else:
            self._set_from_ZZ_pX_both(&(<ntl_ZZ_pX>ntl_ZZ_pX(L, ctx)).x, ctx, absprec - (min_val * self.parent().e()), relprec)
            self._pshift_self(mpz_get_si((<Integer>min_val).value))

    cdef long _check_ZZ_pContext(self, ntl_ZZ_pContext_class ctx) except -1:
        """
        Check that the given ``ntl_ZZ_pContext`` is actually a power
        of the relevant prime.  If so, returns the exponent.

        INPUT:

        - ``ctx`` -- an ``ntl_ZZ_pContext_class``

        OUTPUT:

        - ``val`` -- if ``ctx`` is a context for `p^n`, returns `n`.
          Otherwise, raises an error.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = W(ntl.ZZ_pX([4,1,16],5^2)); z  # indirect doctest
            4 + w + w^2 + 3*w^7 + w^9 + O(w^10)
        """
        cdef ZZ_c leftover
        cdef long val = ZZ_remove(leftover, ctx.p.x, self.prime_pow.pow_ZZ_tmp(1)[0])
        if ZZ_IsOne(leftover):
            return val
        else:
            raise ValueError("context must be a power of the appropriate prime")

    cdef ext_p_list_precs(self, bint pos, long prec):
        """
        Return a list giving a series representation of ``self``.

        - The returned list will consist of:

          + integers (in the Eisenstein case) or

          + lists of integers (in the unramified case).

        - ``self`` can be reconstructed as

          + a sum of elements of the list times powers of the
            uniformiser (in the Eisenstein case), or

          + as a sum of powers of the `p` times polynomials in the
            generator (in the unramified case).

        Note that zeros are truncated from the returned list, so you
        must use the valuation() function to completely recover ``self``.

        INPUT:

        - ``pos`` -- bint; if ``True``, all integers will be in
          the range `[0,p-1]`, otherwise they will be in the range
          `[(1-p)/2, p/2]`

        - ``prec`` -- how many terms to return in the list.  This is
          important since shifting in the Eisenstein case can
          introduce random high order bits.  Thus the process would
          not otherwise necessarily terminate at the right point.

        OUTPUT:

        - ``L`` -- list of integers or list of lists giving the
          series expansion of ``self``

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W(775, 19); y
            w^10 + 4*w^12 + 2*w^14 + w^15 + 2*w^16 + 4*w^17 + w^18 + O(w^19)
            sage: y._ext_p_list(True)  # indirect doctest
            [1, 0, 4, 0, 2, 1, 2, 4, 1]
            sage: y._ext_p_list(False)
            [1, 0, -1, 0, 2, 1, 2, 0, 1]
        """
        ans = []
        cdef ntl_ZZ ZZ_coeff = ntl_ZZ()
        cdef Integer coeff = PY_NEW(Integer)
        cdef Integer zero = Integer(0)
        cdef Integer list_elt
        cdef ZZ_c halfp
        cdef Py_ssize_t i, j
        self.prime_pow.restore_top_context()
        ###ZZ_p_construct(&const_term_holder)
        cdef ntl_ZZ holder = ntl_ZZ()
        cdef ZZ_p_c tmp
        cdef pAdicPrinter_class printer = <pAdicPrinter_class>self.parent()._printer
        cdef ZZ_pX_c shifter = (<ntl_ZZ_pX>self._ntl_rep()).x

        #cdef ntl_ZZ_pContext_class cup = self.prime_pow.get_context(self.prime_pow.prec_cap + (<PowComputer_ZZ_pX_FM_Eis>self.prime_pow).low_length)
        #cdef ntl_ZZ_pX printer = ntl_ZZ_pX([],cup)
        #printer.x = ((<PowComputer_ZZ_pX_FM_Eis>self.prime_pow).low_shifter[0]).val()
        #print(printer)

        if self.prime_pow.e == 1:
            for j from 0 <= j < self.prime_pow.prec_cap:
                ans.append([])
            for i from 0 <= i < self.prime_pow.deg:
                ZZ_coeff.x = ZZ_p_rep(ZZ_pX_coeff(shifter, i))
                ZZ_to_mpz(coeff.value, &ZZ_coeff.x)
                L = printer.base_p_list(coeff, pos)
                for j from 0 <= j < prec:
                    if j < len(L):
                        ans[j].append(L[j])
                    else:
                        ans[j].append(zero)
            for j from 0 <= j < prec:
                while ans[j]:
                    if ans[j][-1] == 0:
                        ans[j].pop()
                    else:
                        break
            zerotest = []
        else:
            halfp = self.prime_pow.pow_ZZ_tmp(1)[0]
            ZZ_DivRem_long(halfp, halfp, 2)
            i = 0
            while True:
                # It's important that one doesn't normalize in between shifting (for capped relative elements):
                # _const_term doesn't normalize and thus we pick up the zeros
                # since we're throwing away leading zeros, it doesn't matter if we start normalized or not.
                for j from 0 <= j < self.prime_pow.e:
                    list_elt = PY_NEW(Integer)
                    if i + j == prec:
                        break
                    ZZ_rem(holder.x, ZZ_p_rep(ZZ_pX_coeff(shifter, j)), self.prime_pow.pow_ZZ_tmp(1)[0])
                    if not pos and not ZZ_IsZero(holder.x) and ZZ_compare(holder.x, halfp) > 0:
                        ZZ_sub(holder.x, self.prime_pow.pow_ZZ_tmp(1)[0], holder.x)
                        ZZ_p_add(tmp, ZZ_to_ZZ_p(holder.x), ZZ_pX_coeff(shifter, j))
                        ZZ_pX_SetCoeff(shifter, j, tmp)
                        ZZ_negate(holder.x, holder.x)
                    ZZ_to_mpz(list_elt.value, &holder.x)
                    ans.append(list_elt)
                i += self.prime_pow.e
                if i >= prec:
                    break
                self.prime_pow.eis_shift(&shifter, &shifter, self.prime_pow.e, self.prime_pow.capdiv(prec - i))
            zerotest = 0
        while ans:
            if ans[-1] == zerotest:
                ans.pop()
            else:
                break
        while ans:
            if ans[0] == zerotest:
                ans.pop(0)
            else:
                break
        return ans

    def norm(self, base=None):
        r"""
        Return the absolute or relative norm of this element.

        .. NOTE::

            This is not the `p`-adic absolute value.  This is a
            field theoretic norm down to a ground ring.  If you want the
            `p`-adic absolute value, use the ``abs()`` function instead.

        If ``base`` is given then ``base`` must be a subfield of the
        parent `L` of ``self``, in which case the norm is the relative
        norm from L to ``base``.

        In all other cases, the norm is the absolute norm down to
        `\QQ_p` or `\ZZ_p`.

        EXAMPLES::

            sage: R = ZpCR(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: ((1+2*w)^5).norm()
            1 + 5^2 + O(5^5)
            sage: ((1+2*w)).norm()^5
            1 + 5^2 + O(5^5)

        TESTS::

           sage: R = ZpCA(5,5)
           sage: S.<x> = ZZ[]
           sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
           sage: W.<w> = R.ext(f)
           sage: ((1+2*w)^5).norm()                                                     # needs sage.geometry.polyhedron
           1 + 5^2 + O(5^5)
           sage: ((1+2*w)).norm()^5                                                     # needs sage.geometry.polyhedron
           1 + 5^2 + O(5^5)
           sage: R = ZpFM(5,5)
           sage: S.<x> = ZZ[]
           sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
           sage: W.<w> = R.ext(f)
           sage: ((1+2*w)^5).norm()                                                     # needs sage.geometry.polyhedron
           1 + 5^2
           sage: ((1+2*w)).norm()^5                                                     # needs sage.geometry.polyhedron
           1 + 5^2

        Check that :issue:`11586` has been resolved::

            sage: R.<x> = QQ[]
            sage: f = x^2 + 3*x + 1
            sage: M.<a> = Qp(7).extension(f)
            sage: M(7).norm()
            7^2 + O(7^22)
            sage: b = 7*a + 35
            sage: b.norm()
            4*7^2 + 7^3 + O(7^22)
            sage: b*b.frobenius()
            4*7^2 + 7^3 + O(7^22)
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
        if self.valuation() == 0:
            return self.parent().ground_ring()(self.matrix_mod_pn().det())
        else:
            if self.parent().e() == 1:
                norm_of_uniformizer = self.parent().ground_ring().uniformizer_pow(self.parent().degree())
            else:
                norm_of_uniformizer = (-1)**self.parent().degree() * self.parent().defining_polynomial()[0]
            return self.parent().ground_ring()(self.unit_part().matrix_mod_pn().det()) * norm_of_uniformizer**self.valuation()

    def trace(self, base = None):
        r"""
        Return the absolute or relative trace of this element.

        If ``base`` is given then ``base`` must be a subfield of the
        parent `L` of ``self``, in which case the norm is the relative
        norm from `L` to ``base``.

        In all other cases, the norm is the absolute norm down to
        `\QQ_p` or `\ZZ_p`.

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

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: R = ZpCA(5,5)
            sage: S.<x> = ZZ[]
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
            sage: R = ZpFM(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = (2+3*w)^7
            sage: b = (6+w^3)^5
            sage: a.trace()
            3*5 + 2*5^2 + 3*5^3 + 2*5^4
            sage: a.trace() + b.trace()
            4*5 + 5^2 + 5^3 + 2*5^4
            sage: (a+b).trace()
            4*5 + 5^2 + 5^3 + 2*5^4

        TESTS:

        We check that :issue:`32072` is resolved::

            sage: F = Qp(2)
            sage: S.<x> = F[]
            sage: L.<w> = F.ext(x^2 - 2)
            sage: L(0, 20).trace()
            O(2^10)
        """
        if base is not None:
            if base is self.parent():
                return self
            else:
                raise NotImplementedError
        if self._is_exact_zero():
            return self.parent().ground_ring()(0)
        elif self._is_inexact_zero():
            return self.parent().ground_ring()(0, (self.valuation() - 1) // self.parent().e() + 1)
        if self.valuation() >= 0:
            return self.parent().ground_ring()(self.matrix_mod_pn().trace())
        else:
            shift = -(self.valuation() // self.parent().e())
            return self.parent().ground_ring()((self * self.parent().prime() ** shift).matrix_mod_pn().trace()) / self.parent().prime()**shift

    def _rational_(self):
        """
        Return a rational approximation of ``self``.

        This does not try to optimize which rational is picked: see
        ``algdep`` for another option.

        EXAMPLES::

            sage: QQ(Qq(125,names='a')(-1/5))  # indirect doctest
            -1/5
        """
        if self.valuation() < 0:
            pk = self.parent().prime()**(-self.ordp()).ceil()
            return (self * pk)._integer_() / pk
        else:
            return QQ(self._integer_())

    def _prime_pow(self):
        """
        Provides access to ``self``'s ``prime_pow``.

        EXAMPLES::

            sage: R = ZpCR(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: w._prime_pow()
            PowComputer_ext for 5, with polynomial [3120 125 3110 75 0 1]
        """
        return self.prime_pow

    cdef int _pshift_self(self, long shift) except -1:
        """
        Multiply this element by ``p^shift``.

        TESTS:

        Check that :issue:`13647` has been fixed::

            sage: # needs sage.libs.flint
            sage: K = ZpCA(3)
            sage: R.<u> = K[]
            sage: L.<u> = K.extension(u^2 + 1)
            sage: L(R.gen())
            u + O(3^20)
            sage: K = ZpFM(3)
            sage: R.<u> = K[]
            sage: L.<u> = K.extension(u^2 + 1)
            sage: L(R.gen())
            u
        """
        if shift != 0:
            raise NotImplementedError


def _test_preprocess_list(R, L):
    r"""
    Given a list of elements convertible to ``ntl_ZZ_p``s, find the
    appropriate absolute precision and return a list of either
    ``ntl_ZZs`` or ``ntl_ZZ_ps``.

    INPUT:

    - ``R`` -- a `p`-adic extension ring

    - ``L`` -- list of rationals, integers, ints, longs,
      ``ntl_ZZ_ps``, ``ntl_ZZs``, ``IntegerMods`` or `p`-adic base
      elements

    OUTPUT:

    - ``LL`` -- if all inputs are integral, a list of ``ntl_ZZs``;
      otherwise, a list of ``ntl_ZZ_ps``, modulo `p^n` which is
      determined by the precision cap of ``R`` and the precisions of
      the elements in ``L``

    - ``min_val`` -- a valuation by which to multiply the elements of
      ``LL`` in order to recover the input elements of ``L``

    - ``ctx`` -- an ``ntl_ZZ_p_Context`` giving the power of `p`
      modulo which the elements in ``LL`` are defined; if ``None``,
      then the elements of ``LL`` are ``ntl_ZZs``.

    EXAMPLES::

        sage: from sage.rings.padics.padic_ZZ_pX_element import _test_preprocess_list
        sage: from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZ_p as ntl_ZZ_p
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [1,2,3])
        ([1, 2, 3], 0, None)
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [10,20,30])
        ([10, 20, 30], 0, None)
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [1/5,2/5,3])
        ([1, 2, 15], -1, NTL modulus 95367431640625)
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [1/5,mod(2,625),3])
        ([1, 10, 15], -1, NTL modulus 3125)
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [1/5,mod(2,625),ntl_ZZ_p(3,25)])
        ([1, 10, 15], -1, NTL modulus 125)
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [1/5,mod(2,625),Zp(5)(5,3)])
        ([1, 10, 25], -1, NTL modulus 625)
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [1/5,mod(2,625),Zp(5)(5,3),0])
        ([1, 10, 25, 0], -1, NTL modulus 625)
        sage: _test_preprocess_list(ZqCA(25,names='a',implementation="NTL"), [1/5,mod(2,625),Zp(5)(5,3),mod(0,3125)])
        ([1, 10, 25, 0], -1, NTL modulus 625)
        sage: x = polygen(ZZ, 'x')
        sage: T.<a> = Qp(5).extension(x^2 - 5)
        sage: _test_preprocess_list(T, [5^-1 + O(5)])
        ([1], -1, NTL modulus 25)
    """
    return preprocess_list(R(0), L)


cdef preprocess_list(pAdicZZpXElement elt, L):
    """
    See the documentation for :func:`_test_preprocess_list`.
    """
    cdef Py_ssize_t i
    cdef ZZ_c tmp
    cdef ntl_ZZ_pContext_class ctx
    cdef ntl_ZZ pshift_z
    cdef Integer pshift_m
    cdef ntl_ZZ py_tmp
    if not isinstance(L, list):
        raise TypeError("L must be a list")
    min_val, min_aprec, total_type = find_val_aprec(elt.prime_pow, L)
    #return "a","b","c"
    if total_type == two:
        # all integers
        return [ntl_ZZ(a) for a in L], zero, None
    if min_val < 0 and not elt.prime_pow.in_field:
        raise ValueError("negative valuation")
    if total_type == one:
        # rationals and integers
        py_tmp = ntl_ZZ.__new__(ntl_ZZ)
        py_tmp.x = elt.prime_pow.pow_ZZ_top()[0]
        ctx = ntl_ZZ_pContext(py_tmp)
    else:
        # integers, rationals and things with finite precision
        # note that min_val will be nonpositive since things with finite precision return nonpositive valuation from get_val_prec
        py_tmp = ntl_ZZ.__new__(ntl_ZZ)
        py_tmp.x = elt.prime_pow.pow_ZZ_tmp(mpz_get_ui((<Integer>(min_aprec - min_val)).value))[0]
        ctx = ntl_ZZ_pContext(py_tmp)
    if min_val < 0:
        pshift_z = ntl_ZZ.__new__(ntl_ZZ)
        pshift_z.x = elt.prime_pow.pow_ZZ_tmp(-mpz_get_si((<Integer>min_val).value))[0]
        pshift_m = elt.prime_pow.pow_Integer(-mpz_get_si((<Integer>min_val).value))
        for i from 0 <= i < len(L):
            if isinstance(L[i], ntl_ZZ):
                L[i] = ntl_ZZ_p(L[i]*pshift_z, ctx)
            elif isinstance(L[i], (Integer, Rational, int)):
                L[i] = ntl_ZZ_p(L[i]*pshift_m, ctx)
            elif isinstance(L[i], pAdicGenericElement) and L[i]._is_base_elt(elt.prime_pow.prime):
                L[i] = ntl_ZZ_p((L[i] >> min_val).lift(), ctx)
            elif isinstance(L[i], IntegerMod_abstract):
                L[i] = ntl_ZZ_p(L[i].lift()*pshift_m, ctx)
            elif (L[i].modulus_context() is not ctx) or min_val != zero:
                L[i] = ntl_ZZ_p(L[i].lift()*pshift_z, ctx)
    elif elt.parent().is_capped_relative() and min_val > 0:
        pshift_z = ntl_ZZ.__new__(ntl_ZZ)
        pshift_z.x = elt.prime_pow.pow_ZZ_tmp(mpz_get_ui((<Integer>min_val).value))[0]
        pshift_m = elt.prime_pow.pow_Integer(mpz_get_ui((<Integer>min_val).value))
        for i from 0 <= i < len(L):
            if isinstance(L[i], ntl_ZZ):
                ZZ_div(tmp, (<ntl_ZZ>L[i]).x, pshift_z.x)
                py_tmp = ntl_ZZ.__new__(ntl_ZZ)
                py_tmp.x = tmp
                L[i] = ntl_ZZ_p(py_tmp, ctx)
            elif isinstance(L[i], (Integer, Rational, int)):
                L[i] = ntl_ZZ_p(L[i]//pshift_m, ctx)
            elif isinstance(L[i], pAdicGenericElement) and L[i]._is_base_elt(elt.prime_pow.prime):
                L[i] = ntl_ZZ_p((L[i] >> min_val).lift(), ctx)
            elif isinstance(L[i], IntegerMod_abstract):
                L[i] = ntl_ZZ_p(L[i].lift()//pshift_m, ctx)
            elif (L[i].modulus_context() is not ctx) or min_val != zero:
                ZZ_div(tmp, (<ntl_ZZ>L[i].lift()).x, pshift_z.x)
                py_tmp = ntl_ZZ.__new__(ntl_ZZ)
                py_tmp.x = tmp
                L[i] = ntl_ZZ_p(py_tmp, ctx)
    else:
        for i from 0 <= i < len(L):
            if isinstance(L[i], (ntl_ZZ, Integer, Rational, int)):
                L[i] = ntl_ZZ_p(L[i], ctx)
            elif (isinstance(L[i], pAdicGenericElement) and L[i]._is_base_elt(elt.prime_pow.prime)) or isinstance(L[i], IntegerMod_abstract) or (L[i].modulus_context() is not ctx):
                L[i] = ntl_ZZ_p(L[i].lift(), ctx)
    return L, min_val, ctx


def _find_val_aprec_test(R, L):
    r"""
    Given a list ``L``, finds the minimum valuation, minimum absolute
    precision and minimum common type of the elements.

    INPUT:

    - ``R`` -- a `p`-adic extension
    - ``L`` -- list of integers, rationals, ``IntegerMods``, etc.

    OUTPUT: ``min_val`` -- the minimum valuation of any element in the list

    - ``min_aprec`` -- the minimum absolute precision of any element
      in the list; if infinite, a predefined constant ``big`` is
      returned instead

    - ``total_type`` --

      * if all elements are integers or ints: 2

      * if all elements are rationals or integers: 1

      * if some elements have finite precision: 0

    EXAMPLES::

        sage: from sage.rings.padics.padic_ZZ_pX_element import _find_val_aprec_test
        sage: from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZ_p as ntl_ZZ_p
        sage: _find_val_aprec_test(Zq(25,names='a',implementation="NTL"), [15, int(75), ntl_ZZ(625)])
        (1, 340282366920938463463374607431768211457, 2)
        sage: _find_val_aprec_test(Zq(25,names='a',implementation="NTL"), [5, int(25), 7/25])
        (-2, 340282366920938463463374607431768211457, 1)
        sage: _find_val_aprec_test(Zq(25,names='a',implementation="NTL"), [mod(4,125), Zp(5)(5,5), ntl_ZZ_p(16,625), 4/125])
        (-3, 3, 0)
        sage: _find_val_aprec_test(Zq(25,names='a',implementation="NTL"), [mod(25,125), Zp(5)(5,5), ntl_ZZ_p(15,625)])
        (0, 3, 0)
    """
    return find_val_aprec(R.prime_pow, L)


cdef find_val_aprec(PowComputer_ext pp, L):
    r"""
    Given a list ``L``, finds the minimum valuation, minimum absolute
    precision and minimum common type of the elements.

    INPUT:

    - ``pp`` -- a PowComputer_ext for the element that this list is
      being initialized into

    - ``L`` -- list of integers, rationals, ``IntegerMods``, etc.

    See the documentation for :func:`_find_val_aprec_test` for more details.
    """
    cdef Py_ssize_t i
    min_val = big
    min_aprec = big
    total_type = two # we begin by defaulting to the list elements being integers
    for i from 0 <= i < len(L):
        cur_val, cur_aprec, cur_type = get_val_prec(pp, L[i])
        #return "a","b","c"
        # proc_type == 0 indicates something with finite precision
        # proc_type == 1 indicates a rational, or something that cannot be coerced to an integer
        # proc_type == 2 indicates an integer, or something that can be coerced to an integer
        # but can be coerced to Z/p^n for any n.
        if cur_aprec < min_aprec:
            min_aprec = cur_aprec
        if cur_val < min_val:
            min_val = cur_val
        if cur_type < total_type:
            total_type = cur_type
    return min_val, min_aprec, total_type


def _test_get_val_prec(R, a):
    """
    Return valuation, absolute precision and type of an input
    element.

    INPUT:

    - ``R`` -- a `p`-adic extension ring to provide a ``PowComputer_ext``

    - ``a`` -- rational, integer, int, long, ``ntl_ZZ_p``,
      ``ntl_ZZ``, ``IntegerMod`` or `p`-adic base element

    OUTPUT:

    - ``val`` -- if ``a`` is exact, ``a.valuation(p)``, otherwise
      ``min(0, a.valuation())``

    - ``aprec`` -- the absolute precision of ``a``; if ``a`` is
      exact, a large predefined constant

    - type --

      * 2 - if ``a`` is an integer, int or long

      * 1 - if ``a`` is a rational

      * 0 - if ``a`` has finite precision

    EXAMPLES::

        sage: from sage.rings.padics.padic_ZZ_pX_element import _test_get_val_prec
        sage: from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZ_p as ntl_ZZ_p
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), 15)
        (1, 340282366920938463463374607431768211457, 2)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), ntl_ZZ(15))
        (1, 340282366920938463463374607431768211457, 2)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), int(15))
        (1, 340282366920938463463374607431768211457, 2)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), 1/15)
        (-1, 340282366920938463463374607431768211457, 1)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), Zp(5)(15,4))
        (0, 4, 0)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), Qp(5)(1/15,4))
        (-1, 4, 0)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), mod(15,625))
        (0, 4, 0)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), ntl_ZZ_p(15,625))
        (0, 4, 0)

    TESTS::

        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), 0)  # indirect doctest
        (340282366920938463463374607431768211457, 340282366920938463463374607431768211457, 2)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), ntl_ZZ(0))
        (340282366920938463463374607431768211457, 340282366920938463463374607431768211457, 2)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), int(0))
        (340282366920938463463374607431768211457, 340282366920938463463374607431768211457, 2)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), 0/1)
        (340282366920938463463374607431768211457, 340282366920938463463374607431768211457, 1)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), Zp(5)(25,4))
        (0, 4, 0)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), Qp(5)(1/25,4))
        (-2, 4, 0)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), Zp(5)(0))
        (340282366920938463463374607431768211457, 340282366920938463463374607431768211457, 1)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), mod(0,625))
        (0, 4, 0)
        sage: _test_get_val_prec(Zq(25,names='a',implementation="NTL"), ntl_ZZ_p(0,625))
        (0, 4, 0)
    """
    return get_val_prec(R.prime_pow, a)


cdef get_val_prec(PowComputer_ext pp, a):
    r"""
    Return valuation, absolute precision and type of an input element.

    INPUT:

    - ``pp`` -- a ``PowComputer_ext``

    - ``a`` -- rational, integer, int, long, ``ntl_ZZ_p``,
      ``ntl_ZZ``, ``IntegerMod`` or `p`-adic base element

    See :func:`_test_get_val_prec` for more details.
    """
    cdef ntl_ZZ py_tmp
    if isinstance(a, Integer):
        if a == 0:
            return (big, big, two)
        return (a.valuation(pp.prime), big, two)
    if isinstance(a, ntl_ZZ):
        if ZZ_IsZero((<ntl_ZZ>a).x):
            return (big, big, two)
        py_tmp = ntl_ZZ.__new__(ntl_ZZ)
        py_tmp.x = pp.pow_ZZ_tmp(1)[0]
        return (Integer(a.valuation(py_tmp)), big, two)
    if isinstance(a, int):
        if a == 0:
            return (big, big, two)
        return (Integer(a).valuation(pp.prime), big, two)
    if isinstance(a, Rational):
        if a == 0:
            return (big, big, one)
        val = a.valuation(pp.prime)
        return (val, big, one)
    if isinstance(a, pAdicGenericElement) and a._is_base_elt(pp.prime):
        if a.parent().prime() == pp.prime:
            if a._is_exact_zero():
                return (big, big, one)
            val = a.valuation()
            return (val if val < zero else zero, a.precision_absolute(), zero)
        else:
            raise TypeError("primes must match")
    cdef mpz_t leftover
    cdef long long_val
    cdef Integer Integer_val
    if isinstance(a, IntegerMod_abstract):
        mpz_init(leftover)
        long_val = mpz_remove(leftover, (<Integer>a.modulus()).value, pp.prime.value)
        if long_val > 0 and mpz_cmp_ui(leftover, 1) == 0:
            mpz_clear(leftover)
            Integer_val = PY_NEW(Integer)
            mpz_set_ui(Integer_val.value, long_val)
            # Since we're guaranteed to be in type 0, we don't care about computing the actual valuation
            return (zero, Integer_val, zero)
        else:
            mpz_clear(leftover)
            raise TypeError("modulus must be a positive power of the appropriate prime")
    cdef ZZ_c leftover_z
    if isinstance(a, ntl_ZZ_p):
        long_val = ZZ_remove(leftover_z, (<ntl_ZZ_p>a).c.p.x, pp.pow_ZZ_tmp(1)[0])
        if long_val > 0 and ZZ_IsOne(leftover_z):
            Integer_val = PY_NEW(Integer)
            mpz_set_ui(Integer_val.value, long_val)
            # Since we're guaranteed to be in type 0, we don't care about computing the actual valuation
            return (zero, Integer_val, zero)
        else:
            print(long_val)
            py_tmp = ntl_ZZ.__new__(ntl_ZZ)
            py_tmp.x = (<ntl_ZZ_p>a).c.p.x
            print(py_tmp)
            py_tmp.x = leftover_z
            print(py_tmp)
            raise TypeError("modulus must be a positive power of the appropriate prime")
    raise TypeError("unsupported type for list element: %s" % type(a))
