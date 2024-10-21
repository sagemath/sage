# sage_setup: distribution = sagemath-ntl
# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
# sage.doctest: needs sage.rings.padics
r"""
`p`-adic Extension Element

A common superclass for all elements of extension rings and field of `\ZZ_p` and
`\QQ_p`.

AUTHORS:

- David Roe (2007): initial version

- Julian Rueth (2012-10-18): added residue
"""
# ****************************************************************************
#       Copyright (C) 2007-2010 David Roe <roed.math@gmail.com>
#                     2012 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p

cdef class pAdicExtElement(pAdicGenericElement):
    cdef int _set_from_list(self, L) except -1:
        """
        Set ``self`` from a list.

        The list should either be uniform in type, or all of the entries
        should be coercible to integers. If any of the entries in ``L``
        is a list, ``L`` will be cast to a ZZ_pEX.

        INPUT:

        - ``L`` -- list
        """
        raise NotImplementedError

    cdef int _set_from_list_rel(self, L, long relprec) except -1:
        raise NotImplementedError

    cdef int _set_from_list_abs(self, L, long absprec) except -1:
        raise NotImplementedError

    cdef int _set_from_list_both(self, L, long absprec, long relprec) except -1:
        raise NotImplementedError

    cdef int _set_from_ZZX(self, ZZX_c poly) except -1:
        """
        Set from a ZZX_c, choosing how to handle based on the
        precision type of ``self.parent()``.

        Fixed modulus elements should override this function.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZX_rel(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZX_abs(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError("_set_from_ZZX should have been overridden")

    cdef int _set_from_ZZX_rel(self, ZZX_c poly, long relprec) except -1:
        """
        Set from a ZZX_c with bounded relative precision.

        Capped relative elements should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZX_both(poly, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZX_abs(self, ZZX_c poly, long absprec) except -1:
        """
        Set from a ZZX_c with bounded absolute precision.

        Capped absolute elements should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZX_both(poly, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, long relprec) except -1:
        """
        Set from a ZZX_c with both absolute and relative precisions bounded.

        This function should be overridden for both capped absolute
        and capped relative elements.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZX(poly)
        else:
            raise RuntimeError("_set_from_ZZX_both should have been overridden")

    cdef int _set_from_ZZ_pX(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx) except -1:
        """
        Set ``self`` from a ZZ_pX defined with context ctx.

        This function should be overridden for fixed modulus elements.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pX_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pX_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError("_set_from_ZZ_pX should have been overridden")

    cdef int _set_from_ZZ_pX_rel(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pX_c with bounded relative precision.

        Capped relative rings should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZ_pX_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pX_abs(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pX_c with bounded absolute precision.

        Capped absolute rings should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZ_pX_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pX_both(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec, long relprec) except -1:
        """
        Set from a ZZ_pX_c with both absolute and relative precision bounded.

        This function should be overridden by both capped absolute and capped relative elements.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pX(poly, ctx)
        else:
            raise RuntimeError("_set_from_ZZ_pX_both should have been overridden")

    cdef int _set_from_ZZ_pE(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx) except -1:
        """
        Set from a ZZ_pE_c.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pE_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pE_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError("_set_from_ZZ_pE should have been overridden")

    cdef int _set_from_ZZ_pE_rel(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pE_c with bounded relative precision.

        Capped relative rings should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZ_pE_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pE_abs(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pE_c with bounded absolute precision.

        Capped absolute elements should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZ_pE_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pE_both(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1:
        """
        Set from a ZZ_pE_c with both absolute and relative precision bounded.

        Capped absolute and capped relative elements should override
        this function.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pE(poly, ctx)
        else:
            raise RuntimeError("_set_from_ZZ_pE_both should have been overridden")

    cdef int _set_from_ZZ_pEX(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx) except -1:
        """
        Set ``self`` from a ZZ_pEX_c.

        Fixed modulus elements should override this function.

        This function is not used internally.
        """
        if self.parent().is_capped_relative():
            self._set_from_ZZ_pEX_rel(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        elif self.parent().is_capped_absolute():
            self._set_from_ZZ_pEX_abs(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap)
        else:
            raise RuntimeError("_set_from_ZZ_pEX should have been overridden")

    cdef int _set_from_ZZ_pEX_rel(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1:
        """
        Set from a ZZ_pEX_c with bounded relative precision.

        Capped relative elements should override this function, so the
        default implementation is for capped absolute.

        This function is not used internally.
        """
        self._set_from_ZZ_pEX_both(poly, ctx, (<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef int _set_from_ZZ_pEX_abs(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1:
        """
        Set from a ZZ_pEX_c with bounded absolute precision.

        Capped absolute elements should override this function, so the
        default implementation is for capped relative.

        This function is not used internally.
        """
        self._set_from_ZZ_pEX_both(poly, ctx, absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef int _set_from_ZZ_pEX_both(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1:
        """
        Set from a ZZ_pEX_c with both absolute and relative precision bounded.

        Capped absolute and capped relative elements should override
        this function.

        This function is not used internally.
        """
        if self.parent().is_fixed_mod():
            self._set_from_ZZ_pEX(poly, ctx)
        else:
            raise RuntimeError("_set_from_ZZ_pEX_both should have been overridden")

    cdef long _check_ZZ_pContext(self, ntl_ZZ_pContext_class ctx) except -1:
        raise NotImplementedError

    cdef long _check_ZZ_pEContext(self, ntl_ZZ_pEContext_class ctx) except -1:
        raise NotImplementedError

    cdef ext_p_list(self, bint pos):
        raise NotImplementedError

    cdef ext_p_list_precs(self, bint pos, long prec):
        raise NotImplementedError

    def _const_term_test(self):
        """
        Return the constant term of a polynomial representing ``self``.

        This function is mainly for troubleshooting, and the meaning
        of the return value will depend on whether ``self`` is capped
        relative or otherwise.

        EXAMPLES::

            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 + 125*x - 5
            sage: W.<w> = R.ext(f)
            sage: a = W(566)
            sage: a._const_term_test()
            566
        """
        cdef ntl_ZZ_p ans = ntl_ZZ_p(modulus=self.parent().prime())
        ans.x = self._const_term()
        return ans

    cdef ZZ_p_c _const_term(self) noexcept:
        raise NotImplementedError

    def _ext_p_list(self, pos):
        r"""
        Return a list of integers (in the Eisenstein case) or a list
        of lists of integers (in the unramified case).

        ``self`` can be reconstructed as a sum of elements of the list times
        powers of the uniformiser (in the Eisenstein case), or as a sum of
        powers of the p times polynomials in the generator (in the
        unramified case).

        Note that zeros are truncated from the returned list, so you must
        use the :func:`valuation()` function to completely recover ``self``.

        INPUT:

        - ``pos`` -- boolean; if ``True``, all integers will be in the range
          `[0,p-1]`, otherwise they will be in the range `[(1-p)/2, p/2]`

        OUTPUT:

        - a list of integers or list of lists giving the
          series expansion of ``self``

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
        """
        return self.ext_p_list(pos)

    def frobenius(self, arithmetic=True):
        r"""
        Return the image of this element under the Frobenius automorphism
        applied to its parent.

        INPUT:

        - ``arithmetic`` -- whether to apply the arithmetic Frobenius (acting
          by raising to the `p`-th power on the residue field). If ``False`` is
          provided, the image of geometric Frobenius (raising to the `(1/p)`-th
          power on the residue field) will be returned instead.

        EXAMPLES::

            sage: R.<a> = Zq(5^4,3)
            sage: a.frobenius()
            (a^3 + a^2 + 3*a) + (3*a + 1)*5 + (2*a^3 + 2*a^2 + 2*a)*5^2 + O(5^3)
            sage: f = R.defining_polynomial()
            sage: f(a)
            O(5^3)
            sage: f(a.frobenius())
            O(5^3)
            sage: for i in range(4): a = a.frobenius()
            sage: a
            a + O(5^3)

            sage: K.<a> = Qq(7^3,4)
            sage: b = (a+1)/7
            sage: c = b.frobenius(); c
            (3*a^2 + 5*a + 1)*7^-1 + (6*a^2 + 6*a + 6) + (4*a^2 + 3*a + 4)*7 + (6*a^2 + a + 6)*7^2 + O(7^3)
            sage: c.frobenius().frobenius()
            (a + 1)*7^-1 + O(7^3)

        An error will be raised if the parent of ``self`` is a ramified extension::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = Qp(5).extension(x^2 - 5)
            sage: a.frobenius()
            Traceback (most recent call last):
            ...
            NotImplementedError: Frobenius automorphism only implemented for unramified extensions
        """
        R = self.parent()
        if R.absolute_e() != 1:
            raise NotImplementedError("Frobenius automorphism only implemented for unramified extensions")
        if self.is_zero(): return self
        L = self.teichmuller_expansion()
        ppow = R.uniformizer_pow(self.valuation())
        if arithmetic:
            exp = R.prime()
        else:
            exp = R.prime()**(R.absolute_degree()-1)
        ans = ppow * L[0]**exp
        for m in range(1,len(L)):
            ppow = ppow << 1
            ans += ppow * L[m]**exp
        return ans

    cpdef bint _is_base_elt(self, p) except -1:
        r"""
        Return ``True`` if this element is an element of Zp or Qp (rather than
        an extension).

        INPUT:

        - ``p`` -- a prime, which is compared with the parent of this element

        EXAMPLES::

            sage: K.<a> = Qq(7^3,4)
            sage: a._is_base_elt(5)
            False
        """
        return False

    def residue(self, absprec=1, field=None, check_prec=True):
        r"""
        Reduces this element modulo `\pi^\mathrm{absprec}`.

        INPUT:

        - ``absprec`` -- nonnegative integer (default: 1)

        - ``field`` -- boolean (default: ``None``); for precision 1, whether to return
          an element of the residue field or a residue ring.  Currently unused.

        - ``check_prec`` -- boolean (default: ``True``); whether to raise an error if this
          element has insufficient precision to determine the reduction.  Errors are never
          raised for fixed-mod or floating-point types.

        OUTPUT: this element reduced modulo `\pi^\mathrm{absprec}`

        If ``absprec`` is zero, then as an element of `\ZZ/(1)`.

        If ``absprec`` is one, then as an element of the residue field.

        .. NOTE::

            Only implemented for ``absprec`` less than or equal to one.

        AUTHORS:

        - Julian Rueth (2012-10-18): initial version

        EXAMPLES:

        Unramified case::

            sage: # needs sage.libs.flint
            sage: R = ZpCA(3,5)
            sage: S.<a> = R[]
            sage: W.<a> = R.extension(a^2 + 9*a + 1)
            sage: (a + 1).residue(1)
            a0 + 1
            sage: a.residue(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: reduction modulo p^n with n>1

        Eisenstein case::

            sage: R = ZpCA(3,5)
            sage: S.<a> = R[]
            sage: W.<a> = R.extension(a^2 + 9*a + 3)
            sage: (a + 1).residue(1)
            1
            sage: a.residue(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: residue() not implemented in extensions for absprec larger than one

        TESTS::

            sage: # needs sage.libs.flint
            sage: K = Qp(3,5)
            sage: S.<a> = R[]
            sage: W.<a> = R.extension(a^2 + 9*a + 1)
            sage: (a/3).residue(0)
            Traceback (most recent call last):
            ...
            ValueError: element must have nonnegative valuation in order to compute residue

            sage: # needs sage.libs.flint
            sage: R = ZpFM(3,5)
            sage: S.<a> = R[]
            sage: W.<a> = R.extension(a^2 + 3)
            sage: W.one().residue(0)
            0
            sage: a.residue(-1)
            Traceback (most recent call last):
            ...
            ValueError: cannot reduce modulo a negative power of the uniformizer
            sage: a.residue(16)
            Traceback (most recent call last):
            ...
            NotImplementedError: residue() not implemented in extensions for absprec larger than one
        """
        if absprec < 0:
            raise ValueError("cannot reduce modulo a negative power of the uniformizer")
        if self.valuation() < 0:
            raise ValueError("element must have nonnegative valuation in order to compute residue")
        R = self.parent()
        if check_prec and (R.is_fixed_mod() or R.is_floating_point()):
            check_prec = False
        if check_prec and absprec > self.precision_absolute():
            from precision_error import PrecisionError
            raise PrecisionError("not enough precision known in order to compute residue")
        if field and absprec != 1:
            raise ValueError("field keyword may only be set at precision 1")

        if absprec == 0:
            from sage.rings.finite_rings.integer_mod import Mod
            return Mod(0,1)
        elif absprec == 1:
            return R.residue_field()(self.expansion(0))
        else:
            raise NotImplementedError("residue() not implemented in extensions for absprec larger than one")
