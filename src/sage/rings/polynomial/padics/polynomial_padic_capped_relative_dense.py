"""
p-adic Capped Relative Dense Polynomials
"""

#*****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.rings.polynomial.polynomial_element_generic
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
import sage.rings.polynomial.polynomial_integer_dense_ntl
import sage.rings.integer
import sage.rings.integer_ring
import sage.rings.padics.misc as misc
import sage.rings.padics.precision_error as precision_error
import sage.rings.fraction_field_element as fraction_field_element
import copy

from sage.libs.all import pari, pari_gen
from sage.libs.ntl.all import ZZX
from sage.rings.infinity import infinity

min = misc.min
ZZ = sage.rings.integer_ring.ZZ
PrecisionError = precision_error.PrecisionError
Integer = sage.rings.integer.Integer
Polynomial_integer_dense = sage.rings.polynomial.polynomial_integer_dense_ntl.Polynomial_integer_dense_ntl
Polynomial_generic_cdv = sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_cdv


class Polynomial_padic_capped_relative_dense(Polynomial_generic_cdv, Polynomial_padic):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct = False, absprec = infinity, relprec = infinity):
        """
        TESTS::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: R([K(13), K(1)])
            (1 + O(13^7))*t + 13 + O(13^8)
            sage: T.<t> = ZZ[]
            sage: R(t + 2)
            (1 + O(13^7))*t + 2 + O(13^7)

        Check that :trac:`13620` has been fixed::

            sage: f = R.zero()
            sage: R(f.dict())
            0

        Check that :trac:`29829` has been fixed::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = x + 5
            sage: S.<y> = PolynomialRing(Qp(5))
            sage: g2 = S(f)
            sage: 25*g2
            (5^2 + O(5^22))*y + 5^3 + O(5^23)
        """
        Polynomial.__init__(self, parent, is_gen=is_gen)
        self._polygon = None
        parentbr = parent.base_ring()
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if construct:
            (self._poly, self._valbase, self._relprecs, self._normalized, self._valaddeds, self._list) = x #the last two of these may be None
            return
        elif is_gen:
            self._poly = PolynomialRing(ZZ, parent.variable_name()).gen()
            self._valbase = 0
            self._valaddeds = [infinity, 0]
            self._relprecs = [infinity, parentbr.precision_cap()]
            self._normalized = True
            self._list = None
            return

        #First we list the types that are turned into Polynomials
        if isinstance(x, ZZX):
            x = Polynomial_integer_dense(PolynomialRing(ZZ, parent.variable_name()), x, construct = True)
        elif isinstance(x, fraction_field_element.FractionFieldElement) and \
               x.denominator() == 1:
            #Currently we ignore precision information in the denominator.  This should be changed eventually
            x = x.numerator()

        #We now coerce various types into lists of coefficients.  There are fast pathways for some types of polynomials
        if isinstance(x, Polynomial):
            if x.parent() is self.parent():
                if absprec is not infinity or relprec is not infinity:
                    x._normalize()
                self._poly = x._poly
                self._valbase = x._valbase
                self._valaddeds = x._valaddeds
                self._relprecs = x._relprecs
                self._normalized = x._normalized
                self._list = x._list
                if absprec is not infinity or relprec is not infinity:
                    self._adjust_prec_info(absprec, relprec)
                return
            elif x.base_ring() is ZZ:
                self._poly = PolynomialRing(ZZ, parent.variable_name())(x)
                self._valbase = Integer(0)
                p = parentbr.prime()
                self._relprecs = [c.valuation(p) + parentbr.precision_cap() for c in x.list()]
                self._comp_valaddeds()
                self._normalized = len(self._valaddeds) == 0 or (min(self._valaddeds) == 0)
                self._list = None
                if absprec is not infinity or relprec is not infinity:
                    self._adjust_prec_info(absprec, relprec)
                return
            else:
                x = [parentbr(a) for a in x.list()]
                check = False
        elif isinstance(x, dict):
            zero = parentbr.zero()
            n = max(x.keys()) if x else 0
            v = [zero] * (n + 1)
            for i, z in x.items():
                v[i] = z
            x = v
        elif isinstance(x, pari_gen):
            x = [parentbr(w) for w in x.list()]
            check = False
        # The default behavior, if we haven't already figured out what
        # the type is, is to assume it coerces into the base_ring as a
        # constant polynomial
        elif not isinstance(x, list):
            x = [x] # constant polynomial

        # In contrast to other polynomials, the zero element is not distinguished
        # by having its list empty. Instead, it has list [0]
        if not x:
            x = [parentbr.zero()]
        if check:
            x = [parentbr(z) for z in x]

        # Remove this -- for p-adics this is terrible, since it kills any non exact zero.
        #if len(x) == 1 and not x[0]:
        #    x = []

        self._list = x
        self._valaddeds = [a.valuation() for a in x]
        self._valbase = sage.rings.padics.misc.min(self._valaddeds)
        if self._valbase is infinity:
            self._valaddeds = []
            self._relprecs = []
            self._poly = PolynomialRing(ZZ, parent.variable_name())()
            self._normalized = True
            if absprec is not infinity or relprec is not infinity:
                self._adjust_prec_info(absprec, relprec)
        else:
            self._valaddeds = [c - self._valbase for c in self._valaddeds]
            self._relprecs = [a.precision_absolute() - self._valbase for a in x]
            self._poly = PolynomialRing(ZZ, parent.variable_name())([a >> self._valbase for a in x])
            self._normalized = True
            if absprec is not infinity or relprec is not infinity:
                self._adjust_prec_info(absprec, relprec)

    def _new_constant_poly(self, a, P):
        """
        Create a new constant polynomial in parent P with value a.

        ASSUMPTION:

        The value a must be an element of the base ring of P. That
        assumption is not verified.

        EXAMPLES::

            sage: R.<t> = Zp(5)[]
            sage: t._new_constant_poly(O(5),R)
            O(5)
        """
        return self.__class__(P, [a], check=False)

    def _normalize(self):
        # Currently slow: need to optimize
        if not self._normalized:
            if self._valaddeds is None:
                self._comp_valaddeds()
            val = sage.rings.padics.misc.min(self._valaddeds)
            prime_pow = self.base_ring().prime_pow
            selflist = self._poly.list()
            if val is infinity:
                pass
            elif val != 0:
                self._relprecs = [max(prec - val,0) for prec in self._relprecs]
                v = [Integer(0) if (e is infinity) else ((c // prime_pow(val)) % prime_pow(e)) for (c,e) in zip(selflist, self._relprecs)]
                self._poly = self._poly.parent()(v, check=False)
                self._valbase += val
                self._valaddeds = [c - val for c in self._valaddeds]
            else:
                self._poly = self._poly.parent()([Integer(0) if (e is infinity) else (c % prime_pow(e)) for (c,e) in zip(selflist, self._relprecs)], check=False)
            self._normalized = True

    def _reduce_poly(self):
        selflist = self._poly.list()
        prime_pow = self.base_ring().prime_pow
        self._poly = self._poly.parent()([Integer(0) if (e is infinity) else (c % prime_pow(e)) for (c, e) in zip(selflist, self._relprecs)], check=False)

    def __reduce__(self):
        """
        For pickling.  This function is here because the relative precisions were getting screwed up for some reason.
        """
        return make_padic_poly, (self.parent(), (self._poly, self._valbase, self._relprecs, self._normalized, self._valaddeds, self._list), 0)

    def _comp_list(self):
        """
        Recomputes the list of coefficients.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = t[:1]
            sage: a._comp_list()
            sage: a
            0
        """
        if self.degree() == -1 and self._valbase == infinity:
            self._list = []
        polylist = self._poly.list()
        polylen = len(polylist)
        self._list = [self.base_ring()(polylist[i], absprec=self._relprecs[i]) << self._valbase for i in range(polylen)] \
                     + [self.base_ring()(0, absprec=self._relprecs[i] + self._valbase) for i in range(polylen, len(self._relprecs))]
        while self._list and self._list[-1]._is_exact_zero():
            self._list.pop()

    def _comp_valaddeds(self):
        self._valaddeds = []
        for i in range(self._poly.degree() + 1):
            tmp = self._poly.list()[i].valuation(self.parent().base_ring().prime())
            if tmp is infinity or tmp > self._relprecs[i]:
                self._valaddeds.append(self._relprecs[i])
            else:
                self._valaddeds.append(tmp)
        for i in range(self._poly.degree() + 1, len(self._relprecs)):
            self._valaddeds.append(self._relprecs[i])

    def _adjust_prec_info(self, absprec=infinity, relprec=infinity):
        r"""
        Assumes that self._poly, self._val and self._relprec are set initially and adjusts self._val and self._relprec to the termwise minimum of absprec and relprec.
        """
        return

#         min = sage.rings.padics.misc.min
#         slen = len(self._relprec)
#         if isinstance(absprec, list):
#             alen = len(absprec)
#         elif absprec is infinity:
#             alen = 0
#             absprec = []
#         else:
#             alen = 1
#         if isinstance(relprec, list):
#             rlen = len(relprec)
#         elif relprec is infinity:
#             rlen = 0
#             relprec = []
#         else:
#             rlen = 1
#         preclen = max(slen, rlen, alen)
#         if not isinstance(absprec, list):
#             absprec = [absprec] * preclen
#         if not isinstance(relprec, list):
#             relprec = [relprec] * preclen
#         vallist = [c.valuation(self.base_ring().prime()) + self._val for c in self._poly.list()] #######
#         vmin = min(vallist)
#         amin = min(absprec)
#         if amin < vmin:
#             vmin = amin
#         if vmin < self._val:
#             vadjust =

#         if not isinstance(absprec, list):
#             self._val = min(vallist + [absprec])
#             absprec = [absprec] * preclen
#         else:
#             self._val = padics.misc.min(vallist + absprec)
#             absprec = absprec + [infinity] * (preclen - len(absprec))
#         if self._val is infinity:
#             self._relprec = []
#             return
#         if not isinstance(relprec, list):
#             relprec = [relprec] * preclen
#         else:
#             relprec = relprec + [parent.base_ring().precision_cap()] * (preclen - len(relprec))
#         self._relprec = [min(a, v + r) - self._val for (a, r, v) in zip(absprec, relprec, vallist)]
#Remember to normalize at the end if self._normalized is true because you need to reduce mod p^n

    def _getprecpoly(self, n):
        one = Integer(1)
        return self._poly.parent()([(0 if (c is infinity) else (one << (n * c))) for c in self._relprecs])

    def _getvalpoly(self, n):
        one = Integer(1)
        if self._valaddeds is None:
            self._comp_valaddeds()
        return self._poly.parent()([(0 if (c is infinity) else (one << (n * c))) for c in self._valaddeds] + \
                                   [(0 if (c is infinity) else (one << (n * c))) for c in self._relprecs[len(self._valaddeds):]])

    def list(self, copy=True):
        """
        Return a list of coefficients of ``self``.

        .. NOTE::

            The length of the list returned may be greater
            than expected since it includes any leading zeros
            that have finite absolute precision.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = 2*t^3 + 169*t - 1
            sage: a
            (2 + O(13^7))*t^3 + (13^2 + O(13^9))*t + 12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + O(13^7)
            sage: a.list()
            [12 + 12*13 + 12*13^2 + 12*13^3 + 12*13^4 + 12*13^5 + 12*13^6 + O(13^7),
             13^2 + O(13^9),
             0,
             2 + O(13^7)]
        """
        if self._list is None:
            self._comp_list()
        if copy:
            return list(self._list)
        else:
            return self._list

    def lift(self):
        """
        Return an integer polynomial congruent to this one modulo the
        precision of each coefficient.

        .. NOTE::

            The lift that is returned will not necessarily be the same
            for polynomials with the same coefficients (i.e. same values
            and precisions): it will depend on how the polynomials are
            created.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = 13^7*t^3 + K(169,4)*t - 13^4
            sage: a.lift()
            62748517*t^3 + 169*t - 28561
        """
        return self.base_ring().prime_pow(self._valbase) * self._poly

    def __getitem__(self, n):
        """
        Returns the coefficient of x^n if `n` is an integer,
        returns the monomials of self of degree in slice `n` if `n` is a slice.

        Return the `n`-th coefficient of ``self``.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = 13^7*t^3 + K(169,4)*t - 13^4
            sage: a[1]
            13^2 + O(13^4)

        Slices can be used to truncate polynomials::

            sage: a[:2]
            (13^2 + O(13^4))*t + 12*13^4 + 12*13^5 + 12*13^6 + 12*13^7 + 12*13^8 + 12*13^9 + 12*13^10 + O(13^11)

        Any other kind of slicing is deprecated or an error, see
        :trac:`18940`::

            sage: a[1:3]
            doctest:warning...:
            DeprecationWarning: polynomial slicing with a start index is deprecated, use list() and slice the resulting list instead
            See http://trac.sagemath.org/18940 for details.
            0*t^2 + (13^2 + O(13^4))*t
            sage: a[1:3:2]
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomial slicing with a step is not defined
        """
        d = len(self._relprecs)  # = degree + 1
        if isinstance(n, slice):
            start, stop, step = n.start, n.stop, n.step
            if step is not None:
                raise NotImplementedError("polynomial slicing with a step is not defined")
            if start is None:
                start = 0
            else:
                if start < 0:
                    start = 0
                from sage.misc.superseded import deprecation
                deprecation(18940, "polynomial slicing with a start index is deprecated, use list() and slice the resulting list instead")
            if stop is None or stop > d:
                stop = d
            values = ([self.base_ring().zero()] * start
                      + [self[i] for i in range(start, stop)])
            return self.parent()(values)

        try:
            n = n.__index__()
        except AttributeError:
            raise TypeError("list indices must be integers, not {0}".format(type(n).__name__))

        if n < 0 or n >= d:
            return self.base_ring().zero()
        if self._list is not None:
            return self._list[n]
        return self.base_ring()(self.base_ring().prime_pow(self._valbase)
            * self._poly[n], absprec = self._valbase + self._relprecs[n])

    def _add_(self, right):
        """
        Return the sum of ``self`` and ``right``.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = t^4 + 17*t^2 + 1
            sage: b = -t^4 + 9*t^2 + 13*t - 1
            sage: c = a + b; c
            O(13^7)*t^4 + (2*13 + O(13^7))*t^2 + (13 + O(13^8))*t + O(13^7)
            sage: c.list()
            [O(13^7), 13 + O(13^8), 2*13 + O(13^7), 0, O(13^7)]
        """
        selfpoly = self._poly
        rightpoly = right._poly
        if self._valbase > right._valbase:
            selfpoly = selfpoly * self.base_ring().prime_pow(self._valbase - right._valbase)
            baseval = right._valbase
        elif self._valbase < right._valbase:
            rightpoly = rightpoly * self.base_ring().prime_pow(right._valbase - self._valbase)
            baseval = self._valbase
        else:
            baseval = self._valbase
        # Currently we don't reduce the coefficients of the answer modulo the appropriate power of p or normalize
        return Polynomial_padic_capped_relative_dense(self.parent(), \
                                                      (selfpoly + rightpoly, \
                                                       baseval, \
                                                       [min(a + self._valbase - baseval, b + right._valbase - baseval) for (a, b) in
                                                              zip(_extend_by_infinity(self._relprecs, max(len(self._relprecs), len(right._relprecs))), \
                                                                  _extend_by_infinity(right._relprecs, max(len(self._relprecs), len(right._relprecs))))], \
                                                       False, None, None), construct = True)

    def _sub_(self, right):
        """
        Return the difference of ``self`` and ``right``.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = t^4 + 17*t^2 + 1
            sage: b = t^4 - 9*t^2 - 13*t + 1
            sage: c = a - b; c
            O(13^7)*t^4 + (2*13 + O(13^7))*t^2 + (13 + O(13^8))*t + O(13^7)
            sage: c.list()
            [O(13^7), 13 + O(13^8), 2*13 + O(13^7), 0, O(13^7)]
        """
        selfpoly = self._poly
        rightpoly = right._poly
        if self._valbase > right._valbase:
            selfpoly = selfpoly * self.base_ring().prime_pow(self._valbase - right._valbase)
            baseval = right._valbase
        elif self._valbase < right._valbase:
            rightpoly = rightpoly * self.base_ring().prime_pow(right._valbase - self._valbase)
            baseval = self._valbase
        else:
            baseval = self._valbase
        # Currently we don't reduce the coefficients of the answer modulo the appropriate power of p or normalize
        return Polynomial_padic_capped_relative_dense(self.parent(), \
                                                      (selfpoly - rightpoly, \
                                                       baseval, \
                                                       [min(a + self._valbase - baseval, b + right._valbase - baseval) for (a, b) in
                                                              zip(_extend_by_infinity(self._relprecs, max(len(self._relprecs), len(right._relprecs))), \
                                                                  _extend_by_infinity(right._relprecs, max(len(self._relprecs), len(right._relprecs))))], \
                                                       False, None, None), construct = True)

    def _mul_(self, right):
        r"""
        Multiplies ``self`` and ``right``.

        ALGORITHM: We use an algorithm thought up by Joe Wetherell to
        find the precisions of the product.  It works as follows:
        Suppose $f = \sum_i a_i x^i$ and $g = \sum_j b_j x^j$. Let $N
        = \max(\deg f, \deg g) + 1$ (in the actual implementation we
        use $N = 2^{\lfloor \log_2\max(\deg f, \deg g)\rfloor + 1}$).
        The valuations and absolute precisions of each coefficient
        contribute to the absolute precision of the kth coefficient of
        the product in the following way: for each $i + j = k$, you
        take the valuation of $a_i$ plus the absolute precision of
        $b_j$, and then take the valuation of $b_j$ plus the absolute
        precision of $a_i$, take the minimum of those two, and then
        take the minimum over all $i$, $j$ summing to $k$.

        You can simulate this as follows. Construct new polynomials of
        degree $N$:

        \begin{align*}
        A &= \sum_i N^{\mbox{valuation of $a_i$}} x^i \\
        B &= \sum_j N^{\mbox{absolute precision of $b_j$}} x^j \\
        C &= \sum_i N^{\mbox{absolute precision of $a_i$}} x^i \\
        D &= \sum_j N^{\mbox{valuation of $b_j$}} x^j \\
        \end{align*}

        Now you compute AB and CD. Because you're representing things
        'N-adically', you don't get any 'overflow', and you can just
        read off what the precisions of the product are. In fact it
        tells you more, it tells you exactly how many terms of each
        combination of valuation modulus contribute to each term of
        the product (though this feature is not currently exposed in
        our implementation.

        Since we're working 'N-adically' we can just consider
        $N^{\infty} = 0$.

        NOTE: The timing of normalization in arithmetic operations
        may very well change as we do more tests on the relative time
        requirements of these operations.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = t^4 + 17*t^2 + 1
            sage: b = -t^4 + 9*t^2 + 13*t - 1
            sage: c = a + b; c
            O(13^7)*t^4 + (2*13 + O(13^7))*t^2 + (13 + O(13^8))*t + O(13^7)
            sage: d = R([K(1,4), K(2, 6), K(1, 5)]); d
            (1 + O(13^5))*t^2 + (2 + O(13^6))*t + 1 + O(13^4)
            sage: e = c * d; e
            O(13^7)*t^6 + O(13^7)*t^5 + (2*13 + O(13^6))*t^4 + (5*13 + O(13^6))*t^3 + (4*13 + O(13^5))*t^2 + (13 + O(13^5))*t + O(13^7)
            sage: e.list()
            [O(13^7),
             13 + O(13^5),
             4*13 + O(13^5),
             5*13 + O(13^6),
             2*13 + O(13^6),
             O(13^7),
             O(13^7)]
        """
        self._normalize()
        right._normalize()
        zzpoly = self._poly * right._poly
        if len(self._relprecs) == 0 or len(right._relprecs) == 0:
            return self.parent()(0)
        n = Integer(len(self._relprecs) + len(right._relprecs) - 1).exact_log(2) + 1
        precpoly1 = self._getprecpoly(n) * right._getvalpoly(n)
        precpoly2 = self._getvalpoly(n) * right._getprecpoly(n)
        # These two will be the same length
        tn = Integer(1) << n
        preclist = [min(a.valuation(tn), b.valuation(tn)) for (a, b) in zip(precpoly1.list(), precpoly2.list())]
        answer = Polynomial_padic_capped_relative_dense(self.parent(), (zzpoly, self._valbase + right._valbase, preclist, False, None, None), construct = True)
        answer._reduce_poly()
        return answer

    def _lmul_(self, right):
        return self._rmul_(right)

    def _rmul_(self, left):
        """
        Return ``self`` multiplied by a constant.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: a = t^4 + K(13,5)*t^2 + 13
            sage: K(13,7) * a
            (13 + O(13^7))*t^4 + (13^2 + O(13^6))*t^2 + 13^2 + O(13^8)
        """
        return None
        # The code below has never been tested and is somehow subtly broken.

        if self._valaddeds is None:
            self._comp_valaddeds()
        if left != 0:
            val, unit = left.val_unit()
            left_rprec = left.precision_relative()
            relprecs = [min(left_rprec + self._valaddeds[i], self._relprecs[i]) for i in range(len(self._relprecs))]
        elif left._is_exact_zero():
            return Polynomial_padic_capped_relative_dense(self.parent(), [])
        else:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly.parent()(0), self._valbase + left.valuation(), self._valaddeds, False, self._valaddeds, None), construct = True)
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly._rmul_(unit), self._valbase + val, relprecs, False, self._valaddeds, None), construct = True)

    def _neg_(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: K = Qp(13,2)
            sage: R.<t> = K[]
            sage: a = t^4 + 13*t^2 + 4
            sage: -a
            (12 + 12*13 + O(13^2))*t^4 + (12*13 + 12*13^2 + O(13^3))*t^2 + 9 + 12*13 + O(13^2)
        """
        return Polynomial_padic_capped_relative_dense(self.parent(), (-self._poly, self._valbase, self._relprecs, False, self._valaddeds, None), construct = True)

    def lshift_coeffs(self, shift, no_list = False):
        """
        Return a new polynomials whose coefficients are multiplied by p^shift.

        EXAMPLES::

            sage: K = Qp(13, 4)
            sage: R.<t> = K[]
            sage: a = t + 52
            sage: a.lshift_coeffs(3)
            (13^3 + O(13^7))*t + 4*13^4 + O(13^8)
        """
        if shift < 0:
            return self.rshift_coeffs(-shift, no_list)
        if no_list or self._list is None:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase + shift, self._relprecs, False, self._valaddeds, None), construct = True)
        else:
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase + shift, self._relprecs, False, self._valaddeds, [c.__lshift__(shift) for c in self._list]), construct = True)

    def rshift_coeffs(self, shift, no_list=False):
        """
        Return a new polynomial whose coefficients are p-adically
        shifted to the right by ``shift``.

        .. NOTE::

            Type ``Qp(5)(0).__rshift__?`` for more information.

        EXAMPLES::

            sage: K = Zp(13, 4)
            sage: R.<t> = K[]
            sage: a = t^2 + K(13,3)*t + 169; a
            (1 + O(13^4))*t^2 + (13 + O(13^3))*t + 13^2 + O(13^6)
            sage: b = a.rshift_coeffs(1); b
            O(13^3)*t^2 + (1 + O(13^2))*t + 13 + O(13^5)
            sage: b.list()
            [13 + O(13^5), 1 + O(13^2), O(13^3)]
            sage: b = a.rshift_coeffs(2); b
            O(13^2)*t^2 + O(13)*t + 1 + O(13^4)
            sage: b.list()
            [1 + O(13^4), O(13), O(13^2)]
        """
        if shift < 0:
            return self.lshift_coeffs(-shift, no_list) # We can't just absorb this into the next if statement because we allow rshift to preserve _normalized
        if self.base_ring().is_field() or shift <= self._valbase:
            if no_list or self._list is None:
                return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase - shift, self._relprecs, self._normalized, self._valaddeds, None), construct = True)
            else:
                return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly, self._valbase - shift, self._relprecs, self._normalized, self._valaddeds, [c.__rshift__(shift) for c in self._list]), construct = True)
        else:
            shift = shift - self._valbase
            fdiv = self.base_ring().prime_pow(shift)
            return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly // fdiv, 0, [0 if a <= shift else a - shift for a in self._relprecs], False, None, None), construct = True)

    #def __floordiv__(self, right):
    #    if is_Polynomial(right) and right.is_constant() and right[0] in self.base_ring():
    #        d = self.base_ring()(right[0])
    #    elif (right in self.base_ring()):
    #        d = self.base_ring()(right)
    #    else:
    #        raise NotImplementedError
    #    return self._rmul_(self.base_ring()(~d.unit_part())).rshift_coeffs(d.valuation())

    def _unsafe_mutate(self, n, value):
        """
        It's a really bad idea to use this function for p-adic
        polynomials.  There are speed issues, and it may not be
        bug-free currently.
        """
        n = int(n)
        value = self.base_ring()(value)
        if self.is_gen():
            raise ValueError("cannot modify generator")
        if n < 0:
            raise IndexError("n must be >= 0")
        if self._valbase is infinity:
            if value._is_exact_zero():
                return
            self._valbase = value.valuation()
            if value != 0:
                self._poly._unsafe_mutate(self, n, value.unit_part().lift())
                self._relprecs = [infinity] * n + [value.precision_relative()]
            else:
                self._relprecs = [infinity] * n + [0]
            self._valaddeds = [infinity] * n + [0]
            zero = self.base_ring()(0)
            self._list = [zero] * n + [value]
            self._normalized = True
        elif value.valuation() >= self._valbase:
            # _valbase and _normalized stay the same
            if value != 0:
                self._poly._unsafe_mutate(self, n, (value.__rshift__(self._valbase)).lift())
            else:
                self._poly._unsafe_mutate(self, n, 0)
            if n < len(self._relprecs):
                self._relprecs[n] = value.precision_absolute() - self._valbase
                if self._valaddeds is not None:
                    self._valaddeds[n] = value.valuation() - self._valbase
                if self._list is not None:
                    self._list[n] = value
            else:
                self._relprecs.extend([infinity] * (n - len(self._relprecs)) + [value.precision_absolute() - self._valbase])
                if self._valaddeds is not None:
                    self._valaddeds.extend([infinity] * (n - len(self._relprecs)) + [value.valuation() - self._valbase])
                if self._list is not None:
                    zero = self.base_ring()(0)
                    self._list.extend([zero] * (n - len(self._relprecs)) + [value])
        else:
            basediff = self._valbase - value.valuation()
            self._valbase = value.valuation()
            if self._valaddeds is not None:
                self._valaddeds = [c + basediff for c in self._valaddeds]
            self._poly = self._poly * self.base_ring().prime_pow(basediff)
            if value != 0:
                self._poly._unsafe_mutate(self, n, value.unit_part().lift())
            else:
                self._poly._unsafe_mutate(self, n, 0)
            if n < len(self._relprecs):
                self._relprecs[n] = value.precision_relative()
            else:
                self._relprecs.extend([infinity] * (n - len(self._relprecs)) + [value.precision_relative()])
            self._normalized = False
            if self._list is not None:
                if n < len(self._list):
                    self._list[n] = value
                else:
                    zero = self._base_ring()(0)
                    self._list.extend([zero] * (n - len(self._list)) + [value])

    def __pari__(self, variable=None):
        """
        Return ``self`` as a Pari object.
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.list()).Polrev(variable)

    def __copy__(self):
        """
        Return a copy of ``self``.
        """
        return Polynomial_padic_capped_relative_dense(self.parent(), (copy.copy(self._poly), self._valbase, copy.copy(self._relprecs), self._normalized, copy.copy(self._valaddeds), copy.copy(self._list)), construct = True)

    def degree(self, secure=False):
        """
        Return the degree of ``self``.

        INPUT:

        - secure  -- a boolean (default: ``False``)

        If ``secure`` is ``True`` and the degree of this polynomial
        is not determined (because the leading coefficient is
        indistinguishable from 0), an error is raised.

        If ``secure`` is ``False``, the returned value is the largest
        $n$ so that the coefficient of $x^n$ does not compare equal
        to $0$.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.degree()
            1
            sage: (f-T).degree()
            0
            sage: (f-T).degree(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: the leading coefficient is indistinguishable from 0

            sage: x = O(3^5)
            sage: li = [3^i * x for i in range(0,5)]; li
            [O(3^5), O(3^6), O(3^7), O(3^8), O(3^9)]
            sage: f = R(li); f
            O(3^9)*T^4 + O(3^8)*T^3 + O(3^7)*T^2 + O(3^6)*T + O(3^5)
            sage: f.degree()
            -1
            sage: f.degree(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: the leading coefficient is indistinguishable from 0
        """
        self._normalize()
        deg = Integer(self._poly.degree())
        if secure and deg < self.prec_degree():
            raise PrecisionError("the leading coefficient is "
                                 "indistinguishable from 0")
        return deg

    def prec_degree(self):
        """
        Return the largest $n$ so that precision information is
        stored about the coefficient of $x^n$.

        Always greater than or equal to degree.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.prec_degree()
            1
        """
        return len(self._relprecs) - 1

    def precision_absolute(self, n = None):
        """
        Return absolute precision information about ``self``.

        INPUT:

        ``self`` -- a p-adic polynomial

        n -- ``None`` or an integer (default ``None``).

        OUTPUT:

        If n == None, returns a list of absolute precisions of
        coefficients.  Otherwise, returns the absolute precision of
        the coefficient of x^n.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.precision_absolute()
            [10, 10]
        """
        if n is None:
            return [c + self._valbase for c in self._relprecs]
        return self._relprecs[n] + self._valbase

    def precision_relative(self, n = None):
        """
        Return relative precision information about ``self``.

        INPUT:

        ``self`` -- a p-adic polynomial

        n -- ``None`` or an integer (default ``None``).

        OUTPUT:

        If n == None, returns a list of relative precisions of
        coefficients.  Otherwise, returns the relative precision of
        the coefficient of x^n.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.precision_relative()
            [10, 10]
        """
        if n is None:
            self._normalize()
            return copy.copy(self._relprecs)
        n = int(n)
        if n < 0 or n >= len(self._relprecs) or self._relprecs[n] is infinity:
            return Integer(0)
        if self._valaddeds is None:
            return self._relprecs[n] - self._poly[n].valuation(self.base_ring().prime())
        else:
            return self._relprecs[n] - self._valaddeds[n]

    def valuation_of_coefficient(self, n=None):
        """
        Return valuation information about ``self``'s coefficients.

        INPUT:

        ``self`` -- a p-adic polynomial

        n -- ``None`` or an integer (default ``None``).

        OUTPUT:

        If n == None, returns a list of valuations of coefficients.  Otherwise,
        returns the valuation of the coefficient of x^n.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.valuation_of_coefficient(1)
            0
        """
        if self._valaddeds is None:
            self._comp_valaddeds()
        if n is None:
            self._normalize()
            return [ c + self._valbase for c in self._valaddeds ]
        n = int(n)
        if n < 0 or n >= len(self._relprecs):
            return infinity
        return self._valbase + self._valaddeds[n]

    def valuation(self, val_of_var=None):
        """
        Return the valuation of ``self``.

        INPUT:

        ``self`` -- a p-adic polynomial

        val_of_var -- ``None`` or a rational (default ``None``).

        OUTPUT:

        If val_of_var == None, returns the largest power of the
        variable dividing self.  Otherwise, returns the valuation of
        ``self`` where the variable is assigned valuation val_of_var

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: R.<T> = K[]
            sage: f = T + 2; f
            (1 + O(3^10))*T + 2 + O(3^10)
            sage: f.valuation()
            0
        """
        if val_of_var is None:
            return self._poly.valuation()
        if self._valaddeds is None:
            self._comp_valaddeds()
        return self._valbase + min([self._valaddeds[i] + val_of_var * i for i in range(len(self._valaddeds))])

    def reverse(self, degree=None):
        """
        Return the reverse of the input polynomial, thought as a polynomial of
        degree ``degree``.

        If `f` is a degree-`d` polynomial, its reverse is `x^d f(1/x)`.

        INPUT:

        - ``degree`` (``None`` or an integer) - if specified, truncate or zero
          pad the list of coefficients to this degree before reversing it.

        EXAMPLES::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: f = t^3 + 4*t; f
            (1 + O(13^7))*t^3 + (4 + O(13^7))*t
            sage: f.reverse()
            0*t^3 + (4 + O(13^7))*t^2 + 1 + O(13^7)
            sage: f.reverse(3)
            0*t^3 + (4 + O(13^7))*t^2 + 1 + O(13^7)
            sage: f.reverse(2)
            0*t^2 + (4 + O(13^7))*t
            sage: f.reverse(4)
            0*t^4 + (4 + O(13^7))*t^3 + (1 + O(13^7))*t
            sage: f.reverse(6)
            0*t^6 + (4 + O(13^7))*t^5 + (1 + O(13^7))*t^3

        TESTS:

        Check that this implementation is compatible with the generic one::

            sage: all(f.reverse(d) == Polynomial.reverse(f, d)
            ....:     for d in [None,0,1,2,3,4,5])
            True
        """
        n = self._poly.degree() if degree is None else degree
        zzlist = self._poly.list()[:(n + 1)] + [0] * (n - self._poly.degree())
        zzlist.reverse()
        relprec = self._relprecs[:(n + 1)] + [infinity] * (n - self.prec_degree())
        relprec.reverse()
        if self._valaddeds is None:
            valadded = None
        else:
            valadded = self._valaddeds[:(n + 1)] + [infinity] * (n - self.prec_degree())
            valadded.reverse()
        if self._list is None:
            L = None
        else:
            L = self._list[:(n + 1)] + [self.base_ring()(0)] * (n - self.prec_degree())
            L.reverse()
        return Polynomial_padic_capped_relative_dense(self.parent(), (self._poly.parent()(zzlist), self._valbase, relprec, self._normalized, valadded, L), construct = True)

    def rescale(self, a):
        r"""
        Return f(a*X)

        .. TODO::

            Need to write this function for integer polynomials before
            this works.

        EXAMPLES::

            sage: K = Zp(13, 5)
            sage: R.<t> = K[]
            sage: f = t^3 + K(13, 3) * t
            sage: f.rescale(2)  # not implemented
        """
        negval = False
        try:
            a = self.base_ring()(a)
        except ValueError as msg:
            if msg == "element has negative valuation.":
                negval = True
            else:
                raise ValueError(msg)
        if negval:
            return self.parent().base_extend(self.base_ring().fraction_field())(self).rescale(a)
        if self.base_ring().is_field() and a.valuation() < 0:
            D = self.prec_degree()
            return a**D * self.reverse(D).rescale(~a).reverse(D)
        aval = a.valuation()
        arprec = a.precision_relative()
        if self._valaddeds is None:
            self._comp_valaddeds()
        valadded = [self._valaddeds[i] + aval * i for i in range(len(self._valaddeds))]
        relprec = [infinity if (self._relprecs[i] is infinity) else (min(self._relprecs[i] - self._valaddeds[i], arprec) + aval * i + self._valaddeds[i]) for i in range(len(self._relprecs))]
        relprec[0] = self._relprecs[0]
        if a == 0:
            zzpoly = self._poly.parent()(0)
        else:
            zzpoly = self._poly.rescale(Integer(a))
        return Polynomial_padic_capped_relative_dense(self.parent(), (zzpoly, self._valbase, relprec, False, valadded, None), construct = True)

    def quo_rem(self, right, secure=False):
        """
        Return the quotient and remainder in division of ``self`` by ``right``.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: R.<T> = K[]
            sage: f = T + 2
            sage: g = T**4 + 3*T+22
            sage: g.quo_rem(f)
            ((1 + O(3^10))*T^3 + (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10))*T^2 + (1 + 3 + O(3^10))*T + 1 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10),
             2 + 3 + 3^3 + O(3^10))

        TESTS:

        Verify that :trac:`15188` has been resolved::

            sage: R.<x> = Qp(3)[]
            sage: x.quo_rem(x)
            (1 + O(3^20), 0)

        """
        return self._quo_rem_list(right, secure=secure)

    def _quo_rem_naive(self, right):
        """
        An implementation of quo_rem that doesn't have good run-time
        or precision characteristics.

        A better one is :meth:`_quo_rem_list`.
        """
        K = self.base_ring().fraction_field()
        f = self.base_extend(K)
        g = right.base_extend(K)
        if g == 0:
            raise ZeroDivisionError("cannot divide by a polynomial "
                                    "indistinguishable from 0")
        x = f.parent().gen()
        quo = f.parent()(0)
        while f.degree() >= g.degree():
            a = f.leading_coefficient() / g.leading_coefficient()
            quo = quo + a * (x ** (f.degree() - g.degree()))
            f = f - a * (x ** (f.degree() - g.degree())) * g
        return (quo, f)

    def _quo_rem_list(self, right, secure):
        """
        An implementation of quo_rem using lists of coefficients.

        Faster than :meth:`_quo_rem_naive`.

        AUTHOR:

        - Xavier Caruso (2013-03)
        """
        if right.is_zero():
            raise ZeroDivisionError("cannot divide by a polynomial "
                                    "indistinguishable from 0")
        a = self.list()
        da = len(a) - 1
        b = right.list()
        db = right.degree(secure=secure)
        inv = ~b[db]
        q = [ ]
        for i in range(da, db - 1, -1):
            c = inv * a[i]
            q.append(c)
            for j in range(db):
                a[j + i - db] -= c * b[j]
        q.reverse()
        K = self.base_ring().fraction_field()
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        parent = PolynomialRing(K, name=self.parent().variable_name())
        return parent(q), parent(a[:db])

    #def gcd(self, right):
    #    raise NotImplementedError

    #def lcm(self, right):
    #    raise NotImplementedError

    #def discriminant(self):
    #    raise NotImplementedError

    def disc(self):
        return self.discriminant()

    #def resultant(self):
    #    raise NotImplementedError

    def newton_polygon(self):
        r"""
        Return the Newton polygon of this polynomial.

        .. NOTE::

            If some coefficients have not enough precision an error is raised.

        OUTPUT:

        - a Newton polygon

        EXAMPLES::

            sage: K = Qp(2, prec=5)
            sage: P.<x> = K[]
            sage: f = x^4 + 2^3*x^3 + 2^13*x^2 + 2^21*x + 2^37
            sage: f.newton_polygon()
            Finite Newton polygon with 4 vertices: (0, 37), (1, 21), (3, 3), (4, 0)

            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10
            sage: f.newton_polygon()
            Finite Newton polygon with 4 vertices: (0, 1), (1, 0), (4, 0), (10, 2)

        Here is an example where the computation fails because precision is
        not sufficient::

            sage: g = f + K(0,0)*t^4; g
            (5^2 + O(5^22))*t^10 + O(5^0)*t^4 + (3 + O(5^20))*t + 5 + O(5^21)
            sage: g.newton_polygon()
            Traceback (most recent call last):
            ...
            PrecisionError: The coefficient of t^4 has not enough precision

        TESTS::

            sage: (5*f).newton_polygon()
            Finite Newton polygon with 4 vertices: (0, 2), (1, 1), (4, 1), (10, 3)

        AUTHOR:

        - Xavier Caruso (2013-03-20)
        """
        if self._valaddeds is None:
            self._comp_valaddeds()
        from sage.geometry.newton_polygon import NewtonPolygon
        valbase = self._valbase
        polygon = NewtonPolygon([(x, val + valbase)
                                 for x, val in enumerate(self._valaddeds)])
        polygon_prec = NewtonPolygon([(x, val + valbase)
                                      for x, val in enumerate(self._relprecs)])
        vertices = polygon.vertices(copy=False)
        vertices_prec = polygon_prec.vertices(copy=False)

        # The two following tests should always fail (i.e. the corresponding errors
        # should never be raised). However, it's probably safer to keep them.
        if vertices[0][0] > vertices_prec[0][0]:
            raise PrecisionError("The constant coefficient has not enough precision")
        if vertices[-1][0] < vertices_prec[-1][0]:
            raise PrecisionError("The leading coefficient has not enough precision")

        for (x, y) in vertices:
            if polygon_prec(x) <= y:
                raise PrecisionError("The coefficient of %s^%s has not enough precision" % (self.parent().variable_name(), x))
        return polygon

    def is_eisenstein(self, secure=False):
        """
        Return ``True`` if this polynomial is an Eisenstein polynomial.

        EXAMPLES::

            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 5*t + t^4
            sage: f.is_eisenstein()
            True

        TESTS::

            sage: f = R([K(5,1),0,0,1]); f
            (1 + O(5^20))*t^3 + O(5)
            sage: f.is_eisenstein()
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision on the constant coefficient

            sage: g = R([5,K(0,0),0,1]); g
            (1 + O(5^20))*t^3 + O(5^0)*t + 5 + O(5^21)
            sage: g.is_eisenstein()
            True
            sage: g.is_eisenstein(secure=True)
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision on the coefficient of t

        AUTHOR:

        - Xavier Caruso (2013-03)
        """
        deg = self.degree()
        if secure and self.prec_degree() > deg:
            raise PrecisionError("The degree of the polynomial is not determined")
        if self._valaddeds is None:
            self._comp_valaddeds()
        compval = 1 - self._valbase
        valaddeds = self._valaddeds
        relprecs = self._relprecs
        if relprecs[0] <= compval:   # not enough precision
            if valaddeds[0] < relprecs[0]:
                return False
            raise PrecisionError("Not enough precision on the constant coefficient")
        else:
            if valaddeds[0] != compval:
                return False
        for i in range(1, deg):
            if relprecs[i] < compval:   # not enough precision
                if valaddeds[i] < relprecs[i]:
                    return False
                if secure:
                    if i == 1:
                        raise PrecisionError("Not enough precision on the coefficient of %s" % self.variable_name())
                    else:
                        raise PrecisionError("Not enough precision on the coefficient of %s^%s" % (self.variable_name(), i))
            else:
                if valaddeds[i] < compval:
                    return False
        if valaddeds[deg] != -self._valbase:
            return False
        return True

    def newton_slopes(self, repetition=True):
        """
        Return a list of the Newton slopes of this polynomial.

        These are the valuations of the roots of this polynomial.

        If ``repetition`` is ``True``, each slope is repeated a number of
        times equal to its multiplicity. Otherwise it appears only one time.

        INPUT:

        - ``repetition`` -- boolean (default ``True``)

        OUTPUT:

        - a list of rationals

        EXAMPLES::

            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10
            sage: f.newton_polygon()
            Finite Newton polygon with 4 vertices: (0, 1), (1, 0), (4, 0),
            (10, 2)
            sage: f.newton_slopes()
            [1, 0, 0, 0, -1/3, -1/3, -1/3, -1/3, -1/3, -1/3]

            sage: f.newton_slopes(repetition=False)
            [1, 0, -1/3]

        AUTHOR:

        - Xavier Caruso (2013-03-20)
        """
        polygon = self.newton_polygon()
        return [-s for s in polygon.slopes(repetition=repetition)]

    def factor_mod(self):
        r"""
        Return the factorization of ``self`` modulo `p`.
        """
        self._normalize()
        if self._valbase < 0:
            raise ValueError("Polynomial does not have integral coefficients")
        elif self._valbase > 0:
            raise ValueError("Factorization of the zero polynomial not defined")
        elif min(self._relprecs) <= 0:
            raise PrecisionError("Polynomial is not known to high enough precision")
        return self._poly.factor_mod(self.base_ring().prime())


def _extend_by_infinity(L, n):
    return L + [infinity] * (n - len(L))


def make_padic_poly(parent, x, version):
    if version == 0:
        return parent(x, construct = True)
    else:
        raise ValueError("unknown pickling version")
