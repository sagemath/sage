# sage_setup: distribution = sagemath-flint
# distutils: libraries = gmp NTL_LIBRARIES
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"Univariate rational functions over prime fields"

from cysignals.signals cimport sig_on, sig_off

from sage.rings.finite_rings.stdint cimport INTEGER_MOD_INT32_LIMIT

from sage.libs.gmp.mpz cimport *
from sage.libs.flint.nmod_poly cimport *
from sage.libs.flint.ulong_extras cimport n_jacobi
from sage.structure.element cimport Element, FieldElement
from sage.rings.integer_ring import ZZ
from sage.rings.fraction_field import FractionField_1poly_field
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_zmod_flint cimport Polynomial_zmod_flint, get_cparent

from sage.structure.richcmp cimport rich_to_bool
from sage.rings.finite_rings.integer_mod cimport mod_inverse_int


class FpT(FractionField_1poly_field):
    r"""
    This class represents the fraction field `\GF{p}(T)` for `2 < p < \sqrt{2^31-1}`.

    EXAMPLES::

        sage: R.<T> = GF(71)[]
        sage: K = FractionField(R); K
        Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 71
        sage: 1-1/T
        (T + 70)/T
        sage: parent(1-1/T) is K
        True
    """
    INTEGER_LIMIT = INTEGER_MOD_INT32_LIMIT

    def __init__(self, R, names=None):  # we include names so that one can use the syntax K.<t> = FpT(GF(5)['t']).  It's actually ignored
        """
        INPUT:

        - ``R`` -- a dense polynomial ring over a finite field of prime order
          `p` with `2 < p < 2^{16}`

        EXAMPLES::

            sage: R.<x> = GF(31)[]
            sage: K = R.fraction_field(); K
            Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 31

        TESTS::

            sage: from sage.rings.fraction_field_FpT import FpT
            sage: FpT(PolynomialRing(GF(37), ['x'], sparse=True))
            Traceback (most recent call last):
            ...
            TypeError: unsupported polynomial ring
        """
        cdef long p = R.base_ring().characteristic()
        assert 2 < p < FpT.INTEGER_LIMIT
        if not issubclass(R.element_class, Polynomial_zmod_flint):
            raise TypeError("unsupported polynomial ring")
        self.p = p
        self.poly_ring = R
        FractionField_1poly_field.__init__(self, R, element_class=FpTElement)
        self._populate_coercion_lists_(coerce_list=[Polyring_FpT_coerce(self), Fp_FpT_coerce(self), ZZ_FpT_coerce(self)])

    def __iter__(self):
        """
        Return an iterator over this fraction field.

        EXAMPLES::

            sage: R.<t> = GF(3)[]; K = R.fraction_field()
            sage: iter(K)
            <sage.rings.fraction_field_FpT.FpT_iter object at ...>
        """
        return self.iter()

    def iter(self, bound=None, start=None):
        """
        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(5)['t'])
            sage: list(R.iter(2))[350:355]
            [(t^2 + t + 1)/(t + 2),
             (t^2 + t + 2)/(t + 2),
             (t^2 + t + 4)/(t + 2),
             (t^2 + 2*t + 1)/(t + 2),
             (t^2 + 2*t + 2)/(t + 2)]
        """
        return FpT_iter(self, bound, start)


cdef class FpTElement(FieldElement):
    """
    An element of an :class:`FpT` fraction field.

    TESTS::

        sage: R.<t> = GF(5)[]
        sage: K = R.fraction_field()
        sage: A.<x> = K[]
        sage: x.divides(x)  # Testing issue #27064
        True
    """

    def __init__(self, parent, numer, denom=1, coerce=True, reduce=True):
        """
        INPUT:

        - ``parent`` -- the Fraction field containing this element
        - ``numer`` -- something that can be converted into the polynomial
          ring, giving the numerator
        - ``denom`` -- something that can be converted into the polynomial
          ring, giving the numerator (default: 1)

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(5)['t'])
            sage: R(7)
            2
        """
        super().__init__(parent)
        if coerce:
            numer = parent.poly_ring(numer)
            denom = parent.poly_ring(denom)
        self.p = parent.p
        nmod_poly_init(self._numer, self.p)
        nmod_poly_init(self._denom, self.p)
        self.initialized = True
        cdef long n
        for n, a in enumerate(numer):
            nmod_poly_set_coeff_ui(self._numer, n, a)
        for n, a in enumerate(denom):
            nmod_poly_set_coeff_ui(self._denom, n, a)
        if reduce:
            normalize(self._numer, self._denom, self.p)

    def __dealloc__(self):
        """
        Deallocation.

        EXAMPLES::

            sage: K = GF(11)['t'].fraction_field()
            sage: t = K.gen()
            sage: del t # indirect doctest
        """
        if self.initialized:
            nmod_poly_clear(self._numer)
            nmod_poly_clear(self._denom)

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: K = GF(11)['t'].fraction_field()
            sage: loads(dumps(K.gen()))
            t
            sage: loads(dumps(1/K.gen()))
            1/t
        """
        return (unpickle_FpT_element,
                (self._parent, self.numer(), self.denom()))

    cdef FpTElement _new_c(self):
        """
        Create a new FpTElement in the same field, leaving the value to be
        initialized.
        """
        cdef FpTElement x = <FpTElement>FpTElement.__new__(FpTElement)
        x._parent = self._parent
        x.p = self.p
        nmod_poly_init_preinv(x._numer, x.p, self._numer.mod.ninv)
        nmod_poly_init_preinv(x._denom, x.p, self._numer.mod.ninv)
        x.initialized = True
        return x

    cdef FpTElement _copy_c(self):
        """
        Create a new FpTElement in the same field, with the same value as
        ``self``.
        """
        cdef FpTElement x = <FpTElement>FpTElement.__new__(FpTElement)
        x._parent = self._parent
        x.p = self.p
        nmod_poly_init2_preinv(x._numer, x.p, self._numer.mod.ninv, self._numer.length)
        nmod_poly_init2_preinv(x._denom, x.p, self._denom.mod.ninv, self._denom.length)
        nmod_poly_set(x._numer, self._numer)
        nmod_poly_set(x._denom, self._denom)
        x.initialized = True
        return x

    def numer(self):
        """
        Return the numerator of this element, as an element of the polynomial ring.

        EXAMPLES::

            sage: K = GF(11)['t'].fraction_field()
            sage: t = K.gen(0); a = (t + 1/t)^3 - 1
            sage: a.numer()
            t^6 + 3*t^4 + 10*t^3 + 3*t^2 + 1
        """
        return self.numerator()

    cpdef numerator(self):
        """
        Return the numerator of this element, as an element of the polynomial ring.

        EXAMPLES::

            sage: K = GF(11)['t'].fraction_field()
            sage: t = K.gen(0); a = (t + 1/t)^3 - 1
            sage: a.numerator()
            t^6 + 3*t^4 + 10*t^3 + 3*t^2 + 1
        """
        cdef Polynomial_zmod_flint res = <Polynomial_zmod_flint>Polynomial_zmod_flint.__new__(Polynomial_zmod_flint)
        nmod_poly_init2_preinv(&res.x, self.p, self._numer.mod.ninv, self._numer.length)
        nmod_poly_set(&res.x, self._numer)
        res._parent = self._parent.poly_ring
        res._cparent = get_cparent(self._parent.poly_ring)
        return res

    def denom(self):
        """
        Return the denominator of this element, as an element of the polynomial ring.

        EXAMPLES::

            sage: K = GF(11)['t'].fraction_field()
            sage: t = K.gen(0); a = (t + 1/t)^3 - 1
            sage: a.denom()
            t^3
        """
        return self.denominator()

    cpdef denominator(self):
        """
        Return the denominator of this element, as an element of the polynomial ring.

        EXAMPLES::

            sage: K = GF(11)['t'].fraction_field()
            sage: t = K.gen(0); a = (t + 1/t)^3 - 1
            sage: a.denominator()
            t^3
        """
        cdef Polynomial_zmod_flint res = <Polynomial_zmod_flint>Polynomial_zmod_flint.__new__(Polynomial_zmod_flint)
        nmod_poly_init2_preinv(&res.x, self.p, self._denom.mod.ninv, self._denom.length)
        nmod_poly_set(&res.x, self._denom)
        res._parent = self._parent.poly_ring
        res._cparent = get_cparent(self._parent.poly_ring)
        return res

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: K = Frac(GF(5)['t'])
            sage: t = K.gen()
            sage: t(3)
            3
            sage: f = t^2/(1-t)
            sage: f(2)
            1
            sage: f(t)
            4*t^2/(t + 4)
            sage: f(t^3)
            4*t^6/(t^3 + 4)
            sage: f((t+1)/t^3)
            (t^2 + 2*t + 1)/(t^6 + 4*t^4 + 4*t^3)
        """
        return self.numer()(*args, **kwds) / self.denom()(*args, **kwds)

    def subs(self, in_dict=None, *args, **kwds):
        """
        EXAMPLES::

            sage: K = Frac(GF(11)['t'])
            sage: t = K.gen()
            sage: f = (t+1)/(t-1)
            sage: f.subs(t=2)
            3
            sage: f.subs(X=2)
            (t + 1)/(t + 10)
        """
        return self.numer().subs(in_dict, *args, **kwds) / self.denom().subs(in_dict, *args, **kwds)

    def valuation(self, v):
        """
        Return the valuation of ``self`` at `v`.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: f = (t+1)^2 * (t^2+t+1) / (t-1)^3
            sage: f.valuation(t+1)
            2
            sage: f.valuation(t-1)
            -3
            sage: f.valuation(t)
            0
        """
        return self.numer().valuation(v) - self.denom().valuation(v)

    def factor(self):
        """
        EXAMPLES::

            sage: K = Frac(GF(5)['t'])
            sage: t = K.gen()
            sage: f = 2 * (t+1) * (t^2+t+1)^2 / (t-1)
            sage: factor(f)
            (2) * (t + 4)^-1 * (t + 1) * (t^2 + t + 1)^2
        """
        return self.numer().factor() / self.denom().factor()

    def _repr_(self):
        """
        Return a string representation of this element.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(17)['t'])
            sage: -t # indirect doctest
            16*t
            sage: 1/t
            1/t
            sage: 1/(t+1)
            1/(t + 1)
            sage: 1-t/t
            0
            sage: (1-t)/t
            (16*t + 1)/t
        """
        if nmod_poly_degree(self._denom) == 0 and nmod_poly_get_coeff_ui(self._denom, 0) == 1:
            return repr(self.numer())
        else:
            numer_s = repr(self.numer())
            denom_s = repr(self.denom())
            if '-' in numer_s or '+' in numer_s:
                numer_s = "(%s)" % numer_s
            if '-' in denom_s or '+' in denom_s:
                denom_s = "(%s)" % denom_s
            return "%s/%s" % (numer_s, denom_s)

    def _latex_(self):
        r"""
        Return a latex representation of this element.

        EXAMPLES::

            sage: K = GF(7)['t'].fraction_field(); t = K.gen(0)
            sage: latex(t^2 + 1) # indirect doctest
            t^{2} + 1
            sage: latex((t + 1)/(t-1))
            \frac{t + 1}{t + 6}
        """
        if nmod_poly_degree(self._denom) == 0 and nmod_poly_get_coeff_ui(self._denom, 0) == 1:
            return self.numer()._latex_()
        else:
            return "\\frac{%s}{%s}" % (self.numer()._latex_(), self.denom()._latex_())

    cpdef _richcmp_(self, other, int op):
        """
        Compare this with another element.

        The ordering is arbitrary, but it is an ordering, and it is
        consistent between runs.  It has nothing to do with the
        algebra structure.

        TESTS::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(7)['t'])
            sage: t == t
            True
            sage: t == -t
            False
            sage: -t == 6*t
            True
            sage: 1/t == 1/t
            True
            sage: 1/t == 1/(t+1)
            False
            sage: 2*t/t == 2
            True
            sage: 2*t/2 == t
            True

            sage: a = (t^3 + 3*t)/(5*t-2); b = (t^2-2)/(t-1)
            sage: b < a
            True
            sage: a < b
            False
            sage: 1/a < b
            True
            sage: b < 1/a
            False

        ::

            sage: K = Frac(GF(5)['t']); t = K.gen()
            sage: t == 1
            False
            sage: t + 1 < t^2
            True
        """
        # They are normalized.
        cdef int j = sage_cmp_nmod_poly_t(self._numer, (<FpTElement>other)._numer)
        if j:
            return rich_to_bool(op, j)
        j = sage_cmp_nmod_poly_t(self._denom, (<FpTElement>other)._denom)
        return rich_to_bool(op, j)

    def __hash__(self):
        """
        Return a hash value for this element.

        TESTS::

            sage: from sage.rings.fraction_field_FpT import *
            sage: K.<t> = FpT(GF(7)['t'])
            sage: hash(K(0))
            0
            sage: hash(K(5))
            5
            sage: set([1, t, 1/t, t, t, 1/t, 1+1/t, t/t])
            {1, 1/t, t, (t + 1)/t}
            sage: a = (t+1)/(t^2-1); hash(a) == hash((a.numer(),a.denom()))
            True
        """
        if self.denom() == 1:
            return hash(self.numer())
        return hash((self.numer(), self.denom()))

    def __neg__(self):
        """
        Negate this element.

        EXAMPLES::

            sage: K = GF(5)['t'].fraction_field(); t = K.gen(0)
            sage: a = (t^2 + 2)/(t-1)
            sage: -a # indirect doctest
            (4*t^2 + 3)/(t + 4)
        """
        cdef FpTElement x = self._copy_c()
        nmod_poly_neg(x._numer, x._numer)
        return x

    def __invert__(self):
        """
        Return the multiplicative inverse of this element.

        EXAMPLES::

            sage: K = GF(5)['t'].fraction_field(); t = K.gen(0)
            sage: a = (t^2 + 2)/(t-1)
            sage: ~a # indirect doctest
            (t + 4)/(t^2 + 2)
        """
        if nmod_poly_degree(self._numer) == -1:
            raise ZeroDivisionError
        cdef FpTElement x = self._copy_c()
        nmod_poly_swap(x._numer, x._denom)
        return x

    cpdef _add_(self, _other):
        """
        Return the sum of this fraction field element and another.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(7)['t'])
            sage: t + t # indirect doctest
            2*t
            sage: (t + 3) + (t + 10)
            2*t + 6
            sage: sum([t] * 7)
            0
            sage: 1/t + t
            (t^2 + 1)/t
            sage: 1/t + 1/t^2
            (t + 1)/t^2
        """
        cdef FpTElement other = <FpTElement>_other
        cdef FpTElement x = self._new_c()
        nmod_poly_mul(x._numer, self._numer, other._denom)
        nmod_poly_mul(x._denom, self._denom, other._numer)  # use x._denom as a temp
        nmod_poly_add(x._numer, x._numer, x._denom)
        nmod_poly_mul(x._denom, self._denom, other._denom)
        normalize(x._numer, x._denom, self.p)
        return x

    cpdef _sub_(self, _other):
        """
        Return the difference of this fraction field element and another.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(7)['t'])
            sage: t - t # indirect doctest
            0
            sage: (t + 3) - (t + 11)
            6
        """
        cdef FpTElement other = <FpTElement>_other
        cdef FpTElement x = self._new_c()
        nmod_poly_mul(x._numer, self._numer, other._denom)
        nmod_poly_mul(x._denom, self._denom, other._numer)  # use x._denom as a temp
        nmod_poly_sub(x._numer, x._numer, x._denom)
        nmod_poly_mul(x._denom, self._denom, other._denom)
        normalize(x._numer, x._denom, self.p)
        return x

    cpdef _mul_(self, _other):
        """
        Return the product of this fraction field element and another.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(7)['t'])
            sage: t * t # indirect doctest
            t^2
            sage: (t + 3) * (t + 10)
            t^2 + 6*t + 2
        """
        cdef FpTElement other = <FpTElement>_other
        cdef FpTElement x = self._new_c()
        nmod_poly_mul(x._numer, self._numer, other._numer)
        nmod_poly_mul(x._denom, self._denom, other._denom)
        normalize(x._numer, x._denom, self.p)
        return x

    cpdef _div_(self, _other):
        """
        Return the quotient of this fraction field element and another.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(5)['t'])
            sage: t / t # indirect doctest
            1
            sage: (t + 3) / (t + 11)
            (t + 3)/(t + 1)
            sage: (t^2 + 2*t + 1) / (t + 1)
            t + 1
        """
        cdef FpTElement other = <FpTElement>_other
        if nmod_poly_degree(other._numer) == -1:
            raise ZeroDivisionError
        cdef FpTElement x = self._new_c()
        nmod_poly_mul(x._numer, self._numer, other._denom)
        nmod_poly_mul(x._denom, self._denom, other._numer)
        normalize(x._numer, x._denom, self.p)
        return x

    cpdef FpTElement next(self):
        """
        Iterate through all polynomials, returning the "next" polynomial after this one.

        The strategy is as follows:

        - We always leave the denominator monic.

        - We progress through the elements with both numerator and denominator monic, and with the denominator less than the numerator.
          For each such, we output all the scalar multiples of it, then all of the scalar multiples of its inverse.

        - So if the leading coefficient of the numerator is less than `p-1`, we scale the numerator to increase it by 1.

        - Otherwise, we consider the multiple with numerator and denominator monic.

          - If the numerator is less than the denominator (lexicographically), we return the inverse of that element.

          - If the numerator is greater than the denominator, we invert, and then increase the numerator (remaining monic) until we either get something relatively prime to the new denominator, or we reach the new denominator.  In this case, we increase the denominator and set the numerator to 1.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(3)['t'])
            sage: a = R(0)
            sage: for _ in range(30):
            ....:     a = a.next()
            ....:     print(a)
            1
            2
            1/t
            2/t
            t
            2*t
            1/(t + 1)
            2/(t + 1)
            t + 1
            2*t + 2
            t/(t + 1)
            2*t/(t + 1)
            (t + 1)/t
            (2*t + 2)/t
            1/(t + 2)
            2/(t + 2)
            t + 2
            2*t + 1
            t/(t + 2)
            2*t/(t + 2)
            (t + 2)/t
            (2*t + 1)/t
            (t + 1)/(t + 2)
            (2*t + 2)/(t + 2)
            (t + 2)/(t + 1)
            (2*t + 1)/(t + 1)
            1/t^2
            2/t^2
            t^2
            2*t^2
        """
        cdef FpTElement next = self._copy_c()
        cdef long a, lead
        cdef nmod_poly_t g
        if nmod_poly_degree(self._numer) == -1:
            # self should be normalized, so denom == 1
            nmod_poly_set_coeff_ui(next._numer, 0, 1)
            return next
        lead = nmod_poly_leading(next._numer)
        if lead < self.p - 1:
            a = mod_inverse_int(lead, self.p)
            # no overflow since self.p < 2^16
            a = a * (lead + 1) % self.p
            nmod_poly_scalar_mul_nmod(next._numer, next._numer, a)
        else:
            a = mod_inverse_int(lead, self.p)
            nmod_poly_scalar_mul_nmod(next._numer, next._numer, a)
            # now both next._numer and next._denom are monic.  We figure out which is lexicographically bigger:
            a = nmod_poly_cmp(next._numer, next._denom)
            if a == 0:
                # next._numer and next._denom are relatively prime, so they're both 1.
                nmod_poly_inc(next._denom, True)
                return next
            nmod_poly_set(next._denom, next._numer)
            nmod_poly_set(next._numer, self._denom)
            if a < 0:
                # since next._numer is smaller, we flip and return the inverse.
                return next
            elif a > 0:
                # since next._numer is bigger, we're in the flipped phase.  We flip back, and increment the numerator (until we reach the denominator).
                nmod_poly_init(g, self.p)
                try:
                    while True:
                        nmod_poly_inc(next._numer, True)
                        if nmod_poly_equal(next._numer, next._denom):
                            # Since we've reached the denominator, we reset the numerator to 1 and increment the denominator.
                            nmod_poly_inc(next._denom, True)
                            nmod_poly_zero(next._numer)
                            nmod_poly_set_coeff_ui(next._numer, 0, 1)
                            break
                        else:
                            # otherwise, we keep incrementing until we have a relatively prime numerator.
                            nmod_poly_gcd(g, next._numer, next._denom)
                            if nmod_poly_is_one(g):
                                break
                finally:
                    nmod_poly_clear(g)
        return next

    cpdef _sqrt_or_None(self):
        """
        Return the square root of ``self``, or ``None``.

        Differs from sqrt() by not raising an exception.

        TESTS::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(17)['t'])
            sage: sqrt(t^2) # indirect doctest
            t
            sage: sqrt(1/t^2)
            1/t
            sage: sqrt((1+t)^2)
            t + 1
            sage: sqrt((1+t)^2 / t^2)
            (t + 1)/t

            sage: sqrt(4 * (1+t)^2 / t^2)
            (2*t + 2)/t

            sage: sqrt(R(0))
            0
            sage: sqrt(R(-1))
            4

            sage: sqrt(t^4)
            t^2
            sage: sqrt(4*t^4/(1+t)^8)
            2*t^2/(t^4 + 4*t^3 + 6*t^2 + 4*t + 1)

            sage: R.<t> = FpT(GF(5)['t'])
            sage: [a for a in R.iter(2) if (a^2).sqrt() not in (a,-a)]
            []
            sage: [a for a in R.iter(2) if a.is_square() and a.sqrt()^2 != a]
            []
        """
        if nmod_poly_is_zero(self._numer):
            return self

        if not nmod_poly_sqrt_check(self._numer) or not nmod_poly_sqrt_check(self._denom):
            return None

        cdef nmod_poly_t numer
        cdef nmod_poly_t denom
        cdef long a
        cdef FpTElement res

        nmod_poly_init(denom, self.p)
        nmod_poly_init(numer, self.p)

        if nmod_poly_sqrt(numer, self._numer) and nmod_poly_sqrt(denom, self._denom):
            # Make denominator monic
            a = nmod_poly_leading(denom)
            if a != 1:
                a = mod_inverse_int(a, self.p)
                nmod_poly_scalar_mul_nmod(numer, numer, a)
                nmod_poly_scalar_mul_nmod(denom, denom, a)
            # Choose numerator with smaller leading coefficient
            a = nmod_poly_leading(numer)
            if a > self.p - a:
                nmod_poly_neg(numer, numer)
            res = self._new_c()
            nmod_poly_swap(numer, res._numer)
            nmod_poly_swap(denom, res._denom)
            return res
        else:
            nmod_poly_clear(numer)
            nmod_poly_clear(denom)
            return None

    cpdef bint is_square(self) noexcept:
        """
        Return ``True`` if this element is the square of another element of the fraction field.

        EXAMPLES::

            sage: K = GF(13)['t'].fraction_field(); t = K.gen()
            sage: t.is_square()
            False
            sage: (1/t^2).is_square()
            True
            sage: K(0).is_square()
            True
        """
        return self._sqrt_or_None() is not None

    def sqrt(self, extend=True, all=False):
        """
        Return the square root of this element.

        INPUT:

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in an extension ring, if necessary. Otherwise, raise a
          :exc:`ValueError` if the square is not in the base ring.

        - ``all`` -- boolean (default: ``False``); if ``True``, return all
          square roots of self, instead of just one

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: K = GF(7)['t'].fraction_field(); t = K.gen(0)
            sage: p = (t + 2)^2/(3*t^3 + 1)^4
            sage: p.sqrt()
            (3*t + 6)/(t^6 + 3*t^3 + 4)
            sage: p.sqrt()^2 == p
            True
        """
        s = self._sqrt_or_None()
        if s is None:
            if extend:
                raise NotImplementedError("function fields not yet implemented")
            else:
                raise ValueError("not a perfect square")
        else:
            if all:
                if not s:
                    return [s]
                else:
                    return [s, -s]
            else:
                return s

    def __pow__(FpTElement self, Py_ssize_t e, dummy):
        r"""
        Return the `e`-th power of this element.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: R.<t> = FpT(GF(7)['t'])
            sage: t^5
            t^5
            sage: t^-5
            1/t^5

            sage: a = (t+1)/(t-1); a
            (t + 1)/(t + 6)
            sage: a^5
            (t^5 + 5*t^4 + 3*t^3 + 3*t^2 + 5*t + 1)/(t^5 + 2*t^4 + 3*t^3 + 4*t^2 + 5*t + 6)
            sage: a^7
            (t^7 + 1)/(t^7 + 6)
            sage: a^14
            (t^14 + 2*t^7 + 1)/(t^14 + 5*t^7 + 1)

            sage: (a^2)^2 == a^4
            True
            sage: a^3 * a^2 == a^5
            True
            sage: a^47 * a^92 == a^(47+92)
            True
        """
        cdef long a
        assert dummy is None
        cdef FpTElement x = self._new_c()
        if e >= 0:
            nmod_poly_pow(x._numer, self._numer, e)
            nmod_poly_pow(x._denom, self._denom, e)
        else:
            nmod_poly_pow(x._denom, self._numer, -e)
            nmod_poly_pow(x._numer, self._denom, -e)
            if nmod_poly_leading(x._denom) != 1:
                a = mod_inverse_int(nmod_poly_leading(x._denom), self.p)
                nmod_poly_scalar_mul_nmod(x._numer, x._numer, a)
                nmod_poly_scalar_mul_nmod(x._denom, x._denom, a)
        return x


cdef class FpT_iter:
    """
    Return a class that iterates over all elements of an FpT.

    EXAMPLES::

        sage: K = GF(3)['t'].fraction_field()
        sage: I = K.iter(1)
        sage: list(I)
        [0,
         1,
         2,
         t,
         t + 1,
         t + 2,
         2*t,
         2*t + 1,
         2*t + 2,
         1/t,
         2/t,
         (t + 1)/t,
         (t + 2)/t,
         (2*t + 1)/t,
         (2*t + 2)/t,
         1/(t + 1),
         2/(t + 1),
         t/(t + 1),
         (t + 2)/(t + 1),
         2*t/(t + 1),
         (2*t + 1)/(t + 1),
         1/(t + 2),
         2/(t + 2),
         t/(t + 2),
         (t + 1)/(t + 2),
         2*t/(t + 2),
         (2*t + 2)/(t + 2)]
    """
    def __init__(self, parent, degree=None, FpTElement start=None):
        """
        INPUT:

        - ``parent`` -- the FpT that we're iterating over

        - ``degree`` -- the maximum degree of the numerator and denominator of
          the elements over which we iterate

        - ``start`` -- (default: 0) the element on which to start

        EXAMPLES::

            sage: K = GF(11)['t'].fraction_field()
            sage: I = K.iter(2) # indirect doctest
            sage: for a in I:
            ....:     if a.denom()[0] == 3 and a.numer()[1] == 2:
            ....:         print(a); break
            2*t/(t + 3)
        """
        # if degree is None:
        #     raise NotImplementedError
        self.parent = parent
        self.cur = start
        self.degree = -2 if degree is None else degree

    def __cinit__(self, parent, *args, **kwds):
        """
        Memory allocation for the temp variable storing the gcd of the numerator and denominator.

        TESTS::

            sage: from sage.rings.fraction_field_FpT import FpT_iter
            sage: K = GF(7)['t'].fraction_field()
            sage: I = FpT_iter(K, 3) # indirect doctest
            sage: I
            <sage.rings.fraction_field_FpT.FpT_iter object at ...>
        """
        nmod_poly_init(self.g, parent.characteristic())

    def __dealloc__(self):
        """
        Deallocating of self.g.

        TESTS::

            sage: from sage.rings.fraction_field_FpT import FpT_iter
            sage: K = GF(7)['t'].fraction_field()
            sage: I = FpT_iter(K, 3)
            sage: del I # indirect doctest
        """
        nmod_poly_clear(self.g)

    def __iter__(self):
        """
        Return this iterator.

        TESTS::

            sage: from sage.rings.fraction_field_FpT import FpT_iter
            sage: K = GF(3)['t'].fraction_field()
            sage: I = FpT_iter(K, 3)
            sage: for a in I: # indirect doctest
            ....:     if a.numer()[1] == 1 and a.denom()[1] == 2 and a.is_square():
            ....:          print(a); break
            (t^2 + t + 1)/(t^2 + 2*t + 1)
        """
        return self

    def __next__(self):
        """
        Return the next element to iterate over.

        This is achieved by iterating over monic denominators, and for each denominator,
        iterating over all numerators relatively prime to the given denominator.

        EXAMPLES::

            sage: from sage.rings.fraction_field_FpT import *
            sage: K.<t> = FpT(GF(3)['t'])
            sage: list(K.iter(1)) # indirect doctest
            [0,
             1,
             2,
             t,
             t + 1,
             t + 2,
             2*t,
             2*t + 1,
             2*t + 2,
             1/t,
             2/t,
             (t + 1)/t,
             (t + 2)/t,
             (2*t + 1)/t,
             (2*t + 2)/t,
             1/(t + 1),
             2/(t + 1),
             t/(t + 1),
             (t + 2)/(t + 1),
             2*t/(t + 1),
             (2*t + 1)/(t + 1),
             1/(t + 2),
             2/(t + 2),
             t/(t + 2),
             (t + 1)/(t + 2),
             2*t/(t + 2),
             (2*t + 2)/(t + 2)]

            sage: len(list(K.iter(3)))
            2187

            sage: K.<t> = FpT(GF(5)['t'])
            sage: L = list(K.iter(3)); len(L)
            78125
            sage: L[:10]
            [0, 1, 2, 3, 4, t, t + 1, t + 2, t + 3, t + 4]
            sage: L[2000]
            (3*t^3 + 3*t^2 + 3*t + 4)/(t + 2)
            sage: L[-1]
            (4*t^3 + 4*t^2 + 4*t + 4)/(t^3 + 4*t^2 + 4*t + 4)
        """
        cdef FpTElement next_
        if self.cur is None:
            self.cur = self.parent(0)
        elif self.degree == -2:
            self.cur = next(self.cur)
        else:
            next_ = self.cur._copy_c()
            sig_on()
            while True:
                nmod_poly_inc(next_._numer, False)
                if nmod_poly_degree(next_._numer) > self.degree:
                    nmod_poly_inc(next_._denom, True)
                    if nmod_poly_degree(next_._denom) > self.degree:
                        sig_off()
                        raise StopIteration
                    nmod_poly_zero(next_._numer)
                    nmod_poly_set_coeff_ui(next_._numer, 0, 1)
                nmod_poly_gcd(self.g, next_._numer, next_._denom)
                if nmod_poly_is_one(self.g):
                    break
            sig_off()
            self.cur = next_
        return self.cur

cdef class Polyring_FpT_coerce(RingHomomorphism):
    """
    This class represents the coercion map from GF(p)[t] to GF(p)(t).

    EXAMPLES::

        sage: R.<t> = GF(5)[]
        sage: K = R.fraction_field()
        sage: f = K.coerce_map_from(R); f
        Ring morphism:
          From: Univariate Polynomial Ring in t over Finite Field of size 5
          To:   Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
        sage: type(f)
        <class 'sage.rings.fraction_field_FpT.Polyring_FpT_coerce'>

    TESTS::

        TestSuite(f).run()
    """
    cdef long p

    def __init__(self, R):
        """
        INPUT:

        - ``R`` -- an FpT

        EXAMPLES::

            sage: R.<t> = GF(next_prime(2000))[]
            sage: K = R.fraction_field() # indirect doctest
        """
        RingHomomorphism.__init__(self, R.ring_of_integers().Hom(R))
        self.p = R.base_ring().characteristic()

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R) # indirect doctest
            sage: f(t^2 + 1)
            t^2 + 1
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['p'] = self.p
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R) # indirect doctest
            sage: f(t^2 + 1)
            t^2 + 1
        """
        self.p = _slots['p']
        RingHomomorphism._update_slots(self, _slots)

    cpdef Element _call_(self, _x):
        """
        Applies the coercion.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(t^2 + 1) # indirect doctest
            t^2 + 1
        """
        cdef Polynomial_zmod_flint x = <Polynomial_zmod_flint?> _x
        cdef FpTElement ans = <FpTElement>FpTElement.__new__(FpTElement)
        ans._parent = self.codomain()
        ans.p = self.p
        nmod_poly_init(ans._numer, ans.p)
        nmod_poly_init(ans._denom, ans.p)
        nmod_poly_set(ans._numer, &x.x)
        nmod_poly_set_coeff_ui(ans._denom, 0, 1)
        ans.initialized = True
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function allows the map to take multiple arguments,
        usually used to specify both numerator and denominator.

        If ``reduce`` is specified as False, then the result won't be
        normalized.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(2*t + 2, t + 3) # indirect doctest
            (2*t + 2)/(t + 3)
            sage: f(2*t + 2, 2)
            t + 1
            sage: f(2*t + 2, 2, reduce=False)
            (2*t + 2)/2

        TESTS:

        Check that :issue:`12217` and :issue:`16811` are fixed::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(t, 0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: fraction has denominator 0
            sage: f(t, GF(5).zero())
            Traceback (most recent call last):
            ...
            ZeroDivisionError: fraction has denominator 0
            sage: f(t, R.zero())
            Traceback (most recent call last):
            ...
            ZeroDivisionError: fraction has denominator 0
        """
        cdef Polynomial_zmod_flint x
        cdef unsigned long r
        try:
            x = <Polynomial_zmod_flint?> _x
        except TypeError:
            raise NotImplementedError('Fraction fields not implemented for this type.')
        cdef FpTElement ans = <FpTElement>FpTElement.__new__(FpTElement)
        ans._parent = self.codomain()
        ans.p = self.p
        nmod_poly_init(ans._numer, ans.p)
        nmod_poly_init(ans._denom, ans.p)
        nmod_poly_set(ans._numer, &x.x)
        if len(args) == 0:
            nmod_poly_set_coeff_ui(ans._denom, 0, 1)  # No need to normalize
        elif len(args) == 1:
            y = args[0]
            if isinstance(y, Integer):
                r = mpz_fdiv_ui((<Integer>y).value, self.p)
                nmod_poly_set_coeff_ui(ans._denom, 0, r)
            else:
                # could use the coerce keyword being set to False to not check this...
                if not (isinstance(y, Element) and y.parent() is self.domain()):
                    # We could special case integers and GF(p) elements here.
                    y = self.domain()(y)
                nmod_poly_set(ans._denom, &((<Polynomial_zmod_flint?>y).x))
            # Normalize the fraction, checking for division by zero
            if nmod_poly_is_zero(ans._denom):
                raise ZeroDivisionError('fraction has denominator 0')
            if kwds.get('reduce', True):
                normalize(ans._numer, ans._denom, ans.p)
        else:
            raise TypeError("FpT only supports two positional arguments")
        ans.initialized = True
        return ans

    def section(self):
        """
        Return the section of this inclusion: the partially defined map from ``GF(p)(t)``
        back to ``GF(p)[t]``, defined on elements with unit denominator.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = f.section(); g
            Section map:
              From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
              To:   Univariate Polynomial Ring in t over Finite Field of size 5
            sage: t = K.gen()
            sage: g(t)
            t
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
        """
        return FpT_Polyring_section(self)

cdef class FpT_Polyring_section(Section):
    """
    This class represents the section from GF(p)(t) back to GF(p)[t].

    EXAMPLES::

        sage: R.<t> = GF(5)[]
        sage: K = R.fraction_field()
        sage: f = R.convert_map_from(K); f
        Section map:
          From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
          To:   Univariate Polynomial Ring in t over Finite Field of size 5
        sage: type(f)
        <class 'sage.rings.fraction_field_FpT.FpT_Polyring_section'>

    .. WARNING::

        Comparison of ``FpT_Polyring_section`` objects is not currently
        implemented. See :issue:`23469`. ::

            sage: fprime = loads(dumps(f))
            sage: fprime == f
            False

            sage: fprime(1+t) == f(1+t)
            True

    TESTS::

        sage: TestSuite(f).run(skip='_test_pickling')
    """
    cdef long p

    def __init__(self, Polyring_FpT_coerce f):
        """
        INPUT:

        - ``f`` -- a Polyring_FpT_coerce homomorphism

        EXAMPLES::

            sage: R.<t> = GF(next_prime(2000))[]
            sage: K = R.fraction_field()
            sage: R(K.gen(0)) # indirect doctest
            t
        """
        self.p = f.p
        Section.__init__(self, f)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = f.section()   # indirect doctest
            sage: t = K.gen()
            sage: g(t^2)
            t^2
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
        """
        slots = Section._extra_slots(self)
        slots['p'] = self.p
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = f.section()   # indirect doctest
            sage: t = K.gen()
            sage: g(t^2)
            t^2
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
        """
        self.p = _slots['p']
        Section._update_slots(self, _slots)

    cpdef Element _call_(self, _x):
        """
        Applies the section.

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = f.section(); g
            Section map:
              From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7
              To:   Univariate Polynomial Ring in t over Finite Field of size 7
            sage: t = K.gen()
            sage: g(t^2) # indirect doctest
            t^2
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
        """
        cdef FpTElement x = <FpTElement?>_x
        cdef Polynomial_zmod_flint ans
        if nmod_poly_degree(x._denom) != 0:
            normalize(x._numer, x._denom, self.p)
            if nmod_poly_degree(x._denom) != 0:
                raise ValueError("not integral")
        ans = Polynomial_zmod_flint.__new__(Polynomial_zmod_flint)
        if nmod_poly_get_coeff_ui(x._denom, 0) != 1:
            normalize(x._numer, x._denom, self.p)
        nmod_poly_init(&ans.x, self.p)
        nmod_poly_set(&ans.x, x._numer)
        ans._parent = self.codomain()
        ans._cparent = get_cparent(ans._parent)
        return ans

cdef class Fp_FpT_coerce(RingHomomorphism):
    """
    This class represents the coercion map from GF(p) to GF(p)(t).

    EXAMPLES::

        sage: R.<t> = GF(5)[]
        sage: K = R.fraction_field()
        sage: f = K.coerce_map_from(GF(5)); f
        Ring morphism:
          From: Finite Field of size 5
          To:   Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
        sage: type(f)
        <class 'sage.rings.fraction_field_FpT.Fp_FpT_coerce'>

    TESTS::

        sage: TestSuite(f).run()
    """
    cdef long p

    def __init__(self, R):
        """
        INPUT:

        - ``R`` -- an FpT

        EXAMPLES::

            sage: R.<t> = GF(next_prime(3000))[]
            sage: K = R.fraction_field() # indirect doctest
        """
        RingHomomorphism.__init__(self, R.base_ring().Hom(R))
        self.p = R.base_ring().characteristic()

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(5))
            sage: g = copy(f)
            sage: g == f
            True
            sage: g(GF(5)(2)) == f(GF(5)(2))
            True
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['p'] = self.p
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(5))
            sage: g = copy(f)
            sage: g == f
            True
            sage: g(GF(5)(2)) == f(GF(5)(2))
            True
        """
        self.p = _slots['p']
        RingHomomorphism._update_slots(self, _slots)

    cpdef Element _call_(self, _x):
        """
        Applies the coercion.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(5))
            sage: f(GF(5)(3)) # indirect doctest
            3
        """
        cdef IntegerMod_int x = <IntegerMod_int?> _x
        cdef FpTElement ans = <FpTElement>FpTElement.__new__(FpTElement)
        ans._parent = self.codomain()
        ans.p = self.p
        nmod_poly_init(ans._numer, ans.p)
        nmod_poly_init(ans._denom, ans.p)
        nmod_poly_set_coeff_ui(ans._numer, 0, x.ivalue)
        nmod_poly_set_coeff_ui(ans._denom, 0, 1)
        ans.initialized = True
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function allows the map to take multiple arguments, usually used to specify both numerator and denominator.

        If ``reduce`` is specified as False, then the result won't be normalized.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(5))
            sage: f(1, t + 3) # indirect doctest
            1/(t + 3)
            sage: f(2, 2*t)
            1/t
            sage: f(2, 2*t, reduce=False)
            2/2*t
        """
        cdef IntegerMod_int x = <IntegerMod_int?> _x
        cdef FpTElement ans = <FpTElement>FpTElement.__new__(FpTElement)
        ans._parent = self.codomain()
        ans.p = self.p
        nmod_poly_init(ans._numer, ans.p)
        nmod_poly_init(ans._denom, ans.p)
        cdef long r
        nmod_poly_set_coeff_ui(ans._numer, 0, x.ivalue)
        if len(args) == 0:
            nmod_poly_set_coeff_ui(ans._denom, 0, 1)
        if len(args) == 1:
            y = args[0]
            if isinstance(y, Integer):
                r = mpz_fdiv_ui((<Integer>y).value, self.p)
                if r == 0:
                    raise ZeroDivisionError
                nmod_poly_set_coeff_ui(ans._denom, 0, r)
            else:
                R = ans._parent.ring_of_integers()
                # could use the coerce keyword being set to False to not check this...
                if not (isinstance(y, Element) and y.parent() is R):
                    # We could special case integers and GF(p) elements here.
                    y = R(y)
                nmod_poly_set(ans._denom, &((<Polynomial_zmod_flint?>y).x))
        else:
            raise ValueError("FpT only supports two positional arguments")
        if 'reduce' not in kwds or kwds['reduce']:
            normalize(ans._numer, ans._denom, ans.p)
        ans.initialized = True
        return ans

    def section(self):
        """
        Return the section of this inclusion: the partially defined map from ``GF(p)(t)``
        back to ``GF(p)``, defined on constant elements.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(5))
            sage: g = f.section(); g
            Section map:
              From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
              To:   Finite Field of size 5
            sage: t = K.gen()
            sage: g(f(1,3,reduce=False))
            2
            sage: g(t)
            Traceback (most recent call last):
            ...
            ValueError: not constant
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
        """
        return FpT_Fp_section(self)

cdef class FpT_Fp_section(Section):
    """
    This class represents the section from GF(p)(t) back to GF(p)[t].

    EXAMPLES::

        sage: R.<t> = GF(5)[]
        sage: K = R.fraction_field()
        sage: f = GF(5).convert_map_from(K); f
        Section map:
          From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
          To:   Finite Field of size 5
        sage: type(f)
        <class 'sage.rings.fraction_field_FpT.FpT_Fp_section'>

    .. WARNING::

        Comparison of ``FpT_Fp_section`` objects is not currently
        implemented. See :issue:`23469`. ::

            sage: fprime = loads(dumps(f))
            sage: fprime == f
            False

            sage: fprime(3) == f(3)
            True

    TESTS::

        sage: TestSuite(f).run(skip='_test_pickling')
    """
    cdef long p

    def __init__(self, Fp_FpT_coerce f):
        """
        INPUT:

        - ``f`` -- an ``Fp_FpT_coerce`` homomorphism

        EXAMPLES::

            sage: R.<t> = GF(next_prime(2000))[]
            sage: K = R.fraction_field()
            sage: GF(next_prime(2000))(K(127)) # indirect doctest
            127
        """
        self.p = f.p
        Section.__init__(self, f)

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(7))
            sage: g = f.section()   # indirect doctest
            sage: t = K.gen()
            sage: g(t^2)
            Traceback (most recent call last):
            ...
            ValueError: not constant
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
            sage: g(K(4))
            4
            sage: g(K(0))
            0
        """
        slots = Section._extra_slots(self)
        slots['p'] = self.p
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(7))
            sage: g = f.section()   # indirect doctest
            sage: t = K.gen()
            sage: g(t^2)
            Traceback (most recent call last):
            ...
            ValueError: not constant
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
            sage: g(K(4))
            4
            sage: g(K(0))
            0
        """
        self.p = _slots['p']
        Section._update_slots(self, _slots)

    cpdef Element _call_(self, _x):
        """
        Applies the section.

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(GF(7))
            sage: g = f.section(); g
            Section map:
              From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 7
              To:   Finite Field of size 7
            sage: t = K.gen()
            sage: g(t^2) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: not constant
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
            sage: g(K(4))
            4
            sage: g(K(0))
            0
        """
        cdef FpTElement x = <FpTElement?>_x
        cdef IntegerMod_int ans
        if nmod_poly_degree(x._denom) != 0 or nmod_poly_degree(x._numer) > 0:
            normalize(x._numer, x._denom, self.p)
            if nmod_poly_degree(x._denom) != 0:
                raise ValueError("not integral")
            if nmod_poly_degree(x._numer) > 0:
                raise ValueError("not constant")
        ans = IntegerMod_int.__new__(IntegerMod_int)
        ans._parent = self.codomain()
        ans._modulus = ans._parent._pyx_order
        if nmod_poly_get_coeff_ui(x._denom, 0) != 1:
            normalize(x._numer, x._denom, self.p)
        ans.ivalue = nmod_poly_get_coeff_ui(x._numer, 0)
        return ans

cdef class ZZ_FpT_coerce(RingHomomorphism):
    """
    This class represents the coercion map from ZZ to GF(p)(t).

    EXAMPLES::

        sage: R.<t> = GF(17)[]
        sage: K = R.fraction_field()
        sage: f = K.coerce_map_from(ZZ); f
        Ring morphism:
          From: Integer Ring
          To:   Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 17
        sage: type(f)
        <class 'sage.rings.fraction_field_FpT.ZZ_FpT_coerce'>

    TESTS::

        sage: TestSuite(f).run()
    """
    cdef long p

    def __init__(self, R):
        """
        INPUT:

        - ``R`` -- an FpT

        EXAMPLES::

            sage: R.<t> = GF(next_prime(3000))[]
            sage: K = R.fraction_field() # indirect doctest
        """
        RingHomomorphism.__init__(self, ZZ.Hom(R))
        self.p = R.base_ring().characteristic()

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g(5) == f(5)
            True
            sage: g(0) == f(0)
            True
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['p'] = self.p
        return slots

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g(5) == f(5)
            True
            sage: g(0) == f(0)
            True
        """
        self.p = _slots['p']
        RingHomomorphism._update_slots(self, _slots)

    cpdef Element _call_(self, _x):
        """
        Applies the coercion.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(ZZ)
            sage: f(3) # indirect doctest
            3
        """
        cdef Integer x = <Integer?> _x
        cdef FpTElement ans = <FpTElement>FpTElement.__new__(FpTElement)
        ans._parent = self.codomain()
        ans.p = self.p
        nmod_poly_init(ans._numer, ans.p)
        nmod_poly_init(ans._denom, ans.p)
        nmod_poly_set_coeff_ui(ans._numer, 0, mpz_fdiv_ui(x.value, self.p))
        nmod_poly_set_coeff_ui(ans._denom, 0, 1)
        ans.initialized = True
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function allows the map to take multiple arguments, usually used to specify both numerator and denominator.

        If ``reduce`` is specified as False, then the result won't be normalized.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(ZZ)
            sage: f(1, t + 3) # indirect doctest
            1/(t + 3)
            sage: f(1,2)
            3
            sage: f(2, 2*t)
            1/t
            sage: f(2, 2*t, reduce=False)
            2/2*t
        """
        cdef Integer x = <Integer?> _x
        cdef FpTElement ans = <FpTElement>FpTElement.__new__(FpTElement)
        ans._parent = self.codomain()
        ans.p = self.p
        nmod_poly_init(ans._numer, ans.p)
        nmod_poly_init(ans._denom, ans.p)
        cdef long r
        nmod_poly_set_coeff_ui(ans._numer, 0, mpz_fdiv_ui(x.value, self.p))
        if len(args) == 0:
            nmod_poly_set_coeff_ui(ans._denom, 0, 1)
        if len(args) == 1:
            y = args[0]
            if isinstance(y, Integer):
                r = mpz_fdiv_ui((<Integer>y).value, self.p)
                if r == 0:
                    raise ZeroDivisionError
                nmod_poly_set_coeff_ui(ans._denom, 0, r)
            else:
                R = ans._parent.ring_of_integers()
                # could use the coerce keyword being set to False to not check this...
                if not (isinstance(y, Element) and y.parent() is R):
                    # We could special case integers and GF(p) elements here.
                    y = R(y)
                nmod_poly_set(ans._denom, &((<Polynomial_zmod_flint?>y).x))
        else:
            raise ValueError("FpT only supports two positional arguments")
        if 'reduce' not in kwds or kwds['reduce']:
            normalize(ans._numer, ans._denom, ans.p)
        ans.initialized = True
        return ans

    def section(self):
        """
        Return the section of this inclusion: the partially defined map from ``GF(p)(t)``
        back to ``ZZ``, defined on constant elements.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(ZZ)
            sage: g = f.section(); g
            Composite map:
              From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
              To:   Integer Ring
              Defn:   Section map:
                      From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
                      To:   Finite Field of size 5
                    then
                      Lifting map:
                      From: Finite Field of size 5
                      To:   Integer Ring
            sage: t = K.gen()
            sage: g(f(1,3,reduce=False))
            2
            sage: g(t)
            Traceback (most recent call last):
            ...
            ValueError: not constant
            sage: g(1/t)
            Traceback (most recent call last):
            ...
            ValueError: not integral
        """
        return ZZ.convert_map_from(self.codomain().base_ring()) * Fp_FpT_coerce(self.codomain()).section()

cdef inline bint normalize(nmod_poly_t numer, nmod_poly_t denom, long p) noexcept:
    """
    Put ``numer`` / ``denom`` into a normal form: denominator monic and sharing no common factor with the numerator.

    The normalized form of 0 is 0/1.

    Return ``True`` if ``numer`` and ``denom`` were changed.
    """
    cdef long a
    cdef bint changed
    if nmod_poly_degree(numer) == -1:
        if nmod_poly_degree(denom) > 0 or nmod_poly_leading(denom) != 1:
            changed = True
        else:
            changed = False
        nmod_poly_truncate(denom, 0)
        nmod_poly_set_coeff_ui(denom, 0, 1)
        return changed
    elif nmod_poly_degree(numer) == 0 or nmod_poly_degree(denom) == 0:
        if nmod_poly_leading(denom) != 1:
            a = mod_inverse_int(nmod_poly_leading(denom), p)
            nmod_poly_scalar_mul_nmod(numer, numer, a)
            nmod_poly_scalar_mul_nmod(denom, denom, a)
            return True
        return False
    cdef nmod_poly_t g
    changed = False
    try:
        nmod_poly_init_preinv(g, p, numer.mod.ninv)
        nmod_poly_gcd(g, numer, denom)
        if nmod_poly_degree(g) != 0:
            # Divide knowing divisible by? Can we get these quotients as a byproduct of the gcd?
            nmod_poly_div(numer, numer, g)
            nmod_poly_div(denom, denom, g)
            changed = True
        if nmod_poly_leading(denom) != 1:
            a = mod_inverse_int(nmod_poly_leading(denom), p)
            nmod_poly_scalar_mul_nmod(numer, numer, a)
            nmod_poly_scalar_mul_nmod(denom, denom, a)
            changed = True
        return changed
    finally:
        nmod_poly_clear(g)


cdef inline unsigned long nmod_poly_leading(nmod_poly_t poly) noexcept:
    """
    Return the leading coefficient of ``poly``.
    """
    return nmod_poly_get_coeff_ui(poly, nmod_poly_degree(poly))


cdef inline void nmod_poly_inc(nmod_poly_t poly, bint monic) noexcept:
    """
    Set poly to the "next" polynomial: this is just counting in base p.

    If monic is ``True`` then will only iterate through monic polynomials.
    """
    cdef long n
    cdef long a
    cdef long p = poly.mod.n
    for n from 0 <= n <= nmod_poly_degree(poly) + 1:
        a = nmod_poly_get_coeff_ui(poly, n) + 1
        if a == p:
            nmod_poly_set_coeff_ui(poly, n, 0)
        else:
            nmod_poly_set_coeff_ui(poly, n, a)
            break
    if monic and a == 2 and n == nmod_poly_degree(poly):
        nmod_poly_set_coeff_ui(poly, n, 0)
        nmod_poly_set_coeff_ui(poly, n + 1, 1)


cdef inline long nmod_poly_cmp(nmod_poly_t a, nmod_poly_t b) noexcept:
    """
    Compare `a` and `b`, returning 0 if they are equal.

    - If the degree of `a` is less than that of `b`, returns `-1`.

    - If the degree of `b` is less than that of `a`, returns `1`.

    - Otherwise, compares `a` and `b` lexicographically, starting at the leading terms.
    """
    cdef long ad = nmod_poly_degree(a)
    cdef long bd = nmod_poly_degree(b)
    if ad < bd:
        return -1
    elif ad > bd:
        return 1
    cdef long d = nmod_poly_degree(a)
    while d >= 0:
        ad = nmod_poly_get_coeff_ui(a, d)
        bd = nmod_poly_get_coeff_ui(b, d)
        if ad < bd:
            return -1
        elif ad > bd:
            return 1
        d -= 1
    return 0


cdef bint nmod_poly_sqrt_check(nmod_poly_t poly) noexcept:
    """
    Quick check to see if ``poly`` could possibly be a square.
    """
    # We could use Sage's jacobi_int which is for 32 bits integers rather
    # than FLINT's n_jacobi which is for longs as the FpT class is crafted
    # for primes 2 < p < 2^16
    return (nmod_poly_degree(poly) % 2 == 0
            and n_jacobi(nmod_poly_leading(poly), poly.mod.n) == 1
            and n_jacobi(nmod_poly_get_coeff_ui(poly, 0), poly.mod.n) != -1)


def unpickle_FpT_element(K, numer, denom):
    """
    Used for pickling.

    TESTS::

        sage: from sage.rings.fraction_field_FpT import unpickle_FpT_element
        sage: R.<t> = GF(13)['t']
        sage: unpickle_FpT_element(Frac(R), t+1, t)
        (t + 1)/t
    """
    return FpTElement(K, numer, denom, coerce=False, reduce=False)


#  Somehow this isn't in FLINT, evidently.  It could be moved
#  elsewhere at some point.
cdef int sage_cmp_nmod_poly_t(nmod_poly_t L, nmod_poly_t R) noexcept:
    """
    Compare two ``nmod_poly_t`` in a Pythonic way, so this returns `-1`, `0`,
    or `1`, and is consistent.
    """
    cdef int j
    cdef Py_ssize_t i

    # First compare the degrees
    j = nmod_poly_degree(L) - nmod_poly_degree(R)
    if j < 0:
        return -1
    elif j > 0:
        return 1

    # Same degree, so compare coefficients, term by term
    for i in range(nmod_poly_degree(L) + 1):
        j = nmod_poly_get_coeff_ui(L, i) - nmod_poly_get_coeff_ui(R, i)
        if j < 0:
            return -1
        elif j > 0:
            return 1

    # Two polynomials are equal
    return 0
