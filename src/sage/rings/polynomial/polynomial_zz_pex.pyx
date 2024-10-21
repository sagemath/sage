# sage_setup: distribution = sagemath-ntl
# sage.doctest: needs sage.libs.ntl sage.rings.finite_rings
# distutils: libraries = NTL_LIBRARIES gmp
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
Univariate Polynomials over GF(p^e) via NTL's ZZ_pEX

AUTHOR:

- Yann Laigle-Chapuy (2010-01) initial implementation
- Lorenz Panny (2023-01): :meth:`minpoly_mod`
- Giacomo Pope (2023-08): :meth:`reverse`, :meth:`inverse_series_trunc`
"""
from cysignals.signals cimport sig_on, sig_off

from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class
from sage.libs.ntl.ZZ_pE cimport ZZ_pE_to_ZZ_pX
from sage.libs.ntl.ZZ_pX cimport ZZ_pX_deg, ZZ_pX_coeff
from sage.libs.ntl.ZZ_p cimport ZZ_p_rep
from sage.libs.ntl.convert cimport ZZ_to_mpz, mpz_to_ZZ

from sage.structure.element import have_same_parent, canonical_coercion

# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

cdef cparent get_cparent(parent) except? NULL:
    if parent is None:
        return NULL
    cdef ntl_ZZ_pEContext_class pec
    try:
        pec = parent._modulus
    except AttributeError:
        return NULL
    return &(pec.ptrs)

# first we include the definitions
include "sage/libs/ntl/ntl_ZZ_pEX_linkage.pxi"

# and then the interface
include "polynomial_template.pxi"

from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE

cdef inline ZZ_pE_c_to_list(ZZ_pE_c x):
    cdef list L = []
    cdef ZZ_pX_c c_pX
    cdef ZZ_p_c c_p
    cdef ZZ_c c_c
    cdef Integer ans

    c_pX = ZZ_pE_to_ZZ_pX(x)
    d = ZZ_pX_deg(c_pX)
    if d>=0:
        for 0 <= j <= d:
            c_p = ZZ_pX_coeff(c_pX, j)
            c_c = ZZ_p_rep(c_p)
            ans = Integer.__new__(Integer)
            ZZ_to_mpz(ans.value, &c_c)
            L.append(ans)
    return L


cdef class Polynomial_ZZ_pEX(Polynomial_template):
    r"""
    Univariate Polynomials over `\GF{p^n}` via NTL's ``ZZ_pEX``.

    EXAMPLES::

        sage: K.<a> = GF(next_prime(2**60)**3)
        sage: R.<x> = PolynomialRing(K, implementation='NTL')
        sage: (x^3 + a*x^2 + 1) * (x + a)
        x^4 + 2*a*x^3 + a^2*x^2 + x + a
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        r"""
        Create a new univariate polynomials over `\GF{p^n}`.

        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: x^2+a
            x^2 + a

        TESTS:

        The following tests against a bug that was fixed in :issue:`9944`.
        With the ring definition above, we now have::

            sage: R([3,'1234'])
            1234*x + 3
            sage: R([3,'12e34'])
            Traceback (most recent call last):
            ...
            TypeError: unable to convert '12e34' to an integer
            sage: R([3,x])
            Traceback (most recent call last):
            ...
            TypeError: x is not a constant polynomial

        Check that NTL contexts are correctly restored and that
        :issue:`9524` has been fixed::

            sage: x = polygen(GF(9, 'a'))
            sage: x = polygen(GF(49, 'a'))
            sage: -x
            6*x
            sage: 5*x
            5*x

        Check that :issue:`11239` is fixed::

            sage: Fq.<a> = GF(2^4); Fqq.<b> = GF(3^7)
            sage: PFq.<x> = Fq[]; PFqq.<y> = Fqq[]
            sage: f = x^3 + (a^3 + 1)*x
            sage: sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX(PFqq, f)
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce from a finite field other than the prime subfield
        """
        cdef ntl_ZZ_pE d
        try:
            if (x.parent() is parent.base_ring()) or (x.parent() == parent.base_ring()):
                Polynomial.__init__(self, parent, is_gen=is_gen)
                (<Polynomial_template>self)._cparent = get_cparent(parent)
                celement_construct(&self.x, (<Polynomial_template>self)._cparent)
                d = parent._modulus.ZZ_pE(list(x.polynomial()))
                ZZ_pEX_SetCoeff(self.x, 0, d.x)
                return
        except AttributeError:
            pass

        if isinstance(x, Polynomial):
            x = x.list()

        if isinstance(x, (list, tuple)):
            Polynomial.__init__(self, parent, is_gen=is_gen)
            (<Polynomial_template>self)._cparent = get_cparent(parent)
            celement_construct(&self.x, (<Polynomial_template>self)._cparent)
            K = parent.base_ring()
            for i,e in enumerate(x):
                # self(x) is supposed to be a conversion,
                # not necessarily a coercion. So, we must
                # not do K.coerce(e) but K(e).
                e = K(e)
                d = parent._modulus.ZZ_pE(list(e.polynomial()))
                ZZ_pEX_SetCoeff(self.x, i, d.x)
            return

        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    cdef get_unsafe(self, Py_ssize_t i):
        r"""
        Return the `i`-th coefficient of ``self``.

        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: f = x^3 + (2*a+1)*x + a
            sage: f[0]
            a
            sage: f[1]
            2*a + 1
            sage: f[2]
            0
            sage: f[:2]
            (2*a + 1)*x + a
            sage: f[:50] == f
            True
        """
        self._parent._modulus.restore()
        cdef ZZ_pE_c c_pE = ZZ_pEX_coeff(self.x, i)
        return self._parent._base(ZZ_pE_c_to_list(c_pE))

    cpdef list list(self, bint copy=True):
        r"""
        Return the list of coefficients.

        EXAMPLES::

            sage: K.<a> = GF(5^3)
            sage: P = PolynomialRing(K, 'x')
            sage: f = P.random_element(100)
            sage: f.list() == [f[i] for i in range(f.degree()+1)]
            True
            sage: P.0.list()
            [0, 1]
        """
        cdef Py_ssize_t i

        self._parent._modulus.restore()

        K = self._parent.base_ring()
        return [K(ZZ_pE_c_to_list(ZZ_pEX_coeff(self.x, i)))
                for i in range(celement_len(&self.x, (<Polynomial_template>self)._cparent))]

    cpdef _lmul_(self, Element left):
        r"""
        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: (2*a+1)*x # indirect doctest
            (2*a + 1)*x
            sage: x*(2*a+1) # indirect doctest
            (2*a + 1)*x
        """
        cdef ntl_ZZ_pE d
        cdef Polynomial_ZZ_pEX r
        r = Polynomial_ZZ_pEX.__new__(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        d = self._parent._modulus.ZZ_pE(list(left.polynomial()))
        ZZ_pEX_mul_ZZ_pE(r.x, self.x, d.x)
        return r

    def __call__(self, *x, **kwds):
        r"""
        Evaluate polynomial at `a`.

        EXAMPLES::

            sage: K.<u> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: P = (x-u)*(x+u+1)
            sage: P(u)
            0
            sage: P(u+1)
            2*u + 2

        TESTS:

        The work around provided in :issue:`10475` is superseded by :issue:`24072`::

            sage: F.<x> = GF(4)
            sage: P.<y> = F[]
            sage: p = y^4 + x*y^3 + y^2 + (x + 1)*y + x + 1
            sage: SR(p)                                                                 # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations

        Check that polynomial evaluation works when using logarithmic
        representation of finite field elements (:issue:`16383`)::

            sage: for i in range(10):
            ....:     F = FiniteField(random_prime(15) ** ZZ.random_element(2, 5), 'a', repr='log')
            ....:     b = F.random_element()
            ....:     P = PolynomialRing(F, 'x')
            ....:     f = P.random_element(8)
            ....:     assert f(b) == sum(c * b^i for i, c in enumerate(f))
        """
        cdef ntl_ZZ_pE _a
        cdef ZZ_pE_c c_b

        K = self._parent.base_ring()

        if kwds:
            if x:
                raise TypeError("%s__call__() takes exactly 1 argument"%type(self))
            try:
                x = [kwds.pop(self.variable_name())]
            except KeyError:
                pass
        if kwds:
            raise TypeError("%s__call__() accepts no named argument except '%s'"%(type(self),self.variable_name()))
        if len(x)!=1:
            raise TypeError("%s__call__() takes exactly 1 positional argument"%type(self))

        a = x[0]
        try:
            if a.parent() is not K:
                a = K.coerce(a)
        except (TypeError, AttributeError, NotImplementedError):
            return Polynomial.__call__(self, a)

        _a = self._parent._modulus.ZZ_pE(list(a.polynomial()))
        ZZ_pEX_eval(c_b, self.x, _a.x)
        return K(ZZ_pE_c_to_list(c_b))

    def resultant(self, other):
        r"""
        Return the resultant of ``self`` and ``other``, which must lie in the same
        polynomial ring.

        INPUT:

        - ``other`` -- a polynomial

        OUTPUT: an element of the base ring of the polynomial ring

        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: f = (x-a)*(x-a**2)*(x+1)
            sage: g = (x-a**3)*(x-a**4)*(x+a)
            sage: r = f.resultant(g)
            sage: r == prod(u - v for (u,eu) in f.roots() for (v,ev) in g.roots())
            True
        """
        cdef ZZ_pE_c r
        self._parent._modulus.restore()

        if other.parent() is not self._parent:
            other = self._parent.coerce(other)

        ZZ_pEX_resultant(r, self.x, (<Polynomial_ZZ_pEX>other).x)

        K = self._parent.base_ring()
        return K(K.polynomial_ring()(ZZ_pE_c_to_list(r)))

    def is_irreducible(self, algorithm='fast_when_false', iter=1):
        r"""
        Return ``True`` precisely when ``self`` is irreducible over its base ring.

        INPUT:

        - ``algorithm`` -- string (default: ``'fast_when_false'``);
          there are 3 available algorithms:
          ``'fast_when_true'``, ``'fast_when_false'``, and ``'probabilistic'``

        - ``iter`` -- (default: 1) if the algorithm is ``'probabilistic'``,
          defines the number of iterations. The error probability is bounded
          by `q^{\text{-iter}}` for polynomials in `\GF{q}[x]`.

        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: P = x^3 + (2-a)*x + 1
            sage: P.is_irreducible(algorithm='fast_when_false')
            True
            sage: P.is_irreducible(algorithm='fast_when_true')
            True
            sage: P.is_irreducible(algorithm='probabilistic')
            True
            sage: Q = (x^2+a)*(x+a^3)
            sage: Q.is_irreducible(algorithm='fast_when_false')
            False
            sage: Q.is_irreducible(algorithm='fast_when_true')
            False
            sage: Q.is_irreducible(algorithm='probabilistic')
            False
        """
        self._parent._modulus.restore()
        if algorithm=="fast_when_false":
            sig_on()
            res = ZZ_pEX_IterIrredTest(self.x)
            sig_off()
        elif algorithm=="fast_when_true":
            sig_on()
            res = ZZ_pEX_DetIrredTest(self.x)
            sig_off()
        elif algorithm=="probabilistic":
            sig_on()
            res = ZZ_pEX_ProbIrredTest(self.x, iter)
            sig_off()
        else:
            raise ValueError("unknown algorithm")
        return res != 0

    def minpoly_mod(self, other):
        r"""
        Compute the minimal polynomial of this polynomial modulo another
        polynomial in the same ring.

        ALGORITHM:

        NTL's ``MinPolyMod()``, which uses Shoup's algorithm [Sho1999]_.

        EXAMPLES::

            sage: R.<x> = GF(101^2)[]
            sage: f = x^17 + x^2 - 1
            sage: (x^2).minpoly_mod(f)
            x^17 + 100*x^2 + 2*x + 100

        TESTS:

        Random testing::

            sage: p = random_prime(50)
            sage: e = randrange(2,10)
            sage: R.<x> = GF((p,e),'a')[]
            sage: d = randrange(1,50)
            sage: f = R.random_element(d)
            sage: g = R.random_element((-1,5*d))
            sage: poly = g.minpoly_mod(f)
            sage: poly(R.quotient(f)(g))
            0
        """
        self._parent._modulus.restore()

        if other.parent() is not self._parent:
            other = self._parent.coerce(other)

        cdef Polynomial_ZZ_pEX r
        r = Polynomial_ZZ_pEX.__new__(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent

        ZZ_pEX_MinPolyMod(r.x, (<Polynomial_ZZ_pEX>(self % other)).x, (<Polynomial_ZZ_pEX>other).x)
        return r

    cpdef _richcmp_(self, other, int op):
        r"""
        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: P1 = (a**2+a+1)*x^2 + a*x + 1
            sage: P2 = (     a+1)*x^2 + a*x + 1
            sage: P1 < P2 # indirect doctests
            False

        TESTS::

            sage: P3 = (a**2+a+1)*x^2 + x + 1
            sage: P4 =                  x + 1
            sage: P1 < P3
            False
            sage: P1 < P4
            False
            sage: P1 > P2
            True
            sage: P1 > P3
            True
            sage: P1 > P4
            True
        """
        return Polynomial._richcmp_(self, other, op)

    def shift(self, int n):
        r"""
        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: f = x^3 + x^2 + 1
            sage: f.shift(1)
            x^4 + x^3 + x
            sage: f.shift(-1)
            x^2 + x
        """
        self._parent._modulus.restore()
        cdef Polynomial_ZZ_pEX r
        r = Polynomial_ZZ_pEX.__new__(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        ZZ_pEX_LeftShift(r.x, self.x, n)
        return r

    def __lshift__(self, int n):
        r"""
        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: f = x^3 + x^2 + 1
            sage: f << 1
            x^4 + x^3 + x
            sage: f << -1
            x^2 + x
        """
        return self.shift(n)

    def __rshift__(self, int n):
        r"""
        EXAMPLES::

            sage: K.<a> = GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K, implementation='NTL')
            sage: f = x^3 + x^2 + 1
            sage: f >> 1
            x^2 + x
            sage: f >> -1
            x^4 + x^3 + x
        """
        return self.shift(-n)

    def reverse(self, degree=None):
        r"""
        Return the polynomial obtained by reversing the coefficients
        of this polynomial.  If degree is set then this function behaves
        as if this polynomial has degree ``degree``.

        EXAMPLES::

            sage: R.<x> = GF(101^2)[]
            sage: f = x^13 + 11*x^10 + 32*x^6 + 4
            sage: f.reverse()
            4*x^13 + 32*x^7 + 11*x^3 + 1
            sage: f.reverse(degree=15)
            4*x^15 + 32*x^9 + 11*x^5 + x^2
            sage: f.reverse(degree=2)
            4*x^2

        TESTS::

            sage: R.<x> = GF(163^2)[]
            sage: f = R([p for p in primes(20)])
            sage: f.reverse()
            2*x^7 + 3*x^6 + 5*x^5 + 7*x^4 + 11*x^3 + 13*x^2 + 17*x + 19
            sage: f.reverse(degree=200)
            2*x^200 + 3*x^199 + 5*x^198 + 7*x^197 + 11*x^196 + 13*x^195 + 17*x^194 + 19*x^193
            sage: f.reverse(degree=0)
            2
            sage: f.reverse(degree=-5)
            Traceback (most recent call last):
            ...
            ValueError: degree argument must be a nonnegative integer, got -5

        Check that this implementation is compatible with the generic one::

            sage: p = R([0,1,0,2])
            sage: all(p.reverse(d) == Polynomial.reverse(p, d)
            ....:     for d in [None, 0, 1, 2, 3, 4])
            True
        """
        self._parent._modulus.restore()

        # Construct output polynomial
        cdef Polynomial_ZZ_pEX r
        r = Polynomial_ZZ_pEX.__new__(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent

        # When a degree has been supplied, ensure it is a valid input
        cdef unsigned long d
        if degree is not None:
            if degree < 0:
                raise ValueError("degree argument must be a nonnegative integer, got %s" % (degree))
            d = degree
            if d != degree:
                raise ValueError("degree argument must be a nonnegative integer, got %s" % (degree))
            ZZ_pEX_reverse_hi(r.x, (<Polynomial_ZZ_pEX> self).x, d)
        else:
            ZZ_pEX_reverse(r.x, (<Polynomial_ZZ_pEX> self).x)
        return r

    def inverse_series_trunc(self, prec):
        r"""
        Compute and return the inverse of ``self`` modulo `x^{prec}`.

        The constant term of ``self`` must be invertible.

        EXAMPLES::

            sage: R.<x> = GF(101^2)[]
            sage: z2 =  R.base_ring().gen()
            sage: f = (3*z2 + 57)*x^3 + (13*z2 + 94)*x^2 + (7*z2 + 2)*x + 66*z2 + 15
            sage: f.inverse_series_trunc(1)
            51*z2 + 92
            sage: f.inverse_series_trunc(2)
            (30*z2 + 30)*x + 51*z2 + 92
            sage: f.inverse_series_trunc(3)
            (42*z2 + 94)*x^2 + (30*z2 + 30)*x + 51*z2 + 92
            sage: f.inverse_series_trunc(4)
            (99*z2 + 96)*x^3 + (42*z2 + 94)*x^2 + (30*z2 + 30)*x + 51*z2 + 92

        TESTS::

            sage: R.<x> = GF(163^2)[]
            sage: f = R([p for p in primes(20)])
            sage: f.inverse_series_trunc(1)
            82
            sage: f.inverse_series_trunc(2)
            40*x + 82
            sage: f.inverse_series_trunc(3)
            61*x^2 + 40*x + 82
            sage: f.inverse_series_trunc(0)
            Traceback (most recent call last):
            ...
            ValueError: the precision must be positive, got 0
            sage: f.inverse_series_trunc(-1)
            Traceback (most recent call last):
            ...
            ValueError: the precision must be positive, got -1
            sage: f = x + x^2 + x^3
            sage: f.inverse_series_trunc(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term 0 is not a unit
        """
        self._parent._modulus.restore()

        # Ensure precision is nonnegative
        if prec <= 0:
            raise ValueError("the precision must be positive, got {}".format(prec))

        # Ensure we can invert the constant term
        const_term = self.get_coeff_c(0)
        if not const_term.is_unit():
            raise ValueError("constant term {} is not a unit".format(const_term))

        # Construct output polynomial
        cdef Polynomial_ZZ_pEX r
        r = Polynomial_ZZ_pEX.__new__(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent

        # Call to NTL for the inverse truncation
        if prec > 0:
            sig_on()
            ZZ_pEX_InvTrunc(r.x, self.x, prec)
            sig_off()
        return r

    def __pow__(self, exp, modulus):
        r"""
        Exponentiation of ``self``.

        If ``modulus`` is not ``None``, the exponentiation is performed
        modulo the polynomial ``modulus``.

        EXAMPLES::

            sage: K.<a> = GF(101^2, 'a', modulus=[1,1,1])
            sage: R.<x> = PolynomialRing(K, implementation="NTL")
            sage: pow(x, 100)
            x^100
            sage: pow(x + 3, 5)
            x^5 + 15*x^4 + 90*x^3 + 68*x^2 + x + 41

        If modulus is not ``None``, performs modular exponentiation::

            sage: K.<a> = GF(101^2, 'a', modulus=[1,1,1])
            sage: R.<x> = PolynomialRing(K, implementation="NTL")
            sage: pow(x, 100, x^2 + x + a)
            (19*a + 64)*x + 30*a + 2
            sage: pow(x, 100 * 101**200, x^2 + x + a)
            (19*a + 64)*x + 30*a + 2

        The modulus can have smaller degree than ``self``::

            sage: K.<a> = GF(101^2, 'a', modulus=[1,1,1])
            sage: R.<x> = PolynomialRing(K, implementation="NTL")
            sage: pow(x^4, 25, x^2 + x + a)
            (19*a + 64)*x + 30*a + 2

        TESTS:

        Canonical coercion should apply::

            sage: xx = GF(101)["x"].gen()
            sage: pow(x+1, 25, 2)
            0
            sage: pow(x + a, 101**2, xx^3 + xx + 1)
            4*x^2 + 44*x + a + 70
            sage: pow(x + a, int(101**2), xx^3 + xx + 1)
            4*x^2 + 44*x + a + 70
            sage: xx = polygen(GF(97))
            sage: _ = pow(x + a, 101**2, xx^3 + xx + 1)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: ...
         """
        exp = Integer(exp)
        if modulus is not None:
            # Handle when modulus is zero
            if modulus.is_zero():
                raise ZeroDivisionError("modulus must be nonzero")

            # Similar to coerce_binop
            if not have_same_parent(self, modulus):
                a, m = canonical_coercion(self, modulus)
                if a is not self:
                    return pow(a, exp, m)
                modulus = m
            self = self % modulus
            if exp > 0 and exp.bit_length() >= 32:
                return (<Polynomial_ZZ_pEX>self)._powmod_bigexp(Integer(exp), modulus)
        return Polynomial_template.__pow__(self, exp, modulus)

    cdef _powmod_bigexp(Polynomial_ZZ_pEX self, Integer exp, Polynomial_ZZ_pEX modulus):
        """
        Modular exponentiation for large exponents.
        """
        self._parent._modulus.restore()
        cdef Polynomial_ZZ_pEX r
        cdef ZZ_c e_ZZ
        cdef ZZ_pEX_c y
        cdef ZZ_pEX_Modulus_c mod

        mpz_to_ZZ(&e_ZZ, exp.value)
        r = Polynomial_ZZ_pEX.__new__(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        ZZ_pEX_Modulus_build(mod, modulus.x)

        sig_on()
        if ZZ_pEX_IsX(self.x):
            ZZ_pEX_PowerXMod_ZZ_pre(r.x, e_ZZ, mod)
        elif ZZ_pEX_deg(self.x) < ZZ_pEX_deg(modulus.x):
            ZZ_pEX_PowerMod_ZZ_pre(r.x, self.x, e_ZZ, mod)
        else:
            ZZ_pEX_rem_pre(y, self.x, mod)
            ZZ_pEX_PowerMod_ZZ_pre(r.x, y, e_ZZ, mod)
        sig_off()
        return r

    def compose_mod(self, other, modulus):
        r"""
        Compute `f(g) \bmod h`.

        To be precise about the order fo compostion, given ``self``, ``other``
        and ``modulus`` as `f(x)`, `g(x)` and `h(x)` compute `f(g(x)) \bmod h(x)`.

        INPUT:

        - ``other`` -- a polynomial `g(x)`
        - ``modulus`` -- a polynomial `h(x)`

        EXAMPLES::

            sage: R.<x> = GF(3**6)[]
            sage: f = R.random_element()
            sage: g = R.random_element()
            sage: g.compose_mod(g, f) == g(g) % f
            True

            sage: F.<z3> = GF(3**6)
            sage: R.<x> = F[]
            sage: f = 2*z3^2*x^2 + (z3 + 1)*x + z3^2 + 2
            sage: g = (z3^2 + 2*z3)*x^2 + (2*z3 + 2)*x + 2*z3^2 + z3 + 2
            sage: h = (2*z3 + 2)*x^2 + (2*z3^2 + 1)*x + 2*z3^2 + z3 + 2
            sage: f.compose_mod(g, h)
            (z3^5 + z3^4 + z3^3 + z3^2 + z3)*x + z3^5 + z3^3 + 2*z3 + 2
            sage: f.compose_mod(g, h) == f(g) % h
            True

        AUTHORS:

        - Giacomo Pope (2024-08) initial implementation
        """
        self._parent._modulus.restore()

        # Ensure all the parents match
        if other.parent() is not self._parent:
            other = self._parent.coerce(other)
        if modulus.parent() is not self._parent:
            modulus = self._parent.coerce(modulus)

        # Create the output polynomial
        cdef Polynomial_ZZ_pEX r
        r = Polynomial_ZZ_pEX.__new__(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent

        # Create ZZ_pEX_Modulus type from modulus input
        cdef ZZ_pEX_Modulus_c mod
        ZZ_pEX_Modulus_build(mod, (<Polynomial_ZZ_pEX>modulus).x)

        # Compute f(g) mod h
        sig_on()
        ZZ_pEX_CompMod(r.x, (<Polynomial_ZZ_pEX>self).x, (<Polynomial_ZZ_pEX>(other % modulus)).x, mod)
        sig_off()

        return r

    # compose_mod is the natural name from the NTL bindings, but polynomial_gf2x
    # has modular_composition as the method name so here we allow both
    modular_composition = compose_mod
