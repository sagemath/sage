"""
Univariate Polynomials over domains and fields

AUTHORS:

- William Stein: first version
- Martin Albrecht: Added singular coercion.
- David Harvey: split off polynomial_integer_dense_ntl.pyx (2007-09)
- Robert Bradshaw: split off polynomial_modn_dense_ntl.pyx (2007-09)

TESTS:

We test coercion in a particularly complicated situation::

    sage: W.<w> = QQ['w']
    sage: WZ.<z> = W['z']
    sage: m = matrix(WZ, 2, 2, [1, z, z, z^2])                                          # needs sage.modules
    sage: a = m.charpoly()                                                              # needs sage.modules
    sage: R.<x> = WZ[]
    sage: R(a)                                                                          # needs sage.modules
    x^2 + (-z^2 - 1)*x
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.rings.polynomial.polynomial_element import Polynomial, Polynomial_generic_dense, Polynomial_generic_dense_inexact
from sage.structure.element import IntegralDomainElement, EuclideanDomainElement

from sage.rings.polynomial.polynomial_singular_interface import Polynomial_singular_repr

try:
    from cypari2.gen import Gen as pari_gen
except ImportError:
    pari_gen = ()

from sage.structure.element import coerce_binop, parent
from sage.structure.factorization import Factorization
from sage.structure.richcmp import richcmp, richcmp_item, rich_to_bool, rich_to_bool_sgn

from sage.rings.infinity import infinity, Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.padics.precision_error import PrecisionError


class Polynomial_generic_sparse(Polynomial):
    """
    A generic sparse polynomial.

    The ``Polynomial_generic_sparse`` class defines functionality for sparse
    polynomials over any base ring. A sparse polynomial is represented using a
    dictionary which maps each exponent to the corresponding coefficient. The
    coefficients must never be zero.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(PolynomialRing(QQ, 'y'), sparse=True)
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_integral_domain_with_category.element_class'>
        sage: loads(f.dumps()) == f
        True

    A more extensive example::

        sage: # needs sage.libs.pari
        sage: A.<T> = PolynomialRing(Integers(5), sparse=True)
        sage: f = T^2 + 1; B = A.quo(f)
        sage: C.<s> = PolynomialRing(B)
        sage: C
        Univariate Polynomial Ring in s over Univariate Quotient Polynomial Ring in Tbar
         over Ring of integers modulo 5 with modulus T^2 + 1
        sage: s + T
        s + Tbar
        sage: (s + T)**2
        s^2 + 2*Tbar*s + 4
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        TESTS::

            sage: PolynomialRing(RIF, 'z', sparse=True)([RIF(-1, 1), RIF(-1,1)])                    # needs sage.rings.real_interval_field
            0.?*z + 0.?
            sage: PolynomialRing(RIF, 'z', sparse=True)((RIF(-1, 1), RIF(-1,1)))                    # needs sage.rings.real_interval_field
            0.?*z + 0.?
            sage: PolynomialRing(CIF, 'z', sparse=True)([CIF(RIF(-1,1), RIF(-1,1)), RIF(-1,1)])     # needs sage.rings.complex_interval_field
            0.?*z + 0.? + 0.?*I
        """
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self.__coeffs = {}
            return
        R = parent.base_ring()
        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                x = x.monomial_coefficients()
            elif x.parent() == R:
                x = {0: x}
            else:
                w = {}
                for n, c in x.monomial_coefficients().items():
                    w[n] = R(c)
                # The following line has been added in github issue #9944.
                # Apparently, the "else" case has never occurred before.
                x = w
        elif isinstance(x, (list, tuple)):
            x = {i: c for i, c in enumerate(x) if c}
        elif isinstance(x, pari_gen):
            y = {}
            for i in range(len(x)):
                y[i] = R(x[i])
            x = y
            check = True
        elif not isinstance(x, dict):
            x = {0: x}   # constant polynomials
        if check:
            self.__coeffs = {}
            for i, z in x.items():
                self.__coeffs[i] = R(z)
        else:
            self.__coeffs = x
        if check:
            self.__normalize()

    def monomial_coefficients(self, copy=None):
        """
        Return a new copy of the dict of the underlying
        elements of ``self``.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: d = f.monomial_coefficients(); d
            {0: 5, 1997: 1, 10000: 7}
            sage: d[0] = 10
            sage: f.monomial_coefficients()
            {0: 5, 1997: 1, 10000: 7}

        ``dict`` is an alias::

            sage: f.dict()
            {0: 5, 1997: 1, 10000: 7}
        """
        return dict(self.__coeffs)

    dict = monomial_coefficients

    def coefficients(self, sparse=True):
        """
        Return the coefficients of the monomials appearing in ``self``.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: f.coefficients()
            [5, 1, 7]

        TESTS:

        Check that all coefficients are in the base ring::

            sage: S.<x> = PolynomialRing(QQ, sparse=True)
            sage: f = x^4
            sage: all(c.parent() is QQ for c in f.coefficients(False))
            True
        """
        if sparse:
            return [self.__coeffs[e] for e in self.exponents()]
        else:
            zero = self.parent().base_ring().zero()
            return [self.__coeffs[i] if i in self.__coeffs else zero
                    for i in range(self.degree() + 1)]

    def exponents(self):
        """
        Return the exponents of the monomials appearing in ``self``.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: f.exponents()
            [0, 1997, 10000]
        """
        return sorted(self.__coeffs)

    def valuation(self, p=None):
        """
        Return the valuation of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<w> = PolynomialRing(GF(9, 'a'), sparse=True)
            sage: f = w^1997 - w^10000
            sage: f.valuation()
            1997
            sage: R(19).valuation()
            0
            sage: R(0).valuation()
            +Infinity
        """
        if not self.__coeffs:
            return infinity

        if p is infinity:
            return -self.degree()

        if p is not None:
            raise NotImplementedError("input p is not support for sparse polynomials")

        return ZZ(min(self.__coeffs))

    def _derivative(self, var=None):
        """
        Return the formal derivative of this polynomial with respect to
        the given variable.

        If ``var`` is ``None`` or is the generator of this ring, the
        derivative is with respect to the generator. Otherwise,
        _derivative(var) is called recursively for each coefficient of
        this polynomial.

        .. SEEALSO:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<w> = PolynomialRing(ZZ, sparse=True)
            sage: f = R(range(9)); f
            8*w^8 + 7*w^7 + 6*w^6 + 5*w^5 + 4*w^4 + 3*w^3 + 2*w^2 + w
            sage: f._derivative()
            64*w^7 + 49*w^6 + 36*w^5 + 25*w^4 + 16*w^3 + 9*w^2 + 4*w + 1
            sage: f._derivative(w)
            64*w^7 + 49*w^6 + 36*w^5 + 25*w^4 + 16*w^3 + 9*w^2 + 4*w + 1

            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<y> = PolynomialRing(R, sparse=True)
            sage: f = x^3*y^4
            sage: f._derivative()
            4*x^3*y^3
            sage: f._derivative(y)
            4*x^3*y^3
            sage: f._derivative(x)
            3*x^2*y^4

        Check that :issue:`28187` is fixed::

            sage: R = PolynomialRing(ZZ, 't', sparse=True)
            sage: t, u = var('t, u')                                                    # needs sage.symbolic
            sage: R.gen()._derivative(t)                                                # needs sage.symbolic
            1
            sage: R.gen()._derivative(u)                                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to u
        """
        P = self.parent()
        if var is not None and var != P.gen():
            try:
                # call _derivative() recursively on coefficients
                return P({n:self.__coeffs[n]._derivative(var)
                            for n in self.__coeffs})
            except AttributeError:
                raise ValueError('cannot differentiate with respect to {}'.format(var))

        # compute formal derivative with respect to generator
        d = {}
        for n, c in self.__coeffs.items():
            d[n-1] = n*c
        if -1 in d:
            del d[-1]
        return P(d)

    def integral(self, var=None):
        """
        Return the integral of this polynomial.

        By default, the integration variable is the variable of the
        polynomial.

        Otherwise, the integration variable is the optional parameter ``var``

        .. NOTE::

            The integral is always chosen so that the constant term is 0.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (1 + 3*x^10 - 2*x^100).integral()
            -2/101*x^101 + 3/11*x^11 + x

        TESTS:

        Check that :issue:`18600` is fixed::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^2^100).integral()
            1/1267650600228229401496703205377*x^1267650600228229401496703205377

        Check the correctness when the base ring is a polynomial ring::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: S.<t> = PolynomialRing(R, sparse=True)
            sage: (x*t+1).integral()
            1/2*x*t^2 + t
            sage: (x*t+1).integral(x)
            1/2*x^2*t + x

        Check the correctness when the base ring is not an integral domain::

            sage: R.<x> = PolynomialRing(Zmod(4), sparse=True)
            sage: (x^4 + 2*x^2  + 3).integral()
            x^5 + 2*x^3 + 3*x
            sage: x.integral()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(2, 4) does not exist
        """
        R = self.parent()
        # TODO:
        # calling the coercion model bin_op is much more accurate than using the
        # true division (which is bypassed by polynomials). But it does not work
        # in all cases!!
        from sage.structure.element import coercion_model as cm
        import operator
        try:
            Q = cm.bin_op(R.one(), ZZ.one(), operator.truediv).parent()
        except TypeError:
            F = (R.base_ring().one()/ZZ.one()).parent()
            Q = R.change_ring(F)

        if var is not None and var != R.gen():
            return Q({k:v.integral(var) for k,v in self.__coeffs.items()}, check=False)

        return Q({ k+1:v/(k+1) for k,v in self.__coeffs.items()}, check=False)

    def _dict_unsafe(self):
        """
        Return unsafe access to the underlying dictionary of coefficients.

        ** DO NOT use this, unless you really really know what you are doing. **

        EXAMPLES::

            sage: R.<w> = PolynomialRing(ZZ, sparse=True)
            sage: f = w^15 - w*3; f
            w^15 - 3*w
            sage: d = f._dict_unsafe(); d
            {1: -3, 15: 1}
            sage: d[1] = 10; f
            w^15 + 10*w
        """
        return self.__coeffs

    def _repr(self, name=None):
        r"""
        EXAMPLES::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: R.<w> = PolynomialRing(CDF, sparse=True)
            sage: f = CDF(1,2) + w^5 - CDF(pi)*w + CDF(e)
            sage: f._repr()   # abs tol 1e-15
            '1.0*w^5 - 3.141592653589793*w + 3.718281828459045 + 2.0*I'
            sage: f._repr(name='z')   # abs tol 1e-15
            '1.0*z^5 - 3.141592653589793*z + 3.718281828459045 + 2.0*I'

        TESTS::

            sage: pol = RIF['x']([0, 0, (-1,1)])                                        # needs sage.rings.real_interval_field
            sage: PolynomialRing(RIF, 'x', sparse=True)(pol)                            # needs sage.rings.real_interval_field
            0.?*x^2

        AUTHOR:

        - David Harvey (2006-08-05), based on Polynomial._repr()
        - Francis Clarke (2008-09-08) improved for 'negative' coefficients
        """
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = sorted(self.__coeffs.items())
        for (n, x) in reversed(coeffs):
            if x:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)" % x
                if n > 1:
                    var = "*%s^%s" % (name,n)
                elif n == 1:
                    var = "*%s" % name
                else:
                    var = ""
                s += "%s%s" % (x,var)
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if s == " ":
            return "0"
        return s[1:]

    def __normalize(self):
        x = self.__coeffs
        D = [n for n, z in x.items() if not z]
        for n in D:
            del x[n]

    def __getitem__(self, n):
        """
        Return the `n`-th coefficient of this polynomial.

        Negative indexes are allowed and always return 0 (so you can
        view the polynomial as embedding Laurent series).

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: R.<w> = PolynomialRing(RDF, sparse=True)
            sage: e = RDF(e)
            sage: f = sum(e^n*w^n for n in range(4)); f   # abs tol 1.1e-14
            20.085536923187664*w^3 + 7.3890560989306495*w^2 + 2.718281828459045*w + 1.0
            sage: f[1]  # abs tol 5e-16
            2.718281828459045
            sage: f[5]
            0.0
            sage: f[-1]
            0.0

            sage: R.<x> = PolynomialRing(RealField(19), sparse=True)                    # needs sage.rings.real_mpfr
            sage: f = (2-3.5*x)^3; f                                                    # needs sage.rings.real_mpfr
            -42.875*x^3 + 73.500*x^2 - 42.000*x + 8.0000

        Using slices, we can truncate polynomials::

            sage: f[:2]                                                                 # needs sage.rings.real_mpfr
            -42.000*x + 8.0000

        Any other kind of slicing is an error, see :issue:`18940`::

            sage: f[1:3]                                                                # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a start is not defined

            sage: f[1:3:2]                                                              # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a step is not defined

        TESTS::

            sage: f["hello"]                                                            # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: list indices must be integers, not str
        """
        if isinstance(n, slice):
            start, stop, step = n.start, n.stop, n.step
            if step is not None:
                raise IndexError("polynomial slicing with a step is not defined")
            if start is not None:
                raise IndexError("polynomial slicing with a start is not defined")
            d = self.degree() + 1
            if stop is None or stop > d:
                stop = d
            v = {key: val for key, val in self.__coeffs.items()
                 if key < stop}
            return self.parent()(v)

        try:
            n = n.__index__()
        except AttributeError:
            raise TypeError("list indices must be integers, not {0}".format(type(n).__name__))
        try:
            return self.__coeffs[n]
        except KeyError:
            return self.base_ring().zero()

    def _unsafe_mutate(self, n, value):
        r"""
        Change the coefficient of `x^n` to value.

        ** NEVER USE THIS ** -- unless you really know what you are doing.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: R.<z> = PolynomialRing(CC, sparse=True)
            sage: f = z^2 + CC.0; f
            1.00000000000000*z^2 + 1.00000000000000*I
            sage: f._unsafe_mutate(0, 10)
            sage: f
            1.00000000000000*z^2 + 10.0000000000000

        Much more nasty::

            sage: z._unsafe_mutate(1, 0)                                                # needs sage.rings.real_mpfr
            sage: z                                                                     # needs sage.rings.real_mpfr
            0
        """
        n = int(n)
        value = self.base_ring()(value)
        x = self.__coeffs
        if n < 0:
            raise IndexError("polynomial coefficient index must be nonnegative")
        if value == 0:
            if n in x:
                del x[n]
        else:
            x[n] = value

    def list(self, copy=True):
        """
        Return a new copy of the list of the underlying
        elements of ``self``.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(Integers(100), sparse=True)
            sage: f = 13*z^5 + 15*z^2 + 17*z
            sage: f.list()
            [0, 17, 15, 0, 0, 13]
        """
        zero = self.base_ring().zero()
        v = [zero] * (self.degree()+1)
        for n, x in self.__coeffs.items():
            v[n] = x
        return v

    def degree(self, gen=None):
        """
        Return the degree of this sparse polynomial.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(ZZ, sparse=True)
            sage: f = 13*z^50000 + 15*z^2 + 17*z
            sage: f.degree()
            50000
        """
        if not self.__coeffs:
            return -1
        return max(self.__coeffs)

    def __floordiv__(self, right):
        """
        Return the quotient upon division (no remainder).

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQbar, sparse=True)
            sage: f = (1+2*x)^3 + 3*x; f
            8*x^3 + 12*x^2 + 9*x + 1
            sage: g = f // (1+2*x); g
            4*x^2 + 4*x + 5/2
            sage: f - g * (1+2*x)
            -3/2
            sage: f.quo_rem(1+2*x)
            (4*x^2 + 4*x + 5/2, -3/2)
        """
        P = self.parent()
        if P is parent(right):
            return self._floordiv_(right)
        d = P.base_ring()(right)
        return self.map_coefficients(lambda c: c // d)

    def _add_(self, right):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(), sparse=True)
            sage: (x^100000 + 2*x^50000) + (4*x^75000 - 2*x^50000 + 3*x)
            x^100000 + 4*x^75000 + 3*x

        AUTHOR:

        - David Harvey (2006-08-05)
        """
        output = dict(self.__coeffs)

        for (index, coeff) in right.__coeffs.items():
            if index in output:
                output[index] += coeff
            else:
                output[index] = coeff

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _neg_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(), sparse=True)
            sage: a = x^10000000; a
            x^10000000
            sage: -a
            -x^10000000
        """
        output = { }
        for (index, coeff) in self.__coeffs.items():
            output[index] = -coeff
        output = self.parent()(output, check=False)
        return output

    def _mul_(self, right):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^100000 - x^50000) * (x^100000 + x^50000)
             x^200000 - x^100000
            sage: (x^100000 - x^50000) * R(0)
             0

        AUTHOR:
        - David Harvey (2006-08-05)
        """
        output = {}

        for (index1, coeff1) in self.__coeffs.items():
            for (index2, coeff2) in right.__coeffs.items():
                product = coeff1 * coeff2
                index = index1 + index2
                if index in output:
                    output[index] += product
                else:
                    output[index] = product

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _rmul_(self, left):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^100000 - x^50000) * (x^100000 + x^50000)
            x^200000 - x^100000
            sage: 7 * (x^100000 - x^50000)   # indirect doctest
            7*x^100000 - 7*x^50000

        AUTHOR:

        - Simon King (2011-03-31)
        """
        output = {}

        for (index, coeff) in self.__coeffs.items():
            output[index] = left * coeff

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _lmul_(self, right):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^100000 - x^50000) * (x^100000 + x^50000)
            x^200000 - x^100000
            sage: (x^100000 - x^50000) * 7   # indirect doctest
            7*x^100000 - 7*x^50000

        AUTHOR:

        - Simon King (2011-03-31)
        """
        output = {}

        for (index, coeff) in self.__coeffs.items():
            output[index] = coeff * right

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _richcmp_(self, other, op):
        """
        Compare this polynomial with other.

        Polynomials are first compared by degree, then in dictionary order
        starting with the coefficient of largest degree.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: 3*x^100 - 12 > 12*x + 5
            True
            sage: 3*x^100 - 12 > 3*x^100 - x^50 + 5
            True
            sage: 3*x^100 - 12 < 3*x^100 - x^50 + 5
            False
            sage: x^100 + x^10 - 1 < x^100 + x^10
            True
            sage: x^100 < x^100 - x^10
            False

        TESTS::

            sage: R.<x> = PolynomialRing(QQ, sparse=True)
            sage: 2*x^2^500 > x^2^500
            True

            sage: Rd = PolynomialRing(ZZ, 'x', sparse=False)
            sage: Rs = PolynomialRing(ZZ, 'x', sparse=True)
            sage: for _ in range(100):
            ....:     pd = Rd.random_element()
            ....:     qd = Rd.random_element()
            ....:     assert bool(pd < qd) == bool(Rs(pd) < Rs(qd))
        """
        d1 = self.degree()
        d2 = other.degree()

        # Special case constant polynomials
        if d1 <= 0 and d2 <= 0:
            return richcmp(self[0], other[0], op)

        # For different degrees, compare the degree
        if d1 != d2:
            return rich_to_bool_sgn(op, d1 - d2)

        degs = set(self.__coeffs) | set(other.__coeffs)
        for i in sorted(degs, reverse=True):
            x = self[i]
            y = other[i]
            res = richcmp_item(x, y, op)
            if res is not NotImplemented:
                return res
        return rich_to_bool(op, 0)

    def shift(self, n):
        r"""
        Return this polynomial multiplied by the power `x^n`.

        If `n` is negative, terms below `x^n` will be discarded. Does
        not change this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^100000 + 2*x + 4
            sage: type(p)
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_integral_domain_with_category.element_class'>
            sage: p.shift(0)
             x^100000 + 2*x + 4
            sage: p.shift(-1)
             x^99999 + 2
            sage: p.shift(-100002)
             0
            sage: p.shift(2)
             x^100002 + 2*x^3 + 4*x^2

        TESTS:

        Check that :issue:`18600` is fixed::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^2^100 - 5
            sage: p.shift(10)
            x^1267650600228229401496703205386 - 5*x^10
            sage: p.shift(-10)
            x^1267650600228229401496703205366
            sage: p.shift(1.5)                                                          # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer

        AUTHOR:
        - David Harvey (2006-08-06)
        """
        n = ZZ(n)
        if n == 0:
            return self
        if n > 0:
            output = {index+n: coeff for index, coeff in self.__coeffs.items()}
            return self.parent()(output, check=False)
        if n < 0:
            output = {index+n:coeff for index, coeff in self.__coeffs.items() if index + n >= 0}
            return self.parent()(output, check=False)

    @coerce_binop
    def quo_rem(self, other):
        """
        Return the quotient and remainder of the Euclidean division of
        ``self`` and ``other``.

        Raises :exc:`ZeroDivisionError` if ``other`` is zero.

        Raises :exc:`ArithmeticError` if ``other`` has a nonunit leading
        coefficient and this causes the Euclidean division to fail.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(ZZ, sparse=True)
            sage: R.<y> = PolynomialRing(P, sparse=True)
            sage: f = R.random_element(10)
            sage: while x.divides(f.leading_coefficient()):
            ....:     f = R.random_element(10)
            sage: g = y^5 + R.random_element(4)
            sage: q, r = f.quo_rem(g)
            sage: f == q*g + r and r.degree() < g.degree()
            True
            sage: g = x*y^5
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ArithmeticError: Division non exact
            (consider coercing to polynomials over the fraction field)
            sage: g = 0
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by zero polynomial

        If the leading coefficient of ``other`` is not a unit, Euclidean division may still work::

            sage: f = -x*y^10 + 2*x*y^7 + y^3 - 2*x^2*y^2 - y
            sage: g = x*y^5
            sage: f.quo_rem(g)
            (-y^5 + 2*y^2, y^3 - 2*x^2*y^2 - y)

        Polynomials over noncommutative rings are also allowed::

            sage: # needs sage.combinat sage.modules
            sage: HH = QuaternionAlgebra(QQ, -1, -1)
            sage: P.<x> = PolynomialRing(HH, sparse=True)
            sage: f = P.random_element(5)
            sage: g = P.random_element((0, 5))
            sage: q, r = f.quo_rem(g)
            sage: f == q*g + r
            True

        TESTS::

            sage: P.<x> = PolynomialRing(ZZ, sparse=True)
            sage: f = x^10 - 4*x^6 - 5
            sage: g = 17*x^22 + x^15 - 3*x^5 + 1
            sage: q, r = g.quo_rem(f)
            sage: g == f*q + r and r.degree() < f.degree()
            True
            sage: zero = P(0)
            sage: zero.quo_rem(f)
            (0, 0)
            sage: Q.<y> = IntegerModRing(14)[]
            sage: f = y^10 - 4*y^6 - 5
            sage: g = 17*y^22 + y^15 - 3*y^5 + 1
            sage: q, r = g.quo_rem(f)
            sage: g == f*q + r and r.degree() < f.degree()
            True
            sage: f += 2*y^10  # 3 is invertible mod 14
            sage: q, r = g.quo_rem(f)
            sage: g == f*q + r and r.degree() < f.degree()
            True

        The following shows that :issue:`16649` is indeed fixed. ::

            sage: P.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (4*x).quo_rem(2*x)
            (2, 0)

        AUTHORS:

        - Bruno Grenet (2014-07-09)
        """
        if other.is_zero():
            raise ZeroDivisionError("Division by zero polynomial")
        if self.is_zero():
            return self, self

        R = self.parent()

        d = other.degree()
        if self.degree() < d:
            return R.zero(), self

        quo = R.zero()
        rem = self

        while rem.degree() >= d:
            try:
                c = R(rem.leading_coefficient() * ~other.leading_coefficient())
            except TypeError:
                raise ArithmeticError("Division non exact (consider coercing to polynomials over the fraction field)")
            e = rem.degree() - d
            quo += c*R.one().shift(e)
            # we know that the leading coefficient of rem vanishes
            # thus we avoid doing a useless computation
            rem = rem[:rem.degree()] - c*other[:d].shift(e)
        return (quo,rem)

    @coerce_binop
    def gcd(self, other, algorithm=None):
        r"""
        Return the gcd of this polynomial and ``other``.

        INPUT:

        - ``other`` -- a polynomial defined over the same ring as this
          polynomial

        ALGORITHM:

        Two algorithms are provided:

        - ``'generic'`` -- uses the generic implementation, which depends on the
          base ring being a UFD or a field
        - ``'dense'`` -- the polynomials are converted to the dense
          representation, their gcd is computed and is converted back to the
          sparse representation

        Default is ``'dense'`` for polynomials over `\ZZ` and ``'generic'`` in the
        other cases.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^6 + 7*x^5 + 8*x^4 + 6*x^3 + 2*x^2 + x + 2
            sage: q = 2*x^4 - x^3 - 2*x^2 - 4*x - 1
            sage: gcd(p, q)
            x^2 + x + 1
            sage: gcd(p, q, algorithm='dense')
            x^2 + x + 1
            sage: gcd(p, q, algorithm='generic')
            x^2 + x + 1
            sage: gcd(p, q, algorithm='foobar')
            Traceback (most recent call last):
            ...
            ValueError: Unknown algorithm 'foobar'

        TESTS:

        Check that :issue:`19676` is fixed::

            sage: S.<y> = R[]
            sage: x.gcd(y)
            1
            sage: (6*x).gcd(9)
            3

        Check that :issue:`36427` is fixed::

            sage: P = PolynomialRing(ZZ, "q", sparse=True)
            sage: q = P.gen()
            sage: 2*q^-100
            2/q^100
            sage: gcd(1, q^100)
            1
            sage: gcd(q^0, q^100)
            1
        """

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        if algorithm is None:
            if self.base_ring() == ZZ:
                algorithm = "dense"
            else:
                algorithm = "generic"
        if algorithm == "dense":
            S = self.parent()
            # FLINT is faster but a bug makes the conversion extremely slow,
            # so NTL is used in those cases where the conversion is too slow. Cf
            # <https://groups.google.com/d/msg/sage-devel/6qhW90dgd1k/Hoq3N7fWe4QJ>
            sd = self.degree()
            od = other.degree()
            if ((sd < 100 or len(self.__coeffs)/sd > .06)
                    and (od < 100 or len(other.__coeffs)/od > .06)):
                implementation = "FLINT"
            else:
                implementation = "NTL"
            D = PolynomialRing(S.base_ring(),'x',implementation=implementation)
            g = D(self).gcd(D(other))
            return S(g)
        elif algorithm == "generic":
            return Polynomial.gcd(self,other)
        else:
            raise ValueError("Unknown algorithm '%s'" % algorithm)

    def reverse(self, degree=None):
        """
        Return this polynomial but with the coefficients reversed.

        If an optional degree argument is given, the coefficient list will be
        truncated or zero padded as necessary and the reverse polynomial will
        have the specified degree.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^4 + 2*x^2^100
            sage: p.reverse()
            x^1267650600228229401496703205372 + 2
            sage: p.reverse(10)
            x^6
        """
        if degree is None:
            degree = self.degree()
        if not isinstance(degree, (int,Integer)):
            raise ValueError("degree argument must be a nonnegative integer, got %s" % degree)
        d = {degree-k: v for k,v in self.__coeffs.items() if degree >= k}
        return self.parent()(d, check=False)

    def truncate(self, n):
        """
        Return the polynomial of degree `< n` equal to ``self`` modulo `x^n`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^11 + x^10 + 1).truncate(11)
            x^10 + 1
            sage: (x^2^500 + x^2^100 + 1).truncate(2^101)
            x^1267650600228229401496703205376 + 1
        """
        return self[:n]

    def number_of_terms(self):
        """
        Return the number of nonzero terms.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^100 - 3*x^10 + 12
            sage: p.number_of_terms()
            3
        """
        return len(self.__coeffs)


class Polynomial_generic_domain(Polynomial, IntegralDomainElement):
    def __init__(self, parent, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

    def is_unit(self):
        r"""
        Return ``True`` if this polynomial is a unit.

        *EXERCISE* (Atiyah-McDonald, Ch 1): Let `A[x]` be a polynomial
        ring in one variable.  Then `f=\sum a_i x^i \in A[x]` is a
        unit if and only if `a_0` is a unit and `a_1,\ldots, a_n` are
        nilpotent.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(ZZ, sparse=True)
            sage: (2 + z^3).is_unit()
            False
            sage: f = -1 + 3*z^3; f
            3*z^3 - 1
            sage: f.is_unit()
            False
            sage: R(-3).is_unit()
            False
            sage: R(-1).is_unit()
            True
            sage: R(0).is_unit()
            False
        """
        if self.degree() > 0:
            return False
        return self[0].is_unit()


class Polynomial_generic_field(Polynomial_singular_repr,
                               Polynomial_generic_domain,
                               EuclideanDomainElement):

    @coerce_binop
    def quo_rem(self, other):
        """
        Return a tuple ``(quotient, remainder)`` where
        ``self = quotient * other + remainder``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<y> = PolynomialRing(QQ)
            sage: K.<t> = NumberField(y^2 - 2)
            sage: P.<x> = PolynomialRing(K)
            sage: x.quo_rem(K(1))
            (x, 0)
            sage: x.xgcd(K(1))
            (1, 0, 1)
        """
        P = self.parent()
        if other.is_zero():
            raise ZeroDivisionError("other must be nonzero")

        # This is algorithm 3.1.1 in Cohen GTM 138
        A = self
        B = other
        R = A
        Q = P.zero()
        while R.degree() >= B.degree():
            aaa = R.leading_coefficient()/B.leading_coefficient()
            diff_deg = R.degree()-B.degree()
            Q += P(aaa).shift(diff_deg)
            # We know that S*B exactly cancels the leading coefficient of R.
            # Thus, we skip the computation of this leading coefficient.
            # For most exact fields, this doesn't matter much; but for
            # inexact fields, the leading coefficient might not end up
            # exactly equal to zero; and for AA/QQbar, verifying that
            # the coefficient is exactly zero triggers exact computation.
            R = R[:R.degree()] - (aaa*B[:B.degree()]).shift(diff_deg)
        return (Q, R)


class Polynomial_generic_sparse_field(Polynomial_generic_sparse, Polynomial_generic_field):
    """
    EXAMPLES::

        sage: R.<x> = PolynomialRing(Frac(RR['t']), sparse=True)
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_field_with_category.element_class'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial_generic_sparse.__init__(self, parent, x, check, is_gen)


class Polynomial_generic_dense_field(Polynomial_generic_dense, Polynomial_generic_field):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial_generic_dense.__init__(self, parent, x, check, is_gen)


##########################################
# Over discrete valuation rings and fields
##########################################

class Polynomial_generic_cdv(Polynomial_generic_domain):
    """
    A generic class for polynomials over complete discrete
    valuation domains and fields.

    AUTHOR:

    - Xavier Caruso (2013-03)
    """
    def newton_slopes(self, repetition=True):
        """
        Return a list of the Newton slopes of this polynomial.

        These are the valuations of the roots of this polynomial.

        If ``repetition`` is ``True``, each slope is repeated a number of
        times equal to its multiplicity. Otherwise it appears only
        one time.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron sage.rings.padics
            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10
            sage: f.newton_polygon()
            Finite Newton polygon with 4 vertices: (0, 1), (1, 0), (4, 0), (10, 2)
            sage: f.newton_slopes()
            [1, 0, 0, 0, -1/3, -1/3, -1/3, -1/3, -1/3, -1/3]
            sage: f.newton_slopes(repetition=False)
            [1, 0, -1/3]

        AUTHOR:

        - Xavier Caruso (2013-03-20)
        """
        polygon = self.newton_polygon()
        return [-s for s in polygon.slopes(repetition=repetition)]

    def newton_polygon(self):
        r"""
        Return a list of vertices of the Newton polygon of this polynomial.

        .. NOTE::

            If some coefficients have not enough precision an error is raised.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron sage.rings.padics
            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10
            sage: f.newton_polygon()
            Finite Newton polygon with 4 vertices: (0, 1), (1, 0), (4, 0), (10, 2)
            sage: g = f + K(0,0)*t^4; g
            (5^2 + O(5^22))*t^10 + O(5^0)*t^4 + (3 + O(5^20))*t + 5 + O(5^21)
            sage: g.newton_polygon()
            Traceback (most recent call last):
            ...
            PrecisionError: The coefficient of t^4 has not enough precision

        TESTS:

        Check that :issue:`22936` is fixed::

            sage: S.<x> = PowerSeriesRing(GF(5))
            sage: R.<y> = S[]
            sage: p = x^2 + y + x*y^2
            sage: p.newton_polygon()                                                    # needs sage.geometry.polyhedron
            Finite Newton polygon with 3 vertices: (0, 2), (1, 0), (2, 1)

        AUTHOR:

        - Xavier Caruso (2013-03-20)
        """
        d = self.degree()
        from sage.geometry.newton_polygon import NewtonPolygon
        polygon = NewtonPolygon([(x, self[x].valuation()) for x in range(d+1)])
        polygon_prec = NewtonPolygon([ (x, self[x].precision_absolute()) for x in range(d+1) ])
        vertices = polygon.vertices(copy=False)
        vertices_prec = polygon_prec.vertices(copy=False)
        if len(vertices_prec) > 0:
            if vertices[0][0] > vertices_prec[0][0]:
                raise PrecisionError("first term with non-infinite valuation must have determined valuation")
            elif vertices[-1][0] < vertices_prec[-1][0]:
                raise PrecisionError("last term with non-infinite valuation must have determined valuation")
            else:
                for (x, y) in vertices:
                    if polygon_prec(x) <= y:
                        raise PrecisionError("The coefficient of %s^%s has not enough precision" % (self.parent().variable_name(), x))
        return polygon

    def hensel_lift(self, a):
        """
        Lift `a` to a root of this polynomial (using
        Newton iteration).

        If `a` is not close enough to a root (so that
        Newton iteration does not converge), an error
        is raised.

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: K = Qp(5, 10)
            sage: P.<x> = PolynomialRing(K)
            sage: f = x^2 + 1
            sage: root = f.hensel_lift(2); root
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: f(root)
            O(5^10)

            sage: g = (x^2 + 1) * (x - 7)                                               # needs sage.rings.padics
            sage: g.hensel_lift(2)  # here, 2 is a multiple root modulo p               # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: a is not close enough to a root of this polynomial

        AUTHOR:

        - Xavier Caruso (2013-03-23)
        """
        selfa = self(a)
        der = self.derivative()
        dera = der(a)
        if selfa.valuation() <= 2 * dera.valuation():
            raise ValueError("a is not close enough to a root of this polynomial")
        # Newton iteration
        # Todo: compute everything up to the adequate precision at each step
        b = ~dera
        while True:
            na = a - selfa * b
            if na == a:
                return a
            a = na
            selfa = self(a)
            dera = der(a)
            b *= 2 - dera * b

    def _factor_of_degree(self, deg):
        """
        Return a factor of ``self`` of degree ``deg``.

        Algorithm is Newton iteration.

        This fails if ``deg`` is not a breakpoint in the Newton
        polygon of ``self``.

        Only for internal use!

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: K = Qp(5)
            sage: R.<x> = K[]
            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10

            sage: g = f._factor_of_degree(4)                                            # needs sage.rings.padics
            sage: (f % g).is_zero()                                                     # needs sage.rings.padics
            True

            sage: g = f._factor_of_degree(3)    # not tested                            # needs sage.rings.padics
            Traceback (most recent call last)
            ...
            KeyboardInterrupt:

        TESTS::

            sage: S.<x> = PowerSeriesRing(GF(5))
            sage: R.<y> = S[]
            sage: p = x^2 + y + x*y^2
            sage: p._factor_of_degree(1)
            (1 + O(x^20))*y + x^2 + x^5 + 2*x^8 + 4*x^14 + 2*x^17 + 2*x^20 + O(x^22)

        AUTHOR:

        - Xavier Caruso (2013-03-20)

        .. TODO::

            Precision is not optimal, and can be improved.
        """
        coeffs = self.list()
        a = coeffs[:deg+1]
        # The leading coefficient need to be known at finite precision
        # in order to ensure that the while loop below terminates
        if a[deg].precision_absolute() is Infinity:
            a[deg] = a[deg].add_bigoh(self.base_ring().default_prec())

        parent = self.parent()
        a = parent(a)
        v = parent.one()
        x = self % a
        while not x.is_zero():
            a += (v * x) % a
            b, x = self.quo_rem(a)
            b %= a
            v = (v * (2 - b * v)) % a
        return a

    def factor_of_slope(self, slope=None):
        """
        INPUT:

        - ``slope`` -- a rational number (default: the first slope
          in the Newton polygon of ``self``)

        OUTPUT:

        The factor of ``self`` corresponding to the slope ``slope`` (i.e.
        the unique monic divisor of ``self`` whose slope is ``slope`` and
        degree is the length of ``slope`` in the Newton polygon).

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron sage.rings.padics
            sage: K = Qp(5)
            sage: R.<x> = K[]
            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10
            sage: f.newton_slopes()
            [1, 0, 0, 0, -1/3, -1/3, -1/3, -1/3, -1/3, -1/3]
            sage: g = f.factor_of_slope(0)
            sage: g.newton_slopes()
            [0, 0, 0]
            sage: (f % g).is_zero()
            True
            sage: h = f.factor_of_slope()
            sage: h.newton_slopes()
            [1]
            sage: (f % h).is_zero()
            True

        If ``slope`` is not a slope of ``self``, the corresponding factor
        is `1`::

            sage: f.factor_of_slope(-1)                                                 # needs sage.geometry.polyhedron sage.rings.padics
            1 + O(5^20)

        AUTHOR:

        - Xavier Caruso (2013-03-20)
        """
        one = self.parent()(1)
        vertices = self.newton_polygon().vertices(copy=False)
        if len(vertices) < 2:
            if slope is Infinity:
                return self.parent().gen() ** self.degree()
            else:
                return one
        if slope is None:
            deg_first = vertices[0][0]
            deg_last = vertices[1][0]
        else:
            (deg_first, y_first) = vertices[0]
            for i in range(1, len(vertices)):
                (deg_last, y_last) = vertices[i]
                slope_cur = (y_first - y_last) / (deg_last - deg_first)
                if slope_cur == slope:
                    break
                elif slope_cur < slope:
                    return one
                deg_first = deg_last
                y_first = y_last
            if slope_cur > slope:
                return one
        if deg_last == self.degree():
            div = self
        else:
            div = self._factor_of_degree(deg_last)
        if deg_first > 0:
            div2 = div._factor_of_degree(deg_first)
            div,_ = div.quo_rem(div2)
        return div.monic()

    def slope_factorization(self):
        """
        Return a factorization of ``self`` into a product of factors
        corresponding to each slope in the Newton polygon.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron sage.rings.padics
            sage: K = Qp(5)
            sage: R.<x> = K[]
            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10
            sage: f.newton_slopes()
            [1, 0, 0, 0, -1/3, -1/3, -1/3, -1/3, -1/3, -1/3]
            sage: F = f.slope_factorization()
            sage: F.prod() == f
            True
            sage: for (f,_) in F:
            ....:     print(f.newton_slopes())
            [-1/3, -1/3, -1/3, -1/3, -1/3, -1/3]
            [0, 0, 0]
            [1]

        TESTS::

            sage: S.<x> = PowerSeriesRing(GF(5))
            sage: R.<y> = S[]
            sage: p = x^2 + y + x*y^2
            sage: p.slope_factorization()                                               # needs sage.geometry.polyhedron
            (x)
            * ((x + O(x^22))*y + 1 + 4*x^3 + 4*x^6 + 3*x^9 + x^15 + 3*x^18 + O(x^21))
            * ((x^-1 + O(x^20))*y + x + x^4 + 2*x^7 + 4*x^13 + 2*x^16 + 2*x^19 + O(x^22))

        AUTHOR:

        - Xavier Caruso (2013-03-20)
        """
        vertices = self.newton_polygon().vertices(copy=False)

        unit = self.leading_coefficient()
        P = ~unit * self

        deg_first = vertices[0][0]
        factors = [ ]
        if deg_first > 0:
            P >>= deg_first
            factors.append((self._parent.gen(), deg_first))
        if len(vertices) > 2:
            for i in range(1, len(vertices)-1):
                deg = vertices[i][0]
                div = P._factor_of_degree(deg-deg_first)
                factors.append((div,1))
                P,_ = P.quo_rem(div)
                deg_first = deg
        if len(vertices) > 1:
            factors.append((P, 1))
        factors.reverse()
        return Factorization(factors, sort=False, unit=unit)

    def _roots(self, secure, minval, hint):
        """
        Return the roots of this polynomial whose valuation is
        at least ``minval``.

        This is a helper method for :meth:`roots`.
        It is not meant to be called directly.

        INPUT:

        - ``secure`` -- boolean; whether we raise an error or
          not in case of multiple roots

        - ``minval`` -- integer

        - ``hint`` -- list or ``None``; if given, it must be the
          list of roots of the residual polynomial of slope ``minval``

        OUTPUT: list of pairs ``(root, multiplicity)``

        TESTS::

            sage: # needs sage.geometry.polyhedron sage.rings.padics
            sage: R = Zp(2)
            sage: S.<x> = R[]
            sage: P = (x-1) * (x-2) * (x-4) * (x-8) * (x-16)
            sage: Q = P^2
            sage: Q.roots(algorithm='sage')  # indirect doctest
            [(2^4 + O(2^14), 2),
             (2^3 + O(2^13), 2),
             (2^2 + O(2^12), 2),
             (2 + O(2^11), 2),
             (1 + O(2^10), 2)]
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        K = self.base_ring()
        Pk = PolynomialRing(K.residue_field(), names='xbar')
        x = self.parent().gen()

        # Trivial cases
        if self.degree() == 0:
            return [ ]
        if self.degree() == 1:
            return [ (-self[0]/self[1], 1) ]

        # We consider the case where zero is a (possibly multiple) root
        i = 0
        while self[i] == 0:
            i += 1
        if secure and i > 1:
            raise PrecisionError("not enough precision to determine the number of roots")
        if i == 0:
            roots = [ ]
            P = self
        else:
            vali = self[i].valuation()
            prec = min((self[j].precision_absolute()-vali) / (i-j) for j in range(i))
            if prec is not Infinity:
                prec = prec.ceil()
            roots = [ (K(0,prec), i) ]
            P = self // self[:i+1]  # we do not shift because we need to track precision here

        # We use Newton polygon and slope factorisation to find roots
        vertices = P.newton_polygon().vertices(copy=False)
        deg = 0
        for i in range(1, len(vertices)):
            deg_left, val_left = vertices[i-1]
            deg_right, val_right = vertices[i]
            slope = (val_right - val_left) / (deg_left - deg_right)
            if slope not in ZZ or slope < minval:
                continue
            if hint is not None and slope == minval:
                rootsbar = hint
                if not rootsbar:
                    continue
            if i < len(vertices) - 1:
                F = P._factor_of_degree(deg_right - deg)
                P = P // F
            else:
                F = P
            if deg < deg_left:
                G = F._factor_of_degree(deg_left - deg)
                F //= G
            deg = deg_right
            val = F[0].valuation()
            if hint is None or slope != minval:
                Fbar = Pk([ F[j] >> (val - j*slope) for j in range(F.degree()+1) ])
                rootsbar = [ r for (r, _) in Fbar.roots() ]
                if not rootsbar:
                    continue
            rbar = rootsbar.pop()
            shift = K(rbar).lift_to_precision() << slope  # probably we should choose a better lift
            roots += [(r+shift, m) for (r, m) in F(x+shift)._roots(secure, slope, [r-rbar for r in rootsbar])]  # recursive call
        return roots


class Polynomial_generic_dense_cdv(Polynomial_generic_dense_inexact, Polynomial_generic_cdv):
    pass


class Polynomial_generic_sparse_cdv(Polynomial_generic_sparse, Polynomial_generic_cdv):
    pass


class Polynomial_generic_cdvr(Polynomial_generic_cdv):
    pass


class Polynomial_generic_dense_cdvr(Polynomial_generic_dense_cdv, Polynomial_generic_cdvr):
    pass


class Polynomial_generic_sparse_cdvr(Polynomial_generic_sparse_cdv, Polynomial_generic_cdvr):
    pass


class Polynomial_generic_cdvf(Polynomial_generic_cdv, Polynomial_generic_field):
    pass


class Polynomial_generic_dense_cdvf(Polynomial_generic_dense_cdv, Polynomial_generic_cdvf):
    pass


class Polynomial_generic_sparse_cdvf(Polynomial_generic_sparse_cdv, Polynomial_generic_cdvf):
    pass


############################################################################
# XXX:  Ensures that the generic polynomials implemented in Sage via PARI  #
# until at least until 4.5.0 unpickle correctly as polynomials implemented #
# via FLINT.                                                               #
try:
    from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
except ImportError:
    pass
else:
    from sage.misc.persist import register_unpickle_override
    register_unpickle_override('sage.rings.polynomial.polynomial_element_generic',
                               'Polynomial_rational_dense', Polynomial_rational_flint)
