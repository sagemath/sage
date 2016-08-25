"""
Univariate Polynomials over domains and fields

AUTHORS:

- William Stein: first version
- Martin Albrecht: Added singular coercion.
- David Harvey: split off polynomial_integer_dense_ntl.pyx (2007-09)
- Robert Bradshaw: split off polynomial_modn_dense_ntl.pyx (2007-09)

TESTS:

We test coercion in a particularly complicated situation::

    sage: W.<w>=QQ['w']
    sage: WZ.<z>=W['z']
    sage: m = matrix(WZ,2,2,[1,z,z,z^2])
    sage: a = m.charpoly()
    sage: R.<x> = WZ[]
    sage: R(a)
    x^2 + (-z^2 - 1)*x
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.rings.polynomial.polynomial_element import Polynomial, Polynomial_generic_dense, Polynomial_generic_dense_inexact
from sage.structure.element import IntegralDomainElement, EuclideanDomainElement

from sage.rings.polynomial.polynomial_singular_interface import Polynomial_singular_repr

from sage.libs.pari.all import pari_gen
from sage.structure.element import coerce_binop

from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.structure.factorization import Factorization


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
        <class 'sage.rings.polynomial.polynomial_element_generic.PolynomialRing_integral_domain_with_category.element_class'>
        sage: loads(f.dumps()) == f
        True

    A more extensive example::

        sage: A.<T> = PolynomialRing(Integers(5),sparse=True) ; f = T^2+1 ; B = A.quo(f)
        sage: C.<s> = PolynomialRing(B)
        sage: C
        Univariate Polynomial Ring in s over Univariate Quotient Polynomial Ring in Tbar over Ring of integers modulo 5 with modulus T^2 + 1
        sage: s + T
        s + Tbar
        sage: (s + T)**2
        s^2 + 2*Tbar*s + 4

    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        TESTS::

            sage: PolynomialRing(RIF, 'z', sparse=True)([RIF(-1, 1), RIF(-1,1)])
            0.?*z + 0.?
            sage: PolynomialRing(CIF, 'z', sparse=True)([CIF(RIF(-1,1), RIF(-1,1)), RIF(-1,1)])
            0.?*z + 0.? + 0.?*I
        """
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self.__coeffs = {}
            return
        R = parent.base_ring()
        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                x = dict(x.dict())
            elif x.parent() == R:
                x = {0:x}
            else:
                w = {}
                for n, c in x.dict().iteritems():
                    w[n] = R(c)
                # The following line has been added in trac ticket #9944.
                # Apparently, the "else" case has never occured before.
                x = w
        elif isinstance(x, list):
            x = dict((i, c) for (i, c) in enumerate(x) if c)
        elif isinstance(x, pari_gen):
            y = {}
            for i in range(len(x)):
                y[i] = R(x[i])
            x = y
            check = True
        elif not isinstance(x, dict):
            x = {0:x}   # constant polynomials
        if check:
            self.__coeffs = {}
            for i, z in x.iteritems():
                self.__coeffs[i] = R(z)
        else:
            self.__coeffs = x
        if check:
            self.__normalize()

    def dict(self):
        """
        Return a new copy of the dict of the underlying
        elements of ``self``.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: d = f.dict(); d
            {0: 5, 1997: 1, 10000: 7}
            sage: d[0] = 10
            sage: f.dict()
            {0: 5, 1997: 1, 10000: 7}
        """
        return dict(self.__coeffs)

    def coefficients(self,sparse=True):
        """
        Return the coefficients of the monomials appearing in ``self``.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: f.coefficients()
            [5, 1, 7]
        """
        if sparse:
          return [c[1] for c in sorted(self.__coeffs.iteritems())]
        else:
          return [self.__coeffs[i] if i in self.__coeffs else 0 for i in xrange(self.degree() + 1)]

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
        keys = self.__coeffs.keys()
        keys.sort()
        return keys

    def valuation(self):
        """
        Return the valuation of ``self``.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(GF(9,'a'), sparse=True)
            sage: f = w^1997 - w^10000
            sage: f.valuation()
            1997
            sage: R(19).valuation()
            0
            sage: R(0).valuation()
            +Infinity
        """
        c = self.__coeffs.keys()
        if len(c) == 0:
            return infinity
        return ZZ(min(self.__coeffs.keys()))

    def _derivative(self, var=None):
        """
        Computes formal derivative of this polynomial with respect to
        the given variable.

        If ``var`` is ``None`` or is the generator of this ring, the
        derivative is with respect to the generator. Otherwise,
        _derivative(var) is called recursively for each coefficient of
        this polynomial.

        .. seealso:: :meth:`.derivative`

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
        """
        P = self.parent()
        if var is not None and var is not P.gen():
            # call _derivative() recursively on coefficients
            return P(dict([(n, c._derivative(var)) \
                                     for (n, c) in self.__coeffs.iteritems()]))

        # compute formal derivative with respect to generator
        d = {}
        for n, c in self.__coeffs.iteritems():
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

        Check that :trac:`18600` is fixed::

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
            ZeroDivisionError: Inverse does not exist.
        """
        R = self.parent()
        # TODO:
        # calling the coercion model bin_op is much more accurate than using the
        # true division (which is bypassed by polynomials). But it does not work
        # in all cases!!
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        import operator
        try:
            Q = cm.bin_op(R.one(), ZZ.one(), operator.div).parent()
        except TypeError:
            F = (R.base_ring().one()/ZZ.one()).parent()
            Q = R.change_ring(F)

        if var is not None and var != R.gen():
            return Q({k:v.integral(var) for k,v in self.__coeffs.iteritems()}, check=False)

        return Q({ k+1:v/(k+1) for k,v in self.__coeffs.iteritems()}, check=False)

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

            sage: R.<w> = PolynomialRing(CDF, sparse=True)
            sage: f = CDF(1,2) + w^5 - CDF(pi)*w + CDF(e)
            sage: f._repr()   # abs tol 1e-15
            '1.0*w^5 - 3.141592653589793*w + 3.718281828459045 + 2.0*I'
            sage: f._repr(name='z')   # abs tol 1e-15
            '1.0*z^5 - 3.141592653589793*z + 3.718281828459045 + 2.0*I'

        TESTS::

            sage: pol = RIF['x']([0, 0, (-1,1)])
            sage: PolynomialRing(RIF, 'x', sparse=True)(pol)
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
        coeffs = sorted(self.__coeffs.iteritems())
        for (n, x) in reversed(coeffs):
            if x:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if s==" ":
            return "0"
        return s[1:]

    def __normalize(self):
        x = self.__coeffs
        D = [n for n, z in x.iteritems() if not z]
        for n in D:
            del x[n]

    def __getitem__(self,n):
        """
        Return the `n`-th coefficient of this polynomial.

        Negative indexes are allowed and always return 0 (so you can
        view the polynomial as embedding Laurent series).

        EXAMPLES::

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
            sage: R.<x> = PolynomialRing(RealField(19), sparse=True)
            sage: f = (2-3.5*x)^3; f
            -42.875*x^3 + 73.500*x^2 - 42.000*x + 8.0000

        Using slices, we can truncate polynomials::

            sage: f[:2]
            -42.000*x + 8.0000

        Any other kind of slicing is deprecated or an error::

            sage: f[1:3]
            doctest:...: DeprecationWarning: polynomial slicing with a start index is deprecated, use list() and slice the resulting list instead
            See http://trac.sagemath.org/18940 for details.
            73.500*x^2 - 42.000*x
            sage: f[1:3:2]
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomial slicing with a step is not defined
            sage: f["hello"]
            Traceback (most recent call last):
            ...
            TypeError: list indices must be integers, not str
        """
        if isinstance(n, slice):
            d = self.degree() + 1
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
            x = self.__coeffs
            v = {k: x[k] for k in x.keys() if start <= k < stop}
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

            sage: R.<z> = PolynomialRing(CC, sparse=True)
            sage: f = z^2 + CC.0; f
            1.00000000000000*z^2 + 1.00000000000000*I
            sage: f._unsafe_mutate(0, 10)
            sage: f
            1.00000000000000*z^2 + 10.0000000000000

        Much more nasty::

            sage: z._unsafe_mutate(1, 0)
            sage: z
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

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of ``self``.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(Integers(100), sparse=True)
            sage: f = 13*z^5 + 15*z^2 + 17*z
            sage: f.list()
            [0, 17, 15, 0, 0, 13]
        """
        zero = self.base_ring()(0)
        v = [zero] * (self.degree()+1)
        for n, x in self.__coeffs.iteritems():
            v[n] = x
        return v

    #def _pari_(self, variable=None):
    #    if variable is None:
    #        return self.__pari
    #    else:
    #        return self.__pari.subst('x',variable)

    def degree(self, gen=None):
        """
        Return the degree of this sparse polynomial.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(ZZ, sparse=True)
            sage: f = 13*z^50000 + 15*z^2 + 17*z
            sage: f.degree()
            50000
        """
        v = self.__coeffs.keys()
        if len(v) == 0:
            return -1
        return max(v)

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

        for (index, coeff) in right.__coeffs.iteritems():
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
        for (index, coeff) in self.__coeffs.iteritems():
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

        for (index1, coeff1) in self.__coeffs.iteritems():
            for (index2, coeff2) in right.__coeffs.iteritems():
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

        for (index, coeff) in self.__coeffs.iteritems():
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

        for (index, coeff) in self.__coeffs.iteritems():
            output[index] = coeff * right

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _cmp_(self, other):
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
            ....:     assert cmp(pd,qd) == cmp(Rs(pd), Rs(qd))
        """
        d1 = self.__coeffs
        keys1 = d1.keys()
        keys1.sort(reverse=True)

        d2 = other.__coeffs
        keys2 = d2.keys()
        keys2.sort(reverse=True)

        zero = self.base_ring().zero()

        if not keys1 and not keys2: return 0
        if not keys1: return -1
        if not keys2: return 1

        c = cmp(keys1[0], keys2[0])
        if c: return c
        c = cmp(d1[keys1[0]],d2[keys2[0]])
        if c: return c

        for k1, k2 in zip(keys1[1:], keys2[1:]):
            c = cmp(k1, k2)
            if c > 0:
                return cmp(d1[k1], zero)
            elif c < 0:
                return cmp(zero, d2[k2])
            c = cmp (d1[k1], d2[k2])
            if c: return c

        n1 = len(keys1)
        n2 = len(keys2)
        c = cmp(n1, n2)
        if c > 0: return cmp(d1[keys1[n2]], zero)
        elif c < 0: return cmp(zero, d2[keys2[n1]])
        return 0

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power `x^n`.

        If `n` is negative, terms below `x^n` will be discarded. Does
        not change this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^100000 + 2*x + 4
            sage: type(p)
            <class 'sage.rings.polynomial.polynomial_element_generic.PolynomialRing_integral_domain_with_category.element_class'>
            sage: p.shift(0)
             x^100000 + 2*x + 4
            sage: p.shift(-1)
             x^99999 + 2
            sage: p.shift(-100002)
             0
            sage: p.shift(2)
             x^100002 + 2*x^3 + 4*x^2

        TESTS:

        Check that :trac:`18600` is fixed::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^2^100 - 5
            sage: p.shift(10)
            x^1267650600228229401496703205386 - 5*x^10
            sage: p.shift(-10)
            x^1267650600228229401496703205366
            sage: p.shift(1.5)
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
            output = {index+n: coeff for index, coeff in self.__coeffs.iteritems()}
            return self.parent()(output, check=False)
        if n < 0:
            output = {index+n:coeff for index, coeff in self.__coeffs.iteritems() if index + n >= 0}
            return self.parent()(output, check=False)

    @coerce_binop
    def quo_rem(self, other):
        """
        Returns the quotient and remainder of the Euclidean division of
        ``self`` and ``other``.

        Raises ZerodivisionError if ``other`` is zero. Raises ArithmeticError
        if ``other`` has a nonunit leading coefficient.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(ZZ,sparse=True)
            sage: R.<y> = PolynomialRing(P,sparse=True)
            sage: f = R.random_element(10)
            sage: g = y^5+R.random_element(4)
            sage: q,r = f.quo_rem(g)
            sage: f == q*g + r and r.degree() < g.degree()
            True
            sage: g = x*y^5
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ArithmeticError: Division non exact (consider coercing to polynomials over the fraction field)
            sage: g = 0
            sage: f.quo_rem(g)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by zero polynomial

        TESTS::

            sage: P.<x> = PolynomialRing(ZZ,sparse=True)
            sage: f = x^10-4*x^6-5
            sage: g = 17*x^22+x^15-3*x^5+1
            sage: q,r = g.quo_rem(f)
            sage: g == f*q + r and r.degree() < f.degree()
            True
            sage: zero = P(0)
            sage: zero.quo_rem(f)
            (0, 0)
            sage: Q.<y> = IntegerModRing(14)[]
            sage: f = y^10-4*y^6-5
            sage: g = 17*y^22+y^15-3*y^5+1
            sage: q,r = g.quo_rem(f)
            sage: g == f*q + r and r.degree() < f.degree()
            True
            sage: f += 2*y^10 # 3 is invertible mod 14
            sage: q,r = g.quo_rem(f)
            sage: g == f*q + r and r.degree() < f.degree()
            True

        The following shows that :trac:`16649` is indeed fixed. ::

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
    def gcd(self,other,algorithm=None):
        """
        Return the gcd of this polynomial and ``other``

        INPUT:

        - ``other`` -- a polynomial defined over the same ring as this
          polynomial.

        ALGORITHM:

        Two algorithms are provided:

        - ``generic``: Uses the generic implementation, which depends on the
          base ring being a UFD or a field.
        - ``dense``: The polynomials are converted to the dense representation,
          their gcd is computed and is converted back to the sparse
          representation.

        Default is ``dense`` for polynomials over ZZ and ``generic`` in the
        other cases.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ,sparse=True)
            sage: p = x^6 + 7*x^5 + 8*x^4 + 6*x^3 + 2*x^2 + x + 2
            sage: q = 2*x^4 - x^3 - 2*x^2 - 4*x - 1
            sage: gcd(p,q)
            x^2 + x + 1
            sage: gcd(p, q, algorithm = "dense")
            x^2 + x + 1
            sage: gcd(p, q, algorithm = "generic")
            x^2 + x + 1
            sage: gcd(p, q, algorithm = "foobar")
            Traceback (most recent call last):
            ...
            ValueError: Unknown algorithm 'foobar'

        TESTS:

        Check that :trac:`19676` is fixed::

            sage: S.<y> = R[]
            sage: x.gcd(y)
            1
            sage: (6*x).gcd(9)
            3
        """

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.arith.all import lcm

        if algorithm is None:
            if self.base_ring() == ZZ:
                algorithm = "dense"
            else:
                algorithm = "generic"
        if algorithm=="dense":
            S = self.parent()
            # FLINT is faster but a bug makes the conversion extremely slow,
            # so NTL is used in those cases where the conversion is too slow. Cf
            # <https://groups.google.com/d/msg/sage-devel/6qhW90dgd1k/Hoq3N7fWe4QJ>
            sd = self.degree()
            od = other.degree()
            if max(sd,od)<100 or \
               min(len(self.__coeffs)/sd, len(other.__coeffs)/od)>.06:
                implementation="FLINT"
            else:
                implementation="NTL"
            D = PolynomialRing(S.base_ring(),'x',implementation=implementation)
            g = D(self).gcd(D(other))
            return S(g)
        elif algorithm=="generic":
            return Polynomial.gcd(self,other)
        else:
            raise ValueError("Unknown algorithm '%s'" % algorithm)

    def reverse(self, degree=None):
        """
        Return this polynomial but with the coefficients reversed.

        If an optional degree argument is given the coefficient list will be
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
            raise ValueError("degree argument must be a nonnegative integer, got %s"%degree)
        d = {degree-k: v for k,v in self.__coeffs.iteritems() if degree >= k}
        return self.parent()(d, check=False)

    def truncate(self, n):
        """
        Return the polynomial of degree `< n` equal to `self` modulo `x^n`.

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

            sage: R.<x> = PolynomialRing(ZZ,sparse=True)
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
        Returns a tuple (quotient, remainder) where
            self = quotient * other + remainder.

        EXAMPLES::

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
            diff_deg=R.degree()-B.degree()
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
        <class 'sage.rings.polynomial.polynomial_element_generic.PolynomialRing_field_with_category.element_class'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
        Polynomial_generic_sparse.__init__(self, parent, x, check, is_gen)


class Polynomial_generic_dense_field(Polynomial_generic_dense, Polynomial_generic_field):
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
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
        Returns a list of the Newton slopes of this polynomial.

        These are the valuations of the roots of this polynomial.

        If ``repetition`` is ``True``, each slope is repeated a number of
        times equal to its multiplicity. Otherwise it appears only
        one time.

        EXAMPLES::

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
        Returns a list of vertices of the Newton polygon of this polynomial.

        .. NOTE::

            If some coefficients have not enough precision an error is raised.

        EXAMPLES::

            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10
            sage: f.newton_polygon()
            Finite Newton polygon with 4 vertices: (0, 1), (1, 0), (4, 0), (10, 2)

            sage: g = f + K(0,0)*t^4; g
            (5^2 + O(5^22))*t^10 + (O(5^0))*t^4 + (3 + O(5^20))*t + (5 + O(5^21))
            sage: g.newton_polygon()
            Traceback (most recent call last):
            ...
            PrecisionError: The coefficient of t^4 has not enough precision

        AUTHOR:

        - Xavier Caruso (2013-03-20)
        """
        d = self.degree()
        from sage.geometry.newton_polygon import NewtonPolygon
        polygon = NewtonPolygon([(x, self[x].valuation()) for x in range(d+1)])
        polygon_prec = NewtonPolygon([ (x, self[x].precision_absolute()) for x in range(d+1) ])
        vertices = polygon.vertices(copy=False)
        vertices_prec = polygon_prec.vertices(copy=False)
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

            sage: K = Qp(5, 10)
            sage: P.<x> = PolynomialRing(K)
            sage: f = x^2 + 1
            sage: root = f.hensel_lift(2); root
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: f(root)
            O(5^10)

            sage: g = (x^2 + 1)*(x - 7)
            sage: g.hensel_lift(2)  # here, 2 is a multiple root modulo p
            Traceback (most recent call last):
            ...
            ValueError: a is not close enough to a root of this polynomial

        AUTHOR:

        - Xavier Caruso (2013-03-23)
        """
        base = self.base_ring()
        selfa = self(a)
        der = self.derivative()
        dera = der(a)
        if selfa.valuation() <= 2 * dera.valuation():
            raise ValueError("a is not close enough to a root of this polynomial")
        # Newton iteration
        # Todo: compute everything up to the adequate precision at each step
        b = ~dera
        while(True):
            na = a - selfa * b
            if na == a: return a
            a = na
            selfa = self(a)
            dera = der(a)
            b *= 2 - dera*b

    def _factor_of_degree(self, deg):
        """
        Return a factor of ``self`` of degree ``deg``.

        Algorithm is Newton iteration.

        This fails if ``deg`` is not a breakpoint in the Newton
        polygon of ``self``.

        Only for internal use!

        EXAMPLES::

            sage: K = Qp(5)
            sage: R.<x> = K[]
            sage: K = Qp(5)
            sage: R.<t> = K[]
            sage: f = 5 + 3*t + t^4 + 25*t^10

            sage: g = f._factor_of_degree(4)
            sage: (f % g).is_zero()
            True

            sage: g = f._factor_of_degree(3)    # not tested
            Traceback (most recent call last)
            ...
            KeyboardInterrupt:

        AUTHOR:

        - Xavier Caruso (2013-03-20)

        TODO:

        Precision is not optimal, and can be improved.
        """
        coeffs = self.list()
        a = self.truncate(deg + 1)
        b = v = self.parent()(1)
        x = self % a
        while(not x.is_zero()):
            a += (v * x) % a
            b, x = self.quo_rem(a)
            b %= a
            v = (v * (2 - b*v)) % a

        return a

    def factor_of_slope(self, slope=None):
        """
        INPUT:

        -  slope -- a rational number (default: the first slope
           in the Newton polygon of ``self``)

        OUTPUT:

        The factor of ``self`` corresponding to the slope ``slope`` (i.e.
        the unique monic divisor of ``self`` whose slope is ``slope`` and
        degree is the length of ``slope`` in the Newton polygon).

        EXAMPLES::

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

            sage: f.factor_of_slope(-1)
            (1 + O(5^20))

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
# XXX:  Ensures that the generic polynomials implemented in SAGE via PARI  #
# until at least until 4.5.0 unpickle correctly as polynomials implemented #
# via FLINT.                                                               #
from sage.structure.sage_object import register_unpickle_override
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint

register_unpickle_override( \
    'sage.rings.polynomial.polynomial_element_generic', \
    'Polynomial_rational_dense', Polynomial_rational_flint)
