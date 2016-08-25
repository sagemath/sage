r"""
Elements of Laurent polynomial rings
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.rings.integer cimport Integer
from sage.structure.element import is_Element, coerce_binop
from sage.misc.latex import latex
import sage.misc.latex
from sage.misc.misc import union
from sage.structure.factorization import Factorization
from sage.misc.derivative import multi_derivative
from sage.rings.polynomial.polynomial_element import Polynomial


cdef class LaurentPolynomial_generic(CommutativeAlgebraElement):
    """
    A generic Laurent polynomial.
    """
    pass


cdef class LaurentPolynomial_univariate(LaurentPolynomial_generic):
    """
    A univariate Laurent polynomial in the form of `t^n \cdot f`
    where `f` is a polynomial in `t`.

    INPUT:

    - ``parent`` -- a Laurent polynomial ring

    - ``f`` -- a polynomial (or something can be coerced to one)

    - ``n`` -- (default: 0) an integer

    AUTHORS:

    - Tom Boothby (2011) copied this class almost verbatim from
      ``laurent_series_ring_element.pyx``, so most of the credit goes to
      William Stein, David Joyner, and Robert Bradshaw
    - Travis Scrimshaw (09-2013): Cleaned-up and added a few extra methods
    """
    def __init__(self, parent, f, n=0):
        r"""
        Create the Laurent polynomial `t^n \cdot f`.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: R([1,2,3])
            1 + 2*q + 3*q^2
            sage: TestSuite(q^-3 + 3*q + 2).run()

        ::

            sage: S.<s> = LaurentPolynomialRing(GF(5))
            sage: T.<t> = PolynomialRing(pAdicRing(5))
            sage: S(t)
            s
            sage: parent(S(t))
            Univariate Laurent Polynomial Ring in s over Finite Field of size 5
            sage: parent(S(t)[1])
            Finite Field of size 5
        """
        CommutativeAlgebraElement.__init__(self, parent)

        if isinstance(f, LaurentPolynomial_univariate):
            n += (<LaurentPolynomial_univariate>f).__n
            if (<LaurentPolynomial_univariate>f).__u._parent is parent.polynomial_ring():
                f = (<LaurentPolynomial_univariate>f).__u
            else:
                f = parent.polynomial_ring()((<LaurentPolynomial_univariate>f).__u)
        elif (not isinstance(f, Polynomial)) or (parent is not f.parent()):
            if isinstance(f, dict):
                v = min(f.keys())
                f = dict((i-v,c) for i,c in f.items())
                n += v
            f = parent.polynomial_ring()(f)

        # self is that t^n * u:
        cdef long val
        self.__u = f
        self.__n = n
        self.__normalize()

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: elt = q^-3 + 2 + q
            sage: loads(dumps(elt)) == elt
            True
        """
        return LaurentPolynomial_univariate, (self._parent, self.__u, self.__n)

    def change_ring(self, R):
        """
        Return a copy of this Laurent polynomial, with coefficients in ``R``.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: a = x^2 + 3*x^3 + 5*x^-1
            sage: a.change_ring(GF(3))
            2*x^-1 + x^2
        """
        return self.parent().change_ring(R)(self)

    def is_unit(self):
        """
        Return ``True`` if this Laurent polynomial is a unit in this ring.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: (2+t).is_unit()
            False
            sage: f = 2*t
            sage: f.is_unit()
            True
            sage: 1/f
            1/2*t^-1
            sage: R(0).is_unit()
            False
            sage: R.<s> = LaurentPolynomialRing(ZZ)
            sage: g = 2*s
            sage: g.is_unit()
            False
            sage: 1/g
            1/2*s^-1

        ALGORITHM: A Laurent polynomial is a unit if and only if its "unit
        part" is a unit.
        """
        return self.__u.is_term() and self.__u.coefficients()[0].is_unit()

    def is_zero(self):
        """
        Return ``1`` if ``self`` is 0, else return ``0``.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x + x + x^2 + 3*x^4
            sage: f.is_zero()
            0
            sage: z = 0*f
            sage: z.is_zero()
            1
        """
        return self.__u.is_zero()

    def __nonzero__(self):
        """
        Check if ``self`` is non-zero.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x + x + x^2 + 3*x^4
            sage: not f
            False
            sage: z = 0*f
            sage: not z
            True
        """
        return not self.__u.is_zero()

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of ``self`` under the morphism defined by
        ``im_gens`` in ``codomain``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: H = Hom(R, QQ)
            sage: mor = H(2)
            sage: mor(t^2 + t^-2)
            17/4
            sage: 4 + 1/4
            17/4
        """
        return codomain(self(im_gens[0]))

    def __normalize(self):
        r"""
        A Laurent series is a pair `(u(t), n)`, where either `u = 0`
        (to some precision) or `u` is a unit. This pair corresponds to
        `t^n \cdot u(t)`.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: elt = t^2 + t^4 # indirect doctest
            sage: elt.polynomial_construction()
            (t^2 + 1, 2)
        """
        from sage.rings.infinity import infinity
        if self.is_zero():
            return
        v = self.__u.valuation()
        if v != 0 and v != infinity:
            self.__n += v
            self.__u = self.__u >> v

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: 2 + (2/3)*t^3
            2 + 2/3*t^3
        """
        if self.is_zero():
            return "0"
        s = " "
        v = self.__u.list()
        valuation = self.__n
        m = len(v)
        X = self._parent.variable_name()
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        first = True
        for n in xrange(m):
            x = v[n]
            e = n + valuation
            x = str(x)
            if x != '0':
                if not first:
                    s += " + "
                if not atomic_repr and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "({})".format(x)
                if e == 1:
                    var = "*{}".format(X)
                elif e == 0:
                    var = ""
                else:
                    var = "*{}^{}".format(X,e)
                s += "{}{}".format(x,var)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        return s[1:]

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = (17/2)*x^-2 + x + x^2 + 3*x^4
            sage: latex(f)
            \frac{\frac{17}{2}}{x^{2}} + x + x^{2} + 3x^{4}

        Verify that :trac:`6656` has been fixed::

            sage: R.<a,b>=PolynomialRing(QQ)
            sage: T.<x>=LaurentPolynomialRing(R)
            sage: y = a*x+b*x
            sage: y._latex_()
            '\\left(a + b\\right)x'
            sage: latex(y)
            \left(a + b\right)x
        """
        if self.is_zero():
            return "0"
        s = " "
        v = self.__u.list()
        valuation = self.__n
        m = len(v)
        X = self._parent.latex_variable_names()[0]
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        first = True
        for n in xrange(m):
            x = v[n]
            e = n + valuation
            x = sage.misc.latex.latex(x)
            if x != '0':
                if not first:
                    s += " + "
                if not atomic_repr and e > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "\\left({}\\right)".format(x)
                if e == 1:
                    var = "|{}".format(X)
                elif e == 0:
                    var = ""
                elif e > 0:
                    var = "|{}^{{{}}}".format(X,e)
                if e >= 0:
                    s += "{}{}".format(x,var)
                else: # negative e
                    if e == -1:
                        s += "\\frac{{{}}}{{{}}}".format(x, X)
                    else:
                        s += "\\frac{{{}}}{{{}^{{{}}}}}".format(x, X,-e)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1|"," ")
        s = s.replace(" -1|", " -")
        s = s.replace("|","")

        return s[1:]

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = -5/t^(10) + t + t^2 - 10/3*t^3
            sage: hash(f) == hash(f)
            True
        """
        return hash(self.__u) ^ self.__n

    def __getitem__(self, i):
        """
        Return the `i`-th coefficient of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = -5/t^(10) + t + t^2 - 10/3*t^3; f
            -5*t^-10 + t + t^2 - 10/3*t^3
            sage: f[-10]
            -5
            sage: f[1]
            1
            sage: f[3]
            -10/3
            sage: f[-9]
            0
            sage: f = -5/t^(10) + 1/3 + t + t^2 - 10/3*t^3; f
            -5*t^-10 + 1/3 + t + t^2 - 10/3*t^3

        Slicing is deprecated::

            sage: f[-10:2]
            doctest:...: DeprecationWarning: polynomial slicing with a start index is deprecated, use list() and slice the resulting list instead
            See http://trac.sagemath.org/18940 for details.
            -5*t^-10 + 1/3 + t
            sage: f[0:]
            1/3 + t + t^2 - 10/3*t^3
            sage: f[:3]
            -5*t^-10 + 1/3 + t + t^2
            sage: f[-14:5:2]
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomial slicing with a step is not defined
        """
        if isinstance(i, slice):
            start = i.start - self.__n if i.start is not None else 0
            stop = i.stop - self.__n if i.stop is not None else self.__u.degree() + 1
            f = self.__u[start:stop:i.step]  # deprecation(18940)
            return LaurentPolynomial_univariate(self._parent, f, self.__n)

        return self.__u[i - self.__n]

    def __iter__(self):
        """
        Iterate through the coefficients from the first nonzero one to the
        last nonzero one.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3; f
            -5*t^-2 + t + t^2 - 10/3*t^3
            sage: for a in f: print(a)
            -5
            0
            0
            1
            1
            -10/3
        """
        return iter(self.__u)

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + 2/x
            sage: g = f._symbolic_(SR); g
            (x^4 + 2)/x
            sage: g(x=2)
            9

            sage: g = SR(f)
            sage: g(x=2)
            9

        The polynomial does not have to be over a field of
        characteristic 0::

            sage: R.<w> = LaurentPolynomialRing(GF(7))
            sage: f = SR(2*w^3 + 1); f
            2*w^3 + 1
            sage: f.variables()
            (w,)
        """
        d = dict([(repr(g), R.var(g)) for g in self.parent().gens()])
        return self.subs(**d)

    def dict(self):
        """
        Return a dictionary representing ``self``.

        EXAMPLES::
            sage: R.<x,y> = ZZ[]
            sage: Q.<t> = LaurentPolynomialRing(R)
            sage: f = (x^3 + y/t^3)^3 + t^2; f
            y^3*t^-9 + 3*x^3*y^2*t^-6 + 3*x^6*y*t^-3 + x^9 + t^2
            sage: f.dict()
            {-9: y^3, -6: 3*x^3*y^2, -3: 3*x^6*y, 0: x^9, 2: 1}
        """
        return dict(zip(self.exponents(), self.coefficients()))

    def coefficients(self):
        """
        Return the nonzero coefficients of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3
            sage: f.coefficients()
            [-5, 1, 1, -10/3]
        """
        return self.__u.coefficients()

    def exponents(self):
        """
        Return the exponents appearing in ``self`` with nonzero coefficients.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3
            sage: f.exponents()
            [-2, 1, 2, 3]
        """
        return [i+self.__n for i in self.__u.exponents()]

    def __setitem__(self, n, value):
        """
        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = t^2 + t^-3
            sage: f[2] = 5
            Traceback (most recent call last):
            ...
            IndexError: Laurent polynomials are immutable
        """
        raise IndexError("Laurent polynomials are immutable")

    def _unsafe_mutate(self, i, value):
        """
        Sage assumes throughout that commutative ring elements are
        immutable. This is relevant for caching, etc. But sometimes you
        need to change a Laurent polynomial and you really know what you're
        doing. That's when this function is for you.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = t^2 + t^-3
            sage: f._unsafe_mutate(2, 3)
            sage: f
            t^-3 + 3*t^2
        """
        j = i - self.__n
        if j >= 0:
            self.__u._unsafe_mutate(j, value)
        else: # off to the left
            if value != 0:
                self.__n = self.__n + j
                R = self._parent.base_ring()
                coeffs = [value] + [R(0) for _ in range(1,-j)] + self.__u.list()
                self.__u = self.__u._parent(coeffs)
        self.__normalize()

    cpdef _add_(self, right_m):
        """
        Add two Laurent polynomials with the same parent.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: t + t
            2*t
            sage: f = 1/t + t^2 + t^3 - 17/3 * t^4
            sage: g = 2/t + t^3
            sage: f + g
            3*t^-1 + t^2 + 2*t^3 - 17/3*t^4
            sage: f + 0
            t^-1 + t^2 + t^3 - 17/3*t^4
            sage: 0 + f
            t^-1 + t^2 + t^3 - 17/3*t^4
            sage: R(0) + R(0)
            0
            sage: t^3 + t^-3
            t^-3 + t^3

        ALGORITHM: Shift the unit parts to align them, then add.
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate>right_m
        cdef long m

        # 1. Special case when one or the other is 0.
        if not right:
            return self
        if not self:
            return right

        # 2. Align the unit parts.
        if self.__n < right.__n:
            m = self.__n
            f1 = self.__u
            f2 = right.__u << right.__n - m
        elif self.__n > right.__n:
            m = right.__n
            f1 = self.__u << self.__n - m
            f2 = right.__u
        else:
            m = self.__n
            f1 = self.__u
            f2 = right.__u
        # 3. Add
        return LaurentPolynomial_univariate(self._parent, f1 + f2, m)

    cpdef _sub_(self, right_m):
        """
        Subtract two Laurent polynomials with the same parent.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: t - t
            0
            sage: t^5 + 2 * t^-5
            2*t^-5 + t^5

        ALGORITHM: Shift the unit parts to align them, then subtract.
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate>right_m
        cdef long m

        # 1. Special case when one or the other is 0.
        if not right:
            return self
        if not self:
            return -right

        # 2. Align the unit parts.
        if self.__n < right.__n:
            m = self.__n
            f1 = self.__u
            f2 = right.__u << right.__n - m
        else:
            m = right.__n
            f1 = self.__u << self.__n - m
            f2 = right.__u
        # 3. Subtract
        return LaurentPolynomial_univariate(self._parent, f1 - f2, m)

    def degree(self):
        """
        Return the degree of this polynomial.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: g = x^2 - x^4
            sage: g.degree()
            4
            sage: g = -10/x^5 + x^2 - x^7
            sage: g.degree()
            7
        """
        return self.__u.degree() + self.__n

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: -(1+t^5)
            -1 - t^5
        """
        return LaurentPolynomial_univariate(self._parent, -self.__u, self.__n)

    cpdef _mul_(self, right_r):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(GF(2))
            sage: f = 1/x^3 + x + x^2 + 3*x^4
            sage: g = 1 - x + x^2 - x^4
            sage: f*g
            x^-3 + x^-2 + x^-1 + x^8
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate>right_r
        return LaurentPolynomial_univariate(self._parent,
                             self.__u * right.__u,
                             self.__n + right.__n)

    cpdef _rmul_(self, RingElement c):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = 1/x^3 + x + x^2 + 3*x^4
            sage: 3 * f
            3*x^-3 + 3*x + 3*x^2 + 9*x^4
        """
        return LaurentPolynomial_univariate(self._parent, self.__u._rmul_(c), self.__n)

    cpdef _lmul_(self, RingElement c):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = 1/x^3 + x + x^2 + 3*x^4
            sage: f * 3
            3*x^-3 + 3*x + 3*x^2 + 9*x^4
        """
        return LaurentPolynomial_univariate(self._parent, self.__u._lmul_(c), self.__n)

    def is_monomial(self):
        """
        Return True if this element is a monomial.  That is, if self is
        `x^n` for some integer `n`.

        EXAMPLES::

            sage: k.<z> = LaurentPolynomialRing(QQ)
            sage: z.is_monomial()
            True
            sage: k(1).is_monomial()
            True
            sage: (z+1).is_monomial()
            False
            sage: (z^-2909).is_monomial()
            True
            sage: (38*z^-2909).is_monomial()
            False
        """

        return self.__u.is_monomial()

    def __pow__(_self, r, dummy):
        """
        EXAMPLES::

            sage: x = LaurentPolynomialRing(QQ,'x').0
            sage: f = x + x^2 + 3*x^4
            sage: g = 1/x^10 - x
            sage: f^3
            x^3 + 3*x^4 + 3*x^5 + 10*x^6 + 18*x^7 + 9*x^8 + 27*x^9 + 27*x^10 + 27*x^12
            sage: g^4
            x^-40 - 4*x^-29 + 6*x^-18 - 4*x^-7 + x^4
        """
        cdef LaurentPolynomial_univariate self = _self
        right = int(r)
        if right != r:
            raise ValueError("exponent must be an integer")
        return LaurentPolynomial_univariate(self._parent, self.__u**right, self.__n*right)

    cpdef _floordiv_(self, rhs):
        """
        Perform division with remainder and return the quotient.

        EXAMPLES::

            sage: L.<x> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + x^-3
            sage: g = x^-1 + x
            sage: f // g
            x^-2 - 1 + x^2
            sage: g * (f // g) == f
            True
            sage: f // 1
            x^-3 + x^3
            sage: 1 // f
            0
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate> rhs
        return LaurentPolynomial_univariate(self._parent,
                                            self.__u // right.__u,
                                            self.__n - right.__n)

    def shift(self, k):
        r"""
        Return this Laurent polynomial multiplied by the power `t^n`.
        Does not change this polynomial.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ['y'])
            sage: f = (t+t^-1)^4; f
            t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
            sage: f.shift(10)
            t^6 + 4*t^8 + 6*t^10 + 4*t^12 + t^14
            sage: f >> 10
            t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6
            sage: f << 4
            1 + 4*t^2 + 6*t^4 + 4*t^6 + t^8
        """
        return LaurentPolynomial_univariate(self._parent, self.__u, self.__n + k)

    def __lshift__(LaurentPolynomial_univariate self, k):
        """
        Return the left shift of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = (t+t^-1)^4; f
            t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
            sage: f << 4
            1 + 4*t^2 + 6*t^4 + 4*t^6 + t^8
        """
        return LaurentPolynomial_univariate(self._parent, self.__u, self.__n + k)

    def __rshift__(LaurentPolynomial_univariate self, k):
        """
        Return the right shift of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = (t+t^-1)^4; f
            t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
            sage: f >> 10
            t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6
        """
        return LaurentPolynomial_univariate(self._parent, self.__u, self.__n - k)

    cpdef _div_(self, rhs):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = x + x^2 + 3*x^4
            sage: g = 1/x^7 - x + x^2 - x^4
            sage: f / x
            1 + x + 3*x^3
            sage: f / g
            (3*x^11 + x^9 + x^8)/(-x^11 + x^9 - x^8 + 1)
            sage: (x^-2 + x)*(x^-2 + 1) / ((x^5 + x^8)*(x + 2))
            (x^2 + 1)/(x^10 + 2*x^9)
            sage: (x^-2 + x)*(x^-2 + 1) / ((x^-5 + x^-8)*(x + 2))
            (x^6 + x^4)/(x + 2)
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate> rhs
        if right.__u.is_zero():
            raise ZeroDivisionError
        return self * ~right

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: i = ~(t^-2); i
            t^2
            sage: i.parent() is R
            True
            sage: i = ~(2*t^2); i
            1/2*t^-2
            sage: i.parent() is R
            True
            sage: i = ~(t^-2 + 2 + t^2); i
            t^2/(t^4 + 2*t^2 + 1)
            sage: i.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        if self.__u.is_constant(): # this has a single term c*x^n
            if self.__u.is_unit():
                return LaurentPolynomial_univariate(self._parent, self.__u.inverse_of_unit(), -self.__n)
            # Enlarge the ring so we can divide by the coefficient
            R = self.base_ring().fraction_field()
            P = self.parent().change_ring(R)
            return LaurentPolynomial_univariate(P, ~R(self.__u), -self.__n)
        P = self.parent().polynomial_ring()
        if self.__n < 0:
            return P.gen()**-self.__n / self.__u
        return P.one() / (P.gen()**self.__n * self.__u)

    def inverse_of_unit(self):
        """
        Return the inverse of ``self`` if a unit.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: (t^-2).inverse_of_unit()
            t^2
            sage: (t + 2).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: element is not a unit
        """
        if self.is_unit():
            return ~self
        raise ArithmeticError("element is not a unit")

    def _fraction_pair(self):
        """
        Return one representation of ``self`` as a pair
        ``(numerator, denominator)``.

        Here both the numerator and the denominator are polynomials.

        This is used for coercion into the fraction field.

        EXAMPLES::

            sage: L.<x> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^-7 + 3*x^3 + 1 + 2*x^4 + x^6
            sage: f._fraction_pair()
            (x^13 + 2*x^11 + 3*x^10 + x^7 + 4, x^7)
        """
        P = self.parent().polynomial_ring()
        numer = self.__u
        denom = P.one()
        if self.__n > 0:
            numer *= P.gen()**self.__n
        elif self.__n < 0:
            denom *= P.gen()**-self.__n
        return (numer, denom)

    def gcd(self, right):
        """
        Return the gcd of ``self`` with ``right`` where the common divisor
        ``d`` makes both ``self`` and ``right`` into polynomials with
        the lowest possible degree.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: t.gcd(2)
            1
            sage: gcd(t^-2 + 1, t^-4 + 3*t^-1)
            t^-4
            sage: gcd((t^-2 + t)*(t + t^-1), (t^5 + t^8)*(1 + t^-2))
            t^-3 + t^-1 + 1 + t^2
        """
        b = <LaurentPolynomial_univariate>self._parent(right)
        return LaurentPolynomial_univariate(self._parent, self.__u.gcd(b.__u), min(self.__n, b.__n))

    @coerce_binop
    def quo_rem(self, right_r):
        """
        Attempts to divide ``self`` by ``right`` and returns a quotient and
        a remainder.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: (t^-3 - t^3).quo_rem(t^-1 - t)
            (t^-2 + 1 + t^2, 0)
            sage: (t^-2 + 3 + t).quo_rem(t^-4)
            (t^2 + 3*t^4 + t^5, 0)
            sage: (t^-2 + 3 + t).quo_rem(t^-4 + t)
            (0, 1 + 3*t^2 + t^3)
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate>right_r
        q,r = self.__u.quo_rem(right.__u)
        ql = LaurentPolynomial_univariate(self._parent, q, self.__n - right.__n)
        rl = LaurentPolynomial_univariate(self._parent, r, 0)
        return (ql, rl)

    cpdef int _cmp_(self, right_r) except -2:
        r"""
        Comparison of ``self`` and ``right_r``.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = x^(-1) + 1 + x
            sage: g = x^(-1) + 1
            sage: f == g
            False

        ::

            sage: f = x^(-1) + 1 + x
            sage: g = x^(-1) + 2
            sage: f == g
            False
            sage: f != g
            True
            sage: f < g
            True
            sage: f <= g
            True
            sage: f > g
            False
            sage: f >= g
            False

        ::

            sage: f = x^(-2) + 1 + x
            sage: g = x^(-1) + 2
            sage: f == g
            False
            sage: f < g
            False
            sage: f > g
            True
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate>right_r

        zero = self.base_ring()(0)

        if not self and not right:
            if self.__n < right.__n:
                return cmp(self.__u[0], zero)
            elif self.__n > right.__n:
                return cmp(zero, right.__u[0])

        # zero pad coefficients on the left, to line them up for comparison
        cdef long n = min(self.__n, right.__n)
        x = [zero] * (self.__n - n) + self.__u.list()
        y = [zero] * (right.__n - n) + right.__u.list()

        # zero pad on right to make the lists the same length
        # (this is necessary since the power series list() function just
        # returns the coefficients of the underlying polynomial, which may
        # have zeroes in the high coefficients)
        if len(x) < len(y):
            x.extend([zero] * (len(y) - len(x)))
        elif len(y) < len(x):
            y.extend([zero] * (len(x) - len(y)))

        return cmp(x,y)

    def valuation(self, p=None):
        """
        Return the valuation of ``self``.

        The valuation of a Laurent polynomial `t^n u` is `n` plus the
        valuation of `u`.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: g = 1 - x + x^2 - x^4
            sage: f.valuation()
            -1
            sage: g.valuation()
            0
        """
        return self.__n + self.__u.valuation(p)

    def truncate(self, n):
        """
        Return a polynomial with degree at most `n-1` whose `j`-th coefficients
        agree with ``self`` for all `j < n`.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x^12 + x^3 + x^5 + x^9
            sage: f.truncate(10)
            x^-12 + x^3 + x^5 + x^9
            sage: f.truncate(5)
            x^-12 + x^3
            sage: f.truncate(-16)
            0
        """
        if n <= self.valuation():
            return self._parent(0)
        return LaurentPolynomial_univariate(self._parent, self.__u.truncate(n-self.__n), self.__n)

    def variable_name(self):
        """
        Return the name of variable of ``self`` as a string.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: f.variable_name()
            'x'
        """
        return self._parent.variable_name()

    def variables(self):
        """
        Return the tuple of variables occuring in this Laurent polynomial.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: f.variables()
            (x,)
            sage: R.one().variables()
            ()
        """
        if self.is_constant():
            return ()
        return self._parent.gens()

    def polynomial_construction(self):
        """
        Return the polynomial and the shift in power used to construct the
        Laurent polynomial `t^n u`.

        OUTPUT:

        A tuple ``(u, n)`` where ``u`` is the underlying polynomial and ``n``
        is the power of the exponent shift.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: f.polynomial_construction()
            (3*x^5 + x^3 + 1, -1)
        """
        return (self.__u, self.__n)

    def is_constant(self):
        """
        Return ``True`` if ``self`` is constant.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: x.is_constant()
            False
            sage: R.one().is_constant()
            True
            sage: (x^-2).is_constant()
            False
            sage: (x^2).is_constant()
            False
            sage: (x^-2 + 2).is_constant()
            False
        """
        return self.__n == 0 and self.__u.is_constant()

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: cf = copy(f)
            sage: cf == f
            True
            sage: cf is not f
            True
        """
        from copy import copy
        return LaurentPolynomial_univariate(self._parent, copy(self.__u), self.__n)

    def derivative(self, *args):
        """
        The formal derivative of this Laurent polynomial, with respect
        to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global :func`derivative()` function for more
        details.

        .. SEEALSO::

           :meth:`_derivative`

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: g = 1/x^10 - x + x^2 - x^4
            sage: g.derivative()
            -10*x^-11 - 1 + 2*x - 4*x^3
            sage: g.derivative(x)
            -10*x^-11 - 1 + 2*x - 4*x^3

        ::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = LaurentPolynomialRing(R)
            sage: f = 2*t/x + (3*t^2 + 6*t)*x
            sage: f.derivative()
            -2*t*x^-2 + (3*t^2 + 6*t)
            sage: f.derivative(x)
            -2*t*x^-2 + (3*t^2 + 6*t)
            sage: f.derivative(t)
            2*x^-1 + (6*t + 6)*x
        """
        return multi_derivative(self, args)

    def _derivative(self, var=None):
        """
        The formal derivative of this Laurent series with respect to ``var``.

        If ``var`` is ``None`` or the generator of this ring, it's the formal
        derivative as expected. Otherwise, ``_derivative(var)`` gets called
        recursively on each coefficient.

        .. SEEALSO::

           :meth:`derivative`

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = x^2 + 3*x^4
            sage: f._derivative()
            2*x + 12*x^3
            sage: f._derivative(x)
            2*x + 12*x^3
            sage: g = 1/x^10 - x + x^2 - x^4
            sage: g._derivative()
            -10*x^-11 - 1 + 2*x - 4*x^3

        Differentiating with respect to something other than the generator
        gets recursed into the base ring::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = LaurentPolynomialRing(R)
            sage: f = 2*t/x + (3*t^2 + 6*t)*x
            sage: f._derivative(t)
            2*x^-1 + (6*t + 6)*x
        """
        if var is not None and var is not self._parent.gen():
            # call _derivative() recursively on coefficients
            u = [coeff._derivative(var) for coeff in self.__u.list()]
            u = self._parent.polynomial_ring()(u)
            return LaurentPolynomial_univariate(self._parent, u, self.__n)

        # compute formal derivative with respect to generator
        if self.is_zero():
            return LaurentPolynomial_univariate(self._parent, 0)
        cdef long m, n = self.__n
        a = self.__u.list()
        v = [(n+m)*a[m] for m from 0 <= m < len(a)]
        u = self._parent.polynomial_ring()(v)
        return LaurentPolynomial_univariate(self._parent, u, n-1)

    def integral(self):
        r"""
        The formal integral of this Laurent series with 0 constant term.

        EXAMPLES:

        The integral may or may not be defined if the base ring
        is not a field.

        ::

            sage: t = LaurentPolynomialRing(ZZ, 't').0
            sage: f = 2*t^-3 + 3*t^2
            sage: f.integral()
            -t^-2 + t^3

        ::

            sage: f = t^3
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: coefficients of integral cannot be coerced into the base ring

        The integral of `1/t` is `\log(t)`, which is not given by a
        Laurent polynomial::

            sage: t = LaurentPolynomialRing(ZZ,'t').0
            sage: f = -1/t^3 - 31/t
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: the integral of is not a Laurent polynomial, since t^-1 has nonzero coefficient

        Another example with just one negative coefficient::

            sage: A.<t> = LaurentPolynomialRing(QQ)
            sage: f = -2*t^(-4)
            sage: f.integral()
            2/3*t^-3
            sage: f.integral().derivative() == f
            True
        """
        cdef long i, n = self.__n
        a = self.__u.list()
        if self[-1] != 0:
            raise ArithmeticError(
                  "the integral of is not a Laurent polynomial, since t^-1 has nonzero coefficient")

        if n < 0:
            v = [a[i]/(n+i+1) for i in range(min(-1-n,len(a)))] + [0]
        else:
            v = []
        v += [a[i]/(n+i+1) for i in range(max(-n,0), len(a))]
        try:
            u = self._parent.polynomial_ring()(v)
        except TypeError:
            raise ArithmeticError("coefficients of integral cannot be coerced into the base ring")
        return LaurentPolynomial_univariate(self._parent, u, n+1)

    def __call__(self, *x):
        """
        Compute value of this Laurent polynomial at ``x``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: f = t^(-2) + t^2
            sage: f(2)
            17/4
            sage: f(-1)
            2
            sage: f(1/3)
            82/9
        """
        if isinstance(x[0], tuple):
            x = x[0]
        return self.__u(x) * (x[0]**self.__n)

    def factor(self):
        """
        Return a Laurent monomial (the unit part of the factorization) and
        a factored polynomial.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: f = 4*t^-7 + 3*t^3 + 2*t^4 + t^-6
            sage: f.factor()
            (t^-7) * (4 + t + 3*t^10 + 2*t^11)
        """
        pf = self.__u.factor()
        u = LaurentPolynomial_univariate(self._parent, pf.unit(), self.__n)

        f = []
        for t in pf:
            d = LaurentPolynomial_univariate(self._parent, t[0], 0)
            if d.is_unit():
                u *= d**t[1]
            else:
                f.append( (d, t[1]) )

        return Factorization(f, unit=u)

    def residue(self):
        """
        Return the residue of ``self``.

        The residue is the coefficient of `t^-1`.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = 3*t^-2 - t^-1 + 3 + t^2
            sage: f.residue()
            -1
            sage: g = -2*t^-2 + 4 + 3*t
            sage: g.residue()
            0
            sage: f.residue().parent()
            Rational Field
        """
        return self[-1]

    def constant_coefficient(self):
        """
        Return the coefficient of the constant term of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = 3*t^-2 - t^-1 + 3 + t^2
            sage: f.constant_coefficient()
            3
            sage: g = -2*t^-2 + t^-1 + 3*t
            sage: g.constant_coefficient()
            0
        """
        return self[0]

cdef class LaurentPolynomial_mpair(LaurentPolynomial_generic):
    """
    Multivariate Laurent polynomials.
    """
    def __init__(self, parent, x, reduce=True):
        """
        Currently, one can only create LaurentPolynomials out of dictionaries
        and elements of the base ring.

        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: f = L({(-1,-1):1}); f
            w^-1*z^-1
            sage: f = L({(1,1):1}); f
            w*z
            sage: f =  L({(-1,-1):1, (1,3):4}); f
            4*w*z^3 + w^-1*z^-1
            sage: L(1/2)
            1/2
        """
        if isinstance(x, LaurentPolynomial_mpair):
            x = x.dict()
        elif isinstance(x, PolyDict):
            x = x.dict()
        if isinstance(x, dict):
            self._mon = ETuple({},int(parent.ngens()))
            for k in x.keys(): # ETuple-ize keys, set _mon
                if not isinstance(k, (tuple, ETuple)) or len(k) != parent.ngens():
                    self._mon = ETuple({}, int(parent.ngens()))
                    break
                if isinstance(k, tuple):
                    a = x[k]
                    del x[k]
                    k = ETuple(k)
                    x[k] = a
                self._mon = self._mon.emin(k) # point-wise min of _mon and k
            if len(self._mon.nonzero_positions()) != 0: # factor out _mon
                D = {}
                for k in x.keys():
                    D[k.esub(self._mon)] = x[k]
                x = D
        else: # since x should coerce into parent, _mon should be (0,...,0)
            self._mon = ETuple({}, int(parent.ngens()))
        self._poly = parent.polynomial_ring()(x)
        CommutativeAlgebraElement.__init__(self, parent)

    def __reduce__(self):
        """
        TESTS::

            sage: R = LaurentPolynomialRing(QQ,2,'x')
            sage: R.<x1,x2> = LaurentPolynomialRing(QQ)
            sage: loads(dumps(x1)) == x1 # indirect doctest
            True
            sage: z = x1/x2
            sage: loads(dumps(z)) == z   # not tested (bug)
            True
        """
        # one should also record the monomial self._mon
        return self._parent, (self._poly,)  # THIS IS WRONG !

    def __hash__(self):
        r"""
        TESTS:

        Test that the hash is non-constant (the hash does not need to be
        deterministic so we leave some slack for collisions)::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: len({hash(w^i*z^j) for i in [-2..2] for j in [-2..2]}) > 20
            True

        Check that :trac:`20490` is fixed::

            sage: R.<a,b> = LaurentPolynomialRing(ZZ)
            sage: p = a*~a
            sage: p._fraction_pair()
            (a, a)
            sage: p == R.one()
            True
            sage: hash(p) == hash(R.one())
            True
        """
        self._normalize()
        return hash(self._poly) ^ hash(self._mon)

    cdef _new_c(self):
        """
        Returns a new Laurent polynomial

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ) # indirect doctest
            sage: x*y
            x*y
        """
        cdef LaurentPolynomial_mpair ans
        ans = LaurentPolynomial_mpair.__new__(LaurentPolynomial_mpair)
        ans._parent = self._parent
        return ans

    def _normalize(self, i = None):
        """
        Removes the common monomials from self._poly and stores them in self._mon

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x*y + 2*y*x^2 + y # indirect doctest
            sage: f.factor() # Notice the y has been factored out.
            (y) * (2*x^2 + x + 1)
        """
        D = self._poly._mpoly_dict_recursive(self.parent().variable_names(), self.parent().base_ring())
        if i is None:
            e = None
            for k in D.keys():
                if e is None:
                    e = k
                else:
                    e = e.emin(k)
            if len(e.nonzero_positions()) > 0:
                self._poly = self._poly // self._poly.parent()({e: 1})
                self._mon = self._mon.eadd(e)
        else:
            e = None
            for k in D.keys():
                if e is None or k[i] < e:
                    e = k[i]
            if e > 0:
                self._poly = self._poly // self._poly.parent().gen(i)
                self._mon = self._mon.eadd_p(e, i)

    def _dict(self):
        """
        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1+3; a
            w^2*z^-1 + 3
            sage: d = a._dict()
            sage: keys = list(sorted(d.keys())); keys
            [(0, 0), (2, -1)]
            sage: d[keys[0]]
            3
            sage: d[keys[1]]
            1

        """
        D = self._poly._mpoly_dict_recursive(self.parent().variable_names(), self.parent().base_ring())
        if len(self._mon.nonzero_positions()) > 0:
            DD = {}
            for k in D.keys():
                DD[k.eadd(self._mon)] = D[k]
            return DD
        else:
            return D

    def _compute_polydict(self):
        """
        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1+3
            sage: a._compute_polydict()
        """
        self._prod = PolyDict(self._dict(), force_etuples = False)

    def is_unit(self):
        """
        Return ``True`` if ``self`` is a unit.

        The ground ring is assumed to be an integral domain.

        This means that the Laurent polynomial is a monomial
        with unit coefficient.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: (x*y/2).is_unit()
            True
            sage: (x + y).is_unit()
            False
            sage: (L.zero()).is_unit()
            False
            sage: (L.one()).is_unit()
            True

            sage: L.<x,y> = LaurentPolynomialRing(ZZ)
            sage: (2*x*y).is_unit()
            False
        """
        coeffs = self.coefficients()
        if len(coeffs) != 1:
            return False
        return coeffs[0].is_unit()

    def _repr_(self):
        """
        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x^2 + x*y/2 + 2*y^-1
            sage: f._repr_()
            'x^2 + 1/2*x*y + 2*y^-1'
        """
        if self._prod is None:
            self._compute_polydict()
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self._prod.poly_repr(self.parent().variable_names(),
                                    atomic_coefficients=atomic, cmpfn=cmpfn)

    def _latex_(self):
        """
        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1+3; a
            w^2*z^-1 + 3
            sage: latex(a)
            w^{2} z^{-1} + 3

        """
        if self._prod is None:
            self._compute_polydict()
        try:
            cmpfn = self.parent().term_order().compare_tuples
        except AttributeError:
            cmpfn = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self._prod.latex(self.parent().variable_names(),
                                atomic_coefficients=atomic, cmpfn=cmpfn)

    def __invert__(LaurentPolynomial_mpair self):
        """
        Return the inverse of ``self``.

        This treats monomials specially so they remain Laurent
        polynomials; the inverse of any other polynomial is an element
        of the rational function field.

        TESTS::

            sage: L.<x,y> = LaurentPolynomialRing(ZZ)
            sage: f = ~x
            sage: parent(f)
            Multivariate Laurent Polynomial Ring in x, y over Integer Ring
            sage: parent(f.coefficients()[0]) is parent(f).base_ring()
            True
            sage: g = ~(2*x)
            sage: parent(g)
            Multivariate Laurent Polynomial Ring in x, y over Rational Field
            sage: parent(g.coefficients()[0]) is parent(g).base_ring()
            True
        """
        cdef dict d = self.dict()
        cdef ETuple e
        if len(d) == 1:
            e, c = d.items()[0]
            e = e.emul(-1)
            P = self.parent()
            try:
                c = c.inverse_of_unit()
            except (AttributeError, ZeroDivisionError, ArithmeticError):
                c = ~c
                if c.parent() is not P.base_ring():
                    P = P.change_ring(c.parent())
            return P({e: c})
        return super(LaurentPolynomial_mpair, self).__invert__()

    def __pow__(LaurentPolynomial_mpair self, n, m):
        """
        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x + y
            sage: f^2
            x^2 + 2*x*y + y^2
            sage: f^(-1)
            1/(x + y)

        TESTS:

        Check that :trac:`2952` is fixed::

            sage: R.<q> = QQ[]
            sage: L.<x,y,z> = LaurentPolynomialRing(R)
            sage: f = (x+y+z^-1)^2
            sage: f.substitute(z=1)
            x^2 + 2*x*y + y^2 + 2*x + 2*y + 1
        """
        cdef LaurentPolynomial_mpair ans
        if n < 0:
            return ~(self ** -n)
        ans = self._new_c()
        ans._poly = self._poly ** n
        ans._mon = self._mon.emul(n)
        return ans

    def __getitem__(self, n):
        r"""
        Return the coefficient of `x^n = x_1^{n_1} \cdots x_k^{n_k}` where
        `n` is a tuple of length `k` and `k` is the number of variables.

        If the number of inputs is not equal to the number of variables, this
        raises a ``TypeError``.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3 + x*z; f
            -x^6 + x*z - 7*x^-2*y^3 + 5*x^-2*y + x^-3*y^2
            sage: f[6,0,0]
            -1
            sage: f[-2,3,0]
            -7
            sage: f[-1,4,2]
            0
            sage: f[1,0,1]
            1
            sage: f[6]
            Traceback (most recent call last):
            ...
            TypeError: Must have exactly 3 inputs
            sage: f[6,0]
            Traceback (most recent call last):
            ...
            TypeError: Must have exactly 3 inputs
            sage: f[6,0,0,0]
            Traceback (most recent call last):
            ...
            TypeError: Must have exactly 3 inputs
        """
        if isinstance(n, slice):
            raise TypeError("Multivariate Laurent polynomials are not iterable")
        if not isinstance(n, tuple) or len(n) != self.parent().ngens():
            raise TypeError("Must have exactly %s inputs"%self.parent().ngens())
        cdef ETuple t = ETuple(n)
        if self._prod is None:
            self._compute_polydict()
        if t not in self._prod.exponents():
            return self.parent().base_ring().zero()
        return self._prod[t]

    def __iter__(self):
        """
        Iterate through all terms by returning a list of the coefficient and
        the corresponding monomial.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: list(f) # indirect doctest
            [(-1, x^6), (1, x^-3*y^2), (5, x^-2*y), (-7, x^-2*y^3)]
        """
        if self._prod is None:
            self._compute_polydict()
        for c, exps in self._prod.list():
            prod = self.parent().one()
            for i in range(len(exps)):
                prod *= self.parent().gens()[i]**exps[i]
            yield (c, prod)

    def monomials(self):
        """
        Return the list of monomials in ``self``.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: f.monomials()
            [x^6, x^-3*y^2, x^-2*y, x^-2*y^3]
        """
        L = []
        if self._prod is None:
            self._compute_polydict()
        for c, exps in self._prod.list():
            prod = self.parent().one()
            for i in range(len(exps)):
                prod *= self.parent().gens()[i]**exps[i]
            L.append(prod)
        return L

    def monomial_coefficient(self, mon):
        """
        Return the coefficient in the base ring of the monomial ``mon`` in
        ``self``, where ``mon`` must have the same parent as ``self``.

        This function contrasts with the function :meth:`coefficient()`
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` - a monomial

        .. SEEALSO::

            For coefficients in a base ring of fewer variables, see
            :meth:`coefficient()`.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: f.monomial_coefficient(x^-2*y^3)
            -7
            sage: f.monomial_coefficient(x^2)
            0
        """
        if mon.parent() != self.parent():
            raise TypeError("Input must have the same parent")
        if self._prod is None:
            self._compute_polydict()
        if (<LaurentPolynomial_mpair>mon)._prod is None:
            mon._compute_polydict()
        return self.parent().base_ring()( self._prod.monomial_coefficient(
                        (<LaurentPolynomial_mpair>mon)._prod.dict()) )

    def constant_coefficient(self):
        """
        Return the constant coefficient of ``self``.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^2 + 5*x*y)*x^-3; f
            -x^6 - 7*x^-2*y^2 + 5*x^-2*y + x^-3*y^2
            sage: f.constant_coefficient()
            0
            sage: f = (x^3 + 2*x^-2*y+y^3)*y^-3; f
            x^3*y^-3 + 1 + 2*x^-2*y^-2
            sage: f.constant_coefficient()
            1
        """
        return self[(0,)*self.parent().ngens()]

    def coefficient(self, mon):
        r"""
        Return the coefficient of ``mon`` in ``self``, where ``mon`` must
        have the same parent as ``self``.

        The coefficient is defined as follows. If `f` is this polynomial, then
        the coefficient `c_m` is sum:

        .. MATH::

            c_m := \sum_T \frac{T}{m}

        where the sum is over terms `T` in `f` that are exactly divisible
        by `m`.

        A monomial `m(x,y)` 'exactly divides' `f(x,y)` if `m(x,y) | f(x,y)`
        and neither `x \cdot m(x,y)` nor `y \cdot m(x,y)` divides `f(x,y)`.

        INPUT:

        - ``mon`` -- a monomial

        OUTPUT:

        Element of the parent of ``self``.

        .. NOTE::

            To get the constant coefficient, call
            :meth:`constant_coefficient()`.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)

        The coefficient returned is an element of the parent of ``self``; in
        this case, ``P``.

        ::

            sage: f = 2 * x * y
            sage: c = f.coefficient(x*y); c
            2
            sage: c.parent()
            Multivariate Laurent Polynomial Ring in x, y over Rational Field

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^2 + 5*x*y)*x^-3; f
            -x^6 - 7*x^-2*y^2 + 5*x^-2*y + x^-3*y^2
            sage: f.coefficient(y)
            5*x^-2
            sage: f.coefficient(y^2)
            -7*x^-2 + x^-3
            sage: f.coefficient(x*y)
            0
            sage: f.coefficient(x^-2)
            -7*y^2 + 5*y
            sage: f.coefficient(x^-2*y^2)
            -7
            sage: f.coefficient(1)
            -x^6 - 7*x^-2*y^2 + 5*x^-2*y + x^-3*y^2
        """
        if mon.parent() is not self.parent():
            mon = self.parent()(mon)
        if self._prod is None:
            self._compute_polydict()
        if (<LaurentPolynomial_mpair>mon)._prod is None:
            mon._compute_polydict()
        return self.parent()(self._prod.coefficient((<LaurentPolynomial_mpair>mon).dict()))

    def coefficients(self):
        """
        Return the nonzero coefficients of this polynomial in a list.
        The returned list is decreasingly ordered by the term ordering
        of ``self.parent()``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ,order='degrevlex')
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.coefficients()
            [4, 3, 2, 1]
            sage: L.<x,y,z> = LaurentPolynomialRing(QQ,order='lex')
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.coefficients()
            [4, 1, 2, 3]
        """
        return self._poly.coefficients()

    def variables(self, sort=True):
        """
        Return a tuple of all variables occurring in self.

        INPUT:

        - ``sort`` -- specifies whether the indices shall be sorted

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.variables()
            (z, y, x)
            sage: f.variables(sort=False) #random
            (y, z, x)
        """
        d = self.dict();
        g = self.parent().gens()
        nvars = len(g)
        vars = []
        for k in d.keys():
            vars = union(vars,k.nonzero_positions())
            if len(vars) == nvars:
                break
        v = [ g[i] for i in vars]
        if sort:
            v.sort()
        return tuple(v)

    def dict(self):
        """
        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: list(sorted(f.dict().iteritems()))
            [((3, 1, 0), 3), ((4, 0, -2), 2), ((6, -7, 0), 1), ((7, 0, -1), 4)]
        """
        if self._prod is None:
            self._compute_polydict()
        return self._prod.dict()

    def _fraction_pair(self):
        """
        Return one representation of ``self`` as a pair
        ``(numerator, denominator)``.

        Here both the numerator and the denominator are polynomials.

        This is used for coercion into the fraction field.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f._fraction_pair()
            (4*x^7*y^7*z + 3*x^3*y^8*z^2 + 2*x^4*y^7 + x^6*z^2, y^7*z^2)
        """
        ring = self.parent().polynomial_ring()
        numer = self._poly
        denom = ring.one()
        var = ring.gens()
        for i, j in enumerate(self._mon):
            if j > 0:
                numer *= var[i] ** j
            else:
                denom *= var[i] ** (-j)
        return (numer, denom)

    cpdef _add_(self, _right):
        """
        Returns the Laurent polynomial self + right.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z
            sage: f + g
            x + y + z + y^-1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        cdef LaurentPolynomial_mpair right = <LaurentPolynomial_mpair>_right
        ans._mon, a, b = self._mon.combine_to_positives(right._mon)
        if len(a.nonzero_positions()) > 0:
            ans._poly = self._poly * self._poly.parent()({a: 1})
        else:
            ans._poly = self._poly
        if len(b.nonzero_positions()) > 0:
            ans._poly += right._poly * self._poly.parent()({b: 1})
        else:
            ans._poly += right._poly
        return ans

    cpdef _sub_(self, _right):
        """
        Returns the Laurent polynomial self - right.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z + x
            sage: f - g
            -y - z + y^-1

        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        cdef LaurentPolynomial_mpair right = <LaurentPolynomial_mpair>_right
        cdef ETuple a, b
        ans._mon, a, b = self._mon.combine_to_positives(right._mon)
        if len(a.nonzero_positions()) > 0:
            ans._poly = self._poly * self._poly.parent()({a: 1})
        else:
            ans._poly = self._poly
        if len(b.nonzero_positions()) > 0:
            ans._poly -= right._poly * self._poly.parent()({b: 1})
        else:
            ans._poly -= right._poly
        return ans

    cpdef _div_(self, rhs):
        """
        Return the division of ``self`` by ``rhs``.

        If the denominator is not a unit,
        the result will be given in the fraction field.

        EXAMPLES::

            sage: R.<s,q,t> = LaurentPolynomialRing(QQ)
            sage: 1/s
            s^-1
            sage: 1/(s*q)
            s^-1*q^-1
            sage: 1/(s+q)
            1/(s + q)
            sage: (1/(s+q)).parent()
            Fraction Field of Multivariate Polynomial Ring in s, q, t over Rational Field
            sage: (1/(s*q)).parent()
            Multivariate Laurent Polynomial Ring in s, q, t over Rational Field
            sage: (s+q)/(q^2*t^(-2))
            s*q^-2*t^2 + q^-1*t^2
        """
        cdef LaurentPolynomial_mpair right = <LaurentPolynomial_mpair> rhs
        if right.is_zero():
            raise ZeroDivisionError
        cdef dict d = right.dict()
        if len(d) == 1:
            return self * ~right
        else:
            return RingElement._div_(self, rhs)

    def is_monomial(self):
        """
        Return True if this element is a monomial.

        EXAMPLES::

            sage: k.<y,z> = LaurentPolynomialRing(QQ)
            sage: z.is_monomial()
            True
            sage: k(1).is_monomial()
            True
            sage: (z+1).is_monomial()
            False
            sage: (z^-2909).is_monomial()
            True
            sage: (38*z^-2909).is_monomial()
            False
        """

        d = self._poly.dict()
        return len(d) == 1 and 1 in d.values()

    cpdef _neg_(self):
        """
        Returns -self.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: -f
            -x - y^-1

        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon
        ans._poly = -self._poly
        return ans

    cpdef _lmul_(self, RingElement right):
        """
        Returns self * right where right is in self's base ring.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: f*(1/2)
            1/2*x + 1/2*y^-1

        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon
        ans._poly = self._poly * right
        return ans

    cpdef _rmul_(self, RingElement left):
        """
        Returns left*self where left is in self's base ring.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: (1/2)*f
            1/2*x + 1/2*y^-1

        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon
        ans._poly = left * self._poly
        return ans

    cpdef _mul_(self, right):
        """
        Return self*right.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z
            sage: f*g
            x*y + x*z + 1 + y^-1*z
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        ans._mon = self._mon.eadd((<LaurentPolynomial_mpair>right)._mon)
        ans._poly = self._poly * (<LaurentPolynomial_mpair>right)._poly
        return ans

    cpdef _floordiv_(self, right):
        """
        Perform division with remainder and return the quotient.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + y^-3
            sage: g = y + x
            sage: f // g
            x^5*y^-3 - x^4*y^-2 + x^3*y^-1

            sage: h = x + y^(-1)
            sage: f // h
            x^2 - x*y^-1 + y^-2
            sage: h * (f // h) == f
            True
            sage: f // 1
            x^3 + y^-3
            sage: 1 // f
            0

        TESTS:

        Check that :trac:`19357` is fixed::

            sage: x // y
            x*y^-1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        self._normalize()
        right._normalize()
        ans._mon = self._mon.esub((<LaurentPolynomial_mpair>right)._mon)
        ans._poly = self._poly // (<LaurentPolynomial_mpair>right)._poly
        return ans

    cpdef int _cmp_(self, right) except -2:
        """
        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + y^-1
            sage: g = y + z
            sage: f == f
            True
            sage: f == g
            False
            sage: f == 2
            False
        """
        if self._prod is None:
            self._compute_polydict()
        if (<LaurentPolynomial_mpair>right)._prod is None:
            right._compute_polydict()
        return cmp(self._prod, (<LaurentPolynomial_mpair>right)._prod)

    def exponents(self):
        """
        Returns a list of the exponents of self.

        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1+3; a
            w^2*z^-1 + 3
            sage: e = a.exponents()
            sage: e.sort(); e
            [(0, 0), (2, -1)]

        """
        return [a.eadd(self._mon) for a in self._poly.exponents()]

    def degree(self,x=None):
        """
        Returns the degree of x in self

        EXAMPLES::

            sage: R.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.degree(x)
            7
            sage: f.degree(y)
            1
            sage: f.degree(z)
            0
        """

        if not x:
            return self._poly.total_degree() + sum(self._mon)

        g = self.parent().gens()
        no_generator_found = True
        for i in range(len(g)):
            if g[i] is x:
                no_generator_found = False
                break
        if no_generator_found:
            raise TypeError("x must be a generator of parent")
        return self._poly.degree(self.parent().polynomial_ring().gens()[i]) + self._mon[i]



    def has_inverse_of(self, i):
        """
        INPUT:

        - ``i`` -- The index of a generator of ``self.parent()``

        OUTPUT:

        Returns True if self contains a monomial including the inverse of
        ``self.parent().gen(i)``, False otherwise.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.has_inverse_of(0)
            False
            sage: f.has_inverse_of(1)
            True
            sage: f.has_inverse_of(2)
            True
        """
        if (not isinstance(i, (int, Integer))) or (i < 0) or (i >= self.parent().ngens()):
            raise TypeError("argument is not the index of a generator")
        if self._mon[i] < 0:
            self._normalize(i)
            if self._mon[i] < 0:
                return True
            return False
        return False

    def has_any_inverse(self):
        """
        Returns True if self contains any monomials with a negative exponent, False otherwise.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.has_any_inverse()
            True
            sage: g = x^2 + y^2
            sage: g.has_any_inverse()
            False
        """
        for m in self._mon.nonzero_values(sort = False):
            if m < 0:
                return True
        return False

    def __call__(self, *x, **kwds):
        """

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + 2*y + 3*z
            sage: f(1,1,1)
            6
            sage: f = x^-1 + y + z
            sage: f(0,1,1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        TESTS::

            sage: f = x + 2*y + 3*z
            sage: f(2)
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match the number of generators in parent.
            sage: f(2,0)
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match the number of generators in parent.
            sage: f( (1,1,1) )
            6
        """
        n = self.parent().ngens()

        if len(kwds) > 0:
            f = self.subs(**kwds)
            if len(x) > 0:
                return f(*x)
            else:
                return f

        cdef int l = len(x)

        if l == 1 and (isinstance(x[0], tuple) or isinstance(x[0], list)):
            x = x[0]
            l = len(x)

        if l != n:
            raise TypeError("number of arguments does not match the number of generators in parent.")

        #Check to make sure that we aren't dividing by zero
        for m in range(n):
            if x[m] == 0:
                if self.has_inverse_of(m):
                    raise ZeroDivisionError

        ans = self._poly(*x)
        if ans != 0:
            for m in self._mon.nonzero_positions():
                ans *= x[m]**self._mon[m]

        return ans

    def subs(self, in_dict=None, **kwds):
        """
        Note that this is a very unsophisticated implementation.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = x + 2*y + 3*z
            sage: f.subs(x=1)
            2*y + 3*z + 1
            sage: f.subs(y=1)
            x + 3*z + 2
            sage: f.subs(z=1)
            x + 2*y + 3
            sage: f.subs(x=1,y=1,z=1)
            6

            sage: f = x^-1
            sage: f.subs(x=2)
            1/2
            sage: f.subs({x:2})
            1/2

            sage: f = x + 2*y + 3*z
            sage: f.subs({x:1,y:1,z:1})
            6
            sage: f.substitute(x=1,y=1,z=1)
            6

        TESTS::

            sage: f = x + 2*y + 3*z
            sage: f(q=10)
            x + 2*y + 3*z

        """
        if in_dict is not None and kwds:
            raise ValueError("you cannot specify both a dictionary and keyword arguments")

        g = self.parent().gens()
        repr_g = [repr(i) for i in g]
        vars = []

        if in_dict is None:
            for i in range(len(g)):
                if repr_g[i] in kwds:
                    vars.append(i)
        else:
            kwds = {}
            for i in range(len(g)):
                if g[i] in in_dict:
                    kwds[ repr(g[i]) ] = in_dict[ g[i] ]
                    vars.append(i)

        d = self._dict()
        out = 0
        for mon in d.keys():
            term = d[mon]
            for i in range(len(mon)):
                if i in vars:
                    term *= kwds[repr_g[i]]**mon[i]
                else:
                    term *= g[i]**mon[i]

            out += term

        return out

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: R.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + y/x
            sage: g = f._symbolic_(SR); g
            x^3 + y/x
            sage: g(x=2,y=2)
            9

            sage: g = SR(f)
            sage: g(x=2,y=2)
            9
        """
        d = dict([(repr(g), R.var(g)) for g in self.parent().gens()])
        return self.subs(**d)

    def derivative(self, *args):
        r"""
        The formal derivative of this Laurent polynomial, with respect
        to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        .. seealso::

           :meth:`_derivative`

        EXAMPLES::

            sage: R = LaurentPolynomialRing(ZZ,'x, y')
            sage: x, y = R.gens()
            sage: t = x**4*y+x*y+y+x**(-1)+y**(-3)
            sage: t.derivative(x, x)
            12*x^2*y + 2*x^-3
            sage: t.derivative(y, 2)
            12*y^-5
        """
        return multi_derivative(self, args)

    # add .diff(), .differentiate() as aliases for .derivative()
    diff = differentiate = derivative

    def _derivative(self, var=None):
        """
        Computes formal derivative of this Laurent polynomial with
        respect to the given variable.

        If var is among the generators of this ring, the derivative
        is with respect to the generator. Otherwise, _derivative(var) is called
        recursively for each coefficient of this polynomial.

        .. seealso:: :meth:`derivative`

        EXAMPLES::

            sage: R = LaurentPolynomialRing(ZZ,'x, y')
            sage: x, y = R.gens()
            sage: t = x**4*y+x*y+y+x**(-1)+y**(-3)
            sage: t._derivative(x)
            4*x^3*y + y - x^-2
            sage: t._derivative(y)
            x^4 + x + 1 - 3*y^-4

            sage: R = LaurentPolynomialRing(QQ['z'],'x')
            sage: z = R.base_ring().gen()
            sage: x = R.gen()
            sage: t = 33*z*x**4+x**(-1)
            sage: t._derivative(z)
            33*x^4
            sage: t._derivative(x)
            -x^-2 + 132*z*x^3
        """
        if var is None:
            raise ValueError("must specify which variable to differentiate "
                             "with respect to")
        P = self.parent()
        gens = list(P.gens())

        # check if var is one of the generators
        try:
            index = gens.index(var)
        except ValueError:
            # call _derivative() recursively on coefficients
            return P({m: c._derivative(var)
                      for (m, c) in self.dict().iteritems()})

        # compute formal derivative with respect to generator
        d = {}
        for m, c in self.dict().iteritems():
            if m[index] != 0:
                new_m = [u for u in m]
                new_m[index] += -1
                d[ETuple(new_m)] = m[index] * c
        return P(d)

    def is_univariate(self):
        """
        Return ``True`` if this is a univariate or constant Laurent polynomial,
        and ``False`` otherwise.

        EXAMPLES::

            sage: R.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = (x^3 + y^-3)*z
            sage: f.is_univariate()
            False
            sage: g = f(1,y,4)
            sage: g.is_univariate()
            True
            sage: R(1).is_univariate()
            True
        """
        return len(self.variables()) < 2

    def univariate_polynomial(self, R=None):
        """
        Returns a univariate polynomial associated to this
        multivariate polynomial.

        INPUT:

        - ``R`` - (default: ``None``) PolynomialRing

        If this polynomial is not in at most one variable, then a
        ``ValueError`` exception is raised.  The new polynomial is over
        the same base ring as the given ``LaurentPolynomial`` and in the
        variable ``x`` if no ring ``R`` is provided.

        EXAMPLES::

            sage: R.<x, y> = LaurentPolynomialRing(ZZ)
            sage: f = 3*x^2 - 2*y^-1 + 7*x^2*y^2 + 5
            sage: f.univariate_polynomial()
            Traceback (most recent call last):
            ...
            TypeError: polynomial must involve at most one variable
            sage: g = f(10,y); g
            700*y^2 + 305 - 2*y^-1
            sage: h = g.univariate_polynomial(); h
            -2*y^-1 + 305 + 700*y^2
            sage: h.parent()
            Univariate Laurent Polynomial Ring in y over Integer Ring
            sage: g.univariate_polynomial(LaurentPolynomialRing(QQ,'z'))
            -2*z^-1 + 305 + 700*z^2

        Here's an example with a constant multivariate polynomial::

            sage: g = R(1)
            sage: h = g.univariate_polynomial(); h
            1
            sage: h.parent()
            Univariate Laurent Polynomial Ring in x over Integer Ring
        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        v = self.variables()
        if len(v) > 1:
            raise TypeError("polynomial must involve at most one variable")
        elif len(v) == 1:
            x = v[0]
            i = self._parent.gens().index(x)
        else:
            x = 'x'
            i = 0

        #construct ring if none
        if R is None:
            R = LaurentPolynomialRing(self.base_ring(),x)

        return R(dict((m[i],c) for m,c in self.dict().items()))

    def factor(self):
        """
        Returns a Laurent monomial (the unit part of the factorization) and a factored multi-polynomial.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.factor()
            (x^3*y^-7*z^-2) * (4*x^4*y^7*z + 3*y^8*z^2 + 2*x*y^7 + x^3*z^2)
        """
        pf = self._poly.factor()
        u = self.parent(pf.unit().dict()) # self.parent won't currently take polynomials

        g = self.parent().gens()
        for i in self._mon.nonzero_positions():
            u *= g[i]**self._mon[i]

        f = []
        for t in pf:
            d = t[0].dict()
            if len(d) == 1:  # monomials are units
                u *= self.parent(d)**t[1]
            else:
                f.append( (self.parent(d),t[1]) )

        return Factorization(f, unit=u)
