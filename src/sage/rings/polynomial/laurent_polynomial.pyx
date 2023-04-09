r"""
Elements of Laurent polynomial rings
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer cimport Integer
from sage.categories.map cimport Map
from sage.structure.element import is_Element, coerce_binop, parent
from sage.structure.factorization import Factorization
from sage.misc.derivative import multi_derivative
from sage.rings.polynomial.polydict cimport monomial_exponent
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.structure.richcmp cimport richcmp, rich_to_bool
from sage.matrix.matrix0 cimport Matrix

cdef class LaurentPolynomial(CommutativeAlgebraElement):
    """
    Base class for Laurent polynomials.
    """
    cdef LaurentPolynomial _new_c(self):
        """
        Return a new Laurent polynomial.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ) # indirect doctest
            sage: x*y
            x*y
        """
        cdef type t = type(self)
        cdef LaurentPolynomial ans
        ans = t.__new__(t)
        ans._parent = self._parent
        return ans

    cpdef _add_(self, other):
        """
        Abstract addition method

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial
            sage: LaurentPolynomial._add_(x, x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    cpdef _mul_(self, other):
        """
        Abstract multiplication method

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial
            sage: LaurentPolynomial._mul_(x, x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    cpdef _floordiv_(self, other):
        """
        Abstract floor division method

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial
            sage: LaurentPolynomial._floordiv_(x, x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _integer_(self, ZZ):
        r"""
        Convert this Laurent polynomial to an integer.

        This is only possible if the Laurent polynomial is constant.

        OUTPUT:

        An integer.

        TESTS::

            sage: L.<a> = LaurentPolynomialRing(QQ)
            sage: L(42)._integer_(ZZ)
            42
            sage: a._integer_(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: L(2/3)._integer_(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: ZZ(L(42))
            42

        ::

            sage: L.<a, b> = LaurentPolynomialRing(QQ)
            sage: L(42)._integer_(ZZ)
            42
            sage: a._integer_(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: L(2/3)._integer_(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: ZZ(L(42))
            42
        """
        if not self.is_constant():
            raise ValueError('{} is not constant'.format(self))
        return ZZ(self.constant_coefficient())

    def _rational_(self):
        r"""
        Convert this Laurent polynomial to a rational.

        This is only possible if the Laurent polynomial is constant.

        OUTPUT:

        A rational.

        TESTS::

            sage: L.<a> = LaurentPolynomialRing(QQ)
            sage: L(42)._rational_()
            42
            sage: a._rational_()
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: QQ(L(2/3))
            2/3

        ::

            sage: L.<a, b> = LaurentPolynomialRing(QQ)
            sage: L(42)._rational_()
            42
            sage: a._rational_()
            Traceback (most recent call last):
            ...
            ValueError: a is not constant
            sage: QQ(L(2/3))
            2/3
        """
        if not self.is_constant():
            raise ValueError('{} is not constant'.format(self))
        from sage.rings.rational_field import QQ
        return QQ(self.constant_coefficient())

    def change_ring(self, R):
        """
        Return a copy of this Laurent polynomial, with coefficients in ``R``.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: a = x^2 + 3*x^3 + 5*x^-1
            sage: a.change_ring(GF(3))
            2*x^-1 + x^2

        Check that :trac:`22277` is fixed::

            sage: R.<x, y> = LaurentPolynomialRing(QQ)
            sage: a = 2*x^2 + 3*x^3 + 4*x^-1
            sage: a.change_ring(GF(3))
            -x^2 + x^-1
        """
        return self._parent.change_ring(R)(self)

    cpdef long number_of_terms(self) except -1:
        """
        Abstract method for number of terms

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial
            sage: LaurentPolynomial.number_of_terms(x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def hamming_weight(self):
        """
        Return the hamming weight of ``self``.

        The hamming weight is number of non-zero coefficients and
        also known as the weight or sparsity.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = x^3 - 1
            sage: f.hamming_weight()
            2
        """
        return self.number_of_terms()

    cpdef dict dict(self):
        """
        Abstract ``dict`` method.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial
            sage: LaurentPolynomial.dict(x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def map_coefficients(self, f, new_base_ring=None):
        """
        Apply ``f`` to the coefficients of ``self``.

        If ``f`` is a :class:`sage.categories.map.Map`, then the resulting
        polynomial will be defined over the codomain of ``f``. Otherwise, the
        resulting polynomial will be over the same ring as ``self``. Set
        ``new_base_ring`` to override this behavior.

        INPUT:

        - ``f`` -- a callable that will be applied to the coefficients of ``self``.

        - ``new_base_ring`` (optional) -- if given, the resulting polynomial
          will be defined over this ring.

        EXAMPLES::

            sage: k.<a> = GF(9)
            sage: R.<x> = LaurentPolynomialRing(k)
            sage: f = x*a + a
            sage: f.map_coefficients(lambda a : a + 1)
            (a + 1) + (a + 1)*x
            sage: R.<x,y> = LaurentPolynomialRing(k, 2)
            sage: f = x*a + 2*x^3*y*a + a
            sage: f.map_coefficients(lambda a : a + 1)
            (2*a + 1)*x^3*y + (a + 1)*x + a + 1

        Examples with different base ring::

            sage: R.<r> = GF(9); S.<s> = GF(81)
            sage: h = Hom(R,S)[0]; h
            Ring morphism:
              From: Finite Field in r of size 3^2
              To:   Finite Field in s of size 3^4
              Defn: r |--> 2*s^3 + 2*s^2 + 1
            sage: T.<X,Y> = LaurentPolynomialRing(R, 2)
            sage: f = r*X+Y
            sage: g = f.map_coefficients(h); g
            (2*s^3 + 2*s^2 + 1)*X + Y
            sage: g.parent()
            Multivariate Laurent Polynomial Ring in X, Y over Finite Field in s of size 3^4
            sage: h = lambda x: x.trace()
            sage: g = f.map_coefficients(h); g
            X - Y
            sage: g.parent()
            Multivariate Laurent Polynomial Ring in X, Y over Finite Field in r of size 3^2
            sage: g = f.map_coefficients(h, new_base_ring=GF(3)); g
            X - Y
            sage: g.parent()
            Multivariate Laurent Polynomial Ring in X, Y over Finite Field of size 3

        """
        R = self.parent()
        if new_base_ring is not None:
            R = R.change_ring(new_base_ring)
        elif isinstance(f, Map):
            R = R.change_ring(f.codomain())
        return R(dict([(k, f(v)) for (k, v) in self.dict().items()]))


cdef class LaurentPolynomial_univariate(LaurentPolynomial):
    r"""
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

        ::

            sage: R({})
            0
        """
        CommutativeAlgebraElement.__init__(self, parent)

        if isinstance(f, LaurentPolynomial_univariate):
            n += (<LaurentPolynomial_univariate>f).__n
            if (<LaurentPolynomial_univariate>f).__u._parent is parent._R:
                f = (<LaurentPolynomial_univariate>f).__u
            else:
                f = parent._R((<LaurentPolynomial_univariate>f).__u)
        elif (not isinstance(f, Polynomial)) or (parent is not f.parent()):
            if isinstance(f, dict):
                v = min(f) if f else 0
                f = {i-v: c for i,c in f.items()}
                n += v
            f = parent._R(f)

        # self is that t^n * u:
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

    def _polynomial_(self, R):
        r"""
        TESTS::

            sage: Lx = LaurentPolynomialRing(QQ, "x")
            sage: Px = PolynomialRing(QQ, "x")
            sage: Pxy = PolynomialRing(QQ, "x,y")
            sage: Paxb = PolynomialRing(QQ, "a,x,b")
            sage: Qx = PolynomialRing(ZZ, "x")
            sage: Rx = PolynomialRing(GF(2), "x")
            sage: p1 = Lx.gen()
            sage: p2 = Lx.zero()
            sage: p3 = Lx.one()
            sage: p4 = Lx.gen()**3 - 3
            sage: p5 = Lx.gen()**3 + 2*Lx.gen()**2
            sage: p6 = Lx.gen() >> 2

            sage: for P,x in [(Px, Px.gen()), (Qx, Qx.gen()), (Rx, Rx.gen()),
            ....:             (Pxy, Pxy.gen(0)), (Paxb, Paxb.gen(1))]:
            ....:     assert P(p1) == x and parent(P(p1)) is P
            ....:     assert P(p2) == P.zero() and parent(P(p2)) is P
            ....:     assert P(p3) == P.one() and parent(P(p3)) is P
            ....:     assert P(p4) == x**3 - 3 and parent(P(p4)) is P
            ....:     assert P(p5) == x**3 + 2*x**2 and parent(P(p5)) is P
            ....:     try: P(p6)
            ....:     except ValueError: pass
            ....:     else: raise RuntimeError

            sage: Pa = ZZ["a"]
            sage: Px = ZZ["x"]
            sage: Pax = ZZ["a,x"]
            sage: Pxa = ZZ["x,a"]
            sage: Pa_x = ZZ["a"]["x"]
            sage: Px_a = ZZ["x"]["a"]
            sage: Lax = LaurentPolynomialRing(Pa, "x")
            sage: Lxa = LaurentPolynomialRing(Px, "a")
            sage: for poly in ["2*a*x^2 - 5*x*a + 3", "a*x^2 - 3*a^3*x"]:
            ....:     assert Pax(Lax(poly)) == Pax(Lxa(poly)) == Pax(poly)
            ....:     assert Pxa(Lax(poly)) == Pxa(Lxa(poly)) == Pxa(poly)
            ....:     assert Pa_x(Lax(poly)) == Pa_x(poly)
            ....:     assert Px_a(Lxa(poly)) == Px_a(poly)
        """
        if self.__n < 0:
            raise ValueError("Laurent polynomial with negative valuation cannot be converted to polynomial")

        if is_PolynomialRing(R):
            return R(self.__u) << self.__n
        elif self.__n == 0:
            return R(self.__u)
        else:
            u = R(self.__u)
            x = R(self.__u._parent.gen())
            return x**self.__n * u

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

    def __bool__(self):
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

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of this element under the morphism defined by
        ``im_gens`` in ``codomain``, where elements of the
        base ring are mapped by ``base_map``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: H = Hom(R, QQ)
            sage: mor = H(2)
            sage: mor(t^2 + t^-2)
            17/4
            sage: 4 + 1/4
            17/4

        You can specify a map on the base ring::

            sage: Zx.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: cc = K.hom([-i])
            sage: R.<t> = LaurentPolynomialRing(K)
            sage: H = Hom(R, R)
            sage: phi = H([t^-2], base_map=cc)
            sage: phi(i*t)
            -i*t^-2
        """
        x = im_gens[0]
        u = self.__u
        if base_map is not None:
            u = u.map_coefficients(base_map)
        return codomain(u(x) * x**self.__n)

    cpdef __normalize(self):
        r"""
        A Laurent series is a pair `(u(t), n)`, where either `u = 0`
        (to some precision) or `u` is a unit. This pair corresponds to
        `t^n \cdot u(t)`.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: elt = t^2 + t^4 # indirect doctest
            sage: elt.polynomial_construction()
            (t^2 + 1, 2)

        Check that :trac:`21272` is fixed::

            sage: (t - t).polynomial_construction()
            (0, 0)
        """
        if self.__u[0]:
            return
        elif self.__u.is_zero():
            self.__n = 0
            return
        # we already caught the infinity and zero cases
        cdef long v = <long> self.__u.valuation()
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

        TESTS::

            sage: L.<lambda2> = LaurentPolynomialRing(QQ)
            sage: latex(L.an_element())
            \lambda_{2}
            sage: L.<y2> = LaurentPolynomialRing(QQ)
            sage: latex(L.an_element())
            y_{2}
        """
        from sage.misc.latex import latex

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
            x = latex(x)
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

        TESTS::

            sage: R = LaurentPolynomialRing(QQ, 't')

            sage: assert hash(R.zero()) == 0
            sage: assert hash(R.one()) == 1
            sage: assert hash(QQ['t'].gen()) == hash(R.gen())

            sage: for _ in range(20):
            ....:     p = QQ.random_element()
            ....:     assert hash(R(p)) == hash(p), "p = {}".format(p)

            sage: S.<t> = QQ[]
            sage: for _ in range(20):
            ....:     p = S.random_element()
            ....:     assert hash(R(p)) == hash(p), "p = {}".format(p)
            ....:     assert hash(R(t*p)) == hash(t*p), "p = {}".format(p)

        Check that :trac:`21272` is fixed::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: hash(R.zero()) == hash(t - t)
            True
        """
        # we reimplement below the hash of polynomials to handle negative
        # degrees
        cdef long result = 0
        cdef long result_mon
        cdef int i,j
        cdef long var_hash_name = hash(self.__u._parent._names[0])
        for i in range(self.__u.degree()+1):
            result_mon = hash(self.__u[i])
            if result_mon:
                j = i + self.__n
                if j > 0:
                    result_mon = (1000003 * result_mon) ^ var_hash_name
                    result_mon = (1000003 * result_mon) ^ j
                elif j < 0:
                    result_mon = (1000003 * result_mon) ^ var_hash_name
                    result_mon = (700005 * result_mon) ^ j
                result += result_mon
        return result

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

        Slicing can be used to truncate Laurent polynomials::

            sage: f[:3]
            -5*t^-10 + 1/3 + t + t^2

        Any other kind of slicing is an error, see :trac:`18940`::

            sage: f[-10:2]
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a start is not defined

            sage: f[-14:5:2]
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a step is not defined
        """
        cdef LaurentPolynomial_univariate ret
        if isinstance(i, slice):
            start, stop, step = i.start, i.stop, i.step
            if start is not None or step is not None:
                self.__u[start:stop:step]  # error out, see issue #18940
            stop = stop - self.__n if stop is not None else self.__u.degree() + 1
            f = self.__u[:stop]
            ret = <LaurentPolynomial_univariate> self._new_c()
            ret.__u = f
            ret.__n = self.__n
            ret.__normalize()
            return ret

        return self.__u[i - self.__n]

    cpdef long number_of_terms(self) except -1:
        """
        Return the number of non-zero coefficients of ``self``.

        Also called weight, hamming weight or sparsity.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = x^3 - 1
            sage: f.number_of_terms()
            2
            sage: R(0).number_of_terms()
            0
            sage: f = (x+1)^100
            sage: f.number_of_terms()
            101

        The method :meth:`hamming_weight` is an alias::

            sage: f.hamming_weight()
            101
        """
        return self.__u.number_of_terms()

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

        Since :trac:`24072` the symbolic ring does not accept positive
        characteristic::

            sage: R.<w> = LaurentPolynomialRing(GF(7))
            sage: SR(2*w^3 + 1)
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations
        """
        d = {repr(g): R.var(g) for g in self._parent.gens()}
        return self.subs(**d)

    cpdef dict dict(self):
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
        cdef dict d = self.__u.dict()
        return {k+self.__n: d[k] for k in d}

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
        return [i + self.__n for i in self.__u.exponents()]

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

    cpdef _unsafe_mutate(self, i, value):
        r"""
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
                coeffs = [value] + [R.zero() for _ in range(1,-j)] + self.__u.list()
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
        cdef LaurentPolynomial_univariate ret

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
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (f1 + f2)
        ret.__n = m
        ret.__normalize()
        return ret

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
        cdef LaurentPolynomial_univariate ret

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
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (f1 - f2)
        ret.__n = m
        ret.__normalize()
        return ret

    def degree(self):
        """
        Return the degree of ``self``.

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
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> -self.__u
        ret.__n = self.__n
        # No need to normalize
        return ret

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
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (self.__u * right.__u)
        ret.__n = self.__n + right.__n
        ret.__normalize()
        return ret

    cpdef _rmul_(self, Element c):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = 1/x^3 + x + x^2 + 3*x^4
            sage: 3 * f
            3*x^-3 + 3*x + 3*x^2 + 9*x^4
        """
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self.__u._rmul_(c)
        ret.__n = self.__n
        ret.__normalize()
        return ret

    cpdef _lmul_(self, Element c):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = 1/x^3 + x + x^2 + 3*x^4
            sage: f * 3
            3*x^-3 + 3*x + 3*x^2 + 9*x^4
        """
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self.__u._lmul_(c)
        ret.__n = self.__n
        ret.__normalize()
        return ret

    def is_monomial(self):
        r"""
        Return ``True`` if ``self`` is a monomial; that is, if ``self``
        is `x^n` for some integer `n`.

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

            sage: R.<x> = LaurentPolynomialRing(Zmod(6))
            sage: x^-2
            x^-2
            sage: (5*x^2)^-4
            x^-8
            sage: (5*x^-4)^-3
            5*x^12
        """
        cdef LaurentPolynomial_univariate self = _self
        cdef long right = r
        if right != r:
            raise ValueError("exponent must be an integer")
        try:
            return self._parent.element_class(self._parent, self.__u**right, self.__n*right)
        except TypeError as err:
            # we need to handle the special case of negative powers and a unit
            if not self.__u.is_constant() or not self.__u.leading_coefficient().is_unit():
                raise
            c = self._parent._R(self.__u.leading_coefficient() ** right)
            return self._parent.element_class(self._parent, c, self.__n*right)

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
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> (self.__u // right.__u)
        ret.__n = self.__n - right.__n
        ret.__normalize()
        return ret

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
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = self.__u
        ret.__n = self.__n + k
        # No need to normalize
        return ret

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
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = self.__u
        ret.__n = self.__n + k
        # No need to normalize
        return ret

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
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = self.__u
        ret.__n = self.__n - k
        # No need to normalize
        return ret

    cpdef _div_(self, rhs):
        """
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = x + x^2 + 3*x^4
            sage: g = 1/x^7 - x + x^2 - x^4
            sage: f / x
            1 + x + 3*x^3
            sage: f / g
            (-3*x^11 - x^9 - x^8)/(x^11 - x^9 + x^8 - 1)
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
        cdef LaurentPolynomial_univariate ret
        if self.__u.is_constant(): # this has a single term c*x^n
            ret = <LaurentPolynomial_univariate> self._new_c()
            if self.__u.is_unit():
                ret.__u = self.__u.inverse_of_unit()
                ret.__n = -self.__n
                ret.__normalize()
                return ret
            # Enlarge the ring so we can divide by the coefficient
            R = self._parent.base_ring().fraction_field()
            P = self._parent.change_ring(R)
            return P.element_class(P, ~R(self.__u), -self.__n)
        P = self._parent._R
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
        P = self._parent._R
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
        b = <LaurentPolynomial_univariate> self._parent(right)
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = self.__u.gcd(b.__u)
        ret.__n = min(self.__n, b.__n)
        ret.__normalize()
        return ret

    @coerce_binop
    def quo_rem(self, other):
        r"""
        Divide ``self`` by ``other`` and return a quotient ``q``
        and a remainder ``r`` such that ``self == q * other + r``.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: (t^-3 - t^3).quo_rem(t^-1 - t)
            (t^-2 + 1 + t^2, 0)
            sage: (t^-2 + 3 + t).quo_rem(t^-4)
            (t^2 + 3*t^4 + t^5, 0)

            sage: num = t^-2 + t
            sage: den = t^-2 + 1
            sage: q, r = num.quo_rem(den)
            sage: num == q * den + r
            True

        TESTS:

        Check that :trac:`34330` is fixed::

            sage: num = t^-2 + 3 + t
            sage: den = t^-4 + t
            sage: q, r = num.quo_rem(den); q, r
            (0, t^-2 + 3 + t)
            sage: num == q * den + r
            True

            sage: num = 2*t^-4 + t^-3 + t^-2 + 2*t + 2*t^2
            sage: q, r = num.quo_rem(den); q, r
            (2 + 2*t, -t^-3 + t^-2)
            sage: num == q * den + r
            True
        """
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate> other
        q, r = self.__u.quo_rem(right.__u)
        cdef LaurentPolynomial_univariate ql, qr
        ql = <LaurentPolynomial_univariate> self._new_c()
        ql.__u = <ModuleElement> q
        ql.__n = self.__n - right.__n
        ql.__normalize()
        qr = <LaurentPolynomial_univariate> self._new_c()
        qr.__u = <ModuleElement> r
        qr.__n = self.__n
        qr.__normalize()
        return ql, qr

    cpdef _richcmp_(self, right_r, int op):
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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate> right_r

        zero = self._parent.base_ring().zero()

        if not self and not right:
            return rich_to_bool(op, 0)

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

        return richcmp(x, y, op)

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
            return self._parent.zero()
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self.__u.truncate(n - self.__n)
        ret.__n = self.__n
        ret.__normalize()
        return ret

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
        Return the tuple of variables occurring in this Laurent polynomial.

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
        Return whether this Laurent polynomial is constant.

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
            sage: R(0).is_constant()
            True
            sage: R(42).is_constant()
            True
            sage: x.is_constant()
            False
            sage: (1/x).is_constant()
            False
        """
        return self.__n == 0 and self.__u.is_constant()


    def is_square(self, root=False):
        r"""
        Return whether this Laurent polynomial is a square.

        If ``root`` is set to ``True`` then return a pair made of the
        boolean answer together with ``None`` or a square root.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)

            sage: R.one().is_square()
            True
            sage: R(2).is_square()
            False

            sage: t.is_square()
            False
            sage: (t**-2).is_square()
            True

        Usage of the ``root`` option::

            sage: p = (1 + t^-1 - 2*t^3)
            sage: p.is_square(root=True)
            (False, None)
            sage: (p**2).is_square(root=True)
            (True, -t^-1 - 1 + 2*t^3)

        The answer is dependent of the base ring::

            sage: S.<u> = LaurentPolynomialRing(QQbar)
            sage: (2 + 4*t + 2*t^2).is_square()
            False
            sage: (2 + 4*u + 2*u^2).is_square()
            True

        TESTS::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: (t - t).is_square(True)
            (True, 0)

            sage: for _ in range(10):
            ....:     p = t ** randint(-15,15) * sum(QQ.random_element() * t**n for n in range(randint(5,10)))
            ....:     ans, r = (p**2).is_square(root=True)
            ....:     assert ans
            ....:     assert r*r == p*p
        """
        cdef LaurentPolynomial_univariate sqrt
        if self.__n % 2:
            return (False, None) if root else False
        elif root:
            ans, r = self.__u.is_square(True)
            if ans:
                sqrt = self._new_c()
                sqrt.__u = r
                sqrt.__n = self.__n // 2
                return (True, sqrt)
            else:
                return (False, None)
        else:
            return self.__u.is_square(False)

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
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = copy(self.__u)
        ret.__n = self.__n
        # No need to normalize
        return ret

    def derivative(self, *args):
        """
        The formal derivative of this Laurent polynomial, with respect
        to variables supplied in args.

        Multiple variables and iteration counts may be supplied. See
        documentation for the global :func:`derivative` function for more
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

        Check that :trac:`28187` is fixed::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: p = 1/x + 1 + x
            sage: x,y = var("x, y")
            sage: p._derivative(x)
            -x^-2 + 1
            sage: p._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        cdef LaurentPolynomial_univariate ret
        if var is not None and var != self._parent.gen():
            try:
                # call _derivative() recursively on coefficients
                u = [coeff._derivative(var) for coeff in self.__u.list(copy=False)]
                ret = <LaurentPolynomial_univariate> self._new_c()
                ret.__u = <ModuleElement> self._parent._R(u)
                ret.__n = self.__n
                ret.__normalize()
                return ret
            except AttributeError:
                raise ValueError('cannot differentiate with respect to {}'.format(var))

        # compute formal derivative with respect to generator
        if self.is_zero():
            return self  # this is already 0
        cdef long m, n = self.__n
        cdef list a = self.__u.list(copy=True)
        for m in range(len(a)):
            a[m] *= n + m
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> self._parent._R(a)
        ret.__n = self.__n - 1
        ret.__normalize()
        return ret

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
        cdef LaurentPolynomial_univariate ret
        if self[-1] != 0:
            raise ArithmeticError("the integral of is not a Laurent polynomial,"
                                  " since t^-1 has nonzero coefficient")

        cdef list a = self.__u.list(copy=False)
        if n < 0:
            v = [a[i]/(n+i+1) for i in range(min(-1-n,len(a)))] + [0]
        else:
            v = []
        v += [a[i]/(n+i+1) for i in range(max(-n,0), len(a))]
        try:
            u = self._parent._R(v)
        except TypeError:
            raise ArithmeticError("coefficients of integral cannot be coerced into the base ring")
        ret = <LaurentPolynomial_univariate> self._new_c()
        ret.__u = <ModuleElement> u
        ret.__n = n + 1
        ret.__normalize()
        return ret

    def __call__(self, *x, **kwds):
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
            sage: f(t=-1)
            2
            sage: f(x=-1)
            t^-2 + t^2
            sage: f()
            t^-2 + t^2
            sage: f(1,2)
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match number of
             variables in parent
        """
        if kwds:
            f = self.subs(**kwds)
            if x: # If there are non-keyword arguments
                return f(*x)
            else:
                return f

        if not x:
            return self
        if len(x) != 1:
            raise TypeError("number of arguments does not match number"
                            " of variables in parent")
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
        cdef LaurentPolynomial_univariate u, d
        pf = self.__u.factor()
        u = <LaurentPolynomial_univariate> self._new_c()
        u.__u = pf.unit()
        u.__n = self.__n
        u.__normalize()

        f = []
        for t in pf:
            d = <LaurentPolynomial_univariate> self._new_c()
            d.__u = t[0]
            d.__n = 0
            d.__normalize()
            if d.is_unit():
                u *= d ** t[1]
            else:
                f.append((d, t[1]))

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
        return self.__u[-1 - self.__n]

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
        return self.__u[-self.__n]


cdef class LaurentPolynomial_mpair(LaurentPolynomial):
    """
    Multivariate Laurent polynomials.
    """
    def __init__(self, parent, x, mon=None, reduce=True):
        """
        Currently, one can only create LaurentPolynomials out of dictionaries
        and elements of the base ring.

        INPUT:

        - ``parent`` -- a SageMath parent

        - ``x`` -- an element or dictionary or anything the underlying
          polynomial ring accepts

        - ``mon`` -- (default: ``None``) a tuple specifying the shift
          in the exponents

        - ``reduce`` -- (default: ``True``) a boolean

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

        TESTS:

        Check that :trac:`19538` is fixed::

            sage: R = LaurentPolynomialRing(QQ,'x2,x0')
            sage: S = LaurentPolynomialRing(QQ,'x',3)
            sage: f = S.coerce_map_from(R)
            sage: f(R.gen(0) + R.gen(1)^2)
            x0^2 + x2
            sage: _.parent()
            Multivariate Laurent Polynomial Ring in x0, x1, x2 over Rational Field

        ::

            sage: from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial_mpair
            sage: LaurentPolynomial_mpair(L, {(1,2): 1/42}, mon=(-3, -3))
            1/42*w^-2*z^-1

        :trac:`22398`::

            sage: LQ = LaurentPolynomialRing(QQ, 'x0, x1, x2, y0, y1, y2, y3, y4, y5')
            sage: LZ = LaurentPolynomialRing(ZZ, 'x0, x1, x2, y0, y1, y2, y3, y4, y5')
            sage: LQ.inject_variables()
            Defining x0, x1, x2, y0, y1, y2, y3, y4, y5
            sage: x2^-1*y0*y1*y2*y3*y4*y5 + x1^-1*x2^-1*y0*y1*y3*y4 + x0^-1 in LZ
            True
            sage: x2^-1*y0*y1*y2*y3*y4*y5 + x1^-1*x2^-1*y0*y1*y3*y4 + x0^-1*x1^-1*y0*y3 + x0^-1 in LZ
            True

        Check that input is not modified::

            sage: LQ.<x,y> = LaurentPolynomialRing(QQ)
            sage: D = {(-1, 1): 1}
            sage: k = tuple(D)[0]
            sage: v = D[k]
            sage: type(k), type(v)
            (<... 'tuple'>, <class 'sage.rings.integer.Integer'>)
            sage: LQ(D)
            x^-1*y
            sage: tuple(D)[0] is k
            True
            sage: D[k] is v
            True
        """
        if isinstance(x, PolyDict):
            x = x.dict()
        if mon is not None:
            if isinstance(mon, ETuple):
                self._mon = mon
            else:
                self._mon = ETuple(mon)
        else:
            if isinstance(x, dict):
                self._mon = ETuple({}, int(parent.ngens()))
                D = {}
                for k, x_k in x.iteritems():  # ETuple-ize keys, set _mon
                    if not isinstance(k, (tuple, ETuple)) or len(k) != parent.ngens():
                        self._mon = ETuple({}, int(parent.ngens()))
                        break
                    if isinstance(k, tuple):
                        k = ETuple(k)
                    D[k] = x_k
                    self._mon = self._mon.emin(k) # point-wise min of _mon and k
                else:
                    x = D
                if not self._mon.is_constant(): # factor out _mon
                    x = {k.esub(self._mon): x_k for k, x_k in x.iteritems()}
            elif (isinstance(x, LaurentPolynomial_mpair) and
                  parent.variable_names() == x.parent().variable_names()):
                self._mon = (<LaurentPolynomial_mpair>x)._mon
                x = (<LaurentPolynomial_mpair>x)._poly
            else: # since x should coerce into parent, _mon should be (0,...,0)
                self._mon = ETuple({}, int(parent.ngens()))
        self._poly = parent._R(x)
        CommutativeAlgebraElement.__init__(self, parent)

    def __reduce__(self):
        """
        TESTS::

            sage: R = LaurentPolynomialRing(QQ,2,'x')
            sage: R.<x1,x2> = LaurentPolynomialRing(QQ)
            sage: loads(dumps(x1)) == x1 # indirect doctest
            True
            sage: z = x1/x2
            sage: loads(dumps(z)) == z
            True
        """
        return self._parent, (self._poly, self._mon)

    def __hash__(self):
        r"""
        TESTS:

        Test that the hash is non-constant (see also :trac:`27914`)::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: len({hash(w^i*z^j) for i in [-2..2] for j in [-2..2]})
            25

        Check that :trac:`20490` is fixed::

            sage: R.<a,b> = LaurentPolynomialRing(ZZ)
            sage: p = a*~a
            sage: p._fraction_pair()
            (a, a)
            sage: p == R.one()
            True
            sage: hash(p)
            1

        Check that :trac:`23864` is fixed (compatibility with integers, rationals
        and polynomial rings)::

            sage: L = LaurentPolynomialRing(QQ, 'x0,x1,x2')
            sage: hash(L.zero())
            0
            sage: hash(L.one())
            1
            sage: hash(-L.one())
            -2
            sage: hash(L(1/2)) == hash(1/2)
            True

            sage: R = PolynomialRing(QQ, 'x0,x1,x2')
            sage: x0,x1,x2 = R.gens()
            sage: hash(x0) == hash(L(x0))
            True
            sage: hash(1 - 7*x0 + x1*x2) == hash(L(1 - 7*x0 + x1*x2))
            True

        Check that :trac:`27914` is fixed::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: Lw = LaurentPolynomialRing(QQ, 'w')
            sage: Lz = LaurentPolynomialRing(QQ, 'z')
            sage: all(hash(w^k) == hash(Lw(w^k))
            ....:     and hash(z^k) == hash(Lz(z^k)) for k in (-5..5))
            True
            sage: p = w^-1 + 2 + w
            sage: hash(p) == hash(Lw(p))
            True
        """
        # we reimplement the hash from multipolynomial to handle negative exponents
        # (see multi_polynomial.pyx)
        cdef long result = 0
        cdef long exponent
        cdef list var_name_hash = [hash(v) for v in self._parent.variable_names()]
        cdef int p
        cdef int n = len(var_name_hash)
        cdef long c_hash
        for m, c in self._poly.iterator_exp_coeff():
            c_hash = hash(c)
            if c_hash != 0:
                for p in range(n):
                    exponent = m[p] + self._mon[p]
                    if exponent > 0:
                        c_hash = (1000003 * c_hash) ^ var_name_hash[p]
                        c_hash = (1000003 * c_hash) ^ exponent
                    elif exponent < 0:
                        c_hash = (1000003 * c_hash) ^ var_name_hash[p]
                        c_hash = (700005 * c_hash) ^ exponent
                result += c_hash

        return result

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the morphism defined by
        ``im_gens`` in ``codomain``.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(ZZ)
            sage: M.<u,v> = LaurentPolynomialRing(ZZ)
            sage: phi = L.hom([u,v])
            sage: phi(x^2*~y -5*y**3)            # indirect doctest
            -5*v^3 + u^2*v^-1

        TESTS:

        check compatibility with  :trac:`26105`::

            sage: F.<t> = GF(4)
            sage: LF.<a,b> = LaurentPolynomialRing(F)
            sage: rho = LF.hom([b,a], base_map=F.frobenius_endomorphism())
            sage: s = t*~a + b +~t*(b**-3)*a**2; rs = rho(s); rs
            a + (t + 1)*b^-1 + t*a^-3*b^2
            sage: s == rho(rs)
            True
        """
        p = self._poly
        m = self._mon
        if base_map is not None:
            p = p.map_coefficients(base_map)
        from sage.misc.misc_c import prod
        return codomain(p(im_gens) * prod(ig**m[im_gens.index(ig)] for ig in im_gens))

    cdef _normalize(self, i=None):
        r"""
        Remove the common monomials from ``self._poly`` and store
        them in ``self._mon``.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x*y + 2*y*x^2 + y  # indirect doctest
            sage: f.factor() # Notice the y has been factored out.
            (y) * (2*x^2 + x + 1)

        Check that :trac:`23864` has been fixed::

            sage: hash(L.zero())
            0
        """
        if not self._poly:
            self._mon = ETuple({}, int(self._parent.ngens()))
            return

        #cdef dict D = <dict> self._poly._mpoly_dict_recursive(
        #                                <tuple> self._parent.variable_names(),
        #                                self._parent.base_ring()
        #                                )
        cdef dict D = <dict> self._poly.dict()

        cdef ETuple e
        if i is None:
            e = None
            for k in D:
                if e is None:
                    e = <ETuple> k
                else:
                    e = e.emin(k)
            if not e.is_constant():
                self._poly = <ModuleElement> (self._poly // self._poly._parent({e: 1}))
                self._mon = self._mon.eadd(e)
        else:
            e = None
            for k in D:
                if e is None or k[i] < e:
                    e = <ETuple> k[i]
            if e > 0:
                self._poly = <ModuleElement> (self._poly // self._poly._parent.gen(i))
                self._mon = self._mon.eadd_p(e, i)

    cdef _compute_polydict(self):
        """
        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1 +3
            sage: a.dict()  # indirect doctest
            {(0, 0): 3, (2, -1): 1}
        """
        #cdef dict D = <dict> self._poly._mpoly_dict_recursive(self._parent.variable_names(),
        #                                                      self._parent.base_ring())
        cdef dict D = <dict> self._poly.dict()
        cdef dict DD
        if self._mon.is_constant():
            self._prod = PolyDict(D)
            return
        DD = {}
        for k in D:
            DD[k.eadd(self._mon)] = D[k]
        self._prod = PolyDict(DD)

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
            key = self.parent().term_order().sortkey
        except AttributeError:
            key = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self._prod.poly_repr(self.parent().variable_names(),
                                    atomic_coefficients=atomic, sortkey=key)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1+3; a
            w^2*z^-1 + 3
            sage: latex(a)
            w^{2} z^{-1} + 3

        TESTS::

            sage: L.<lambda2, y2> = LaurentPolynomialRing(QQ)
            sage: latex(1/lambda2 + y2^(-3))
            \lambda_{2}^{-1} + y_{2}^{-3}
        """
        if self._prod is None:
            self._compute_polydict()
        try:
            key = self.parent().term_order().sortkey
        except AttributeError:
            key = None
        atomic = self.parent().base_ring()._repr_option('element_is_atomic')
        return self._prod.latex(self.parent().latex_variable_names(),
                                atomic_coefficients=atomic, sortkey=key)

    cpdef long number_of_terms(self) except -1:
        """
        Return the number of non-zero coefficients of ``self``.

        Also called weight, hamming weight or sparsity.

        EXAMPLES::

            sage: R.<x, y> = LaurentPolynomialRing(ZZ)
            sage: f = x^3 - y
            sage: f.number_of_terms()
            2
            sage: R(0).number_of_terms()
            0
            sage: f = (x+1/y)^100
            sage: f.number_of_terms()
            101

        The method :meth:`hamming_weight` is an alias::

            sage: f.hamming_weight()
            101
        """
        return self._poly.number_of_terms()

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
        cdef ETuple e
        if self._poly.is_term():
            (e, c), = self.dict().items()
            e = e.emul(-1)
            P = self._parent
            try:
                c = c.inverse_of_unit()
            except (AttributeError, ZeroDivisionError, ArithmeticError):
                c = ~c
                if c.parent() is not P.base_ring():
                    P = P.change_ring(c.parent())
            return P({e: c})
        return super().__invert__()

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
            TypeError: must have exactly 3 inputs
            sage: f[6,0]
            Traceback (most recent call last):
            ...
            TypeError: must have exactly 3 inputs
            sage: f[6,0,0,0]
            Traceback (most recent call last):
            ...
            TypeError: must have exactly 3 inputs
        """
        if isinstance(n, slice):
            raise TypeError("multivariate Laurent polynomials are not iterable")
        if not isinstance(n, tuple) or len(n) != self._parent.ngens():
            raise TypeError("must have exactly %s inputs" %
                            self.parent().ngens())
        cdef ETuple t = ETuple(n)
        if self._prod is None:
            self._compute_polydict()
        try:
            return self._prod[t]
        except KeyError:
            return self._parent.base_ring().zero()

    def __iter__(self):
        """
        Iterate through all terms by returning a list of the coefficient and
        the corresponding monomial.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: sorted(f) # indirect doctest
            [(-7, x^-2*y^3), (-1, x^6), (1, x^-3*y^2), (5, x^-2*y)]
        """
        P = self._parent
        one = P._R.one()
        if self._mon.is_constant():
            for exp, coeff in self._poly.iterator_exp_coeff():
                yield (coeff, P.element_class(P, one, exp))
        else:
            for exp, coeff in self._poly.iterator_exp_coeff():
                yield (coeff, P.element_class(P, one, exp.eadd(self._mon)))

    def iterator_exp_coeff(self):
        """
        Iterate over ``self`` as pairs of (ETuple, coefficient).

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: list(f.iterator_exp_coeff())
            [((6, 0), -1), ((-2, 3), -7), ((-2, 1), 5), ((-3, 2), 1)]
        """
        for exp, coeff in self._poly.iterator_exp_coeff():
            yield (exp.eadd(self._mon), coeff)

    def monomials(self):
        """
        Return the list of monomials in ``self``.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = (y^2 - x^9 - 7*x*y^3 + 5*x*y)*x^-3
            sage: sorted(f.monomials())
            [x^-3*y^2, x^-2*y, x^-2*y^3, x^6]
        """
        return [mon for coeff, mon in self]

    def monomial_coefficient(self, mon):
        """
        Return the coefficient in the base ring of the monomial ``mon`` in
        ``self``, where ``mon`` must have the same parent as ``self``.

        This function contrasts with the function :meth:`coefficient()`
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` -- a monomial

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

        TESTS::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = y^2 * x^-2
            sage: f.monomial_coefficient(x + y)
            Traceback (most recent call last):
            ...
            ValueError: not a monomial
        """
        if parent(mon) != self._parent:
            raise TypeError("input must have the same parent")
        cdef LaurentPolynomial_mpair m = <LaurentPolynomial_mpair> mon
        if m._prod is None:
            m._compute_polydict()
        if self._prod is None:
            self._compute_polydict()
        exp = monomial_exponent(m._prod)
        zero = self._parent.base_ring().zero()
        return self._prod.get(exp, zero)

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
        return self[(0,)*self._parent.ngens()]

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
        this case, ``P``. ::

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
        if mon.parent() is not self._parent:
            mon = self._parent(mon)
        cdef LaurentPolynomial_mpair m = <LaurentPolynomial_mpair> mon
        if self._prod is None:
            self._compute_polydict()
        if m._prod is None:
            m._compute_polydict()
        return self._parent(self._prod.coefficient(m.dict()))

    def coefficients(self):
        """
        Return the nonzero coefficients of ``self`` in a list.

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
        Return a tuple of all variables occurring in ``self``.

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
        cdef dict d = self.dict()
        cdef tuple g = self._parent.gens()
        cdef Py_ssize_t nvars = len(g)
        cdef set vars = set()
        for k in d:
            vars.update(k.nonzero_positions())
            if len(vars) == nvars:
                break
        cdef list v = [g[i] for i in vars]
        if sort:
            v.sort()
        return tuple(v)

    cpdef dict dict(self):
        """
        Return ``self`` represented as a ``dict``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: sorted(f.dict().items())
            [((3, 1, 0), 3), ((4, 0, -2), 2), ((6, -7, 0), 1), ((7, 0, -1), 4)]
        """
        if self._prod is None:
            self._compute_polydict()
        return <dict> self._prod.dict()

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
        ring = self._parent._R
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
        Return the Laurent polynomial ``self + right``.

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
        if not a.is_constant():
            ans._poly = self._poly * self._poly._parent({a: 1})
        else:
            ans._poly = self._poly
        if not b.is_constant():
            ans._poly += right._poly * self._poly._parent({b: 1})
        else:
            ans._poly += right._poly
        return ans

    cpdef _sub_(self, _right):
        """
        Return the Laurent polynomial ``self - right``.

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
        if not a.is_constant():
            ans._poly = self._poly * self._poly._parent({a: 1})
        else:
            ans._poly = self._poly
        if not b.is_constant():
            ans._poly -= right._poly * self._poly._parent({b: 1})
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
        if right._poly.is_term():
            return self * ~right
        else:
            return RingElement._div_(self, rhs)

    def is_monomial(self):
        """
        Return ``True`` if ``self`` is a monomial.

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
        return self._poly.is_monomial()

    cpdef _neg_(self):
        """
        Return ``-self``.

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

    cpdef _lmul_(self, Element right):
        """
        Return ``self * right`` where ``right`` is in ``self``'s base ring.

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

    cpdef _rmul_(self, Element left):
        """
        Return ``left * self`` where ``left`` is in ``self``'s base ring.

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
        Return ``self * right``.

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

        Check that :trac:`21999` is fixed::

            sage: L.<a,b> = LaurentPolynomialRing(QQbar)
            sage: (a+a*b) // a
            b + 1
        """
        cdef LaurentPolynomial_mpair ans = self._new_c()
        cdef LaurentPolynomial_mpair rightl = <LaurentPolynomial_mpair> right
        self._normalize()
        rightl._normalize()
        ans._mon = self._mon.esub(rightl._mon)
        ans._poly = self._poly // rightl._poly
        return ans

    @coerce_binop
    def quo_rem(self, right):
        """
        Divide this Laurent polynomial by ``right`` and return a quotient and
        a remainder.

        INPUT:

        - ``right`` -- a Laurent polynomial

        OUTPUT:

        A pair of Laurent polynomials.

        EXAMPLES::

            sage: R.<s, t> = LaurentPolynomialRing(QQ)
            sage: (s^2-t^2).quo_rem(s-t)
            (s + t, 0)
            sage: (s^-2-t^2).quo_rem(s-t)
            (s + t, -s^2 + s^-2)
            sage: (s^-2-t^2).quo_rem(s^-1-t)
            (t + s^-1, 0)

        TESTS:

        Verify that :trac:`31257` is fixed::

            sage: R.<x,y> = LaurentPolynomialRing(QQ)
            sage: q, r = (1/x).quo_rem(y)
            sage: q, r
            (x^-1*y^-1, 0)
            sage: q*y + r == 1/x
            True
            sage: q,r = (x^-2 - y^2).quo_rem(x - y)
            sage: q*(x - y) + r == x^-2 - y^2
            True
        """
        # make copies of self and right so that the input can be normalized
        # without affecting the objects that were passed to the method
        cdef LaurentPolynomial_mpair selfl = self._new_c()
        selfl._poly = self._poly
        selfl._mon = self._mon
        cdef LaurentPolynomial_mpair rightl = self._new_c()
        rightl._poly = (<LaurentPolynomial_mpair> right)._poly
        rightl._mon = (<LaurentPolynomial_mpair> right)._mon

        selfl._normalize()
        rightl._normalize()
        q, r = selfl._poly.quo_rem(rightl._poly)
        ql = LaurentPolynomial_mpair(self._parent, q,
                                     mon=selfl._mon.esub(rightl._mon))
        rl = LaurentPolynomial_mpair(self._parent, r,
                                     mon=selfl._mon)
        ql._normalize()
        rl._normalize()
        return (ql, rl)

    cpdef _richcmp_(self, right, int op):
        """
        Compare two polynomials in a `LaurentPolynomialRing` based on the term
        order from the parent ring.  If the parent ring does not specify a term
        order then only comparison by equality is supported.

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
        if (<LaurentPolynomial_mpair> right)._prod is None:
            (<LaurentPolynomial_mpair> right)._compute_polydict()

        try:
            sortkey = self._parent.term_order().sortkey
        except AttributeError:
            sortkey = None

        return self._prod.rich_compare((<LaurentPolynomial_mpair>right)._prod,
                                       op, sortkey)

    def exponents(self):
        """
        Return a list of the exponents of ``self``.

        EXAMPLES::

            sage: L.<w,z> = LaurentPolynomialRing(QQ)
            sage: a = w^2*z^-1+3; a
            w^2*z^-1 + 3
            sage: e = a.exponents()
            sage: e.sort(); e
            [(0, 0), (2, -1)]

        """
        return [a.eadd(self._mon) for a in self._poly.exponents()]

    def degree(self, x=None):
        """
        Return the degree of ``x`` in ``self``.

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

        cdef tuple g = <tuple> self._parent.gens()
        cdef Py_ssize_t i
        cdef bint no_generator_found = True
        for i in range(len(g)):
            if g[i] is x:
                no_generator_found = False
                break
        if no_generator_found:
            raise TypeError("x must be a generator of parent")
        return self._poly.degree(self._parent._R.gens()[i]) + self._mon[i]

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
        if (not isinstance(i, (int, Integer))) or (i < 0) or (i >= self._parent.ngens()):
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
        for m in self._mon.nonzero_values(sort=False):
            if m < 0:
                return True
        return False

    def __call__(self, *x, **kwds):
        """
        Compute value of ``self`` at ``x``.

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
            TypeError: number of arguments does not match the number of generators in parent
            sage: f(2,0)
            Traceback (most recent call last):
            ...
            TypeError: number of arguments does not match the number of generators in parent
            sage: f( (1,1,1) )
            6
        """
        if kwds:
            f = self.subs(**kwds)
            if x: # More than 1 non-keyword argument
                return f(*x)
            else:
                return f

        cdef int l = len(x)

        if l == 1 and isinstance(x[0], (tuple, list)):
            x = x[0]
            l = len(x)

        if l != self._parent.ngens():
            raise TypeError("number of arguments does not match the number"
                            " of generators in parent")

        #Check to make sure that we aren't dividing by zero
        cdef Py_ssize_t m
        for m in range(l):
            if x[m] == 0:
                if self.has_inverse_of(m):
                    raise ZeroDivisionError

        ans = self._poly(*x)
        if ans:
            for m in self._mon.nonzero_positions():
                ans *= x[m]**self._mon[m]

        return ans

    def subs(self, in_dict=None, **kwds):
        """
        Substitute some variables in this Laurent polynomial.

        Variable/value pairs for the substitution may be given
        as a dictionary or via keyword-value pairs. If both are
        present, the latter take precedence.

        INPUT:

        - ``in_dict`` -- dictionary (optional)

        - ``**kwargs`` -- keyword arguments

        OUTPUT:

        A Laurent polynomial.

        EXAMPLES::

            sage: L.<x, y, z> = LaurentPolynomialRing(QQ)
            sage: f = x + 2*y + 3*z
            sage: f.subs(x=1)
            2*y + 3*z + 1
            sage: f.subs(y=1)
            x + 3*z + 2
            sage: f.subs(z=1)
            x + 2*y + 3
            sage: f.subs(x=1, y=1, z=1)
            6

            sage: f = x^-1
            sage: f.subs(x=2)
            1/2
            sage: f.subs({x: 2})
            1/2

            sage: f = x + 2*y + 3*z
            sage: f.subs({x: 1, y: 1, z: 1})
            6
            sage: f.substitute(x=1, y=1, z=1)
            6

        TESTS::

            sage: f = x + 2*y + 3*z
            sage: f(q=10)
            x + 2*y + 3*z

            sage: x.subs({x: 2}, x=1)
            1
        """
        cdef list variables = list(self._parent.gens())
        cdef Py_ssize_t i
        for i in range(len(variables)):
            if str(variables[i]) in kwds:
                variables[i] = kwds[str(variables[i])]
            elif in_dict and variables[i] in in_dict:
                variables[i] = in_dict[variables[i]]
        return self(tuple(variables))

    def is_constant(self):
        r"""
        Return whether this Laurent polynomial is constant.

        EXAMPLES::

            sage: L.<a, b> = LaurentPolynomialRing(QQ)
            sage: L(0).is_constant()
            True
            sage: L(42).is_constant()
            True
            sage: a.is_constant()
            False
            sage: (1/b).is_constant()
            False
        """
        return (self._mon == ETuple({}, int(self._parent.ngens())) and
                self._poly.is_constant())

    def _symbolic_(self, R):
        """
        EXAMPLES::

            sage: R.<x,y> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + y/x
            sage: g = f._symbolic_(SR); g
            (x^4 + y)/x
            sage: g(x=2,y=2)
            9

            sage: g = SR(f)
            sage: g(x=2,y=2)
            9
        """
        d = {repr(g): R.var(g) for g in self._parent.gens()}
        return self.subs(**d)

    def derivative(self, *args):
        r"""
        The formal derivative of this Laurent polynomial, with respect
        to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global :func:`derivative` function for more
        details.

        .. SEEALSO::

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
        is with respect to the generator. Otherwise, ``_derivative(var)`` is called
        recursively for each coefficient of this polynomial.

        .. SEEALSO:: :meth:`derivative`

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
        P = self._parent
        cdef list gens = list(P.gens())

        # check if var is one of the generators
        try:
            index = gens.index(var)
        except ValueError:
            # call _derivative() recursively on coefficients
            return P({m: c._derivative(var)
                      for (m, c) in self.dict().iteritems()})

        # compute formal derivative with respect to generator
        cdef dict d = {}
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

        - ``R`` - (default: ``None``) a univariate Laurent polynomial ring

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
            x = str(x)
        else:
            x = 'x'
            i = 0

        #construct ring if none
        if R is None:
            R = LaurentPolynomialRing(self.base_ring(), x)

        return R({m[i]: c for m,c in self.dict().iteritems()})

    def factor(self):
        """
        Returns a Laurent monomial (the unit part of the factorization) and a factored multi-polynomial.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: f = 4*x^7*z^-1 + 3*x^3*y + 2*x^4*z^-2 + x^6*y^-7
            sage: f.factor()
            (x^3*y^-7*z^-2) * (4*x^4*y^7*z + 3*y^8*z^2 + 2*x*y^7 + x^3*z^2)

        TESTS:

        Tests for :trac:`29173`::

            sage: L.<a, b> = LaurentPolynomialRing(ZZ, 'a, b')
            sage: (a*b + a + b + 1).factor()
            (b + 1) * (a + 1)
            sage: ((a^-1)*(a*b + a + b + 1)).factor()
            (a^-1) * (b + 1) * (a + 1)
            sage: L(-12).factor()
            -1 * 2^2 * 3
        """
        pf = self._poly.factor()

        if self._poly.degree() == 0:
            # Factorization is broken for polynomials, see
            # https://github.com/sagemath/sage/issues/20214
            return pf

        u = self.parent(pf.unit())

        cdef tuple g = <tuple> self._parent.gens()
        for i in self._mon.nonzero_positions():
            u *= g[i] ** self._mon[i]

        cdef list f = []
        cdef dict d
        for t in pf:
            d = <dict> (t[0].dict())
            if len(d) == 1:  # monomials are units
                u *= self.parent(d) ** t[1]
            else:
                f.append((self.parent(d), t[1]))

        return Factorization(f, unit=u)

    def is_square(self, root=False):
        r"""
        Test whether this Laurent polynomial is a square.

        INPUT:

        - ``root`` - boolean (default ``False``) - if set to ``True``
          then return a pair ``(True, sqrt)`` with ``sqrt`` a square
          root of this Laurent polynomial when it exists or
          ``(False, None)``.

        EXAMPLES::

            sage: L.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: p = (1 + x*y + z^-3)
            sage: (p**2).is_square()
            True
            sage: (p**2).is_square(root=True)
            (True, x*y + 1 + z^-3)

            sage: x.is_square()
            False
            sage: x.is_square(root=True)
            (False, None)

            sage: (x**-4 * (1 + z)).is_square(root=False)
            False
            sage: (x**-4 * (1 + z)).is_square(root=True)
            (False, None)
        """
        self._normalize()
        if not self._mon.is_multiple_of(2):
            return (False, None) if root else False

        cdef LaurentPolynomial_mpair ans

        if not root:
            return self._poly.is_square(root=False)
        else:
            (pans, root) = self._poly.is_square(root=True)
            if not pans:
                return (False, None)

            mon = self._mon.escalar_div(2)
            ans = self._new_c()
            ans._mon = mon
            ans._poly = root
            return (True, ans)

    cpdef rescale_vars(self, dict d, h=None, new_ring=None):
        r"""
        Rescale variables in a Laurent polynomial.

        INPUT:

        - ``d`` -- a ``dict`` whose keys are the generator indices
          and values are the coefficients; so a pair ``(i, v)``
          means `x_i \mapsto v x_i`
        - ``h`` -- (optional) a map to be applied to coefficients
          done after rescaling
        - ``new_ring`` -- (optional) a new ring to map the result into

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = x^-2*y + x*y^-2
            sage: p.rescale_vars({0: 2, 1: 3})
            2/9*x*y^-2 + 3/4*x^-2*y
            sage: F = GF(2)
            sage: p.rescale_vars({0: 3, 1: 7}, new_ring=L.change_ring(F))
            x*y^-2 + x^-2*y

        Test for :trac:`30331`::

            sage: F.<z> = CyclotomicField(3)
            sage: p.rescale_vars({0: 2, 1: z}, new_ring=L.change_ring(F))
            2*z*x*y^-2 + 1/4*z*x^-2*y
        """
        cdef int i
        cdef dict df
        cdef ETuple v
        cdef LaurentPolynomial_mpair ans

        if self._prod is None:
            self._compute_polydict()

        df = dict(self._prod.__repn)  # This makes a copy for us to manipulate
        if new_ring is None:
            R = self._parent._base
        else:
            R = new_ring._base
        if h is None:
            for v in df:
                val = df[v]
                for i in d:
                    val *= d[i]**v[i]
                df[v] = val
        else:
            for v in df:
                val = df[v]
                for i in d:
                    val *= d[i]**v[i]
                df[v] = R(h(val))

        ans = <LaurentPolynomial_mpair> self._new_c()
        ans._prod = PolyDict(df)
        ans._mon = self._mon
        if new_ring is None:
            S = self._poly._parent
        else:
            S = self._poly._parent.change_ring(R)
        ans._poly = <MPolynomial> S({v.esub(ans._mon): df[v] for v in df})
        if new_ring is not None:
            return new_ring(ans)
        return ans

    cpdef toric_coordinate_change(self, M, h=None, new_ring=None):
        r"""
        Apply a matrix to the exponents in a Laurent polynomial.

        For efficiency, we implement this directly, rather than as a substitution.

        The optional argument ``h`` is a map to be applied to coefficients.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = 2*x^2 + y - x*y
            sage: p.toric_coordinate_change(Matrix([[1,-3],[1,1]]))
            2*x^2*y^2 - x^-2*y^2 + x^-3*y
            sage: F = GF(2)
            sage: p.toric_coordinate_change(Matrix([[1,-3],[1,1]]), new_ring=L.change_ring(F))
            x^-2*y^2 + x^-3*y

        """
        cdef int n, i, j, x
        cdef dict d, dr
        cdef ETuple v
        cdef LaurentPolynomial_mpair ans
        cdef list L, mon, exp
        cdef Matrix mat = M

        n = self._parent.ngens()
        if mat.dimensions() != (n, n):
            raise ValueError("the matrix M must be a {k} x {k} matrix".format(k=n))

        if not self:
            if new_ring is None:
                return self._parent.zero()
            else:
                return new_ring.zero()

        if self._prod is None:
            self._compute_polydict()

        d = self._prod.__repn
        dr = {}
        mon = [0] * n
        for v in d:
            # Make a copy of mon as this might be faster than creating the data from scratch.
            # We will set every entry, so no need to clear the data.
            exp = list(mon)
            for j in range(n):
                x = 0
                for i in range(n):
                    if not mat.get_is_zero_unsafe(j, i):
                        x += (<int> v[i]) * int(mat.get_unsafe(j, i))
                if x < (<int> mon[j]):
                    mon[j] = x
                exp[j] = x
            dr[ETuple(exp)] = d[v]

        if h is not None:
            for v in dr:
                dr[v] = self._parent._base(h(dr[v]))

        ans = <LaurentPolynomial_mpair> self._new_c()
        ans._prod = PolyDict(dr)
        ans._mon = ETuple(mon)
        ans._poly = <MPolynomial> self._poly._parent({v.esub(ans._mon): dr[v] for v in dr})
        if new_ring is not None:
            return new_ring(ans)
        return ans

    cpdef toric_substitute(self, v, v1, a, h=None, new_ring=None):
        r"""
        Perform a single-variable substitution up to a toric coordinate change.

        The optional argument ``h`` is a map to be applied to coefficients.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = x + y
            sage: p.toric_substitute((2,3), (-1,1), 2)
            1/2*x^3*y^3 + 2*x^-2*y^-2
            sage: F = GF(5)
            sage: p.toric_substitute((2,3), (-1,1), 2, new_ring=L.change_ring(F))
            3*x^3*y^3 + 2*x^-2*y^-2

        TESTS:

        Tests for :trac:`30331`::

            sage: L.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: p = x + y
            sage: F.<z> = CyclotomicField(3)
            sage: p.toric_substitute((2,3), (-1,1), z, new_ring=L.change_ring(F))
            (-z - 1)*x^3*y^3 + z*x^-2*y^-2

            sage: P.<x> = LaurentPolynomialRing(QQ, 1)
            sage: u = x - 1
            sage: v = u.toric_substitute((-1,), (-1,), 1)
            sage: v.is_zero()
            True
        """
        cdef dict d, dr
        cdef ETuple ve, v1e, w, w1, mon
        cdef LaurentPolynomial_mpair ans
        cdef int t

        if self._prod is None:
            self._compute_polydict()

        d = self._prod.__repn
        dr = {}
        ve = ETuple(v)
        v1e = ETuple(v1)
        mon = self._mon
        if h is not None:
            d = dict(d)  # Make a copy so we can manipulate it
            for w in d:
                d[w] = h(d[w])
        for w in d:
            x = d[w]
            t = w.dotprod(v1e)
            w1 = w.eadd_scaled(ve, -t)
            if w1 in dr:
                dr[w1] += x * a**t
            else:
                dr[w1] = x * a**t
            mon = mon.emin(w1)
        for v in tuple(dr.keys()):
            if not dr[v]:
                del dr[v]

        if new_ring is None:
            S = self._poly._parent
        else:
            S = self._poly._parent.change_ring(new_ring._base)
        ans = <LaurentPolynomial_mpair> self._new_c()
        ans._prod = PolyDict(dr)
        ans._mon = mon
        ans._poly = <MPolynomial> S({v.esub(ans._mon): dr[v] for v in dr})
        if new_ring is not None:
            return new_ring(ans)
        return ans
