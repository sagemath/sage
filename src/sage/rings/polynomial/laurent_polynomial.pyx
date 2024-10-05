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
from sage.categories.map cimport Map
from sage.structure.element import coerce_binop, parent
from sage.structure.factorization import Factorization
from sage.misc.derivative import multi_derivative
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.structure.richcmp cimport richcmp, rich_to_bool
from sage.rings.infinity import minus_infinity


cdef class LaurentPolynomial(CommutativeAlgebraElement):
    """
    Base class for Laurent polynomials.
    """
    cdef LaurentPolynomial _new_c(self):
        """
        Return a new Laurent polynomial.

        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)  # indirect doctest               # needs sage.modules
            sage: x*y                                                                   # needs sage.modules
            x*y
        """
        cdef type t = type(self)
        cdef LaurentPolynomial ans
        ans = t.__new__(t)
        ans._parent = self._parent
        return ans

    cpdef _add_(self, other):
        """
        Abstract addition method.

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
        Abstract multiplication method.

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
        Abstract floor division method.

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

        OUTPUT: integer

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

            sage: # needs sage.modules
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

        OUTPUT: a rational

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

            sage: # needs sage.modules
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

        Check that :issue:`22277` is fixed::

            sage: # needs sage.modules
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

        The hamming weight is number of nonzero coefficients and
        also known as the weight or sparsity.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: f = x^3 - 1
            sage: f.hamming_weight()
            2
        """
        return self.number_of_terms()

    cpdef dict monomial_coefficients(self):
        """
        Abstract ``dict`` method.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial
            sage: LaurentPolynomial.monomial_coefficients(x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    dict = monomial_coefficients

    def map_coefficients(self, f, new_base_ring=None):
        """
        Apply ``f`` to the coefficients of ``self``.

        If ``f`` is a :class:`sage.categories.map.Map`, then the resulting
        polynomial will be defined over the codomain of ``f``. Otherwise, the
        resulting polynomial will be over the same ring as ``self``. Set
        ``new_base_ring`` to override this behavior.

        INPUT:

        - ``f`` -- a callable that will be applied to the coefficients of ``self``

        - ``new_base_ring`` -- (optional) if given, the resulting polynomial
          will be defined over this ring

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(9)
            sage: R.<x> = LaurentPolynomialRing(k)
            sage: f = x*a + a
            sage: f.map_coefficients(lambda a: a + 1)
            (a + 1) + (a + 1)*x
            sage: R.<x,y> = LaurentPolynomialRing(k, 2)                                 # needs sage.modules
            sage: f = x*a + 2*x^3*y*a + a                                               # needs sage.modules
            sage: f.map_coefficients(lambda a: a + 1)                                   # needs sage.modules
            (2*a + 1)*x^3*y + (a + 1)*x + a + 1

        Examples with different base ring::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: R.<r> = GF(9); S.<s> = GF(81)
            sage: h = Hom(R, S)[0]; h
            Ring morphism:
              From: Finite Field in r of size 3^2
              To:   Finite Field in s of size 3^4
              Defn: r |--> 2*s^3 + 2*s^2 + 1
            sage: T.<X,Y> = LaurentPolynomialRing(R, 2)
            sage: f = r*X + Y
            sage: g = f.map_coefficients(h); g
            (2*s^3 + 2*s^2 + 1)*X + Y
            sage: g.parent()
            Multivariate Laurent Polynomial Ring in X, Y
             over Finite Field in s of size 3^4
            sage: h = lambda x: x.trace()
            sage: g = f.map_coefficients(h); g
            X - Y
            sage: g.parent()
            Multivariate Laurent Polynomial Ring in X, Y
             over Finite Field in r of size 3^2
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
        return R(dict([(k, f(v))
                       for k, v in self.monomial_coefficients().items()]))


cdef class LaurentPolynomial_univariate(LaurentPolynomial):
    r"""
    A univariate Laurent polynomial in the form of `t^n \cdot f`
    where `f` is a polynomial in `t`.

    INPUT:

    - ``parent`` -- a Laurent polynomial ring

    - ``f`` -- a polynomial (or something that can be coerced to one)

    - ``n`` -- integer (default: 0)

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

            sage: # needs sage.rings.padics
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
            n += (< LaurentPolynomial_univariate > f).__n
            if (< LaurentPolynomial_univariate > f).__u._parent is parent._R:
                f = (< LaurentPolynomial_univariate > f).__u
            else:
                f = parent._R((< LaurentPolynomial_univariate > f).__u)
        elif (not isinstance(f, Polynomial)) or (parent is not f.parent()):
            if isinstance(f, dict):
                v = min(f) if f else 0
                f = {i-v: c for i, c in f.items()}
                n += v
            f = parent._R(f)

        # self is that t^n * u:
        self.__u = f
        self.__n = n
        self._normalize()

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

            sage: Pxes = [(Px, Px.gen()), (Qx, Qx.gen()),
            ....:         (Pxy, Pxy.gen(0)), (Paxb, Paxb.gen(1))]
            sage: Pxes += [(Rx, Rx.gen())]
            sage: for P, x in Pxes:
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

        if isinstance(R, PolynomialRing_general):
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
            sage: (2 + t).is_unit()
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
        Check if ``self`` is nonzero.

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

            sage: # needs sage.rings.number_field
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

    cpdef _normalize(self):
        r"""
        A Laurent series is a pair `(u(t), n)`, where either `u = 0`
        (to some precision) or `u` is a unit. This pair corresponds to
        `t^n \cdot u(t)`.

        EXAMPLES::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: elt = t^2 + t^4  # indirect doctest
            sage: elt.polynomial_construction()
            (t^2 + 1, 2)

        Check that :issue:`21272` is fixed::

            sage: (t - t).polynomial_construction()
            (0, 0)
        """
        if self.__u[0]:
            return
        elif self.__u.is_zero():
            self.__n = 0
            return
        # we already caught the infinity and zero cases
        cdef long v = <long > self.__u.valuation()
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
        for n in range(m):
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
                    var = "*{}^{}".format(X, e)
                s += "{}{}".format(x, var)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*", " ")
        s = s.replace(" -1*", " -")
        return s[1:]

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = (17/2)*x^-2 + x + x^2 + 3*x^4
            sage: latex(f)
            \frac{\frac{17}{2}}{x^{2}} + x + x^{2} + 3x^{4}

        Verify that :issue:`6656` has been fixed::

            sage: R.<a,b> = PolynomialRing(QQ)
            sage: T.<x> = LaurentPolynomialRing(R)
            sage: y = a*x + b*x
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
        for n in range(m):
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
                    var = "|{}^{{{}}}".format(X, e)
                if e >= 0:
                    s += "{}{}".format(x, var)
                else:  # negative e
                    if e == -1:
                        s += "\\frac{{{}}}{{{}}}".format(x, X)
                    else:
                        s += "\\frac{{{}}}{{{}^{{{}}}}}".format(x, X, -e)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1|", " ")
        s = s.replace(" -1|", " -")
        s = s.replace("|", "")

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

        Check that :issue:`21272` is fixed::

            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: hash(R.zero()) == hash(t - t)
            True
        """
        # we reimplement below the hash of polynomials to handle negative
        # degrees
        cdef long result = 0
        cdef long result_mon
        cdef int i, j
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

        Any other kind of slicing is an error, see :issue:`18940`::

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
            ret = <LaurentPolynomial_univariate > self._new_c()
            ret.__u = f
            ret.__n = self.__n
            ret._normalize()
            return ret

        return self.__u[i - self.__n]

    cpdef long number_of_terms(self) except -1:
        """
        Return the number of nonzero coefficients of ``self``.

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

            sage: # needs sage.symbolic
            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = x^3 + 2/x
            sage: g = f._symbolic_(SR); g
            (x^4 + 2)/x
            sage: g(x=2)
            9
            sage: g = SR(f)
            sage: g(x=2)
            9

        Since :issue:`24072` the symbolic ring does not accept positive
        characteristic::

            sage: R.<w> = LaurentPolynomialRing(GF(7))
            sage: SR(2*w^3 + 1)                                                         # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: positive characteristic not allowed in symbolic computations
        """
        d = {repr(g): R.var(g) for g in self._parent.gens()}
        return self.subs(**d)

    cpdef dict monomial_coefficients(self):
        """
        Return a dictionary representing ``self``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: Q.<t> = LaurentPolynomialRing(R)
            sage: f = (x^3 + y/t^3)^3 + t^2; f
            y^3*t^-9 + 3*x^3*y^2*t^-6 + 3*x^6*y*t^-3 + x^9 + t^2
            sage: f.monomial_coefficients()
            {-9: y^3, -6: 3*x^3*y^2, -3: 3*x^6*y, 0: x^9, 2: 1}

        ``dict`` is an alias::

            sage: f.dict()
            {-9: y^3, -6: 3*x^3*y^2, -3: 3*x^6*y, 0: x^9, 2: 1}
        """
        cdef dict d = self.__u.monomial_coefficients()
        return {k + self.__n: d[k] for k in d}

    dict = monomial_coefficients

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
        else:  # off to the left
            if value != 0:
                self.__n = self.__n + j
                R = self._parent.base_ring()
                coeffs = [value] + [R.zero() for _ in range(1, -j)] + self.__u.list()
                self.__u = self.__u._parent(coeffs)
        self._normalize()

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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate > right_m
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
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > (f1 + f2)
        ret.__n = m
        ret._normalize()
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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate > right_m
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
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > (f1 - f2)
        ret.__n = m
        ret._normalize()
        return ret

    def degree(self):
        r"""
        Return the degree of ``self``.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: g = x^2 - x^4
            sage: g.degree()
            4
            sage: g = -10/x^5 + x^2 - x^7
            sage: g.degree()
            7

        The zero polynomial is defined to have degree `-\infty`::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: R.zero().degree()
            -Infinity
        """
        # The zero polynomial is defined to have degree -Infinity
        if self.is_zero():
            return minus_infinity
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
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > -self.__u
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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate > right_r
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > (self.__u * right.__u)
        ret.__n = self.__n + right.__n
        ret._normalize()
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
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > self.__u._rmul_(c)
        ret.__n = self.__n
        ret._normalize()
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
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > self.__u._lmul_(c)
        ret.__n = self.__n
        ret._normalize()
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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate > rhs
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > (self.__u // right.__u)
        ret.__n = self.__n - right.__n
        ret._normalize()
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
        ret = <LaurentPolynomial_univariate > self._new_c()
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
        ret = <LaurentPolynomial_univariate > self._new_c()
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
        ret = <LaurentPolynomial_univariate > self._new_c()
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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate > rhs
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
        if self.__u.is_constant():  # this has a single term c*x^n
            ret = <LaurentPolynomial_univariate > self._new_c()
            if self.__u.is_unit():
                ret.__u = self.__u.inverse_of_unit()
                ret.__n = -self.__n
                ret._normalize()
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

    @coerce_binop
    def xgcd(self, other):
        r"""
        Extended :meth:`gcd` for univariate Laurent polynomial rings over a field.

        OUTPUT:

        A triple ``(g, p, q)`` such that ``g`` is the :meth:`gcd` of
        ``self`` (`= a`) and ``other`` (`= b`), and ``p`` and ``q`` are
        cofactors satisfying the Bezout identity

        .. MATH::

            g = p \cdot a + q \cdot b.

        EXAMPLES::

            sage: S.<t> = LaurentPolynomialRing(QQ)
            sage: a = t^-2 + 1
            sage: b = t^-3 + 1
            sage: g, p, q = a.xgcd(b); (g, p, q)
            (t^-3, 1/2*t^-1 - 1/2 - 1/2*t, 1/2 + 1/2*t)
            sage: g == p * a + q * b
            True
            sage: g == a.gcd(b)
            True
            sage: t.xgcd(t)
            (t, 0, 1)
            sage: t.xgcd(5)
            (1, 0, 1/5)
        """
        cdef LaurentPolynomial_univariate elt = other
        cdef LaurentPolynomial_univariate ret_gcd, ret_p, ret_q
        cdef long n = min(self.__n, elt.__n)

        h, p, q = self.__u.xgcd(elt.__u)

        ret_gcd = <LaurentPolynomial_univariate> self._new_c()
        ret_gcd.__u = h
        ret_gcd.__n = n
        ret_gcd._normalize()

        ret_p = <LaurentPolynomial_univariate> self._new_c()
        ret_p.__u = p
        ret_p.__n = n - self.__n
        ret_p._normalize()

        ret_q = <LaurentPolynomial_univariate> self._new_c()
        ret_q.__u = q
        ret_q.__n = n - elt.__n
        ret_q._normalize()

        return (ret_gcd, ret_p, ret_q)

    def inverse_mod(a, m):
        """
        Invert the polynomial ``a`` with respect to ``m``, or raise a :exc:`ValueError`
        if no such inverse exists.

        The parameter ``m`` may be either a single polynomial or an ideal
        (for consistency with :meth:`inverse_mod` in other rings).

        ALGORITHM: Solve the system `as + mt = 1`, returning `s` as the inverse
        of `a` mod `m`.

        EXAMPLES::

            sage: S.<t> = LaurentPolynomialRing(QQ)
            sage: f = inverse_mod(t^-2 + 1, t^-3 + 1); f
            1/2*t^2 - 1/2*t^3 - 1/2*t^4
            sage: f * (t^-2 + 1) + (1/2*t^4 + 1/2*t^3) * (t^-3 + 1)
            1
        """
        from sage.rings.ideal import Ideal_generic
        if isinstance(m, Ideal_generic):
            v = m.gens_reduced()
            if len(v) > 1:
                raise NotImplementedError("only inversion modulo principal ideals implemented")
            m = v[0]
        if m.degree() == 1 and m[1].is_unit():
            # a(x) mod (x-r) = a(r)
            r = -m[0]
            if not m[1].is_one():
                r *= m.base_ring()(~m[1])
            u = a(r)
            if u.is_unit():
                return a.parent()(~u)
        g, s, _ = a.xgcd(m)
        if g == 1:
            return s
        elif g.is_unit():
            return g.inverse_of_unit() * s
        raise ValueError("impossible inverse modulo")

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
        b = <LaurentPolynomial_univariate > self._parent(right)
        cdef LaurentPolynomial_univariate ret
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = self.__u.gcd(b.__u)
        ret.__n = min(self.__n, b.__n)
        ret._normalize()
        return ret

    def euclidean_degree(self):
        r"""
        Return the degree of ``self`` as an element of an Euclidean domain.

        This is the Euclidean degree of the underlying polynomial.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: (x^-5 + x^2).euclidean_degree()
            7

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: (x^-5 + x^2).euclidean_degree()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.__u.euclidean_degree()

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

        Check that :issue:`34330` is fixed::

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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate > other
        q, r = self.__u.quo_rem(right.__u)
        cdef LaurentPolynomial_univariate ql, qr
        ql = <LaurentPolynomial_univariate > self._new_c()
        ql.__u = <ModuleElement > q
        ql.__n = self.__n - right.__n
        ql._normalize()
        qr = <LaurentPolynomial_univariate > self._new_c()
        qr.__u = <ModuleElement > r
        qr.__n = self.__n
        qr._normalize()
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
        cdef LaurentPolynomial_univariate right = <LaurentPolynomial_univariate > right_r

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
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > self.__u.truncate(n - self.__n)
        ret.__n = self.__n
        ret._normalize()
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

    def monomial_reduction(self):
        """
        Return the decomposition as a polynomial and a power of the variable.
        Constructed for compatibility with the multivariate case.

        OUTPUT:

        A tuple ``(u, t^n)`` where ``u`` is the underlying polynomial and ``n``
        is the power of the exponent shift.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4
            sage: f.monomial_reduction()
            (3*x^5 + x^3 + 1, x^-1)
        """
        return (self.__u, self._parent.gen(0) ** self.__n)

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

            sage: # needs sage.rings.number_field
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
        ret = <LaurentPolynomial_univariate > self._new_c()
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

        Check that :issue:`28187` is fixed::

            sage: # needs sage.symbolic
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
                ret = <LaurentPolynomial_univariate > self._new_c()
                ret.__u = <ModuleElement > self._parent._R(u)
                ret.__n = self.__n
                ret._normalize()
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
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > self._parent._R(a)
        ret.__n = self.__n - 1
        ret._normalize()
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
            v = [a[i]/(n+i+1) for i in range(min(-1-n, len(a)))] + [0]
        else:
            v = []
        v += [a[i]/(n+i+1) for i in range(max(-n, 0), len(a))]
        try:
            u = self._parent._R(v)
        except TypeError:
            raise ArithmeticError("coefficients of integral cannot be coerced into the base ring")
        ret = <LaurentPolynomial_univariate > self._new_c()
        ret.__u = <ModuleElement > u
        ret.__n = n + 1
        ret._normalize()
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
            if x:  # If there are non-keyword arguments
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
            sage: f.factor()                                                            # needs sage.libs.pari
            (t^-7) * (4 + t + 3*t^10 + 2*t^11)
        """
        cdef LaurentPolynomial_univariate u, d
        pf = self.__u.factor()
        u = <LaurentPolynomial_univariate > self._new_c()
        u.__u = pf.unit()
        u.__n = self.__n
        u._normalize()

        f = []
        for t in pf:
            d = <LaurentPolynomial_univariate > self._new_c()
            d.__u = t[0]
            d.__n = 0
            d._normalize()
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

    def _as_extended_polynomial(self):
        """
        This Laurent polynomial seen as a polynomial in twice as many variables,
        where half of the variables are the inverses of the other half.

        EXAMPLES::

            sage: L.<t> = LaurentPolynomialRing(QQ)
            sage: f = t-t^-2
            sage: f._as_extended_polynomial()
            -tinv^2 + t
            sage: _.parent()
            Multivariate Polynomial Ring in t, tinv over Rational Field
        """
        dres = {}
        for e, c in self.monomial_coefficients().items():
            if e > 0:
                dres[(e, 0)] = c
            else:
                dres[(0, -e)] = c
        return self.parent()._extended_ring(dres)

    @coerce_binop
    def divides(self, other):
        r"""
        Return ``True`` if ``self`` divides ``other``.

        EXAMPLES::

            sage: R.<x> = LaurentPolynomialRing(ZZ)
            sage: (2*x**-1 + 1).divides(4*x**-2 - 1)
            True
            sage: (2*x + 1).divides(4*x**2 + 1)
            False
            sage: (2*x + x**-1).divides(R(0))
            True
            sage: R(0).divides(2*x ** -1 + 1)
            False
            sage: R(0).divides(R(0))
            True
            sage: R.<x> = LaurentPolynomialRing(Zmod(6))
            sage: p = 4*x + 3*x^-1
            sage: q = 5*x^2 + x + 2*x^-2
            sage: p.divides(q)
            False

            sage: R.<x,y> = GF(2)[]
            sage: S.<z> = LaurentPolynomialRing(R)
            sage: p = (x+y+1) * z**-1 + x*y
            sage: q = (y^2-x^2) * z**-2 + z + x-y
            sage: p.divides(q), p.divides(p*q)                                          # needs sage.libs.singular
            (False, True)
        """
        p = self.polynomial_construction()[0]
        q = other.polynomial_construction()[0]
        return p.divides(q)
