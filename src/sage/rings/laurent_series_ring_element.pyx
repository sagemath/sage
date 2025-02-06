"""
Laurent Series

EXAMPLES::

    sage: R.<t> = LaurentSeriesRing(GF(7), 't'); R
    Laurent Series Ring in t over Finite Field of size 7
    sage: f = 1/(1-t+O(t^10)); f
    1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)

Laurent series are immutable::

    sage: f[2]
    1
    sage: f[2] = 5
    Traceback (most recent call last):
    ...
    IndexError: Laurent series are immutable

We compute with a Laurent series over the complex mpfr numbers.

::

    sage: K.<q> = Frac(CC[['q']]); K                                                    # needs sage.rings.real_mpfr
    Laurent Series Ring in q over Complex Field with 53 bits of precision
    sage: q                                                                             # needs sage.rings.real_mpfr
    1.00000000000000*q

Saving and loading.

::

    sage: loads(q.dumps()) == q                                                         # needs sage.rings.real_mpfr
    True
    sage: loads(K.dumps()) == K                                                         # needs sage.rings.real_mpfr
    True

IMPLEMENTATION: Laurent series in Sage are represented internally
as a power of the variable times the unit part (which need not be a
unit - it's a polynomial with nonzero constant term). The zero
Laurent series has unit part 0.

AUTHORS:

- William Stein: original version

- David Joyner (2006-01-22): added examples

- Robert Bradshaw (2007-04): optimizations, shifting

- Robert Bradshaw: Cython version
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.infinity import infinity

from sage.rings.rational_field import QQ
import sage.misc.latex
from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial_univariate
from sage.rings.power_series_ring_element cimport PowerSeries
from sage.structure.element cimport Element, AlgebraElement
from sage.structure.richcmp cimport richcmp_not_equal, rich_to_bool
from sage.misc.derivative import multi_derivative


def is_LaurentSeries(x):
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38266,
                       "The function is_LaurentSeries is deprecated; "
                       "use 'isinstance(..., LaurentSeries)' instead.")
    return isinstance(x, LaurentSeries)


cdef class LaurentSeries(AlgebraElement):
    r"""
    A Laurent Series.

    We consider a Laurent series of the form `t^n \cdot f` where `f` is a
    power series.

    INPUT:

    - ``parent`` -- a Laurent series ring

    - ``f`` -- a power series (or something can be coerced
      to one); note that ``f`` does *not* have to be a unit

    - ``n`` -- (default: 0) integer
    """
    def __init__(self, parent, f, n=0):
        r"""
        Initialize ``self``.

        OUTPUT: a Laurent series

        EXAMPLES::

            sage: R.<q> = LaurentSeriesRing(ZZ)
            sage: R([1,2,3])
            1 + 2*q + 3*q^2
            sage: R([1,2,3], -5)
            q^-5 + 2*q^-4 + 3*q^-3

        ::

            sage: # needs sage.rings.finite_rings sage.rings.padics
            sage: S.<s> = LaurentSeriesRing(GF(5))
            sage: T.<t> = PowerSeriesRing(pAdicRing(5))
            sage: S(t)
            s
            sage: parent(S(t))
            Laurent Series Ring in s over Finite Field of size 5
            sage: parent(S(t)[1])
            Finite Field of size 5
        """
        AlgebraElement.__init__(self, parent)

        if isinstance(f, LaurentSeries):
            n += (<LaurentSeries>f).__n
            if (<LaurentSeries>f).__u._parent is parent._power_series_ring:
                f = (<LaurentSeries>f).__u
            else:
                f = parent._power_series_ring((<LaurentSeries>f).__u)
        elif isinstance(f, LaurentPolynomial_univariate):
            f = f(parent.gen())
        elif isinstance(f, dict):
            ## Sanitize input to make sure all exponents are nonnegative,
            ## adjusting n to match.
            n1 = min(f.keys())
            if n1 < 0:
                f = {e - n1: c for e, c in f.items()}
                n += n1
            f = parent._power_series_ring(f)
        elif not isinstance(f, PowerSeries):
            f = parent._power_series_ring(f)
        ## now this is a power series, over a different ring ...
        ## requires that power series rings with same vars over the
        ## same parent are unique.
        elif parent is not f.parent():
            f = parent._power_series_ring(f)

        # self is that t^n * u:
        if not f:
            if n is infinity:
                self.__n = 0
                self.__u = parent._power_series_ring.zero()
            else:
                self.__n = n
                self.__u = f
        else:
            val = f.valuation()
            if val is infinity:
                self.__n = 0
                self.__u = f
            elif val == 0:
                self.__n = n    # power of the variable
                self.__u = f    # unit part
            else:
                self.__n = n + val
                self.__u = f >> val

    def __reduce__(self):
        return self._parent, (self.__u, self.__n)

    def change_ring(self, R):
        """
        Change the base ring of ``self``.

        EXAMPLES::

            sage: R.<q> = LaurentSeriesRing(ZZ)
            sage: p = R([1,2,3]); p
            1 + 2*q + 3*q^2
            sage: p.change_ring(GF(2))
            1 + q^2
        """
        return self._parent.change_ring(R)(self)

    def is_unit(self):
        """
        Return ``True`` if this is Laurent series is a unit in this ring.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: (2 + t).is_unit()
            True
            sage: f = 2 + t^2 + O(t^10); f.is_unit()
            True
            sage: 1/f
            1/2 - 1/4*t^2 + 1/8*t^4 - 1/16*t^6 + 1/32*t^8 + O(t^10)
            sage: R(0).is_unit()
            False
            sage: R.<s> = LaurentSeriesRing(ZZ)
            sage: f = 2 + s^2 + O(s^10)
            sage: f.is_unit()
            False
            sage: 1/f
            Traceback (most recent call last):
            ...
            ValueError: constant term 2 is not a unit

        ALGORITHM: A Laurent series is a unit if and only if its "unit
        part" is a unit.
        """
        return self.__u.is_unit()

    def is_zero(self):
        """
        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x + x + x^2 + 3*x^4 + O(x^7)
            sage: f.is_zero()
            0
            sage: z = 0*f
            sage: z.is_zero()
            1
        """
        return self.__u.is_zero()

    def is_monomial(self):
        """
        Return ``True`` if this element is a monomial.  That is, if ``self`` is
        `x^n` for some integer `n`.

        EXAMPLES::

            sage: k.<z> = LaurentSeriesRing(QQ, 'z')
            sage: (30*z).is_monomial()
            False
            sage: k(1).is_monomial()
            True
            sage: (z+1).is_monomial()
            False
            sage: (z^-2909).is_monomial()
            True
            sage: (3*z^-2909).is_monomial()
            False
        """
        return self.__u.is_monomial()

    def __bool__(self):
        """
        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: bool(t)
            True
            sage: bool(1/t)
            True
            sage: bool(2 + t)
            True
            sage: bool(1/(1-t))
            True
            sage: bool(1 + O(t^3))
            True
            sage: bool(t + O(t^3))
            True
            sage: bool(O(t^3))
            False
            sage: bool(O(t^-3))
            False
            sage: bool(R.zero())
            False
        """
        return bool(self.__u)

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of this series under the map that sends the generators of
        the parent to ``im_gens``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: Zx.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: R.<t> = LaurentSeriesRing(K)
            sage: z = t^-1 + i*t
            sage: z._im_gens_(R, [t^2])
            t^-2 + i*t^2

        The argument ``base_map`` is not yet supported, because it isn't over power series::

            sage: cc = K.hom([i])                                                       # needs sage.rings.number_field
            sage: z._im_gens_(R, [t^2], base_map=cc)                                    # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        x = im_gens[0]
        return codomain(self.__u._im_gens_(codomain, im_gens, base_map=base_map) * x**self.__n)

    cdef _normalize(self):
        r"""
        A Laurent series is a pair (u(t), n), where either u=0 (to some
        precision) or u is a unit. This pair corresponds to
        `t^n\cdot u(t)`.
        """
        if self.is_zero():
            return
        v = self.__u.valuation()
        if v == 0:
            return
        self.__n += v
        self.__u = self.__u.valuation_zero_part()

    def _repr_(self):
        """
        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: (2 + (2/3)*t^3).__repr__()
            '2 + 2/3*t^3'
        """
        if self.is_zero():
            if self.prec() is infinity:
                return "0"
            else:
                return "O(%s^%s)"%(self._parent.variable_name(),self.prec())
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
                    x = "(%s)"%x
                if e == 1:
                    var = "*%s"%X
                elif e == 0:
                    var = ""
                else:
                    var = "*%s^%s"%(X,e)
                s += "%s%s"%(x,var)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if self.prec() == 0:
            bigoh = "O(1)"
        elif self.prec() == 1:
            bigoh = "O(%s)"%self._parent.variable_name()
        else:
            bigoh = "O(%s^%s)"%(self._parent.variable_name(),self.prec())
        if self.prec() != infinity:
            if s == " ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def verschiebung(self, n):
        r"""
        Return the ``n``-th Verschiebung of ``self``.

        If `f = \sum a_m x^m` then this function returns `\sum a_m x^{mn}`.

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: f = -1/x + 1 + 2*x^2 + 5*x^5
            sage: f.V(2)
            -x^-2 + 1 + 2*x^4 + 5*x^10
            sage: f.V(-1)
            5*x^-5 + 2*x^-2 + 1 - x
            sage: h = f.add_bigoh(7)
            sage: h.V(2)
            -x^-2 + 1 + 2*x^4 + 5*x^10 + O(x^14)
            sage: h.V(-2)
            Traceback (most recent call last):
            ...
            ValueError: For finite precision only positive arguments allowed

        TESTS::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: f = x
            sage: f.V(3)
            x^3
            sage: f.V(-3)
            x^-3
            sage: g = 2*x^(-1) + 3 + 5*x
            sage: g.V(-1)
            5*x^-1 + 3 + 2*x
        """
        if n == 0:
            raise ValueError('n must be nonzero')

        if n < 0:
            if not self.prec() is infinity:
                raise ValueError('For finite precision only positive arguments allowed')

            exponents = [e * n for e in self.exponents()]
            u = min(exponents)
            exponents = [e - u for e in exponents]
            coefficients = self.coefficients()
            zero = self.base_ring().zero()
            w = [zero] * (max(exponents) + 1)
            for i in range(len(exponents)):
                e = exponents[i]
                c = coefficients[i]
                w[e] = c
            l = LaurentSeries(self._parent, w, u)
        else:
            __u = self.__u.V(n)
            __n = <long>self.__n * n
            l = LaurentSeries(self._parent, __u, __n)
        return l

    V = verschiebung

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = (17/2)*x^-2 + x + x^2 + 3*x^4 + O(x^7)
            sage: latex(f)
            \frac{\frac{17}{2}}{x^{2}} + x + x^{2} + 3x^{4} + O(x^{7})

        Verify that :issue:`6656` has been fixed::

            sage: R.<a,b> = PolynomialRing(QQ)
            sage: T.<x> = LaurentSeriesRing(R)
            sage: y = a*x + b*x
            sage: y._latex_()
            '\\left(a + b\\right)x'
            sage: latex(y)
            \left(a + b\right)x
        """
        if self.is_zero():
            if self.prec() is infinity:
                return "0"
            else:
                return "0 + \\cdots"
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
            x = sage.misc.latex.latex(x)
            if x != '0':
                if not first:
                    s += " + "
                if not atomic_repr and e > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if e == 1:
                    var = "|%s"%X
                elif e == 0:
                    var = ""
                elif e > 0:
                    var = "|%s^{%s}"%(X,e)
                if e >= 0:
                    s += "%s%s"%(x,var)
                else: # negative e
                    if e == -1:
                        s += "\\frac{%s}{%s}"%(x, X)
                    else:
                        s += "\\frac{%s}{%s^{%s}}"%(x, X,-e)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1|"," ")
        s = s.replace(" -1|", " -")
        s = s.replace("|","")
        pr = self.prec()
        if pr != infinity:
            if pr == 0:
                bigoh = "O(1)"
            elif pr == 1:
                bigoh = "O(%s)"%(X,)
            else:
                bigoh = "O(%s^{%s})"%(X,pr)
            if s == " ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def __hash__(self):
        return hash(self.__u) ^ self.__n

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
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
            sage: f = -5/t^(10) + 1/3 + t + t^2 - 10/3*t^3 + O(t^5); f
            -5*t^-10 + 1/3 + t + t^2 - 10/3*t^3 + O(t^5)

        Slicing can be used to truncate, keeping the same precision::

            sage: f[:2]
            -5*t^-10 + 1/3 + t + O(t^5)

        Any other kind of slicing is an error, see :issue:`18940`::

            sage: f[-10:2:2]
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a step is not defined

            sage: f[0:]
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a start is not defined
        """
        if isinstance(i, slice):
            start, stop, step = i.start, i.stop, i.step
            if step is not None:
                raise IndexError("polynomial slicing with a step is not defined")
            if start is not None:
                raise IndexError("polynomial slicing with a start is not defined")
            if stop > self.__u.degree() or stop is None:
                stop = self.__u.degree()
            f = self.__u[:stop - self.__n]
            return type(self)(self._parent, f, self.__n)

        return self.__u[i - self.__n]

    def __iter__(self):
        """
        Iterate through the coefficients from the first nonzero one to the
        last nonzero one.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
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

    def list(self):
        """
        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3
            sage: f.list()
            [-5, 0, 0, 1, 1, -10/3]
        """
        return self.__u.list()

    def coefficients(self):
        """
        Return the nonzero coefficients of ``self``.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3
            sage: f.coefficients()
            [-5, 1, 1, -10/3]
        """
        zero = self._parent.base_ring().zero()
        return [c for c in self.list() if c != zero]

    def residue(self):
        r"""
        Return the residue of ``self``.

        Consider the Laurent series

        .. MATH::

            f = \sum_{n \in \ZZ} a_n t^n
            = \cdots + \frac{a_{-2}}{t^2} + \frac{a_{-1}}{t} + a_0
            + a_1 t + a_2 t^2 + \cdots,

        then the residue of `f` is `a_{-1}`. Alternatively this is the
        coefficient of `1/t`.

        EXAMPLES::

            sage: t = LaurentSeriesRing(ZZ,'t').gen()
            sage: f = 1/t**2 + 2/t + 3 + 4*t
            sage: f.residue()
            2
            sage: f = t + t**2
            sage: f.residue()
            0
            sage: f.residue().parent()
            Integer Ring
        """
        return self[-1]

    def exponents(self):
        """
        Return the exponents appearing in ``self`` with nonzero coefficients.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3
            sage: f.exponents()
            [-2, 1, 2, 3]
        """
        zero = self._parent.base_ring().zero()
        v = self.valuation()
        return [i+v for i,val in enumerate(self.list()) if val != zero]

    def laurent_polynomial(self):
        """
        Return the corresponding Laurent polynomial.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = t^-3 + t + 7*t^2 + O(t^5)
            sage: g = f.laurent_polynomial(); g
            t^-3 + t + 7*t^2
            sage: g.parent()
            Univariate Laurent Polynomial Ring in t over Rational Field
        """
        R = self._parent.laurent_polynomial_ring()
        return R(self.__u.polynomial()) * R.gen()**(self.__n)

    def lift_to_precision(self, absprec=None):
        """
        Return a congruent Laurent series with absolute precision at least
        ``absprec``.

        INPUT:

        - ``absprec`` -- integer or ``None`` (default: ``None``); the
          absolute precision of the result. If ``None``, lifts to an exact
          element.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A.<t> = LaurentSeriesRing(GF(5))
            sage: x = t^(-1) + t^2 + O(t^5)
            sage: x.lift_to_precision(10)
            t^-1 + t^2 + O(t^10)
            sage: x.lift_to_precision()
            t^-1 + t^2
        """
        if absprec is not None and absprec <= self.precision_absolute():
            return self

        exact = self._parent(0) if self.is_zero() else self._parent(self.list()) << self.__n
        if absprec is None:
            return exact
        else:
            return exact.add_bigoh(absprec)

    def __setitem__(self, n, value):
        """
        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = t^2 + t^3 + O(t^10)
            sage: f[2] = 5
            Traceback (most recent call last):
            ...
            IndexError: Laurent series are immutable
        """
        raise IndexError("Laurent series are immutable")

    def _unsafe_mutate(self, i, value):
        """
        Never use this unless you really know what you are doing.

        .. WARNING::

           This could easily introduce subtle bugs, since Sage assumes
           everywhere that Laurent series are immutable. It's OK to use
           this if you really know what you're doing.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = t^2 + t^3 + O(t^10); f
            t^2 + t^3 + O(t^10)
            sage: f._unsafe_mutate(2, -3)
            sage: f
            -3*t^2 + t^3 + O(t^10)
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
        self._normalize()

    cpdef _add_(self, right_m):
        """
        Add two power series with the same parent.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: t + t
            2*t
            sage: f = 1/t + t^2 + t^3 - 17/3 * t^4 + O(t^5)
            sage: g = 1/(1-t + O(t^7)); g
            1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + O(t^7)
            sage: f + g
            t^-1 + 1 + t + 2*t^2 + 2*t^3 - 14/3*t^4 + O(t^5)
            sage: f + 0
            t^-1 + t^2 + t^3 - 17/3*t^4 + O(t^5)
            sage: 0 + f
            t^-1 + t^2 + t^3 - 17/3*t^4 + O(t^5)
            sage: R(0) + R(0)
            0
            sage: (t^3 + O(t^10)) + (t^-3 +O(t^9))
            t^-3 + t^3 + O(t^9)

        ALGORITHM: Shift the unit parts to align them, then add.

        TESTS:

        Verify that :issue:`35860` is fixed::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: sqrt(t^2) + t^-1
            t^-1 + t
        """
        cdef LaurentSeries right = <LaurentSeries>right_m
        cdef long m

        # 1. Special case when one or the other is 0.
        if not right:
            return self.add_bigoh(right.prec())
        if not self:
            return right.add_bigoh(self.prec())

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
        return type(self)(self._parent, f1 + f2, m)

    cpdef _sub_(self, right_m):
        """
        Subtract two power series with the same parent.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: t - t
            0
            sage: t^5 + 2 * t^-5
            2*t^-5 + t^5

        ALGORITHM: Shift the unit parts to align them, then subtract.
        """
        cdef LaurentSeries right = <LaurentSeries>right_m
        cdef long m

        # 1. Special case when one or the other is 0.
        if not right:
            return self.add_bigoh(right.prec())
        if not self:
            return -right.add_bigoh(self.prec())

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
        return type(self)(self._parent, f1 - f2, m)

    def add_bigoh(self, prec):
        """
        Return the truncated series at chosen precision ``prec``.

        The methods :meth:`add_bigoh` and :meth:`O` are synonyms.

        INPUT:

        - ``prec`` -- the precision of the series as an integer

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = t^2 + t^3 + O(t^10); f
            t^2 + t^3 + O(t^10)
            sage: f.add_bigoh(5)
            t^2 + t^3 + O(t^5)

        TESTS:

        Check that :issue:`28239` is fixed::

            sage: (t^(-2)).add_bigoh(-1)
            t^-2 + O(t^-1)
            sage: (t^(-2)).add_bigoh(-2)
            O(t^-2)
            sage: (t^(-2)).add_bigoh(-3)
            O(t^-3)
        """
        if prec is infinity or prec >= self.prec():
            return self
        P = self._parent
        if not self or prec < self.__n:
            return type(self)(P, P._power_series_ring(0, prec=0), prec)
        u = self.__u.add_bigoh(prec - self.__n)
        return type(self)(P, u, self.__n)

    O = add_bigoh

    def degree(self):
        """
        Return the degree of a polynomial equivalent to this power series
        modulo big oh of the precision.

        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: g = x^2 - x^4 + O(x^8)
            sage: g.degree()
            4
            sage: g = -10/x^5 + x^2 - x^4 + O(x^8)
            sage: g.degree()
            4
            sage: (x^-2 + O(x^0)).degree()
            -2
        """
        return self.__u.degree() + self.__n

    def __neg__(self):
        """
        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: -(1+t^5)
            -1 - t^5
            sage: -(1/(1+t+O(t^5)))
            -1 + t - t^2 + t^3 - t^4 + O(t^5)
        """
        return type(self)(self._parent, -self.__u, self.__n)

    cpdef _mul_(self, right_r):
        """
        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x^3 + x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1 - x + x^2 - x^4 + O(x^8)
            sage: f*g
            x^-3 - x^-2 + x^-1 + 4*x^4 + O(x^5)
        """
        cdef LaurentSeries right = <LaurentSeries>right_r
        return type(self)(self._parent,
                          self.__u * right.__u,
                          self.__n + right.__n)

    cpdef _rmul_(self, Element c):
        return type(self)(self._parent, self.__u._rmul_(c), self.__n)

    cpdef _lmul_(self, Element c):
        return type(self)(self._parent, self.__u._lmul_(c), self.__n)

    def __pow__(_self, r, dummy):
        """
        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: f^7
            x^7 + 7*x^8 + 21*x^9 + 56*x^10 + 161*x^11 + 336*x^12 + O(x^13)
            sage: g^7
            x^-70 - 7*x^-59 + 7*x^-58 - 7*x^-56 + O(x^-52)
            sage: g^(1/2)
            x^-5 - 1/2*x^6 + 1/2*x^7 - 1/2*x^9 + O(x^13)
            sage: g^(1/5)
            x^-2 - 1/5*x^9 + 1/5*x^10 - 1/5*x^12 + O(x^16)
            sage: g^(2/5)
            x^-4 - 2/5*x^7 + 2/5*x^8 - 2/5*x^10 + O(x^14)
            sage: h = x^2 + 2*x^4 + x^6
            sage: h^(1/2)
            x + x^3
        """
        cdef LaurentSeries self = _self

        try:
            right = QQ.coerce(r)
        except TypeError:
            raise ValueError("exponent must be a rational number")

        if right.denominator() == 1:
            right = right.numerator()
            return type(self)(self._parent, self.__u**right, self.__n*right)

        if self.is_zero():
            return self._parent(0).O(self.prec()*right)

        d = right.denominator()
        n = right.numerator()

        val = self.valuation()

        if val % d:
            raise ValueError("power series valuation would be fractional")

        u = self.valuation_zero_part().nth_root(d)

        s = type(self)(self._parent, u, val // d)

        return s**n

    def shift(self, k):
        r"""
        Return this Laurent series multiplied by the power `t^n`.
        Does not change this series.

        .. NOTE::

           Despite the fact that higher order terms are printed to the
           right in a power series, right shifting decreases the
           powers of `t`, while left shifting increases
           them. This is to be consistent with polynomials, integers,
           etc.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ['y'])
            sage: f = (t+t^-1)^4; f
            t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
            sage: f.shift(10)
            t^6 + 4*t^8 + 6*t^10 + 4*t^12 + t^14
            sage: f >> 10
            t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6
            sage: t << 4
            t^5
            sage: t + O(t^3) >> 4
            t^-3 + O(t^-1)

        AUTHORS:

        - Robert Bradshaw (2007-04-18)
        """
        return type(self)(self._parent, self.__u, self.__n + k)

    def __lshift__(LaurentSeries self, k):
        return type(self)(self._parent, self.__u, self.__n + k)

    def __rshift__(LaurentSeries self, k):
        return type(self)(self._parent, self.__u, self.__n - k)

    def truncate(self, long n):
        r"""
        Return the Laurent series of degree ` < n` which is
        equivalent to ``self`` modulo `x^n`.

        EXAMPLES::

            sage: A.<x> = LaurentSeriesRing(ZZ)
            sage: f = 1/(1-x)
            sage: f
            1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + x^11
             + x^12 + x^13 + x^14 + x^15 + x^16 + x^17 + x^18 + x^19 + O(x^20)
            sage: f.truncate(10)
            1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9
        """
        if n <= self.__n:
            return self._parent.zero()
        else:
            return type(self)(self._parent, self.__u.truncate(n - self.__n), self.__n)

    def truncate_laurentseries(self, long n):
        r"""
        Replace any terms of degree >= n by big oh.

        EXAMPLES::

            sage: A.<x> = LaurentSeriesRing(ZZ)
            sage: f = 1/(1-x)
            sage: f
            1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + x^11
             + x^12 + x^13 + x^14 + x^15 + x^16 + x^17 + x^18 + x^19 + O(x^20)
            sage: f.truncate_laurentseries(10)
            1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)
        """
        if n <= self.__n:
            return self._parent.zero()
        else:
            return type(self)(self._parent, self.__u.truncate_powerseries(n - self.__n), self.__n)

    def truncate_neg(self, long n):
        r"""
        Return the Laurent series equivalent to ``self`` except without any
        degree ``n`` terms.

        This is equivalent to::

            ``self - self.truncate(n)``

        EXAMPLES::

            sage: A.<t> = LaurentSeriesRing(ZZ)
            sage: f = 1/(1-t)
            sage: f.truncate_neg(15)
            t^15 + t^16 + t^17 + t^18 + t^19 + O(t^20)
        """
        return type(self)(self._parent, self.__u >> (n - self.__n), n)

    cpdef _div_(self, right_r):
        """
        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
            sage: f/x
            1 + x + 3*x^3 + O(x^6)
            sage: f/g
            x^8 + x^9 + 3*x^11 + O(x^14)
        """
        cdef LaurentSeries right = <LaurentSeries>right_r
        if right.__u.is_zero():
            raise ZeroDivisionError
        try:
            return type(self)(self._parent,
                              self.__u / right.__u,
                              self.__n - right.__n)
        except TypeError as msg:
            # todo: this could also make something in the formal fraction field.
            raise ArithmeticError("division not defined")

    def common_prec(self, other):
        r"""
        Return the minimum precision of ``self`` and ``other``.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)

        ::

            sage: f = t^(-1) + t + t^2 + O(t^3)
            sage: g = t + t^3 + t^4 + O(t^4)
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

        ::

            sage: f = t + t^2 + O(t^3)
            sage: g = t^(-3) + t^2
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

        ::

            sage: f = t + t^2
            sage: g = t^2
            sage: f.common_prec(g)
            +Infinity

        ::

            sage: f = t^(-3) + O(t^(-2))
            sage: g = t^(-5) + O(t^(-1))
            sage: f.common_prec(g)
            -2

        ::

            sage: f = O(t^2)
            sage: g = O(t^5)
            sage: f.common_prec(g)
            2
        """
        return min(self.prec(), other.prec())

    def common_valuation(self, other):
        """
        Return the minimum valuation of ``self`` and ``other``.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)

        ::

            sage: f = t^(-1) + t + t^2 + O(t^3)
            sage: g = t + t^3 + t^4 + O(t^4)
            sage: f.common_valuation(g)
            -1
            sage: g.common_valuation(f)
            -1

        ::

            sage: f = t + t^2 + O(t^3)
            sage: g = t^(-3) + t^2
            sage: f.common_valuation(g)
            -3
            sage: g.common_valuation(f)
            -3

        ::

            sage: f = t + t^2
            sage: g = t^2
            sage: f.common_valuation(g)
            1

        ::

            sage: f = t^(-3) + O(t^(-2))
            sage: g = t^(-5) + O(t^(-1))
            sage: f.common_valuation(g)
            -5

        ::

            sage: f = O(t^2)
            sage: g = O(t^5)
            sage: f.common_valuation(g)
            +Infinity
        """
        return min(self.valuation(), other.valuation())

    cpdef _richcmp_(self, right_r, int op):
        r"""
        Comparison of ``self`` and ``right``.

        We say two approximate Laurent series are equal, if they agree for
        all coefficients up to the *minimum* of the precisions of each.
        Comparison is done in dictionary order from lowest degree to
        highest degree coefficients. This is different than polynomials,
        but consistent with the idea that the variable of a Laurent
        series is considered to be "very small".

        See :meth:`power_series_ring_element._richcmp_` for more
        information.

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: f = x^(-1) + 1 + x + O(x^2)
            sage: g = x^(-1) + 1 + O(x)
            sage: f == g
            True

        ::

            sage: f = x^(-1) + 1 + x + O(x^2)
            sage: g = x^(-1) + 2 + O(x)
            sage: f == g
            False
            sage: f < g
            True
            sage: f > g
            False

        ::

            sage: f = x^(-2) + 1 + x + O(x^2)
            sage: g = x^(-1) + 2 + O(x)
            sage: f == g
            False
            sage: f < g
            False
            sage: f > g
            True

        Check that :issue:`19664` is fixed::

            sage: R.<x> = LaurentSeriesRing(RR)
            sage: x^(10^9) > 0
            True
        """
        cdef LaurentSeries right = <LaurentSeries>right_r

        val = self.common_valuation(right)
        if val is infinity:
            return rich_to_bool(op, 0)  # Both arguments are zero

        cdef long deg = max(self.degree(), right.degree())
        prec = self.common_prec(right)
        if deg >= prec:
            deg = prec - 1

        cdef long i
        for i in range(val, deg + 1):
            li = self[i]
            ri = right[i]
            if li != ri:
                return richcmp_not_equal(li, ri, op)
        return rich_to_bool(op, 0)

    def valuation_zero_part(self):
        """
        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: f/x
            1 + x + 3*x^3 + O(x^6)
            sage: f.valuation_zero_part()
            1 + x + 3*x^3 + O(x^6)
            sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
            sage: g.valuation_zero_part()
            1 - x^8 + x^9 - x^11 + O(x^15)
        """
        return self.__u

    def valuation(self):
        """
        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: f = 1/x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1 - x + x^2 - x^4 + O(x^8)
            sage: f.valuation()
            -1
            sage: g.valuation()
            0

        Note that the valuation of an element undistinguishable from
        zero is infinite::

            sage: h = f - f; h
            O(x^7)
            sage: h.valuation()
            +Infinity

        TESTS:

        The valuation of the zero element is ``+Infinity``
        (see :issue:`15088`)::

            sage: zero = R(0)
            sage: zero.valuation()
            +Infinity
        """
        if self.is_zero():
            return infinity
        return self.__n

    def variable(self):
        """
        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x + x^2 + 3*x^4 + O(x^7)
            sage: f.variable()
            'x'
        """
        return self._parent.variable_name()

    def prec(self):
        """
        This function returns the n so that the Laurent series is of the
        form (stuff) + `O(t^n)`. It doesn't matter how many
        negative powers appear in the expansion. In particular, prec could
        be negative.

        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f.prec()
            7
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: g.prec()
            8
        """
        return self.__u.prec() + self.__n

    def precision_absolute(self):
        """
        Return the absolute precision of this series.

        By definition, the absolute precision of
        `...+O(x^r)` is `r`.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: (t^2 + O(t^3)).precision_absolute()
            3
            sage: (1 - t^2 + O(t^100)).precision_absolute()
            100
        """
        return self.prec()

    def precision_relative(self):
        """
        Return the relative precision of this series, that
        is the difference between its absolute precision
        and its valuation.

        By convention, the relative precision of `0` (or
        `O(x^r)` for any `r`) is `0`.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: (t^2 + O(t^3)).precision_relative()
            1
            sage: (1 - t^2 + O(t^100)).precision_relative()
            100
            sage: O(t^4).precision_relative()
            0
        """
        if self.is_zero():
            return 0
        else:
            return self.prec() - self.valuation()

    def __copy__(self):
        return type(self)(self._parent, self.__u.__copy__(), self.__n)

    def reverse(self, precision=None):
        """
        Return the reverse of f, i.e., the series g such that g(f(x)) = x.
        Given an optional argument ``precision``, return the reverse with given
        precision (note that the reverse can have precision at most
        ``f.prec()``).  If ``f`` has infinite precision, and the argument
        ``precision`` is not given, then the precision of the reverse defaults
        to the default precision of ``f.parent()``.

        Note that this is only possible if the valuation of ``self`` is exactly
        1.

        The implementation depends on the underlying power series element
        implementing a reverse method.

        EXAMPLES::

            sage: R.<x> = Frac(QQ[['x']])
            sage: f = 2*x + 3*x^2 - x^4 + O(x^5)
            sage: g = f.reverse()
            sage: g
            1/2*x - 3/8*x^2 + 9/16*x^3 - 131/128*x^4 + O(x^5)
            sage: f(g)
            x + O(x^5)
            sage: g(f)
            x + O(x^5)

            sage: A.<t> = LaurentSeriesRing(ZZ)
            sage: a = t - t^2 - 2*t^4 + t^5 + O(t^6)
            sage: b = a.reverse(); b
            t + t^2 + 2*t^3 + 7*t^4 + 25*t^5 + O(t^6)
            sage: a(b)
            t + O(t^6)
            sage: b(a)
            t + O(t^6)

            sage: B.<b,c> = ZZ[ ]
            sage: A.<t> = LaurentSeriesRing(B)
            sage: f = t + b*t^2 + c*t^3 + O(t^4)
            sage: g = f.reverse(); g
            t - b*t^2 + (2*b^2 - c)*t^3 + O(t^4)
            sage: f(g)
            t + O(t^4)
            sage: g(f)
            t + O(t^4)

            sage: A.<t> = PowerSeriesRing(ZZ)
            sage: B.<s> = LaurentSeriesRing(A)
            sage: f = (1 - 3*t + 4*t^3 + O(t^4))*s + (2 + t + t^2 + O(t^3))*s^2 + O(s^3)
            sage: set_verbose(1)
            sage: g = f.reverse(); g
            verbose 1 (<module>) passing to pari failed; trying Lagrange inversion
            (1 + 3*t + 9*t^2 + 23*t^3 + O(t^4))*s + (-2 - 19*t - 118*t^2 + O(t^3))*s^2 + O(s^3)
            sage: set_verbose(0)
            sage: f(g) == g(f) == s
            True

        If the leading coefficient is not a unit, we pass to its fraction
        field if possible::

            sage: A.<t> = LaurentSeriesRing(ZZ)
            sage: a = 2*t - 4*t^2 + t^4 - t^5 + O(t^6)
            sage: a.reverse()
            1/2*t + 1/2*t^2 + t^3 + 79/32*t^4 + 437/64*t^5 + O(t^6)

            sage: B.<b> = PolynomialRing(ZZ)
            sage: A.<t> = LaurentSeriesRing(B)
            sage: f = 2*b*t + b*t^2 + 3*b^2*t^3 + O(t^4)
            sage: g = f.reverse(); g
            1/(2*b)*t - 1/(8*b^2)*t^2 + ((-3*b + 1)/(16*b^3))*t^3 + O(t^4)
            sage: f(g)
            t + O(t^4)
            sage: g(f)
            t + O(t^4)

        We can handle some base rings of positive characteristic::

            sage: A8.<t> = LaurentSeriesRing(Zmod(8))
            sage: a = t - 15*t^2 - 2*t^4 + t^5 + O(t^6)
            sage: b = a.reverse(); b
            t + 7*t^2 + 2*t^3 + 5*t^4 + t^5 + O(t^6)
            sage: a(b)
            t + O(t^6)
            sage: b(a)
            t + O(t^6)

        The optional argument ``precision`` sets the precision of the output::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: f = 2*x + 3*x^2 - 7*x^3 + x^4 + O(x^5)
            sage: g = f.reverse(precision=3); g
            1/2*x - 3/8*x^2 + O(x^3)
            sage: f(g)
            x + O(x^3)
            sage: g(f)
            x + O(x^3)

        If the input series has infinite precision, the precision of the
        output is automatically set to the default precision of the parent
        ring::

            sage: R.<x> = LaurentSeriesRing(QQ, default_prec=20)
            sage: (x - x^2).reverse()  # get some Catalan numbers
            x + x^2 + 2*x^3 + 5*x^4 + 14*x^5 + 42*x^6 + 132*x^7 + 429*x^8 + 1430*x^9
             + 4862*x^10 + 16796*x^11 + 58786*x^12 + 208012*x^13 + 742900*x^14
             + 2674440*x^15 + 9694845*x^16 + 35357670*x^17 + 129644790*x^18
             + 477638700*x^19 + O(x^20)
            sage: (x - x^2).reverse(precision=3)
            x + x^2 + O(x^3)

        TESTS::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: f = 1 + 2*x + 3*x^2 - x^4 + O(x^5)
            sage: f.reverse()
            Traceback (most recent call last):
            ...
            ValueError: Series must have valuation one for reversion.
        """
        val = self.valuation()
        if val != 1:
            raise ValueError("Series must have valuation one for reversion.")
        u = self.valuation_zero_part()
        u = u.parent().gen(0) * u

        rev = u.reverse(precision=precision)

        if rev.parent() == u.parent():
            return self._parent(rev)
        else:
            P = self._parent.change_ring(rev.parent().base_ring())
            return P(rev)

    def derivative(self, *args):
        """
        The formal derivative of this Laurent series, with respect to
        variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        .. SEEALSO::

           :meth:`_derivative`

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: g.derivative()
            -10*x^-11 - 1 + 2*x - 4*x^3 + O(x^7)
            sage: g.derivative(x)
            -10*x^-11 - 1 + 2*x - 4*x^3 + O(x^7)

        ::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = LaurentSeriesRing(R)
            sage: f = 2*t/x + (3*t^2 + 6*t)*x + O(x^2)
            sage: f.derivative()
            -2*t*x^-2 + (3*t^2 + 6*t) + O(x)
            sage: f.derivative(x)
            -2*t*x^-2 + (3*t^2 + 6*t) + O(x)
            sage: f.derivative(t)
            2*x^-1 + (6*t + 6)*x + O(x^2)
        """
        return multi_derivative(self, args)

    def _derivative(self, var=None):
        """
        The formal derivative of this Laurent series with respect to var.

        If var is None or the generator of this ring, it's the formal
        derivative as expected. Otherwise, _derivative(var) gets called
        recursively on each coefficient.

        .. SEEALSO::

           :meth:`derivative`

        EXAMPLES::

            sage: x = Frac(QQ[['x']]).0
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f._derivative()
            2*x + 12*x^3 + O(x^6)
            sage: f._derivative(x)
            2*x + 12*x^3 + O(x^6)
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: g._derivative()
            -10*x^-11 - 1 + 2*x - 4*x^3 + O(x^7)

        Differentiating with respect to something other than the generator
        gets recursed into the base ring::

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = LaurentSeriesRing(R)
            sage: f = 2*t/x + (3*t^2 + 6*t)*x + O(x^2)
            sage: f._derivative(t)
            2*x^-1 + (6*t + 6)*x + O(x^2)

        TESTS::

            sage: y = var('y')                                                          # needs sage.symbolic
            sage: f.derivative(y)                                                       # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var != self._parent.gen():
            try:
                # call _derivative() recursively on coefficients
                u = [coeff._derivative(var) for coeff in self.__u.list()]
                u = self._parent._power_series_ring(u, self.__u.prec())
                return type(self)(self._parent, u, self.__n)
            except AttributeError:
                raise ValueError("cannot differentiate with respect to {}".format(var))

        # compute formal derivative with respect to generator
        if self.is_zero():
            return type(self)(self._parent, 0, self.__u.prec() - 1)
        cdef long m, n = self.__n
        a = self.__u.list()
        v = [(n+m)*a[m] for m from 0 <= m < len(a)]
        u = self._parent._power_series_ring(v, self.__u.prec())
        return type(self)(self._parent, u, n-1)

    def integral(self):
        r"""
        The formal integral of this Laurent series with 0 constant term.

        EXAMPLES: The integral may or may not be defined if the base ring
        is not a field.

        ::

            sage: t = LaurentSeriesRing(ZZ, 't').0
            sage: f = 2*t^-3 + 3*t^2 + O(t^4)
            sage: f.integral()
            -t^-2 + t^3 + O(t^5)

        ::

            sage: f = t^3
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: Coefficients of integral cannot be coerced into the base ring

        The integral of 1/t is `\log(t)`, which is not given by a
        Laurent series::

            sage: t = Frac(QQ[['t']]).0
            sage: f = -1/t^3 - 31/t + O(t^3)
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: The integral of is not a Laurent series, since t^-1 has nonzero coefficient.

        Another example with just one negative coefficient::

            sage: A.<t> = QQ[[]]
            sage: f = -2*t^(-4) + O(t^8)
            sage: f.integral()
            2/3*t^-3 + O(t^9)
            sage: f.integral().derivative() == f
            True
        """
        cdef long i, n = self.__n
        a = self.__u.list()
        if self[-1] != 0:
            raise ArithmeticError("The integral of is not a Laurent series, since t^-1 has nonzero coefficient.")

        if n < 0:
            v = [a[i]/(n+i+1) for i in range(min(-1-n,len(a)))] + [0]
        else:
            v = []
        v += [a[i]/(n+i+1) for i in range(max(-n,0), len(a))]
        try:
            u = self._parent._power_series_ring(v, self.__u.prec())
        except TypeError:
            raise ArithmeticError("Coefficients of integral cannot be coerced into the base ring")
        return type(self)(self._parent, u, n+1)

    def nth_root(self, long n, prec=None):
        r"""
        Return the ``n``-th root of this Laurent power series.

        INPUT:

        - ``n`` -- integer

        - ``prec`` -- integer (optional); precision of the result. Though, if
          this series has finite precision, then the result cannot have larger
          precision.

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: (x^-2 + 1 + x).nth_root(2)
            x^-1 + 1/2*x + 1/2*x^2 - ... - 19437/65536*x^18 + O(x^19)
            sage: (x^-2 + 1 + x).nth_root(2)**2
            x^-2 + 1 + x + O(x^18)

            sage: # needs sage.modular
            sage: j = j_invariant_qexp()
            sage: q = j.parent().gen()
            sage: j(q^3).nth_root(3)
            q^-1 + 248*q^2 + 4124*q^5 + ... + O(q^29)
            sage: (j(q^2) - 1728).nth_root(2)
            q^-1 - 492*q - 22590*q^3 - ... + O(q^19)
        """
        if prec is None:
            prec = self.prec()
            if prec is infinity:
                prec = self.parent().default_prec()
        else:
            prec = min(self.prec(), prec)

        if n <= 0:
            raise ValueError('n must be positive')

        i = self.valuation()
        if i % n:
            raise ValueError('valuation must be divisible by n')

        q = self.__u.nth_root(n, prec)
        return type(self)(self._parent, q + self.parent()(0).O(prec), i // n)

    def power_series(self):
        """
        Convert this Laurent series to a power series.

        An error is raised if the Laurent series has a term (or an error
        term `O(x^k)`) whose exponent is negative.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: f = 1/(1-t+O(t^10)); f.parent()
            Laurent Series Ring in t over Integer Ring
            sage: g = f.power_series(); g
            1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)
            sage: parent(g)
            Power Series Ring in t over Integer Ring
            sage: f = 3/t^2 +  t^2 + t^3 + O(t^10)
            sage: f.power_series()
            Traceback (most recent call last):
            ...
            TypeError: self is not a power series

        TESTS:

        Check whether a polynomial over a Laurent series ring is contained in the
        polynomial ring over the power series ring (see :issue:`19459`):

            sage: # needs sage.rings.finite_rings
            sage: L.<t> = LaurentSeriesRing(GF(2))
            sage: R.<x,y> = PolynomialRing(L)
            sage: S.<x,y> = PolynomialRing(L._power_series_ring)
            sage: t**(-1)*x*y in S
            False

        There used to be an issue with non-canonical representations of zero,
        see :issue:`31383`::

            sage: S.<x> = PowerSeriesRing(QQ)
            sage: L = Frac(S)
            sage: s = L(O(x^2))
            sage: (s*x^(-1)).power_series()
            O(x^1)
            sage: (s*x^(-2)).power_series()
            O(x^0)
            sage: (s*x^(-3)).power_series()
            Traceback (most recent call last):
            ...
            TypeError: self is not a power series

        Test for :issue:`32440`::

            sage: L.<x> = LaurentSeriesRing(QQ, implementation='pari')                  # needs sage.libs.pari
            sage: (x + O(x^3)).power_series()                                           # needs sage.libs.pari
            x + O(x^3)
        """
        if self.__n < 0:
            if self.__u.is_zero() and self.__u.prec() >= - self.__n:
                return self.__u >> (- self.__n)
            else:
                raise TypeError("self is not a power series")
        return self.__u << self.__n

    def inverse(self):
        """
        Return the inverse of self, i.e., self^(-1).

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: t.inverse()
            t^-1
            sage: (1-t).inverse()
            1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + ...
        """
        return ~self

    def __call__(self, *x, **kwds):
        """
        Compute value of this Laurent series at x.

        EXAMPLES::

            sage: P.<x, y> = ZZ[]
            sage: R.<t> = LaurentSeriesRing(P)
            sage: f = x*t^-2 + y*t^2 + O(t^8)
            sage: f(t^3)
            x*t^-6 + y*t^6 + O(t^24)
            sage: f(t=t^3)
            x*t^-6 + y*t^6 + O(t^24)
            sage: f(t + O(t^5))
            x*t^-2 + O(t^2)
            sage: f(y=x)
            x*t^-2 + x*t^2 + O(t^8)
            sage: f(t^3, x=2, y=x + x^2)
            2*t^-6 + (x^2 + x)*t^6 + O(t^24)
            sage: f(t^3, 2, x+x^2)
            2*t^-6 + (x^2 + x)*t^6 + O(t^24)
            sage: f(x=2, t=t^3, y=x + x^2)
            2*t^-6 + (x^2 + x)*t^6 + O(t^24)
            sage: f(2, x + x^2, t=t^3)
            Traceback (most recent call last):
            ...
            ValueError: must not specify t keyword and positional argument

        It is only possible to substitute elements of positive valuation::

            sage: f(t^-2)
            Traceback (most recent call last):
            ...
            ValueError: Can only substitute elements of positive valuation

        Test for :issue:`23928`::

            sage: R.<x> = LaurentSeriesRing(QQ, implementation='pari')                  # needs sage.libs.pari
            sage: f = x.add_bigoh(7)                                                    # needs sage.libs.pari
            sage: f(x)                                                                  # needs sage.libs.pari
            x + O(x^7)
            """
        if len(kwds) >= 1:
            name = self.parent().variable_name()
            if name in kwds: # a keyword specifies the Laurent series generator
                if x:
                    raise ValueError("must not specify %s keyword and positional argument" % name)
                a = self(kwds[name])
                del kwds[name]
                try:
                    return a(**kwds)
                except TypeError:
                    return a
            elif x:       # both keywords and positional arguments
                a = self(*x)
                try:
                    return a(**kwds)
                except TypeError:
                    return a
            else:                  # keywords but no positional arguments
                return self.__u(**kwds)*(self.parent().gen()**self.__n)

        if len(x) == 0:
            return self

        if isinstance(x[0], tuple):
            x = x[0]

        return self.__u(*x)*(x[0]**self.__n)

    def __pari__(self):
        """
        Convert ``self`` to a PARI object.

        TESTS::

            sage: L.<x> = LaurentSeriesRing(QQ)
            sage: f = x + 1/x + O(x^2); f
            x^-1 + x + O(x^2)
            sage: f.__pari__()                                                          # needs sage.libs.pari
            x^-1 + x + O(x^2)

        Check that :issue:`32437` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: F.<u> = GF(257^2)
            sage: R.<t> = LaurentSeriesRing(F)
            sage: g = t + O(t^99)
            sage: f = u*t + O(t^99)
            sage: g(f)  # indirect doctest                                              # needs sage.libs.pari
            u*t + O(t^99)
        """
        f = self.__u
        x = f.parent().gen()
        return f.__pari__() * x.__pari__()**self.__n
