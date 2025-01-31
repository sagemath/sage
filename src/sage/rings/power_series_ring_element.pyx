"""
Power Series

Sage provides an implementation of dense and sparse power series
over any Sage base ring. This is the base class of the implementations
of univariate and multivariate power series ring elements in Sage
(see also :doc:`power_series_poly`, :doc:`multi_power_series_ring_element`).

AUTHORS:

- William Stein
- David Harvey (2006-09-11): added solve_linear_de() method
- Robert Bradshaw (2007-04): sqrt, rmul, lmul, shifting
- Robert Bradshaw (2007-04): Cython version
- Simon King (2012-08): use category and coercion framework, :issue:`13412`

EXAMPLES::

    sage: R.<x> = PowerSeriesRing(ZZ)
    sage: TestSuite(R).run()
    sage: R([1,2,3])
    1 + 2*x + 3*x^2
    sage: R([1,2,3], 10)
    1 + 2*x + 3*x^2 + O(x^10)
    sage: f = 1 + 2*x - 3*x^3 + O(x^4); f
    1 + 2*x - 3*x^3 + O(x^4)
    sage: f^10
    1 + 20*x + 180*x^2 + 930*x^3 + O(x^4)
    sage: g = 1/f; g
    1 - 2*x + 4*x^2 - 5*x^3 + O(x^4)
    sage: g * f
    1 + O(x^4)

In Python (as opposed to Sage) create the power series ring and
its generator as follows::

    sage: R = PowerSeriesRing(ZZ, 'x')
    sage: x = R.gen()
    sage: parent(x)
    Power Series Ring in x over Integer Ring

EXAMPLES:

This example illustrates that coercion for power
series rings is consistent with coercion for polynomial rings.

::

    sage: poly_ring1.<gen1> = PolynomialRing(QQ)
    sage: poly_ring2.<gen2> = PolynomialRing(QQ)
    sage: huge_ring.<x> = PolynomialRing(poly_ring1)

The generator of the first ring gets coerced in as itself, since it
is the base ring.

::

    sage: huge_ring(gen1)
    gen1

The generator of the second ring gets mapped via the natural map
sending one generator to the other.

::

    sage: huge_ring(gen2)
    x

With power series the behavior is the same.

::

    sage: power_ring1.<gen1> = PowerSeriesRing(QQ)
    sage: power_ring2.<gen2> = PowerSeriesRing(QQ)
    sage: huge_power_ring.<x> = PowerSeriesRing(power_ring1)
    sage: huge_power_ring(gen1)
    gen1
    sage: huge_power_ring(gen2)
    x
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

from cpython.object cimport Py_EQ, Py_NE
from sage.rings.infinity import infinity, InfinityElement

from sage.rings.rational_field import QQ

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
import sage.misc.misc
import sage.arith.all as arith
import sage.misc.latex
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

from sage.categories.fields import Fields
_Fields = Fields()

from sage.misc.derivative import multi_derivative
from sage.structure.element cimport AlgebraElement, RingElement
from sage.structure.richcmp cimport richcmp


def is_PowerSeries(x):
    """
    Return ``True`` if ``x`` is an instance of a univariate
    or multivariate power series.

    EXAMPLES::

        sage: R.<x> = PowerSeriesRing(ZZ)
        sage: from sage.rings.power_series_ring_element import is_PowerSeries
        sage: is_PowerSeries(1 + x^2)
        doctest:warning...
        DeprecationWarning: The function is_PowerSeries is deprecated; use 'isinstance(..., PowerSeries)' instead.
        See https://github.com/sagemath/sage/issues/38266 for details.
        True
        sage: is_PowerSeries(x - x)
        True
        sage: is_PowerSeries(0)
        False
        sage: var('x')                                                                  # needs sage.symbolic
        x
        sage: is_PowerSeries(1 + x^2)                                                   # needs sage.symbolic
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38266,
                       "The function is_PowerSeries is deprecated; "
                       "use 'isinstance(..., PowerSeries)' instead.")
    return isinstance(x, PowerSeries)


cdef class PowerSeries(AlgebraElement):
    """
    A power series. Base class of univariate and
    multivariate power series. The following methods
    are available with both types of objects.
    """
    def __init__(self, parent, prec, is_gen=False):
        """
        Initialize a power series. Not for public use.
        It gets called by the ``PowerSeries_poly`` and
        ``MPowerSeries`` constructors.

        EXAMPLES::

            sage: PowerSeriesRing(CC, 'q')                                              # needs sage.rings.real_mpfr
            Power Series Ring in q over Complex Field with 53 bits of precision
            sage: T = PowerSeriesRing(GF(3), 5, 't'); T
            Multivariate Power Series Ring in t0, t1, t2, t3, t4
             over Finite Field of size 3
        """
        AlgebraElement.__init__(self, parent)
        self._is_gen = is_gen
        self._prec = prec

    def __hash__(self):
        """
        Compute a hash of ``self``.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: (1+x^2).__hash__()     # random
            15360174650385709
        """
        return hash(self.polynomial())

    def __reduce__(self):
        """
        EXAMPLES::

            sage: K.<t> = PowerSeriesRing(QQ, 5)
            sage: f = 1 + t - t^3 + O(t^12)
            sage: loads(dumps(f)) == f
            True
        """
        return make_element_from_parent_v0, (self._parent, self.polynomial(), self.prec())

    def is_sparse(self):
        """
        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t.is_sparse()
            False
            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: t.is_sparse()
            True
        """
        return self._parent.is_sparse()

    def is_dense(self):
        """
        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t.is_dense()
            True
            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: t.is_dense()
            False
        """
        return self._parent.is_dense()

    def is_gen(self):
        """
        Return ``True`` if this is the generator (the variable) of the power
        series ring.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: t.is_gen()
            True
            sage: (1 + 2*t).is_gen()
            False

        Note that this only returns ``True`` on the actual generator, not on
        something that happens to be equal to it.

        ::

            sage: (1*t).is_gen()
            False
            sage: 1*t == t
            True
        """
        return bool(self._is_gen)

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of this series under the map that sends the
        generators to ``im_gens``. This is used internally for computing
        homomorphisms.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: f = 1 + t + t^2
            sage: f._im_gens_(ZZ, [3])
            13
        """
        if base_map is None:
            return codomain(self(im_gens[0]))
        else:
            raise NotImplementedError

    cpdef base_extend(self, R):
        """
        Return a copy of this power series but with coefficients in R.

        The following coercion uses ``base_extend`` implicitly::

            sage: R.<t> = ZZ[['t']]
            sage: (t - t^2) * Mod(1, 3)
            t + 2*t^2
        """
        S = self._parent.base_extend(R)
        return S(self)

    def change_ring(self, R):
        """
        Change if possible the coefficients of ``self`` to lie in R.

        EXAMPLES::

            sage: R.<T> = QQ[[]]; R
            Power Series Ring in T over Rational Field
            sage: f = 1 - 1/2*T + 1/3*T^2 + O(T^3)
            sage: f.base_extend(GF(5))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
            sage: f.change_ring(GF(5))
            1 + 2*T + 2*T^2 + O(T^3)
            sage: f.change_ring(GF(3))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(0, 3) does not exist

        We can only change the ring if there is a ``__call__`` coercion
        defined. The following succeeds because ``ZZ(K(4))`` is defined.

        ::

            sage: K.<a> = NumberField(cyclotomic_polynomial(3), 'a')                    # needs sage.rings.number_field
            sage: R.<t> = K[['t']]                                                      # needs sage.rings.number_field
            sage: (4*t).change_ring(ZZ)                                                 # needs sage.rings.number_field
            4*t

        This does not succeed because ``ZZ(K(a+1))`` is not defined.

        ::

            sage: K.<a> = NumberField(cyclotomic_polynomial(3), 'a')                    # needs sage.rings.number_field
            sage: R.<t> = K[['t']]                                                      # needs sage.rings.number_field
            sage: ((a+1)*t).change_ring(ZZ)                                             # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce a + 1 to an integer
        """
        S = self._parent.change_ring(R)
        return S(self)

    cpdef _richcmp_(self, right, int op):
        r"""
        Comparison of ``self`` and ``right``.

        We say two approximate power series are equal if they agree for
        all coefficients up to the *minimum* of the precisions of each.
        Thus, e.g., `f = 1 + q + O(q^2)` is equal to `g = 1 + O(q)`.

        This is how PARI defines equality of power series, but not how
        Magma defines equality. (Magma would declare `f` and `g` unequal.)
        The PARI/Sage convention is consistent with the idea that
        comparison should take place after coercing both elements into
        a common parent.  Hence, in the above example `f` is truncated
        to `f + O(q)`, which is considered to be equal to `g`, even
        though the coefficients of `q` are unknown for both series in
        that comparison.

        Comparison is done in dictionary order from lowest degree to
        highest degree coefficients.  This is different than polynomial
        comparison.

        EXAMPLES::

            sage: R.<q> = ZZ[[ ]]; R
            Power Series Ring in q over Integer Ring
            sage: f=1+q+O(q^2); g = 1+O(q)
            sage: f == g
            True
            sage: 1 - 2*q + q^2 +O(q^3) == 1 - 2*q^2 + q^2 + O(q^4)
            False

        ::

            sage: R.<t> = ZZ[[]]
            sage: 1 + t^2 < 2 - t
            True
            sage: f = 1 + t + t^7 - 5*t^10
            sage: g = 1 + t + t^7 - 5*t^10 + O(t^15)
            sage: f == f
            True
            sage: f < g
            False
            sage: f == g
            True

        TESTS:

        :issue:`9457` is fixed::

            sage: A.<t> = PowerSeriesRing(ZZ)
            sage: g = t + t^3 + t^5 + O(t^6); g
            t + t^3 + t^5 + O(t^6)
            sage: [g == g.add_bigoh(i) for i in range(7)]
            [True, True, True, True, True, True, True]
            sage: A(g.polynomial()) == g
            True

            sage: f = t + t^2 + O(t^10)
            sage: f == f.truncate()
            True
        """
        if op == Py_EQ:
            return not (self-right)
        if op == Py_NE:
            return bool(self-right)
        prec = self.common_prec(right)
        x = self.list()
        y = right.list()
        if not (prec is infinity):
            x += [0]*(prec - len(x)) # self.list() does not include trailing zeroes
            x = x[:prec] # truncate x to common prec
            y += [0]*(prec - len(y))
            y = y[:prec]
        return richcmp(x, y, op)

    def __call__(self, x):
        """
        Implementations *MUST* override this in the derived class.

        EXAMPLES::

            sage: from sage.rings.power_series_ring_element import PowerSeries
            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: PowerSeries.__call__(1+x^2,x)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def coefficients(self):
        """
        Return the nonzero coefficients of ``self``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t + t^2 - 10/3*t^3
            sage: f.coefficients()
            [1, 1, -10/3]
        """
        zero = self.parent().base_ring().zero()
        return [c for c in self.list() if c != zero]

    def exponents(self):
        """
        Return the exponents appearing in ``self`` with nonzero coefficients.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t + t^2 - 10/3*t^3
            sage: f.exponents()
            [1, 2, 3]
        """
        zero = self.parent().base_ring().zero()
        l = self.list()
        return [i for i in range(len(l)) if l[i] != zero]

    def list(self):
        """
        See this method in derived classes:

        - :meth:`sage.rings.power_series_poly.PowerSeries_poly.list`,

        - :meth:`sage.rings.multi_power_series_ring_element.MPowerSeries.list`

        Implementations *MUST* override this in the derived class.

        EXAMPLES::

            sage: from sage.rings.power_series_ring_element import PowerSeries
            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: PowerSeries.list(1 + x^2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def polynomial(self):
        """
        See this method in derived classes:

        - :meth:`sage.rings.power_series_poly.PowerSeries_poly.polynomial`,

        - :meth:`sage.rings.multi_power_series_ring_element.MPowerSeries.polynomial`

        Implementations *MUST* override this in the derived class.

        EXAMPLES::

            sage: from sage.rings.power_series_ring_element import PowerSeries
            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: PowerSeries.polynomial(1 + x^2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def lift_to_precision(self, absprec=None):
        """
        Return a congruent power series with absolute precision at least
        ``absprec``.

        INPUT:

        - ``absprec`` -- integer or ``None`` (default: ``None``); the
          absolute precision of the result. If ``None``, lifts to an exact
          element.

        EXAMPLES::

            sage: A.<t> = PowerSeriesRing(GF(5))
            sage: x = t + t^2 + O(t^5)
            sage: x.lift_to_precision(10)
            t + t^2 + O(t^10)
            sage: x.lift_to_precision()
            t + t^2
        """
        if absprec is not None and absprec <= self.precision_absolute():
            return self

        exact = self._parent(self.list())
        if absprec is None:
            return exact
        else:
            return exact.add_bigoh(absprec)

    def __copy__(self):
        """
        Return this power series. Power series are immutable so copy can
        safely just return the same polynomial.

        EXAMPLES::

            sage: R.<q> = ZZ[[ ]]; R
            Power Series Ring in q over Integer Ring
            sage: f = 1 + 3*q + O(q^10)
            sage: copy(f) is f       # !!! ok since power series are immutable.
            True
        """
        return self

    def base_ring(self):
        """
        Return the base ring that this power series is defined over.

        EXAMPLES::

            sage: R.<t> = GF(49,'alpha')[[]]                                            # needs sage.rings.finite_rings
            sage: (t^2 + O(t^3)).base_ring()                                            # needs sage.rings.finite_rings
            Finite Field in alpha of size 7^2
        """
        return self._parent.base_ring()

    def padded_list(self, n=None):
        """
        Return a list of coefficients of ``self`` up to (but not including)
        `q^n`.

        Includes 0s in the list on the right so that the list has length `n`.

        INPUT:

        - ``n`` -- (optional) an integer that is at least 0. If ``n`` is
          not given, it will be taken to be the precision of ``self``,
          unless this is ``+Infinity``, in which case we just return
          ``self.list()``.

        EXAMPLES::

            sage: R.<q> = PowerSeriesRing(QQ)
            sage: f = 1 - 17*q + 13*q^2 + 10*q^4 + O(q^7)
            sage: f.list()
            [1, -17, 13, 0, 10]
            sage: f.padded_list(7)
            [1, -17, 13, 0, 10, 0, 0]
            sage: f.padded_list(10)
            [1, -17, 13, 0, 10, 0, 0, 0, 0, 0]
            sage: f.padded_list(3)
            [1, -17, 13]
            sage: f.padded_list()
            [1, -17, 13, 0, 10, 0, 0]
            sage: g = 1 - 17*q + 13*q^2 + 10*q^4
            sage: g.list()
            [1, -17, 13, 0, 10]
            sage: g.padded_list()
            [1, -17, 13, 0, 10]
            sage: g.padded_list(10)
            [1, -17, 13, 0, 10, 0, 0, 0, 0, 0]
        """
        if n is None:
            if self._prec is infinity:
                return self.list()
            else:
                n = self._prec
        if n < 0:
            raise ValueError("n must be at least 0")
        v = self.list()
        if len(v) < n:
            z = self._parent.base_ring()(0)
            return v + [z]*(n - len(v))
        else:
            return v[:int(n)]

    def prec(self):
        """
        The precision of `...+O(x^r)` is by definition
        `r`.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: (t^2 + O(t^3)).prec()
            3
            sage: (1 - t^2 + O(t^100)).prec()
            100
        """
        return self._prec

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

    def _repr_(self):
        """
        Return the string representation of this power series.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: (t^2 + O(t^3))._repr_()
            't^2 + O(t^3)'

        ::

            sage: R.<t> = QQ[[]]
            sage: 1 / (1+2*t +O(t^5))
            1 - 2*t + 4*t^2 - 8*t^3 + 16*t^4 + O(t^5)

        ::

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: 1 / (1+2*t +O(t^5))
            1 - 2*t + 4*t^2 - 8*t^3 + 16*t^4 + O(t^5)
            sage: -13/2 * t^3  + 5*t^5 + O(t^10)
            -13/2*t^3 + 5*t^5 + O(t^10)
        """
        if self.is_zero():
            if self.prec() is infinity:
                return "0"
            else:
                return "O(%s^%s)"%(self._parent.variable_name(),self.prec())

        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        X = self._parent.variable_name()

        s = " "
        if self.is_sparse():
            f = self.polynomial()
            m = f.degree() + 1
            d = f._dict_unsafe()
            coeffs = sorted(d.items())
            for (n, x) in coeffs:
                x = repr(x)
                if x != '0':
                    if s != ' ':
                        s += " + "
                    if not atomic_repr and n > 0 and (x.find("+") != -1 or x.find("-") != -1):
                        x = "(%s)"%x
                    if n > 1:
                        var = "*%s^%s"%(X,n)
                    elif n==1:
                        var = "*%s"%X
                    else:
                        var = ""
                    s += "%s%s"%(x,var)
        else:
            v = self.list()
            m = len(v)
            first = True
            for n in range(m):
                x = v[n]
                x = repr(x)
                if x != '0':
                    if not first:
                        s += " + "
                    if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                        x = "(%s)"%x
                    if n > 1:
                        var = "*%s^%s"%(X,n)
                    elif n==1:
                        var = "*%s"%X
                    else:
                        var = ""
                    s += "%s%s"%(x,var)
                    first = False
        # end

        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if not (self._prec is infinity):
            if self._prec == 0:
                bigoh = "O(1)"
            elif self._prec == 1:
                bigoh = "O(%s)"%self._parent.variable_name()
            else:
                bigoh = "O(%s^%s)"%(self._parent.variable_name(),self._prec)
            if s==" ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def _latex_(self):
        r"""
        Return the latex representation of this power series.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: f = -1/2 * t + 2/3*t^2 - 9/7 * t^15 + O(t^20); f
            -1/2*t + 2/3*t^2 - 9/7*t^15 + O(t^20)
            sage: latex(f)
            -\frac{1}{2} t + \frac{2}{3} t^{2} - \frac{9}{7} t^{15} + O(t^{20})

        Check that :issue:`26606` is fixed::

            sage: R.<beta> = QQ[]
            sage: S.<x> = R[[]]
            sage: latex(beta*x)
            \beta x
        """
        if self.is_zero():
            if self.prec() is infinity:
                return "0"
            else:
                return "0 + \\cdots"
        s = " "
        v = self.list()
        m = len(v)
        X = self._parent.latex_variable_names()[0]
        atomic_repr = self._parent.base_ring()._repr_option('element_is_atomic')
        first = True
        for n in range(m):
            x = v[n]
            x = sage.misc.latex.latex(x)
            if x != '0':
                if not first:
                    s += " + "
                if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if n > 1:
                    var = "%s^{%s}"%(X,n)
                elif n==1:
                    var = "%s"%X
                else:
                    var = ""
                if n > 0:
                    s += "%s| %s"%(x,var)
                else:
                    s += repr(x)
                first = False

        s = s.replace(" + -", " - ")
        s = s.replace(" -1|", " -")
        s = s.replace(" 1|"," ")
        s = s.replace("|","")
        if not (self._prec is infinity):
            if self._prec == 0:
                bigoh = "O(1)"
            elif self._prec == 1:
                bigoh = "O(%s)"%(X,)
            else:
                bigoh = "O(%s^{%s})"%(X,self._prec)
            if s == " ":
                return bigoh
            s += " + %s"%bigoh
        return s.lstrip(" ")

    def truncate(self, prec=infinity):
        """
        The polynomial obtained from power series by truncation.

        EXAMPLES::

            sage: R.<I> = GF(2)[[]]
            sage: f = 1/(1+I+O(I^8)); f
            1 + I + I^2 + I^3 + I^4 + I^5 + I^6 + I^7 + O(I^8)
            sage: f.truncate(5)
            I^4 + I^3 + I^2 + I + 1
        """
        if prec is infinity:
            prec = self._prec
        a = self.list()
        v = [a[i] for i in range(min(prec, len(a)))]
        return self._parent._poly_ring()(v)

    cdef _inplace_truncate(self, long prec):
        return self.truncate(prec)

    def add_bigoh(self, prec):
        r"""
        Return the power series of precision at most ``prec`` got by adding
        `O(q^\text{prec})` to `f`, where `q` is the variable.

        EXAMPLES::

            sage: R.<A> = RDF[[]]
            sage: f = (1+A+O(A^5))^5; f
            1.0 + 5.0*A + 10.0*A^2 + 10.0*A^3 + 5.0*A^4 + O(A^5)
            sage: f.add_bigoh(3)
            1.0 + 5.0*A + 10.0*A^2 + O(A^3)
            sage: f.add_bigoh(5)
            1.0 + 5.0*A + 10.0*A^2 + 10.0*A^3 + 5.0*A^4 + O(A^5)
        """
        if prec is infinity or prec > self.prec():
            return self
        a = self.list()
        v = [a[i] for i in range(min(prec, len(a)))]
        return self._parent(v, prec)

    def __getitem__(self, n):
        r"""
        Return the coefficient of `t^n` in this power series, where
        `t` is the indeterminate of the power series ring.

        If `n` is negative return 0. If `n` is beyond the precision, raise an
        :exc:`IndexError`.

        EXAMPLES::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: R.<m> = CDF[[]]
            sage: f = CDF(pi)^2 + m^3 + CDF(e)*m^4 + O(m^10); f   # abs tol 5e-16
            9.869604401089358 + 0.0*m + 0.0*m^2 + 1.0*m^3 + 2.718281828459045*m^4 + O(m^10)
            sage: f[-5]
            0.0
            sage: f[0]
            9.869604401089358
            sage: f[4]   # abs tol 5e-16
            2.718281828459045
            sage: f[9]
            0.0
            sage: f[10]
            Traceback (most recent call last):
            ...
            IndexError: coefficient not known
            sage: f[1000]
            Traceback (most recent call last):
            ...
            IndexError: coefficient not known
        """
        if n<0:
            return self.base_ring()(0)
        c = self.list()
        if n >= len(c):
            if self._prec > n:
                return self.base_ring()(0)
            else:
                raise IndexError("coefficient not known")
        return c[n]

    def common_prec(self, f):
        r"""
        Return minimum precision of `f` and ``self``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)

        ::

            sage: f = t + t^2 + O(t^3)
            sage: g = t + t^3 + t^4 + O(t^4)
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

        ::

            sage: f = t + t^2 + O(t^3)
            sage: g = t^2
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

        ::

            sage: f = t + t^2
            sage: f = t^2
            sage: f.common_prec(g)
            +Infinity
        """
        if self.prec() is infinity:
            return f.prec()
        elif f.prec() is infinity:
            return self.prec()
        return min(self.prec(), f.prec())

    cdef common_prec_c(self, PowerSeries f):
        if self._prec is infinity:
            return f._prec
        elif f._prec is infinity:
            return self._prec
        elif self._prec < f._prec:
            return self._prec
        else:
            return f._prec

    def _mul_prec(self, RingElement right_r):
        cdef PowerSeries right = <PowerSeries>right_r
        sp = self._prec
        rp = right._prec
        if sp is infinity:
            if rp is infinity:
                prec = infinity
            else:
                prec = rp + self.valuation()
        else:  # sp != infinity
            if rp is infinity:
                prec = sp + right.valuation()
            else:
                prec = min(rp + self.valuation(), sp + right.valuation())
        # endif
        return prec

    def __bool__(self):
        """
        Return ``True`` if this power series is not equal to 0.

        EXAMPLES::

            sage: R.<q> = ZZ[[ ]]; R
            Power Series Ring in q over Integer Ring
            sage: f = 1 + 3*q + O(q^10)
            sage: f.is_zero()
            False
            sage: (0 + O(q^2)).is_zero()
            True
            sage: R(0).is_zero()
            True
            sage: (0 + O(q^1000)).is_zero()
            True
        """
        return not not self.polynomial()

    def is_unit(self):
        """
        Return ``True`` if this power series is invertible.

        A power series is invertible precisely when the
        constant term is invertible.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: (-1 + t - t^5).is_unit()
            True
            sage: (3 + t - t^5).is_unit()
            False
            sage: O(t^0).is_unit()
            False

        AUTHORS:

        - David Harvey (2006-09-03)
        """
        # Return False for O(t^0) (the only power-series-ring element with precision_absolute == 0).
        if self.precision_absolute() == 0:
            return False
        return self[0].is_unit()

    def inverse(self):
        """
        Return the inverse of self, i.e., self^(-1).

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: t.inverse()
            t^-1
            sage: type(_)
            <class 'sage.rings.laurent_series_ring_element.LaurentSeries'>
            sage: (1-t).inverse()
            1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + ...
        """
        return ~self

    def valuation_zero_part(self):
        r"""
        Factor ``self`` as `q^n \cdot (a_0 + a_1 q + \cdots)` with
        `a_0` nonzero. Then this function returns
        `a_0 + a_1 q + \cdots` .

        .. NOTE::

           This valuation zero part need not be a unit if, e.g.,
           `a_0` is not invertible in the base ring.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: ((1/3)*t^5*(17-2/3*t^3)).valuation_zero_part()
            17/3 - 2/9*t^3

        In this example the valuation 0 part is not a unit::

            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: u = (-2*t^5*(17-t^3)).valuation_zero_part(); u
            -34 + 2*t^3
            sage: u.is_unit()
            False
            sage: u.valuation()
            0
        """
        if self.is_zero():
            raise ValueError("power series has no valuation 0 part")
        n = self.valuation()
        if n == 0:
            return self
        elif self.is_dense():
            v = self.list()[int(n):]
        else:
            n = int(n)
            v = {}
            for k, x in self.monomial_coefficients().items():
                if k >= n:
                    v[k-n] = x
        return self._parent(v, self.prec()-n)

    cpdef _div_(self, denom_r):
        """
        EXAMPLES::

            sage: k.<t> = QQ[[]]
            sage: t/t
            1
            sage: (t/(t^3 + 1)) * (t^3 + 1)
            t + O(t^21)
            sage: (t^5/(t^2 - 2)) * (t^2 -2 )
            t^5 + O(t^25)

        TESTS:

        The following tests against bugs that were fixed in :issue:`8972`::

            sage: P.<t> = ZZ[]
            sage: R.<x> = P[[]]
            sage: 1/(t*x)
            1/t*x^-1
            sage: R.<x> = ZZ[[]]
            sage: (1/x).parent()
            Laurent Series Ring in x over Rational Field
            sage: F = FractionField(R)
            sage: 1/x in F
            True
            sage: (1/(2*x)).parent()
            Laurent Series Ring in x over Rational Field
        """
        denom = <PowerSeries>denom_r
        if denom.is_zero():
            raise ZeroDivisionError("Can't divide by something indistinguishable from 0")
        u = denom.valuation_zero_part()
        inv = ~u  # inverse

        v = denom.valuation()
        if v > self.valuation():
            try:
                R = self._parent.fraction_field()
            except (TypeError, NotImplementedError):  # no integral domain
                R = self._parent.laurent_series_ring()
            return R(self) / R(denom)

        # Algorithm: Cancel common factors of q from top and bottom,
        # then invert the denominator.  We do the cancellation first
        # because we can only invert a unit (and remain in the ring
        # of power series).

        if v > 0:
            num = self >> v
        else:
            num = self
        return num*inv

    def __mod__(self, other):
        """
        EXAMPLES::

            sage: R.<T> = Qp(7)[[]]                                                     # needs sage.rings.padics
            sage: f = (48*67 + 46*67^2)*T + (1 + 42*67 + 5*67^3)*T^2 + O(T^3)           # needs sage.rings.padics
            sage: f % 67                                                                # needs sage.rings.padics
            T^2 + O(T^3)
        """
        from sage.rings.power_series_ring import PowerSeriesRing
        if isinstance(other, (int, Integer)):
            return PowerSeriesRing(IntegerModRing(other), self.variable())(self)
        raise NotImplementedError("Mod on power series ring elements not defined except modulo an integer.")

    def __pow__(self, r, dummy):
        """
        EXAMPLES::

            sage: x = QQ[['x']].0
            sage: f = x^2 + x^4 + O(x^6)
            sage: f^(1/2)
            x + 1/2*x^3 + O(x^5)
            sage: f^7
            x^14 + 7*x^16 + O(x^18)
            sage: f^(7/2)
            x^7 + 7/2*x^9 + O(x^11)
            sage: h = x^2 + 2*x^4 + x^6
            sage: h^(1/2)
            x + x^3
            sage: O(x^4)^(1/2)
            O(x^2)
        """
        try:
            right = QQ.coerce(r)
        except TypeError:
            raise ValueError("exponent must be a rational number")

        if right.denominator() == 1:
            right = right.numerator()
            return super().__pow__(right, dummy)

        if self.is_zero():
            return self.parent()(0).O((self.prec()*right).floor())

        d = right.denominator()
        n = right.numerator()

        return self.nth_root(d)**n

    def shift(self, n):
        r"""
        Return this power series multiplied by the power `t^n`.

        If `n` is negative, terms below `t^{-n}` are discarded.

        This power series is left unchanged.

        .. NOTE::

           Despite the fact that higher order terms are printed to the
           right in a power series, right shifting decreases the
           powers of `t`, while left shifting increases them.
           This is to be consistent with polynomials, integers, etc.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ['y'], 't', 5)
            sage: f = ~(1+t); f
            1 - t + t^2 - t^3 + t^4 + O(t^5)
            sage: f.shift(3)
            t^3 - t^4 + t^5 - t^6 + t^7 + O(t^8)
            sage: f >> 2
            1 - t + t^2 + O(t^3)
            sage: f << 10
            t^10 - t^11 + t^12 - t^13 + t^14 + O(t^15)
            sage: t << 29
            t^30

        AUTHORS:

        - Robert Bradshaw (2007-04-18)
        """
        if not n:
            return self
        prec = max(0, self.prec() + Integer(n))
        return self._parent(self.polynomial().shift(n), prec)

    def __lshift__(self, n):
        """
        Left-shift this power series by `n`, i.e., multiply by `t^n`.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = exp(x) + O(x^7); f
            1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + 1/720*x^6 + O(x^7)
            sage: f << 2
            x^2 + x^3 + 1/2*x^4 + 1/6*x^5 + 1/24*x^6 + 1/120*x^7 + 1/720*x^8 + O(x^9)
            sage: (f << 99) >> 99
            1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + 1/720*x^6 + O(x^7)
        """
        return self.shift(n)

    def __rshift__(self, n):
        """
        Right-shift this power series by `n`, i.e., divide by `t^n`.

        Terms below `t^n` are discarded.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = exp(x) + O(x^7)
            sage: f >> 3
            1/6 + 1/24*x + 1/120*x^2 + 1/720*x^3 + O(x^4)
            sage: f >> 7
            O(x^0)
            sage: f >> 99
            O(x^0)
            sage: (f >> 99) << 99
            O(x^99)
        """
        return self.shift(-n)

    def is_monomial(self):
        """
        Return ``True`` if this element is a monomial.  That is, if ``self`` is
        `x^n` for some nonnegative integer `n`.

        EXAMPLES::

            sage: k.<z> = PowerSeriesRing(QQ, 'z')
            sage: z.is_monomial()
            True
            sage: k(1).is_monomial()
            True
            sage: (z+1).is_monomial()
            False
            sage: (z^2909).is_monomial()
            True
            sage: (3*z^2909).is_monomial()
            False
        """

        return self.polynomial().is_monomial()

    def map_coefficients(self, f, new_base_ring=None):
        r"""
        Return the series obtained by applying ``f`` to the nonzero
        coefficients of ``self``.

        If ``f`` is a :class:`sage.categories.map.Map`, then the resulting
        series will be defined over the codomain of ``f``. Otherwise, the
        resulting polynomial will be over the same ring as ``self``. Set
        ``new_base_ring`` to override this behaviour.

        INPUT:

        - ``f`` -- a callable that will be applied to the coefficients of ``self``

        - ``new_base_ring`` -- (optional) if given, the resulting polynomial
          will be defined over this ring

        EXAMPLES::

            sage: R.<x> = SR[[]]                                                        # needs sage.symbolic
            sage: f = (1+I)*x^2 + 3*x - I                                               # needs sage.symbolic
            sage: f.map_coefficients(lambda z: z.conjugate())                           # needs sage.symbolic
            I + 3*x + (-I + 1)*x^2
            sage: R.<x> = ZZ[[]]
            sage: f = x^2 + 2
            sage: f.map_coefficients(lambda a: a + 42)
            44 + 43*x^2

        Examples with different base ring::

            sage: R.<x> = ZZ[[]]
            sage: k = GF(2)
            sage: residue = lambda x: k(x)
            sage: f = 4*x^2+x+3
            sage: g = f.map_coefficients(residue); g
            1 + x
            sage: g.parent()
            Power Series Ring in x over Integer Ring
            sage: g = f.map_coefficients(residue, new_base_ring=k); g
            1 + x
            sage: g.parent()
            Power Series Ring in x over Finite Field of size 2
            sage: residue = k.coerce_map_from(ZZ)
            sage: g = f.map_coefficients(residue); g
            1 + x
            sage: g.parent()
            Power Series Ring in x over Finite Field of size 2

        Tests other implementations::

            sage: # needs sage.libs.pari
            sage: R.<q> = PowerSeriesRing(GF(11), implementation='pari')
            sage: f = q - q^3 + O(q^10)
            sage: f.map_coefficients(lambda c: c - 2)
            10*q + 8*q^3 + O(q^10)
        """
        pol = self.polynomial()
        res = pol.map_coefficients(f, new_base_ring)
        if res.base_ring() != pol.base_ring():
            return self.parent().change_ring(res.base_ring())(res, self.prec())
        else:
            return self.parent()(res, self.prec())

    def jacobi_continued_fraction(self):
        r"""
        Return the Jacobi continued fraction of ``self``.

        The J-fraction or Jacobi continued fraction of a power series
        is a continued fraction expansion with steps of size two. We use
        the following convention

        .. MATH::

            1 / (1 + A_0 t + B_0 t^2 / (1 + A_1 t + B_1 t^2 / (1 + \cdots)))

        OUTPUT:

        tuple of pairs `(A_n, B_n)` for `n \geq 0`

        The expansion is done as long as possible given the precision.
        Whenever the expansion is not well-defined, because it would
        require to divide by zero, an exception is raised.

        See section 2.7 of [Kra1999det]_ for the close relationship
        of this kind of expansion with Hankel determinants and
        orthogonal polynomials.

        EXAMPLES::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: s = sum(factorial(k) * t**k for k in range(12)).O(12)
            sage: s.jacobi_continued_fraction()
            ((-1, -1), (-3, -4), (-5, -9), (-7, -16), (-9, -25))

        Another example::

            sage: (log(1+t)/t).jacobi_continued_fraction()
            ((1/2, -1/12),
             (1/2, -1/15),
             (1/2, -9/140),
             (1/2, -4/63),
             (1/2, -25/396),
             (1/2, -9/143),
             (1/2, -49/780),
             (1/2, -16/255),
             (1/2, -81/1292))

        TESTS::

             sage: (t).jacobi_continued_fraction()
             Traceback (most recent call last):
             ...
             ValueError: vanishing constant term, no expansion
             sage: (1/(1+3*t)).jacobi_continued_fraction()
             Traceback (most recent call last):
             ...
             ValueError: vanishing term, no further expansion
        """
        t = self.parent().gen()
        if self[0] == 0:
            raise ValueError('vanishing constant term, no expansion')
        serie = self / self[0]
        resu = []
        while serie.prec() >= 3:
            u = serie.inverse()
            A, B = u[1], u[2]
            resu.append((A, B))
            if B == 0:
                raise ValueError('vanishing term, no further expansion')
            serie = (u - 1 - A * t) / (B * t ** 2)
        return tuple(resu)

    def super_fraction(self, delta):
        """
        Return the super delta continued fraction of ``self``.

        This is a fraction of the following shape::

                                    v0 x^k0
             --------------------------------------------------------
                                        v1 x^(k0 + k1 + delta)
              U1(x)   -   -------------------------------------------
                                               v2 x^(k1 +k2 + delta)
                                U2(x)    -  -------------------------
                                                     U3(x) - ...

        where each `U_j(x) = 1 + u_j(x)x`.

        INPUT:

        - ``delta`` -- positive integer

        OUTPUT:

        The list of  ( vj, kj, U(j+1) ), j \geq 0.

        REFERENCES:

        - [Han2016]_

        EXAMPLES::

            sage: deg = 30
            sage: PS = PowerSeriesRing(QQ, 'q', default_prec=deg+1)
            sage: q = PS.gen()
            sage: F = prod([(1+q**k).add_bigoh(deg+1) for k in range(1,deg)])
            sage: F.super_fraction(2)
            [(1, 0, -q + 1),
             (1, 1, q + 1),
             (-1, 2, -q^3 + q^2 - q + 1),
             (1, 1, q^2 + q + 1),
             (-1, 0, -q + 1),
             (-1, 1, q^2 + q + 1),
             (-1, 0, -q + 1),
             (1, 1, 3*q^2 + 2*q + 1),
             (-4, 0, -q + 1)]

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: s = sum(factorial(k) * t**k for k in range(12)).O(12)
            sage: s.super_fraction(2)
            [(1, 0, -t + 1),
            (1, 0, -3*t + 1),
            (4, 0, -5*t + 1),
            (9, 0, -7*t + 1),
            (16, 0, -9*t + 1),
            (25, 0, -11*t + 1)]
        """
        ring = self.parent()
        q = ring.gen()
        Gi, Gj = ring.one(), self

        deg = self.prec()

        vkU = []

        di = Gi.valuation()
        ci = Gi[di]

        while True:
            dj = Gj.valuation()
            cj = Gj[dj]
            k, v = dj - di, cj / ci
            c = v * q**k
            gi = Gi.add_bigoh(dj + delta)
            gj = Gj.add_bigoh(k + dj + delta)
            U = (c * gi / gj).truncate()
            Gk = (U * Gj - Gi * c) >> (k + delta)
            deg = deg - 2 * k - delta
            if deg < 0:
                break
            vkU.append((v, k, U))
            if deg == 0:
                break
            if Gk.degree() == -1:
                break
            di, ci, Gi, Gj = dj, cj, Gj, Gk

        return vkU

    def stieltjes_continued_fraction(self):
        r"""
        Return the Stieltjes continued fraction of ``self``.

        The S-fraction or Stieltjes continued fraction of a power series
        is a continued fraction expansion with steps of size one. We use
        the following convention

        .. MATH::

            1 / (1 - A_1 t / (1 - A_2 t / (1 - A_3 t / (1 - \cdots))))

        OUTPUT: `A_n` for `n \geq 1`

        The expansion is done as long as possible given the precision.
        Whenever the expansion is not well-defined, because it would
        require to divide by zero, an exception is raised.

        EXAMPLES::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: s = sum(catalan_number(k) * t**k for k in range(12)).O(12)
            sage: s.stieltjes_continued_fraction()
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

        Another example::

            sage: (exp(t)).stieltjes_continued_fraction()
            (1,
             -1/2,
             1/6,
             -1/6,
             1/10,
             -1/10,
             1/14,
             -1/14,
             1/18,
             -1/18,
             1/22,
             -1/22,
             1/26,
             -1/26,
             1/30,
             -1/30,
             1/34,
             -1/34,
             1/38)

        TESTS::

             sage: (t).stieltjes_continued_fraction()
             Traceback (most recent call last):
             ...
             ValueError: vanishing constant term, no expansion
             sage: (1/(1+3*t)).stieltjes_continued_fraction()
             Traceback (most recent call last):
             ...
             ValueError: vanishing term, no further expansion
        """
        t = self.parent().gen()
        if self[0] == 0:
            raise ValueError('vanishing constant term, no expansion')
        serie = self / self[0]
        resu = []
        while serie.prec() >= 2:
            u = serie.inverse()
            A = -u[1]
            resu.append(A)
            if A == 0:
                raise ValueError('vanishing term, no further expansion')
            serie = (u - 1) / (-A * t)
        return tuple(resu)

    def is_square(self):
        """
        Return ``True`` if this function has a square root in this ring, e.g.,
        there is an element `y` in ``self.parent()``
        such that `y^2` equals ``self``.

        ALGORITHM: If the base ring is a field, this is true whenever the
        power series has even valuation and the leading coefficient is a
        perfect square.

        For an integral domain, it attempts the square root in the
        fraction field and tests whether or not the result lies in the
        original ring.

        EXAMPLES::

            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: (1+t).is_square()
            True
            sage: (2+t).is_square()
            False
            sage: (2+t.change_ring(RR)).is_square()
            True
            sage: t.is_square()
            False
            sage: K.<t> = PowerSeriesRing(ZZ, 't', 5)
            sage: (1+t).is_square()
            False
            sage: f = (1+t)^100
            sage: f.is_square()
            True
        """
        val = self.valuation()
        if val is not infinity and val % 2 == 1:
            return False
        elif not self[val].is_square():
            return False
        elif self.base_ring() in _Fields:
            return True
        else:
            try:
                self.parent()(self.sqrt())
                return True
            except TypeError:
                return False

    def sqrt(self, prec=None, extend=False, all=False, name=None):
        r"""
        Return a square root of ``self``.

        INPUT:

          - ``prec`` -- integer (default: ``None``); if not ``None`` and the
            series has infinite precision, truncates series at precision ``prec``

          - ``extend`` -- boolean (default: ``False``); if ``True``, return a
            square root in an extension ring, if necessary. Otherwise, raise
            a :exc:`ValueError` if the square root is not in the base power series
            ring. For example, if ``extend`` is ``True``, the square root of a
            power series with odd degree leading coefficient is
            defined as an element of a formal extension ring.

          - ``name`` -- string; if ``extend`` is ``True``, you must also
            specify the print name of the formal square root

          - ``all`` -- boolean (default: ``False``); if ``True``, return all
            square roots of ``self``, instead of just one

        ALGORITHM: Newton's method

        .. MATH::

           x_{i+1} = \frac{1}{2}( x_i + \mathrm{self}/x_i )

        EXAMPLES::

            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: sqrt(t^2)
            t
            sage: sqrt(1 + t)
            1 + 1/2*t - 1/8*t^2 + 1/16*t^3 - 5/128*t^4 + O(t^5)
            sage: sqrt(4 + t)
            2 + 1/4*t - 1/64*t^2 + 1/512*t^3 - 5/16384*t^4 + O(t^5)
            sage: u = sqrt(2 + t, prec=2, extend=True, name = 'alpha'); u
            alpha
            sage: u^2
            2 + t
            sage: u.parent()
            Univariate Quotient Polynomial Ring in alpha
             over Power Series Ring in t over Rational Field
             with modulus x^2 - 2 - t
            sage: K.<t> = PowerSeriesRing(QQ, 't', 50)
            sage: sqrt(1 + 2*t + t^2)
            1 + t
            sage: sqrt(t^2 + 2*t^4 + t^6)
            t + t^3
            sage: sqrt(1 + t + t^2 + 7*t^3)^2
            1 + t + t^2 + 7*t^3 + O(t^50)
            sage: sqrt(K(0))
            0
            sage: sqrt(t^2)
            t

        ::

            sage: # needs sage.rings.complex_double
            sage: K.<t> = PowerSeriesRing(CDF, 5)
            sage: v = sqrt(-1 + t + t^3, all=True); v
            [1.0*I - 0.5*I*t - 0.125*I*t^2 - 0.5625*I*t^3 - 0.2890625*I*t^4 + O(t^5),
             -1.0*I + 0.5*I*t + 0.125*I*t^2 + 0.5625*I*t^3 + 0.2890625*I*t^4 + O(t^5)]
            sage: [a^2 for a in v]
            [-1.0 + 1.0*t + 0.0*t^2 + 1.0*t^3 + O(t^5), -1.0 + 1.0*t + 0.0*t^2 + 1.0*t^3 + O(t^5)]

        A formal square root::

            sage: K.<t> = PowerSeriesRing(QQ, 5)
            sage: f = 2*t + t^3 + O(t^4)
            sage: s = f.sqrt(extend=True, name='sqrtf'); s
            sqrtf
            sage: s^2
            2*t + t^3 + O(t^4)
            sage: parent(s)
            Univariate Quotient Polynomial Ring in sqrtf
             over Power Series Ring in t over Rational Field
             with modulus x^2 - 2*t - t^3 + O(t^4)

        TESTS::

            sage: R.<x> = QQ[[]]
            sage: (x^10/2).sqrt()
            Traceback (most recent call last):
            ...
            ValueError: unable to take the square root of 1/2

        Check :issue:`30655`::

            sage: t = polygen(QQ, 't')
            sage: x = t.parent()[['x']].0
            sage: W = (t*x + 1 - x).O(3)
            sage: W.sqrt()
            1 + (1/2*t - 1/2)*x + (-1/8*t^2 + 1/4*t - 1/8)*x^2 + O(x^3)

        AUTHORS:

        - Robert Bradshaw

        - William Stein
        """
        if self.is_zero():
            ans = self._parent(0).O(self.prec()/2)
            if all:
                return [ans]
            else:
                return ans

        if all and not self.base_ring().is_integral_domain():
            raise NotImplementedError('all roots not implemented over a non-integral domain')

        formal_sqrt = False
        u = self.valuation_zero_part()
        # TODO, fix underlying element sqrt()
        try:
            try:
                s = u[0].sqrt(extend=False)
            except TypeError:
                s = u[0].sqrt()
        except ValueError:
            formal_sqrt = True
        if self.degree() == 0:
            if not formal_sqrt:
                a = self.parent()([s], self.prec())
                if all:
                    return [a, -a]
                else:
                    return a

        val = self.valuation()

        if formal_sqrt or val % 2:
            if extend:
                if name is None:
                    raise ValueError("the square root generates an extension, so you must specify the name of the square root")
                R = self._parent['x']
                S = R.quotient(R([-self, 0, 1]), names=name)
                a = S.gen()
                if all:
                    if not self.base_ring().is_integral_domain():
                        raise NotImplementedError('all roots not implemented over a non-integral domain')
                    return [a, -a]
                else:
                    return a
            elif formal_sqrt:
                raise ValueError("unable to take the square root of %s" % u[0])
            else:
                raise ValueError("power series does not have a square root since it has odd valuation.")

        pr = self.prec()
        if pr == infinity:
            test_exact = True
            if prec is None:
                pr = self._parent.default_prec()
            else:
                pr = prec
        else:
            test_exact = False
        prec = pr

        R = s.parent()
        a = self.valuation_zero_part()
        P = self._parent
        if not P.base_ring().has_coerce_map_from(R):
            a = a.change_ring(R)
        half = ~R(2)

        s = a.parent()([s])
        for cur_prec in sage.misc.misc.newton_method_sizes(prec)[1:]:
            (<PowerSeries>s)._prec = cur_prec
            s = half * (s + a/s)

        ans = s
        if val != 0:
            ans *= P.gen(0) ** (val // 2)
        if test_exact and ans.degree() < prec/2:
            if ans*ans == self:
                (<PowerSeries>ans)._prec = infinity

        if all:
            return [ans, -ans]  # since over an integral domain
        else:
            return ans

    def square_root(self):
        """
        Return the square root of ``self`` in this ring. If this cannot be done,
        then an error will be raised.

        This function succeeds if and only if
        ``self``. :meth:`.is_square`

        EXAMPLES::

            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: (1 + t).square_root()
            1 + 1/2*t - 1/8*t^2 + 1/16*t^3 - 5/128*t^4 + O(t^5)
            sage: (2 + t).square_root()
            Traceback (most recent call last):
            ...
            ValueError: Square root does not live in this ring.
            sage: (2 + t.change_ring(RR)).square_root()                                 # needs sage.rings.real_mpfr
            1.41421356237309 + 0.353553390593274*t - 0.0441941738241592*t^2
             + 0.0110485434560398*t^3 - 0.00345266983001244*t^4 + O(t^5)
            sage: t.square_root()
            Traceback (most recent call last):
            ...
            ValueError: Square root not defined for power series of odd valuation.
            sage: K.<t> = PowerSeriesRing(ZZ, 't', 5)
            sage: f = (1+t)^20
            sage: f.square_root()
            1 + 10*t + 45*t^2 + 120*t^3 + 210*t^4 + O(t^5)
            sage: f = 1 + t
            sage: f.square_root()
            Traceback (most recent call last):
            ...
            ValueError: Square root does not live in this ring.

        AUTHORS:

        - Robert Bradshaw
        """
        val = self.valuation()
        if val is not infinity and val % 2 == 1:
            raise ValueError("Square root not defined for power series of odd valuation.")
        elif not self[val].is_square():
            raise ValueError("Square root does not live in this ring.")
        elif self.base_ring() in _Fields:
            return self.sqrt()
        else:
            try:
                return self.parent()(self.sqrt())
            except TypeError:
                raise ValueError("Square root does not live in this ring.")

    def nth_root(self, n, prec=None):
        r"""
        Return the ``n``-th root of this power series.

        INPUT:

        - ``n`` -- integer

        - ``prec`` -- integer (optional); precision of the result. Though, if
          this series has finite precision, then the result cannot have larger
          precision.

        EXAMPLES::

            sage: R.<x> = QQ[[]]
            sage: (1+x).nth_root(5)
            1 + 1/5*x - 2/25*x^2 + ... + 12039376311816/2384185791015625*x^19 + O(x^20)

            sage: (1 + x + O(x^5)).nth_root(5)
            1 + 1/5*x - 2/25*x^2 + 6/125*x^3 - 21/625*x^4 + O(x^5)

        Check that the results are consistent with taking log and exponential::

            sage: R.<x> = PowerSeriesRing(QQ, default_prec=100)
            sage: p = (1 + 2*x - x^4)**200
            sage: p1 = p.nth_root(1000, prec=100)
            sage: p2 = (p.log()/1000).exp()
            sage: p1.prec() == p2.prec() == 100
            True
            sage: p1.polynomial() == p2.polynomial()
            True

        Positive characteristic::

            sage: R.<u> = GF(3)[[]]
            sage: p = 1 + 2 * u^2
            sage: p.nth_root(4)
            1 + 2*u^2 + u^6 + 2*u^8 + u^12 + 2*u^14 + O(u^20)
            sage: p.nth_root(4)**4
            1 + 2*u^2 + O(u^20)

        TESTS:

        Check that exact roots show infinite precision::

            sage: ((1+x)^5).nth_root(5)
            1 + x

        Check precision on `O(x^r)`::

            sage: O(x^4).nth_root(2)
            O(x^2)
            sage: O(x^4).nth_root(3)
            O(x^1)
            sage: O(x^4).nth_root(4)
            O(x^1)

        Check precision on higher valuation series::

            sage: (x^5+x^6+O(x^7)).nth_root(5)
            x + 1/5*x^2 + O(x^3)
        """

        val = self.valuation()

        if self.is_zero():
            if val is infinity:
                return self
            else:
                return self.parent()(0).O(val // n)

        if val is not infinity and val % n != 0:
            raise ValueError("power series valuation is not a multiple of %s" % n)

        maxprec = (val // n) + self.precision_relative()

        if prec is None:
            prec = maxprec
            if prec == infinity:
                prec = self.parent().default_prec()
        else:
            prec = min(maxprec, prec)

        p = self.polynomial()
        q = p._nth_root_series(n, prec)
        ans = self.parent()(q)
        if not (self.prec() == infinity and q.degree() * n <= prec and q**n == p):
            ans = ans.add_bigoh(prec)
        return ans

    def cos(self, prec=infinity):
        r"""
        Apply cos to the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT: a new power series

        EXAMPLES:

        For one variable::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: f = (t + t**2).O(4)
            sage: cos(f)
            1 - 1/2*t^2 - t^3 + O(t^4)

        For several variables::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(3)
            sage: cos(f)
            1 - 1/2*a^2 - a*b - 1/2*b^2 + O(a, b)^3
            sage: f.cos()
            1 - 1/2*a^2 - a*b - 1/2*b^2 + O(a, b)^3
            sage: f.cos(prec=2)
            1 + O(a, b)^2

        If the power series has a nonzero constant coefficient `c`,
        one raises an error::

            sage: g = 2+f
            sage: cos(g)
            Traceback (most recent call last):
            ...
            ValueError: can only apply cos to formal power series with zero constant term

        If no precision is specified, the default precision is used::

            sage: T.default_prec()
            12
            sage: cos(a)
            1 - 1/2*a^2 + 1/24*a^4 - 1/720*a^6 + 1/40320*a^8 - 1/3628800*a^10 + O(a, b)^12
            sage: a.cos(prec=5)
            1 - 1/2*a^2 + 1/24*a^4 + O(a, b)^5
            sage: cos(a + T.O(5))
            1 - 1/2*a^2 + 1/24*a^4 + O(a, b)^5

        TESTS::

            sage: cos(a^2 + T.O(5))
            1 - 1/2*a^4 + O(a, b)^5
        """
        R = self.parent()

        c = self[0]
        if not c.is_zero():
            raise ValueError('can only apply cos to formal power '
                             'series with zero constant term')
        x = self
        val = x.valuation()
        assert(val >= 1)

        prec = min(prec, self.prec())
        if isinstance(prec, InfinityElement):
            prec = R.default_prec()
        n_inv_factorial = R.base_ring().one()
        x_pow_n = R.one()
        x2 = x ** 2
        cos_x = R.one().add_bigoh(prec)
        for n in range(2, prec // val + 1, 2):
            x_pow_n = (x_pow_n * x2).add_bigoh(prec)
            n_inv_factorial /= - n * (n - 1)
            cos_x += x_pow_n * n_inv_factorial
        result_bg = cos_x

        if result_bg.base_ring() is not self.base_ring():
            R = R.change_ring(self.base_ring().fraction_field())
        return R(result_bg, prec=prec)

    def sin(self, prec=infinity):
        r"""
        Apply sin to the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT: a new power series

        EXAMPLES:

        For one variable::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: f = (t + t**2).O(4)
            sage: sin(f)
            t + t^2 - 1/6*t^3 + O(t^4)

        For several variables::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(3)
            sage: sin(f)
            a + b + a*b + O(a, b)^3
            sage: f.sin()
            a + b + a*b + O(a, b)^3
            sage: f.sin(prec=2)
            a + b + O(a, b)^2

        If the power series has a nonzero constant coefficient `c`,
        one raises an error::

            sage: g = 2+f
            sage: sin(g)
            Traceback (most recent call last):
            ...
            ValueError: can only apply sin to formal power series with zero constant term

        If no precision is specified, the default precision is used::

            sage: T.default_prec()
            12
            sage: sin(a)
            a - 1/6*a^3 + 1/120*a^5 - 1/5040*a^7 + 1/362880*a^9 - 1/39916800*a^11 + O(a, b)^12
            sage: a.sin(prec=5)
            a - 1/6*a^3 + O(a, b)^5
            sage: sin(a + T.O(5))
            a - 1/6*a^3 + O(a, b)^5

        TESTS::

            sage: sin(a^2 + T.O(5))
            a^2 + O(a, b)^5
        """
        R = self.parent()

        c = self[0]
        if not c.is_zero():
            raise ValueError('can only apply sin to formal power '
                             'series with zero constant term')
        val = self.valuation()
        assert(val >= 1)

        x = self

        prec = min(prec, self.prec())
        if isinstance(prec, InfinityElement):
            prec = R.default_prec()
        n_inv_factorial = R.base_ring().one()
        x_pow_n = x
        x2 = x ** 2
        sin_x = x.add_bigoh(prec)
        for n in range(3, prec // val + 1, 2):
            x_pow_n = (x_pow_n * x2).add_bigoh(prec)
            n_inv_factorial /= - n * (n - 1)
            sin_x += x_pow_n * n_inv_factorial
        result_bg = sin_x

        if result_bg.base_ring() is not self.base_ring():
            R = R.change_ring(self.base_ring().fraction_field())
        return R(result_bg, prec=prec)

    def tan(self, prec=infinity):
        r"""
        Apply tan to the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT: a new power series

        EXAMPLES:

        For one variable::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: f = (t + t**2).O(4)
            sage: tan(f)
            t + t^2 + 1/3*t^3 + O(t^4)

        For several variables::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(3)
            sage: tan(f)
            a + b + a*b + O(a, b)^3
            sage: f.tan()
            a + b + a*b + O(a, b)^3
            sage: f.tan(prec=2)
            a + b + O(a, b)^2

        If the power series has a nonzero constant coefficient `c`,
        one raises an error::

            sage: g = 2 + f
            sage: tan(g)
            Traceback (most recent call last):
            ...
            ValueError: can only apply tan to formal power series with zero constant term

        If no precision is specified, the default precision is used::

            sage: T.default_prec()
            12
            sage: tan(a)
            a + 1/3*a^3 + 2/15*a^5 + 17/315*a^7 + 62/2835*a^9 + 1382/155925*a^11 + O(a, b)^12
            sage: a.tan(prec=5)
            a + 1/3*a^3 + O(a, b)^5
            sage: tan(a + T.O(5))
            a + 1/3*a^3 + O(a, b)^5

        TESTS::

            sage: tan(a^2 + T.O(5))
            a^2 + O(a, b)^5
        """
        if not self[0].is_zero():
            raise ValueError('can only apply tan to formal power '
                             'series with zero constant term')
        assert(self.valuation() >= 1)
        return self.sin(prec) / self.cos(prec)

    def sinh(self, prec=infinity):
        r"""
        Apply sinh to the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT: a new power series

        EXAMPLES:

        For one variable::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: f = (t + t**2).O(4)
            sage: sinh(f)
            t + t^2 + 1/6*t^3 + O(t^4)

        For several variables::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(3)
            sage: sinh(f)
            a + b + a*b + O(a, b)^3
            sage: f.sinh()
            a + b + a*b + O(a, b)^3
            sage: f.sinh(prec=2)
            a + b + O(a, b)^2

        If the power series has a nonzero constant coefficient `c`,
        one raises an error::

            sage: g = 2 + f
            sage: sinh(g)
            Traceback (most recent call last):
            ...
            ValueError: can only apply sinh to formal power series with zero
            constant term

        If no precision is specified, the default precision is used::

            sage: T.default_prec()
            12
            sage: sinh(a)
            a + 1/6*a^3 + 1/120*a^5 + 1/5040*a^7 + 1/362880*a^9 +
             1/39916800*a^11 + O(a, b)^12
            sage: a.sinh(prec=5)
            a + 1/6*a^3 + O(a, b)^5
            sage: sinh(a + T.O(5))
            a + 1/6*a^3 + O(a, b)^5

        TESTS::

            sage: sinh(a^2 + T.O(5))
            a^2 + O(a, b)^5
        """
        R = self.parent()

        c = self[0]
        if not c.is_zero():
            raise ValueError('can only apply sinh to formal power '
                             'series with zero constant term')
        val = self.valuation()
        assert(val >= 1)

        x = self

        prec = min(prec, self.prec())
        if isinstance(prec, InfinityElement):
            prec = R.default_prec()
        n_inv_factorial = R.base_ring().one()
        x_pow_n = x
        x2 = x ** 2
        sinh_x = x.add_bigoh(prec)
        for n in range(3, prec // val + 1, 2):
            x_pow_n = (x_pow_n * x2).add_bigoh(prec)
            n_inv_factorial /= n * (n - 1)
            sinh_x += x_pow_n * n_inv_factorial
        result_bg = sinh_x

        if result_bg.base_ring() is not self.base_ring():
            R = R.change_ring(self.base_ring().fraction_field())
        return R(result_bg, prec=prec)

    def cosh(self, prec=infinity):
        r"""
        Apply cosh to the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT: a new power series

        EXAMPLES:

        For one variable::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: f = (t + t**2).O(4)
            sage: cosh(f)
            1 + 1/2*t^2 + t^3 + O(t^4)

        For several variables::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(3)
            sage: cosh(f)
            1 + 1/2*a^2 + a*b + 1/2*b^2 + O(a, b)^3
            sage: f.cosh()
            1 + 1/2*a^2 + a*b + 1/2*b^2 + O(a, b)^3
            sage: f.cosh(prec=2)
            1 + O(a, b)^2

        If the power series has a nonzero constant coefficient `c`,
        one raises an error::

            sage: g = 2 + f
            sage: cosh(g)
            Traceback (most recent call last):
            ...
            ValueError: can only apply cosh to formal power series with zero
            constant term

        If no precision is specified, the default precision is used::

            sage: T.default_prec()
            12
            sage: cosh(a)
            1 + 1/2*a^2 + 1/24*a^4 + 1/720*a^6 + 1/40320*a^8 + 1/3628800*a^10 +
             O(a, b)^12
            sage: a.cosh(prec=5)
            1 + 1/2*a^2 + 1/24*a^4 + O(a, b)^5
            sage: cosh(a + T.O(5))
            1 + 1/2*a^2 + 1/24*a^4 + O(a, b)^5

        TESTS::

            sage: cosh(a^2 + T.O(5))
            1 + 1/2*a^4 + O(a, b)^5
        """
        R = self.parent()

        c = self[0]
        if not c.is_zero():
            raise ValueError('can only apply cosh to formal power '
                             'series with zero constant term')
        x = self
        val = x.valuation()
        assert(val >= 1)

        prec = min(prec, self.prec())
        if isinstance(prec, InfinityElement):
            prec = R.default_prec()
        n_inv_factorial = R.base_ring().one()
        x_pow_n = R.one()
        x2 = x ** 2
        cosh_x = R.one().add_bigoh(prec)
        for n in range(2, prec // val + 1, 2):
            x_pow_n = (x_pow_n * x2).add_bigoh(prec)
            n_inv_factorial /= n * (n - 1)
            cosh_x += x_pow_n * n_inv_factorial
        result_bg = cosh_x

        if result_bg.base_ring() is not self.base_ring():
            R = R.change_ring(self.base_ring().fraction_field())
        return R(result_bg, prec=prec)

    def tanh(self, prec=infinity):
        r"""
        Apply tanh to the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT: a new power series

        EXAMPLES:

        For one variable::

            sage: t = PowerSeriesRing(QQ, 't').gen()
            sage: f = (t + t**2).O(4)
            sage: tanh(f)
            t + t^2 - 1/3*t^3 + O(t^4)

        For several variables::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(3)
            sage: tanh(f)
            a + b + a*b + O(a, b)^3
            sage: f.tanh()
            a + b + a*b + O(a, b)^3
            sage: f.tanh(prec=2)
            a + b + O(a, b)^2

        If the power series has a nonzero constant coefficient `c`,
        one raises an error::

            sage: g = 2 + f
            sage: tanh(g)
            Traceback (most recent call last):
            ...
            ValueError: can only apply tanh to formal power series with zero
             constant term

        If no precision is specified, the default precision is used::

            sage: T.default_prec()
            12
            sage: tanh(a)
            a - 1/3*a^3 + 2/15*a^5 - 17/315*a^7 + 62/2835*a^9 -
            1382/155925*a^11 + O(a, b)^12
            sage: a.tanh(prec=5)
            a - 1/3*a^3 + O(a, b)^5
            sage: tanh(a + T.O(5))
            a - 1/3*a^3 + O(a, b)^5

        TESTS::

            sage: tanh(a^2 + T.O(5))
            a^2 + O(a, b)^5
        """
        if not self[0].is_zero():
            raise ValueError('can only apply tanh to formal power '
                             'series with zero constant term')
        assert(self.valuation() >= 1)
        return self.sinh(prec) / self.cosh(prec)

    def O(self, prec):
        r"""
        Return this series plus `O(x^\text{prec})`. Does not change
        ``self``.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: p = 1 + x^2 + x^10; p
            1 + x^2 + x^10
            sage: p.O(15)
            1 + x^2 + x^10 + O(x^15)
            sage: p.O(5)
            1 + x^2 + O(x^5)
            sage: p.O(-5)
            Traceback (most recent call last):
            ...
            ValueError: prec (= -5) must be nonnegative
        """
        if prec is infinity or prec >= self.prec():
            return self
        coeffs = self[:prec]
        return self._parent(coeffs, prec)

    def solve_linear_de(self, prec=infinity, b=None, f0=None):
        r"""
        Obtain a power series solution to an inhomogeneous linear
        differential equation of the form:

        .. MATH::

              f'(t) = a(t) f(t) + b(t).

        INPUT:

        - ``self`` -- the power series `a(t)`

        - ``b`` -- the power series `b(t)` (default: zero)

        - ``f0`` -- the constant term of `f` ("initial condition")
          (default: 1)

        - ``prec`` -- desired precision of result (this will be
          reduced if either a or b have less precision available)

        OUTPUT: the power series `f`, to indicated precision

        ALGORITHM: A divide-and-conquer strategy; see the source code.
        Running time is approximately `M(n) \log n`, where
        `M(n)` is the time required for a polynomial multiplication
        of length `n` over the coefficient ring. (If you're working
        over something like `\QQ`, running time analysis can be a
        little complicated because the coefficients tend to explode.)

        .. NOTE::

           - If the coefficient ring is a field of characteristic
             zero, then the solution will exist and is unique.

           - For other coefficient rings, things are more
             complicated. A solution may not exist, and if it does it
             may not be unique. Generally, by the time the `n`-th term
             has been computed, the algorithm will have attempted
             divisions by `n!` in the coefficient ring. So if
             your coefficient ring has enough 'precision', and if your
             coefficient ring can perform divisions even when the
             answer is not unique, and if you know in advance that a
             solution exists, then this function will find a solution
             (otherwise it will probably crash).

        AUTHORS:

        - David Harvey (2006-09-11): factored functionality out from
          exp() function, cleaned up precision tests a bit

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)

        ::

            sage: a = 2 - 3*t + 4*t^2 + O(t^10)
            sage: b = 3 - 4*t^2 + O(t^7)
            sage: f = a.solve_linear_de(prec=5, b=b, f0=3/5)
            sage: f
             3/5 + 21/5*t + 33/10*t^2 - 38/15*t^3 + 11/24*t^4 + O(t^5)
            sage: f.derivative() - a*f - b
             O(t^4)

        ::

            sage: a = 2 - 3*t + 4*t^2
            sage: b = b = 3 - 4*t^2
            sage: f = a.solve_linear_de(b=b, f0=3/5)
            Traceback (most recent call last):
            ...
            ValueError: cannot solve differential equation to infinite precision

        ::

            sage: a.solve_linear_de(prec=5, b=b, f0=3/5)
             3/5 + 21/5*t + 33/10*t^2 - 38/15*t^3 + 11/24*t^4 + O(t^5)
        """
        if b is None:
            b = self._parent(0)

        if f0 is None:
            f0 = 1
        f0 = self.base_ring()(f0)

        # reduce precision to whatever is available from a and b
        prec = min(prec, self.prec() + 1, b.prec() + 1)
        if prec is infinity:
            raise ValueError("cannot solve differential equation to infinite precision")
        prec = int(prec)
        if prec == 0:
            return self._parent(0, 0)
        if prec < 0:
            raise ValueError("prec (=%s) must be a nonnegative integer" % prec)

        base_ring = self._parent.base_ring()
        R = PolynomialRing(base_ring, self._parent.variable_name())

        a_list = self.list()
        b_list = [base_ring(0)]
        b_list.extend(b.list())
        while len(b_list) < prec:
            b_list.append(base_ring(0))
        f = _solve_linear_de(R, 0, prec, a_list, b_list, f0)
        return self._parent(f, prec)

    def exp(self, prec=None):
        r"""
        Return exp of this power series to the indicated precision.

        INPUT:

        - ``prec`` -- integer (default: ``self.parent().default_prec``)

        ALGORITHM: See :meth:`.solve_linear_de`.

        .. NOTE::

           - Screwy things can happen if the coefficient ring is not a
             field of characteristic zero. See :meth:`.solve_linear_de`.

        AUTHORS:

        - David Harvey (2006-09-08): rewrote to use simplest possible "lazy" algorithm.

        - David Harvey (2006-09-10): rewrote to use divide-and-conquer
          strategy.

        - David Harvey (2006-09-11): factored functionality out to
          solve_linear_de().

        - Sourav Sen Gupta, David Harvey (2008-11): handle constant
          term

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)

        Check that `\exp(t)` is, well, `\exp(t)`::

            sage: (t + O(t^10)).exp()
            1 + t + 1/2*t^2 + 1/6*t^3 + 1/24*t^4 + 1/120*t^5 + 1/720*t^6
              + 1/5040*t^7 + 1/40320*t^8 + 1/362880*t^9 + O(t^10)

        Check that `\exp(\log(1+t))` is `1+t`::

            sage: (sum([-(-t)^n/n for n in range(1, 10)]) + O(t^10)).exp()
            1 + t + O(t^10)

        Check that `\exp(2t + t^2 - t^5)` is whatever it is::

            sage: (2*t + t^2 - t^5 + O(t^10)).exp()
            1 + 2*t + 3*t^2 + 10/3*t^3 + 19/6*t^4 + 8/5*t^5 - 7/90*t^6
              - 538/315*t^7 - 425/168*t^8 - 30629/11340*t^9 + O(t^10)

        Check requesting lower precision::

            sage: (t + t^2 - t^5 + O(t^10)).exp(5)
            1 + t + 3/2*t^2 + 7/6*t^3 + 25/24*t^4 + O(t^5)

        Can't get more precision than the input::

            sage: (t + t^2 + O(t^3)).exp(10)
            1 + t + 3/2*t^2 + O(t^3)

        Check some boundary cases::

            sage: (t + O(t^2)).exp(1)
            1 + O(t)
            sage: (t + O(t^2)).exp(0)
            O(t^0)

        Handle nonzero constant term (fixes :issue:`4477`)::

            sage: # needs sage.rings.real_mpfr
            sage: R.<x> = PowerSeriesRing(RR)
            sage: (1 + x + x^2 + O(x^3)).exp()
            2.71828182845905 + 2.71828182845905*x + 4.07742274268857*x^2 + O(x^3)

        ::

            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: (1 + x + O(x^2)).exp()                                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ArithmeticError: exponential of constant term does not belong
            to coefficient ring (consider working in a larger ring)

        ::

            sage: R.<x> = PowerSeriesRing(GF(5))
            sage: (1 + x + O(x^2)).exp()
            Traceback (most recent call last):
            ...
            ArithmeticError: constant term of power series does not support exponentiation
        """
        if prec is None:
            prec = self._parent.default_prec()

        t = self.derivative().solve_linear_de(prec)

        if not self[0].is_zero():
            try:
                C = self[0].exp()
            except AttributeError:
                raise ArithmeticError("constant term of power series does not support exponentiation")

            if C.parent() is not self.base_ring():
                raise ArithmeticError("exponential of constant term does not belong to coefficient ring (consider working in a larger ring)")

            t = C * t

        return t

    def log(self, prec=None):
        r"""
        Return log of this power series to the indicated precision.

        This works only if the constant term of the power series is 1
        or the base ring can take the logarithm of the constant coefficient.

        INPUT:

        - ``prec`` -- integer (default: ``self.parent().default_prec()``)

        ALGORITHM: See :meth:`solve_linear_de()`.

        .. WARNING::

            Screwy things can happen if the coefficient ring is not a
            field of characteristic zero. See :meth:`solve_linear_de()`.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)
            sage: (1 + t + O(t^10)).log()
            t - 1/2*t^2 + 1/3*t^3 - 1/4*t^4 + 1/5*t^5 - 1/6*t^6 + 1/7*t^7
              - 1/8*t^8 + 1/9*t^9 + O(t^10)

            sage: t.exp().log()
            t + O(t^10)

            sage: (1 + t).log().exp()
            1 + t + O(t^10)

            sage: (-1 + t + O(t^10)).log()
            Traceback (most recent call last):
            ...
            ArithmeticError: constant term of power series is not 1

            sage: # needs sage.rings.real_mpfr
            sage: R.<t> = PowerSeriesRing(RR)
            sage: (2 + t).log().exp()
            2.00000000000000 + 1.00000000000000*t + O(t^20)
        """
        if prec is None:
            prec = self._parent.default_prec()

        zero = self.parent().zero()
        const_off = zero

        if not self[0].is_one():
            if self.base_ring() in _Fields and self[0] > 0:
                const_off = self[0].log()
            else:
                raise ArithmeticError("constant term of power series is not 1")

        t = zero.solve_linear_de(prec,b=self.derivative()/self,f0=0)
        return t + const_off

    def V(self, n):
        r"""
        If `f = \sum a_m x^m`, then this function returns
        `\sum a_m x^{nm}`.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: p = 1 + x^2 + x^10; p
            1 + x^2 + x^10
            sage: p.V(3)
            1 + x^6 + x^30
            sage: (p + O(x^20)).V(3)
            1 + x^6 + x^30 + O(x^60)
        """
        v = self.list()
        m = 0
        w = []
        zero = self.base_ring()(0)
        for i in range(len(v)*n):
            if i%n != 0:
                w.append(zero)
            else:
                w.append(v[m])
                m += 1
        return self._parent(w, self.prec()*n)

    def valuation(self):
        """
        Return the valuation of this power series.

        This is equal to the valuation of the underlying polynomial.

        EXAMPLES:

        Sparse examples::

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: f = t^100000 + O(t^10000000)
            sage: f.valuation()
            100000
            sage: R(0).valuation()
            +Infinity

        Dense examples::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: f = 17*t^100 + O(t^110)
            sage: f.valuation()
            100
            sage: t.valuation()
            1
        """
        return self.polynomial().valuation()

    def variable(self):
        """
        Return a string with the name of the variable
        of this power series.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(Rationals())
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f.variable()
            'x'

        AUTHORS:

        - David Harvey (2006-08-08)
        """
        return self._parent.variable_name()

    def degree(self):
        """
        Return the degree of this power series, which is by definition the
        degree of the underlying polynomial.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: f = t^100000 + O(t^10000000)
            sage: f.degree()
            100000
        """
        return self.polynomial().degree()

    def derivative(self, *args):
        r"""
        The formal derivative of this power series, with respect to
        variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        .. SEEALSO::

           :meth:`_derivative`

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(QQ)
            sage: g = -x + x^2/2 - x^4 + O(x^6)
            sage: g.derivative()
            -1 + x - 4*x^3 + O(x^5)
            sage: g.derivative(x)
            -1 + x - 4*x^3 + O(x^5)
            sage: g.derivative(x, x)
            1 - 12*x^2 + O(x^4)
            sage: g.derivative(x, 2)
            1 - 12*x^2 + O(x^4)
        """
        return multi_derivative(self, args)

    def __setitem__(self, n, value):
        """
        Called when an attempt is made to change a power series.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: f = 3 - t^3 + O(t^5)
            sage: f[1] = 5
            Traceback (most recent call last):
            ...
            IndexError: power series are immutable
        """
        raise IndexError("power series are immutable")

    def laurent_series(self):
        """
        Return the Laurent series associated to this power series, i.e.,
        this series considered as a Laurent series.

        EXAMPLES::

            sage: k.<w> = QQ[[]]
            sage: f = 1 + 17*w + 15*w^3 + O(w^5)
            sage: parent(f)
            Power Series Ring in w over Rational Field
            sage: g = f.laurent_series(); g
            1 + 17*w + 15*w^3 + O(w^5)
        """
        return self._parent.laurent_series_ring()(self)

    def egf_to_ogf(self):
        r"""
        Return the ordinary generating function power series,
        assuming ``self`` is an exponential generating function power series.

        This is a formal Laplace transform.

        This function is known as :pari:`serlaplace` in PARI/GP.

        .. SEEALSO:: :meth:`ogf_to_egf` for the inverse method.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t + t^2/factorial(2) + 2*t^3/factorial(3)
            sage: f.egf_to_ogf()
            t + t^2 + 2*t^3
        """
        return self.parent()([self[i] * arith.factorial(i) for i in range(self.degree() + 1)])

    def ogf_to_egf(self):
        r"""
        Return the exponential generating function power series,
        assuming ``self`` is an ordinary generating function power series.

        This is a formal Borel transform.

        This can also be computed as ``serconvol(f,exp(t))`` in PARI/GP.

        .. SEEALSO:: :meth:`egf_to_ogf` for the inverse method.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t + t^2 + 2*t^3
            sage: f.ogf_to_egf()
            t + 1/2*t^2 + 1/3*t^3
        """
        return self.parent()([self[i] / arith.factorial(i) for i in range(self.degree() + 1)])

    def __pari__(self):
        """
        Return a PARI representation of this series.

        There are currently limits to the possible base rings over which this
        function works.  See the documentation for
        :meth:`~sage.rings.polynomial.polynomial_element.Polynomial.__pari__`

        EXAMPLES::

            sage: k.<w> = QQ[[]]
            sage: f = 1 + 17*w + 15*w^3 + O(w^5)
            sage: pari(f)  # indirect doctest                                           # needs sage.libs.pari
            1 + 17*w + 15*w^3 + O(w^5)
            sage: pari(1 - 19*w + w^5)  # indirect doctest                              # needs sage.libs.pari
            w^5 - 19*w + 1
            sage: R.<x> = Zmod(6)[[]]
            sage: pari(1 + x + 8*x^3 + O(x^8))  # indirect doctest                      # needs sage.libs.pari
            Mod(1, 6) + Mod(1, 6)*x + Mod(2, 6)*x^3 + O(x^8)

        TESTS::

            sage: pari(1 + O(x^1))                                                      # needs sage.libs.pari
            Mod(1, 6) + O(x)
            sage: pari(O(x^1))                                                          # needs sage.libs.pari
            O(x)
            sage: pari(O(x^0))                                                          # needs sage.libs.pari
            O(x^0)
        """
        n = self.prec()
        s = self.truncate().__pari__()  # PARI polynomial
        if n is not infinity:
            v = s.variable()
            s = s.Ser(v, n - s.valuation(v) if s else n)
        return s


def _solve_linear_de(R, N, L, a, b, f0):
    r"""
    Internal function used by ``PowerSeries.solve_linear_de()``.

    INPUT:

    - ``R`` -- a PolynomialRing

    - ``N`` -- integer >= 0

    - ``L`` -- integer >= 1

    - ``a`` -- list of coefficients of `a`, any length, all coefficients should
      belong to base ring of `R`

    - ``b`` -- list of coefficients of `b`, length at least `L` (only first `L`
      coefficients are used), all coefficients should belong to base ring of `R`

    - ``f0`` -- constant term of `f` (only used if ``N == 0``), should belong
      to base ring of `R`

    OUTPUT: list of coefficients of `f` (length exactly `L`), where `f` is a
    solution to the linear inhomogeneous differential equation:

    .. MATH::

         (t^N f)'  =  a t^N f  +  t^{N-1} b  +  O(t^{N+L-1}).


    The constant term of `f` is taken to be f0 if `N == 0`, and otherwise is
    determined by `N f_0 = b_0`.

    ALGORITHM: The algorithm implemented here is inspired by the "fast
    lazy multiplication algorithm" described in the paper "Lazy
    multiplication of formal power series" by J van der Hoeven (1997
    ISSAC conference).

    Our situation is a bit simpler than the one described there,
    because only one of the series being multiplied needs the lazy
    treatment; the other one is known fully in advance. I have
    reformulated the algorithm as an explicit divide-and-conquer
    recursion. Possibly it is slower than the one described by van der
    Hoeven, by a constant factor, but it seems easier to implement.
    Also, because we repeatedly split in half starting at the top
    level, instead of repeatedly doubling in size from the bottom
    level, it's easier to avoid problems with "overshooting" in the
    last iteration.

    The idea is to split the problem into two instances with
    `L` about half the size. Take `L'` to be the
    ceiling of `(L/2)`. First recursively find `g`
    modulo `t^{L'}` such that

    .. MATH::

         (t^N g)'  =  a t^N g  +  t^{N-1} b  +  O(t^{N+L'-1}).

    Next we want to find `h` modulo `t^{L-L'}` such
    that `f = g + t^{L'} h` is a solution of the original
    equation. We can find a suitable `h` by recursively solving
    another differential equation of the same form, namely

    .. MATH::

         (t^{N+L'} h)'  =  a t^{N+L'} h  +  c t^{N+L'-1} + O(t^{N+L'-1}),


    where `c` is determined by

    .. MATH::

         (t^N g)' - a t^N g - t^{N-1} b  =  -c t^{N+L'-1} + O(t^{N+L-1}).


    Once `g` is known, `c` can be recovered from this
    relation by computing the coefficients of `t^j` of the
    product `a g` for `j` in the range
    `L'-1 \leq j < L-1`.

    At the bottom of the recursion we have `L = 1`, in which
    case the solution is simply given by `f_0 = b_0/N` (or by
    the supplied initial condition `f_0` when
    `N == 0`).

    The algorithm has to do one multiplication of length `N`,
    two of length `N/2`, four of length `N/4`, etc,
    hence giving runtime `O(M(N) \log N)`.

    AUTHORS:

    - David Harvey (2006-09-11): factored out of PowerSeries.exp().
    """
    if L == 1:
        # base case
        if N == 0:
            return [f0]
        else:
            return [b[0] / N]

    L2 = (L + 1) >> 1    # ceil(L/2)

    g = _solve_linear_de(R, N, L2, a, b, f0)

    term1 = R(g)  # we must not have check=False, since otherwise [..., 0, 0] is not stripped
    term2 = R(a[:L])  # , check=False)
    product = (term1 * term2).list()

    # todo: perhaps next loop could be made more efficient
    c = b[L2:L]
    for j in range(L2 - 1, min(L-1, len(product))):
        c[j - (L2-1)] = c[j - (L2-1)] + product[j]

    h = _solve_linear_de(R, N + L2, L - L2, a, c, f0)

    g.extend(h)
    return g


def make_powerseries_poly_v0(parent,  f, prec, is_gen):
    # This is only used to unpickle old pickles. The new pickling
    # works differently!
    from sage.rings import power_series_poly
    return power_series_poly.PowerSeries_poly(parent, f, prec, 0, is_gen)


def make_element_from_parent_v0(parent, *args):
    # This is only used to unpickle old pickles. The new pickling
    # works differently!
    return parent(*args)
