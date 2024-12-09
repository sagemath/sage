r"""
Multivariate Power Series

Construct and manipulate multivariate power series (in finitely many
variables) over a given commutative ring. Multivariate power series
are implemented with total-degree precision.

EXAMPLES:

Power series arithmetic, tracking precision::

    sage: R.<s,t> = PowerSeriesRing(ZZ); R
    Multivariate Power Series Ring in s, t over Integer Ring

    sage: f = 1 + s + 3*s^2; f
    1 + s + 3*s^2
    sage: g = t^2*s + 3*t^2*s^2 + R.O(5); g
    s*t^2 + 3*s^2*t^2 + O(s, t)^5
    sage: g = t^2*s + 3*t^2*s^2 + O(s, t)^5; g
    s*t^2 + 3*s^2*t^2 + O(s, t)^5
    sage: f = f.O(7); f
    1 + s + 3*s^2 + O(s, t)^7
    sage: f += s; f
    1 + 2*s + 3*s^2 + O(s, t)^7
    sage: f*g
    s*t^2 + 5*s^2*t^2 + O(s, t)^5
    sage: (f-1)*g
    2*s^2*t^2 + 9*s^3*t^2 + O(s, t)^6
    sage: f*g - g
    2*s^2*t^2 + O(s, t)^5

    sage: f *= s; f
    s + 2*s^2 + 3*s^3 + O(s, t)^8
    sage: f%2
    s + s^3 + O(s, t)^8
    sage: (f%2).parent()
    Multivariate Power Series Ring in s, t over Ring of integers modulo 2

As with univariate power series, comparison of `f` and `g` is
done up to the minimum precision of `f` and `g`::

    sage: f = 1 + t + s + s*t + R.O(3); f
    1 + s + t + s*t + O(s, t)^3
    sage: g = s^2 + 2*s^4 - s^5 + s^2*t^3 + R.O(6); g
    s^2 + 2*s^4 - s^5 + s^2*t^3 + O(s, t)^6
    sage: f == g
    False
    sage: g == g.add_bigoh(3)
    True
    sage: f < g
    False
    sage: f > g
    True

Calling::

    sage: f = s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + R.O(5); f
    s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + O(s, t)^5
    sage: f(t, s)
    s*t + t^2 + s*t^2 + t^3 + 3*s*t^3 + 3*t^4 + O(s, t)^5
    sage: f(t^2, s^2)
    s^2*t^2 + t^4 + s^2*t^4 + t^6 + 3*s^2*t^6 + 3*t^8 + O(s, t)^10

Substitution is defined only for elements of positive valuation, unless `f`
has infinite precision::

    sage: f(t^2, s^2 + 1)
    Traceback (most recent call last):
    ...
    TypeError: Substitution defined only for elements of positive valuation,
    unless self has infinite precision.

    sage: g = f.truncate()
    sage: g(t^2, s^2 + 1)
    t^2 + s^2*t^2 + 2*t^4 + s^2*t^4 + 4*t^6 + 3*s^2*t^6 + 3*t^8
    sage: g(t^2, (s^2+1).O(3))
    t^2 + s^2*t^2 + 2*t^4 + O(s, t)^5

0 has valuation ``+Infinity``::

    sage: f(t^2, 0)
    t^4 + t^6 + 3*t^8 + O(s, t)^10
    sage: f(t^2, s^2 + s)
    s*t^2 + s^2*t^2 + t^4 + O(s, t)^5

Substitution of power series with finite precision works too::

    sage: f(s.O(2), t)
    s^2 + s*t + O(s, t)^3
    sage: f(f, f)
    2*s^4 + 4*s^3*t + 2*s^2*t^2 + 4*s^5 + 8*s^4*t + 4*s^3*t^2 + 16*s^6 +
    34*s^5*t + 20*s^4*t^2 + 2*s^3*t^3 + O(s, t)^7
    sage: t(f, f)
    s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + O(s, t)^5
    sage: t(0, f) == s(f, 0)
    True

The ``subs`` syntax works as expected::

    sage: r0 = -t^2 - s*t^3 - 2*t^6 + s^7 + s^5*t^2 + R.O(10)
    sage: r1 = s^4 - s*t^4 + s^6*t - 4*s^2*t^5 - 6*s^3*t^5 + R.O(10)
    sage: r2 = 2*s^3*t^2 - 2*s*t^4 - 2*s^3*t^4 + s*t^7 + R.O(10)
    sage: r0.subs({t: r2, s: r1})
    -4*s^6*t^4 + 8*s^4*t^6 - 4*s^2*t^8 + 8*s^6*t^6 - 8*s^4*t^8 - 4*s^4*t^9
    + 4*s^2*t^11 - 4*s^6*t^8 + O(s, t)^15
    sage: r0.subs({t: r2, s: r1}) == r0(r1, r2)
    True

Construct ring homomorphisms from one power series ring to another::

    sage: A.<a,b> = PowerSeriesRing(QQ)
    sage: X.<x,y> = PowerSeriesRing(QQ)

    sage: phi = Hom(A,X)([x,2*y]); phi
    Ring morphism:
      From: Multivariate Power Series Ring in a, b over Rational Field
      To:   Multivariate Power Series Ring in x, y over Rational Field
      Defn: a |--> x
            b |--> 2*y

    sage: phi(a+b+3*a*b^2 + A.O(5))
    x + 2*y + 12*x*y^2 + O(x, y)^5

Multiplicative inversion of power series::

    sage: h = 1 + s + t + s*t + s^2*t^2 + 3*s^4 + 3*s^3*t + R.O(5)
    sage: k = h^-1; k
    1 - s - t + s^2 + s*t + t^2 - s^3 - s^2*t - s*t^2 - t^3 - 2*s^4 -
    2*s^3*t + s*t^3 + t^4 + O(s, t)^5
    sage: h*k
    1 + O(s, t)^5

    sage: f = 1 - 5*s^29 - 5*s^28*t + 4*s^18*t^35 + \
    ....: 4*s^17*t^36 - s^45*t^25 - s^44*t^26 + s^7*t^83 + \
    ....: s^6*t^84 + R.O(101)
    sage: h = ~f; h
    1 + 5*s^29 + 5*s^28*t - 4*s^18*t^35 - 4*s^17*t^36 + 25*s^58 + 50*s^57*t
    + 25*s^56*t^2 + s^45*t^25 + s^44*t^26 - 40*s^47*t^35 - 80*s^46*t^36
    - 40*s^45*t^37 + 125*s^87 + 375*s^86*t + 375*s^85*t^2 + 125*s^84*t^3
    - s^7*t^83 - s^6*t^84 + 10*s^74*t^25 + 20*s^73*t^26 + 10*s^72*t^27
    + O(s, t)^101
    sage: h*f
    1 + O(s, t)^101

AUTHORS:

- Niles Johnson (07/2010): initial code
- Simon King (08/2012): Use category and coercion framework, :issue:`13412`
"""
# ****************************************************************************
#       Copyright (C) 2010 Niles Johnson <nilesj@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.infinity import infinity, InfinityElement
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.power_series_ring_element import PowerSeries
from sage.structure.richcmp import richcmp

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')
lazy_import('sage.rings.multi_power_series_ring', 'MPowerSeriesRing_generic')
lazy_import('sage.rings.power_series_ring', 'PowerSeriesRing_generic')


def is_MPowerSeries(f):
    """
    Return ``True`` if ``f`` is a multivariate power series.

    TESTS::

        sage: from sage.rings.power_series_ring_element import is_PowerSeries
        sage: from sage.rings.multi_power_series_ring_element import is_MPowerSeries
        sage: M = PowerSeriesRing(ZZ,4,'v')
        sage: is_PowerSeries(M.random_element(10))
        doctest:warning...
        DeprecationWarning: The function is_PowerSeries is deprecated; use 'isinstance(..., PowerSeries)' instead.
        See https://github.com/sagemath/sage/issues/38266 for details.
        True
        sage: is_MPowerSeries(M.random_element(10))
        doctest:warning...
        DeprecationWarning: The function is_MPowerSeries is deprecated; use 'isinstance(..., MPowerSeries)' instead.
        See https://github.com/sagemath/sage/issues/38266 for details.
        True
        sage: T.<v> = PowerSeriesRing(RR)
        sage: is_MPowerSeries(1 - v + v^2 +O(v^3))
        False
        sage: is_PowerSeries(1 - v + v^2 +O(v^3))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38266,
                "The function is_MPowerSeries is deprecated; "
                "use 'isinstance(..., MPowerSeries)' instead.")
    return isinstance(f, MPowerSeries)


class MPowerSeries(PowerSeries):
    ### methods from PowerSeries that we *don't* override:
    #
    # __hash__ : works just fine
    #
    # __reduce__ : don't really understand this
    #
    # is_sparse : works just fine
    #
    # is_dense : works just fine
    #
    # is_gen : works just fine
    #
    # base_extend : works just fine
    #
    # change_ring : works just fine
    #
    # _cmp_ : don't understand this
    #
    # __copy__ : works just fine
    #
    # base_ring : works just fine
    #
    # common_prec : works just fine
    #
    # common_prec_c : seems fine
    #
    # _mul_prec : works just fine
    #
    # __bool__ : works just fine
    #
    """
    Multivariate power series; these are the elements of Multivariate Power
    Series Rings.

    INPUT:

    - ``parent`` -- a multivariate power series

    - ``x`` -- the element (default: 0).  This can be another
      :class:`MPowerSeries` object, or an element of one of the following:

      - the background univariate power series ring
      - the foreground polynomial ring
      - a ring that coerces to one of the above two

    - ``prec`` -- (default: ``infinity``) the precision

    - ``is_gen`` -- boolean (default: ``False``); whether this element is one
      of the generators

    - ``check`` -- boolean (default: ``False``); needed by univariate power
      series class

    EXAMPLES:

    Construct multivariate power series from generators::

        sage: S.<s,t> = PowerSeriesRing(ZZ)
        sage: f = s + 4*t + 3*s*t
        sage: f in S
        True
        sage: f = f.add_bigoh(4); f
        s + 4*t + 3*s*t + O(s, t)^4
        sage: g = 1 + s + t - s*t + S.O(5); g
        1 + s + t - s*t + O(s, t)^5

        sage: T = PowerSeriesRing(GF(3),5,'t'); T
        Multivariate Power Series Ring in t0, t1, t2, t3, t4
         over Finite Field of size 3
        sage: t = T.gens()
        sage: w = t[0] - 2*t[1]*t[3] + 5*t[4]^3 - t[0]^3*t[2]^2; w
        t0 + t1*t3 - t4^3 - t0^3*t2^2
        sage: w = w.add_bigoh(5); w
        t0 + t1*t3 - t4^3 + O(t0, t1, t2, t3, t4)^5
        sage: w in T
        True

        sage: w = t[0] - 2*t[0]*t[2] + 5*t[4]^3 - t[0]^3*t[2]^2 + T.O(6)
        sage: w
        t0 + t0*t2 - t4^3 - t0^3*t2^2 + O(t0, t1, t2, t3, t4)^6


    Get random elements::

        sage: S.random_element(4)   # random
        -2*t + t^2 - 12*s^3 + O(s, t)^4

        sage: T.random_element(10)  # random
        -t1^2*t3^2*t4^2 + t1^5*t3^3*t4 + O(t0, t1, t2, t3, t4)^10


    Convert elements from polynomial rings::

        sage: # needs sage.rings.finite_rings
        sage: R = PolynomialRing(ZZ, 5, T.variable_names())
        sage: t = R.gens()
        sage: r = -t[2]*t[3] + t[3]^2 + t[4]^2
        sage: T(r)
        -t2*t3 + t3^2 + t4^2
        sage: r.parent()
        Multivariate Polynomial Ring in t0, t1, t2, t3, t4 over Integer Ring
        sage: r in T
        True
    """

    def __init__(self, parent, x=0, prec=infinity, is_gen=False, check=False):
        """
        Input ``x`` can be an :class:`MPowerSeries`, or an element of

            - the background univariate power series ring
            - the foreground polynomial ring
            - a ring that coerces to one of the above two

        EXAMPLES::

            sage: R.<s,t> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in s, t over Integer Ring
            sage: f = 1 + t + s + s*t + R.O(3)
            sage: g = (1/2) * f; g
            1/2 + 1/2*s + 1/2*t + 1/2*s*t + O(s, t)^3
            sage: g.parent()
            Multivariate Power Series Ring in s, t over Rational Field
            sage: g = (1/2)*f; g
            1/2 + 1/2*s + 1/2*t + 1/2*s*t + O(s, t)^3
            sage: g.parent()
            Multivariate Power Series Ring in s, t over Rational Field

            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x - 3,'a')                                            # needs sage.rings.number_field
            sage: g = K.random_element()*f                                              # needs sage.rings.number_field
            sage: g.parent()                                                            # needs sage.rings.number_field
            Multivariate Power Series Ring in s, t over
             Number Field in a with defining polynomial x - 3

        TESTS::

            sage: S.<s,t> = PowerSeriesRing(ZZ)
            sage: f = s + 4*t + 3*s*t
            sage: f in S
            True
            sage: f = f.add_bigoh(4); f
            s + 4*t + 3*s*t + O(s, t)^4
            sage: g = 1 + s + t - s*t + S.O(5); g
            1 + s + t - s*t + O(s, t)^5

            sage: B.<s, t> = PowerSeriesRing(QQ)
            sage: C.<z> = PowerSeriesRing(QQ)
            sage: B(z)
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce input to polynomial ring.

            sage: D.<s> = PowerSeriesRing(QQ)
            sage: s.parent() is D
            True
            sage: B(s) in B
            True
            sage: d = D.random_element(20)
            sage: b = B(d) # test coercion from univariate power series ring
            sage: b in B
            True
        """
        PowerSeries.__init__(self, parent, prec, is_gen=is_gen)
        self._PowerSeries__is_gen = is_gen

        try:
            prec = min(prec, x.prec()) # use precision of input, if defined
        except AttributeError:
            pass

        # set the correct background value, depending on what type of input x is
        try:
            xparent = x.parent() # 'int' types have no parent
        except AttributeError:
            xparent = None

        # test whether x coerces to background univariate
        # power series ring of parent
        if isinstance(xparent, (PowerSeriesRing_generic, MPowerSeriesRing_generic, LazyPowerSeriesRing)):
            # x is either a multivariate or univariate power series
            #
            # test whether x coerces directly to designated parent
            if isinstance(x, MPowerSeries):
                try:
                    self._bg_value = parent._bg_ps_ring(x._bg_value)
                except TypeError:
                    raise TypeError("Unable to coerce into background ring.")

            # test whether x coerces to background univariate
            # power series ring of parent
            elif xparent == parent._bg_ps_ring():
                self._bg_value = x
            elif parent._bg_ps_ring().has_coerce_map_from(xparent):
                # previous test may fail if precision or term orderings of
                # base rings do not match
                self._bg_value = parent._bg_ps_ring(x)
            else:
                # x is a univariate power series, but not from the
                # background power series ring
                #
                # convert x to a polynomial and send to background
                # ring of parent
                x = x.polynomial()
                self._bg_value = parent._send_to_bg(x).add_bigoh(prec)

        # test whether x coerces to underlying polynomial ring of parent
        elif isinstance(xparent, PolynomialRing_general):
            self._bg_value = parent._send_to_bg(x).add_bigoh(prec)

        else:
            try:
                x = parent._poly_ring(x)
                #self._value = x
                self._bg_value = parent._send_to_bg(x).add_bigoh(prec)
            except (TypeError, AttributeError):
                raise TypeError("Input does not coerce to any of the "
                                "expected rings.")

        self._go_to_fg = parent._send_to_fg
        self._prec = self._bg_value.prec()

        # self._parent is used a lot by the class PowerSeries
        self._parent = self.parent()

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: K.<s,t> = PowerSeriesRing(QQ)
            sage: f = 1 + t - s + s*t - s*t^3 + K.O(12)
            sage: loads(dumps(f)) == f
            True
        """
        return self.__class__, (self._parent,self._bg_value,self._prec)

    def __call__(self, *x, **kwds):
        """
        Evaluate ``self``.

        EXAMPLES::

            sage: R.<s,t> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in s, t over Integer Ring
            sage: f = s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + R.O(5); f
            s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + O(s, t)^5
            sage: f(t,s)
            s*t + t^2 + s*t^2 + t^3 + 3*s*t^3 + 3*t^4 + O(s, t)^5

            sage: f(t,0)
            t^2 + t^3 + 3*t^4 + O(s, t)^5
            sage: f(t,2)
            Traceback (most recent call last):
            ...
            TypeError: Substitution defined only for elements of positive
            valuation, unless self has infinite precision.

            sage: f.truncate()(t,2)
            2*t + 3*t^2 + 7*t^3 + 3*t^4

        Checking that :issue:`15059` is fixed::

            sage: M.<u,v> = PowerSeriesRing(GF(5))
            sage: s = M.hom([u, u+v])
            sage: s(M.one())
            1

        Since :issue:`26105` you can specify a map on the base ring::

            sage: # needs sage.rings.number_field
            sage: Zx.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: cc = K.hom([-i])
            sage: R.<s,t> = PowerSeriesRing(K)
            sage: f = s^2 + i*s*t + (3+4*i)*s^3 + R.O(4); f
            s^2 + i*s*t + (4*i + 3)*s^3 + O(s, t)^4
            sage: f(t, s, base_map=cc)
            (-i)*s*t + t^2 + (-4*i + 3)*t^3 + O(s, t)^4
        """
        if len(x) != self.parent().ngens():
            raise ValueError("Number of arguments does not match number of variables in parent.")
        if kwds:
            return self._subs_formal(*x, **kwds)

        sub_dict = {}
        valn_list = []
        for i in range(len(x)):
            try:
                xi = self.parent(x[i])
            except (AttributeError, TypeError):
                # Input does not coerce to parent ring of self
                # attempt formal substitution
                return self._subs_formal(*x, **kwds)
            if xi.valuation() == 0 and self.prec() is not infinity:
                raise TypeError("Substitution defined only for elements of positive valuation, unless self has infinite precision.")
            elif xi.valuation() > 0:
                sub_dict[self.parent()._poly_ring().gens()[i]] = xi.add_bigoh(xi.valuation() * self.prec())
                valn_list.append(xi.valuation())
            else:
                sub_dict[self.parent()._poly_ring().gens()[i]] = xi
        if self.prec() is infinity:
            newprec = infinity
        else:
            newprec = self.prec()*min(valn_list)
        return self.parent()(self._value().subs(sub_dict)).add_bigoh(newprec)

    def _subs_formal(self, *x, **kwds):
        """
        Substitution of inputs as variables of ``self``.  This is formal
        in the sense that the inputs do not need to be elements of
        same multivariate power series ring as ``self``.  They can be any
        objects which support addition and multiplication with
        each other and with the coefficients of ``self``.  If ``self`` has
        finite precision, the inputs must also support an ``add_bigoh``
        method.

        TESTS::

            sage: B.<s, t> = PowerSeriesRing(QQ)
            sage: C.<z> = PowerSeriesRing(QQ)
            sage: s(z,z)
            z

            sage: f = -2/33*s*t^2 - 1/5*t^5 - s^5*t + s^2*t^4
            sage: f(z,z)  # indirect doctest
            -2/33*z^3 - 1/5*z^5
            sage: f(z,1)  # indirect doctest
            -1/5 - 2/33*z + z^2 - z^5
            sage: RF = RealField(10)                                                    # needs sage.rings.real_mpfr
            sage: f(z, RF(1))  # indirect doctest                                       # needs sage.rings.real_mpfr
            -0.20 - 0.061*z + 1.0*z^2 - 0.00*z^3 - 0.00*z^4 - 1.0*z^5

            sage: m = matrix(QQ,[[1,0,1],[0,2,1],[-1,0,0]]); m                          # needs sage.modules
            [ 1  0  1]
            [ 0  2  1]
            [-1  0  0]
            sage: f(m,m)  # indirect doctest                                            # needs sage.modules
            [     2/33         0       1/5]
            [   131/55 -1136/165    -24/11]
            [     -1/5         0   -23/165]
            sage: f(m,m) == -2/33*m^3 - 1/5*m^5  # indirect doctest                     # needs sage.modules
            True

            sage: f = f.add_bigoh(10)
            sage: f(z,z)
            -2/33*z^3 - 1/5*z^5 + O(z^10)
            sage: f(m,m)                                                                # needs sage.modules
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'
            object has no attribute 'add_bigoh'
        """
        from sage.misc.misc_c import prod

        if len(x) == 1 and isinstance(x[0], (list, tuple)):
            x = x[0]
        n = self.parent().ngens()
        if len(x) != n:
            raise ValueError("Input must be of correct length.")
        if n == 0:
            return self

        y = 0
        base_map = kwds.get('base_map')
        if base_map is None:
            base_map = lambda t: t
        for m, c in self.monomial_coefficients().items():
            y += base_map(c)*prod([x[i]**m[i] for i in range(n) if m[i] != 0])
        if self.prec() == infinity:
            return y
        else:
            return y.add_bigoh(self.prec())

    def _value(self):
        """
        Return the value of ``self`` in the foreground polynomial ring.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(GF(5)); R
            Multivariate Power Series Ring in a, b, c over Finite Field of
            size 5
            sage: f = 1 + a + b - a*b + R.O(3); f
            1 + a + b - a*b + O(a, b, c)^3
            sage: f._value()
            1 + a + b - a*b
            sage: f._value().parent()
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 5
        """
        return self._go_to_fg(self._bg_value)

    def _repr_(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: B.<s,t,v> = PowerSeriesRing(QQ)
            sage: e = 1 + s - s*t + t*v/2 - 2*s*t*v/8 + B.O(4)
            sage: e._repr_()
            '1 + s - s*t + 1/2*t*v - 1/4*s*t*v + O(s, t, v)^4'
        """
        if self._prec == infinity:
            return "%s" % self._value()
        return "%(val)s + O(%(gens)s)^%(prec)s" \
               % {'val':self._value(),
                 'gens':', '.join(str(g) for g in self.parent().gens()),
                 'prec':self._prec}

    def _latex_(self):
        """
        Return latex representation of this multivariate power series.

        EXAMPLES::

            sage: M = PowerSeriesRing(GF(5),3,'t'); M
            Multivariate Power Series Ring in t0, t1, t2 over Finite Field of size 5
            sage: t = M.gens()
            sage: f = (-t[0]^4*t[1]^3*t[2]^4 - 2*t[0]*t[1]^4*t[2]^7
            ....:      + 2*t[1]*t[2]^12 + 2*t[0]^7*t[1]^5*t[2]^2 + M.O(15))
            sage: f
            -t0^4*t1^3*t2^4 - 2*t0*t1^4*t2^7 + 2*t1*t2^12 + 2*t0^7*t1^5*t2^2
            + O(t0, t1, t2)^15
            sage: f._latex_()
            '-t_{0}^{4} t_{1}^{3} t_{2}^{4} + 3 t_{0} t_{1}^{4} t_{2}^{7} +
            2 t_{1} t_{2}^{12} + 2 t_{0}^{7} t_{1}^{5} t_{2}^{2}
            + O(t_{0}, t_{1}, t_{2})^{15}'

        TESTS:

        Check that :issue:`25156` is fixed::

            sage: R.<x1,y1> = PowerSeriesRing(QQ, ('x', 'y'))
            sage: element = 1 + y1^10 + x1^5
            sage: element._latex_()
            '1 + x_{1}^{5} + y_{1}^{10}'
        """
        if self._prec == infinity:
            return "%s" % self._value()._latex_()
        return "%(val)s + O(%(gens)s)^{%(prec)s}" \
               % {'val':self._value()._latex_(),
                 'gens':', '.join(g._latex_() for g in self.parent().gens()),
                 'prec':self._prec}

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of this series under the map that sends the
        generators to ``im_gens``. This is used internally for computing
        homomorphisms.

        EXAMPLES::

            sage: A.<a,b> = PowerSeriesRing(QQ)
            sage: X.<x,y> = PowerSeriesRing(QQ)
            sage: phi = Hom(A,X)([x,2*y])
            sage: phi = Hom(A,X)([x,2*y]); phi
            Ring morphism:
              From: Multivariate Power Series Ring in a, b over Rational Field
              To:   Multivariate Power Series Ring in x, y over Rational Field
              Defn: a |--> x
                    b |--> 2*y
            sage: phi(a+b+3*a*b^2 + A.O(5))  # indirect doctest
            x + 2*y + 12*x*y^2 + O(x, y)^5
        """
        if base_map is None:
            # __call__ might be faster if codomain coerces into the base ring
            return codomain(self(*im_gens))
        else:
            return codomain(self._subs_formal(*im_gens, base_map=base_map))

    def __getitem__(self, n):
        """
        Return summand of total degree ``n``.

        TESTS::

            sage: R.<a,b> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b + R.O(4)
            sage: f[0]
            1
            sage: f[2]
            -a*b
            sage: f[3]
            0
            sage: f[4]
            Traceback (most recent call last):
            ...
            IndexError: Cannot return terms of total degree greater than or
            equal to precision of self.
        """
        if n >= self.prec():
            raise IndexError("Cannot return terms of total degree greater than or equal to precision of self.")
        return self.parent(self._bg_value[n])

    def __invert__(self):
        """
        Return multiplicative inverse of this multivariate power series.

        Currently implemented only if constant coefficient is a unit in the
        base ring.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b - b*c - a*c + R.O(4)
            sage: ~f
            1 - a - b + a^2 + 3*a*b + a*c + b^2 + b*c - a^3 - 5*a^2*b
            - 2*a^2*c - 5*a*b^2 - 4*a*b*c - b^3 - 2*b^2*c + O(a, b, c)^4
        """
        if self.valuation() == 0:
            return self.parent(~self._bg_value)
        else:
            raise NotImplementedError("Multiplicative inverse of multivariate power series currently implemented only if constant coefficient is a unit.")

    ## comparisons
    def _richcmp_(self, other, op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(GF(5)); R
            Multivariate Power Series Ring in a, b, c over Finite Field of size 5
            sage: f = a + b + c + a^2*c
            sage: f == f^2
            False
            sage: f = f.truncate()
            sage: f == f.O(4)
            True

        Ordering is determined by underlying polynomial ring::

            sage: a > b
            True
            sage: a > a^2
            True
            sage: b > a^2
            True
            sage: (f^2).O(3)
            a^2 + 2*a*b + 2*a*c + b^2 + 2*b*c + c^2 + O(a, b, c)^3
            sage: f < f^2
            False
            sage: f > f^2
            True
            sage: f < 2*f
            True
        """
        return richcmp(self._bg_value, other._bg_value, op)

    ## arithmetic
    def _add_(left, right):
        """
        Add ``left`` to ``right``.

        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f0 = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: f1 = -6*b*c^3 - 4*a^2*b*c^2 + a^6*b^2*c - 2*a^3*b^3*c^3 + R.O(10)
            sage: g = f0 + f1; g #indirect doctest
            -6*b*c^3 - 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 + a^6*b^2*c
            - 14*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == f0.polynomial() + f1.polynomial()
            True
        """
        f = left._bg_value + right._bg_value
        return MPowerSeries(left.parent(), f, prec=f.prec())

    def _sub_(left, right):
        """
        Subtract ``right`` from ``left``.

        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f0 = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: f1 = -6*b*c^3 - 4*a^2*b*c^2 + a^6*b^2*c - 2*a^3*b^3*c^3 + R.O(10)
            sage: g = f0 - f1; g #indirect doctest
            6*b*c^3 + 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 - a^6*b^2*c
            - 10*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == f0.polynomial() - f1.polynomial()
            True
        """
        f = left._bg_value - right._bg_value
        return MPowerSeries(left.parent(), f, prec=f.prec())

    def _mul_(left, right):
        """
        Multiply ``left`` and ``right``.

        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f0 = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: f1 = -6*b*c^3 - 4*a^2*b*c^2 + a^6*b^2*c - 2*a^3*b^3*c^3 + R.O(10)
            sage: g = f0*f1; g #indirect doctest
            6*a^3*b^2*c^5 + 4*a^5*b^2*c^4 - 6*a^2*b^3*c^7 - 4*a^4*b^3*c^6
            + 72*a^3*b^4*c^6 + O(a, b, c)^14
            sage: g in R
            True

        The power series product and polynomial product agree up to
        total degree < precision of `g`::

            sage: diff = g.polynomial() - f0.polynomial() * f1.polynomial()
            sage: all(S >= g.prec() for S in [sum(e) for e in diff.exponents()])
            True
        """
        f = left._bg_value * right._bg_value
        return MPowerSeries(left.parent(), f, prec=f.prec())

    def _lmul_(self, c):
        """
        Multiply ``self`` with ``c`` on the left.

        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: g = 3*f; g #indirect doctest
            -3*a^3*b*c^2 + 3*a^2*b^2*c^4 - 36*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == 3 * (f.polynomial())
            True
            sage: g = f*5; g #indirect doctest
            -5*a^3*b*c^2 + 5*a^2*b^2*c^4 - 60*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == (f.polynomial()) * 5
            True
        """
        f = c * self._bg_value
        return MPowerSeries(self.parent(), f, prec=f.prec())

    def trailing_monomial(self):
        """
        Return the trailing monomial of ``self``.

        This is defined here as the lowest term of the underlying polynomial.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b + R.O(3)
            sage: f.trailing_monomial()
            1
            sage: f = a^2*b^3*f; f
            a^2*b^3 + a^3*b^3 + a^2*b^4 - a^3*b^4 + O(a, b, c)^8
            sage: f.trailing_monomial()
            a^2*b^3

        TESTS::

            sage: (f-f).trailing_monomial()
            0
        """
        return self.polynomial().lt()

    def quo_rem(self, other, precision=None):
        r"""
        Return the pair of quotient and remainder for the increasing power
        division of ``self`` by ``other``.

        If `a` and `b` are two elements of a power series ring
        `R[[x_1, x_2, \cdots, x_n]]` such that the trailing term of
        `b` is invertible in `R`, then the pair of quotient and
        remainder for the increasing power division of `a` by `b` is
        the unique pair `(u, v) \in R[[x_1, x_2, \cdots, x_n]] \times
        R[x_1, x_2, \cdots, x_n]` such that `a = bu + v` and such that
        no monomial appearing in `v` divides the trailing monomial
        (:meth:`trailing_monomial`) of `b`. Note that this depends on
        the order of the variables.

        This method returns both quotient and remainder as power series,
        even though in mathematics, the remainder for the increasing
        power division of two power series is a polynomial. This is
        because Sage's power series come with a precision, and that
        precision is not always sufficient to determine the remainder
        completely. Disregarding this issue, the :meth:`polynomial`
        method can be used to recast the remainder as an actual
        polynomial.

        INPUT:

        - ``other`` -- an element of the same power series ring as
          ``self`` such that the trailing term of ``other`` is
          invertible in ``self`` (this is automatically satisfied
          if the base ring is a field, unless ``other`` is zero)

        - ``precision`` -- (default: the default precision of the
          parent of ``self``) nonnegative integer, determining the
          precision to be cast on the resulting quotient and
          remainder if both ``self`` and ``other`` have infinite
          precision (ignored otherwise); note that the resulting
          precision might be lower than this integer

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b + R.O(3)
            sage: g = 1 + 2*a - 3*a*b + R.O(3)
            sage: q, r = f.quo_rem(g); q, r
            (1 - a + b + 2*a^2 + O(a, b, c)^3, 0 + O(a, b, c)^3)
            sage: f == q*g + r
            True
            sage: q, r = (a*f).quo_rem(g); q, r
            (a - a^2 + a*b + 2*a^3 + O(a, b, c)^4, 0 + O(a, b, c)^4)
            sage: a*f == q*g + r
            True
            sage: q, r = (a*f).quo_rem(a*g); q, r
            (1 - a + b + 2*a^2 + O(a, b, c)^3, 0 + O(a, b, c)^4)
            sage: a*f == q*(a*g) + r
            True
            sage: q, r = (a*f).quo_rem(b*g); q, r
            (a - 3*a^2 + O(a, b, c)^3, a + a^2 + O(a, b, c)^4)
            sage: a*f == q*(b*g) + r
            True

        Trying to divide two polynomials, we run into the issue that
        there is no natural setting for the precision of the quotient
        and remainder (and if we wouldn't set a precision, the
        algorithm would never terminate). Here, default precision
        comes to our help::

            sage: # needs sage.libs.singular
            sage: (1 + a^3).quo_rem(a + a^2)
            (a^2 - a^3 + a^4 - a^5 + a^6 - a^7 + a^8 - a^9 + a^10 + O(a, b, c)^11,
             1 + O(a, b, c)^12)
            sage: (1 + a^3 + a*b).quo_rem(b + c)
            (a + O(a, b, c)^11, 1 - a*c + a^3 + O(a, b, c)^12)
            sage: (1 + a^3 + a*b).quo_rem(b + c, precision=17)
            (a + O(a, b, c)^16, 1 - a*c + a^3 + O(a, b, c)^17)
            sage: (a^2 + b^2 + c^2).quo_rem(a + b + c)
            (a - b - c + O(a, b, c)^11, 2*b^2 + 2*b*c + 2*c^2 + O(a, b, c)^12)
            sage: (a^2 + b^2 + c^2).quo_rem(1/(1+a+b+c))
            (a^2 + b^2 + c^2 + a^3 + a^2*b + a^2*c + a*b^2 + a*c^2
               + b^3 + b^2*c + b*c^2 + c^3 + O(a, b, c)^14,
             0)
            sage: (a^2 + b^2 + c^2).quo_rem(a/(1+a+b+c))
            (a + a^2 + a*b + a*c + O(a, b, c)^13, b^2 + c^2)
            sage: (1 + a + a^15).quo_rem(a^2)
            (0 + O(a, b, c)^10, 1 + a + O(a, b, c)^12)
            sage: (1 + a + a^15).quo_rem(a^2, precision=15)
            (0 + O(a, b, c)^13, 1 + a + O(a, b, c)^15)
            sage: (1 + a + a^15).quo_rem(a^2, precision=16)
            (a^13 + O(a, b, c)^14, 1 + a + O(a, b, c)^16)

        Illustrating the dependency on the ordering of variables::

            sage: # needs sage.libs.singular
            sage: (1 + a + b).quo_rem(b + c)
            (1 + O(a, b, c)^11, 1 + a - c + O(a, b, c)^12)
            sage: (1 + b + c).quo_rem(c + a)
            (0 + O(a, b, c)^11, 1 + b + c + O(a, b, c)^12)
            sage: (1 + c + a).quo_rem(a + b)
            (1 + O(a, b, c)^11, 1 - b + c + O(a, b, c)^12)

        TESTS::

            sage: (f).quo_rem(R.zero())                                                 # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ZeroDivisionError

            sage: (f).quo_rem(R.zero().add_bigoh(2))                                    # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        Coercion is applied on ``other``::

            sage: (a + b).quo_rem(1)                                                    # needs sage.libs.singular
            (a + b + O(a, b, c)^12, 0 + O(a, b, c)^12)

            sage: R.<a,b,c> = PowerSeriesRing(QQ)
            sage: R(3).quo_rem(2)
            (3/2 + O(a, b, c)^12, 0 + O(a, b, c)^12)
        """
        parent = self.parent()
        if other.parent() is not parent:
            other = self.parent(other)
        other_tt = other.trailing_monomial()
        if not other_tt:
            raise ZeroDivisionError()
        self_prec = self.prec()
        if self_prec == infinity and other.prec() == infinity:
            if precision is None:
                precision = parent.default_prec()
            self = self.add_bigoh(precision)
            self_prec = self.prec()
        rem = parent.zero().add_bigoh(self_prec)
        quo = parent.zero().add_bigoh(self_prec-other.valuation())
        while self:
            # Loop invariants:
            # ``(the original value of self) - self == quo * other + rem``
            # and
            # ``(quo * other).prec() <= self.prec().
            # (``other`` doesn't change throughout the loop.)
            # The loop terminates because:
            # (1) every step increases ``self_tt``;
            # (2) either ``self`` has finite precision, or ``self`` is a
            #     polynomial and ``other`` has infinite precision (in
            #     which case either ``self`` will run out of nonzero
            #     coefficients after sufficiently many iterations of the
            #     if-case, or ``self``'s precision gets reduced to finite
            #     in one iteration of the else-case).
            # These show that at the end we have
            # ``(the original value of self) == quo * other + rem``
            # up to the minimum of the precision of either side of this
            # equality and the precision of self.
            self_tt = self.trailing_monomial()
            #assert self_tt
            if not other_tt.divides(self_tt):
                self -= self_tt
                rem += self_tt
            else:
                d = self_tt//other_tt
                self -= d * other
                quo += d
                quo = quo.add_bigoh(self.prec()-other_tt.degree())
        return quo, rem

    def _div_(self, denom_r):
        r"""
        Division in the ring of power series.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b + R.O(3)
            sage: g = 1/f; g #indirect doctest
            1 - a - b + a^2 + 3*a*b + b^2 + O(a, b, c)^3
            sage: g in R
            True
            sage: g == ~f
            True

        When possible, division by non-units also works::

            sage: a/(a*f)                                                               # needs sage.libs.singular
            1 - a - b + a^2 + 3*a*b + b^2 + O(a, b, c)^3

            sage: a/(R.zero())
            Traceback (most recent call last):
            ZeroDivisionError

            sage: (a*f)/f
            a + O(a, b, c)^4
            sage: f/(a*f)                                                               # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ValueError: not divisible

        An example where one loses precision::

            sage: ((1+a)*f - f) / a*f                                                   # needs sage.libs.singular
            1 + 2*a + 2*b + O(a, b, c)^2

        TESTS::

            sage: ((a+b)*f) / f == (a+b)
            True
            sage: ((a+b)*f) / (a+b) == f                                                # needs sage.libs.singular
            True
        """
        if denom_r.is_unit(): # faster if denom_r is a unit
            return self.parent(self._bg_value * ~denom_r._bg_value)
        quo, rem = self.quo_rem(denom_r)
        if rem:
            raise ValueError("not divisible")
        else:
            return quo

    def __mod__(self, other):
        """
        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: g = f % 2; g
            a^3*b*c^2 + a^2*b^2*c^4 + O(a, b, c)^10
            sage: g in R
            False
            sage: g in R.base_extend(Zmod(2))
            True
            sage: g.polynomial() == f.polynomial() % 2                                  # needs sage.libs.singular
            True
        """
        if isinstance(other, (int, Integer)):
            return self.change_ring(Zmod(other))
        raise NotImplementedError("Mod on multivariate power series ring elements not defined except modulo an integer.")

    def monomial_coefficients(self, copy=None):
        """
        Return underlying dictionary with keys the exponents and values the
        coefficients of this power series.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'t',sparse=True); M
            Sparse Multivariate Power Series Ring in t0, t1, t2, t3 over
            Rational Field

            sage: M.inject_variables()
            Defining t0, t1, t2, t3

            sage: m = 2/3*t0*t1^15*t3^48 - t0^15*t1^21*t2^28*t3^5
            sage: m2 = 1/2*t0^12*t1^29*t2^46*t3^6 - 1/4*t0^39*t1^5*t2^23*t3^30 + M.O(100)
            sage: s = m + m2
            sage: s.monomial_coefficients()
            {(1, 15, 0, 48): 2/3,
             (12, 29, 46, 6): 1/2,
             (15, 21, 28, 5): -1,
             (39, 5, 23, 30): -1/4}

        ``dict`` is an alias::

            sage: s.dict()
            {(1, 15, 0, 48): 2/3,
             (12, 29, 46, 6): 1/2,
             (15, 21, 28, 5): -1,
             (39, 5, 23, 30): -1/4}
        """
        out_dict = {}
        for j in self._bg_value.coefficients():
            out_dict.update(j.monomial_coefficients())
        return out_dict

    dict = monomial_coefficients

    def polynomial(self):
        """
        Return the underlying polynomial of ``self`` as an element of
        the underlying multivariate polynomial ring (the "foreground
        polynomial ring").

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'t'); M
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational
            Field
            sage: t = M.gens()
            sage: f = 1/2*t[0]^3*t[1]^3*t[2]^2 + 2/3*t[0]*t[2]^6*t[3] \
            - t[0]^3*t[1]^3*t[3]^3 - 1/4*t[0]*t[1]*t[2]^7 + M.O(10)
            sage: f
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7 + O(t0, t1, t2, t3)^10

            sage: f.polynomial()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7

            sage: f.polynomial().parent()
            Multivariate Polynomial Ring in t0, t1, t2, t3 over Rational Field

        Contrast with :meth:`truncate`::

            sage: f.truncate()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3 - 1/4*t0*t1*t2^7
            sage: f.truncate().parent()
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational Field
        """
        return self._value()

    def variables(self):
        """
        Return tuple of variables occurring in ``self``.

        EXAMPLES::

            sage: T = PowerSeriesRing(GF(3),5,'t'); T
            Multivariate Power Series Ring in t0, t1, t2, t3, t4 over
            Finite Field of size 3
            sage: t = T.gens()
            sage: w = t[0] - 2*t[0]*t[2] + 5*t[4]^3 - t[0]^3*t[2]^2 + T.O(6)
            sage: w
            t0 + t0*t2 - t4^3 - t0^3*t2^2 + O(t0, t1, t2, t3, t4)^6
            sage: w.variables()
            (t0, t2, t4)
        """
        return tuple(self.parent(v) for v in self._value().variables())

    def monomials(self):
        """
        Return a list of monomials of ``self``.

        These are the keys of the dict returned by :meth:`coefficients`.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b, c over Integer Ring
            sage: f = 1 + a + b - a*b - b*c - a*c + R.O(4)
            sage: sorted(f.monomials())
            [b*c, a*c, a*b, b, a, 1]
            sage: f = 1 + 2*a + 7*b - 2*a*b - 4*b*c - 13*a*c + R.O(4)
            sage: sorted(f.monomials())
            [b*c, a*c, a*b, b, a, 1]
            sage: f = R.zero()
            sage: f.monomials()
            []
        """
        return list(self.coefficients())

    def coefficients(self):
        """
        Return a dict of monomials and coefficients.

        EXAMPLES::

            sage: R.<s,t> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in s, t over Integer Ring
            sage: f = 1 + t + s + s*t + R.O(3)
            sage: f.coefficients()
            {s*t: 1, t: 1, s: 1, 1: 1}
            sage: (f^2).coefficients()
            {t^2: 1, s*t: 4, s^2: 1, t: 2, s: 2, 1: 1}

            sage: g = f^2 + f - 2; g
            3*s + 3*t + s^2 + 5*s*t + t^2 + O(s, t)^3
            sage: cd = g.coefficients()
            sage: g2 = sum(k*v for (k,v) in cd.items()); g2
            3*s + 3*t + s^2 + 5*s*t + t^2
            sage: g2 == g.truncate()
            True
        """
        if self.is_sparse():
            return self.monomial_coefficients()
        tmp = {}
        for j in self._bg_value.coefficients():
            for m in j.monomials():
                tmp[self.parent(m)] = j.monomial_coefficient(self.parent()._poly_ring(m))
        return tmp

    def constant_coefficient(self):
        """
        Return constant coefficient of ``self``.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b, c over Integer Ring
            sage: f = 3 + a + b - a*b - b*c - a*c + R.O(4)
            sage: f.constant_coefficient()
            3
            sage: f.constant_coefficient().parent()
            Integer Ring
        """
        return self.base_ring()(self._bg_value[0])

    def exponents(self):
        """
        Return a list of tuples which hold the exponents of each monomial
        of ``self``.

        EXAMPLES::

            sage: H = QQ[['x,y']]
            sage: (x,y) = H.gens()
            sage: h = -y^2 - x*y^3 - 6/5*y^6 - x^7 + 2*x^5*y^2 + H.O(10)
            sage: h
            -y^2 - x*y^3 - 6/5*y^6 - x^7 + 2*x^5*y^2 + O(x, y)^10
            sage: h.exponents()
            [(0, 2), (1, 3), (0, 6), (7, 0), (5, 2)]
        """
        exp_list = []
        for m in self._bg_value.coefficients():
            exp_list += m.exponents()
        return exp_list

    def V(self, n):
        r"""
        If

        .. MATH::

            f = \sum a_{m_0, \ldots, m_k} x_0^{m_0} \cdots x_k^{m_k},

        then this function returns

        .. MATH::

            \sum a_{m_0, \ldots, m_k} x_0^{n m_0} \cdots x_k^{n m_k}.

        The total-degree precision of the output is ``n`` times the precision
        of ``self``.

        EXAMPLES::

            sage: H = QQ[['x,y,z']]
            sage: (x,y,z) = H.gens()
            sage: h = -x*y^4*z^7 - 1/4*y*z^12 + 1/2*x^7*y^5*z^2 \
            + 2/3*y^6*z^8 + H.O(15)
            sage: h.V(3)
            -x^3*y^12*z^21 - 1/4*y^3*z^36 + 1/2*x^21*y^15*z^6 + 2/3*y^18*z^24 + O(x, y, z)^45
        """
        cd = self.coefficients()
        Vs = sum(v * k**n for k, v in cd.items())
        return Vs.add_bigoh(self.prec()*n)

    def prec(self):
        """
        Return precision of ``self``.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b, c over Integer Ring
            sage: f = 3 + a + b - a*b - b*c - a*c + R.O(4)
            sage: f.prec()
            4
            sage: f.truncate().prec()
            +Infinity
        """
        return self._prec

    def add_bigoh(self, prec):
        """
        Return a multivariate power series of precision ``prec``
        obtained by truncating ``self`` at precision ``prec``.

        This is the same as :meth:`O`.

        EXAMPLES::

            sage: B.<x,y> = PowerSeriesRing(QQ); B
            Multivariate Power Series Ring in x, y over Rational Field
            sage: r = 1 - x*y + x^2
            sage: r.add_bigoh(4)
            1 + x^2 - x*y + O(x, y)^4
            sage: r.add_bigoh(2)
            1 + O(x, y)^2

        Note that this does not change ``self``::

            sage: r
            1 + x^2 - x*y
        """
        return self.parent(self._bg_value.add_bigoh(prec))

    def O(self, prec):
        """
        Return a multivariate power series of precision ``prec``
        obtained by truncating ``self`` at precision ``prec``.

        This is the same as :meth:`add_bigoh`.

        EXAMPLES::

            sage: B.<x,y> = PowerSeriesRing(QQ); B
            Multivariate Power Series Ring in x, y over Rational Field
            sage: r = 1 - x*y + x^2
            sage: r.O(4)
            1 + x^2 - x*y + O(x, y)^4
            sage: r.O(2)
            1 + O(x, y)^2

        Note that this does not change ``self``::

            sage: r
            1 + x^2 - x*y
        """
        return self.parent(self._bg_value.O(prec))

    def truncate(self, prec=infinity):
        """
        Return infinite precision multivariate power series formed by
        truncating ``self`` at precision ``prec``.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'t'); M
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational Field
            sage: t = M.gens()
            sage: f = 1/2*t[0]^3*t[1]^3*t[2]^2 + 2/3*t[0]*t[2]^6*t[3] \
            - t[0]^3*t[1]^3*t[3]^3 - 1/4*t[0]*t[1]*t[2]^7 + M.O(10)
            sage: f
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7 + O(t0, t1, t2, t3)^10

            sage: f.truncate()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7
            sage: f.truncate().parent()
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational Field

        Contrast with polynomial::

            sage: f.polynomial()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3 - 1/4*t0*t1*t2^7
            sage: f.polynomial().parent()
            Multivariate Polynomial Ring in t0, t1, t2, t3 over Rational Field
        """
        return self.parent((self.O(prec))._value())

    def valuation(self):
        r"""
        Return the valuation of ``self``.

        The valuation of a power series `f` is the highest nonnegative
        integer `k` less or equal to the precision of `f` and such
        that the coefficient of `f` before each term of degree `< k` is
        zero. (If such an integer does not exist, then the valuation is
        the precision of `f` itself.)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<a,b> = PowerSeriesRing(GF(4949717)); R
            Multivariate Power Series Ring in a, b
             over Finite Field of size 4949717
            sage: f = a^2 + a*b + a^3 + R.O(9)
            sage: f.valuation()
            2
            sage: g = 1 + a + a^3
            sage: g.valuation()
            0
            sage: R.zero().valuation()
            +Infinity
        """
        try:
            return self._bg_value.valuation()
        except (TypeError, AttributeError):
            if self._bg_value == 0:
                return infinity

            # at this stage, self is probably a nonzero
            # element of the base ring
            for a in range(len(self._bg_value.list())):
                if self._bg_value.list()[a] != 0:
                    return a

    def is_nilpotent(self):
        """
        Return ``True`` if ``self`` is nilpotent. This occurs if

        - ``self`` has finite precision and positive valuation, or
        - ``self`` is constant and nilpotent in base ring.

        Otherwise, return ``False``.

        .. WARNING::

            This is so far just a sufficient condition, so don't trust
            a ``False`` output to be legit!

        .. TODO::

            What should we do about this method? Is nilpotency of a
            power series even decidable (assuming a nilpotency oracle
            in the base ring)? And I am not sure that returning
            ``True`` just because the series has finite precision and
            zero constant term is a good idea.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(Zmod(8)); R
            Multivariate Power Series Ring in a, b, c over Ring of integers modulo 8
            sage: f = a + b + c + a^2*c
            sage: f.is_nilpotent()
            False
            sage: f = f.O(4); f
            a + b + c + a^2*c + O(a, b, c)^4
            sage: f.is_nilpotent()
            True

            sage: g = R(2)
            sage: g.is_nilpotent()
            True
            sage: (g.O(4)).is_nilpotent()
            True

            sage: S = R.change_ring(QQ)
            sage: S(g).is_nilpotent()
            False
            sage: S(g.O(4)).is_nilpotent()
            False
        """
        if self.prec() < infinity and self.valuation() > 0:
            return True
        elif self == self.constant_coefficient() and \
           self.base_ring()(self.constant_coefficient()).is_nilpotent():
            return True
        else:
            return False

    def degree(self):
        """
        Return degree of underlying polynomial of ``self``.

        EXAMPLES::

            sage: B.<x,y> = PowerSeriesRing(QQ)
            sage: B
            Multivariate Power Series Ring in x, y over Rational Field
            sage: r = 1 - x*y + x^2
            sage: r = r.add_bigoh(4); r
            1 + x^2 - x*y + O(x, y)^4
            sage: r.degree()
            2
        """
        return self._value().degree()

    def is_unit(self):
        """
        A multivariate power series is a unit if and only if its constant
        coefficient is a unit.

        EXAMPLES::

            sage: R.<a,b> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b over Integer Ring
            sage: f = 2 + a^2 + a*b + a^3 + R.O(9)
            sage: f.is_unit()
            False
            sage: f.base_extend(QQ).is_unit()
            True
            sage: (O(a,b)^0).is_unit()
            False
        """
        # Return False for 0 + O(a, b)^0, which is the only element with precision_absolute == 0.
        if self.precision_absolute() == 0:
            return False
        return self._bg_value[0].is_unit()

    ###
    ### the following could be implemented, but aren't
    ###

    def padded_list(self):
        """
        Method from univariate power series not yet implemented.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.padded_list()
            Traceback (most recent call last):
            ...
            NotImplementedError: padded_list
        """
        raise NotImplementedError("padded_list")

    def is_square(self):
        """
        Method from univariate power series not yet implemented.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.is_square()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_square
        """
        raise NotImplementedError("is_square")

    def square_root(self):
        """
        Method from univariate power series not yet implemented.
        Depends on square root method for multivariate polynomials.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.square_root()
            Traceback (most recent call last):
            ...
            NotImplementedError: square_root
        """
        raise NotImplementedError("square_root")

    sqrt = square_root

    def derivative(self, *args):
        """
        The formal derivative of this power series, with respect to
        variables supplied in ``args``.

        EXAMPLES::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a^2*b + T.O(5)
            sage: f.derivative(a)
            1 + 2*a*b + O(a, b)^4
            sage: f.derivative(a,2)
            2*b + O(a, b)^3
            sage: f.derivative(a,a)
            2*b + O(a, b)^3
            sage: f.derivative([a,a])
            2*b + O(a, b)^3
            sage: f.derivative(a,5)
            0 + O(a, b)^0
            sage: f.derivative(a,6)
            0 + O(a, b)^0
        """
        from sage.misc.derivative import derivative_parse
        R = self.parent()
        variables = [ x.polynomial() for x in derivative_parse(args) ]
        deriv = self.polynomial().derivative(variables)
        new_prec = max(self.prec()-len(variables), 0)
        return R(deriv) + R.O(new_prec)

    def integral(self, *args):
        """
        The formal integral of this multivariate power series, with respect to
        variables supplied in ``args``.

        The variable sequence ``args`` can contain both variables and
        counts; for the syntax, see
        :meth:`~sage.misc.derivative.derivative_parse`.

        EXAMPLES::

            sage: T.<a,b> = PowerSeriesRing(QQ, 2)
            sage: f = a + b + a^2*b + T.O(5)
            sage: f.integral(a, 2)
            1/6*a^3 + 1/2*a^2*b + 1/12*a^4*b + O(a, b)^7
            sage: f.integral(a, b)
            1/2*a^2*b + 1/2*a*b^2 + 1/6*a^3*b^2 + O(a, b)^7
            sage: f.integral(a, 5)
            1/720*a^6 + 1/120*a^5*b + 1/2520*a^7*b + O(a, b)^10

        Only integration with respect to variables works::

            sage: f.integral(a + b)
            Traceback (most recent call last):
            ...
            ValueError: a + b is not a variable

        .. warning:: Coefficient division.

            If the base ring is not a field (e.g. `ZZ`), or if it has a
            nonzero characteristic, (e.g. `ZZ/3ZZ`), integration is not
            always possible while staying with the same base ring. In the
            first case, Sage will report that it has not been able to
            coerce some coefficient to the base ring::

                sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
                sage: f = a + T.O(5)
                sage: f.integral(a)
                Traceback (most recent call last):
                ...
                TypeError: no conversion of this rational to integer

            One can get the correct result by changing the base ring first::

                sage: f.change_ring(QQ).integral(a)
                1/2*a^2 + O(a, b)^6

            However, a correct result is returned even without base change
            if the denominator cancels::

                sage: f = 2*b + T.O(5)
                sage: f.integral(b)
                b^2 + O(a, b)^6

            In nonzero characteristic, Sage will report that a zero division
            occurred ::

                sage: T.<a,b> = PowerSeriesRing(Zmod(3), 2)
                sage: (a^3).integral(a)
                a^4
                sage: (a^2).integral(a)
                Traceback (most recent call last):
                ...
                ZeroDivisionError: inverse of Mod(0, 3) does not exist
        """
        from sage.misc.derivative import derivative_parse
        res = self
        for v in derivative_parse(args):
            res = res._integral(v)
        return res

    def _integral(self, xx):
        """
        Formal integral for multivariate power series.

        INPUT:

        - ``xx`` -- a generator of the power series ring (the
          one with respect to which to integrate)

        EXAMPLES::

            sage: T.<a,b> = PowerSeriesRing(QQ,2)
            sage: f = a + b + a^2*b + T.O(5)
            sage: f._integral(a)
            1/2*a^2 + a*b + 1/3*a^3*b + O(a, b)^6
            sage: f._integral(b)
            a*b + 1/2*b^2 + 1/2*a^2*b^2 + O(a, b)^6

        TESTS:

        We try to recognize variables even if they are not recognized as
        generators of the rings::

            sage: T.<a,b> = PowerSeriesRing(QQ,2)
            sage: a.is_gen()
            True
            sage: (a+0).is_gen()
            False
            sage: (a+b).integral(a+0)
            1/2*a^2 + a*b

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: aa = a.change_ring(Zmod(5))
            sage: aa.is_gen()
            False
            sage: aa.integral(aa)
            3*a^2
            sage: aa.integral(a)
            3*a^2
        """
        P = self.parent()
        R = P.base_ring()
        xx = P(xx)
        if not xx.is_gen():
            for g in P.gens():  # try to find a generator equal to xx
                if g == xx:
                    xx = g
                    break
            else:
                raise ValueError("%s is not a variable" % xx)
        xxe = xx.exponents()[0]
        pos = [i for i, c in enumerate(xxe) if c != 0][0]  # get the position of the variable
        res = {mon.eadd(xxe): R(co / (mon[pos]+1))
               for mon, co in self.monomial_coefficients().items()}
        return P( res ).add_bigoh(self.prec()+1)

    def ogf(self):
        """
        Method from univariate power series not yet implemented.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.ogf()
            Traceback (most recent call last):
            ...
            NotImplementedError: ogf
        """
        raise NotImplementedError("ogf")

    def egf(self):
        """
        Method from univariate power series not yet implemented.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.egf()
            Traceback (most recent call last):
            ...
            NotImplementedError: egf
        """
        raise NotImplementedError("egf")

    def __pari__(self):
        """
        Method from univariate power series not yet implemented.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.__pari__()
            Traceback (most recent call last):
            ...
            NotImplementedError: __pari__
        """
        raise NotImplementedError("__pari__")

    ###
    ### the following don't make sense for multivariable power series
    ###
    def list(self):
        """
        Doesn't make sense for multivariate power series.
        Multivariate polynomials don't have list of coefficients either.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.list()
            Traceback (most recent call last):
            ...
            NotImplementedError: Multivariate power series do not have list
            of coefficients; use 'coefficients' to get a dict of coefficients.
        """
        #return [self.parent(c) for c in self._bg_value.list()]
        raise NotImplementedError("Multivariate power series do not have list of coefficients; use 'coefficients' to get a dict of coefficients.")

    def variable(self):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.variable()
            Traceback (most recent call last):
            ...
            NotImplementedError: variable not defined for multivariate power
            series; use 'variables' instead.
        """
        raise NotImplementedError("variable not defined for multivariate power series; use 'variables' instead.")

    def shift(self, n):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.shift(3)
            Traceback (most recent call last):
            ...
            NotImplementedError: shift not defined for multivariate power series.
        """
        raise NotImplementedError("shift not defined for multivariate power series.")

    def __lshift__(self, n):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.__lshift__(3)
            Traceback (most recent call last):
            ...
            NotImplementedError: __lshift__ not defined for multivariate power series.
        """
        raise NotImplementedError("__lshift__ not defined for multivariate power series.")

    def __rshift__(self, n):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.__rshift__(3)
            Traceback (most recent call last):
            ...
            NotImplementedError: __rshift__ not defined for multivariate power series.
        """
        raise NotImplementedError("__rshift__ not defined for multivariate power series.")

    def valuation_zero_part(self):
        """
        Doesn't make sense for multivariate power series;
        valuation zero with respect to which variable?

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.valuation_zero_part()
            Traceback (most recent call last):
            ...
            NotImplementedError: valuation_zero_part not defined for multivariate
            power series; perhaps 'constant_coefficient' is what you want.
        """
        raise NotImplementedError("valuation_zero_part not defined for multivariate power series; perhaps 'constant_coefficient' is what you want.")

    def solve_linear_de(self, prec=infinity, b=None, f0=None):
        """
        Not implemented for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.solve_linear_de()
            Traceback (most recent call last):
            ...
            NotImplementedError: solve_linear_de not defined for multivariate power series.
        """
        raise NotImplementedError("solve_linear_de not defined for multivariate power series.")

    def exp(self, prec=infinity):
        r"""
        Exponentiate the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT:

        The exponentiated multivariate power series as a new
        multivariate power series.

        EXAMPLES::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = a + b + a*b + T.O(3)
            sage: exp(f)
            1 + a + b + 1/2*a^2 + 2*a*b + 1/2*b^2 + O(a, b)^3
            sage: f.exp()
            1 + a + b + 1/2*a^2 + 2*a*b + 1/2*b^2 + O(a, b)^3
            sage: f.exp(prec=2)
            1 + a + b + O(a, b)^2
            sage: log(exp(f)) - f
            0 + O(a, b)^3

        If the power series has a constant coefficient `c` and
        `\exp(c)` is transcendental, then `\exp(f)` would have to be a
        power series over the :class:`~sage.symbolic.ring.SymbolicRing`. These
        are not yet implemented and therefore such cases raise an error::

            sage: g = 2 + f
            sage: exp(g)                                                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Symbolic Ring' and
            'Power Series Ring in Tbg over Multivariate Polynomial Ring in a, b
            over Rational Field'

        Another workaround for this limitation is to change base ring
        to one which is closed under exponentiation, such as `\RR` or `\CC`::

            sage: exp(g.change_ring(RDF))
            7.38905609... + 7.38905609...*a + 7.38905609...*b + 3.69452804...*a^2 +
            14.7781121...*a*b + 3.69452804...*b^2 + O(a, b)^3

        If no precision is specified, the default precision is used::

            sage: T.default_prec()
            12
            sage: exp(a)
            1 + a + 1/2*a^2 + 1/6*a^3 + 1/24*a^4 + 1/120*a^5 + 1/720*a^6 + 1/5040*a^7 +
            1/40320*a^8 + 1/362880*a^9 + 1/3628800*a^10 + 1/39916800*a^11 + O(a, b)^12
            sage: a.exp(prec=5)
            1 + a + 1/2*a^2 + 1/6*a^3 + 1/24*a^4 + O(a, b)^5
            sage: exp(a + T.O(5))
            1 + a + 1/2*a^2 + 1/6*a^3 + 1/24*a^4 + O(a, b)^5

        TESTS::

            sage: exp(a^2 + T.O(5))
            1 + a^2 + 1/2*a^4 + O(a, b)^5
        """
        R = self.parent()
        Rbg = R._bg_power_series_ring

        c = self.constant_coefficient()
        if not c:
            exp_c = self.base_ring().one()
        else:
            from sage.functions.log import exp
            exp_c = exp(c)
        x = self._bg_value - c
        if x.is_zero():
            return exp_c
        val = x.valuation()
        assert (val >= 1)

        prec = min(prec, self.prec())
        if isinstance(prec, InfinityElement):
            prec = R.default_prec()
        n_inv_factorial = R.base_ring().one()
        x_pow_n = Rbg.one()
        exp_x = Rbg.one().add_bigoh(prec)
        for n in range(1,prec//val+1):
            x_pow_n = (x_pow_n * x).add_bigoh(prec)
            n_inv_factorial /= n
            exp_x += x_pow_n * n_inv_factorial
        result_bg = exp_c*exp_x

        if result_bg.base_ring() is not self.base_ring():
            R = R.change_ring(self.base_ring().fraction_field())
        return R(result_bg, prec=prec)

    def log(self, prec=infinity):
        r"""
        Return the logarithm of the formal power series.

        INPUT:

        - ``prec`` -- integer or ``infinity``; the degree to truncate
          the result to

        OUTPUT:

        The logarithm of the multivariate power series as a new
        multivariate power series.

        EXAMPLES::

            sage: T.<a,b> = PowerSeriesRing(ZZ, 2)
            sage: f = 1 + a + b + a*b + T.O(5)
            sage: f.log()
            a + b - 1/2*a^2 - 1/2*b^2 + 1/3*a^3 + 1/3*b^3 - 1/4*a^4 - 1/4*b^4 + O(a, b)^5
            sage: log(f)
            a + b - 1/2*a^2 - 1/2*b^2 + 1/3*a^3 + 1/3*b^3 - 1/4*a^4 - 1/4*b^4 + O(a, b)^5
            sage: exp(log(f)) - f
            0 + O(a, b)^5

        If the power series has a constant coefficient `c` and
        `\exp(c)` is transcendental, then `\exp(f)` would have to be a
        power series over the :class:`~sage.symbolic.ring.SymbolicRing`. These
        are not yet implemented and therefore such cases raise an error::

            sage: g = 2 + f
            sage: log(g)                                                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for -: 'Symbolic Ring' and 'Power
            Series Ring in Tbg over Multivariate Polynomial Ring in a, b over Rational Field'

        Another workaround for this limitation is to change base ring
        to one which is closed under exponentiation, such as `\RR` or `\CC`::

            sage: log(g.change_ring(RDF))
            1.09861228... + 0.333333333...*a + 0.333333333...*b - 0.0555555555...*a^2
            + 0.222222222...*a*b - 0.0555555555...*b^2 + 0.0123456790...*a^3
            - 0.0740740740...*a^2*b - 0.0740740740...*a*b^2 + 0.0123456790...*b^3
            - 0.00308641975...*a^4 + 0.0246913580...*a^3*b + 0.0246913580...*a*b^3
            - 0.00308641975...*b^4 + O(a, b)^5

        TESTS::

            sage: (1+a).log(prec=10).exp()
            1 + a + O(a, b)^10
            sage: a.exp(prec=10).log()
            a + O(a, b)^10

            sage: log(1+a)
            a - 1/2*a^2 + 1/3*a^3 - 1/4*a^4 + 1/5*a^5 - 1/6*a^6 + 1/7*a^7
            - 1/8*a^8 + 1/9*a^9 - 1/10*a^10 + 1/11*a^11 + O(a, b)^12
            sage: -log(1-a+T.O(5))
            a + 1/2*a^2 + 1/3*a^3 + 1/4*a^4 + O(a, b)^5
            sage: a.log(prec=10)
            Traceback (most recent call last):
            ...
            ValueError: Can only take formal power series for nonzero constant term.
        """
        R = self.parent()
        Rbg = R._bg_power_series_ring

        c = self.constant_coefficient()
        if c.is_zero():
            raise ValueError('Can only take formal power series for nonzero constant term.')
        if c.is_one():
            log_c = self.base_ring().zero()
        else:
            from sage.functions.log import log
            log_c = log(c)
        x = 1 - self._bg_value/c
        if x.is_zero():
            return log_c
        val = x.valuation()
        assert (val >= 1)

        prec = min(prec, self.prec())
        if isinstance(prec, InfinityElement):
            prec = R.default_prec()
        x_pow_n = Rbg.one()
        log_x = Rbg.zero().add_bigoh(prec)
        for n in range(1,prec//val+1):
            x_pow_n = (x_pow_n * x).add_bigoh(prec)
            log_x += x_pow_n / n
        result_bg = log_c - log_x

        if result_bg.base_ring() is not self.base_ring():
            R = R.change_ring(self.base_ring().fraction_field())
        return R(result_bg, prec=prec)

    def laurent_series(self):
        """
        Not implemented for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.laurent_series()
            Traceback (most recent call last):
            ...
            NotImplementedError: laurent_series not defined for multivariate power series.
        """
        raise NotImplementedError("laurent_series not defined for multivariate power series.")


class MO:
    """
    Object representing a zero element with given precision.

    EXAMPLES::

        sage: R.<u,v> = QQ[[]]
        sage: m = O(u, v)
        sage: m^4
        0 + O(u, v)^4
        sage: m^1
        0 + O(u, v)^1

        sage: T.<a,b,c> = PowerSeriesRing(ZZ, 3)
        sage: z = O(a, b, c)
        sage: z^1
        0 + O(a, b, c)^1
        sage: 1 + a + z^1
        1 + O(a, b, c)^1

        sage: w = 1 + a + O(a, b, c)^2; w
        1 + a + O(a, b, c)^2
        sage: w^2
        1 + 2*a + O(a, b, c)^2
    """
    def __init__(self, x):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R.<u,v> = QQ[[]]
            sage: m = O(u, v)
        """
        self._vars = x

    def __pow__(self, prec):
        """
        Raise ``self`` to the given precision ``prec``.

        EXAMPLES::

            sage: R.<u,v> = QQ[[]]
            sage: m = O(u, v)
            sage: m^4
            0 + O(u, v)^4
        """
        parent = self._vars[0].parent()
        if self._vars != parent.gens():
            raise NotImplementedError
        return self._vars[0].parent()(0,prec)
