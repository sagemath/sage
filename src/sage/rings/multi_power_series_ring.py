r"""
Multivariate Power Series Rings

Construct a multivariate power series ring (in finitely many variables)
over a given (commutative) base ring.

EXAMPLES:

Construct rings and elements::

    sage: R.<t,u,v> = PowerSeriesRing(QQ); R
    Multivariate Power Series Ring in t, u, v over Rational Field
    sage: TestSuite(R).run()
    sage: p = -t + 1/2*t^3*u - 1/4*t^4*u + 2/3*v^5 + R.O(6); p
    -t + 1/2*t^3*u - 1/4*t^4*u + 2/3*v^5 + O(t, u, v)^6
    sage: p in R
    True

    sage: g = 1 + v + 3*u*t^2 - 2*v^2*t^2; g
    1 + v + 3*t^2*u - 2*t^2*v^2
    sage: g in R
    True

Add big O as with single variable power series::

    sage: g.add_bigoh(3)
    1 + v + O(t, u, v)^3
    sage: g = g.O(5); g
    1 + v + 3*t^2*u - 2*t^2*v^2 + O(t, u, v)^5

Sage keeps track of total-degree precision::

    sage: f = (g-1)^2 - g + 1; f
    -v + v^2 - 3*t^2*u + 6*t^2*u*v + 2*t^2*v^2 + O(t, u, v)^5
    sage: f in R
    True
    sage: f.prec()
    5
    sage: ((g-1-v)^2).prec()
    8

Construct multivariate power series rings over various base rings.

::

    sage: M = PowerSeriesRing(QQ, 4, 'k'); M
    Multivariate Power Series Ring in k0, k1, k2, k3 over Rational Field
    sage: loads(dumps(M)) is M
    True
    sage: TestSuite(M).run()

    sage: H = PowerSeriesRing(PolynomialRing(ZZ, 3, 'z'), 4, 'f'); H
    Multivariate Power Series Ring in f0, f1, f2, f3
     over Multivariate Polynomial Ring in z0, z1, z2 over Integer Ring
    sage: TestSuite(H).run()
    sage: loads(dumps(H)) is H
    True

    sage: z = H.base_ring().gens()
    sage: f = H.gens()
    sage: h = 4*z[1]^2 + 2*z[0]*z[2] + z[1]*z[2] + z[2]^2 \
    ....: + (-z[2]^2 - 2*z[0] + z[2])*f[0]*f[2] \
    ....: + (-22*z[0]^2 + 2*z[1]^2 - z[0]*z[2] + z[2]^2 - 1955*z[2])*f[1]*f[2] \
    ....: + (-z[0]*z[1] - 2*z[1]^2)*f[2]*f[3] \
    ....: + (2*z[0]*z[1] + z[1]*z[2] - z[2]^2 - z[1] + 3*z[2])*f[3]^2 \
    ....: + H.O(3)
    sage: h in H
    True
    sage: h
    4*z1^2 + 2*z0*z2 + z1*z2 + z2^2 + (-z2^2 - 2*z0 + z2)*f0*f2
    + (-22*z0^2 + 2*z1^2 - z0*z2 + z2^2 - 1955*z2)*f1*f2
    + (-z0*z1 - 2*z1^2)*f2*f3 + (2*z0*z1 + z1*z2 - z2^2 - z1 + 3*z2)*f3^2
    + O(f0, f1, f2, f3)^3


- Use angle-bracket notation::

    sage: # needs sage.rings.finite_rings
    sage: S.<x,y> = PowerSeriesRing(GF(65537)); S
    Multivariate Power Series Ring in x, y over Finite Field of size 65537
    sage: s = -30077*x + 9485*x*y - 6260*y^3 + 12870*x^2*y^2 - 20289*y^4 + S.O(5); s
    -30077*x + 9485*x*y - 6260*y^3 + 12870*x^2*y^2 - 20289*y^4 + O(x, y)^5
    sage: s in S
    True
    sage: TestSuite(S).run()
    sage: loads(dumps(S)) is S
    True

- Use double square bracket notation::

    sage: ZZ[['s,t,u']]
    Multivariate Power Series Ring in s, t, u over Integer Ring
    sage: GF(127931)[['x,y']]                                                           # needs sage.rings.finite_rings
    Multivariate Power Series Ring in x, y over Finite Field of size 127931

Variable ordering determines how series are displayed.

::

    sage: T.<a,b> = PowerSeriesRing(ZZ,order='deglex'); T
    Multivariate Power Series Ring in a, b over Integer Ring
    sage: TestSuite(T).run()
    sage: loads(dumps(T)) is T
    True
    sage: T.term_order()
    Degree lexicographic term order
    sage: p = - 2*b^6 + a^5*b^2 + a^7 - b^2 - a*b^3 + T.O(9); p
    a^7 + a^5*b^2 - 2*b^6 - a*b^3 - b^2 + O(a, b)^9

    sage: U = PowerSeriesRing(ZZ,'a,b',order='negdeglex'); U
    Multivariate Power Series Ring in a, b over Integer Ring
    sage: U.term_order()
    Negative degree lexicographic term order
    sage: U(p)
    -b^2 - a*b^3 - 2*b^6 + a^7 + a^5*b^2 + O(a, b)^9

Change from one base ring to another::

    sage: R.<t,u,v> = PowerSeriesRing(QQ); R
    Multivariate Power Series Ring in t, u, v over Rational Field
    sage: R.base_extend(RR)                                                             # needs sage.rings.real_mpfr
    Multivariate Power Series Ring in t, u, v
     over Real Field with 53 bits of precision
    sage: R.change_ring(IntegerModRing(10))
    Multivariate Power Series Ring in t, u, v
     over Ring of integers modulo 10

    sage: S = PowerSeriesRing(GF(65537),2,'x,y'); S                                     # needs sage.rings.finite_rings
    Multivariate Power Series Ring in x, y over Finite Field of size 65537
    sage: S.change_ring(GF(5))                                                          # needs sage.rings.finite_rings
    Multivariate Power Series Ring in x, y over Finite Field of size 5

Coercion from polynomial ring::

    sage: R.<t,u,v> = PowerSeriesRing(QQ); R
    Multivariate Power Series Ring in t, u, v over Rational Field
    sage: A = PolynomialRing(ZZ,3,'t,u,v')
    sage: g = A.gens()
    sage: a = 2*g[0]*g[2] - 2*g[0] - 2; a
    2*t*v - 2*t - 2
    sage: R(a)
    -2 - 2*t + 2*t*v
    sage: R(a).O(4)
    -2 - 2*t + 2*t*v + O(t, u, v)^4
    sage: a.parent()
    Multivariate Polynomial Ring in t, u, v over Integer Ring
    sage: a in R
    True

Coercion from polynomial ring in subset of variables::

    sage: R.<t,u,v> = PowerSeriesRing(QQ); R
    Multivariate Power Series Ring in t, u, v over Rational Field
    sage: A = PolynomialRing(QQ,2,'t,v')
    sage: g = A.gens()
    sage: a = -2*g[0]*g[1] - 1/27*g[1]^2 + g[0] - 1/2*g[1]; a
    -2*t*v - 1/27*v^2 + t - 1/2*v
    sage: a in R
    True

Coercion from symbolic ring::

    sage: # needs sage.symbolic
    sage: x,y = var('x,y')
    sage: S = PowerSeriesRing(GF(11),2,'x,y'); S
    Multivariate Power Series Ring in x, y over Finite Field of size 11
    sage: type(x)
    <class 'sage.symbolic.expression.Expression'>
    sage: type(S(x))
    <class 'sage.rings.multi_power_series_ring.MPowerSeriesRing_generic_with_category.element_class'>
    sage: f = S(2/7 -100*x^2 + 1/3*x*y + y^2).O(3); f
    5 - x^2 + 4*x*y + y^2 + O(x, y)^3
    sage: f.parent()
    Multivariate Power Series Ring in x, y over Finite Field of size 11
    sage: f.parent() == S
    True

The implementation of the multivariate power series ring uses a combination
of multivariate polynomials and univariate power series. Namely, in order
to construct the multivariate power series ring `R[[x_1, x_2, \cdots, x_n]]`,
we consider the univariate power series ring `S[[T]]` over the multivariate
polynomial ring `S := R[x_1, x_2, \cdots, x_n]`, and in it we take the
subring formed by all power series whose `i`-th coefficient has degree `i`
for all `i \geq 0`. This subring is isomorphic to
`R[[x_1, x_2, \cdots, x_n]]`. This is how `R[[x_1, x_2, \cdots, x_n]]` is
implemented in this class. The ring `S` is called the foreground polynomial
ring, and the ring `S[[T]]` is called the background univariate power
series ring.

AUTHORS:

- Niles Johnson (2010-07): initial code
- Simon King (2012-08, 2013-02): Use category and coercion framework, :issue:`13412` and :issue:`14084`
"""

#*****************************************************************************
#       Copyright (C) 2010 Niles Johnson <nilesj@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.misc.latex as latex

from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import infinity
from sage.rings.multi_power_series_ring_element import MPowerSeries
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.power_series_ring import PowerSeriesRing, PowerSeriesRing_generic
from sage.structure.nonexact import Nonexact
from sage.structure.parent import Parent

from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()

from sage.categories.integral_domains import IntegralDomains
_IntegralDomains = IntegralDomains()

try:
    from sage.rings.laurent_series_ring import LaurentSeriesRing
except ImportError:
    LaurentSeriesRing = ()

lazy_import('sage.rings.lazy_series_ring', ('LazyPowerSeriesRing', 'LazyLaurentSeriesRing'))


def is_MPowerSeriesRing(x):
    """
    Return ``True`` if input is a multivariate power series ring.

    TESTS::

        sage: from sage.rings.power_series_ring import is_PowerSeriesRing
        sage: from sage.rings.multi_power_series_ring import is_MPowerSeriesRing
        sage: M = PowerSeriesRing(ZZ, 4, 'v')
        sage: is_PowerSeriesRing(M)
        doctest:warning...
        DeprecationWarning: The function is_PowerSeriesRing is deprecated;
        use 'isinstance(..., (PowerSeriesRing_generic, LazyPowerSeriesRing) and ....ngens() == 1)' instead.
        See https://github.com/sagemath/sage/issues/38290 for details.
        False
        sage: is_MPowerSeriesRing(M)
        doctest:warning...
        DeprecationWarning: The function is_MPowerSeriesRing is deprecated;
        use 'isinstance(..., (MPowerSeriesRing_generic, LazyPowerSeriesRing))' instead.
        See https://github.com/sagemath/sage/issues/38290 for details.
        True
        sage: T = PowerSeriesRing(RR, 'v')
        sage: is_PowerSeriesRing(T)
        True
        sage: is_MPowerSeriesRing(T)
        False
        sage: L = LazyPowerSeriesRing(QQ, 'x')
        sage: is_MPowerSeriesRing(L)
        True
        sage: L = LazyPowerSeriesRing(QQ, 'x, y')
        sage: is_MPowerSeriesRing(L)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38290,
                "The function is_MPowerSeriesRing is deprecated; "
                "use 'isinstance(..., (MPowerSeriesRing_generic, LazyPowerSeriesRing))' instead.")
    return isinstance(x, (MPowerSeriesRing_generic, LazyPowerSeriesRing))


class MPowerSeriesRing_generic(PowerSeriesRing_generic, Nonexact):
    r"""
    A multivariate power series ring.  This class is implemented as a
    single variable power series ring in the variable ``T`` over a
    multivariable polynomial ring in the specified generators.  Each
    generator ``g`` of the multivariable polynomial ring (called the
    "foreground ring") is mapped to ``g*T`` in the single variable power series
    ring (called the "background ring").  The background power series ring
    is used to do arithmetic and track total-degree precision.  The
    foreground polynomial ring is used to display elements.

    For usage and examples, see above, and :meth:`PowerSeriesRing`.
    """
    ### methods from PowerSeriesRing_generic that we *don't* override:
    #
    # variable_names_recursive : works just fine
    #
    # __contains__ : works just fine
    #
    # base_extend : works just fine
    #
    # is_exact : works just fine
    #
    # random_element : works just fine
    #
    # is_field : works just fine
    #
    # is_finite : works just fine
    #
    # __setitem__ : works just fine
    #
    #
    #### notes
    #
    # sparse setting may not be implemented completely
    Element = MPowerSeries

    @staticmethod
    def __classcall__(cls, base_ring, num_gens, name_list,
                 order='negdeglex', default_prec=10, sparse=False):
        """
        Preprocessing of arguments: The term order can be given as string
        or as a :class:`~sage.rings.polynomial.term_order.TermOrder` instance.

        TESTS::

            sage: P1 = PowerSeriesRing(QQ, ['f0','f1','f2','f3'], order=TermOrder('degrevlex'))
            sage: P2 = PowerSeriesRing(QQ,4,'f', order='degrevlex')
            sage: P1 is P2   # indirect doctest
            True
        """
        order = TermOrder(order, num_gens)
        return super().__classcall__(cls, base_ring, num_gens, name_list,
                                     order, default_prec, sparse)

    def __init__(self, base_ring, num_gens, name_list,
                 order='negdeglex', default_prec=10, sparse=False):
        """
        Initialize a multivariate power series ring.  See PowerSeriesRing
        for complete documentation.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``num_gens`` -- number of generators

        - ``name_list`` -- list of indeterminate names or a single name
            If a single name is given, indeterminates will be this name
            followed by a number from 0 to num_gens - 1.  If a list is
            given, these will be the indeterminate names and the length
            of the list must be equal to num_gens.

        - ``order`` -- ordering of variables; default is
          negative degree lexicographic

        - ``default_prec`` -- (default: 10) the default total-degree precision
          for elements

        - ``sparse`` -- whether or not the power series are sparse; the
          underlying polynomial ring is always sparse

        EXAMPLES::

            sage: R.<t,u,v> = PowerSeriesRing(QQ)
            sage: g = 1 + v + 3*u*t^2 - 2*v^2*t^2
            sage: g = g.add_bigoh(5); g
            1 + v + 3*t^2*u - 2*t^2*v^2 + O(t, u, v)^5
            sage: g in R
            True

        TESTS:

        By :issue:`14084`, the multi-variate power series ring belongs to the
        category of integral domains, if the base ring does::

            sage: P = ZZ[['x','y']]
            sage: P.category()
            Category of integral domains
            sage: TestSuite(P).run()

        Otherwise, it belongs to the category of commutative rings::

            sage: P = Integers(15)[['x','y']]
            sage: P.category()
            Category of commutative rings
            sage: TestSuite(P).run()
        """
        self._term_order = order
        if not base_ring.is_commutative():
            raise TypeError("Base ring must be a commutative ring.")
        n = int(num_gens)
        if n < 0:
            raise ValueError("Multivariate Polynomial Rings must have more than 0 variables.")
        self._ngens = n
        self._has_singular = False #cannot convert to Singular by default
        # Multivariate power series rings inherit from power series rings. But
        # apparently we can not call their initialisation. Instead, initialise
        # CommutativeRing and Nonexact:
        Parent.__init__(self, base=base_ring, names=name_list,
                        category=_IntegralDomains if base_ring in
                        _IntegralDomains else _CommutativeRings)
        Nonexact.__init__(self, default_prec)

        # underlying polynomial ring in which to represent elements
        self._poly_ring_ = PolynomialRing(base_ring, self.variable_names(), order=order)
        # because sometimes PowerSeriesRing_generic calls self.__poly_ring
        self._PowerSeriesRing_generic__poly_ring = self._poly_ring()

        # background univariate power series ring
        self._bg_power_series_ring = PowerSeriesRing(self._poly_ring_, 'Tbg', sparse=sparse, default_prec=default_prec)
        self._bg_indeterminate = self._bg_power_series_ring.gen()

        self._is_sparse = sparse
        self._params = (base_ring, num_gens, name_list,
                         order, default_prec, sparse)
        self._populate_coercion_lists_()

    def _repr_(self):
        """
        Print out a multivariate power series ring.

        EXAMPLES::

            sage: R.<x,y> = PowerSeriesRing(GF(17))
            sage: R  #indirect doctest
            Multivariate Power Series Ring in x, y over Finite Field of size 17
            sage: R.rename('my multivariate power series ring')
            sage: R
            my multivariate power series ring
        """
        if self.ngens() == 0:
            generators_rep = "no variables"
        else:
            generators_rep = ", ".join(self.variable_names())

        s = "Multivariate Power Series Ring in %s over %s" % (generators_rep, self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def _latex_(self):
        """
        Return latex representation of power series ring.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'v'); M
            Multivariate Power Series Ring in v0, v1, v2, v3 over Rational Field
            sage: M._latex_()
            '\\Bold{Q}[[v_{0}, v_{1}, v_{2}, v_{3}]]'
        """
        generators_latex = ", ".join(self.latex_variable_names())
        return "%s[[%s]]" % (latex.latex(self.base_ring()), generators_latex)

    def is_integral_domain(self, proof=False):
        """
        Return ``True`` if the base ring is an integral domain; otherwise
        return False.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'v'); M
            Multivariate Power Series Ring in v0, v1, v2, v3 over Rational Field
            sage: M.is_integral_domain()
            True
        """
        return self.base_ring().is_integral_domain()

    def is_noetherian(self, proof=False):
        """
        Power series over a Noetherian ring are Noetherian.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'v'); M
            Multivariate Power Series Ring in v0, v1, v2, v3 over Rational Field
            sage: M.is_noetherian()
            True

            sage: W = PowerSeriesRing(InfinitePolynomialRing(ZZ,'a'),2,'x,y')
            sage: W.is_noetherian()
            False
        """
        return self.base_ring().is_noetherian()

    def term_order(self):
        """
        Print term ordering of ``self``.  Term orderings are implemented by the
        TermOrder class.

        EXAMPLES::

            sage: M.<x,y,z> = PowerSeriesRing(ZZ,3)
            sage: M.term_order()
            Negative degree lexicographic term order
            sage: m = y*z^12 - y^6*z^8 - x^7*y^5*z^2 + x*y^2*z + M.O(15); m
            x*y^2*z + y*z^12 - x^7*y^5*z^2 - y^6*z^8 + O(x, y, z)^15

            sage: N = PowerSeriesRing(ZZ,3,'x,y,z', order='deglex')
            sage: N.term_order()
            Degree lexicographic term order
            sage: N(m)
            -x^7*y^5*z^2 - y^6*z^8 + y*z^12 + x*y^2*z + O(x, y, z)^15
        """
        return self._term_order

    def characteristic(self):
        """
        Return characteristic of base ring, which is characteristic of ``self``.

        EXAMPLES::

            sage: H = PowerSeriesRing(GF(65537),4,'f'); H                               # needs sage.rings.finite_rings
            Multivariate Power Series Ring in f0, f1, f2, f3 over
            Finite Field of size 65537
            sage: H.characteristic()                                                    # needs sage.rings.finite_rings
            65537
        """
        return self.base_ring().characteristic()

    def construction(self):
        r"""
        Return a functor `F` and base ring `R` such that ``F(R) == self``.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ, 4, 'f'); M
            Multivariate Power Series Ring in f0, f1, f2, f3 over Rational Field

            sage: (c,R) = M.construction(); (c,R)
            (Completion[('f0', 'f1', 'f2', 'f3'), prec=12],
             Multivariate Polynomial Ring in f0, f1, f2, f3 over Rational Field)
            sage: c
            Completion[('f0', 'f1', 'f2', 'f3'), prec=12]
            sage: c(R)
            Multivariate Power Series Ring in f0, f1, f2, f3 over Rational Field
            sage: c(R) == M
            True

        TESTS::

            sage: M2 = PowerSeriesRing(QQ,4,'f', sparse=True)
            sage: M == M2
            False
            sage: c,R = M2.construction()
            sage: c(R)==M2
            True
            sage: M3 = PowerSeriesRing(QQ,4,'f', order='degrevlex')
            sage: M3 == M
            False
            sage: M3 == M2
            False
            sage: c,R = M3.construction()
            sage: c(R)==M3
            True
        """
        from sage.categories.pushout import CompletionFunctor
        extras = {'order':self.term_order(), 'num_gens':self.ngens()}
        if self.is_sparse():
            extras['sparse'] = True
        return (CompletionFunctor(self._names, self.default_prec(),
                                  extras=extras),
                self._poly_ring())

    def change_ring(self, R):
        """
        Return the power series ring over `R` in the same variable as ``self``.
        This function ignores the question of whether the base ring of self
        is or can extend to the base ring of `R`; for the latter, use
        ``base_extend``.

        EXAMPLES::

            sage: R.<t,u,v> = PowerSeriesRing(QQ); R
            Multivariate Power Series Ring in t, u, v over Rational Field
            sage: R.base_extend(RR)                                                     # needs sage.rings.real_mpfr
            Multivariate Power Series Ring in t, u, v over Real Field with
            53 bits of precision
            sage: R.change_ring(IntegerModRing(10))
            Multivariate Power Series Ring in t, u, v over Ring of integers
            modulo 10
            sage: R.base_extend(IntegerModRing(10))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined


            sage: S = PowerSeriesRing(GF(65537),2,'x,y'); S                             # needs sage.rings.finite_rings
            Multivariate Power Series Ring in x, y over Finite Field of size
            65537
            sage: S.change_ring(GF(5))                                                  # needs sage.rings.finite_rings
            Multivariate Power Series Ring in x, y over Finite Field of size 5
        """
        return PowerSeriesRing(R, names=self.variable_names(), default_prec=self.default_prec())

    def remove_var(self, *var):
        """
        Remove given variable or sequence of variables from ``self``.

        EXAMPLES::

            sage: A.<s,t,u> = PowerSeriesRing(ZZ)
            sage: A.remove_var(t)
            Multivariate Power Series Ring in s, u over Integer Ring
            sage: A.remove_var(s,t)
            Power Series Ring in u over Integer Ring


            sage: M = PowerSeriesRing(GF(5),5,'t'); M
            Multivariate Power Series Ring in t0, t1, t2, t3, t4
             over Finite Field of size 5
            sage: M.remove_var(M.gens()[3])
            Multivariate Power Series Ring in t0, t1, t2, t4
             over Finite Field of size 5

        Removing all variables results in the base ring::

            sage: M.remove_var(*M.gens())
            Finite Field of size 5
        """
        vars = list(self.variable_names())
        for v in var:
            vars.remove(str(v))
        if len(vars) == 0:
            return self.base_ring()
        return PowerSeriesRing(self.base_ring(), names=vars)

    ## this is defined in PowerSeriesRing_generic
    # def __call__(self, f, prec=infinity):
    #     """
    #     Coerce object to this multivariate power series ring.
    #     """
    #     return

    def _coerce_impl(self, f):
        """
        Return the canonical coercion of ``f`` into this multivariate power
        series ring, if one is defined, or raise a :exc:`TypeError`.

        The rings that canonically coerce to this multivariate power series
        ring are:

            - this ring itself

            - a polynomial or power series ring in the same variables or a
              subset of these variables (possibly empty), over any base
              ring that canonically coerces into the base ring of this ring

        EXAMPLES::

            sage: R.<t,u,v> = PowerSeriesRing(QQ); R
            Multivariate Power Series Ring in t, u, v over Rational Field
            sage: S1.<t,v> = PolynomialRing(ZZ); S1
            Multivariate Polynomial Ring in t, v over Integer Ring
            sage: f1 = -t*v + 2*v^2 + v; f1
            -t*v + 2*v^2 + v
            sage: R(f1)
            v - t*v + 2*v^2
            sage: S2.<u,v> = PowerSeriesRing(ZZ); S2
            Multivariate Power Series Ring in u, v over Integer Ring
            sage: f2 = -2*v^2 + 5*u*v^2 + S2.O(6); f2
            -2*v^2 + 5*u*v^2 + O(u, v)^6
            sage: R(f2)
            -2*v^2 + 5*u*v^2 + O(t, u, v)^6

            sage: R2 = R.change_ring(GF(2))
            sage: R2(f1)
            v + t*v
            sage: R2(f2)
            u*v^2 + O(t, u, v)^6

        TESTS::

            sage: R.<t,u,v> = PowerSeriesRing(QQ)
            sage: S1.<t,v> = PolynomialRing(ZZ)
            sage: f1 = S1.random_element()
            sage: g1 = R._coerce_impl(f1)
            sage: f1.parent() == R
            False
            sage: g1.parent() == R
            True
        """
        P = f.parent()
        if isinstance(P, (PolynomialRing_general, MPolynomialRing_base,
                          PowerSeriesRing_generic, MPowerSeriesRing_generic,
                          LazyPowerSeriesRing)):
            if set(P.variable_names()).issubset(set(self.variable_names())):
                if self.has_coerce_map_from(P.base_ring()):
                    return self(f)
        else:
            return self(self.base_ring().coerce(f))

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        Replacement for method of PowerSeriesRing_generic.

        To be valid, a homomorphism must send generators to elements of
        positive valuation or to nilpotent elements.

        Note that the method is_nilpotent doesn't (as of sage 4.4) seem to
        be defined for obvious examples (matrices, quotients of polynomial
        rings).

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(Zmod(8)); R
            Multivariate Power Series Ring in a, b, c
             over Ring of integers modulo 8
            sage: M = PowerSeriesRing(ZZ,3,'x,y,z')
            sage: M._is_valid_homomorphism_(R,[a,c,b])
            True

            sage: M._is_valid_homomorphism_(R,[0,c,b])
            True

        2 is nilpotent in `ZZ/8`, but 3 is not::

            sage: M._is_valid_homomorphism_(R,[2,c,b])
            True
            sage: M._is_valid_homomorphism_(R,[3,c,b])
            False

        Over `ZZ`, 2 is not nilpotent::

            sage: S = R.change_ring(ZZ); S
            Multivariate Power Series Ring in a, b, c over Integer Ring
            sage: M._is_valid_homomorphism_(S,[a,c,b])
            True
            sage: M._is_valid_homomorphism_(S,[0,c,b])
            True
            sage: M._is_valid_homomorphism_(S,[2,c,b])
            False

            sage: g = [S.random_element(10)*v for v in S.gens()]
            sage: M._is_valid_homomorphism_(S,g)
            True

        You must either give a base map or there must be a coercion
        from the base ring to the codomain::

            sage: # needs sage.rings.number_field
            sage: T.<t> = ZZ[]
            sage: K.<i> = NumberField(t^2 + 1)
            sage: Q8.<z> = CyclotomicField(8)
            sage: X.<x> = PowerSeriesRing(Q8)
            sage: M.<a,b,c> = PowerSeriesRing(K)
            sage: M._is_valid_homomorphism_(X, [x,x,x+x^2])  # no coercion
            False
            sage: M._is_valid_homomorphism_(X, [x,x,x+x^2], base_map=K.hom([z^2]))
            True
        """
        try:
            im_gens = [codomain(v) for v in im_gens]
        except TypeError:
            raise TypeError("The given generator images do not coerce to codomain.")

        if len(im_gens) is not self.ngens():
            raise ValueError("You must specify the image of each generator.")
        if base_map is None and not codomain.has_coerce_map_from(self.base_ring()):
            return False
        if all(v == 0 for v in im_gens):
            return True

        if isinstance(codomain, (PowerSeriesRing_generic, MPowerSeriesRing_generic, LazyPowerSeriesRing,
                                 LaurentSeriesRing, LazyLaurentSeriesRing)):
            try:
                B = all(v.valuation() > 0 or v.is_nilpotent() for v in im_gens)
            except NotImplementedError:
                B = all(v.valuation() > 0 for v in im_gens)
            return B
        try:
            return all(v.is_nilpotent() for v in im_gens)
        except NotImplementedError:
            pass
        return False

    def _coerce_map_from_(self, P):
        """
        The rings that canonically coerce to this multivariate power series
        ring are:

        - this ring itself

        - a polynomial or power series ring in the same variables or a
          subset of these variables (possibly empty), over any base
          ring that canonically coerces into this ring

        - any ring that coerces into the foreground polynomial ring of this ring

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A = GF(17)[['x','y']]
            sage: A.has_coerce_map_from(ZZ)
            True
            sage: A.has_coerce_map_from(ZZ['x'])
            True
            sage: A.has_coerce_map_from(ZZ['y','x'])
            True
            sage: A.has_coerce_map_from(ZZ[['x']])
            True
            sage: A.has_coerce_map_from(ZZ[['y','x']])
            True
            sage: A.has_coerce_map_from(ZZ['x','z'])
            False
            sage: A.has_coerce_map_from(GF(3)['x','y'])
            False
            sage: A.has_coerce_map_from(Frac(ZZ['y','x']))
            False

        TESTS::

            sage: M = PowerSeriesRing(ZZ, 3, 'x,y,z')
            sage: M._coerce_map_from_(M)
            True
            sage: M._coerce_map_from_(M.remove_var(x))                                  # needs sage.symbolic
            True
            sage: M._coerce_map_from_(PowerSeriesRing(ZZ,x))                            # needs sage.symbolic
            True
            sage: M._coerce_map_from_(PolynomialRing(ZZ,'x,z'))
            True
            sage: M._coerce_map_from_(PolynomialRing(ZZ,0,''))
            True
            sage: M._coerce_map_from_(ZZ)
            True

            sage: M._coerce_map_from_(Zmod(13))
            False
            sage: M._coerce_map_from_(PolynomialRing(ZZ,2,'x,t'))
            False
            sage: M._coerce_map_from_(PolynomialRing(Zmod(11),2,'x,y'))
            False

            sage: P = PolynomialRing(ZZ,3,'z')
            sage: H = PowerSeriesRing(P,4,'f'); H
            Multivariate Power Series Ring in f0, f1, f2, f3 over
             Multivariate Polynomial Ring in z0, z1, z2 over Integer Ring
            sage: H._coerce_map_from_(P)
            True
            sage: H._coerce_map_from_(P.remove_var(P.gen(1)))
            True
            sage: H._coerce_map_from_(PolynomialRing(ZZ,'z2,f0'))
            True

            sage: L.<x,y> = LazyPowerSeriesRing(QQ)
            sage: R = PowerSeriesRing(QQ, names=('x','y','z'))
            sage: R.has_coerce_map_from(L)
            True
        """
        if isinstance(P, (MPolynomialRing_base, MPowerSeriesRing_generic, LazyPowerSeriesRing,
                          PolynomialRing_general, PowerSeriesRing_generic)):
            if set(P.variable_names()).issubset(set(self.variable_names())):
                if self.has_coerce_map_from(P.base_ring()):
                    return True

        return self._poly_ring().has_coerce_map_from(P)

    def _element_constructor_(self, f, prec=None):
        """
        TESTS::

            sage: M = PowerSeriesRing(ZZ,5,'t')
            sage: t = M.gens()
            sage: m = -2*t[0]*t[3]^6*t[4] - 12*t[0]^2*t[3]*t[4]^6 + t[1]*t[2]*t[3]^4*t[4]^3 + M.O(10)
            sage: M._element_constructor_(m)
            -2*t0*t3^6*t4 - 12*t0^2*t3*t4^6 + t1*t2*t3^4*t4^3 +
            O(t0, t1, t2, t3, t4)^10
            sage: R = PolynomialRing(ZZ,5,'t')
            sage: t = R.gens()
            sage: p = -4*t[0]*t[4] + t[1]^2 + t[1]*t[2] - 6*t[2]*t[4] - t[3]*t[4]
            sage: M._element_constructor_(p)
            -4*t0*t4 + t1^2 + t1*t2 - 6*t2*t4 - t3*t4
            sage: p.parent()
            Multivariate Polynomial Ring in t0, t1, t2, t3, t4 over Integer Ring
            sage: M._element_constructor_(p).parent()
            Multivariate Power Series Ring in t0, t1, t2, t3, t4 over
            Integer Ring

            sage: L.<x,y> = LazyPowerSeriesRing(QQ)
            sage: R = PowerSeriesRing(QQ, names=('x','y','z'))
            sage: R(1/(1-x-y), prec=3)
            1 + x + y + x^2 + 2*x*y + y^2 + O(x, y, z)^3
            sage: R(x + y^2)
            x + y^2
        """
        if prec is None:
            try:
                prec = f.prec()
            except AttributeError:
                prec = infinity
        from sage.rings.lazy_series import LazyPowerSeries
        if isinstance(f, LazyPowerSeries):
            if prec is infinity:
                try:
                    f = f.polynomial()
                except ValueError:
                    f = f.add_bigoh(self.default_prec())
            else:
                f = f.add_bigoh(prec)
        return self.element_class(parent=self, x=f, prec=prec)

    def laurent_series_ring(self):
        """
        Laurent series not yet implemented for multivariate power series rings.

        TESTS::

            sage: M = PowerSeriesRing(ZZ,3,'x,y,z')
            sage: M.laurent_series_ring()
            Traceback (most recent call last):
            ...
            NotImplementedError: Laurent series not implemented for
            multivariate power series.
        """
        raise NotImplementedError("Laurent series not implemented for multivariate power series.")

    def _poly_ring(self, x=None):
        """
        Return the underlying polynomial ring used to represent elements of
        this power series ring.  If given an input x, returns x coerced
        into this polynomial ring.

        EXAMPLES::

            sage: R.<t,u> = PowerSeriesRing(QQ)
            sage: R._poly_ring()
            Multivariate Polynomial Ring in t, u over Rational Field
            sage: R._poly_ring(2).parent()
            Multivariate Polynomial Ring in t, u over Rational Field
        """
        if x is None:
            return self._poly_ring_
        else:
            return self._poly_ring_(x)

    def _mpoly_ring(self, x=None):
        """
        Same as ``_poly_ring``.

        TESTS::

            sage: R.<t,u> = PowerSeriesRing(QQ)
            sage: R._mpoly_ring()
            Multivariate Polynomial Ring in t, u over Rational Field
            sage: R._mpoly_ring(2).parent()
            Multivariate Polynomial Ring in t, u over Rational Field
        """
        return self._poly_ring(x)

    def _bg_ps_ring(self, x=None):
        """
        Return the background univariate power series ring.  If given an
        input x, returns x coerced into this power series ring.

        EXAMPLES::

            sage: R.<t,u> = PowerSeriesRing(QQ)
            sage: R._bg_ps_ring()
            Power Series Ring in Tbg over Multivariate Polynomial Ring in
            t, u over Rational Field
            sage: R._bg_ps_ring(4).parent() == R
            False
        """
        if x is None:
            return self._bg_power_series_ring
        else:
            return self._bg_power_series_ring(x)

    def is_sparse(self):
        """
        Check whether ``self`` is sparse.

        EXAMPLES::

            sage: M = PowerSeriesRing(ZZ, 3, 's,t,u'); M
            Multivariate Power Series Ring in s, t, u over Integer Ring
            sage: M.is_sparse()
            False
            sage: N = PowerSeriesRing(ZZ, 3, 's,t,u', sparse=True); N
            Sparse Multivariate Power Series Ring in s, t, u over Integer Ring
            sage: N.is_sparse()
            True
        """
        return self._is_sparse

    def is_dense(self):
        """
        Is ``self`` dense? (opposite of sparse)

        EXAMPLES::

            sage: M = PowerSeriesRing(ZZ, 3, 's,t,u'); M
            Multivariate Power Series Ring in s, t, u over Integer Ring
            sage: M.is_dense()
            True
            sage: N = PowerSeriesRing(ZZ, 3, 's,t,u', sparse=True); N
            Sparse Multivariate Power Series Ring in s, t, u over Integer Ring
            sage: N.is_dense()
            False
        """
        return not self.is_sparse()

    def gen(self, n=0):
        """
        Return the `n`-th generator of ``self``.

        EXAMPLES::

            sage: M = PowerSeriesRing(ZZ, 10, 'v')
            sage: M.gen(6)
            v6
        """
        if n < 0 or n >= self._ngens:
            raise ValueError("Generator not defined.")
        #return self(self._poly_ring().gens()[int(n)])
        return self.element_class(parent=self, x=self._poly_ring().gens()[int(n)], is_gen=True)

    def ngens(self):
        """
        Return number of generators of ``self``.

        EXAMPLES::

            sage: M = PowerSeriesRing(ZZ, 10, 'v')
            sage: M.ngens()
            10
        """
        return self._ngens

    def gens(self) -> tuple:
        """
        Return the generators of this ring.

        EXAMPLES::

            sage: M = PowerSeriesRing(ZZ, 3, 'v')
            sage: M.gens()
            (v0, v1, v2)
        """
        return tuple(self.gen(i) for i in range(self._ngens))

    def prec_ideal(self):
        """
        Return the ideal which determines precision; this is the ideal
        generated by all of the generators of our background polynomial
        ring.

        EXAMPLES::

            sage: A.<s,t,u> = PowerSeriesRing(ZZ)
            sage: A.prec_ideal()
            Ideal (s, t, u) of
             Multivariate Polynomial Ring in s, t, u over Integer Ring
        """
        return self._poly_ring().ideal(self._poly_ring().gens())

    def bigoh(self, prec):
        """
        Return big oh with precision ``prec``.  The function ``O`` does the same thing.

        EXAMPLES::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2); T
            Multivariate Power Series Ring in a, b over Integer Ring
            sage: T.bigoh(10)
            0 + O(a, b)^10
            sage: T.O(10)
            0 + O(a, b)^10
        """
        return self.zero().O(prec)

    def O(self, prec):
        """
        Return big oh with precision ``prec``.  This function is an alias for ``bigoh``.

        EXAMPLES::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2); T
            Multivariate Power Series Ring in a, b over Integer Ring
            sage: T.O(10)
            0 + O(a, b)^10
            sage: T.bigoh(10)
            0 + O(a, b)^10
        """
        return self.bigoh(prec)

    def _send_to_bg(self, f):
        """
        Send an element of the foreground polynomial ring to the background
        power series ring.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'f')
            sage: f = M._poly_ring().gens()
            sage: bg = M._send_to_bg((f[0] + f[2] + 2)**2); bg
            4 + (4*f0 + 4*f2)*Tbg + (f0^2 + 2*f0*f2 + f2^2)*Tbg^2

            sage: M._send_to_bg(bg)
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce input to polynomial ring.
        """
        try:
            f = self._poly_ring(f)
        except TypeError:
            raise TypeError("Cannot coerce input to polynomial ring.")
        return self._bg_ps_ring(f.homogeneous_components())

    def _send_to_fg(self, f):
        """
        Send an element of the background univariate power series ring to
        the foreground multivariate polynomial ring.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'f')
            sage: f = M._poly_ring().gens()
            sage: bg = M._send_to_bg((f[0] + f[2] + 2)**2); bg
            4 + (4*f0 + 4*f2)*Tbg + (f0^2 + 2*f0*f2 + f2^2)*Tbg^2
            sage: bg.parent()
            Power Series Ring in Tbg over Multivariate Polynomial Ring in f0, f1,
            f2, f3 over Rational Field
            sage: fg = M._send_to_fg(bg); fg
            4 + 4*f0 + 4*f2 + f0^2 + 2*f0*f2 + f2^2
            sage: fg.parent()
            Multivariate Polynomial Ring in f0, f1, f2, f3 over Rational Field
            sage: fg = M._send_to_fg(bg.add_bigoh(3)); fg
            4 + 4*f0 + 4*f2 + f0^2 + 2*f0*f2 + f2^2
            sage: fg = M._send_to_fg(bg.add_bigoh(2)); fg
            4 + 4*f0 + 4*f2
        """
        return self._poly_ring_.sum(f.polynomial().coefficients())


def unpickle_multi_power_series_ring_v0(base_ring, num_gens, names, order, default_prec, sparse):
    """
    Unpickle (deserialize) a multivariate power series ring according
    to the given inputs.

    EXAMPLES::

        sage: P.<x,y> = PowerSeriesRing(QQ)
        sage: loads(dumps(P)) == P # indirect doctest
        True
    """
    return PowerSeriesRing(base_ring, num_gens=num_gens, names=names, order=order, default_prec=default_prec, sparse=sparse)
