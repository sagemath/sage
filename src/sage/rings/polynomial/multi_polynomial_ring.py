r"""
Multivariate Polynomial Rings over Generic Rings

Sage implements multivariate polynomial rings through several
backends. This generic implementation uses the classes ``PolyDict``
and ``ETuple`` to construct a dictionary with exponent tuples as keys
and coefficients as values.

AUTHORS:

- David Joyner and William Stein

- Kiran S. Kedlaya (2006-02-12): added Macaulay2 analogues of Singular
  features

- Martin Albrecht (2006-04-21): reorganize class hierarchy for singular
  rep

- Martin Albrecht (2007-04-20): reorganized class hierarchy to support
  Pyrex implementations

- Robert Bradshaw (2007-08-15): Coercions from rings in a subset of
  the variables.

EXAMPLES:

We construct the Frobenius morphism on `\GF{5}[x,y,z]` over
`\GF{5}`::

    sage: R.<x,y,z> = GF(5)[]
    sage: frob = R.hom([x^5, y^5, z^5])
    sage: frob(x^2 + 2*y - z^4)
    -z^20 + x^10 + 2*y^5
    sage: frob((x + 2*y)^3)                                                             # needs sage.rings.finite_rings
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15
    sage: (x^5 + 2*y^5)^3                                                               # needs sage.rings.finite_rings
    x^15 + x^10*y^5 + 2*x^5*y^10 - 2*y^15

We make a polynomial ring in one variable over a polynomial ring in
two variables::

    sage: R.<x, y> = PolynomialRing(QQ, 2)
    sage: S.<t> = PowerSeriesRing(R)
    sage: t*(x+y)
    (x + y)*t

TESTS::

    sage: PolynomialRing(GF(5), 3, 'xyz').objgens()
    (Multivariate Polynomial Ring in x, y, z over Finite Field of size 5,
    (x, y, z))
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.ring import IntegralDomain
import sage.rings.fraction_field_element as fraction_field_element

from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base, is_MPolynomialRing
from sage.rings.polynomial.polynomial_singular_interface import PolynomialRing_singular_repr
from sage.rings.polynomial.polydict import PolyDict, ETuple
from sage.rings.polynomial.term_order import TermOrder

import sage.interfaces.abc

try:
    from sage.libs.pari.all import pari_gen
except ImportError:
    pari_gen = ()

from sage.structure.element import Element


class MPolynomialRing_macaulay2_repr:
    """
    A mixin class for polynomial rings that support conversion to Macaulay2.
    """
    def _macaulay2_init_(self, macaulay2=None):
        """
        EXAMPLES::

            sage: PolynomialRing(QQ, 'x', 2, implementation='generic')._macaulay2_init_()   # optional - macaulay2
            'sage...[symbol x0,symbol  x1, MonomialSize=>16, MonomialOrder=>GRevLex]'
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default
        return macaulay2._macaulay2_input_ring(self.base_ring(), self.gens(),
                                               self.term_order().macaulay2_str())


class MPolynomialRing_polydict(MPolynomialRing_macaulay2_repr, PolynomialRing_singular_repr, MPolynomialRing_base):
    """
    Multivariable polynomial ring.

    EXAMPLES::

        sage: R = PolynomialRing(Integers(12), 'x', 5); R
        Multivariate Polynomial Ring in x0, x1, x2, x3, x4 over Ring of integers modulo 12
        sage: loads(R.dumps()) == R
        True
    """
    from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict as Element_hidden
    # should be just Element, once polynomial use new coercion framework

    def __init__(self, base_ring, n, names, order):
        from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular
        order = TermOrder(order, n)
        # MPolynomialRing_base.__init__() normally initialises the base ring,
        # but it also needs the generators to construct a coercion map from the
        # base ring, and the base ring must be set to initialise the generators.
        # We set the base ring manually to break this circular dependency.
        self._base = base_ring
        # Construct the generators
        v = [0] * n
        one = base_ring.one()
        self._gens = []
        for i in range(n):
            v[i] = 1  # int's!
            self._gens.append(self.Element_hidden(self, {tuple(v): one}))
            v[i] = 0
        self._gens = tuple(self._gens)
        self._zero_tuple = tuple(v)
        MPolynomialRing_base.__init__(self, base_ring, n, names, order)
        self._has_singular = can_convert_to_singular(self)

    def _monomial_order_function(self):
        return self.__monomial_order_function

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: R = PolynomialRing(Integers(10), 'x', 4)
            sage: loads(R.dumps()) == R
            True
        """
        if not isinstance(other, MPolynomialRing_base):
            return False
        return ((self.base_ring(), self.ngens(),
                self.variable_names(), self.term_order()) ==
                (other.base_ring(), other.ngens(),
                 other.variable_names(), other.term_order()))

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: R = PolynomialRing(Integers(8), 'x', 3)
            sage: loads(R.dumps()) != R
            False
        """
        return not (self == other)

    def __hash__(self):
        """
        Compute the hash of ``self``.

        EXAMPLES::

            sage: h = hash(PolynomialRing(Integers(8), 'x', 3))
        """
        return hash((self.base_ring(), self.ngens(),
                     self.variable_names(), self.term_order()))

    def __call__(self, x=0, check=True):
        """
        Convert ``x`` to an element of this multivariate polynomial ring,
        possibly non-canonically.

        EXAMPLES:

        We create a Macaulay2 multivariate polynomial via ideal
        arithmetic, then convert it into R.

        ::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: I = R.ideal([x^3 + y, y])
            sage: S = I._macaulay2_()                                    # optional - macaulay2
            sage: T = S*S*S                                              # optional - macaulay2
            sage: U = T.gens().entries().flatten()                       # optional - macaulay2
            sage: f = U[2]; f                                            # optional - macaulay2
             6      3 2    3
            x y + 2x y  + y
            sage: R(f.external_string())                                 # optional - macaulay2
            x^6*y + 2*x^3*y^2 + y^3

        Some other subtle conversions. We create polynomial rings in 2
        variables over the rationals, integers, and a finite field.

        ::

            sage: R.<x,y> = QQ[]
            sage: S.<x,y> = ZZ[]
            sage: T.<x,y> = GF(7)[]

        We convert from integer polynomials to rational polynomials,
        and back::

            sage: f = R(S.0^2 - 4*S.1^3); f
            -4*y^3 + x^2
            sage: parent(f)
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: parent(S(f))
            Multivariate Polynomial Ring in x, y over Integer Ring

        We convert from polynomials over the finite field.

        ::

            sage: f = R(T.0^2 - 4*T.1^3); f
            3*y^3 + x^2
            sage: parent(f)
            Multivariate Polynomial Ring in x, y over Rational Field

        We dump and load the polynomial ring S::

            sage: S2 = loads(dumps(S))
            sage: S2 == S
            True

        Coerce works and gets the right parent.

        ::

            sage: parent(S2.coerce(S.0)) is S2
            True

        Conversion to reduce modulo a prime between rings with different
        variable names::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = PolynomialRing(GF(7),2)
            sage: f = x^2 + 2/3*y^3
            sage: S(f)
            3*b^3 + a^2

        Conversion from symbolic variables::

            sage: # needs sage.symbolic
            sage: x,y,z = var('x,y,z')
            sage: R = QQ['x,y,z']
            sage: type(x)
            <class 'sage.symbolic.expression.Expression'>
            sage: type(R(x))
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f = R(x^3 + y^3 - z^3); f
            x^3 + y^3 - z^3
            sage: type(f)
            <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: parent(f)
            Multivariate Polynomial Ring in x, y, z over Rational Field

        A more complicated symbolic and computational mix. Behind the
        scenes Singular and Maxima are doing the real work.

        ::

            sage: # needs sage.symbolic
            sage: R = QQ['x,y,z']
            sage: f = (x^3 + y^3 - z^3)^10; f
            (x^3 + y^3 - z^3)^10
            sage: g = R(f); parent(g)
            Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: (f - g).expand()
            0

        It intelligently handles conversions from polynomial rings in a subset
        of the variables too.

        ::

            sage: R = GF(5)['x,y,z']
            sage: S = ZZ['y']
            sage: R(7*S.0)
            2*y
            sage: T = ZZ['x,z']
            sage: R(2*T.0 + 6*T.1 + T.0*T.1^2)
            x*z^2 + 2*x + z

        ::

            sage: R = QQ['t,x,y,z']
            sage: S.<x> = ZZ['x']
            sage: T.<z> = S['z']
            sage: T
            Univariate Polynomial Ring in z over Univariate Polynomial Ring in x over Integer Ring
            sage: f = (x+3*z+5)^2; f
            9*z^2 + (6*x + 30)*z + x^2 + 10*x + 25
            sage: R(f)
            x^2 + 6*x*z + 9*z^2 + 10*x + 30*z + 25

        Arithmetic with a constant from a base ring::

            sage: R.<u,v> = QQ[]
            sage: S.<x,y> = R[]
            sage: u^3*x^2 + v*y
            u^3*x^2 + v*y

        Stacked polynomial rings convert into constants if possible. First,
        the univariate case::

            sage: R.<x> = QQ[]
            sage: S.<u,v> = R[]
            sage: S(u + 2)
            u + 2
            sage: S(u + 2).degree()
            1
            sage: S(x + 3)
            x + 3
            sage: S(x + 3).degree()
            0

        Second, the multivariate case::

            sage: R.<x,y> = QQ[]
            sage: S.<u,v> = R[]
            sage: S(x + 2*y)
            x + 2*y
            sage: S(u + 2*v)
            u + 2*v

        Conversion from strings::

            sage: R.<x,y> = QQ[]
            sage: R('x+(1/2)*y^2')
            1/2*y^2 + x
            sage: S.<u,v> = ZZ[]
            sage: S('u^2 + u*v + v^2')
            u^2 + u*v + v^2

        Foreign polynomial rings convert into the highest ring; the point
        here is that an element of T could convert to an element of R or an
        element of S; it is anticipated that an element of T is more likely
        to be "the right thing" and is historically consistent.

        ::

            sage: R.<x,y> = QQ[]
            sage: S.<u,v> = R[]
            sage: T.<a,b> = QQ[]
            sage: S(a + b)
            u + v

        TESTS:

        Check if we still allow nonsense (see :issue:`7951`)::

            sage: P = PolynomialRing(QQ, 0, '')
            sage: P('pi')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert pi to a rational

        Check that it is possible to convert strings to iterated
        polynomial rings (see :issue:`13327`)::

            sage: Rm = QQ["a"]["b, c"]
            sage: Rm("a*b")
            a*b
            sage: parent(_) is Rm
            True

        Check that conversion from PARI works correctly (see
        :issue:`17974`)::

            sage: A.<a> = PolynomialRing(QQ)
            sage: B.<d,e> = PolynomialRing(A)
            sage: f = pari(a*d)                                                         # needs sage.libs.pari
            sage: B(f)
            a*d
            sage: f = pari(a*d - (a+1)*d*e^3 + a*d^2)                                   # needs sage.libs.pari
            sage: B(f)
            (-a - 1)*d*e^3 + a*d^2 + a*d

            sage: A.<a,b> = PolynomialRing(QQ)
            sage: B.<d,e> = PolynomialRing(A)
            sage: f = pari(a*d)                                                         # needs sage.libs.pari
            sage: B(f)
            a*d

        It is possible to convert `f` into `B` by using ``f.sage()``,
        but this requires specifying a ``locals`` argument::

            sage: f
            d*a
            sage: f.sage()
            Traceback (most recent call last):
            ...
            NameError: name 'd' is not defined
            sage: f.sage(locals={'a': a, 'd': d})
            a*d

        Check that :issue:`21999` is fixed::

            sage: R = QQbar['s,t']                                                      # needs sage.rings.number_field
            sage: type(R({(1,2): 3}).coefficients()[0])
            <class 'sage.rings.qqbar.AlgebraicNumber'>
        """
        from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
        import sage.rings.polynomial.polynomial_element as polynomial_element

        # handle constants that coerce into self.base_ring() first, if possible
        if isinstance(x, Element) and x.parent() is self.base_ring():
            # A Constant multi-polynomial
            return self({self._zero_tuple: x})

        try:
            y = self.base_ring().coerce(x)
            return MPolynomial_polydict(self, {self._zero_tuple: y})
        except TypeError:
            pass

        from .multi_polynomial import MPolynomial_libsingular

        if isinstance(x, MPolynomial_polydict):
            P = x.parent()

            if P is self:
                return x

            if P == self:
                return MPolynomial_polydict(self, x.element().dict())

            if self.base_ring().has_coerce_map_from(P):
                # it might be in the base ring (i.e. a poly ring over a poly ring)
                c = self.base_ring()(x)
                return MPolynomial_polydict(self, {self._zero_tuple: c})

            if len(P.variable_names()) == len(self.variable_names()):
                # Map the variables in some crazy way (but in order,
                # of course).  This is here since R(blah) is supposed
                # to be "make an element of R if at all possible with
                # no guarantees that this is mathematically solid."
                K = self.base_ring()
                D = x.element().dict()
                for i, a in D.items():
                    D[i] = K(a)
                return MPolynomial_polydict(self, D)

            if (set(P.variable_names()).issubset(set(self.variable_names()))
                and self.base_ring().has_coerce_map_from(P.base_ring())):
                # If the named variables are a superset of the input, map the variables by name
                return MPolynomial_polydict(self, self._extract_polydict(x))

            return MPolynomial_polydict(self,
                                        x._mpoly_dict_recursive(self.variable_names(),
                                                                self.base_ring()))

        elif isinstance(x, MPolynomial_libsingular):
            P = x.parent()
            if P == self:
                return MPolynomial_polydict(self, x.monomial_coefficients())

            if self.base_ring().has_coerce_map_from(P):
                # it might be in the base ring (i.e. a poly ring over a poly ring)
                c = self.base_ring()(x)
                return MPolynomial_polydict(self, {self._zero_tuple: c})

            if len(P.variable_names()) == len(self.variable_names()):
                # Map the variables in some crazy way (but in order,
                # of course).  This is here since R(blah) is supposed
                # to be "make an element of R if at all possible with
                # no guarantees that this is mathematically solid."
                K = self.base_ring()
                D = x.monomial_coefficients()
                for i, a in D.items():
                    D[i] = K(a)
                return MPolynomial_polydict(self, D)

            if (set(P.variable_names()).issubset(set(self.variable_names()))
                and self.base_ring().has_coerce_map_from(P.base_ring())):
                # If the named variables are a superset of the input, map the variables by name
                return MPolynomial_polydict(self, self._extract_polydict(x))

            return MPolynomial_polydict(self,
                                        x._mpoly_dict_recursive(self.variable_names(),
                                                                self.base_ring()))

        if isinstance(x, polynomial_element.Polynomial):
            return MPolynomial_polydict(self,
                                        x._mpoly_dict_recursive(self.variable_names(),
                                                                self.base_ring()))

        if isinstance(x, PolyDict):
            return MPolynomial_polydict(self, x)

        if isinstance(x, dict):
            K = self.base_ring()
            return MPolynomial_polydict(self, {i: K(a) for i, a in x.items()})

        if (isinstance(x, fraction_field_element.FractionFieldElement)
            and x.parent().ring() == self):
            if x.denominator() == 1:
                return x.numerator()

            raise TypeError("unable to coerce since the denominator is not 1")

        if isinstance(x, sage.interfaces.abc.SingularElement) and self._has_singular:
            self._singular_().set_ring()
            try:
                return x.sage_poly(self)
            except TypeError:
                raise TypeError("unable to coerce singular object")

        if hasattr(x, '_polynomial_'):
            return x._polynomial_(self)

        if isinstance(x, str):
            from sage.misc.sage_eval import sage_eval
            try:
                x = sage_eval(x, self.gens_dict_recursive())
            except NameError:
                raise TypeError("unable to evaluate {!r} in {}".format(x, self))
            return self(x)

        if isinstance(x, sage.interfaces.abc.Macaulay2Element):
            try:
                s = x.sage_polystring()
                if len(s) == 0:
                    raise TypeError
                # NOTE: It's CRUCIAL to use the eval command as follows,
                # i.e., with the gen dict as the third arg and the second
                # empty.  Otherwise pickling won't work after calls to this eval!!!
                # This took a while to figure out!
                return self(eval(s, {}, self.gens_dict()))
            except (AttributeError, TypeError, NameError, SyntaxError):
                raise TypeError("Unable to coerce macaulay2 object")
            return MPolynomial_polydict(self, x)

        if isinstance(x, pari_gen) and x.type() == 't_POL':
            # This recursive approach is needed because PARI
            # represents multivariate polynomials as iterated
            # univariate polynomials.  Below, v is the variable
            # with highest priority, and the x[i] are expressions
            # in the remaining variables.
            d = x.poldegree()
            if d.type() == 't_INFINITY':
                return self.zero()
            v = self.gens_dict_recursive()[str(x.variable())]
            return sum(self(x[i]) * v ** i for i in range(d + 1))

        if isinstance(x, dict):
            return MPolynomial_polydict(self, x)

        c = self.base_ring()(x)
        return MPolynomial_polydict(self, {self._zero_tuple: c})

# The following methods are handy for implementing Groebner
# basis algorithms. They do only superficial type/sanity checks
# and should be called carefully.

    def monomial_quotient(self, f, g, coeff=False):
        r"""
        Return ``f/g``, where both ``f`` and ``g`` are treated as monomials.

        Coefficients are ignored by default.

        INPUT:

        - ``f`` -- monomial

        - ``g`` -- monomial

        - ``coeff`` -- divide coefficients as well (default: ``False``)

        OUTPUT: monomial

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(QQ, 3, order='degrevlex')
            sage: P.monomial_quotient(3/2*x*y, x)
            y

        ::

            sage: P.monomial_quotient(3/2*x*y, 2*x, coeff=True)
            3/4*y

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: R.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_quotient(x*y, x)
            y

        ::

            sage: P.monomial_quotient(x*y, R.gen())
            y

        ::

            sage: P.monomial_quotient(P(0), P(1))
            0

        ::

            sage: P.monomial_quotient(P(1), P(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        ::

            sage: P.monomial_quotient(P(3/2), P(2/3), coeff=True)
            9/4

        ::

            sage: P.monomial_quotient(x, y) # Note the wrong result
            x*y^-1

        ::

            sage: P.monomial_quotient(x, P(1))
            x

        .. NOTE::

            Assumes that the head term of f is a multiple of the head
            term of g and return the multiplicant m. If this rule is
            violated, funny things may happen.
        """
        from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict

        if not f:
            return f
        if not g:
            raise ZeroDivisionError

        fd = f.monomial_coefficients()
        gd = g.monomial_coefficients()

        if not coeff:
            f = next(iter(fd))
            g = next(iter(gd))
            coeff = self.base_ring().one()
        else:
            f, cf = next(iter(fd.items()))
            g, cg = next(iter(gd.items()))
            coeff = self.base_ring()(cf / cg)

        res = f.esub(g)

        return MPolynomial_polydict(self, PolyDict({res: coeff}))

    def monomial_lcm(self, f, g):
        """
        LCM for monomials. Coefficients are ignored.

        INPUT:

        - ``f`` -- monomial

        - ``g`` -- monomial

        OUTPUT: monomial

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_lcm(3/2*x*y, x)
            x*y

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: R.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_lcm(x*y, R.gen())
            x*y

        ::

            sage: P.monomial_lcm(P(3/2), P(2/3))
            1

        ::

            sage: P.monomial_lcm(x, P(1))
            x
        """
        one = self.base_ring().one()

        f, = f.monomial_coefficients()
        g, = g.monomial_coefficients()

        length = len(f)

        res = {i: max(f[i], g[i])
               for i in f.common_nonzero_positions(g)}

        return self(PolyDict({ETuple(res, length): one}))

    def monomial_reduce(self, f, G):
        r"""
        Try to find a ``g`` in ``G`` where ``g.lm()`` divides ``f``.

        If found, ``(flt,g)`` is returned, ``(0,0)`` otherwise, where
        ``flt`` is ``f/g.lm()``. It is assumed that ``G`` is iterable and contains
        ONLY elements in this ring.

        INPUT:

        - ``f`` -- monomial

        - ``G`` -- list/set of mpolynomials

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: f = x*y^2
            sage: G = [3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, P(1/2)]
            sage: P.monomial_reduce(f,G)
            (y, 1/4*x*y + 2/7)

        ::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(Zmod(23432),3, order='degrevlex')
            sage: f = x*y^2
            sage: G = [3*x^3 + y^2 + 2, 4*x*y + 7, P(2)]
            sage: P.monomial_reduce(f,G)
            (y, 4*x*y + 7)

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z>=MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: f = x*y^2
            sage: G = [3/2*x^3 + y^2 + 1/2, 1/4*x*y + 2/7, P(1/2)]

        ::

            sage: P.monomial_reduce(P(0),G)
            (0, 0)

        ::

            sage: P.monomial_reduce(f,[P(0)])
            (0, 0)
        """
        if not f:
            return 0, 0
        for g in G:
            t = g.lm()
            try:
                if self.monomial_divides(t, f):
                    return self.monomial_quotient(f, t), g
            except ZeroDivisionError:
                return 0, 0
        return 0, 0

    def monomial_divides(self, a, b):
        """
        Return ``False`` if ``a`` does not divide ``b`` and ``True`` otherwise.

        INPUT:

        - ``a`` -- monomial

        - ``b`` -- monomial

        OUTPUT: boolean

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(ZZ,3, order='degrevlex')
            sage: P.monomial_divides(x*y*z, x^3*y^2*z^4)
            True
            sage: P.monomial_divides(x^3*y^2*z^4, x*y*z)
            False

        TESTS::

            sage: P.<x,y,z> = PolynomialRing(ZZ,3, order='degrevlex')
            sage: P.monomial_divides(P(1), P(0))
            True
            sage: P.monomial_divides(P(1), x)
            True
        """
        if not b:
            return True
        if not a:
            raise ZeroDivisionError

        a, = a.monomial_coefficients()
        b, = b.monomial_coefficients()

        return all(b[i] >= a[i]
                   for i in b.common_nonzero_positions(a))

    def monomial_pairwise_prime(self, h, g):
        r"""
        Return ``True`` if ``h`` and ``g`` are pairwise prime.

        Both are treated as monomials.

        INPUT:

        - ``h`` -- monomial

        - ``g`` -- monomial

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_pairwise_prime(x^2*z^3, y^4)
            True

        ::

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, 3/4*y^3)
            False

        TESTS::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: Q.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_pairwise_prime(x^2*z^3, Q('y^4'))
            True

        ::

            sage: P.monomial_pairwise_prime(1/2*x^3*y^2, Q(0))
            True

        ::

            sage: P.monomial_pairwise_prime(P(1/2),x)
            False
        """
        not_g = not g
        not_h = not h
        if not_g and not_h:  # GCD(0,0) = 0
            return False
        if not_g or not_h:  # GCD(x,0) = 1
            return True
        return self.monomial_lcm(g, h) == g * h

    def monomial_all_divisors(self, t):
        r"""
        Return a list of all monomials that divide ``t``, coefficients are
        ignored.

        INPUT:

        - ``t`` -- a monomial

        OUTPUT: list of monomials

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict_domain
            sage: P.<x,y,z> = MPolynomialRing_polydict_domain(QQ,3, order='degrevlex')
            sage: P.monomial_all_divisors(x^2*z^3)
            [x, x^2, z, x*z, x^2*z, z^2, x*z^2, x^2*z^2, z^3, x*z^3, x^2*z^3]

        ALGORITHM: addwithcarry idea by Toon Segers
        """

        def addwithcarry(tempvector, maxvector, pos):
            if tempvector[pos] < maxvector[pos]:
                tempvector[pos] += 1
            else:
                tempvector[pos] = 0
                tempvector = addwithcarry(tempvector, maxvector, pos + 1)
            return tempvector

        if not t.is_monomial():
            raise TypeError("only monomials are supported")

        R = self
        one = self.base_ring().one()
        M = list()

        v, = t.monomial_coefficients()
        maxvector = list(v)

        tempvector = [0] * len(maxvector)

        pos = 0

        while tempvector != maxvector:
            tempvector = addwithcarry(list(tempvector), maxvector, pos)
            M.append(R(PolyDict({ETuple(tempvector): one})))
        return M

    def sum(self, terms):
        r"""
        Return a sum of elements of this multipolynomial ring.

        This is method is much faster than the Python builtin :func:`sum`.

        EXAMPLES::

            sage: R = QQ['x']
            sage: S = R['y, z']
            sage: x = R.gen()
            sage: y, z = S.gens()
            sage: S.sum([x*y, 2*x^2*z - 2*x*y, 1 + y + z])
            (-x + 1)*y + (2*x^2 + 1)*z + 1

        Comparison with builtin :func:`sum`::

            sage: sum([x*y, 2*x^2*z - 2*x*y, 1 + y + z])
            (-x + 1)*y + (2*x^2 + 1)*z + 1
        """
        elt = PolyDict({}, check=False)
        for t in terms:
            elt += self(t).element()
        # NOTE: here we should be using self.element_class but polynomial rings are not complient
        # with categories...
        from sage.rings.polynomial.multi_polynomial_element import MPolynomial_polydict
        return MPolynomial_polydict(self, elt)


class MPolynomialRing_polydict_domain(IntegralDomain,
                                      MPolynomialRing_polydict):
    def __init__(self, base_ring, n, names, order):
        order = TermOrder(order, n)
        MPolynomialRing_polydict.__init__(self, base_ring, n, names, order)

    def is_integral_domain(self, proof=True):
        return True

    def is_field(self, proof=True):
        if self.ngens() == 0:
            return self.base_ring().is_field(proof)
        return False
