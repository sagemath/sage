r"""
Completion of polynomial rings and their fraction fields
"""

# ***************************************************************************
#    Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************


from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element

from sage.categories.fields import Fields

from sage.rings.morphism import RingHomomorphism
from sage.rings.infinity import Infinity

from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.fraction_field import FractionField_1poly_field
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing

from sage.rings.ring_extension import RingExtension_generic
from sage.rings.ring_extension_element import RingExtensionElement


class CompletionToPowerSeries(RingHomomorphism):
    r"""
    Conversion morphism from a completion to the
    underlying power series ring.

    TESTS::

        sage: A.<x> = QQ[]
        sage: Ap = A.completion(x - 1)
        sage: S.<u> = Ap.power_series_ring()
        sage: f = S.convert_map_from(Ap)
        sage: type(f)
        <class 'sage.rings.completion_polynomial_ring.CompletionToPowerSeries'>

        sage: # TestSuite(f).run()
    """
    def _call_(self, x):
        r"""
        Return the image of ``x`` by this morphism.

        TESTS::

            sage: A.<x> = QQ[]
            sage: Ap = A.completion(x - 1)
            sage: S.<u> = Ap.power_series_ring()
            sage: S(Ap(x))  # indirect doctest
            1 + u
        """
        return self.codomain()(x.backend(force=True))


class CompletionPolynomial(RingExtensionElement):
    r"""
    An element in the completion of a polynomial ring
    or a field of rational functions.

    TESTS::

        sage: A.<x> = QQ[]
        sage: Ap = A.completion(x - 1)
        sage: u = Ap.random_element()
        sage: type(u)
        <class 'sage.rings.completion_polynomial_ring.CompletionPolynomialRing_with_category.element_class'>
    """
    def __init__(self, parent, f):
        if isinstance(f, Element):
            R = f.parent()
            ring = parent._ring
            integer_ring = parent._integer_ring
            if ring.has_coerce_map_from(R):
                f = ring(f)
                f = f(parent._gen)
            #elif integer_ring.has_coerce_map_from(R):
            #    f = integer_ring(f).backend(force=True)
        return super().__init__(parent, f)

    def _repr_(self):
        r"""
        Return a string representation of this element.

        TESTS::

            sage: A.<x> = QQ[]
            sage: Ap = A.completion(x - 1)
            sage: y = Ap(x)  # indirect doctest
            sage: y
            1 + (x - 1)
            sage: y.add_bigoh(5)  # indirect doctest
            1 + (x - 1) + O((x - 1)^5)

        ::

            sage: Ainf = A.completion(infinity)
            sage: Ainf.uniformizer()  # indirect doctest
            x^-1
        """
        # Uniformizer
        S = self.parent()
        u = S._p
        if u is Infinity:
            step = -1
            unif = S._ring.variable_name()
        else:
            step = 1
            if u._is_atomic():
                unif = str(u)
            else:
                unif = "(%s)" % u

        # Bigoh
        prec = self.precision_absolute()
        if prec is Infinity:
            prec = self.parent()._default_prec
            bigoh = "..."
        elif prec == 0:
            bigoh = "O(1)"
        elif step*prec == 1:
            bigoh = "O(%s)" % unif
        else:
            bigoh = "O(%s^%s)" % (unif, step*prec)

        E = self.expansion(include_final_zeroes=False)
        is_exact = False
        terms = []
        if S._integer_ring is S:
            start = 0
        else:
            start = min(prec, self.valuation())
            if start is Infinity:
                return "0"
        if prec is Infinity:
            prec = start + 10  # ???
        for e in range(step*start, step*prec, step):
            try:
                coeff = next(E)
            except StopIteration:
                is_exact = True
                break
            if coeff.is_zero():
                continue
            if coeff._is_atomic():
                coeff = str(coeff)
            else:
                coeff = "(%s)" % coeff
            if e == 0:
                term = coeff
            elif e == 1:
                term = "%s*%s" % (coeff, unif)
            else:
                term = "%s*%s^%s" % (coeff, unif, e)
            terms.append(term)
        if len(terms) == 0 and bigoh == "...":
            terms = ["0"]
        if not is_exact:
            terms.append(bigoh)
        s = " " + " + ".join(terms)
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        return s[1:]

    def expansion(self, include_final_zeroes=True):
        r"""
        Return a generator producing the list of coefficients
        of this element, when written as a series in `p`.

        If the parent of this element does not contain elements
        of negative valuation, the expansion starts at `0`;
        otherwise, it starts at the valuation of the elment.

        INPUT:

        - ``include_final_zeroes`` (default: ``True``) : a boolean;
          if ``False``, stop the iterator as soon as all the next
          coefficients are all known to be zero

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: Ap = A.completion(x - 1)
            sage: y = Ap(x)  # indirect doctest
            sage: y
            1 + (x - 1)
            sage: E = y.expansion()
            sage: next(E)
            1
            sage: next(E)
            1
            sage: next(E)
            0

        Using ``include_final_zeroes=False`` stops the iterator after
        the two first `1`::

            sage: E = y.expansion(include_final_zeroes=False)
            sage: next(E)
            1
            sage: next(E)
            1
            sage: next(E)
            Traceback (most recent call last):
            ...
            StopIteration

        We underline that, for a nonexact element the iterator always
        stops at the precision::

            sage: z = y.add_bigoh(3)
            sage: z
            1 + (x - 1) + O((x - 1)^3)
            sage: E = z.expansion(include_final_zeroes=False)
            sage: next(E)
            1
            sage: next(E)
            1
            sage: next(E)
            0
            sage: next(E)
            Traceback (most recent call last):
            ...
            StopIteration

        Over the completion of a field of rational functions, the
        iterator always starts at the first nonzero coefficient
        (correspoding to the valuation of the element)::

            sage: Kp = Ap.fraction_field()
            sage: u = Kp.uniformizer()
            sage: u
            (x - 1)
            sage: E = u.expansion()
            sage: next(E)
            1
        """
        S = self.parent()
        A = S._base
        p = S._p
        prec = self.precision_absolute()
        if p is Infinity:
            elt = self.backend(force=True)
            n = elt.valuation()
            while n < prec:
                if (not include_final_zeroes
                and elt.precision_absolute() is Infinity and elt.is_zero()):
                     break
                coeff = elt[n]
                yield coeff
                elt -= elt.parent()(coeff) << n
                n += 1
        else:
            if S._integer_ring is S:
                n = 0
            else:
                n = self.valuation()
                if n < 0:
                    self *= p ** (-n)
                    prec -= n
                    n = 0
                self = S._integer_ring(self)
            try:
                f = self.polynomial()
                current_prec = prec
            except ValueError:
                current_prec = n + 20
                f = self.polynomial(current_prec)
            if current_prec is not Infinity:
                include_final_zeroes = True
            f //= p ** n
            while n < prec:
                if n >= current_prec:
                    current_prec *= 2
                    f = self.add_bigoh(current_prec).polynomial()
                    f //= p ** n
                if not include_final_zeroes and f.is_zero():
                     break
                f, coeff = f.quo_rem(p)
                yield coeff
                n += 1

    def teichmuller_lift(self):
        r"""
        Return the Teichmüller representative of this element.

        EXAMPLES::

            sage: A.<x> = GF(5)[]
            sage: Ap = A.completion(x^2 + x + 1, prec=5)
            sage: a = Ap(x).teichmuller_lift()
            sage: a
            x + (4*x + 2)*(x^2 + x + 1) + (4*x + 2)*(x^2 + x + 1)^2 + ...
            sage: a^2 + a + 1
            0
        """
        elt = self.backend(force=True)
        lift = elt.parent()(elt[0])
        return self.parent()(lift)

    def valuation(self):
        r"""
        Return the valuation of this element.

        EXAMPLES::

            sage: A.<x> = GF(5)[]
            sage: Ap = A.completion(x^2 + x + 1)
            sage: Ap(x).valuation()
            0
            sage: u = Ap(x^2 + x + 1)
            sage: u.valuation()
            1
            sage: (1/u).valuation()
            -1
        """
        return self.backend(force=True).valuation()

    def lift_to_precision(self, prec):
        r"""
        Return this element lifted to precision ``prec``.

        INPUT:

        - ``prec`` -- an integer

        EXAMPLES::

            sage: A.<x> = GF(5)[]
            sage: Ap = A.completion(x^2 + x + 1)
            sage: y = Ap(x).add_bigoh(10)
            sage: y
            x + O((x^2 + x + 1)^10)
            sage: y.lift_to_precision(20)
            x + O((x^2 + x + 1)^20)

        When ``prec`` is less than the absolute precision of
        the element, the same element is returned without any
        change on the precision::

            sage: y.lift_to_precision(5)
            x + O((x^2 + x + 1)^10)
        """
        elt = self.backend(force=True)
        elt = elt.lift_to_precision(prec)
        return self.parent()(elt)

    def polynomial(self, prec=None):
        r"""
        Return a polynomial (in the underlying polynomial ring)
        which is indistinguishable from ``self`` at precision ``prec``.

        INPUT:

        - ``prec`` (optional) -- an integer; if not given, defaults
          to the absolute precision of this element

        EXAMPLES::

            sage: A.<x> = GF(5)[]
            sage: Ap = A.completion(x - 2, prec=5)
            sage: y = 1 / Ap(1 - x)
            sage: y
            4 + (x + 3) + 4*(x + 3)^2 + (x + 3)^3 + 4*(x + 3)^4 + O((x + 3)^5)
            sage: y.polynomial()
            4*x^4 + 4*x^3 + 4*x^2 + 4*x + 4

        When called with a nonintegral element, this method raises an error::

            sage: z = 1 / Ap(2 - x)
            sage: z.polynomial()
            Traceback (most recent call last):
            ...
            ValueError: this element of negative valuation cannot be approximated by a polynomial

        TESTS::

            sage: Kq = A.completion(infinity)
            sage: Kq(1/x).polynomial()
            Traceback (most recent call last):
            ...
            ValueError: approximation by a polynomial does not make sense for a completion at infinity
        """
        S = self.parent()
        A = S._A
        p = S._p
        if p is Infinity:
            raise ValueError("approximation by a polynomial does not make sense for a completion at infinity")
        backend = self.backend(force=True)
        if backend.valuation() < 0:
            raise ValueError("this element of negative valuation cannot be approximated by a polynomial")
        try:
            backend = backend.power_series()
        except AttributeError:
            pass
        d = p.degree()
        k = S.residue_field()
        if prec is not None:
            backend = backend.add_bigoh(prec)
        prec = backend.precision_absolute()
        f = backend.polynomial().change_ring(k)
        g = f(f.parent().gen() - k.gen())
        gs = [A([c[i] for c in g.list()]) for i in range(d)]
        if prec is Infinity:
            for i in range(1, d):
                if gs[i]:
                    raise ValueError("this element is exact and does not lie in the polynomial ring")
            return gs[0]
        xbar = S._xbar(prec)
        modulus = p ** prec
        res = gs[0]
        xbari = A.one()
        for i in range(1, d):
            xbari = (xbari * xbar) % modulus
            res += xbari * gs[i]
        return res % modulus

    def shift(self, n):
        raise NotImplementedError


class CompletionPolynomialRing(UniqueRepresentation, RingExtension_generic):
    r"""
    A class for completions of polynomial rings and their fraction fields.
    """
    Element = CompletionPolynomial

    def __classcall_private__(cls, ring, p, default_prec=20, sparse=False):
        r"""
        Normalize the parameters and call the appropriate constructor.

        INPUT:

        - ``ring`` -- the underlying polynomial ring or field

        - ``p`` -- a generator of the ideal at which the completion is
          done; it could also be ``Infinity``

        - ``default_pres`` (default: ``20``) -- the default precision
          of the completion

        - ``sparse`` (default: ``False``) -- a boolean

        TESTS::

            sage: A.<x> = ZZ[]
            sage: A.completion(x^2 + x + 1)
            Completion of Univariate Polynomial Ring in x over Integer Ring at x^2 + x + 1
            sage: A.completion(2*x - 1)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of p must be invertible

        ::

            sage: B.<y> = QQ[]
            sage: B.completion(x^2 + x^3)
            Traceback (most recent call last):
            ...
            ValueError: p must be a squarefree polynomial

        When passing a variable name, a standard ring of power series is
        constructed::

            sage: A.completion('x')
            Power Series Ring in x over Integer Ring
            sage: A.completion(x)
            Completion of Univariate Polynomial Ring in x over Integer Ring at x
        """
        if isinstance(ring, PolynomialRing_generic):
            if isinstance(p, str):
                return PowerSeriesRing(ring.base_ring(),
                                       default_prec=default_prec,
                                       names=p, sparse=sparse)
            A = ring
        elif isinstance(ring, FractionField_1poly_field):
            if isinstance(p, str):
                return LaurentSeriesRing(ring.base_ring(),
                                         default_prec=default_prec,
                                         names=p, sparse=sparse)
            A = ring.ring()
        else:
            raise NotImplementedError("not a polynomial ring or a rational function field")
        if p is Infinity:
            ring = ring.fraction_field()
        else:
            p = A(p)
            if not p.leading_coefficient().is_unit():
                raise NotImplementedError("the leading coefficient of p must be invertible")
            p = A(p.monic())
            try:
                if not p.is_squarefree():
                    raise ValueError("p must be a squarefree polynomial")
            except (AttributeError, NotImplementedError):
                pass
        return cls.__classcall__(cls, ring, p, default_prec, sparse)

    def __init__(self, ring, p, default_prec, sparse):
        r"""
        Initialize this ring.

        INPUT:

        - ``ring`` -- the underlying polynomial ring or field

        - ``p`` -- a generator of the ideal at which the completion is
          done; it could also be ``Infinity``

        - ``default_pres`` -- the default precision of the completion

        - ``sparse`` -- a boolean

        TESTS::

            sage: A.<x> = QQ[]
            sage: Ap = A.completion(x - 1)
            sage: type(Ap)
            <class 'sage.rings.completion_polynomial_ring.CompletionPolynomialRing_with_category'>
            sage: TestSuite(Ap).run()

        ::

            sage: K = Frac(A)
            sage: Kp = K.completion(x - 1)
            sage: type(Kp)
            <class 'sage.rings.completion_polynomial_ring.CompletionPolynomialRing_with_category'>
            sage: TestSuite(Kp).run()
        """
        A = self._ring = ring
        if not isinstance(A, PolynomialRing_generic):
            A = A.ring()
        self._A = A
        base = A.base_ring()
        self._p = p
        self._default_prec = default_prec
        # Construct backend
        if p is Infinity:
            self._residue_ring = k = base
            self._integer_ring = None
            name = "u_%s" % id(self)
            backend = LaurentSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
            x = backend.gen().inverse()
        else:
            self._residue_ring = k = A.quotient(p)
            self._xbar_approx = A.gen()
            self._xbar_modulus = p
            self._xbar_prec = 1
            if isinstance(ring, PolynomialRing_generic):
                self._integer_ring = self
                name = "u_%s" % id(self)
                backend = PowerSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
            else:
                self._integer_ring = CompletionPolynomialRing(A, p, default_prec, sparse)
                name = "u_%s" % id(self._integer_ring)
                backend = LaurentSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
            x = backend.gen() + k(A.gen())
        super().__init__(backend.coerce_map_from(k) * k.coerce_map_from(base), category=backend.category())
        # Set generator
        self._gen = self(x)
        # Set coercions
        self.register_coercion(ring)
        if A is not ring:
            self.register_coercion(A)
        if self._integer_ring is not None:
            self.register_coercion(self._integer_ring)
        if p is not Infinity and p.is_gen():
            S = PowerSeriesRing(base, A.variable_name(), default_prec=default_prec, sparse=sparse)
            self.register_coercion(S)
            if self._integer_ring is not self:
                S = LaurentSeriesRing(base, A.variable_name(), default_prec=default_prec, sparse=sparse)
                self.register_coercion(S)

    def _repr_(self):
        r"""
        Return a string representation of this ring.

        EXAMPLES::

            sage: A.<x> = GF(7)[]
            sage: A.completion(x)  # indirect doctest
            Completion of Univariate Polynomial Ring in x over Finite Field of size 7 at x
            sage: A.completion(infinity)  # indirect doctest
            Completion of Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 7 at infinity

        ::

            sage: K = Frac(A)
            sage: K.completion(x^2 + 2*x + 2)  # indirect doctest
            Completion of Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 7 at x^2 + 2*x + 2
        """
        if self._p is Infinity:
            return "Completion of %s at infinity" % self._ring
        else:
            return "Completion of %s at %s" % (self._ring, self._p)

    def construction(self):
        r"""
        Return the functorial construction of this ring.

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: Ap = A.completion(x^2 - 1)
            sage: F, X = Ap.construction()
            sage: F
            Completion[x^2 - 1, prec=20]
            sage: X
            Univariate Polynomial Ring in x over Rational Field
            sage: F(X) is Ap
            True

        TESTS::

            sage: R = PolynomialRing(QQ, 'y', sparse=True)
            sage: Kinf = R.completion(infinity)
            sage: F, X = Kinf.construction()
            sage: F(X) is Kinf
            True
        """
        from sage.categories.pushout import CompletionFunctor
        extras = {
            'names': [str(x) for x in self._defining_names()],
            'sparse': self.is_sparse()
        }
        return CompletionFunctor(self._p, self._default_prec, extras), self._ring

    @cached_method
    def uniformizer(self):
        r"""
        Return a uniformizer of this ring.

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: Ap = A.completion(x^2 - 1)
            sage: Ap.uniformizer()
            (x^2 - 1)
        """
        if self._p is Infinity:
            return self._gen.inverse()
        else:
            return self(self._p)

    def gen(self):
        r"""
        Return the generator of this ring.

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: Ap = A.completion(x^2 - 1)
            sage: Ap.gen()
            x
        """
        return self._gen

    def _xbar(self, prec):
        r"""
        Return an approximation at precision ``prec`` of the
        Teichmüller representative of `x`, the generator of
        this ring.

        This method is only for internal use.

        EXAMPLES::

            sage: A.<x> = GF(7)[]
            sage: Ap = A.completion(x^2 + 2*x + 2)
            sage: Ap._xbar(5)
            4*x^15 + 4*x^14 + x^8 + x + 6

        The method returns nothing in case of a completion
        at infinity::

            sage: Ap = A.completion(infinity)
            sage: Ap._xbar(5)
        """
        if self._p is Infinity:
            return
        # We solve the equation p(xbar) = 0 in A_p
        p = self._p
        pp = p.derivative()
        while self._xbar_prec < prec:
            self._xbar_modulus *= self._xbar_modulus
            self._xbar_prec *= 2
            u = p(self._xbar_approx)
            _, v, _ = pp(self._xbar_approx).xgcd(self._xbar_modulus)
            self._xbar_approx -= u * v
            self._xbar_approx %= self._xbar_modulus
        return self._xbar_approx

    def residue_ring(self):
        r"""
        Return the residue ring of this completion.

        EXAMPLES::

            sage: A.<x> = GF(7)[]
            sage: Ap = A.completion(x^2 + 2*x + 2)
            sage: Ap.residue_ring()
            Univariate Quotient Polynomial Ring in xbar over Finite Field of size 7 with modulus x^2 + 2*x + 2
        """
        return self._residue_ring

    residue_field = residue_ring

    def power_series_ring(self, names=None, sparse=None):
        r"""
        Return a power series ring, which is isomorphic to
        the integer ring of this completion.

        INPUT:

        - ``names`` -- a string, the variable name

        - ``sparse`` (optional) -- a boolean; if not given,
          use the sparcity of this ring

        EXAMPLES::

            sage: A.<x> = GF(7)[]
            sage: Ap = A.completion(x^2 + 2*x + 2)
            sage: S.<u> = Ap.power_series_ring()
            sage: S
            Power Series Ring in u over Univariate Quotient Polynomial Ring in xbar over Finite Field of size 7 with modulus x^2 + 2*x + 2

        A conversion from ``Ap`` to ``S`` is set::

            sage: xp = Ap.gen()
            sage: S(xp)
            xbar + u

        ::

            sage: Kp = Frac(Ap)
            sage: T.<u> = Kp.power_series_ring()
            sage: T
            Power Series Ring in u over Univariate Quotient Polynomial Ring in xbar over Finite Field of size 7 with modulus x^2 + 2*x + 2
            sage: T is S
            True

        .. SEEALSO::

            :meth:`laurent_series_ring`
        """
        if isinstance(names, (list, tuple)):
            names = names[0]
        if sparse is None:
            sparse = self.is_sparse()
        S = PowerSeriesRing(self._residue_ring, names, default_prec=self._default_prec, sparse=sparse)
        S.register_conversion(CompletionToPowerSeries(self.Hom(S)))
        return S

    def laurent_series_ring(self, names=None, sparse=None):
        r"""
        Return a Laurent series ring, which is isomorphic to
        the fraction field of this completion.

        INPUT:

        - ``names`` -- a string, the variable name

        - ``sparse`` (optional) -- a boolean; if not given,
          use the sparcity of this ring

        EXAMPLES::

            sage: A.<x> = GF(7)[]
            sage: Ap = A.completion(x^2 + 2*x + 2)
            sage: S.<u> = Ap.laurent_series_ring()
            sage: S
            Laurent Series Ring in u over Univariate Quotient Polynomial Ring in xbar over Finite Field of size 7 with modulus x^2 + 2*x + 2

        .. SEEALSO::

            :meth:`power_series_ring`
        """
        if isinstance(names, (list, tuple)):
            names = names[0]
        if sparse is None:
            sparse = self.is_sparse()
        S = LaurentSeriesRing(self._residue_ring, names, default_prec=self._default_prec, sparse=sparse)
        S.register_conversion(CompletionToPowerSeries(self.Hom(S)))
        return S

    def integer_ring(self):
        r"""
        Return the integer ring of this completion.

        EXAMPLES::

            sage: A.<x> = GF(7)[]
            sage: K = Frac(A)
            sage: Kp = K.completion(x^2 + 2*x + 2)
            sage: Kp
            Completion of Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 7 at x^2 + 2*x + 2
            sage: Kp.integer_ring()
            Completion of Univariate Polynomial Ring in x over Finite Field of size 7 at x^2 + 2*x + 2
        """
        if self._p is Infinity:
            raise NotImplementedError
        return self._integer_ring

    def is_integral_domain(self):
        return self._residue_ring.is_integral_domain()

    def fraction_field(self, permissive=False):
        r"""
        Return the fraction field of this completion.

        - ``permissive`` (default: ``False``) -- a boolean;
          if ``True`` and this ring is a domain, return
          instead the completion at the same place of the
          fraction field of the underlying polynomial ring

        EXAMPLES::

            sage: A.<x> = GF(7)[]
            sage: Ap = A.completion(x^2 + 2*x + 2)
            sage: Ap
            Completion of Univariate Polynomial Ring in x over Finite Field of size 7 at x^2 + 2*x + 2
            sage: Ap.fraction_field()
            Completion of Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 7 at x^2 + 2*x + 2

        Trying to apply this method with a completion at a nonprime ideal
        produces an error::

            sage: Aq = A.completion(x^2 + x + 1)
            sage: Aq.fraction_field()
            Traceback (most recent call last):
            ...
            ValueError: this ring is not an integral domain

        Nonetheless, if the flag ``permissive`` is set to ``True``,
        the localization at all nonzero divisors is returned::

            sage: Aq.fraction_field(permissive=True)
            Completion of Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 7 at x^2 + x + 1
        """
        if not (permissive or self.is_integral_domain()):
            raise ValueError("this ring is not an integral domain")
        field = self._ring.fraction_field()
        if field is self._ring:
            return self
        else:
            return CompletionPolynomialRing(field, self._p,
                       default_prec = self._default_prec,
                       sparse = self.is_sparse())
