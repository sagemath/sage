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
        <class 'sage.rings.completion.CompletionToPowerSeries'>

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
        <class 'sage.rings.completion.CompletionPolynomialRing_with_category.element_class'>
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
                f = self.add_bigoh(current_prec).polynomial()
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
        elt = self.backend(force=True)
        elt = elt.lift_to_precision(prec)
        return self.parent()(elt)

    def polynomial(self):
        S = self.parent()
        A = S._ring
        p = S._p
        d = p.degree()
        k = S.residue_field()
        f = self.backend(force=True).polynomial().change_ring(k)
        g = f(f.parent().gen() - k.gen())
        gs = [A([c[i] for c in g.list()]) for i in range(d)]
        prec = self.precision_absolute()
        if prec is Infinity:
            for i in range(1, d):
                if gs[i]:
                    raise ValueError("exact element which is not in the polynomial ring")
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

    TESTS::

        sage: A.<x> = QQ[]
        sage: Ap = A.completion(x - 1)
        sage: type(Ap)
        <class 'sage.rings.completion.CompletionPolynomialRing_with_category'>

        sage: TestSuite(Ap).run()

    """
    Element = CompletionPolynomial

    def __classcall_private__(cls, ring, p, default_prec=20, sparse=False):
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
            raise ValueError("not a polynomial ring or a rational function field")
        if p is not Infinity:
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
        A = self._ring = ring
        if not isinstance(A, PolynomialRing_generic):
            A = A.ring()
        self._A = A
        base = A.base_ring()
        self._p = p
        self._default_prec = default_prec
        # Construct backend
        if p is Infinity:
            self._ring = ring.fraction_field()
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

    def _repr_(self):
        if self._p is Infinity:
            return "Completion of %s at infinity" % self._ring
        else:
            return "Completion of %s at %s" % (self._ring, self._p)

    def construction(self):
        from sage.categories.pushout import CompletionFunctor
        extras = {
            'names': [str(x) for x in self._defining_names()],
            'sparse': self.is_sparse()
        }
        return CompletionFunctor(self._p, self._default_prec, extras), self._ring

    @cached_method
    def uniformizer(self):
        if self._p is Infinity:
            return self._gen.inverse()
        else:
            return self(self._p)

    def gen(self):
        return self._gen

    def _xbar(self, prec):
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
        return self._residue_ring

    residue_field = residue_ring

    def power_series_ring(self, names=None, sparse=None):
        if isinstance(names, (list, tuple)):
            names = names[0]
        if sparse is None:
            sparse = self.is_sparse()
        S = PowerSeriesRing(self._residue_ring, names, default_prec=self._default_prec, sparse=sparse)
        S.register_conversion(CompletionToPowerSeries(self.Hom(S)))
        return S

    def laurent_series_ring(self, names=None, sparse=None):
        if isinstance(names, (list, tuple)):
            names = names[0]
        if sparse is None:
            sparse = self.is_sparse()
        S = LaurentSeriesRing(self._residue_ring, names, default_prec=self._default_prec, sparse=sparse)
        S.register_conversion(CompletionToPowerSeries(self.Hom(S)))
        return S

    def integer_ring(self):
        if self._p is Infinity:
            raise NotImplementedError
        return self._integer_ring

    def fraction_field(self):
        field = self._ring.fraction_field()
        if field is self._ring:
            return self
        else:
            return CompletionPolynomialRing(field, self._p,
                       default_prec = self._default_prec,
                       sparse = self.is_sparse())
