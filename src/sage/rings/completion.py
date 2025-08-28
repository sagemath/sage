from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation

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
    def __init__(self, parent):
        super().__init__(parent)

    def _call_(self, x):
        return self.codomain()(x.backend(force=True))


class CompletionPolynomial(RingExtensionElement):
    def _repr_(self):
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
        # TODO: improve performance
        S = self.parent()
        A = S._base
        incl = S._incl
        u = self.parent()._p
        elt = self.backend(force=True)
        prec = self.precision_absolute()
        n = 0
        if u is Infinity:
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
            if S._integer_ring is not S:
                val = elt.valuation()
                if val is Infinity:
                    return
                if val < 0:
                    elt *= incl(u) ** (-val)
                else:
                    n = val
            v = S._p.derivative()(S._xbar).inverse()
            vv = v**n
            uu = u**n
            while n < prec:
                if (not include_final_zeroes
                and elt.precision_absolute() is Infinity and elt.is_zero()):
                     break
                coeff = (vv * elt[n]).lift()
                yield coeff
                elt -= incl(uu * coeff)
                uu *= u
                vv *= v
                n += 1

    def teichmuller(self):
        elt = self.backend(force=True)
        lift = elt.parent()(elt[0])
        return self.parent()(lift)

    def valuation(self):
        return self.backend(force=True).valuation()

    def shift(self, n):
        raise NotImplementedError


class CompletionPolynomialRing(UniqueRepresentation, RingExtension_generic):
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
            self._xbar = xbar = k(A.gen())
            if isinstance(ring, PolynomialRing_generic):
                self._integer_ring = self
                name = "u_%s" % id(self)
                backend = PowerSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
            else:
                self._integer_ring = CompletionPolynomialRing(A, p, default_prec, sparse)
                name = "u_%s" % id(self._integer_ring)
                backend = LaurentSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
            x = backend.gen() + xbar
        super().__init__(backend.coerce_map_from(k) * k.coerce_map_from(base), category=backend.category())
        # Set generator
        self._gen = self(x)
        # Set coercions
        self._incl = ring.hom([x])
        self.register_coercion(A.Hom(self)(self._incl))
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

    def residue_ring(self):
        return self._residue_ring

    residue_field = residue_ring

    def power_series(self, names=None, sparse=None):
        if isinstance(names, (list, tuple)):
            names = names[0]
        if sparse is None:
            sparse = self.is_sparse()
        base = self._base.base_ring()
        if self._integer_ring is self:
            S = PowerSeriesRing(self._residue_ring, names, default_prec=self._default_prec, sparse=sparse)
        else:
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
