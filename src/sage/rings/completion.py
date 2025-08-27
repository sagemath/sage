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
        if u._is_atomic():
            unif = str(u)
        else:
            unif = "(%s)" % u

        # Bigoh
        prec = self.precision_absolute()
        if prec is Infinity:
            prec = self.parent().default_prec()
            bigoh = "..."
        elif prec == 0:
            bigoh = "O(1)"
        elif prec == 1:
            bigoh = "O(%s)" % unif
        else:
            bigoh = "O(%s^%s)" % (unif, prec)

        E = self.expansion(include_final_zeroes=False)
        is_exact = False
        terms = []
        if S._integer_ring is S:
            start = 0
        else:
            start = self.valuation()
            if start is Infinity:
                return "0"
        for i in range(start, prec):
            try:
                coeff = next(E)
            except StopIteration:
                is_exact = True
                break
            if coeff.is_zero():
                continue
            if coeff.is_one():
                if i == 0:
                    term = "1"
                elif i == 1:
                    term = unif
                else:
                    term = "%s^%s" % (unif, i)
            else:
                if coeff._is_atomic():
                    coeff = str(coeff)
                else:
                    coeff = "(%s)" % coeff
                if i == 0:
                    term = coeff
                elif i == 1:
                    term = "%s*%s" % (coeff, unif)
                else:
                    term = "%s*%s^%s" % (coeff, unif, i)
            terms.append(term)
        if len(terms) == 0 and bigoh == "...":
            terms = ["0"]
        if not is_exact:
            terms.append(bigoh)
        return " + ".join(terms)

    def expansion(self, include_final_zeroes=True):
        # TODO: improve performance
        S = self.parent()
        A = S._base
        incl = S._incl
        u = self.parent()._p
        v = S._p.derivative()(S._xbar).inverse()
        elt = self.backend(force=True)
        prec = self.precision_absolute()
        n = 0
        if S._integer_ring is not S:
            val = elt.valuation()
            if val is Infinity:
                return
            if val < 0:
                elt *= incl(u) ** (-val)
            else:
                n = val
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
        if not A.base_ring() in Fields():
            raise NotImplementedError
        p = A(p)
        try:
            p = p.squarefree_part()
        except (AttributeError, NotImplementedError):
            pass
        p = p.monic()
        return cls.__classcall__(cls, ring, p, default_prec, sparse)

    def __init__(self, ring, p, default_prec, sparse):
        A = self._ring = ring
        if not isinstance(A, PolynomialRing_generic):
            A = A.ring()
        base = A.base_ring()
        self._p = p
        # Construct backend
        self._residue_ring = k = A.quotient(p)
        self._xbar = xbar = k(A.gen())
        name = "u_%s" % hash(self._p)
        if isinstance(ring, PolynomialRing_generic):
            self._integer_ring = self
            name = "u_%s" % id(self)
            backend = PowerSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
        else:
            self._integer_ring = CompletionPolynomialRing(A, p, default_prec, sparse)
            name = "u_%s" % id(self._integer_ring)
            backend = LaurentSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
        super().__init__(backend.coerce_map_from(k) * k.coerce_map_from(base))
        # Set generator
        x = backend.gen() + xbar
        self._gen = self(x)
        # Set coercions
        self._incl = ring.hom([x])
        self.register_coercion(A.Hom(self)(self._incl))
        self.register_coercion(self._integer_ring)

    def _repr_(self):
        s = "Completion of %s at %s" % (self._ring, self._p)
        return s

    def construction(self):
        from sage.categories.pushout import CompletionFunctor
        extras = {
            'names': [str(x) for x in self._defining_names()],
            'sparse': self.is_sparse()
        }
        return CompletionFunctor(self._p, self.default_prec(), extras), self._ring

    @cached_method
    def uniformizer(self):
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
            S = PowerSeriesRing(self._residue_ring, names, sparse=sparse)
        else:
            S = LaurentSeriesRing(self._residue_ring, names, sparse=sparse)
        S.register_conversion(CompletionToPowerSeries(self.Hom(S)))
        return S

    def integer_ring(self):
        return self._integer_ring

    def fraction_field(self):
        field = self._ring.fraction_field()
        if field is self._ring:
            return self
        else:
            return CompletionPolynomialRing(field, self._p,
                       default_prec = self.default_prec(),
                       sparse = self.is_sparse())
