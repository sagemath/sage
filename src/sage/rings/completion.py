from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.fields import Fields

from sage.rings.morphism import RingHomomorphism
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.power_series_ring import PowerSeriesRing

from sage.rings.ring_extension import RingExtension_generic
from sage.rings.ring_extension_element import RingExtensionElement


class CompletionToPowerSeries(RingHomomorphism):
    def __init__(self, parent):
        super().__init__(parent)

    def _call_(self, x):
        return self.codomain()(x.backend(force=True))


class CompletionPolynomial(RingExtensionElement):
    def _repr_(self):
        prec = self.precision_absolute()

        # Uniformizer
        u = self.parent()._place
        if u._is_atomic():
            unif = str(u)
        else:
            unif = "(%s)" % u

        # Bigoh
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
        for i in range(prec):
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
        S = self.parent()
        A = S._base
        incl = S._incl
        u = self.parent()._place
        v = S._place.derivative()(S._xbar).inverse()
        vv = v.parent().one()
        uu = A.one()
        elt = self.backend(force=True)
        prec = self.precision_absolute()
        n = 0
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


class CompletionPolynomialRing(UniqueRepresentation, RingExtension_generic):
    Element = CompletionPolynomial

    def __classcall_private__(cls, A, p, default_prec=20, sparse=False):
        if not isinstance(A, PolynomialRing_generic):
            raise ValueError("not a polynomial ring")
        if isinstance(p, str):
            return PowerSeriesRing(A.base_ring(), default_prec=default_prec,
                                   names=p, sparse=sparse)
        if not A.base_ring() in Fields():
            raise NotImplementedError
        try:
            p = p.squarefree_part()
        except (AttributeError, NotImplementedError):
            pass
        p = p.monic()
        return cls.__classcall__(cls, A, p, default_prec, sparse)

    def __init__(self, A, p, default_prec, sparse):
        self._A = A
        base = A.base_ring()
        self._place = p
        self._residue_ring = k = A.quotient(p)
        self._xbar = xbar = k(A.gen())
        name = "u_%s" % id(self)
        backend = PowerSeriesRing(k, name, default_prec=default_prec, sparse=sparse)
        gen = backend.gen() + xbar
        self._incl = A.hom([gen])
        self._base_morphism = self._incl * A.coerce_map_from(base)
        super().__init__(self._base_morphism)
        self._gen = self(gen)
        coerce = A.Hom(self)(self._incl)
        self.register_coercion(coerce)

    def _repr_(self):
        s = "Completion of %s at %s" % (self._A, self._place)
        return s

    def construction(self):
        from sage.categories.pushout import CompletionFunctor
        extras = {
            'names': [str(x) for x in self._defining_names()],
            'sparse': self.is_sparse()
        }
        return CompletionFunctor(self._place, self.default_prec(), extras), self._A

    @cached_method
    def uniformizer(self):
        return self(self._place)

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
        S = PowerSeriesRing(self._residue_ring, names, sparse=sparse)
        S.register_conversion(CompletionToPowerSeries(self.Hom(S)))
        return S
