from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.fields import Fields

from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.power_series_ring import PowerSeriesRing

from sage.rings.ring_extension import RingExtension_generic
from sage.rings.ring_extension_element import RingExtensionElement


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

        S = self.parent()
        A = S._base
        incl = S._incl
        v = S._place.derivative()(S._a).inverse()
        w = v.parent().one()
        uu = A.one()
        elt = self.backend(force=True)
        terms = []
        is_exact = False
        for i in range(prec):
            if elt.precision_absolute() is Infinity and elt.is_zero():
                is_exact = True
                break
            coeff = A((w * elt[i]).polynomial())
            elt -= incl(uu * coeff)
            uu *= u
            w *= v
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

    def teichmuller(self):
        elt = self.backend(force=True)
        lift = elt.parent()(elt[0])
        return self.parent()(lift)


class CompletionPolynomialRing(UniqueRepresentation, RingExtension_generic):
    Element = CompletionPolynomial

    def __classcall_private__(cls, A, p, default_prec=20, sparse=False):
        if not isinstance(A, PolynomialRing_generic):
            raise ValueError("not a polynomial ring")
        if p == A.gen():
            return PowerSeriesRing(A.base_ring(), default_prec=default_prec,
                                   names=A.variable_name(), sparse=sparse)
        if not A.base_ring() in Fields():
            raise NotImplementedError
        if not p.is_irreducible():
            raise ValueError("the place must be an irreducible polynomial")
        p = p.monic()
        return cls.__classcall__(cls, A, p, default_prec, sparse)

    def __init__(self, A, p, default_prec, sparse):
        base = A.base_ring()
        self._place = p
        if p.degree() > 1:
            self._extension = True
            k = base.extension(p, names='a')
            a = k.gen()
        else:
            self._extension = False
            k = base
            a = -p[0]
        backend = PowerSeriesRing(k, 'u', sparse=sparse)
        self._gen = backend.gen() + a
        self._incl = A.hom([self._gen])
        super().__init__(self._incl)
        self._a = a

    def _repr_(self):
        A = self._base
        s = "Completion of %s at %s" % (A, self._place)
        return s

    def construction(self):
        from sage.categories.pushout import CompletionFunctor
        extras = {
            'names': [str(x) for x in self._defining_names()],
            'sparse': self.is_sparse()
        }
        return CompletionFunctor(self._place, self.default_prec(), extras), self._poly_ring()

    @cached_method
    def uniformizer(self):
        return self(self._place)

    def gen(self):
        return self._gen

    def residue_field(self, names=None):
        base = self._base.base_ring()
        if self._extension:
            if names is None:
                raise ValueError("you must give a variable name")
            return base.extension(self._place, names=names)
        else:
            return base

    def power_series(self, sparse=False, names=None):
        if not isinstance(names, (list, tuple)):
            names = [names]
        base = self._base.base_ring()
        if self._extension:
            if len(names) < 2:
                raise ValueError("you must give variable names for the uniformizer and the generator of the residue field")
            k = self.residue_field(names[1])
        else:
            if len(names) < 1:
                raise ValueError("you must give variable names for the uniformizer")
            k = self.residue_field()
        return PowerSeriesRing(k, names[0], sparse=sparse)
