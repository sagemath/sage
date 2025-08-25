from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent

from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.power_series_ring import PowerSeriesRing, PowerSeriesRing_generic
from sage.rings.power_series_ring_element import PowerSeries


class CompletionPolynomialRing(PowerSeriesRing_generic):
    def __classcall_private__(cls, A, p, default_prec=20, names=None, sparse=False):
        if not isinstance(A, PolynomialRing_generic):
            raise ValueError("not a polynomial ring")
        if p == A.gen():
            return PowerSeriesRing(A.base_ring(), default_prec=default_prec,
                                   names=A.variable_name(), sparse=sparse)
        lc = p.leading_coefficient()
        if not lc.is_unit():
            raise NotImplementedError("leading coefficient must be a unit")
        p = lc.inverse_of_unit() * p
        if not isinstance(names, (list, tuple)):
            names = (names,)
        if p.degree() == 0:
            raise ValueError("place must be an irreducible polynomial")
        if p.degree() == 1:
            if len(names) != 1:
                raise ValueError("you should provide a variable name")
        else:
            if len(names) != 2:
                raise ValueError("you should provide two variable names")
        return cls.__classcall__(cls, A, p, default_prec, names)

    def __init__(self, A, p, default_prec, names):
        name_uniformizer = names[0]
        name_residue = None
        if len(names) > 1:
            name_residue = names[1]
        base = A.base_ring()
        self._polynomial_ring = A
        self._place = p
        if p.degree() > 1:
            k = base.extension(p, names=name_residue)
            a = k.gen()
        else:
            k = base
            a = -p[0]
        super().__init__(k, name_uniformizer, default_prec)
        self._a = a
        self._inject = A.hom([self.gen() + a])
        self.register_coercion(self._inject)

    def _repr_(self):
        A = self._polynomial_ring
        s = "Completion of %s at %s:\n" % (A, self._place)
        s += "Power Series Ring in %s = %s - %s over %s" % (self.gen(), A.gen(), self._a, self.residue_field())
        return s

    def construction(self):
        from sage.categories.pushout import CompletionFunctor
        extras = {
            'names': [str(x) for x in self._defining_names()],
            'sparse': self.is_sparse()
        }
        return CompletionFunctor(self._place, self.default_prec(), extras), self._poly_ring()

    def _defining_names(self):
        if self._place.degree() > 1:
            return (self.gen(), self(self.base_ring().gen()))
        else:
            return (self.gen())
