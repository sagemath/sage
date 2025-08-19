r"""
Morphisms attached to polynomial rings.
"""

# *****************************************************************************
#        Copyright (C) 2025 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.category_object import normalize_names
from sage.categories.rings import Rings

from sage.rings.morphism import RingHomomorphism
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing


class MorphismToCompletion(RingHomomorphism):
    r"""
    Morphisms from a polynomial ring (or its fraction field)
    to the completion at one place.

    TESTS::

        sage: A.<t> = GF(5)[]
        sage: B.<u> = A.completion(t+1)
        sage: f = B.coerce_map_from(A)
        sage: type(f)
        <class 'sage.rings.polynomial.morphism.MorphismToCompletion'>

        sage: TestSuite(f).run(skip='_test_category')
    """
    def __init__(self, domain, place, prec, name, residue_name):
        r"""
        Initialize this morphism.

        INPUT:

        - ``domain`` -- a polynomial ring of its fraction field

        - ``place`` -- an irreducible polynomial or ``Infinity``

        - ``prec`` -- an integer or ``Infinity``

        - ``name`` -- a string, the variable name of the uniformizer

        - ``residue_name`` -- a string, the variable name of the
          generator of the residue ring

        TESTS::

            sage: A.<t> = ZZ[]
            sage: B.<u> = A.completion(2*t + 1)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the place is not a unit

        ::

            sage: A.<t> = QQ[]
            sage: B.<u,a> = A.completion(x^2 + 2*x + 1)
            Traceback (most recent call last):
            ...
            ValueError: place must be Infinity or an irreducible polynomial
        """
        if domain.is_field():
            ring = domain.ring()
            SeriesRing = LaurentSeriesRing
        else:
            ring = domain
            SeriesRing = PowerSeriesRing
        k = base = ring.base_ring()
        x = ring.gen()
        if place is Infinity:
            pass
        elif place in ring:
            place = ring(place)
            if place.leading_coefficient().is_unit():
                place = ring(place.monic())
                if not place.is_irreducible():
                    raise ValueError("place must be Infinity or an irreducible polynomial")
            else:
                raise NotImplementedError("the leading coefficient of the place is not a unit")
        else:
            raise ValueError("place must be Infinity or an irreducible polynomial")
        self._place = place
        if name is None:
            raise ValueError("you must specify a variable name")
        name = normalize_names(1, name)
        if place is Infinity:
            codomain = LaurentSeriesRing(base, names=name, default_prec=prec)
            image = codomain.one() >> 1
        elif place.degree() == 1:
            codomain = SeriesRing(base, names=name, default_prec=prec)
            image = codomain.gen() - place[0]
        else:
            if residue_name is None:
                raise ValueError("you must specify a variable name for the residue field")
            residue_name = normalize_names(1, residue_name)
            k = base.extension(place, names=residue_name)
            codomain = SeriesRing(k, names=name, default_prec=prec)
            image = codomain.gen() + k.gen()
        parent = domain.Hom(codomain, category=Rings())
        RingHomomorphism.__init__(self, parent)
        self._image = image
        self._k = k
        self._q = k.cardinality()

    def _repr_type(self):
        r"""
        Return a string that describes the type of this morphism.

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: B.<u> = A.completion(x + 1)
            sage: f = B.coerce_map_from(A)
            sage: f  # indirect doctest
            Completion morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Power Series Ring in u over Rational Field
        """
        return "Completion"

    def place(self):
        r"""
        Return the place at which we completed.

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: B.<u> = A.completion(x + 1)
            sage: f = B.coerce_map_from(A)
            sage: f.place()
            x + 1
        """
        return self._place

    def _call_(self, P):
        r"""
        Return the image of ``P`` under this morphism.

        EXAMPLES::

            sage: A.<x> = QQ[]
            sage: B.<u> = A.completion(x + 1)
            sage: f = B.coerce_map_from(A)
            sage: f(x)  # indirect doctest
            -1 + u
        """
        return P(self._image)
