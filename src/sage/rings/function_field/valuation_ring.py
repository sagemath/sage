# sage.doctest: needs sage.rings.finite_rings
# sage.doctest: needs sage.rings.function_field
r"""
Valuation rings of function fields

A valuation ring of a function field is associated with a place of the
function field. The valuation ring consists of all elements of the function
field that have nonnegative valuation at the place.

EXAMPLES::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
    sage: p = L.places_finite()[0]
    sage: p
    Place (x, x*y)
    sage: R = p.valuation_ring()
    sage: R
    Valuation ring at Place (x, x*y)
    sage: R.place() == p
    True

Thus any nonzero element or its inverse of the function field lies in the
valuation ring, as shown in the following example::

    sage: f = y/(1+y)
    sage: f in R
    True
    sage: f not in R
    False
    sage: f.valuation(p)
    0

The residue field at the place is defined as the quotient ring of the valuation
ring modulo its unique maximal ideal. The method :meth:`residue_field()` of the
valuation ring returns an extension field of the constant base field, isomorphic
to the residue field, along with lifting and evaluation homomorphisms::

    sage: k,phi,psi = R.residue_field()
    sage: k
    Finite Field of size 2
    sage: phi
    Ring morphism:
      From: Finite Field of size 2
      To:   Valuation ring at Place (x, x*y)
    sage: psi
    Ring morphism:
      From: Valuation ring at Place (x, x*y)
      To:   Finite Field of size 2
    sage: psi(f) in k
    True

AUTHORS:

- Kwankyu Lee (2017-04-30): initial version
"""

# ****************************************************************************
#       Copyright (C) 2016-2019 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.categories.homset import Hom
from sage.categories.rings import Rings


class FunctionFieldValuationRing(UniqueRepresentation, Parent):
    """
    Base class for valuation rings of function fields.

    INPUT:

    - ``field`` -- function field

    - ``place`` -- place of the function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
        sage: p = L.places_finite()[0]
        sage: p.valuation_ring()
        Valuation ring at Place (x, x*y)
    """
    def __init__(self, field, place, category=None):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R = p.valuation_ring()
            sage: TestSuite(R).run()
        """
        Parent.__init__(self, category=Rings().or_subcategory(category).Infinite(), facade=field)

        self._field = field
        self._place = place

    def _element_constructor_(self, x):
        """
        Construct an element of the function field belonging to the
        valuation ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R = p.valuation_ring()
            sage: y in R
            False
            sage: 1/y in R
            True
            sage: x + y in R
            False
            sage: 1/(x + y) in R
            True
        """
        x = self._field(x)
        if x.valuation(self._place) >= 0:
            return x
        else:
            raise TypeError

    def _repr_(self):
        """
        Return the string representation of the valuation ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: p.valuation_ring()
            Valuation ring at Place (x, x*y)
        """
        return 'Valuation ring at {}'.format(self._place)

    def place(self):
        """
        Return the place associated with the valuation ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R = p.valuation_ring()
            sage: p == R.place()
            True
        """
        return self._place

    @cached_method
    def residue_field(self, name=None):
        """
        Return the residue field of the valuation ring together with
        the maps from and to it.

        INPUT:

        - ``name`` -- string; name of the generator of the field

        OUTPUT:

        - a field isomorphic to the residue field

        - a ring homomorphism from the valuation ring to the field

        - a ring homomorphism from the field to the valuation ring

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R = p.valuation_ring()
            sage: k, fr_k, to_k = R.residue_field()
            sage: k
            Finite Field of size 2
            sage: fr_k
            Ring morphism:
              From: Finite Field of size 2
              To:   Valuation ring at Place (x, x*y)
            sage: to_k
            Ring morphism:
              From: Valuation ring at Place (x, x*y)
              To:   Finite Field of size 2
            sage: to_k(1/y)
            0
            sage: to_k(y/(1+y))
            1
        """
        from .maps import FunctionFieldRingMorphism as morphism

        k, from_k, to_k = self._place._residue_field(name=name)
        mor_from_k = morphism(Hom(k,self), from_k)
        mor_to_k = morphism(Hom(self,k), to_k)
        return k, mor_from_k, mor_to_k
