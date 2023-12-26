# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.rings.finite_rings       (because all doctests use finite fields)
"""
Places of function fields: rational
"""

# ****************************************************************************
#       Copyright (C) 2016-2022 Kwankyu Lee <ekwankyu@gmail.com>
#                     2019      Brent Baccala
#                     2021      Jonathan Kliem
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from .place import FunctionFieldPlace


class FunctionFieldPlace_rational(FunctionFieldPlace):
    """
    Places of rational function fields.
    """
    def degree(self):
        """
        Return the degree of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: i = O.ideal(x^2 + x + 1)
            sage: p = i.place()
            sage: p.degree()
            2
        """
        if self.is_infinite_place():
            return 1
        else:
            return self._prime.gen().numerator().degree()

    def is_infinite_place(self):
        """
        Return ``True`` if the place is at infinite.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: F.places()
            [Place (1/x), Place (x), Place (x + 1)]
            sage: [p.is_infinite_place() for p in F.places()]
            [True, False, False]
        """
        F = self.function_field()
        return self.prime_ideal().ring() == F.maximal_order_infinite()

    def local_uniformizer(self):
        """
        Return a local uniformizer of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: F.places()
            [Place (1/x), Place (x), Place (x + 1)]
            sage: [p.local_uniformizer() for p in F.places()]
            [1/x, x, x + 1]
        """
        return self.prime_ideal().gen()

    def residue_field(self, name=None):
        """
        Return the residue field of the place.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: p = O.ideal(x^2 + x + 1).place()
            sage: k, fr_k, to_k = p.residue_field()                                     # needs sage.rings.function_field
            sage: k                                                                     # needs sage.rings.function_field
            Finite Field in z2 of size 2^2
            sage: fr_k                                                                  # needs sage.rings.function_field
            Ring morphism:
              From: Finite Field in z2 of size 2^2
              To:   Valuation ring at Place (x^2 + x + 1)
            sage: to_k                                                                  # needs sage.rings.function_field
            Ring morphism:
              From: Valuation ring at Place (x^2 + x + 1)
              To:   Finite Field in z2 of size 2^2
        """
        return self.valuation_ring().residue_field(name=name)

    def _residue_field(self, name=None):
        """
        Return the residue field of the place along with the maps from
        and to it.

        INPUT:

        - ``name`` -- string; name of the generator of the residue field

        EXAMPLES::

            sage: # needs sage.modules
            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: i = O.ideal(x^2 + x + 1)
            sage: p = i.place()
            sage: R, fr, to = p._residue_field()
            sage: R
            Finite Field in z2 of size 2^2
            sage: [fr(e) for e in R.list()]
            [0, x, x + 1, 1]
            sage: to(x*(x+1)) == to(x) * to(x+1)
            True
        """
        F = self.function_field()
        prime = self.prime_ideal()

        if self.is_infinite_place():
            K = F.constant_base_field()

            def from_K(e):
                return F(e)

            def to_K(f):
                n = f.numerator()
                d = f.denominator()

                n_deg = n.degree()
                d_deg = d.degree()

                if n_deg < d_deg:
                    return K(0)
                elif n_deg == d_deg:
                    return n.lc() / d.lc()
                else:
                    raise TypeError("not in the valuation ring")
        else:
            O = F.maximal_order()
            K, from_K, _to_K = O._residue_field(prime, name=name)

            def to_K(f):
                if f in O: # f.denominator() is 1
                    return _to_K(f.numerator())
                else:
                    d = F(f.denominator())
                    n = d * f

                    nv = prime.valuation(O.ideal(n))
                    dv = prime.valuation(O.ideal(d))

                    if nv > dv:
                        return K(0)
                    elif dv > nv:
                        raise TypeError("not in the valuation ring")

                    s = ~prime.gen()
                    rd = d * s**dv # in O but not in prime
                    rn = n * s**nv # in O but not in prime
                    return to_K(rn) / to_K(rd)

        return K, from_K, to_K

    def valuation_ring(self):
        """
        Return the valuation ring at the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)                                # needs sage.rings.function_field
            sage: p = L.places_finite()[0]                                              # needs sage.rings.function_field
            sage: p.valuation_ring()                                                    # needs sage.rings.function_field
            Valuation ring at Place (x, x*y)
        """
        from .valuation_ring import FunctionFieldValuationRing

        return FunctionFieldValuationRing(self.function_field(), self)
