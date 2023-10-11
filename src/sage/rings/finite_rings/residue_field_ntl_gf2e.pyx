r"""
Finite residue fields (NTL implementation)
"""

# *****************************************************************************
#       Copyright (C) 2007-2019 David Roe <roed@math.harvard.edu>
#                     2007      William Stein <wstein@gmail.com>
#                     2008      John Cremona
#                     2008      Robert Bradshaw
#                     2009      Nick Alexander
#                     2010      Robert L. Miller
#                     2010-2013 Simon King
#                     2010-2017 Jeroen Demeyer
#                     2012      Travis Scrimshaw
#                     2016-2021 Frédéric Chapoton
#                     2021-2022 Antonio Rojas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.finite_rings.finite_field_ntl_gf2e import FiniteField_ntl_gf2e
from sage.rings.finite_rings.residue_field import ResidueField_generic, ResidueFieldHomomorphism_global, ReductionMap


class ResidueFiniteField_ntl_gf2e(ResidueField_generic, FiniteField_ntl_gf2e):
    """
    The class representing residue fields with order a power of 2.

    When the order is less than `2^16`, givaro is used by default instead.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^3 - 7)
        sage: P = K.ideal(29).factor()[0][0]
        sage: k = K.residue_field(P)
        sage: k.degree()
        2
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(c)
        sage: b*c^2
        7
        sage: b*c
        13*abar + 5

        sage: R.<t> = GF(2)[]; P = R.ideal(t^19 + t^5 + t^2 + t + 1)
        sage: k.<a> = R.residue_field(P); type(k)
        <class 'sage.rings.finite_rings.residue_field_ntl_gf2e.ResidueFiniteField_ntl_gf2e_with_category'>
        sage: k(1/t)
        a^18 + a^4 + a + 1
        sage: k(1/t)*t
        1
    """
    # we change the order for consistency with FiniteField_ntl_gf2e's __cinit__
    def __init__(self, q, name, modulus, repr, p, to_vs, to_order, PB):
        r"""
        INPUT:

        - ``p`` -- the prime ideal defining this residue field

        - ``q`` -- the order of this residue field

        - ``name`` -- the name of the generator of this extension

        - ``modulus`` -- the polynomial modulus for this extension

        - ``to_vs`` -- the map from the number field (or function field) to
          the appropriate vector space (over `\QQ` or `F_p(t)`)

        - ``to_order`` -- the map from a lattice in that vector space to the
          maximal order

        - ``PB`` -- a matrix used in defining the reduction and lifting maps

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4 + 3*x^2 - 17)                                 # needs sage.rings.number_field
            sage: P = K.ideal(61).factor()[0][0]                                        # needs sage.rings.number_field
            sage: k = K.residue_field(P)                                                # needs sage.rings.number_field

            sage: R.<t> = GF(3)[]; P = R.ideal(t^4 - t^3 + t + 1); k.<a> = P.residue_field(); type(k)
            <class 'sage.rings.finite_rings.residue_field_givaro.ResidueFiniteField_givaro_with_category'>
            sage: a^5
            a^3 + 2*a^2 + a + 2
        """
        ResidueField_generic.__init__(self, p)
        FiniteField_ntl_gf2e.__init__(self, q, name, modulus, repr)
        K = OK = p.ring()
        if OK.is_field():
            OK = OK.ring_of_integers()
        else:
            K = K.fraction_field()
        if PB is None:
            PBinv = None
        else:
            PBinv = PB**(-1)
        self._populate_coercion_lists_(coerce_list=[self.base_ring(),
                                                    ResidueFieldHomomorphism_global(OK, self, to_vs, to_order, PB, PBinv)],
                                       convert_list=[ReductionMap(K, self, to_vs, to_order, PB, PBinv)])

    def _element_constructor_(self, x):
        """
        INPUT:

        - ``x`` -- Something to cast into ``self``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4 + 3*x^2 - 17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: k(77*a^7 + 4)
            2*abar + 4
            sage: V = k.vector_space(map=False); v = V([3,-2])
            sage: type(k.convert_map_from(V))
            <class 'sage.structure.coerce_maps.DefaultConvertMap_unique'>
            sage: k(v) # indirect doctest
            59*abar + 3

            sage: R.<t> = GF(3)[]; P = R.ideal(t^4 - t^3 + t + 1); k.<a> = P.residue_field()
            sage: V = k.vector_space(map=False); v = V([0,1,2,3])
            sage: k(v) # indirect doctest
            2*a^2 + a
        """
        try:
            return FiniteField_ntl_gf2e._element_constructor_(self, x)
        except TypeError:
            return ResidueField_generic._element_constructor_(self, x)
