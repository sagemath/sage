r"""
Finite residue fields (PARI implementation)
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

from sage.rings.finite_rings.finite_field_pari_ffelt import FiniteField_pari_ffelt
from sage.rings.finite_rings.residue_field import ResidueField_generic, ResidueFieldHomomorphism_global, ReductionMap


class ResidueFiniteField_pari_ffelt(ResidueField_generic, FiniteField_pari_ffelt):
    """
    The class representing residue fields of number fields that have non-prime
    order at least `2^16`.

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^3 - 7)
        sage: P = K.ideal(923478923).factor()[0][0]
        sage: k = K.residue_field(P)
        sage: k.degree()
        2
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(c)
        sage: b + c
        2*abar
        sage: b*c
        664346875*abar + 535606347
        sage: k.base_ring()
        Finite Field of size 923478923

        sage: R.<t> = GF(5)[]
        sage: P = R.ideal(4*t^12 + 3*t^11 + 4*t^10 + t^9 + t^8
        ....:             + 3*t^7 + 2*t^6 + 3*t^4 + t^3 + 3*t^2 + 2)
        sage: k.<a> = P.residue_field()
        sage: type(k)
        <class 'sage.rings.finite_rings.residue_field_pari_ffelt.ResidueFiniteField_pari_ffelt_with_category'>
        sage: k(1/t)
        3*a^11 + a^10 + 3*a^9 + 2*a^8 + 2*a^7 + a^6 + 4*a^5 + a^3 + 2*a^2 + a
    """

    def __init__(self, p, characteristic, name, modulus, to_vs, to_order, PB):
        """
        Initialize ``self``.

        EXAMPLES:

        We create a residue field with implementation ``pari_ffelt``::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 7)
            sage: P = K.ideal(923478923).factor()[0][0]
            sage: type(P.residue_field())
            <class 'sage.rings.finite_rings.residue_field_pari_ffelt.ResidueFiniteField_pari_ffelt_with_category'>
        """
        ResidueField_generic.__init__(self, p)
        FiniteField_pari_ffelt.__init__(self, characteristic, modulus, name)
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
        Coerce ``x`` into ``self``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K.<aa> = NumberField(x^3 - 2)
            sage: P = K.factor(10007)[0][0]
            sage: P.residue_class_degree()
            2
            sage: ff.<alpha> = P.residue_field(); ff
            Residue field in alpha of Fractional ideal (-12*aa^2 + 189*aa - 475)
            sage: type(ff)
            <class 'sage.rings.finite_rings.residue_field_pari_ffelt.ResidueFiniteField_pari_ffelt_with_category'>
            sage: ff(alpha^2 + 1)
            7521*alpha + 4131
            sage: ff(17/3)
            6677
            sage: V = ff.vector_space(map=False); v = V([3,-2])                         # needs sage.modules
            sage: type(ff.convert_map_from(V))                                          # needs sage.modules
            <class 'sage.structure.coerce_maps.DefaultConvertMap_unique'>
            sage: ff(v)  # indirect doctest                                             # needs sage.modules
            10005*alpha + 3

            sage: R.<t> = GF(5)[]; P = R.ideal(4*t^12 + 3*t^11 + 4*t^10 + t^9 + t^8 + 3*t^7 + 2*t^6 + 3*t^4 + t^3 + 3*t^2 + 2)
            sage: k.<a> = P.residue_field()
            sage: V = k.vector_space(map=False); v = V([1,2,3,4,5,6,7,8,9,0,1,2]); k(v)  # indirect doctest             # needs sage.modules
            2*a^11 + a^10 + 4*a^8 + 3*a^7 + 2*a^6 + a^5 + 4*a^3 + 3*a^2 + 2*a + 1
        """
        try:
            return self.element_class(self, x)
        except TypeError:
            return ResidueField_generic._element_constructor_(self, x)
