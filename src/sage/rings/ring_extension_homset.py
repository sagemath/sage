r"""
Homset between extensions of rings

AUTHOR:

- Xavier Caruso (2019)
"""

#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.rings.homset import RingHomset_generic
from sage.rings.ring_extension_morphism import RingExtensionHomomorphism

class RingExtensionHomset(RingHomset_generic):
    r"""
    A generic class for homsets between ring extensions.

    TESTS::

        sage: K = GF(5^2).over()
        sage: L = GF(5^8).over(K)
        sage: H = Hom(K,L)
        sage: H
        Set of Homomorphisms from Field in z2 with defining polynomial x^2 + 4*x + 2 over its base to Field in z8 with defining polynomial x^4 + (3 - z2)*x + z2 over its base

        sage: type(H)
        <... 'sage.rings.ring_extension_homset.RingExtensionHomset_with_category'>
    """
    def __call__(self, *args, **kwargs):
        r"""
        Return the morphism in this parent defined by the
        given parameters.

        TESTS::

            sage: K.<a> = GF(5^2).over()
            sage: L.<b> = GF(5^4).over(K)
            sage: Hom(L,L)([b^5, a^5])
            Ring endomorphism of Field in b with defining polynomial x^2 + (3 - a)*x + a over its base
              Defn: b |--> (2 + a) + 2*b
                    with map on base ring:
                    a |--> 1 - a

            sage: Hom(K,L)(GF(5^4).coerce_map_from(GF(5^2)))
            Ring morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Field in b with defining polynomial x^2 + (3 - a)*x + a over its base
              Defn: a |--> a
        """
        return RingExtensionHomomorphism(self, *args, **kwargs)
