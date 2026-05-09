r"""
Rational point sets on a Jacobian of a hyperelliptic curve (ramified case)

AUTHORS:

- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2025 Sabrina Kunzweiler <sabrina.kunzweiler@math.u-bordeaux.fr>
#                     2025 Gareth Ma <grhkm21@gmail.com>
#                     2025 Giacomo Pope <giacomopope@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.hyperelliptic_curves.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves.jacobian_morphism import (
    MumfordDivisorClassFieldRamified,
)


class HyperellipticJacobianHomsetRamified(HyperellipticJacobianHomset):
    Element = MumfordDivisorClassFieldRamified

    def __init__(self, Y, X, **kwds):
        r"""
        Create the Jacobian Hom-set of a hyperelliptic curve with
        precisely one rational points at infinity.

        TESTS::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(x^5 + 2*x^2 + 1)
            sage: assert H.is_ramified()
            sage: JK = Jacobian(H)(GF(7))
            sage: type(JK)
            <class 'sage.schemes.hyperelliptic_curves.jacobian_g2_homset_ramified.HyperellipticJacobianHomsetRamified_g2_with_category'>
            sage: TestSuite(JK).run(skip='_test_elements')
        """
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldRamified
