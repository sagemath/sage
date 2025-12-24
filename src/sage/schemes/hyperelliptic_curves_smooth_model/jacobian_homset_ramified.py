r"""
Rational point sets on a Jacobian of a hyperelliptic curve (ramified case)

AUTHORS:

- Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
"""

# ****************************************************************************
#       Copyright (C) 2025 Sabrina Kunzweiler, Gareth Ma, Giacomo Pope
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism import (
    MumfordDivisorClassFieldRamified,
)


class HyperellipticJacobianHomsetRamified(HyperellipticJacobianHomset):
    def __init__(self, Y, X, **kwds):
        """
        Create the Jacobian Hom-set of a hyperelliptic curve with
        precisely one rational points at infinity.

        TESTS::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + 2*x^2 + 1)
            sage: assert H.is_ramified()
            sage: JK = Jacobian(H)(GF(7))
            sage: type(JK)
            <class 'sage.schemes.hyperelliptic_curves_smooth_model.jacobian_g2_homset_ramified.HyperellipticJacobianHomsetRamified_g2_with_category'>
        """
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldRamified
