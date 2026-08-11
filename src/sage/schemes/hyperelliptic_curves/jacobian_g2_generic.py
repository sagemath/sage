r"""
Jacobians of genus-2 hyperelliptic curves

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

from sage.schemes.hyperelliptic_curves import (
    jacobian_g2_homset_inert,
    jacobian_g2_homset_ramified,
    jacobian_g2_homset_split,
)
from sage.schemes.hyperelliptic_curves.jacobian_generic import (
    HyperellipticJacobian_generic,
)
from sage.schemes.hyperelliptic_curves.kummer_surface import KummerSurface


class HyperellipticJacobian_g2_generic(HyperellipticJacobian_generic):
    r"""
    Special class to handle optimisations for jacobian computations
    in genus two
    """

    def _point_homset(self, *args, **kwds):
        r"""
        Create the point Hom-set of the Jacobian of a genus-2 curve.

        TODO: make a constructor for this??

        TESTS::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(-x^6 + 15*x^4 - 75*x^2 -56, x^3 + x)
            sage: J = Jacobian(H)(QQ)
            sage: type(J) # indirect doctest
            <class 'sage.schemes.hyperelliptic_curves.jacobian_g2_homset_inert.HyperellipticJacobianHomsetInert_g2_with_category'>
        """
        H = self.curve()
        if H.is_ramified():
            return jacobian_g2_homset_ramified.HyperellipticJacobianHomsetRamified_g2(
                *args, **kwds
            )
        if H.is_split():
            return jacobian_g2_homset_split.HyperellipticJacobianHomsetSplit_g2(
                *args, **kwds
            )
        return jacobian_g2_homset_inert.HyperellipticJacobianHomsetInert_g2(
            *args, **kwds
        )

    def kummer_surface(self):
        r"""
        Construct the Kummer surface from the Jacobian of a genus-2 curve.

        INPUT: ``jacobian`` -- the Jacobian of a genus-2 curve

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurve(x**5 + x)
            sage: J = Jacobian(H)
            sage: K = KummerSurface(J); K
            Kummer Surface of Jacobian of Hyperelliptic Curve over Finite Field of size 13 defined by y^2 = x^5 + x.
            The defining equation is X0^4 - 4*X0*X1^2*X2 + 2*X0^2*X2^2 + X2^4 - 2*X0^2*X1*X3 - 2*X1*X2^2*X3 + X1^2*X3^2 - 4*X0*X2*X3^2
        """
        return KummerSurface(self)
