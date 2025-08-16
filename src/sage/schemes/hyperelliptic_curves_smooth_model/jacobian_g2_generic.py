"""
Jacobians of genus-2 curves
"""

from sage.schemes.hyperelliptic_curves_smooth_model import (
    jacobian_g2_homset_inert,
    jacobian_g2_homset_ramified,
    jacobian_g2_homset_split,
)
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_generic import (
    HyperellipticJacobian_generic,
)


class HyperellipticJacobian_g2_generic(HyperellipticJacobian_generic):
    """
    Special class to handle optimisations for jacobian computations
    in genus two
    """

    def _point_homset(self, *args, **kwds):
        """
        Create the point Hom-set of the Jacobian of a genus-2 curve.

        TODO: make a constructor for this??

        TESTS::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(-x^6 + 15*x^4 - 75*x^2 -56, x^3 + x)
            sage: J = Jacobian(H)(QQ)
            sage: type(J) # indirect doctest
            <class 'sage.schemes.hyperelliptic_curves_smooth_model.jacobian_g2_homset_inert.HyperellipticJacobianHomsetInert_g2_with_category'>
        """
        H = self.curve()
        if H.is_ramified():
            return jacobian_g2_homset_ramified.HyperellipticJacobianHomsetRamified_g2(
                *args, **kwds
            )
        elif H.is_split():
            return jacobian_g2_homset_split.HyperellipticJacobianHomsetSplit_g2(
                *args, **kwds
            )
        return jacobian_g2_homset_inert.HyperellipticJacobianHomsetInert_g2(
            *args, **kwds
        )
