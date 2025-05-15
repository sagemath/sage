"""
Rational point sets on a Jacobian of a hyperelliptic curve (inert case)
"""
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism import (
    MumfordDivisorClassFieldInert,
)


class HyperellipticJacobianHomsetInert(HyperellipticJacobianHomset):
    def __init__(self, Y, X, **kwds):
        """
        Create the Jacobian Hom-set of a hyperelliptic curve without
        rational points at infinity.

        TESTS::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(3*x^6 + 2*x^2 + 1)
            sage: assert H.is_inert()
            sage: JK = Jacobian(H)(GF(7))
            sage: type(JK)
            <class 'sage.schemes.hyperelliptic_curves_smooth_model.jacobian_g2_homset_inert.HyperellipticJacobianHomsetInert_g2_with_category'>
        """
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldInert

    def zero(self, check=True):
        """
        Return the zero element of the Jacobian.

        The Mumford presentation of the zero element is given by
        `(1, 0 : g/2)`, `g` is the genus of the hyperelliptic curve.

        NOTE: We require that the genus is even if the hyperelliptic
        curve is inert.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(2*x^6 + 1)
            sage: H.is_inert()
            True
            sage: J = H.jacobian()
            sage: J.zero()
            (1, 0 : 1)

            sage: H = HyperellipticCurveSmoothModel(3*x^10 + 1)
            sage: J = H.jacobian()
            sage: J.zero()
            (1, 0 : 2)
        """
        g = self.curve().genus()
        if g % 2:
            raise ValueError(
                "unable to perform arithmetic for inert models of odd genus"
            )
        R = self.curve().polynomial_ring()
        return self._morphism_element(self, R.one(), R.zero(), check=check)
