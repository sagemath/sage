from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism import (
    MumfordDivisorClassFieldInert,
)


class HyperellipticJacobianHomsetInert(HyperellipticJacobianHomset):
    def __init__(self, Y, X, **kwds):
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldInert

    def zero(self, check=True):
        """
        Return the zero element of the Jacobian
        """
        g = self.curve().genus()
        if g % 2:
            raise ValueError(
                "unable to perform arithmetic for inert models of odd genus"
            )
        R = self.curve().polynomial_ring()
        return self._morphism_element(self, R.one(), R.zero(), check=check)
