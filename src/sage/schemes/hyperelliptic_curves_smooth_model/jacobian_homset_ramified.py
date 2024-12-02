from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_generic import (
    HyperellipticJacobianHomset,
)
from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism import (
    MumfordDivisorClassFieldRamified,
)


class HyperellipticJacobianHomsetRamified(HyperellipticJacobianHomset):
    def __init__(self, Y, X, **kwds):
        super().__init__(Y, X, **kwds)
        self._morphism_element = MumfordDivisorClassFieldRamified
