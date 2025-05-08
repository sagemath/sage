from sage.schemes.curves.curve import Curve_generic
from sage.schemes.weighted_projective.weighted_projective_space import WeightedProjectiveSpace_ring


class WeightedProjectiveCurve(Curve_generic):
    """
    Curves in weighted projective spaces.

    EXAMPLES:

    We construct a hyperelliptic curve manually::

        sage: P.<x, y, z> = WeightedProjectiveSpace([1, 3, 1], QQ)
        sage: C = P.curve(y^2 - x^5 * z - 3 * x^2 * z^4 - 2 * z^6); C
    """
    def __init__(self, A, X, *kwargs):
        if not isinstance(A, WeightedProjectiveSpace_ring):
            raise TypeError(f"A(={A}) is not a weighted projective space")
        super().__init__(A, X, *kwargs)

    def curve(self):
        return
