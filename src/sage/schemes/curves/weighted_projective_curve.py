# sage.doctest: needs sage.libs.singular
r"""
Weighted projective curves

Weighted projective curves in Sage are curves in a weighted projective space or
a weighted projective plane.

EXAMPLES:

For now, only curves in weighted projective plane is supported::

    sage: WP.<x, y, z> = WeightedProjectiveSpace([1, 3, 1], QQ)
    sage: C1 = WP.curve(y^2 - x^5 * z - 3 * x^2 * z^4 - 2 * z^6); C1
    Weighted Projective Curve over Rational Field defined by y^2 - x^5*z - 3*x^2*z^4 - 2*z^6
    sage: C2 = Curve(y^2 - x^5 * z - 3 * x^2 * z^4 - 2 * z^6, WP); C2
    Weighted Projective Curve over Rational Field defined by y^2 - x^5*z - 3*x^2*z^4 - 2*z^6
    sage: C1 == C2
    True
"""

from sage.schemes.curves.curve import Curve_generic
from sage.schemes.weighted_projective.weighted_projective_space import WeightedProjectiveSpace_ring


class WeightedProjectiveCurve(Curve_generic):
    """
    Curves in weighted projective spaces.

    EXAMPLES:

    We construct a hyperelliptic curve manually::

        sage: WP.<x, y, z> = WeightedProjectiveSpace([1, 3, 1], QQ)
        sage: C = Curve(y^2 - x^5 * z - 3 * x^2 * z^4 - 2 * z^6, WP); C
        Weighted Projective Curve over Rational Field defined by y^2 - x^5*z - 3*x^2*z^4 - 2*z^6
    """
    def __init__(self, A, X, *kwargs):
        if not isinstance(A, WeightedProjectiveSpace_ring):
            raise TypeError(f"A(={A}) is not a weighted projective space")
        super().__init__(A, X, *kwargs)

    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: WP.<x,y,z> = WeightedProjectiveSpace([1, 3, 1], QQ)
            sage: C = Curve(y^2 - x^5 * z - 3 * x^2 * z^4 - 2 * z^6, WP); C
            Weighted Projective Curve over Rational Field defined by y^2 - x^5*z - 3*x^2*z^4 - 2*z^6
            sage: C._repr_type()
            'Weighted Projective'
        """
        return "Weighted Projective"

    def curve(self):
        return self
