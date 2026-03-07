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

AUTHORS:

- Gareth Ma (2025)
"""

# ****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@gmail.com>
#  Copyright (C) 2025 Gareth Ma <grhkm21@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

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

    def _repr_type(self) -> str:
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

    def projective_curve(self):
        r"""
        Return this weighted projective curve as a projective curve.

        A weighted homogeneous polynomial `f(x_1, \ldots, x_n)`, where `x_i` has
        weight `w_i`, can be viewed as an unweighted homogeneous polynomial
        `f(y_1^{w_1}, \ldots, y_n^{w_n})`. This correspondence extends to
        varieties.

        .. TODO:

            Implement homsets for weighted projective spaces and implement this
            as a ``projective_embedding`` method instead.

        EXAMPLES::

            sage: WP = WeightedProjectiveSpace([1, 3, 1], QQ, "x, y, z")
            sage: x, y, z = WP.gens()
            sage: C = WP.curve(y^2 - (x^5*z + 3*x^2*z^4 - 2*x*z^5 + 4*z^6)); C
            Weighted Projective Curve over Rational Field defined by y^2 - x^5*z - 3*x^2*z^4 + 2*x*z^5 - 4*z^6
            sage: C.projective_curve()
            Projective Plane Curve over Rational Field defined by y^6 - x^5*z - 3*x^2*z^4 + 2*x*z^5 - 4*z^6
        """
        from sage.schemes.projective.projective_space import ProjectiveSpace

        WP = self.ambient_space()
        PP = ProjectiveSpace(WP.dimension_relative(), WP.base_ring(), WP.variable_names())
        PP_ring = PP.coordinate_ring()
        subs_dict = {name: var**weight for (name, var), weight in
                     zip(WP.gens_dict().items(), WP.weights())}

        wp_polys = self.defining_polynomials()
        pp_polys = [PP_ring(poly.subs(**subs_dict)) for poly in wp_polys]

        return PP.curve(pp_polys)
