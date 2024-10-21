# sage_setup: distribution = sagemath-schemes
"""
Plane quartic curves over a general ring

These are generic genus 3 curves, as distinct from hyperelliptic curves of genus 3.

EXAMPLES::

    sage: PP.<X,Y,Z> = ProjectiveSpace(2, QQ)
    sage: f = X^4 + Y^4 + Z^4 - 3*X*Y*Z*(X+Y+Z)
    sage: C = QuarticCurve(f); C
    Quartic Curve over Rational Field
     defined by X^4 + Y^4 - 3*X^2*Y*Z - 3*X*Y^2*Z - 3*X*Y*Z^2 + Z^4
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#*****************************************************************************


import sage.schemes.curves.projective_curve as projective_curve


def is_QuarticCurve(C):
    """
    Check whether ``C`` is a Quartic Curve.

    EXAMPLES::

        sage: from sage.schemes.plane_quartics.quartic_generic import is_QuarticCurve
        sage: x,y,z = PolynomialRing(QQ, ['x','y','z']).gens()
        sage: Q = QuarticCurve(x**4 + y**4 + z**4)
        sage: is_QuarticCurve(Q)
        doctest:warning...
        DeprecationWarning: The function is_QuarticCurve is deprecated; use 'isinstance(..., QuarticCurve_generic)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_QuarticCurve is deprecated; use 'isinstance(..., QuarticCurve_generic)' instead.")
    return isinstance(C, QuarticCurve_generic)


class QuarticCurve_generic(projective_curve.ProjectivePlaneCurve):
    # DRK: Note that we should check whether the curve is

    def _repr_type(self):
        """
        Return the representation of ``self``.

        EXAMPLES::

            sage: x,y,z = PolynomialRing(QQ, ['x','y','z']).gens()
            sage: Q = QuarticCurve(x**4 + y**4 + z**4)
            sage: Q._repr_type()
            'Quartic'
        """
        return "Quartic"

    def genus(self):
        """
        Return the genus of ``self``.

        EXAMPLES::

            sage: x,y,z = PolynomialRing(QQ, ['x','y','z']).gens()
            sage: Q = QuarticCurve(x**4 + y**4 + z**4)
            sage: Q.genus()
            3
        """
        return 3
