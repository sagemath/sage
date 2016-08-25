"""
General curve constructors.

AUTHORS:

- William Stein (2005-11-13)

- David Kohel (2006-01)
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField

from sage.structure.all import Sequence

from sage.schemes.affine.affine_space import is_AffineSpace
from sage.schemes.generic.ambient_space import is_AmbientSpace
from sage.schemes.generic.algebraic_scheme import is_AlgebraicScheme
from sage.schemes.projective.projective_space import is_ProjectiveSpace

from sage.schemes.affine.all import AffineSpace

from sage.schemes.projective.all import ProjectiveSpace


from .projective_curve import (ProjectivePlaneCurve,
                              ProjectiveCurve,
                              ProjectivePlaneCurve_finite_field,
                              ProjectivePlaneCurve_prime_finite_field)

from .affine_curve import (AffinePlaneCurve,
                          AffineCurve,
                          AffinePlaneCurve_finite_field,
                          AffinePlaneCurve_prime_finite_field)

from sage.schemes.plane_conics.constructor import Conic

def Curve(F, A=None):
    """
    Return the plane or space curve defined by ``F``, where
    ``F`` can be either a multivariate polynomial, a list or
    tuple of polynomials, or an algebraic scheme.

    If no ambient space is passed in for ``A``, and if ``F`` is not
    an algebraic scheme, a new ambient space is constructed.

    Also not specifying an ambient space will cause the curve to be defined
    in either affine or projective space based on properties of ``F``. In
    particular, if ``F`` contains a nonhomogenous polynomial, the curve is
    affine, and if ``F`` consists of homogenous polynomials, then the curve
    is projective.

    INPUT:

    - ``F`` -- a multivariate polynomial, or a list or tuple of polynomials, or an algebraic scheme.

    - ``A`` -- (default: None) an ambient space in which to create the curve.

    EXAMPLE: A projective plane curve

    ::

        sage: x,y,z = QQ['x,y,z'].gens()
        sage: C = Curve(x^3 + y^3 + z^3); C
        Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
        sage: C.genus()
        1

    EXAMPLE: Affine plane curves

    ::

        sage: x,y = GF(7)['x,y'].gens()
        sage: C = Curve(y^2 + x^3 + x^10); C
        Affine Plane Curve over Finite Field of size 7 defined by x^10 + x^3 + y^2
        sage: C.genus()
        0
        sage: x, y = QQ['x,y'].gens()
        sage: Curve(x^3 + y^3 + 1)
        Affine Plane Curve over Rational Field defined by x^3 + y^3 + 1

    EXAMPLE: A projective space curve

    ::

        sage: x,y,z,w = QQ['x,y,z,w'].gens()
        sage: C = Curve([x^3 + y^3 - z^3 - w^3, x^5 - y*z^4]); C
        Projective Curve over Rational Field defined by x^3 + y^3 - z^3 - w^3, x^5 - y*z^4
        sage: C.genus()
        13

    EXAMPLE: An affine space curve

    ::

        sage: x,y,z = QQ['x,y,z'].gens()
        sage: C = Curve([y^2 + x^3 + x^10 + z^7,  x^2 + y^2]); C
        Affine Curve over Rational Field defined by x^10 + z^7 + x^3 + y^2, x^2 + y^2
        sage: C.genus()
        47

    EXAMPLE: We can also make non-reduced non-irreducible curves.

    ::

        sage: x,y,z = QQ['x,y,z'].gens()
        sage: Curve((x-y)*(x+y))
        Projective Conic Curve over Rational Field defined by x^2 - y^2
        sage: Curve((x-y)^2*(x+y)^2)
        Projective Plane Curve over Rational Field defined by x^4 - 2*x^2*y^2 + y^4

    EXAMPLE: A union of curves is a curve.

    ::

        sage: x,y,z = QQ['x,y,z'].gens()
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: D = Curve(x^4 + y^4 + z^4)
        sage: C.union(D)
        Projective Plane Curve over Rational Field defined by
        x^7 + x^4*y^3 + x^3*y^4 + y^7 + x^4*z^3 + y^4*z^3 + x^3*z^4 + y^3*z^4 + z^7

    The intersection is not a curve, though it is a scheme.

    ::

        sage: X = C.intersection(D); X
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
         x^3 + y^3 + z^3,
         x^4 + y^4 + z^4

    Note that the intersection has dimension `0`.

    ::

        sage: X.dimension()
        0
        sage: I = X.defining_ideal(); I
        Ideal (x^3 + y^3 + z^3, x^4 + y^4 + z^4) of Multivariate Polynomial Ring in x, y, z over Rational Field

    EXAMPLE: In three variables, the defining equation must be
    homogeneous.

    If the parent polynomial ring is in three variables, then the
    defining ideal must be homogeneous.

    ::

        sage: x,y,z = QQ['x,y,z'].gens()
        sage: Curve(x^2+y^2)
        Projective Conic Curve over Rational Field defined by x^2 + y^2
        sage: Curve(x^2+y^2+z)
        Traceback (most recent call last):
        ...
        TypeError: x^2 + y^2 + z is not a homogeneous polynomial

    The defining polynomial must always be nonzero::

        sage: P1.<x,y> = ProjectiveSpace(1,GF(5))
        sage: Curve(0*x)
        Traceback (most recent call last):
        ...
        ValueError: defining polynomial of curve must be nonzero

    ::

        sage: A.<x,y,z> = AffineSpace(QQ, 3)
        sage: C = Curve([y - x^2, z - x^3], A)
        sage: A == C.ambient_space()
        True
    """
    if not A is None:
        if not isinstance(F, (list, tuple)):
            return Curve([F], A)
        if not is_AmbientSpace(A):
            raise TypeError("A (=%s) must be either an affine or projective space"%A)
        if not all([f.parent() == A.coordinate_ring() for f in F]):
            raise TypeError("F (=%s) must be a list or tuple of polynomials of the coordinate ring of " \
            "A (=%s)"%(F, A))
        n = A.dimension_relative()
        if n < 2:
            raise TypeError("A (=%s) must be either an affine or projective space of dimension > 1"%A)
        # there is no dimension check when initializing a plane curve, so check here that F consists
        # of a single nonconstant polynomial
        if n == 2:
            if len(F) != 1 or F[0] == 0 or not is_MPolynomial(F[0]):
                raise TypeError("F (=%s) must consist of a single nonconstant polynomial to define a plane curve"%(F,))
        if is_AffineSpace(A):
            if n > 2:
                return AffineCurve(A, F)
            k = A.base_ring()
            if is_FiniteField(k):
                if k.is_prime_field():
                    return AffinePlaneCurve_prime_finite_field(A, F[0])
                return AffinePlaneCurve_finite_field(A, F[0])
            return AffinePlaneCurve(A, F[0])
        elif is_ProjectiveSpace(A):
            if not all([f.is_homogeneous() for f in F]):
                raise TypeError("polynomials defining a curve in a projective space must be homogeneous")
            if n > 2:
                return ProjectiveCurve(A, F)
            k = A.base_ring()
            if is_FiniteField(k):
                if k.is_prime_field():
                    return ProjectivePlaneCurve_prime_finite_field(A, F[0])
                return ProjectivePlaneCurve_finite_field(A, F[0])
            return ProjectivePlaneCurve(A, F[0])

    if is_AlgebraicScheme(F):
        return Curve(F.defining_polynomials(), F.ambient_space())

    if isinstance(F, (list, tuple)):
        if len(F) == 1:
            return Curve(F[0])
        F = Sequence(F)
        P = F.universe()
        if not is_MPolynomialRing(P):
            raise TypeError("universe of F must be a multivariate polynomial ring")

        for f in F:
            if not f.is_homogeneous():
                A = AffineSpace(P.ngens(), P.base_ring())
                A._coordinate_ring = P
                return AffineCurve(A, F)

        A = ProjectiveSpace(P.ngens()-1, P.base_ring())
        A._coordinate_ring = P
        return ProjectiveCurve(A, F)

    if not is_MPolynomial(F):
        raise TypeError("F (=%s) must be a multivariate polynomial"%F)

    P = F.parent()
    k = F.base_ring()
    if F.parent().ngens() == 2:
        if F == 0:
            raise ValueError("defining polynomial of curve must be nonzero")
        A2 = AffineSpace(2, P.base_ring())
        A2._coordinate_ring = P

        if is_FiniteField(k):
            if k.is_prime_field():
                return AffinePlaneCurve_prime_finite_field(A2, F)
            else:
                return AffinePlaneCurve_finite_field(A2, F)
        else:
            return AffinePlaneCurve(A2, F)

    elif F.parent().ngens() == 3:
        if F == 0:
            raise ValueError("defining polynomial of curve must be nonzero")
        P2 = ProjectiveSpace(2, P.base_ring())
        P2._coordinate_ring = P

        if F.total_degree() == 2 and k.is_field():
            return Conic(F)

        if is_FiniteField(k):
            if k.is_prime_field():
                return ProjectivePlaneCurve_prime_finite_field(P2, F)
            else:
                return ProjectivePlaneCurve_finite_field(P2, F)
        else:
            return ProjectivePlaneCurve(P2, F)


    else:

        raise TypeError("Number of variables of F (=%s) must be 2 or 3"%F)


