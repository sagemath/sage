# sage_setup: distribution = sagemath-schemes
r"""
Plane conic constructor

AUTHORS:

- Marco Streng (2010-07-20)

- Nick Alexander (2008-01-08)
"""
# ****************************************************************************
#       Copyright (C) 2008 Nick Alexander <ncalexander@gmail.com>
#       Copyright (C) 2009/2010 Marco Streng <marco.streng@gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.categories.integral_domains import IntegralDomains
from sage.rings.rational_field import RationalField
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.fraction_field import FractionField_generic

from sage.rings.number_field.number_field_base import NumberField
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.affine.affine_point import SchemeMorphism_point_affine
from sage.schemes.projective.projective_point import SchemeMorphism_point_projective_field
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.structure.all import Sequence
from sage.structure.element import Matrix

from .con_field import ProjectiveConic_field
from .con_finite_field import ProjectiveConic_finite_field
from .con_number_field import ProjectiveConic_number_field
from .con_rational_field import ProjectiveConic_rational_field
from .con_rational_function_field import ProjectiveConic_rational_function_field


def Conic(base_field, F=None, names=None, unique=True):
    r"""
    Return the plane projective conic curve defined by ``F``
    over ``base_field``.

    The input form ``Conic(F, names=None)`` is also accepted,
    in which case the fraction field of the base ring of ``F``
    is used as base field.

    INPUT:

    - ``base_field`` -- the base field of the conic

    - ``names`` -- list, tuple, or comma separated string
      of three variable names specifying the names
      of the coordinate functions of the ambient
      space `\Bold{P}^3`. If not specified or read
      off from ``F``, then this defaults to ``'x,y,z'``.

    - ``F`` -- a polynomial, list, matrix, ternary quadratic form,
      or list or tuple of 5 points in the plane.

                   If ``F`` is a polynomial or quadratic form,
                   then the output is the curve in the projective plane
                   defined by ``F = 0``.

                   If ``F`` is a polynomial, then it must be a polynomial
                   of degree at most 2 in 2 variables, or a homogeneous
                   polynomial in of degree 2 in 3 variables.

                   If ``F`` is a matrix, then the output is the zero locus
                   of `(x,y,z) F (x,y,z)^t`.

                   If ``F`` is a list of coefficients, then it has
                   length 3 or 6 and gives the coefficients of
                   the monomials `x^2, y^2, z^2` or all 6 monomials
                   `x^2, xy, xz, y^2, yz, z^2` in lexicographic order.

                   If ``F`` is a list of 5 points in the plane, then the output
                   is a conic through those points.

    - ``unique`` -- used only if ``F`` is a list of points in the plane;
      if the conic through the points is not unique, then
      raise :exc:`ValueError` if and only if ``unique`` is ``True``

    OUTPUT: a plane projective conic curve defined by ``F`` over a field

    EXAMPLES:

    Conic curves given by polynomials ::

        sage: X,Y,Z = QQ['X,Y,Z'].gens()
        sage: Conic(X^2 - X*Y + Y^2 - Z^2)
        Projective Conic Curve over Rational Field defined by X^2 - X*Y + Y^2 - Z^2
        sage: x,y = GF(7)['x,y'].gens()
        sage: Conic(x^2 - x + 2*y^2 - 3, 'U,V,W')
        Projective Conic Curve over Finite Field of size 7
         defined by U^2 + 2*V^2 - U*W - 3*W^2

    Conic curves given by matrices ::

        sage: Conic(matrix(QQ, [[1, 2, 0], [4, 0, 0], [7, 0, 9]]), 'x,y,z')
        Projective Conic Curve over Rational Field defined by x^2 + 6*x*y + 7*x*z + 9*z^2

        sage: x,y,z = GF(11)['x,y,z'].gens()
        sage: C = Conic(x^2 + y^2 - 2*z^2); C
        Projective Conic Curve over Finite Field of size 11 defined by x^2 + y^2 - 2*z^2
        sage: Conic(C.symmetric_matrix(), 'x,y,z')
        Projective Conic Curve over Finite Field of size 11 defined by x^2 + y^2 - 2*z^2

    Conics given by coefficients ::

        sage: Conic(QQ, [1,2,3])
        Projective Conic Curve over Rational Field defined by x^2 + 2*y^2 + 3*z^2
        sage: Conic(GF(7), [1,2,3,4,5,6], 'X')
        Projective Conic Curve over Finite Field of size 7
        defined by X0^2 + 2*X0*X1 - 3*X1^2 + 3*X0*X2 - 2*X1*X2 - X2^2

    The conic through a set of points ::

        sage: C = Conic(QQ, [[10,2],[3,4],[-7,6],[7,8],[9,10]]); C
        Projective Conic Curve over Rational Field
         defined by x^2 + 13/4*x*y - 17/4*y^2 - 35/2*x*z + 91/4*y*z - 37/2*z^2
        sage: C.rational_point()
        (10 : 2 : 1)
        sage: C.point([3,4])
        (3 : 4 : 1)

        sage: a = AffineSpace(GF(13), 2)
        sage: Conic([a([x,x^2]) for x in range(5)])
        Projective Conic Curve over Finite Field of size 13 defined by x^2 - y*z
    """
    if not (base_field is None or base_field in IntegralDomains()):
        if names is None:
            names = F
        F = base_field
        base_field = None
    if isinstance(F, (list, tuple)):
        if len(F) == 1:
            return Conic(base_field, F[0], names)
        if names is None:
            names = 'x,y,z'
        if len(F) == 5:
            L = []
            for f in F:
                if isinstance(f, SchemeMorphism_point_affine):
                    C = Sequence(f, universe=base_field)
                    if len(C) != 2:
                        raise TypeError("points in F (=%s) must be planar" % F)
                    C.append(1)
                elif isinstance(f, SchemeMorphism_point_projective_field):
                    C = Sequence(f, universe=base_field)
                elif isinstance(f, (list, tuple)):
                    C = Sequence(f, universe=base_field)
                    if len(C) == 2:
                        C.append(1)
                else:
                    raise TypeError("F (=%s) must be a sequence of planar "
                                    "points" % F)
                if len(C) != 3:
                    raise TypeError("points in F (=%s) must be planar" % F)
                P = C.universe()
                if P not in IntegralDomains():
                    raise TypeError("coordinates of points in F (=%s) must "
                                    "be in an integral domain" % F)
                L.append(Sequence([C[0]**2, C[0] * C[1],
                                   C[0] * C[2], C[1]**2,
                                   C[1] * C[2], C[2]**2], P.fraction_field()))
            M = matrix(L)
            if unique and M.rank() != 5:
                raise ValueError("points in F (=%s) do not define a unique "
                                 "conic" % F)
            con = Conic(base_field, Sequence(M.right_kernel().gen()), names)
            con.point(F[0])
            return con
        F = Sequence(F, universe=base_field)
        base_field = F.universe().fraction_field()
        temp_ring = PolynomialRing(base_field, 3, names)
        x, y, z = temp_ring.gens()
        if len(F) == 3:
            return Conic(F[0] * x**2 + F[1] * y**2 + F[2] * z**2)
        if len(F) == 6:
            return Conic(F[0] * x**2 + F[1] * x * y + F[2] * x * z +
                         F[3] * y**2 + F[4] * y * z + F[5] * z**2)
        raise TypeError("F (=%s) must be a sequence of 3 or 6"
                        "coefficients" % F)

    from sage.quadratic_forms.quadratic_form import QuadraticForm

    if isinstance(F, QuadraticForm):
        F = F.matrix()
    if isinstance(F, Matrix) and F.is_square() and F.ncols() == 3:
        if names is None:
            names = 'x,y,z'
        temp_ring = PolynomialRing(F.base_ring(), 3, names)
        F = vector(temp_ring.gens()) * F * vector(temp_ring.gens())

    if not isinstance(F, MPolynomial):
        raise TypeError("F (=%s) must be a three-variable polynomial or "
                        "a sequence of points or coefficients" % F)

    if F.total_degree() != 2:
        raise TypeError("F (=%s) must have degree 2" % F)

    if base_field is None:
        base_field = F.base_ring()
    if base_field not in IntegralDomains():
        raise ValueError(f"Base field (={base_field}) must be a field")
    base_field = base_field.fraction_field()
    if names is None:
        names = F.parent().variable_names()
    pol_ring = PolynomialRing(base_field, 3, names)

    if F.parent().ngens() == 2:
        x, y, z = pol_ring.gens()
        F = pol_ring(F(x / z, y / z) * z**2)

    if F == 0:
        raise ValueError("F must be nonzero over base field %s" % base_field)

    if F.total_degree() != 2:
        raise TypeError("F (=%s) must have degree 2 over base field %s" %
                        (F, base_field))

    if F.parent().ngens() == 3:
        P2 = ProjectiveSpace(2, base_field, names)
        if isinstance(base_field, FiniteField):
            return ProjectiveConic_finite_field(P2, F)
        if isinstance(base_field, RationalField):
            return ProjectiveConic_rational_field(P2, F)
        if isinstance(base_field, NumberField):
            return ProjectiveConic_number_field(P2, F)
        if isinstance(base_field, FractionField_generic) and (isinstance(base_field.ring(), PolynomialRing_general) or isinstance(base_field.ring(), MPolynomialRing_base)):
            return ProjectiveConic_rational_function_field(P2, F)

        return ProjectiveConic_field(P2, F)

    raise TypeError("Number of variables of F (=%s) must be 2 or 3" % F)
