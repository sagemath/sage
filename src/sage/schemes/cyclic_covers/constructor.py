# sage_setup: distribution = sagemath-schemes
"""
Cyclic cover curve constructor
"""

# *****************************************************************************
#  Copyright (C) 2018 Edgar Costa <edgarcosta@math.dartmouth.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.schemes.affine.affine_space import AffineSpace

from .cycliccover_generic import CyclicCover_generic
from .cycliccover_finite_field import CyclicCover_finite_field


def CyclicCover(r, f, names=None, check_smooth=True):
    r"""
    Return the cyclic cover of the projective line given by `y^r = f`, for
    a univariate polynomial `f`.

    INPUT:

    - ``r`` -- the order of the cover

    - ``f`` -- univariate polynomial if not given, then it defaults to 0

    - ``names`` -- (default: ``["x","y"]``) names for the coordinate functions

    - ``check_squarefree`` -- boolean (default: ``True``); test if
      the input defines a unramified cover of the projective line

    .. WARNING::

        When setting ``check_smooth=False`` or using a base ring that is
        not a field, the output curves are not to be trusted. For example, the
        output of ``is_singular`` or ``is_smooth`` only tests smoothness over
        the field of fractions.

    .. NOTE::

        The words "cyclic cover" are usually used for covers of degree
        greater than two.
        We usually refer to smooth double covers of the projective line as
        "hyperelliptic curves" or "elliptic curves" if the genus is one.
        We allow such cases in this implementation, but we highly recommend
        to use the more specific constructors/classes HyperellipticCurve and
        EllipticCurve for a wider range of tools.

    EXAMPLES:

    Basic examples::

        sage: R.<x> = QQ[]
        sage: CyclicCover(2, x^5 + x + 1)
        Cyclic Cover of P^1 over Rational Field defined by y^2 = x^5 + x + 1
        sage: CyclicCover(3, x^5 + x + 1)
        Cyclic Cover of P^1 over Rational Field defined by y^3 = x^5 + x + 1
        sage: CyclicCover(5, x^5 + x + 1)
        Cyclic Cover of P^1 over Rational Field defined by y^5 = x^5 + x + 1
        sage: CyclicCover(15, x^9 + x + 1)
        Cyclic Cover of P^1 over Rational Field defined by y^15 = x^9 + x + 1

        sage: k.<a> = GF(9); R.<x> = k[]                                                # needs sage.rings.finite_rings
        sage: CyclicCover(5, x^9 + x + 1)                                               # needs sage.rings.finite_rings
        Cyclic Cover of P^1 over Finite Field in a of size 3^2
         defined by y^5 = x^9 + x + 1
        sage: CyclicCover(15, x^9 + x + 1)                                              # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: As the characteristic divides the order of the cover,
        this model is not smooth.

    We can change the names of the variables in the output::

        sage: k.<a> = GF(9); R.<x> = k[]                                                # needs sage.rings.finite_rings
        sage: CyclicCover(5, x^9 + x + 1, names=["A","B"])                              # needs sage.rings.finite_rings
        Cyclic Cover of P^1 over Finite Field in a of size 3^2
         defined by B^5 = A^9 + A + 1

    Double roots::

        sage: P.<x> = GF(7)[]
        sage: CyclicCover(2, (x^3-x+2)^2*(x^6-1))
        Traceback (most recent call last):
        ...
        ValueError: Not a smooth Cyclic Cover of P^1: singularity in the provided affine patch.

        sage: CyclicCover(2, (x^3-x+2)^2*(x^6-1), check_smooth=False)                   # needs sage.rings.finite_rings
        Cyclic Cover of P^1 over Finite Field of size 7
         defined by y^2 = x^12 - 2*x^10 - 3*x^9 + x^8 + 3*x^7 + 3*x^6
                          + 2*x^4 + 3*x^3 - x^2 - 3*x + 3


    Input with integer coefficients creates objects with the integers
    as base ring, but only checks smoothness over `\QQ`, not over Spec(`\ZZ`).
    In other words, it is checked that the discriminant is nonzero, but it is
    not checked whether the discriminant is a unit in `\ZZ^*`::

        sage: R.<x> = ZZ[]
        sage: CyclicCover(5, (x^3-x+2)*(x^6-1))
        Cyclic Cover of P^1 over Integer Ring defined by y^5 = x^9 - x^7 + 2*x^6 - x^3 + x - 2
    """
    if not isinstance(f, Polynomial):
        raise TypeError("Arguments f (= %s) must be a polynomial" % (f,))
    P = f.parent()
    f = P(f)
    if check_smooth:
        if P(r) == 0:
            raise ValueError(
                "As the characteristic divides the order of the cover, "
                "this model is not smooth."
            )

        try:
            smooth = f.is_squarefree()
        except NotImplementedError as err:
            raise NotImplementedError(
                str(err) + "Use " "check_smooth=False to skip this check."
            )
        if not smooth:
            raise ValueError(
                "Not a smooth Cyclic Cover of P^1: "
                "singularity in the provided affine patch."
            )
    R = P.base_ring()
    if names is None:
        names = ["x", "y"]
    A2 = AffineSpace(2, R, names=names)

    if isinstance(R, FiniteField):
        return CyclicCover_finite_field(A2, r, f, names=names)
    else:
        return CyclicCover_generic(A2, r, f, names=names)
