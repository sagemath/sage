# sage.doctest: needs sage.schemes
r"""
Hypersurfaces in affine and projective space

AUTHORS:

- William Stein <wstein@gmail.com> (2005-12-08)
- David Kohel <kohel@maths.usyd.edu.au> (2005-12-08)
- Alex Ghitza <aghitza@alum.mit.edu> (2009-04-17)
"""

#*****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine
from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective


def is_Hypersurface(self):
    """
    Return ``True`` if ``self`` is a hypersurface, i.e. an object of the type
    :class:`ProjectiveHypersurface` or :class:`AffineHypersurface`.

    EXAMPLES::

        sage: from sage.schemes.generic.hypersurface import is_Hypersurface
        sage: R.<x, y, z> = ZZ[]
        sage: H = ProjectiveHypersurface(x*z + y^2)
        sage: is_Hypersurface(H)
        doctest:warning...
        DeprecationWarning: The function is_Hypersurface is deprecated; use 'isinstance(..., (ProjectiveHypersurface, AffineHypersurface))' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        True

    ::

        sage: H = AffineHypersurface(x*z + y^2)
        sage: is_Hypersurface(H)
        True

    ::

        sage: H = ProjectiveSpace(QQ, 5)
        sage: is_Hypersurface(H)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_Hypersurface is deprecated; use 'isinstance(..., (ProjectiveHypersurface, AffineHypersurface))' instead.")
    return isinstance(self, (ProjectiveHypersurface, AffineHypersurface))


class ProjectiveHypersurface(AlgebraicScheme_subscheme_projective):
    """
    The projective hypersurface defined by the given polynomial.

    EXAMPLES::

        sage: P.<x, y, z> = ProjectiveSpace(ZZ, 2)
        sage: ProjectiveHypersurface(x - y, P)
        Projective hypersurface defined by x - y
         in Projective Space of dimension 2 over Integer Ring

    ::

        sage: R.<x, y, z> = QQ[]
        sage: ProjectiveHypersurface(x - y)
        Projective hypersurface defined by x - y
         in Projective Space of dimension 2 over Rational Field
    """

    def __init__(self, poly, ambient=None):
        """
        Return the projective hypersurface in the space ambient
        defined by the polynomial ``poly``.

        If ambient is not given, it will be constructed based on ``poly``.

        EXAMPLES::

            sage: P.<x, y, z> = ProjectiveSpace(ZZ, 2)
            sage: ProjectiveHypersurface(x - y, P)
            Projective hypersurface defined by x - y in Projective Space of dimension 2 over Integer Ring

        ::

            sage: R.<x, y, z> = QQ[]
            sage: ProjectiveHypersurface(x - y)
            Projective hypersurface defined by x - y in Projective Space of dimension 2 over Rational Field

        TESTS::

            sage: H = ProjectiveHypersurface(x - y)
            sage: H == loads(dumps(H))
            True
        """
        if not isinstance(poly, MPolynomial):
            raise TypeError("Defining polynomial (=%s) must be a multivariate polynomial." % poly)
        if not poly.is_homogeneous():
            raise TypeError("Defining polynomial (=%s) must be homogeneous." % poly)
        if ambient is None:
            R = poly.parent()
            from sage.schemes.projective.projective_space import ProjectiveSpace
            ambient = ProjectiveSpace(R.base_ring(), R.ngens()-1)
            ambient._coordinate_ring = R
        AlgebraicScheme_subscheme_projective.__init__(self, ambient, [poly])

    def _repr_(self):
        """
        Return a string representation of this projective
        hypersurface.

        EXAMPLES::

            sage: R.<x, y, z> = ZZ[]
            sage: H = ProjectiveHypersurface(x*z + y^2)
            sage: H
            Projective hypersurface defined by y^2 + x*z
             in Projective Space of dimension 2 over Integer Ring
            sage: H._repr_()
            'Projective hypersurface defined by y^2 + x*z in Projective Space of dimension 2 over Integer Ring'
        """
        return "Projective hypersurface defined by %s in %s" % (
            self.defining_polynomial(), self.ambient_space())

    def defining_polynomial(self):
        """
        Return the polynomial equation that cuts out this projective
        hypersurface.

        EXAMPLES::

            sage: R.<x, y, z> = ZZ[]
            sage: H = ProjectiveHypersurface(x*z + y^2)
            sage: H.defining_polynomial()
            y^2 + x*z
        """
        return self.defining_polynomials()[0]


class AffineHypersurface(AlgebraicScheme_subscheme_affine):
    """
    The affine hypersurface defined by the given polynomial.

    EXAMPLES::

        sage: A.<x, y, z> = AffineSpace(ZZ, 3)
        sage: AffineHypersurface(x*y - z^3, A)
        Affine hypersurface defined by -z^3 + x*y
         in Affine Space of dimension 3 over Integer Ring

    ::

        sage: A.<x, y, z> = QQ[]
        sage: AffineHypersurface(x*y - z^3)
        Affine hypersurface defined by -z^3 + x*y
         in Affine Space of dimension 3 over Rational Field
    """
    def __init__(self, poly, ambient=None):
        """
        Return the affine hypersurface in the space ambient
        defined by the polynomial poly.

        If ambient is not given, it will be constructed based on
        poly.

        EXAMPLES::

            sage: A.<x, y, z> = AffineSpace(ZZ, 3)
            sage: AffineHypersurface(x*y - z^3, A)
            Affine hypersurface defined by -z^3 + x*y
             in Affine Space of dimension 3 over Integer Ring

        ::

            sage: A.<x, y, z> = QQ[]
            sage: AffineHypersurface(x*y - z^3)
            Affine hypersurface defined by -z^3 + x*y
             in Affine Space of dimension 3 over Rational Field

        TESTS::

            sage: H = AffineHypersurface(x*y - z^3)
            sage: H == loads(dumps(H))
            True
        """
        if not isinstance(poly, MPolynomial):
            raise TypeError("Defining polynomial (= %s) must be a multivariate polynomial" % poly)
        if ambient is None:
            R = poly.parent()
            from sage.schemes.affine.affine_space import AffineSpace
            ambient = AffineSpace(R.base_ring(), R.ngens())
            ambient._coordinate_ring = R
        AlgebraicScheme_subscheme_affine.__init__(self, ambient, [poly])

    def _repr_(self):
        """
        Return a string representation of this affine
        hypersurface.

        EXAMPLES::

            sage: R.<x, y, z> = ZZ[]
            sage: H = AffineHypersurface(x*z + y^2)
            sage: H
            Affine hypersurface defined by y^2 + x*z
             in Affine Space of dimension 3 over Integer Ring
            sage: H._repr_()
            'Affine hypersurface defined by y^2 + x*z in Affine Space of dimension 3 over Integer Ring'
        """
        return "Affine hypersurface defined by %s in %s" % (
            self.defining_polynomial(), self.ambient_space())

    def defining_polynomial(self):
        """
        Return the polynomial equation that cuts out this affine
        hypersurface.

        EXAMPLES::

            sage: R.<x, y, z> = ZZ[]
            sage: H = AffineHypersurface(x*z + y^2)
            sage: H.defining_polynomial()
            y^2 + x*z
        """
        return self.defining_polynomials()[0]
