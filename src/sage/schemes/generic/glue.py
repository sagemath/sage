"""
Scheme obtained by gluing two other schemes
"""

#*******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#*******************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.schemes.generic.scheme import Scheme

lazy_import('sage.schemes.generic.morphism', 'SchemeMorphism')


class GluedScheme(Scheme):
    r"""
    INPUT:

    - ``f`` -- open immersion from a scheme `U` to a scheme `X`

    - ``g`` -- open immersion from `U` to a scheme `Y`

    OUTPUT: the scheme obtained by gluing `X` and `Y` along the open set `U`

    .. NOTE::

       Checking that `f` and `g` are open
       immersions is not implemented.

    EXAMPLES::

        sage: # needs sage.libs.singular
        sage: R.<x, y> = QQ[]
        sage: S.<xbar, ybar> = R.quotient(x * y - 1)
        sage: Rx = QQ["x"]
        sage: Ry = QQ["y"]
        sage: phi_x = Rx.hom([xbar])
        sage: phi_y = Ry.hom([ybar])
        sage: Sx = Schemes()(phi_x)
        sage: Sy = Schemes()(phi_y)
        sage: Sx.glue_along_domains(Sy)
        Scheme obtained by gluing X and Y along U, where
          X: Spectrum of Univariate Polynomial Ring in x over Rational Field
          Y: Spectrum of Univariate Polynomial Ring in y over Rational Field
          U: Spectrum of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y - 1)
    """
    def __init__(self, f, g, check=True):
        if check:
            if not isinstance(f, SchemeMorphism):
                raise TypeError("f (=%s) must be a scheme morphism" % f)
            if not isinstance(g, SchemeMorphism):
                raise TypeError("g (=%s) must be a scheme morphism" % g)
            if f.domain() != g.domain():
                raise ValueError("f (=%s) and g (=%s) must have the same domain" % (f,g))
        self.__f = f
        self.__g = g

    def gluing_maps(self):
        r"""
        Return the gluing maps of this glued scheme, i.e. the maps `f` and `g`.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<x, y> = QQ[]
            sage: S.<xbar, ybar> = R.quotient(x * y - 1)
            sage: Rx = QQ["x"]
            sage: Ry = QQ["y"]
            sage: phi_x = Rx.hom([xbar])
            sage: phi_y = Ry.hom([ybar])
            sage: Sx = Schemes()(phi_x)
            sage: Sy = Schemes()(phi_y)
            sage: Sxy = Sx.glue_along_domains(Sy)
            sage: Sxy.gluing_maps() == (Sx, Sy)
            True
        """
        return self.__f, self.__g

    def _repr_(self):
        return "Scheme obtained by gluing X and Y along U, where\n  X: %s\n  Y: %s\n  U: %s" % (
            self.__f.codomain(), self.__g.codomain(), self.__f.domain())
