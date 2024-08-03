"""
Points on schemes
"""

#*******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.structure.element import Element
from sage.structure.richcmp import richcmp

########################################################
# Base class for points on a scheme, either topological
# or defined by a morphism.
########################################################

class SchemePoint(Element):
    """
    Base class for points on a scheme, either topological or defined
    by a morphism.
    """
    def __init__(self, S, parent=None):
        """
        INPUT:

        - ``S`` -- a scheme

        - ``parent`` -- the parent in which to construct this point

        TESTS::

            sage: from sage.schemes.generic.point import SchemePoint
            sage: S = Spec(ZZ)
            sage: P = SchemePoint(S); P
            Point on Spectrum of Integer Ring
        """
        Element.__init__(self, parent)
        self.__S = S

    def scheme(self):
        """
        Return the scheme on which ``self`` is a point.

        EXAMPLES::

            sage: from sage.schemes.generic.point import SchemePoint
            sage: S = Spec(ZZ)
            sage: P = SchemePoint(S)
            sage: P.scheme()
            Spectrum of Integer Ring
        """
        return self.__S

    def _repr_(self):
        """
        Return a string representation of this generic scheme point.

        TESTS::

            sage: from sage.schemes.generic.point import SchemePoint
            sage: S = Spec(ZZ)
            sage: P = SchemePoint(S); P
            Point on Spectrum of Integer Ring
            sage: P._repr_()
            'Point on Spectrum of Integer Ring'
        """
        return "Point on %s" % self.__S

########################################################
# Topological points on a scheme
########################################################

def is_SchemeTopologicalPoint(x):
    from sage.misc.superseded import deprecation
    deprecation(38296,
                "The function is_SchemeTopologicalPoint is deprecated; "
                "use 'isinstance(..., SchemeTopologicalPoint)' instead.")
    return isinstance(x, SchemeTopologicalPoint)

class SchemeTopologicalPoint(SchemePoint):
    """
    Base class for topological points on schemes.
    """
    def __init__(self, S):
        """
        INPUT:

        - ``S`` -- a scheme

        TESTS:

        The parent of a topological point is the scheme on which it
        lies (see :issue:`7946`)::

            sage: R = Zmod(8)
            sage: S = Spec(R)
            sage: x = S(R.ideal(2))
            sage: isinstance(x, sage.schemes.generic.point.SchemeTopologicalPoint)
            True
            sage: x.parent() is S
            True
        """
        SchemePoint.__init__(self, S, parent=S)


class SchemeTopologicalPoint_affine_open(SchemeTopologicalPoint):
    def __init__(self, u, x):
        """
        INPUT:

        - ``u`` -- morphism with domain an affine scheme `U`

        - ``x`` -- topological point on `U`
        """
        SchemeTopologicalPoint.__init__(self, u.codomain())
        self.__u = u
        self.__x = x

    def _repr_(self):
        return "Point on %s defined by x in U, where:\n  U: %s\n  x: %s" % (
            self.scheme(), self.embedding_of_affine_open().domain(),
            self.point_on_affine())

    def point_on_affine(self):
        """
        Return the scheme point on the affine open `U`.
        """
        return self.__x

    def affine_open(self):
        """
        Return the affine open subset `U`.
        """
        return self.__u.domain()

    def embedding_of_affine_open(self):
        """
        Return the embedding from the affine open subset `U` into this scheme.
        """
        return self.__u


class SchemeTopologicalPoint_prime_ideal(SchemeTopologicalPoint):
    def __init__(self, S, P, check=False):
        """
        INPUT:

        - ``S`` -- an affine scheme

        - ``P`` -- a prime ideal of the coordinate ring of `S`, or
          anything that can be converted into such an ideal

        TESTS::

            sage: from sage.schemes.generic.point import SchemeTopologicalPoint_prime_ideal
            sage: S = Spec(ZZ)
            sage: P = SchemeTopologicalPoint_prime_ideal(S, 3); P
            Point on Spectrum of Integer Ring defined by the Principal ideal (3) of Integer Ring
            sage: SchemeTopologicalPoint_prime_ideal(S, 6, check=True)
            Traceback (most recent call last):
            ...
            ValueError: The argument Principal ideal (6) of Integer Ring must be a prime ideal of Integer Ring
            sage: SchemeTopologicalPoint_prime_ideal(S, ZZ.ideal(7))
            Point on Spectrum of Integer Ring defined by the Principal ideal (7) of Integer Ring

        We define a parabola in the projective plane as a point
        corresponding to a prime ideal::

            sage: P2.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: SchemeTopologicalPoint_prime_ideal(P2, y*z - x^2)
            Point on Projective Space of dimension 2 over Rational Field defined by
             the Ideal (-x^2 + y*z) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        R = S.coordinate_ring()
        from sage.rings.ideal import Ideal_generic
        if not isinstance(P, Ideal_generic):
            P = R.ideal(P)
        elif P.ring() is not R:
            P = R.ideal(P.gens())
        # ideally we would have check=True by default, but
        # unfortunately is_prime() is only implemented in a small
        # number of cases
        if check and not P.is_prime():
            raise ValueError("The argument %s must be a prime ideal of %s" % (P, R))
        SchemeTopologicalPoint.__init__(self, S)
        self.__P = P

    def _repr_(self):
        """
        Return a string representation of this scheme point.

        TESTS::

            sage: from sage.schemes.generic.point import SchemeTopologicalPoint_prime_ideal
            sage: P2.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: pt = SchemeTopologicalPoint_prime_ideal(P2, y*z - x^2); pt
            Point on Projective Space of dimension 2 over Rational Field defined by
             the Ideal (-x^2 + y*z) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: pt._repr_()
            'Point on Projective Space of dimension 2 over Rational Field defined by the Ideal (-x^2 + y*z) of Multivariate Polynomial Ring in x, y, z over Rational Field'
        """
        return "Point on %s defined by the %s" % (self.scheme(),
                                                self.prime_ideal())

    def prime_ideal(self):
        """
        Return the prime ideal that defines this scheme point.

        EXAMPLES::

            sage: from sage.schemes.generic.point import SchemeTopologicalPoint_prime_ideal
            sage: P2.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: pt = SchemeTopologicalPoint_prime_ideal(P2, y*z - x^2)
            sage: pt.prime_ideal()
            Ideal (-x^2 + y*z) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.__P

    def _richcmp_(self, other, op):
        """
        Compare ``self`` to ``other``.

        TESTS::

            sage: S = Spec(ZZ)
            sage: x = S(ZZ.ideal(5))
            sage: y = S(ZZ.ideal(7))
            sage: x == y
            False
        """
        return richcmp(self.__P, other.__P, op)
