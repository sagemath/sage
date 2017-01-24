r"""
Set of homomorphisms between two affine schemes

For schemes `X` and `Y`, this module implements the set of morphisms
`Hom(X,Y)`. This is done by :class:`SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a
scheme. Recall that the `K`-rational points of a scheme `X` over `k`
can be identified with the set of morphisms `Spec(K) \to X`. In Sage
the rational points are implemented by such scheme morphisms. This is
done by :class:`SchemeHomset_points` and its subclasses.

.. note::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- William Stein (2006): initial version.
"""


#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.all import ZZ
from sage.rings.rational_field import is_RationalField
from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
import sage.schemes.generic.homset

#*******************************************************************
# Affine varieties
#*******************************************************************
class SchemeHomset_points_spec(sage.schemes.generic.homset.SchemeHomset_generic):
    """
    Set of rational points of an affine variety.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_spec
        sage: SchemeHomset_points_spec(Spec(QQ), Spec(QQ))
        Set of rational points of Spectrum of Rational Field
    """

    def _element_constructor_(self, *args, **kwds):
        """
        The element constructor.

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: ring_hom = QQ.hom((1,), QQ);  ring_hom
            Ring endomorphism of Rational Field
              Defn: 1 |--> 1
            sage: H = X.Hom(X)
            sage: H(ring_hom)
            Affine Scheme endomorphism of Spectrum of Rational Field
              Defn: Ring endomorphism of Rational Field
                      Defn: 1 |--> 1

        TESTS::

            sage: H._element_constructor_(ring_hom)
            Affine Scheme endomorphism of Spectrum of Rational Field
              Defn: Ring endomorphism of Rational Field
                      Defn: 1 |--> 1
        """
        return sage.schemes.generic.homset.SchemeHomset_generic._element_constructor_(self, *args, **kwds)

    def _repr_(self):
        """
        Return a string representation of a homset.

        OUTPUT: A string.

        EXAMPLES::

            sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_spec
            sage: S = SchemeHomset_points_spec(Spec(QQ), Spec(QQ))
            sage: S._repr_()
            'Set of rational points of Spectrum of Rational Field'
        """
        return 'Set of rational points of '+str(self.codomain())



#*******************************************************************
# Affine varieties
#*******************************************************************
class SchemeHomset_points_affine(sage.schemes.generic.homset.SchemeHomset_points):
    """
    Set of rational points of an affine variety.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_affine
        sage: SchemeHomset_points_affine(Spec(QQ), AffineSpace(ZZ,2))
        Set of rational points of Affine Space of dimension 2 over Rational Field
    """

    def points(self, B=0):
        r"""
        Return some or all rational points of an affine scheme.

        INPUT:

        - ``B`` -- integer (optional, default: 0). The bound for the
          height of the coordinates.

        OUTPUT:

        - If the base ring is a finite field: all points of the scheme,
          given by coordinate tuples.

        - If the base ring is `\QQ` or `\ZZ`: the subset of points whose
          coordinates have height ``B`` or less.

        EXAMPLES: The bug reported at #11526 is fixed::

            sage: A2 = AffineSpace(ZZ, 2)
            sage: F = GF(3)
            sage: A2(F).points()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

            sage: R = ZZ
            sage: A.<x,y> = R[]
            sage: I = A.ideal(x^2-y^2-1)
            sage: V = AffineSpace(R, 2)
            sage: X = V.subscheme(I)
            sage: M = X(R)
            sage: M.points(1)
            [(-1, 0), (1, 0)]

        ::

            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: len(A(K).points(9))
            361

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: E = A.subscheme([x^2 + y^2 - 1, y^2 - x^3 + x^2 + x - 1])
            sage: E(A.base_ring()).points()
            [(-1, 0), (0, -1), (0, 1), (1, 0)]
        """
        X = self.codomain()

        from sage.schemes.affine.affine_space import is_AffineSpace
        if not is_AffineSpace(X) and X.base_ring() in Fields():
            # Then X must be a subscheme
            dim_ideal = X.defining_ideal().dimension()
            if dim_ideal < 0: # no points
                return []
            if dim_ideal == 0: # if X zero-dimensional
                N = len(X.ambient_space().gens())
                S = X.defining_polynomials()[0].parent()
                R = PolynomialRing(S.base_ring(), 's', N, order='lex')
                phi = S.hom(R.gens(),R)
                J = R.ideal([phi(t) for t in X.defining_polynomials()])
                D = J.variety()
                points = []
                for d in D:
                    P = [d[t] for t in R.gens()]
                    points.append(X(P))
                points.sort()
                return points
        R = self.value_ring()
        if is_RationalField(R) or R == ZZ:
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.affine.affine_rational_point import enum_affine_rational_field
            return enum_affine_rational_field(self,B)
        if R in NumberFields():
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.affine.affine_rational_point import enum_affine_number_field
            return enum_affine_number_field(self,B)
        elif is_FiniteField(R):
            from sage.schemes.affine.affine_rational_point import enum_affine_finite_field
            return enum_affine_finite_field(self)
        else:
            raise TypeError("unable to enumerate points over %s"%R)
