r"""
Set of homomorphisms between two projective schemes

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

- Volker Braun (2011-08-11): significant improvement and refactoring.

- Ben Hutz (June 2012): added support for projective ring
"""


#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ
from sage.schemes.generic.homset import SchemeHomset_points

from sage.rings.rational_field import is_RationalField
from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from copy import copy

#*******************************************************************
# Projective varieties
#*******************************************************************
class SchemeHomset_points_projective_field(SchemeHomset_points):
    """
    Set of rational points of a projective variety over a field.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_homset import SchemeHomset_points_projective_field
        sage: SchemeHomset_points_projective_field(Spec(QQ), ProjectiveSpace(QQ,2))
        Set of rational points of Projective Space of dimension 2 over Rational Field
    """
    def points(self, B=0, prec=53):
        """
        Return some or all rational points of a projective scheme.

        INPUT:

        - ``B`` - integer (optional, default=0). The bound for the
          coordinates.

        - ``prec`` - he precision to use to compute the elements of bounded height for number fields.

        OUTPUT:

        A list of points. Over a finite field, all points are
        returned. Over an infinite field, all points satisfying the
        bound are returned.

        .. WARNING::

           In the current implementation, the output of the [Doyle-Krumm] algorithm
           cannot be guaranteed to be correct due to the necessity of floating point
           computations. In some cases, the default 53-bit precision is
           considerably lower than would be required for the algorithm to
           generate correct output.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: P(QQ).points(4)
            [(-4 : 1), (-3 : 1), (-2 : 1), (-3/2 : 1), (-4/3 : 1), (-1 : 1),
            (-3/4 : 1), (-2/3 : 1), (-1/2 : 1), (-1/3 : 1), (-1/4 : 1), (0 : 1),
            (1/4 : 1), (1/3 : 1), (1/2 : 1), (2/3 : 1), (3/4 : 1), (1 : 0), (1 : 1),
            (4/3 : 1), (3/2 : 1), (2 : 1), (3 : 1), (4 : 1)]

        ::

            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)
            sage: P.<x,y,z> = ProjectiveSpace(K,2)
            sage: len(P(K).points(1.8))
            381

        ::

            sage: P1 = ProjectiveSpace(GF(2),1)
            sage: F.<a> = GF(4,'a')
            sage: P1(F).points()
            [(0 : 1), (1 : 0), (1 : 1), (a : 1), (a + 1 : 1)]

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
            sage: E = P.subscheme([(y^3-y*z^2) - (x^3-x*z^2),(y^3-y*z^2) + (x^3-x*z^2)])
            sage: E(P.base_ring()).points()
            [(-1 : -1 : 1), (-1 : 0 : 1), (-1 : 1 : 1), (0 : -1 : 1), (0 : 0 : 1), (0 : 1 : 1),
            (1 : -1 : 1), (1 : 0 : 1), (1 : 1 : 1)]
        """
        X = self.codomain()
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if not is_ProjectiveSpace(X) and X.base_ring() in Fields():
            #Then it must be a subscheme
            dim_ideal = X.defining_ideal().dimension()
            if dim_ideal < 1: # no points
                return []
            if dim_ideal == 1: # if X zero-dimensional
                rat_points = set()
                PS = X.ambient_space()
                N = PS.dimension_relative()
                BR = X.base_ring()
                #need a lexicographic ordering for elimination
                R = PolynomialRing(BR, N + 1, PS.gens(), order='lex')
                I = R.ideal(X.defining_polynomials())
                I0 = R.ideal(0)
                #Determine the points through elimination
                #This is much faster than using the I.variety() function on each affine chart.
                for k in range(N + 1):
                    #create the elimination ideal for the kth affine patch
                    G = I.substitute({R.gen(k):1}).groebner_basis()
                    if G != [1]:
                        P = {}
                        #keep track that we know the kth coordinate is 1
                        P.update({R.gen(k):1})
                        points = [P]
                        #work backwards from solving each equation for the possible
                        #values of the next coordinate
                        for i in range(len(G) - 1, -1, -1):
                            new_points = []
                            good = 0
                            for P in points:
                                #substitute in our dictionary entry that has the values
                                #of coordinates known so far. This results in a single
                                #variable polynomial (by elimination)
                                L = G[i].substitute(P)
                                if L != 0:
                                    L = L.factor()
                                    #the linear factors give the possible rational values of
                                    #this coordinate
                                    for pol, pow in L:
                                        if pol.degree() == 1 and len(pol.variables()) == 1:
                                            good = 1
                                            r = pol.variables()[0]
                                            varindex = R.gens().index(r)
                                            #add this coordinates information to
                                            #each dictionary entry
                                            P.update({R.gen(varindex):-pol.constant_coefficient() / pol.monomial_coefficient(r)})
                                            new_points.append(copy(P))
                                else:
                                    new_points.append(P)
                                    good = 1
                            if good:
                                points = new_points
                        #the dictionary entries now have values for all coordinates
                        #they are the rational solutions to the equations
                        #make them into projective points
                        for i in range(len(points)):
                            if len(points[i]) == N + 1 and I.subs(points[i]) == I0:
                                S = X([points[i][R.gen(j)] for j in range(N + 1)])
                                S.normalize_coordinates()
                                rat_points.add(S)
                rat_points = sorted(rat_points)
                return rat_points
        R = self.value_ring()
        if is_RationalField(R):
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
            return enum_projective_rational_field(self,B)
        elif R in NumberFields():
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.projective.projective_rational_point import enum_projective_number_field
            return enum_projective_number_field(self,B, prec=prec)
        elif is_FiniteField(R):
            from sage.schemes.projective.projective_rational_point import enum_projective_finite_field
            return enum_projective_finite_field(self.extended_codomain())
        else:
            raise TypeError("unable to enumerate points over %s"%R)

class SchemeHomset_points_projective_ring(SchemeHomset_points):
    """
    Set of rational points of a projective variety over a commutative ring.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_homset import SchemeHomset_points_projective_ring
        sage: SchemeHomset_points_projective_ring(Spec(ZZ), ProjectiveSpace(ZZ,2))
        Set of rational points of Projective Space of dimension 2 over Integer Ring
    """

    def points(self, B=0):
        """
        Return some or all rational points of a projective scheme.

        INPUT:

        - ``B`` -- integer (optional, default=0). The bound for the
          coordinates.

        EXAMPLES::

            sage: from sage.schemes.projective.projective_homset import SchemeHomset_points_projective_ring
            sage: H = SchemeHomset_points_projective_ring(Spec(ZZ), ProjectiveSpace(ZZ,2))
            sage: H.points(3)
            [(0 : 0 : 1), (0 : 1 : -3), (0 : 1 : -2), (0 : 1 : -1), (0 : 1 : 0), (0
            : 1 : 1), (0 : 1 : 2), (0 : 1 : 3), (0 : 2 : -3), (0 : 2 : -1), (0 : 2 :
            1), (0 : 2 : 3), (0 : 3 : -2), (0 : 3 : -1), (0 : 3 : 1), (0 : 3 : 2),
            (1 : -3 : -3), (1 : -3 : -2), (1 : -3 : -1), (1 : -3 : 0), (1 : -3 : 1),
            (1 : -3 : 2), (1 : -3 : 3), (1 : -2 : -3), (1 : -2 : -2), (1 : -2 : -1),
            (1 : -2 : 0), (1 : -2 : 1), (1 : -2 : 2), (1 : -2 : 3), (1 : -1 : -3),
            (1 : -1 : -2), (1 : -1 : -1), (1 : -1 : 0), (1 : -1 : 1), (1 : -1 : 2),
            (1 : -1 : 3), (1 : 0 : -3), (1 : 0 : -2), (1 : 0 : -1), (1 : 0 : 0), (1
            : 0 : 1), (1 : 0 : 2), (1 : 0 : 3), (1 : 1 : -3), (1 : 1 : -2), (1 : 1 :
            -1), (1 : 1 : 0), (1 : 1 : 1), (1 : 1 : 2), (1 : 1 : 3), (1 : 2 : -3),
            (1 : 2 : -2), (1 : 2 : -1), (1 : 2 : 0), (1 : 2 : 1), (1 : 2 : 2), (1 :
            2 : 3), (1 : 3 : -3), (1 : 3 : -2), (1 : 3 : -1), (1 : 3 : 0), (1 : 3 :
            1), (1 : 3 : 2), (1 : 3 : 3), (2 : -3 : -3), (2 : -3 : -2), (2 : -3 :
            -1), (2 : -3 : 0), (2 : -3 : 1), (2 : -3 : 2), (2 : -3 : 3), (2 : -2 :
            -3), (2 : -2 : -1), (2 : -2 : 1), (2 : -2 : 3), (2 : -1 : -3), (2 : -1 :
            -2), (2 : -1 : -1), (2 : -1 : 0), (2 : -1 : 1), (2 : -1 : 2), (2 : -1 :
            3), (2 : 0 : -3), (2 : 0 : -1), (2 : 0 : 1), (2 : 0 : 3), (2 : 1 : -3),
            (2 : 1 : -2), (2 : 1 : -1), (2 : 1 : 0), (2 : 1 : 1), (2 : 1 : 2), (2 :
            1 : 3), (2 : 2 : -3), (2 : 2 : -1), (2 : 2 : 1), (2 : 2 : 3), (2 : 3 :
            -3), (2 : 3 : -2), (2 : 3 : -1), (2 : 3 : 0), (2 : 3 : 1), (2 : 3 : 2),
            (2 : 3 : 3), (3 : -3 : -2), (3 : -3 : -1), (3 : -3 : 1), (3 : -3 : 2),
            (3 : -2 : -3), (3 : -2 : -2), (3 : -2 : -1), (3 : -2 : 0), (3 : -2 : 1),
            (3 : -2 : 2), (3 : -2 : 3), (3 : -1 : -3), (3 : -1 : -2), (3 : -1 : -1),
            (3 : -1 : 0), (3 : -1 : 1), (3 : -1 : 2), (3 : -1 : 3), (3 : 0 : -2), (3
            : 0 : -1), (3 : 0 : 1), (3 : 0 : 2), (3 : 1 : -3), (3 : 1 : -2), (3 : 1
            : -1), (3 : 1 : 0), (3 : 1 : 1), (3 : 1 : 2), (3 : 1 : 3), (3 : 2 : -3),
            (3 : 2 : -2), (3 : 2 : -1), (3 : 2 : 0), (3 : 2 : 1), (3 : 2 : 2), (3 :
            2 : 3), (3 : 3 : -2), (3 : 3 : -1), (3 : 3 : 1), (3 : 3 : 2)]
        """
        R = self.value_ring()
        if R == ZZ:
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
            return enum_projective_rational_field(self,B)
        else:
            raise TypeError("unable to enumerate points over %s"%R)


#*******************************************************************
# Abelian varieties
#*******************************************************************
class SchemeHomset_points_abelian_variety_field(SchemeHomset_points_projective_field):
    r"""
    Set of rational points of an Abelian variety.

    INPUT:

    See :class:`SchemeHomset_generic`.

    TESTS:

    The bug reported at :trac:`1785` is fixed::

        sage: K.<a> = NumberField(x^2 + x - (3^3-3))
        sage: E = EllipticCurve('37a')
        sage: X = E(K)
        sage: X
        Abelian group of points on Elliptic Curve defined by
        y^2 + y = x^3 + (-1)*x over Number Field in a with
        defining polynomial x^2 + x - 24
        sage: P = X([3,a])
        sage: P
        (3 : a : 1)
        sage: P in E
        False
        sage: P in E.base_extend(K)
        True
        sage: P in X.codomain()
        False
        sage: P in X.extended_codomain()
        True
    """

    def _element_constructor_(self, *v, **kwds):
        """
        The element constructor.

        INPUT:

        - ``v`` -- anything that determines a scheme morphism in the
          Hom-set.

        OUTPUT:

        The scheme morphism determined by ``v``.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: X = E(QQ)
            sage: P = X([0,1,0]);  P
            (0 : 1 : 0)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_number_field'>

        TESTS::

            sage: X._element_constructor_([0,1,0])
            (0 : 1 : 0)
        """
        if len(v) == 1:
            v = v[0]
        return self.codomain()._point(self.extended_codomain(), v, **kwds)

    def _repr_(self):
        """
        Return a string representation of this homset.

        OUTPUT:

        String.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: X = E(QQ)
            sage: X._repr_()
            'Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field'
        """
        s = 'Abelian group of points on ' + str(self.extended_codomain())
        return s

    def base_extend(self, R):
        """
        Extend the base ring.

        This is currently not implemented except for the trivial case
        ``R==ZZ``.

        INPUT:

        - ``R`` -- a ring.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: Hom = E.point_homset();  Hom
            Abelian group of points on Elliptic Curve defined
            by y^2 + y = x^3 - x over Rational Field
            sage: Hom.base_ring()
            Integer Ring
            sage: Hom.base_extend(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Abelian variety point sets are not
            implemented as modules over rings other than ZZ
        """
        if R is not ZZ:
            raise NotImplementedError('Abelian variety point sets are not '
                            'implemented as modules over rings other than ZZ')
        return self


from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.schemes.generic.homset',
                           'SchemeHomsetModule_abelian_variety_coordinates_field',
                           SchemeHomset_points_abelian_variety_field)
