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

- Ben Hutz (2018): add numerical point support
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

from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.rings.cc import CC
from sage.schemes.generic.homset import SchemeHomset_points

from sage.misc.verbose import verbose

from sage.rings.rational_field import is_RationalField
from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
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
    def points(self, **kwds):
        """
        Return some or all rational points of a projective scheme.

        For dimension 0 subschemes points are determined through a groebner
        basis calculation. For schemes or subschemes with dimension greater than 1
        points are determined through enumeration up to the specified bound.

        INPUT:

        kwds:

        - ``bound`` - real number (optional, default=0). The bound for the coordinates for
          subschemes with dimension at least 1.

        - ``precision`` - integer (optional, default=53). The precision to use to
          compute the elements of bounded height for number fields.

        - ``point_tolerance`` - positive real number (optional, default=10^(-10)).
          For numerically inexact fields, two points are considered the same
          if their coordinates are within tolerance.

        - ``zero_tolerance`` - positive real number (optional, default=10^(-10)).
          For numerically inexact fields, points are on the subscheme if they
          satisfy the equations to within tolerance.

        - ``tolerance`` - a rational number in (0,1] used in doyle-krumm algorithm-4
          for enumeration over number fields.

        OUTPUT:

        - a list of rational points of a projective scheme

        .. WARNING::
        
            For numerically inexact fields such as ComplexField or RealField the
            list of points returned is very likely to be incomplete. It may also
            contain repeated points due to tolerances.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: P(QQ).points(bound=4)
            [(-4 : 1), (-3 : 1), (-2 : 1), (-3/2 : 1), (-4/3 : 1), (-1 : 1),
            (-3/4 : 1), (-2/3 : 1), (-1/2 : 1), (-1/3 : 1), (-1/4 : 1), (0 : 1),
            (1/4 : 1), (1/3 : 1), (1/2 : 1), (2/3 : 1), (3/4 : 1), (1 : 0), (1 : 1),
            (4/3 : 1), (3/2 : 1), (2 : 1), (3 : 1), (4 : 1)]

        ::

            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)
            sage: P.<x,y,z> = ProjectiveSpace(K,2)
            sage: len(P(K).points(bound=1.8))
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

        ::

            sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
            sage: E = P.subscheme([y^3 - x^3 - x*z^2, x*y*z])
            sage: L=E(P.base_ring()).points(); sorted(L, key=str)
            verbose 0 (...: projective_homset.py, points) Warning: computations in the numerical fields are inexact;points may be computed partially or incorrectly.
            [(-0.500000000000000 + 0.866025403784439*I : 1.00000000000000 : 0.000000000000000),
            (-0.500000000000000 - 0.866025403784439*I : 1.00000000000000 : 0.000000000000000),
            (-1.00000000000000*I : 0.000000000000000 : 1.00000000000000),
            (0.000000000000000 : 0.000000000000000 : 1.00000000000000),
            (1.00000000000000 : 1.00000000000000 : 0.000000000000000),
            (1.00000000000000*I : 0.000000000000000 : 1.00000000000000)]
            sage: L[0].codomain()
            Projective Space of dimension 2 over Complex Field with 53 bits of precision

        ::

            sage: P.<x,y,z> = ProjectiveSpace(CDF, 2)
            sage: E = P.subscheme([y^2 + x^2 + z^2, x*y*z])
            sage: len(E(P.base_ring()).points())
            verbose 0 (...: projective_homset.py, points) Warning: computations in the numerical fields are inexact;points may be computed partially or incorrectly.
            6
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        X = self.codomain()
        if not is_ProjectiveSpace(X) and X.base_ring() in Fields():
            if hasattr(X.base_ring(), 'precision'):
                numerical = True
                verbose("Warning: computations in the numerical fields are inexact;points may be computed partially or incorrectly.", level=0)
                pt_tol = RR(kwds.pop('point_tolerance', 10**(-10)))
                zero_tol = RR(kwds.pop('zero_tolerance', 10**(-10)))
                if pt_tol <= 0 or zero_tol <= 0:
                    raise ValueError("tolerance must be positive")
            else:
                numerical = False
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
                R = PolynomialRing(BR, N + 1, PS.variable_names(), order='lex')
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
                                if R(L).degree() > 0:
                                    if numerical:
                                        for pol in L.univariate_polynomial().roots(multiplicities=False):
                                            good = 1
                                            r = L.variables()[0]
                                            varindex = R.gens().index(r)
                                            P.update({R.gen(varindex):pol})
                                            new_points.append(copy(P))
                                    else:
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
                            if numerical:
                                if len(points[i]) == N + 1:
                                    S = PS([points[i][R.gen(j)] for j in range(N + 1)])
                                    S.normalize_coordinates()
                                    if all(g(list(S)) < zero_tol for g in X.defining_polynomials()):
                                        rat_points.add(S)
                            else:
                                if len(points[i]) == N + 1 and I.subs(points[i]) == I0:
                                    S = X([points[i][R.gen(j)] for j in range(N + 1)])
                                    S.normalize_coordinates()
                                    rat_points.add(S)

                # remove duplicate element using tolerance
                if numerical:
                    dupl_points = list(rat_points)
                    for i in range(len(dupl_points)):
                        u = dupl_points[i]
                        for j in range(i+1, len(dupl_points)):
                            v = dupl_points[j]
                            if all((u[k] - v[k]).abs() < pt_tol
                                   for k in range(len(u))):
                                rat_points.remove(u)
                                break

                rat_points = sorted(rat_points)
                return rat_points
        R = self.value_ring()
        B = kwds.pop('bound', 0)
        tol = kwds.pop('tolerance', 1e-2)
        prec = kwds.pop('precision', 53)
        if is_RationalField(R):
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            if isinstance(X, AlgebraicScheme_subscheme): # sieve should only be called for subschemes
                from sage.schemes.projective.projective_rational_point import sieve
                return sieve(X, B)
            else:
                from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
                return enum_projective_rational_field(self, B)
        elif R in NumberFields():
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.projective.projective_rational_point import enum_projective_number_field
            return enum_projective_number_field(self, bound=B, tolerance=tol, precision=prec)
        elif is_FiniteField(R):
            from sage.schemes.projective.projective_rational_point import enum_projective_finite_field
            return enum_projective_finite_field(self.extended_codomain())
        else:
            raise TypeError("unable to enumerate points over %s"%R)

    def numerical_points(self, F=None, **kwds):
        """
        Return some or all numerical approximations of rational points of a projective scheme.

        This is for dimension 0 subschemes only and the points are determined
        through a groebner calculation over the base ring and then numerically
        approximating the roots of the resulting polynomials. If the base ring
        is a number field, the embedding into ``F`` must be known.

        INPUT:

        ``F`` - numerical ring

        kwds:

        - ``point_tolerance`` - positive real number (optional, default=10^(-10)).
          For numerically inexact fields, two points are considered the same
          if their coordinates are within tolerance.

        - ``zero_tolerance`` - positive real number (optional, default=10^(-10)).
          For numerically inexact fields, points are on the subscheme if they
          satisfy the equations to within tolerance.

        OUTPUT: A list of points in the ambient space.

        .. WARNING::

           For numerically inexact fields the list of points returned may contain repeated
           or be missing points due to tolerance.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: E = P.subscheme([y^3 - x^3 - x*z^2, x*y*z])
            sage: L = E(QQ).numerical_points(F=RR); L
            [(0.000000000000000 : 0.000000000000000 : 1.00000000000000),
             (1.00000000000000 : 1.00000000000000 : 0.000000000000000)]
            sage: L[0].codomain()
            Projective Space of dimension 2 over Real Field with 53 bits of precision

        ::

            sage: S.<a> = QQ[]
            sage: K.<v> = NumberField(a^5 - 7, embedding=CC((7)**(1/5)))
            sage: P.<x,y,z> = ProjectiveSpace(K,2)
            sage: X = P.subscheme([x^2 - v^2*z^2, y-v*z])
            sage: len(X(K).numerical_points(F=CDF))
            2

        ::

            sage: P.<x1, x2, x3> = ProjectiveSpace(QQ, 2)
            sage: E = P.subscheme([3000*x1^50 + 9875643*x2^2*x3^48 + 12334545*x2^50, x1 + x2])
            sage: len(E(P.base_ring()).numerical_points(F=CDF, zero_tolerance=1e-6))
            49

        TESTS::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: E = P.subscheme([y^3 - x^3 - x*z^2, x*y*z])
            sage: E(QQ).numerical_points(F=CDF, point_tolerance=-1)
            Traceback (most recent call last):
            ...
            ValueError: tolerance must be positive

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: E = P.subscheme([y^3 - x^3 - x*z^2, x*y*z])
            sage: E(QQ).numerical_points(F=CC, zero_tolerance=-1)
            Traceback (most recent call last):
            ...
            ValueError: tolerance must be positive

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: E = P.subscheme([y^3 - x^3 - x*z^2, x*y*z])
            sage: E(QQ).numerical_points(F=QQbar)
            Traceback (most recent call last):
            ...
            TypeError: F must be a numerical field
        """
        from sage.schemes.projective.projective_space import is_ProjectiveSpace
        if F is None:
            F = CC
        if F not in Fields() or not hasattr(F, 'precision'):
            raise TypeError('F must be a numerical field')
        X = self.codomain()
        if X.base_ring() not in NumberFields():
            raise TypeError('base ring must be a number field')

        PP = X.ambient_space().change_ring(F)
        if not is_ProjectiveSpace(X) and X.base_ring() in Fields():
            #Then it must be a subscheme
            dim_ideal = X.defining_ideal().dimension()
            if dim_ideal < 1: # no points
                return []
            if dim_ideal == 1: # if X zero-dimensional
                pt_tol = RR(kwds.pop('point_tolerance', 10**(-10)))
                zero_tol = RR(kwds.pop('zero_tolerance', 10**(-10)))
                if pt_tol <= 0 or zero_tol <= 0:
                    raise ValueError("tolerance must be positive")
                rat_points = set()
                PS = X.ambient_space()
                N = PS.dimension_relative()
                BR = X.base_ring()
                #need a lexicographic ordering for elimination
                R = PolynomialRing(BR, N + 1, PS.variable_names(), order='lex')
                RF = R.change_ring(F)
                I = R.ideal(X.defining_polynomials())
                #Determine the points through elimination
                #This is much faster than using the I.variety() function on each affine chart.
                for k in range(N + 1):
                    #create the elimination ideal for the kth affine patch
                    G = I.substitute({R.gen(k):1}).groebner_basis()
                    G = [RF(g) for g in G]
                    if G != [1]:
                        P = {}
                        #keep track that we know the kth coordinate is 1
                        P.update({RF.gen(k):1})
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
                                if len(RF(L).variables())==1:
                                    for pol in L.univariate_polynomial().roots(ring=F, multiplicities=False):
                                        r = L.variables()[0]
                                        varindex = RF.gens().index(r)
                                        P.update({RF.gen(varindex):pol})
                                        new_points.append(copy(P))
                                        good = 1
                                else:
                                    new_points.append(P)
                                    good = 1
                            if good:
                                points = new_points
                        #the dictionary entries now have values for all coordinates
                        #they are approximate solutions to the equations
                        #make them into projective points
                        polys = [g.change_ring(F) for g in X.defining_polynomials()]
                        for i in range(len(points)):
                            if len(points[i]) == N + 1:
                                S = PP([points[i][RF.gen(j)] for j in range(N + 1)])
                                S.normalize_coordinates()
                                if all(g(list(S)) < zero_tol for g in polys):
                                    rat_points.add(S)
                        # remove duplicate element using tolerance
                        #since they are normalized we can just compare coefficients
                        dupl_points = list(rat_points)
                        for i in range(len(dupl_points)):
                            u = dupl_points[i]
                            for j in range(i+1, len(dupl_points)):
                                v = dupl_points[j]
                                if all((u[k] - v[k]).abs() < pt_tol
                                       for k in range(len(u))):
                                    rat_points.remove(u)
                                    break

                rat_points = sorted(rat_points)
                return rat_points
            raise NotImplementedError('numerical approximation of points only for dimension 0 subschemes')

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

    Check for :trac:`11982`::

        sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
        sage: d = 7
        sage: C = Curve(x^3 + y^3 - d*z^3)
        sage: E = EllipticCurve([0,-432*d^2])
        sage: transformation = [(36*d*z-y)/(72*d),(36*d*z+y)/(72*d),x/(12*d)]
        sage: phi = E.hom(transformation, C); phi
        Scheme morphism:
          From: Elliptic Curve defined by y^2 = x^3 - 21168 over Rational Field
          To:   Projective Plane Curve over Rational Field defined by x^3 + y^3 - 7*z^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (-1/504*y + 1/2*z : 1/504*y + 1/2*z : 1/84*x)
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


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.schemes.generic.homset',
                           'SchemeHomsetModule_abelian_variety_coordinates_field',
                           SchemeHomset_points_abelian_variety_field)
