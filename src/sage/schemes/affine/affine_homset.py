r"""
Set of homomorphisms between two affine schemes

For schemes `X` and `Y`, this module implements the set of morphisms
`Hom(X,Y)`. This is done by :class:`SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a
scheme. Recall that the `K`-rational points of a scheme `X` over `k`
can be identified with the set of morphisms `Spec(K) \to X`. In Sage
the rational points are implemented by such scheme morphisms. This is
done by :class:`SchemeHomset_points` and its subclasses.

.. NOTE::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- William Stein (2006): initial version.

- Ben Hutz (2018): add numerical point support
"""

# *****************************************************************************
#        Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from copy import copy

from sage.misc.verbose import verbose
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import RationalField
from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.generic.homset import SchemeHomset_points, SchemeHomset_generic


# *******************************************************************
#  Affine varieties
# *******************************************************************

class SchemeHomset_points_spec(SchemeHomset_generic):
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
        return super()._element_constructor_(*args, **kwds)

    def _repr_(self):
        """
        Return a string representation of a homset.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_spec
            sage: S = SchemeHomset_points_spec(Spec(QQ), Spec(QQ))
            sage: S._repr_()
            'Set of rational points of Spectrum of Rational Field'
        """
        return 'Set of rational points of {}'.format(self.codomain())


class SchemeHomset_polynomial_affine_space(SchemeHomset_generic):
    """
    Set of morphisms between affine spaces defined by polynomials.

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(2, QQ)
        sage: Hom(A, A)
        Set of morphisms
          From: Affine Space of dimension 2 over Rational Field
          To:   Affine Space of dimension 2 over Rational Field
    """
    def identity(self):
        """
        The identity morphism of this homset.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(2, QQ)
            sage: I = A.identity_morphism()
            sage: I.parent()
            Set of morphisms
              From: Affine Space of dimension 2 over Rational Field
              To:   Affine Space of dimension 2 over Rational Field
            sage: _.identity() == I
            True
        """
        if self.is_endomorphism_set():
            from sage.schemes.generic.morphism import SchemeMorphism_polynomial_id
            return SchemeMorphism_polynomial_id(self.domain())
        raise TypeError("identity map is only defined for endomorphisms")


# *******************************************************************
#  Affine varieties
# *******************************************************************

class SchemeHomset_points_affine(SchemeHomset_points):
    """
    Set of rational points of an affine variety.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_affine
        sage: SchemeHomset_points_affine(Spec(QQ), AffineSpace(ZZ,2))
        Set of rational points of Affine Space of dimension 2 over Rational Field
    """

    def points(self, **kwds):
        r"""
        Return some or all rational points of an affine scheme.

        For dimension 0 subschemes points are determined through a groebner
        basis calculation. For schemes or subschemes with dimension greater than 1
        points are determined through enumeration up to the specified bound.

        Over a finite field, all points are returned. Over an infinite field, all points satisfying the bound
        are returned. For a zero-dimensional subscheme, all points are returned regardless of whether the field
        is infinite or not.

        For number fields, this uses the
        Doyle-Krumm algorithm 4 (algorithm 5 for imaginary quadratic) for
        computing algebraic numbers up to a given height [DK2013]_.

        The algorithm requires floating point arithmetic, so the user is
        allowed to specify the precision for such calculations.
        Additionally, due to floating point issues, points
        slightly larger than the bound may be returned. This can be controlled
        by lowering the tolerance.


        INPUT: keyword arguments:

        - ``bound`` -- real number (default: 0). The bound for the
          height of the coordinates. Only used for subschemes with
          dimension at least 1.

        - ``zero_tolerance`` -- positive real number (default: 10^(-10)).
          For numerically inexact fields, points are on the subscheme if they
          satisfy the equations to within tolerance.

        - ``tolerance`` -- a rational number in (0,1] used in Doyle-Krumm
          algorithm-4 for enumeration over number fields

        - ``precision`` -- the precision to use for computing the elements of
          bounded height of number fields

        OUTPUT: list of rational points of a affine scheme

        .. WARNING::

           For numerically inexact fields such as ComplexField or RealField the
           list of points returned is very likely to be incomplete. It may also
           contain repeated points due to tolerance.

        EXAMPLES: The bug reported at #11526 is fixed::

            sage: A2 = AffineSpace(ZZ, 2)
            sage: F = GF(3)
            sage: A2(F).points()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

        ::

            sage: A.<x,y> = ZZ[]
            sage: I = A.ideal(x^2 - y^2 - 1)
            sage: V = AffineSpace(ZZ, 2)
            sage: X = V.subscheme(I)
            sage: M = X(ZZ)
            sage: M.points(bound=1)
            [(-1, 0), (1, 0)]

        ::

            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)                                          # needs sage.rings.number_field
            sage: A.<x,y> = AffineSpace(K, 2)                                           # needs sage.rings.number_field
            sage: len(A(K).points(bound=2))                                             # needs sage.rings.number_field
            1849

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: E = A.subscheme([x^2 + y^2 - 1, y^2 - x^3 + x^2 + x - 1])
            sage: E(A.base_ring()).points()                                             # needs sage.libs.singular
            [(-1, 0), (0, -1), (0, 1), (1, 0)]

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)                                          # needs sage.rings.real_mpfr
            sage: E = A.subscheme([y^3 - x^3 - x^2, x*y])
            sage: E(A.base_ring()).points()                                             # needs sage.libs.singular sage.rings.real_mpfr
            verbose 0 (...: affine_homset.py, points)
            Warning: computations in the numerical fields are inexact;points
            may be computed partially or incorrectly.
            [(-1.00000000000000, 0.000000000000000),
             (0.000000000000000, 0.000000000000000)]

        ::

            sage: A.<x1,x2> = AffineSpace(CDF, 2)                                       # needs sage.rings.complex_double
            sage: E = A.subscheme([x1^2 + x2^2 + x1*x2, x1 + x2])                       # needs sage.libs.singular sage.rings.complex_double
            sage: E(A.base_ring()).points()                                             # needs sage.libs.singular sage.rings.complex_double
            verbose 0 (...: affine_homset.py, points)
            Warning: computations in the numerical fields are inexact;points
            may be computed partially or incorrectly.
            [(0.0, 0.0)]
        """
        from sage.schemes.affine.affine_space import AffineSpace_generic

        X = self.codomain()
        if not isinstance(X, AffineSpace_generic) and X.base_ring() in Fields():
            if hasattr(X.base_ring(), 'precision'):
                numerical = True
                verbose("Warning: computations in the numerical fields are inexact;points may be computed partially or incorrectly.", level=0)
                from sage.rings.real_mpfr import RR
                zero_tol = RR(kwds.pop('zero_tolerance', 10**(-10)))
                if zero_tol <= 0:
                    raise ValueError("tolerance must be positive")
            else:
                numerical = False
            # Then X must be a subscheme
            dim_ideal = X.defining_ideal().dimension()
            if dim_ideal < 0:  # no points
                return []
            if dim_ideal == 0:  # if X zero-dimensional
                rat_points = []
                AS = X.ambient_space()
                N = AS.dimension_relative()
                BR = X.base_ring()
                # need a lexicographic ordering for elimination
                R = PolynomialRing(BR, N, AS.gens(), order='lex')
                I = R.ideal(X.defining_polynomials())
                I0 = R.ideal(0)
                # Determine the points through elimination
                # This is much faster than using the I.variety() function on each affine chart.
                G = I.groebner_basis()
                if G != [1]:
                    P = {}
                    points = [P]
                    # work backwards from solving each equation for the possible
                    # values of the next coordinate
                    for i in range(len(G) - 1, -1, -1):
                        new_points = []
                        good = 0
                        for P in points:
                            # substitute in our dictionary entry that has the values
                            # of coordinates known so far. This results in a single
                            # variable polynomial (by elimination)
                            L = G[i].substitute(P)
                            if R(L).degree() > 0:
                                if numerical:
                                    for pol in L.univariate_polynomial().roots(multiplicities=False):
                                        r = L.variables()[0]
                                        varindex = R.gens().index(r)
                                        P.update({R.gen(varindex): pol})
                                        new_points.append(copy(P))
                                        good = 1
                                else:
                                    L = L.factor()
                                # the linear factors give the possible rational values of
                                # this coordinate
                                    for pol, pow in L:
                                        if pol.degree() == 1 and len(pol.variables()) == 1:
                                            good = 1
                                            r = pol.variables()[0]
                                            varindex = R.gens().index(r)
                                            # add this coordinates information to
                                            # each dictionary entry
                                            P.update({R.gen(varindex):
                                                      -pol.constant_coefficient() / pol.monomial_coefficient(r)})
                                            new_points.append(copy(P))
                            else:
                                new_points.append(P)
                                good = 1
                        if good:
                            points = new_points
                    # the dictionary entries now have values for all coordinates
                    # they are the rational solutions to the equations
                    # make them into affine points
                    for i in range(len(points)):
                        if numerical:
                            if len(points[i]) == N:
                                S = AS([points[i][R.gen(j)] for j in range(N)])
                                if all(g(list(S)) < zero_tol for g in X.defining_polynomials()):
                                    rat_points.append(S)
                        else:
                            if len(points[i]) == N and I.subs(points[i]) == I0:
                                S = X([points[i][R.gen(j)] for j in range(N)])
                                rat_points.append(S)

                rat_points = sorted(rat_points)
                return rat_points
        R = self.value_ring()
        B = kwds.pop('bound', 0)
        tol = kwds.pop('tolerance', 1e-2)
        prec = kwds.pop('precision', 53)
        if isinstance(R, RationalField) or R == ZZ:
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified" % B)
            from sage.schemes.affine.affine_rational_point import enum_affine_rational_field
            return enum_affine_rational_field(self, B)
        if R in NumberFields():
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified" % B)
            from sage.schemes.affine.affine_rational_point import enum_affine_number_field
            return enum_affine_number_field(self, bound=B, tolerance=tol, precision=prec)
        elif isinstance(R, FiniteField):
            from sage.schemes.affine.affine_rational_point import enum_affine_finite_field
            return enum_affine_finite_field(self)
        else:
            raise TypeError("unable to enumerate points over %s" % R)

    def numerical_points(self, F=None, **kwds):
        """
        Return some or all numerical approximations of rational points of an affine scheme.

        This is for dimension 0 subschemes only and the points are determined
        through a groebner calculation over the base ring and then numerically
        approximating the roots of the resulting polynomials. If the base ring
        is a number field, the embedding into ``F`` must be known.

        INPUT:

        - ``F`` -- numerical ring

        kwds:

        - ``zero_tolerance`` -- positive real number (default: 10^(-10)).
          For numerically inexact fields, points are on the subscheme if they
          satisfy the equations to within tolerance.

        OUTPUT: list of points in the ambient space

        .. WARNING::

           For numerically inexact fields the list of points returned may contain repeated
           or be missing points due to tolerance.

        EXAMPLES::

            sage: # needs sage.libs.singular sage.rings.number_field
            sage: K.<v> = QuadraticField(3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: X = A.subscheme([x^3 - v^2*y, y - v*x^2 + 3])
            sage: L = X(K).numerical_points(F=RR); L  # abs tol 1e-14
            [(-1.18738247880014, -0.558021142104134),
             (1.57693558184861, 1.30713548084184),
             (4.80659931965815, 37.0162574656220)]
            sage: L[0].codomain()
            Affine Space of dimension 2 over Real Field with 53 bits of precision

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme([y^2 - x^2 - 3*x, x^2 - 10*y])
            sage: len(X(QQ).numerical_points(F=ComplexField(100)))                      # needs sage.libs.singular
            4

        ::

            sage: A.<x1, x2> = AffineSpace(QQ, 2)
            sage: E = A.subscheme([30*x1^100 + 1000*x2^2 + 2000*x1*x2 + 1, x1 + x2])
            sage: len(E(A.base_ring()).numerical_points(F=CDF, zero_tolerance=1e-9))    # needs sage.libs.singular
            100

        TESTS::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme([y^2 - x^2 - 3*x, x^2 - 10*y])
            sage: X(QQ).numerical_points(F=QQ)
            Traceback (most recent call last):
            ...
            TypeError: F must be a numerical field

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme([y^2 - x^2 - 3*x, x^2 - 10*y])
            sage: X(QQ).numerical_points(F=CC, zero_tolerance=-1)                       # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            ValueError: tolerance must be positive
        """
        from sage.schemes.affine.affine_space import AffineSpace_generic
        if F is None:
            from sage.rings.cc import CC as F
        if F not in Fields() or not hasattr(F, 'precision'):
            raise TypeError('F must be a numerical field')
        X = self.codomain()
        if X.base_ring() not in NumberFields():
            raise TypeError('base ring must be a number field')

        AA = X.ambient_space().change_ring(F)
        if not isinstance(X, AffineSpace_generic) and X.base_ring() in Fields():
            # Then X must be a subscheme
            dim_ideal = X.defining_ideal().dimension()
            if dim_ideal != 0:  # no points
                return []
        else:
            return []

        # if X zero-dimensional
        from sage.rings.real_mpfr import RR
        zero_tol = RR(kwds.pop('zero_tolerance', 10**(-10)))
        if zero_tol <= 0:
            raise ValueError("tolerance must be positive")
        rat_points = []
        PS = X.ambient_space()
        N = PS.dimension_relative()
        BR = X.base_ring()
        # need a lexicographic ordering for elimination
        R = PolynomialRing(BR, N, PS.gens(), order='lex')
        RF = R.change_ring(F)
        I = R.ideal(X.defining_polynomials())
        # Determine the points through elimination This is much faster
        # than using the I.variety() function on each affine chart.
        G = I.groebner_basis()
        G = [RF(g) for g in G]
        if G != [1]:
            P = {}
            points = [P]
            # work backwards from solving each equation for the possible
            # values of the next coordinate
            for g in reversed(G):
                new_points = []
                good = False
                for P in points:
                    # substitute in our dictionary entry that has the
                    # values of coordinates known so far. This results
                    # in a single variable polynomial (by elimination)
                    L = g.substitute(P)
                    if len(RF(L).variables()) == 1:
                        r = L.variables()[0]
                        var = RF.gen(RF.gens().index(r))

                        for pol in L.univariate_polynomial().roots(ring=F,
                                multiplicities=False):
                            P[var] = pol
                            new_points.append(copy(P))
                            good = True
                    else:
                        new_points.append(P)
                        good = True
                if good:
                    points = new_points
            # the dictionary entries now have values for all
            # coordinates they are the rational solutions to the
            # equations make them into affine points
            polys = [g.change_ring(F) for g in X.defining_polynomials()]
            for P in points:
                if len(P) == N:
                    S = AA([P[R.gen(j)] for j in range(N)])
                    if all(g(list(S)) < zero_tol for g in polys):
                        rat_points.append(S)

        rat_points = sorted(rat_points)
        return rat_points
