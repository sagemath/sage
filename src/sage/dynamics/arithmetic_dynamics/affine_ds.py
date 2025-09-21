r"""
Dynamical systems on affine schemes

An endomorphism of an affine scheme or subscheme determined by polynomials
or rational functions.

The main constructor function is given by ``DynamicalSystem_affine``.
The constructor function can take polynomials, rational functions, or
morphisms from which to construct a dynamical system. If the domain is
not specified, it is constructed. However, if you plan on working with
points or subvarieties in the domain, it recommended to specify the
domain.

The initialization checks are always performed by the constructor functions.
It is possible, but not recommended, to skip these checks by calling the
class initialization directly.

AUTHORS:

- David Kohel, William Stein

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (2017) relocate code and create new class
"""

# ****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.fields import Fields
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sage.matrix.constructor import identity_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import typecall
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.fraction_field import FractionField
from sage.rings.fraction_field import FractionField_generic
from sage.rings.quotient_ring import QuotientRing_nc
from sage.schemes.affine.affine_morphism import SchemeMorphism_polynomial_affine_space
from sage.schemes.affine.affine_morphism import SchemeMorphism_polynomial_affine_space_field
from sage.schemes.affine.affine_morphism import SchemeMorphism_polynomial_affine_space_finite_field
from sage.schemes.affine.affine_space import AffineSpace_generic
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.structure.element import get_coercion_model

import sage.rings.abc


class DynamicalSystem_affine(SchemeMorphism_polynomial_affine_space,
                             DynamicalSystem):
    r"""
    An endomorphism of affine schemes determined by rational functions.

    .. WARNING::

        You should not create objects of this class directly because
        no type or consistency checking is performed. The preferred
        method to construct such dynamical systems is to use
        :func:`~sage.dynamics.arithmetic_dynamics.generic_ds.DynamicalSystem_affine`
        function.

    INPUT:

    - ``morphism_or_polys`` -- a SchemeMorphism, a polynomial, a
      rational function, or a list or tuple of polynomials or rational
      functions

    - ``domain`` -- (optional) affine space or subscheme of such;
      the following combinations of ``morphism_or_polys`` and
      ``domain`` are meaningful:

      * ``morphism_or_polys`` is a SchemeMorphism; ``domain`` is
        ignored in this case

      * ``morphism_or_polys`` is a list of polynomials or rational
        functions that define a rational endomorphism of ``domain``

      * ``morphism_or_polys`` is a list of polynomials or rational
        functions and ``domain`` is unspecified; ``domain`` is then
        taken to be the affine space of appropriate dimension over the
        common base ring, if one exists, of the elements of ``morphism_or_polys``

      * ``morphism_or_polys`` is a single polynomial or rational
        function; ``domain`` is ignored and assumed to be the
        1-dimensional affine space over the base ring of
        ``morphism_or_polys``

    OUTPUT: :class:`DynamicalSystem_affine`

    EXAMPLES::

        sage: A3.<x,y,z> = AffineSpace(QQ, 3)
        sage: DynamicalSystem_affine([x, y, 1])
        Dynamical System of Affine Space of dimension 3 over Rational Field
              Defn: Defined on coordinates by sending (x, y, z) to
                    (x, y, 1)

    ::

        sage: R.<x,y> = QQ[]
        sage: DynamicalSystem_affine([x/y, y^2 + 1])
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (x/y, y^2 + 1)

    ::

        sage: R.<t> = ZZ[]
        sage: DynamicalSystem_affine(t^2 - 1)
        Dynamical System of Affine Space of dimension 1 over Integer Ring
          Defn: Defined on coordinates by sending (t) to
                (t^2 - 1)

    ::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: X = A.subscheme([x-y^2])
        sage: DynamicalSystem_affine([9/4*x^2, 3/2*y], domain=X)
        Dynamical System of Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -y^2 + x
              Defn: Defined on coordinates by sending (x, y) to
                    (9/4*x^2, 3/2*y)

    ::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: H = End(A)
        sage: f = H([x^2, y^2])
        sage: DynamicalSystem_affine(f)
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (x^2, y^2)

    Notice that `ZZ` becomes `QQ` since the function is rational::

        sage: A.<x,y> = AffineSpace(ZZ, 2)
        sage: DynamicalSystem_affine([3*x^2/(5*y), y^2/(2*x^2)])
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (3*x^2/(5*y), y^2/(2*x^2))

    ::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: DynamicalSystem_affine([3/2*x^2, y^2])
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (3/2*x^2, y^2)

    If you pass in quotient ring elements, they are reduced::

        sage: # needs sage.libs.singular
        sage: A.<x,y,z> = AffineSpace(QQ, 3)
        sage: X = A.subscheme([x - y])
        sage: u,v,w = X.coordinate_ring().gens()
        sage: DynamicalSystem_affine([u, v, u+v], domain=X)
        Dynamical System of Closed subscheme of Affine Space of dimension 3
        over Rational Field defined by:
          x - y
          Defn: Defined on coordinates by sending (x, y, z) to
                (y, y, 2*y)

    ::

        sage: # needs sage.libs.singular
        sage: R.<t> = PolynomialRing(QQ)
        sage: A.<x,y,z> = AffineSpace(R, 3)
        sage: X = A.subscheme(x^2 - y^2)
        sage: H = End(X)
        sage: f = H([x^2/(t*y), t*y^2, x*z])
        sage: DynamicalSystem_affine(f)
        Dynamical System of Closed subscheme of Affine Space of dimension 3
        over Univariate Polynomial Ring in t over Rational Field defined by:
          x^2 - y^2
          Defn: Defined on coordinates by sending (x, y, z) to
                (x^2/(t*y), t*y^2, x*z)

    ::

        sage: x = var('x')                                                              # needs sage.symbolic
        sage: DynamicalSystem_affine(x^2 + 1)                                           # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: symbolic ring cannot be the base ring
    """

    @staticmethod
    def __classcall_private__(cls, morphism_or_polys, domain=None):
        r"""
        Return the appropriate dynamical system on an affine scheme.

        TESTS::

            sage: A.<x> = AffineSpace(ZZ,1)
            sage: A1.<z> = AffineSpace(CC,1)
            sage: H = End(A1)
            sage: f2 = H([z^2+1])
            sage: f = DynamicalSystem_affine(f2, A)
            sage: f.domain() is A
            False

        ::

            sage: P1.<x,y> = ProjectiveSpace(QQ,1)
            sage: DynamicalSystem_affine([y, 2*x], domain=P1)
            Traceback (most recent call last):
            ...
            ValueError: "domain" must be an affine scheme
            sage: H = End(P1)
            sage: DynamicalSystem_affine(H([y, 2*x]))
            Traceback (most recent call last):
            ...
            ValueError: "domain" must be an affine scheme

        ::

            sage: R.<x,y,z> = QQ[]
            sage: f = DynamicalSystem_affine([x+y+z, y*z])
            Traceback (most recent call last):
            ...
            ValueError: number of polys does not match dimension of Affine Space of dimension 3 over Rational Field

        ::

            sage: A.<x,y> = AffineSpace(QQ,2)
            sage: f = DynamicalSystem_affine([CC.0*x^2, y^2], domain=A)
            Traceback (most recent call last):
            ...
            TypeError: coefficients of polynomial not in Rational Field
        """
        if isinstance(morphism_or_polys, SchemeMorphism_polynomial):
            morphism = morphism_or_polys
            R = morphism.base_ring()
            polys = list(morphism)
            domain = morphism.domain()
            if not isinstance(domain, AffineSpace_generic) and not isinstance(domain, AlgebraicScheme_subscheme_affine):
                raise ValueError('"domain" must be an affine scheme')
            if domain != morphism_or_polys.codomain():
                raise ValueError('domain and codomain do not agree')
            if R not in Fields():
                return typecall(cls, polys, domain)
            if isinstance(R, FiniteField):
                return DynamicalSystem_affine_finite_field(polys, domain)
            return DynamicalSystem_affine_field(polys, domain)
        elif isinstance(morphism_or_polys, (list, tuple)):
            polys = list(morphism_or_polys)
        else:
            polys = [morphism_or_polys]

        PR = get_coercion_model().common_parent(*polys)
        fraction_field = any(isinstance(poly.parent(), FractionField_generic) for poly in polys)
        if fraction_field:
            K = PR.base_ring().fraction_field()
            # Replace base ring with its fraction field
            PR = PR.ring().change_ring(K).fraction_field()
            polys = [PR(poly) for poly in polys]
        else:
            quotient_ring = any(isinstance(poly.parent(), QuotientRing_nc) for poly in polys)
            # If any of the list entries lies in a quotient ring, we try
            # to lift all entries to a common polynomial ring.
            if quotient_ring:
                polys = [PR(poly).lift() for poly in polys]
            else:
                polys = [PR(poly) for poly in polys]
        if domain is None:
            if isinstance(PR, sage.rings.abc.SymbolicRing):
                raise TypeError("symbolic ring cannot be the base ring")
            if fraction_field:
                PR = PR.ring()
            domain = AffineSpace(PR)
        else:
            # Check if we can coerce the given polynomials over the given domain
            PR = domain.ambient_space().coordinate_ring()
            try:
                if fraction_field:
                    PR = PR.fraction_field()
                polys = [PR(poly) for poly in polys]
            except TypeError:
                raise TypeError('coefficients of polynomial not in {}'.format(domain.base_ring()))
        if len(polys) != domain.ambient_space().coordinate_ring().ngens():
            raise ValueError(f'number of polys does not match dimension of {domain}')
        R = domain.base_ring()
        if isinstance(R, sage.rings.abc.SymbolicRing):
            raise TypeError("symbolic ring cannot be the base ring")
        if not isinstance(domain, AffineSpace_generic) and not isinstance(domain, AlgebraicScheme_subscheme_affine):
            raise ValueError('"domain" must be an affine scheme')

        if R not in Fields():
            return typecall(cls, polys, domain)
        if isinstance(R, FiniteField):
            return DynamicalSystem_affine_finite_field(polys, domain)
        return DynamicalSystem_affine_field(polys, domain)

    def __init__(self, polys_or_rat_fncts, domain):
        r"""
        The Python constructor.

        See :class:`DynamicalSystem` for details.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: DynamicalSystem_affine([3/5*x^2, y^2/(2*x^2)], domain=A)
            Dynamical System of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (3/5*x^2, y^2/(2*x^2))
        """
        L = polys_or_rat_fncts
        # Next attribute needed for _fast_eval and _fastpolys
        R = L[0].base_ring()
        self._is_prime_finite_field = isinstance(R, FiniteField) and R.degree() == 1
        DynamicalSystem.__init__(self, L, domain)

    def __copy__(self):
        r"""
        Return a copy of this dynamical system.

        OUTPUT: :class:`DynamicalSystem_affine`

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ,2)
            sage: f = DynamicalSystem_affine([3/5*x^2,6*y^2])
            sage: g = copy(f)
            sage: f == g
            True
            sage: f is g
            False
        """
        return DynamicalSystem_affine(self._polys, self.domain())

    def homogenize(self, n):
        r"""
        Return the homogenization of this dynamical system.

        If its domain is a subscheme, the domain of the homogenized map
        is the projective embedding of the domain. The domain and codomain
        can be homogenized at different coordinates: ``n[0]`` for the
        domain and ``n[1]`` for the codomain.

        INPUT:

        - ``n`` -- tuple of nonnegative integers; if `n` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: :class:`DynamicalSystem_projective`

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: f = DynamicalSystem_affine([(x^2-2)/x^5, y^2])
            sage: f.homogenize(2)
            Dynamical System of Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (x0^2*x2^5 - 2*x2^7 : x0^5*x1^2 : x0^5*x2^2)

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)
            sage: f = DynamicalSystem_affine([(x^2-2)/(x*y), y^2-x])
            sage: f.homogenize((2, 0))
            Dynamical System of Projective Space of dimension 2 over Complex Field with 53 bits of precision
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (x0*x1*x2^2 : x0^2*x2^2 + (-2.00000000000000)*x2^4 : x0*x1^3 - x0^2*x1*x2)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: f = DynamicalSystem_affine([9*y^2, 3*y], domain=X)
            sage: f.homogenize(2)
            Dynamical System of Closed subscheme of Projective Space
            of dimension 2 over Integer Ring defined by:
                x1^2 - x0*x2
                Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                      (9*x1^2 : 3*x1*x2 : x2^2)

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([x^2-1])
            sage: f.homogenize((1, 0))
            Dynamical System of Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x0 : x1) to
                    (x1^2 : x0^2 - x1^2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<a> = PolynomialRing(QQbar)
            sage: A.<x,y> = AffineSpace(R, 2)
            sage: f = DynamicalSystem_affine([QQbar(sqrt(2))*x*y, a*x^2])               # needs sage.symbolic
            sage: f.homogenize(2)                                                       # needs sage.symbolic
            Dynamical System of Projective Space of dimension 2 over Univariate
            Polynomial Ring in a over Algebraic Field
              Defn: Defined on coordinates by sending (x0 : x1 : x2) to
                    (1.414213562373095?*x0*x1 : a*x0^2 : x2^2)

        ::

            sage: P.<x,y,z> = AffineSpace(QQ, 3)
            sage: f = DynamicalSystem_affine([x^2 - 2*x*y + z*x, z^2 -y^2 , 5*z*y])
            sage: f.homogenize(2).dehomogenize(2) == f
            True
        """
        F = self.as_scheme_morphism().homogenize(n)
        return F.as_dynamical_system()

    def dynatomic_polynomial(self, period):
        r"""
        Compute the (affine) dynatomic polynomial of a dynamical system
        `f: \mathbb{A}^1 \to \mathbb{A}^1`.

        The dynatomic polynomial is the analog of the cyclotomic polynomial
        and its roots are the points of formal period `n`.

        ALGORITHM:

        Homogenize to a map `f: \mathbb{P}^1 \to \mathbb{P}^1` and compute
        the dynatomic polynomial there. Then, dehomogenize.

        INPUT:

        - ``period`` -- positive integer or a list/tuple `[m,n]`,
          where `m` is the preperiod and `n` is the period

        OUTPUT:

        If possible, a single variable polynomial in the coordinate ring
        of the polynomial. Otherwise a fraction field element of the
        coordinate ring of the polynomial.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([x^2+y^2, y^2])
            sage: f.dynatomic_polynomial(2)
            Traceback (most recent call last):
            ...
            TypeError: does not make sense in dimension >1

        ::

            sage: A.<x> = AffineSpace(ZZ, 1)
            sage: f = DynamicalSystem_affine([(x^2+1)/x])
            sage: f.dynatomic_polynomial(4)                                             # needs sage.libs.pari
            2*x^12 + 18*x^10 + 57*x^8 + 79*x^6 + 48*x^4 + 12*x^2 + 1

        ::

            sage: A.<x> = AffineSpace(CC, 1)
            sage: f = DynamicalSystem_affine([(x^2+1)/(3*x)])
            sage: f.dynatomic_polynomial(3)
            13.0000000000000*x^6 + 117.000000000000*x^4 + 78.0000000000000*x^2 +
            1.00000000000000

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([x^2-10/9])
            sage: f.dynatomic_polynomial([2, 1])
            531441*x^4 - 649539*x^2 - 524880

        ::

            sage: A.<x> = AffineSpace(CC, 1)
            sage: f = DynamicalSystem_affine([x^2+CC.0])
            sage: f.dynatomic_polynomial(2)
            x^2 + x + 1.00000000000000 + 1.00000000000000*I

        ::

            sage: K.<c> = FunctionField(QQ)
            sage: A.<x> = AffineSpace(K, 1)
            sage: f = DynamicalSystem_affine([x^2 + c])
            sage: f.dynatomic_polynomial(4)
            x^12 + 6*c*x^10 + x^9 + (15*c^2 + 3*c)*x^8 + 4*c*x^7 + (20*c^3 + 12*c^2 + 1)*x^6
            + (6*c^2 + 2*c)*x^5 + (15*c^4 + 18*c^3 + 3*c^2 + 4*c)*x^4 + (4*c^3 + 4*c^2 + 1)*x^3
            + (6*c^5 + 12*c^4 + 6*c^3 + 5*c^2 + c)*x^2 + (c^4 + 2*c^3 + c^2 + 2*c)*x
            + c^6 + 3*c^5 + 3*c^4 + 3*c^3 + 2*c^2 + 1

        ::

            sage: A.<z> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([z^2+3/z+1/7])
            sage: f.dynatomic_polynomial(1).parent()
            Multivariate Polynomial Ring in z over Rational Field

        ::

            sage: R.<c> = QQ[]
            sage: A.<z> = AffineSpace(R,1)
            sage: f = DynamicalSystem_affine([z^2 + c])
            sage: f.dynatomic_polynomial([1,1])
            z^2 + z + c

        ::

            sage: A.<x> = AffineSpace(CC,1)
            sage: F = DynamicalSystem_affine([1/2*x^2 + CC(sqrt(3))])                   # needs sage.symbolic
            sage: F.dynatomic_polynomial([1,1])                                         # needs sage.symbolic
            (0.125000000000000*x^4 + 0.366025403784439*x^2 + 1.50000000000000)/(0.500000000000000*x^2 - x + 1.73205080756888)

        TESTS::

            sage: R.<c> = QQ[]
            sage: Pc.<x,y> = ProjectiveSpace(R, 1)
            sage: G = DynamicalSystem_projective([(1/2*c + 1/2)*x^2 + (-2*c)*x*y + 2*c*y^2 , \
            ....: (1/4*c + 1/2)*x^2 + (-c - 1)*x*y + (c + 1)*y^2])
            sage: G.dehomogenize(1).dynatomic_polynomial(2)
            (1/4*c + 1/4)*x^2 + (-c - 1/2)*x + c + 1
        """
        from sage.schemes.affine.affine_space import AffineSpace_generic
        if not isinstance(self.domain(), AffineSpace_generic):
            raise NotImplementedError("not implemented for subschemes")
        if self.domain().dimension_relative() > 1:
            raise TypeError("does not make sense in dimension >1")
        G = self.homogenize(1)
        F = G.dynatomic_polynomial(period)
        T = G.domain().coordinate_ring()
        S = self.domain().coordinate_ring()
        if isinstance(F.parent(), sage.rings.abc.SymbolicRing):
            from sage.symbolic.ring import var
            u = var(self.domain().coordinate_ring().variable_name())
            return F.subs({F.variables()[0]: u, F.variables()[1]: 1})
        elif T(F.denominator()).degree() == 0:
            R = F.parent()
            phi = R.hom([S.gen(0), 1], S)
            return phi(F)
        else:
            R = F.numerator().parent()
            phi = R.hom([S.gen(0), 1], S)
            return phi(F.numerator()) / phi(F.denominator())

    def nth_iterate_map(self, n):
        r"""
        Return the `n`-th iterate of ``self``.

        ALGORITHM:

        Uses a form of successive squaring to reducing computations.

        .. TODO::

            This could be improved.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: a dynamical system of affine space

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: f = DynamicalSystem_affine([(x^2-2)/(2*y), y^2-3*x])
            sage: f.nth_iterate_map(2)
            Dynamical System of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    ((x^4 - 4*x^2 - 8*y^2 + 4)/(8*y^4 - 24*x*y^2), (2*y^5 - 12*x*y^3
            + 18*x^2*y - 3*x^2 + 6)/(2*y))

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([(3*x^2-2)/(x)])
            sage: f.nth_iterate_map(3)
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x) to
                    ((2187*x^8 - 6174*x^6 + 6300*x^4 - 2744*x^2 + 432)/(81*x^7 -
            168*x^5 + 112*x^3 - 24*x))

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: f = DynamicalSystem_affine([9*x^2, 3*y], domain=X)
            sage: f.nth_iterate_map(2)
            Dynamical System of Closed subscheme of Affine Space of dimension 2
            over Integer Ring defined by:
              -y^2 + x
              Defn: Defined on coordinates by sending (x, y) to
                    (729*x^4, 9*y)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: f = DynamicalSystem_affine([3/5*x^2, y^2/(2*x^2)])
            sage: f.nth_iterate_map(2)
            Dynamical System of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (27/125*x^4, y^4/(72/25*x^8))
        """
        F = list(self._polys)
        R = F[0].parent()
        Coord_ring = self.codomain().ambient_space().coordinate_ring()
        D = Integer(n).digits(2)
        PHI = list(Coord_ring.gens())
        for i in range(len(D)):
            for k in range(D[i]):
                PHI = [poly(F) for poly in PHI]
            if i != len(D) - 1:  # avoid extra iterate
                F = [R(poly(F)) for poly in F]  # 'square'
        return DynamicalSystem_affine(PHI, domain=self.domain())

    def nth_iterate(self, P, n):
        r"""
        Return the `n`-th iterate of the point `P` by this dynamical system.

        INPUT:

        - ``P`` -- a point in the map's domain

        - ``n`` -- positive integer

        OUTPUT: a point in the map's codomain

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([(x-2*y^2)/x, 3*x*y])
            sage: f.nth_iterate(A(9, 3), 3)
            (-104975/13123, -9566667)

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: f = DynamicalSystem_affine([9*y^2, 3*y], domain=X)
            sage: f.nth_iterate(X(9, 3), 4)
            (59049, 243)

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: A.<x,y> = AffineSpace(FractionField(R), 2)
            sage: f = DynamicalSystem_affine([(x-t*y^2)/x, t*x*y])
            sage: f.nth_iterate(A(1, t), 3)
            ((-t^16 + 3*t^13 - 3*t^10 + t^7 + t^5 + t^3 - 1)/(t^5 + t^3 - 1), -t^9 - t^7 + t^4)

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([x^2-y^2])
            sage: f = DynamicalSystem_affine([x^2, y^2, x+y], domain=X)
            sage: f.nth_iterate_map(2)
            Dynamical System of Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
                  x^2 - y^2
                  Defn: Defined on coordinates by sending (x, y, z) to
                        (x^4, y^4, x^2 + y^2)
        """
        n = int(n)
        if n == 0:
            return P
        Q = P
        for i in range(n):
            Q = self(Q)
        return Q

    def orbit(self, P, n):
        r"""
        Return the orbit of `P` by the dynamical system.

        Let `F` be this dynamical system. If `n` is an integer
        return `[P, F(P), \ldots, F^n(P)]`. If `n` is a list or
        tuple `n = [m,k]` return `[F^m(P), \ldots, F^k(P)]`.

        INPUT:

        - ``P`` -- a point in the map's domain

        - ``n`` -- nonnegative integer or list or tuple of
          two nonnegative integers

        OUTPUT: list of points in the map's codomain

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([(x-2*y^2)/x, 3*x*y])
            sage: f.orbit(A(9, 3), 3)
            [(9, 3), (-1, 81), (13123, -243), (-104975/13123, -9566667)]

        ::

            sage: A.<x> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([(x-2)/x])
            sage: f.orbit(A(1/2), [1, 3])
            [(-3), (5/3), (-1/5)]

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: X = A.subscheme([x-y^2])
            sage: f = DynamicalSystem_affine([9*y^2, 3*y], domain=X)
            sage: f.orbit(X(9, 3), (0, 4))
            [(9, 3), (81, 9), (729, 27), (6561, 81), (59049, 243)]

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: A.<x,y> = AffineSpace(FractionField(R), 2)
            sage: f = DynamicalSystem_affine([(x-t*y^2)/x, t*x*y])
            sage: f.orbit(A(1, t), 3)
            [(1, t),
             (-t^3 + 1, t^2),
             ((t^5 + t^3 - 1)/(t^3 - 1), -t^6 + t^3),
             ((-t^16 + 3*t^13 - 3*t^10 + t^7 + t^5 + t^3 - 1)/(t^5 + t^3 - 1), -t^9 - t^7 + t^4)]
        """
        Q = P
        if isinstance(n, (list, tuple)):
            bounds = list(n)
        else:
            bounds = [0, n]
        for i in range(1, bounds[0] + 1):
            Q = self(Q)
        orb = [Q]
        for i in range(bounds[0] + 1, bounds[1] + 1):
            Q = self(Q)
            orb.append(Q)
        return orb

    def multiplier(self, P, n, check=True):
        r"""
        Return the multiplier of the point `P` of period `n` by this
        dynamical system.

        INPUT:

        - ``P`` -- a point on domain of the map

        - ``n`` -- positive integer, the period of `P`

        - ``check`` -- boolean (default: ``True``); verify that `P`
          has period `n`

        OUTPUT: a square matrix of size ``self.codomain().dimension_relative()``
        in the ``base_ring`` of the map

        EXAMPLES::

            sage: P.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([x^2, y^2])
            sage: f.multiplier(P([1, 1]), 1)
            [2 0]
            [0 2]

        ::

            sage: P.<x,y,z> = AffineSpace(QQ, 3)
            sage: f = DynamicalSystem_affine([x, y^2, z^2 - y])
            sage: f.multiplier(P([1/2, 1, 0]), 2)
            [1 0 0]
            [0 4 0]
            [0 0 0]

        ::

            sage: P.<x> = AffineSpace(CC, 1)
            sage: f = DynamicalSystem_affine([x^2 + 1/2])
            sage: f.multiplier(P([0.5 + 0.5*I]), 1)                                     # needs sage.symbolic
            [1.00000000000000 + 1.00000000000000*I]

        ::

            sage: R.<t> = PolynomialRing(CC, 1)
            sage: P.<x> = AffineSpace(R, 1)
            sage: f = DynamicalSystem_affine([x^2 - t^2 + t])
            sage: f.multiplier(P([-t + 1]), 1)
            [(-2.00000000000000)*t + 2.00000000000000]

        ::

            sage: P.<x,y> = AffineSpace(QQ, 2)
            sage: X = P.subscheme([x^2 - y^2])
            sage: f = DynamicalSystem_affine([x^2, y^2], domain=X)
            sage: f.multiplier(X([1, 1]), 1)
            [2 0]
            [0 2]
        """
        if check:
            if self.nth_iterate(P, n) != P:
                raise ValueError("%s is not periodic of period %s" % (P, n))
            if n < 1:
                raise ValueError("period must be a positive integer")
        N = self.domain().ambient_space().dimension_relative()
        l = identity_matrix(FractionField(self.codomain().base_ring()), N, N)
        Q = P
        J = self.jacobian()
        for i in range(n):
            R = self(Q)
            l = J(tuple(Q)) * l  # chain rule matrix multiplication
            Q = R
        return l

    def conjugate(self, M):
        r"""
        Conjugate this dynamical system by ``M``, i.e. `M^{-1} \circ f \circ M`.

        If possible the new map will be defined over the same space.
        Otherwise, will try to coerce to the base ring of ``M``.

        INPUT:

        - ``M`` -- a square invertible matrix

        OUTPUT: an affine dynamical system

        EXAMPLES::

            sage: A.<t> = AffineSpace(QQ, 1)
            sage: f = DynamicalSystem_affine([t^2 + 1])
            sage: f.conjugate(matrix([[1,2], [0,1]]))
            Dynamical System of Affine Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (t) to
                    (t^2 + 4*t + 3)

        ::

            sage: A.<x,y> = AffineSpace(ZZ,2)
            sage: f = DynamicalSystem_affine([x^3 + y^3, y^2])
            sage: f.conjugate(matrix([[1,2,3], [0,1,2], [0,0,1]]))
            Dynamical System of Affine Space of dimension 2 over Integer Ring
              Defn: Defined on coordinates by sending (x, y) to
                    (x^3 + 6*x^2*y + 12*x*y^2 + 9*y^3
                      + 9*x^2 + 36*x*y + 40*y^2 + 27*x + 58*y + 28, y^2 + 4*y + 2)

        ::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<i> = NumberField(x^2 + 1)
            sage: A.<x> = AffineSpace(ZZ,1)
            sage: f = DynamicalSystem_affine([x^3 + 2*x^2 + 3])
            sage: f.conjugate(matrix([[i,i], [0,-i]]))
            Dynamical System of Affine Space of dimension 1 over Integer Ring
              Defn: Defined on coordinates by sending (x) to
                    (x^3 + x^2 - x - 5)
        """
        d = self.codomain().ngens()
        f = self.homogenize(d).conjugate(M)
        return f.dehomogenize(d)

    def degree(self):
        r"""
        Return the degree of the affine dynamical system.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<c> = QuadraticField(7)
            sage: A.<x,y,z> = AffineSpace(R, 3)
            sage: f = DynamicalSystem_affine([x^2 + y^5 + c, x^11, z^19])
            sage: f.degree()
            19

        ::

            sage: R.<c> = QQ[]
            sage: A.<x> = AffineSpace(R, 1)
            sage: f = DynamicalSystem_affine([x^4])
            sage: f.degree()
            4

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([x, y/(1 + x^2)])
            sage: f.degree()
            2
        """
        return self.as_scheme_morphism().degree()


class DynamicalSystem_affine_field(DynamicalSystem_affine,
                                   SchemeMorphism_polynomial_affine_space_field):
    @cached_method
    def weil_restriction(self):
        r"""
        Compute the Weil restriction of this morphism over some extension field.

        If the field is a finite field, then this computes
        the Weil restriction to the prime subfield.

        A Weil restriction of scalars - denoted `Res_{L/k}` - is a
        functor which, for any finite extension of fields `L/k` and
        any algebraic variety `X` over `L`, produces another
        corresponding variety `Res_{L/k}(X)`, defined over `k`. It is
        useful for reducing questions about varieties over large
        fields to questions about more complicated varieties over
        smaller fields. Since it is a functor it also applied to morphisms.
        In particular, the functor applied to a morphism gives the equivalent
        morphism from the Weil restriction of the domain to the Weil restriction
        of the codomain.

        OUTPUT:

        Scheme morphism on the Weil restrictions of the domain
        and codomain of the map.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K.<v> = QuadraticField(5)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: f = DynamicalSystem_affine([x^2 - y^2, y^2])
            sage: f.weil_restriction()
            Dynamical System of Affine Space of dimension 4 over Rational Field
              Defn: Defined on coordinates by sending (z0, z1, z2, z3) to
                    (z0^2 + 5*z1^2 - z2^2 - 5*z3^2, 2*z0*z1 - 2*z2*z3, z2^2 + 5*z3^2, 2*z2*z3)

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = QuadraticField(5)
            sage: PS.<x,y> = AffineSpace(K, 2)
            sage: f = DynamicalSystem_affine([x, y])
            sage: F = f.weil_restriction()
            sage: P = PS(2, 1)
            sage: Q = P.weil_restriction()
            sage: f(P).weil_restriction() == F(Q)
            True
        """
        F = self.as_scheme_morphism().weil_restriction()
        return F.as_dynamical_system()

    def reduce_base_field(self):
        """
        Return this map defined over the field of definition of the coefficients.

        The base field of the map could be strictly larger than
        the field where all of the coefficients are defined. This function
        reduces the base field to the minimal possible. This can be done when
        the base ring is a number field, ``QQbar``, a finite field, or algebraic
        closure of a finite field.

        OUTPUT: a dynamical system

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<t> = GF(5^2)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: f = DynamicalSystem_affine([x^2 + 3*y^2, 3*y^2])
            sage: f.reduce_base_field()
            Dynamical System of Affine Space of dimension 2 over Finite Field of size 5
              Defn: Defined on coordinates by sending (x, y) to
                    (x^2 - 2*y^2, -2*y^2)

        ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: A.<x,y> = AffineSpace(QQbar, 2)
            sage: f = DynamicalSystem_affine([x^2 + QQbar(sqrt(3))*y^2,
            ....:                             QQbar(sqrt(-1))*y^2])
            sage: f.reduce_base_field()
            Dynamical System of Affine Space of dimension 2 over Number Field in a
             with defining polynomial y^4 - y^2 + 1
             with a = -0.866025403784439? + 0.50000000000000000?*I
              Defn: Defined on coordinates by sending (x, y) to
                    (x^2 + (a^3 - 2*a)*y^2, (a^3)*y^2)

        ::

            sage: # needs sage.rings.number_field
            sage: K.<v> = CyclotomicField(5)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: f = DynamicalSystem_affine([(3*x^2 + y) / (5*x), (y^2+1) / (x+y)])
            sage: f.reduce_base_field()
            Dynamical System of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    ((3*x^2 + y)/(5*x), (5*y^2 + 5)/(5*x + 5*y))
        """
        return self.as_scheme_morphism().reduce_base_field().as_dynamical_system()


class DynamicalSystem_affine_finite_field(DynamicalSystem_affine_field,
                                    SchemeMorphism_polynomial_affine_space_finite_field):

    def orbit_structure(self, P):
        r"""
        Every point is preperiodic over a finite field.

        This function returns the pair `[m,n]` where `m` is the
        preperiod and `n` is the period of the point `P` by this map.

        INPUT:

        - ``P`` -- a point in the map's domain

        OUTPUT: list `[m, n]` of integers

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(13), 2)
            sage: f = DynamicalSystem_affine([x^2 - 1, y^2])
            sage: f.orbit_structure(A(2, 3))
            [1, 6]

        ::

            sage: # needs sage.rings.finite_rings
            sage: A.<x,y,z> = AffineSpace(GF(49, 't'), 3)
            sage: f = DynamicalSystem_affine([x^2 - z, x - y + z, y^2 - x^2])
            sage: f.orbit_structure(A(1, 1, 2))
            [7, 6]
        """
        orbit = []
        index = 1
        Q = P
        while Q not in orbit:
            orbit.append(Q)
            Q = self(Q)
            index += 1
        I = orbit.index(Q)
        return [I, index - I - 1]

    def cyclegraph(self):
        r"""
        Return the digraph of all orbits of this morphism mod `p`.

        For subschemes, only points on the subscheme whose
        image are also on the subscheme are in the digraph.

        OUTPUT: a digraph

        EXAMPLES::

            sage: P.<x,y> = AffineSpace(GF(5), 2)
            sage: f = DynamicalSystem_affine([x^2 - y, x*y + 1])
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 25 vertices

        ::

            sage: # needs sage.rings.finite_rings
            sage: P.<x> = AffineSpace(GF(3^3, 't'), 1)
            sage: f = DynamicalSystem_affine([x^2 - 1])
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 27 vertices

        ::

            sage: P.<x,y> = AffineSpace(GF(7), 2)
            sage: X = P.subscheme(x - y)
            sage: f = DynamicalSystem_affine([x^2, y^2], domain=X)
            sage: f.cyclegraph()                                                        # needs sage.graphs
            Looped digraph on 7 vertices
        """
        V = []
        E = []
        from sage.schemes.affine.affine_space import AffineSpace_generic
        if isinstance(self.domain(), AffineSpace_generic):
            for P in self.domain():
                V.append(str(P))
                Q = self(P)
                E.append([str(Q)])
        else:
            X = self.domain()
            for P in X.ambient_space():
                try:
                    XP = X.point(P)
                    V.append(str(XP))
                    Q = self(XP)
                    E.append([str(Q)])
                except TypeError:  # not on the scheme
                    pass
        from sage.graphs.digraph import DiGraph
        return DiGraph(dict(zip(V, E)), loops=True)
