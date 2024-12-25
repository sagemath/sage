r"""
Subschemes of projective space

AUTHORS:

- David Kohel (2005): initial version.
- William Stein (2005): initial version.
- Volker Braun (2010-12-24): documentation of schemes and
  refactoring. Added coordinate neighborhoods and is_smooth()
- Ben Hutz (2013) refactoring
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Ben Hutz <bn4941@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import binomial
from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
from sage.schemes.projective.projective_morphism import SchemeMorphism_polynomial_projective_subscheme_field

lazy_import('sage.dynamics.arithmetic_dynamics.generic_ds', 'DynamicalSystem')
lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.schemes.elliptic_curves.ell_generic', 'EllipticCurve_generic', as_='EllipticCurve')


class AlgebraicScheme_subscheme_projective(AlgebraicScheme_subscheme):
    r"""
    Construct an algebraic subscheme of projective space.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`~sage.schemes.projective.projective_space.ProjectiveSpace_field.subscheme`
        method of :class:`projective space
        <sage.schemes.projective.projective_space.ProjectiveSpace_field>`.

    INPUT:

    - ``A`` -- ambient :class:`projective space
      <sage.schemes.projective.projective_space.ProjectiveSpace_field>`

    - ``polynomials`` -- single polynomial, ideal or iterable of
      defining homogeneous polynomials

    EXAMPLES::

        sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
        sage: P.subscheme([x^2 - y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z

    TESTS::

        sage: from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective
        sage: AlgebraicScheme_subscheme_projective(P, [x^2-y*z])
        Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          x^2 - y*z
    """
    def point(self, v, check=True):
        """
        Create a point on this projective subscheme.

        INPUT:

        - ``v`` -- anything that defines a point

        - ``check`` -- boolean (default: ``True``); whether
          to check the defining data for consistency

        OUTPUT: a point of the subscheme

        EXAMPLES::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P2.subscheme([x - y, y - z])
            sage: X.point([1,1,1])
            (1 : 1 : 1)

        ::

            sage: P2.<x,y> = ProjectiveSpace(QQ, 1)
            sage: X = P2.subscheme([y])
            sage: X.point(infinity)
            (1 : 0)

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: X = P.subscheme(x^2 + 2*y^2)
            sage: X.point(infinity)
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [1, 0] do not define a point on Closed subscheme
            of Projective Space of dimension 1 over Rational Field defined by:
              x^2 + 2*y^2
        """
        from sage.rings.infinity import infinity
        if v is infinity or\
           (isinstance(v, (list, tuple)) and len(v) == 1 and v[0] is infinity):
            if self.ambient_space().dimension_relative() > 1:
                raise ValueError("%s not well defined in dimension > 1" % v)
            v = [1, 0]
        # todo: update elliptic curve stuff to take point_homset as argument
        from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
        if isinstance(self, EllipticCurve_generic):
            try:
                return self._point(self.point_homset(), v, check=check)
            except AttributeError:  # legacy code without point_homset
                return self._point(self, v, check=check)

        return self.point_homset()(v, check=check)

    def _morphism(self, *args, **kwds):
        r"""
        Construct a morphism determined by action on points of ``self``.

        For internal use only.

        INPUT:

        - same as for
          :class:`~sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space`.

        OUTPUT: :class:`~sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space`

        TESTS::

            sage: P1.<x,y> = ProjectiveSpace(1, QQ)
            sage: P2 = ProjectiveSpace(2, QQ)
            sage: H12 = P1.Hom(P2)
            sage: H12([x^2, x*y, y^2])    # indirect doctest
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 : x*y : y^2)
            sage: P1._morphism(H12, [x^2, x*y, y^2])
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 : x*y : y^2)
        """
        return self.ambient_space()._morphism(*args, **kwds)

    def dimension(self):
        """
        Return the dimension of the projective algebraic subscheme.

        OUTPUT: integer

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: P2.<x,y,z> = ProjectiveSpace(2, QQ)
            sage: P2.subscheme([]).dimension()
            2
            sage: P2.subscheme([x]).dimension()
            1
            sage: P2.subscheme([x^5]).dimension()
            1
            sage: P2.subscheme([x^2 + y^2 - z^2]).dimension()
            1
            sage: P2.subscheme([x*(x-z), y*(y-z)]).dimension()
            0

        Something less obvious::

            sage: P3.<x,y,z,w,t> = ProjectiveSpace(4, QQ)
            sage: X = P3.subscheme([x^2, x^2*y^2 + z^2*t^2,
            ....:                   z^2 - w^2, 10*x^2 + w^2 - z^2]); X
            Closed subscheme of Projective Space of dimension 4 over Rational Field defined by:
              x^2,
              x^2*y^2 + z^2*t^2,
              z^2 - w^2,
              10*x^2 - z^2 + w^2
            sage: X.dimension()                                                         # needs sage.libs.singular
            1
        """
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self.defining_ideal().dimension() - 1
            return self.__dimension

    def affine_patch(self, i, AA=None):
        r"""
        Return the `i`-th affine patch of this projective scheme.

        This is the intersection with this `i`-th affine patch of
        its ambient space.

        INPUT:

        - ``i`` -- integer between 0 and dimension of ``self``, inclusive

        - ``AA`` -- (default: ``None``) ambient affine space, this
          is constructed if it is not given

        OUTPUT:

        An affine algebraic scheme with fixed
        :meth:`embedding_morphism` equal to the default
        :meth:`projective_embedding` map`.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: PP = ProjectiveSpace(2, QQ, names='X,Y,Z')
            sage: X,Y,Z = PP.gens()
            sage: C = PP.subscheme(X^3*Y + Y^3*Z + Z^3*X)
            sage: U = C.affine_patch(0); U
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              Y^3*Z + Z^3 + Y
            sage: U.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: Y^3*Z + Z^3 + Y
              To:   Closed subscheme of Projective Space of dimension 2 over Rational Field
                    defined by: X^3*Y + Y^3*Z + X*Z^3
              Defn: Defined on coordinates by sending (Y, Z) to (1 : Y : Z)
            sage: U.projective_embedding() is U.embedding_morphism()
            True

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: X = A.subscheme([x - y*z])
            sage: Y = X.projective_embedding(1).codomain()                              # needs sage.libs.singular
            sage: Y.affine_patch(1, A).ambient_space() == A                             # needs sage.libs.singular
            True

        ::

            sage: P.<u,v,w> = ProjectiveSpace(2, ZZ)
            sage: S = P.subscheme([u^2 - v*w])
            sage: A.<x, y> = AffineSpace(2, ZZ)
            sage: S.affine_patch(1, A)                                                  # needs sage.libs.singular
            Closed subscheme of Affine Space of dimension 2 over Integer Ring defined by:
              x^2 - y
        """
        i = int(i)   # implicit type checking
        PP = self.ambient_space()
        n = PP.dimension_relative()
        if i < 0 or i > n:
            raise ValueError(f"Argument i (= {i}) must be between 0 and {n}.")
        try:
            A = self.__affine_patches[i]
            # assume that if you've passed in a new ambient affine space
            # you want to override the existing patch
            if AA is None or A.ambient_space() == AA:
                return self.__affine_patches[i]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        if AA is None:
            AA = PP.affine_patch(i)
        elif AA.dimension_relative() != n:
            raise ValueError("Affine Space must be of the dimension %s" % (n))
        phi = AA.projective_embedding(i, PP)
        polys = self.defining_polynomials()
        xi = phi.defining_polynomials()
        U = AA.subscheme([f(xi) for f in polys])
        phi = U.projective_embedding(i, PP)
        U._embedding_morphism = phi
        self.__affine_patches[i] = U
        return U

    def _best_affine_patch(self, point):
        r"""
        Return the best affine patch of the ambient projective space.

        The "best" affine patch is where you end up dividing by the
        homogeneous coordinate with the largest absolute
        value. Division by small numbers is numerically unstable.

        INPUT:

        - ``point`` -- a point of the algebraic subscheme

        OUTPUT: integer. The index of the patch. See :meth:`affine_patch`

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: S = P.subscheme(x + 2*y + 3*z)
            sage: S._best_affine_patch(P.point([0,-3,2]))
            1
            sage: S._best_affine_patch([0,-3,2])
            1

        TESTS::

            sage: F = GF(3)
            sage: P.<x,y,z> = ProjectiveSpace(F, 2)
            sage: S._best_affine_patch([0,1,2])
            2
        """
        point = list(point)
        try:
            abs_point = [abs(c) for c in point]
        except ArithmeticError:
            # our base ring does not know abs
            abs_point = point
        # find best patch
        i_max = 0
        p_max = abs_point[i_max]
        for i in range(1, len(point)):
            if abs_point[i] > p_max:
                i_max = i
                p_max = abs_point[i_max]
        return i_max

    def neighborhood(self, point):
        r"""
        Return an affine algebraic subscheme isomorphic to a
        neighborhood of the ``point``.

        INPUT:

        - ``point`` -- a point of the projective subscheme

        OUTPUT:

        An affine algebraic scheme (polynomial equations in affine
        space) ``result`` such that

        * :meth:`embedding_morphism
          <AlgebraicScheme.embedding_morphism>` is an isomorphism to a
          neighborhood of ``point``

        * :meth:`embedding_center <AlgebraicScheme.embedding_center>`
          is mapped to ``point``.

        EXAMPLES::

            sage: P.<x,y,z>= ProjectiveSpace(QQ, 2)
            sage: S = P.subscheme(x + 2*y + 3*z)
            sage: s = S.point([0,-3,2]); s
            (0 : -3/2 : 1)
            sage: patch = S.neighborhood(s); patch
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x + 3*z
            sage: patch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field
                    defined by: x + 3*z
              To:   Closed subscheme of Projective Space of dimension 2 over Rational Field
                    defined by: x + 2*y + 3*z
              Defn: Defined on coordinates by sending (x, z) to (x : -3/2 : z + 1)
            sage: patch.embedding_center()
            (0, 0)
            sage: patch.embedding_morphism()([0,0])
            (0 : -3/2 : 1)
            sage: patch.embedding_morphism()(patch.embedding_center())
            (0 : -3/2 : 1)
        """
        point = list(point)
        self._check_satisfies_equations(point)
        PP = self.ambient_space()
        n = PP.dimension()
        i = self._best_affine_patch(point)

        patch_cover = PP.affine_patch(i)
        R = patch_cover.coordinate_ring()

        phi = list(point)
        for j in range(i):
            phi[j] += R.gen(j)
        for j in range(i, n):
            phi[j + 1] += R.gen(j)

        pullback_polys = [f(phi) for f in self.defining_polynomials()]
        return patch_cover.subscheme(pullback_polys, embedding_center=[0] * n,
                                     embedding_codomain=self,
                                     embedding_images=phi)

    def is_smooth(self, point=None) -> bool:
        r"""
        Test whether the algebraic subscheme is smooth.

        INPUT:

        - ``point`` -- a point or ``None`` (default); the point to
          test smoothness at

        OUTPUT:

        boolean; if no point was specified, returns whether the
        algebraic subscheme is smooth everywhere. Otherwise,
        smoothness at the specified point is tested.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: P2.<x,y,z> = ProjectiveSpace(2, QQ)
            sage: cuspidal_curve = P2.subscheme([y^2*z - x^3])
            sage: cuspidal_curve
            Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
              -x^3 + y^2*z
            sage: cuspidal_curve.is_smooth([1,1,1])
            True
            sage: cuspidal_curve.is_smooth([0,0,1])
            False
            sage: cuspidal_curve.is_smooth()
            False
            sage: P2.subscheme([y^2*z - x^3 + z^3 + 1/10*x*y*z]).is_smooth()
            True

        TESTS::

            sage: H = P2.subscheme(x)                                                   # needs sage.libs.singular
            sage: H.is_smooth()  # one of the few cases where the cone over the subvariety is smooth                    # needs sage.libs.singular
            True
        """
        if point is not None:
            self._check_satisfies_equations(point)
            R = self.ambient_space().coordinate_ring()
            point_subs = dict(zip(R.gens(), point))
            Jac = self.Jacobian().subs(point_subs)
            return not Jac.is_zero()

        # testing smoothness everywhere tends to be expensive
        try:
            return self._smooth
        except AttributeError:
            pass
        sing_dim = self.Jacobian().dimension()
        # We really test the affine cone here; the origin is always a
        # singular point:
        self._smooth = (sing_dim <= 0)
        return self._smooth

    def orbit(self, f, N) -> list:
        r"""
        Return the orbit of this scheme by ``f``.

        If `N` is an integer it returns `[self,f(self),\ldots,f^N(self)]`.
        If `N` is a list or tuple `N=[m,k]` it returns `[f^m(self),\ldots,f^k(self)`].

        INPUT:

        - ``f`` -- a :class:`DynamicalSystem_projective` with ``self`` in ``f.domain()``

        - ``N`` -- nonnegative integer or list or tuple of two nonnegative integers

        OUTPUT: list of projective subschemes

        EXAMPLES::

            sage: # needs sage.libs.singular sage.schemes
            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: f = DynamicalSystem_projective([(x-2*y)^2, (x-2*z)^2,
            ....:                                 (x-2*w)^2, x^2])
            sage: f.orbit(P.subscheme([x]), 5)
            [Closed subscheme of Projective Space of dimension 3 over Rational Field
              defined by: x,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
              defined by: w,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
              defined by: z - w,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
              defined by: y - z,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
              defined by: x - y,
             Closed subscheme of Projective Space of dimension 3 over Rational Field
              defined by: x - w]

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
            sage: H = Hom(PS, P1)
            sage: f = H([x^2, y^2])
            sage: X = PS.subscheme([x - y])
            sage: X.orbit(f, 2)
            Traceback (most recent call last):
            ...
            TypeError: map must be a dynamical system for iteration

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])                       # needs sage.schemes
            sage: X = PS.subscheme([x - y])
            sage: X.orbit(f, [-1,2])                                                    # needs sage.schemes
            Traceback (most recent call last):
            ...
            TypeError: orbit bounds must be nonnegative
        """
        if not isinstance(f, DynamicalSystem):
            raise TypeError("map must be a dynamical system for iteration")
        if not isinstance(N, (list, tuple)):
            N = [0, N]
        N[0] = ZZ(N[0])
        N[1] = ZZ(N[1])
        if N[0] < 0 or N[1] < 0:
            raise TypeError("orbit bounds must be nonnegative")
        if N[0] > N[1]:
            return []

        Q = self
        for i in range(1, N[0] + 1):
            Q = f(Q)
        Orb = [Q]

        for i in range(N[0] + 1, N[1] + 1):
            Q = f(Q)
            Orb.append(Q)
        return Orb

    def nth_iterate(self, f, n):
        r"""
        The `n`-th forward image of this scheme by the map ``f``.

        INPUT:

        - ``f`` -- a :class:`DynamicalSystem_projective` with ``self`` in ``f.domain()``

        - ``n`` -- positive integer

        OUTPUT:

        - A subscheme in ``f.codomain()``

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: f = DynamicalSystem_projective([y^2, z^2, x^2, w^2])                  # needs sage.schemes
            sage: f.nth_iterate(P.subscheme([x - w, y - z]), 3)                         # needs sage.libs.singular sage.schemes
            Closed subscheme of Projective Space of dimension 3 over Rational Field
             defined by:
              y - z,
              x - w

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])                       # needs sage.schemes
            sage: X = PS.subscheme([x - y])
            sage: X.nth_iterate(f, -2)                                                  # needs sage.libs.singular sage.schemes
            Traceback (most recent call last):
            ...
            TypeError: must be a forward orbit

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: P2.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(PS, P2)
            sage: f = H([x^2, y^2, z^2])
            sage: X = PS.subscheme([x - y])
            sage: X.nth_iterate(f, 2)
            Traceback (most recent call last):
            ...
            TypeError: map must be a dynamical system for iteration

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])                       # needs sage.schemes
            sage: X = PS.subscheme([x - y])
            sage: X.nth_iterate(f, 2.5)                                                 # needs sage.schemes
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
        """
        n = ZZ(n)
        if n < 0:
            raise TypeError("must be a forward orbit")
        return self.orbit(f, [n, n + 1])[0]

    def _forward_image(self, f, check=True):
        r"""
        Compute the forward image of this subscheme by the morphism ``f``.

        The forward image is computed through elimination and ``f`` must be
        a morphism for this to be well defined.
        In particular, let `X = V(h_1,\ldots, h_t)` and define the ideal
        `I = (h_1,\ldots,h_t,y_0-f_0(\bar{x}), \ldots, y_n-f_n(\bar{x}))`.
        Then the elimination ideal `I_{n+1} = I \cap K[y_0,\ldots,y_n]` is a homogeneous
        ideal and `self(X) = V(I_{n+1})`.

        INPUT:

        - ``f`` -- a map whose domain contains ``self``

        - ``check`` -- boolean, if ``False`` no input checking is done

        OUTPUT: a subscheme in the codomain of ``f``

        EXAMPLES::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = End(PS)
            sage: f = H([x^2, y^2 - 2*z^2, z^2])
            sage: X = PS.subscheme(y - 2*z)
            sage: X._forward_image(f)                                                   # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field
             defined by:
              y - 2*z

        ::

            sage: set_verbose(None)
            sage: PS.<x,y,z,w> = ProjectiveSpace(ZZ, 3)
            sage: H = End(PS)
            sage: f = H([y^2, x^2, w^2, z^2])
            sage: X = PS.subscheme([z^2 + y*w, x - w])
            sage: f(X)                                                                  # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 3 over Integer Ring
             defined by:
              y - z,
              x*z - w^2

        ::

            sage: # needs sage.rings.real_mpfr
            sage: PS.<x,y,z,w> = ProjectiveSpace(CC, 3)
            sage: H = End(PS)
            sage: f = H([x^2 + y^2, y^2, z^2-y^2, w^2])
            sage: X = PS.subscheme([z - 2*w])
            sage: f(X)                                                                  # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 3 over Complex Field
             with 53 bits of precision defined by:
              y + z + (-4.00000000000000)*w

        ::

            sage: R.<t> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
            sage: H = End(P)
            sage: f = H([x^2 + 2*y*z, t^2*y^2, z^2])
            sage: f([t^2*y - z])                                                        # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Fraction Field
             of Univariate Polynomial Ring in t over Rational Field defined by:
              y - 1/(t^2)*z

        ::

            sage: # needs sage.rings.padics
            sage: set_verbose(-1)
            sage: PS.<x,y,z> = ProjectiveSpace(Qp(3), 2)
            sage: H = End(PS)
            sage: f = H([x^2, 2*y^2, z^2])
            sage: X = PS.subscheme([2*x - y, z])
            sage: f(X)                                                                  # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over 3-adic Field
             with capped relative precision 20 defined by:
              z,
              x + (1 + 3^2 + 3^4 + 3^6 + 3^8 + 3^10 + 3^12 + 3^14 + 3^16 + 3^18 + O(3^20))*y

        ::

            sage: R.<y0,y1,y2,y3> = PolynomialRing(QQ)
            sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
            sage: H = End(P)
            sage: f = H([y0*x^2 + y1*z^2, y2*y^2 + y3*z^2, z^2])
            sage: X = P.subscheme(x*z)
            sage: X._forward_image(f)                                                   # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Fraction Field
             of Multivariate Polynomial Ring in y0, y1, y2, y3 over Rational Field
             defined by:
              x*z + (-y1)*z^2

            ::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P5.<z0,z1,z2,z3,z4,z5> = ProjectiveSpace(QQ, 5)
            sage: H = Hom(P2, P5)
            sage: f = H([x^2, x*y, x*z, y^2, y*z, z^2])  # Veronese map
            sage: X = P2.subscheme([])
            sage: f(X)                                                                  # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 5 over Rational Field
             defined by:
              -z4^2 + z3*z5,
              -z2*z4 + z1*z5,
              -z2*z3 + z1*z4,
              -z2^2 + z0*z5,
              -z1*z2 + z0*z4,
              -z1^2 + z0*z3

            ::

            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P3.<u,v,w,t> = ProjectiveSpace(QQ, 3)
            sage: H = Hom(P2, P3)
            sage: X = P2.subscheme([x - y, x - z])
            sage: f = H([x^2, y^2, z^2, x*y])
            sage: f(X)                                                                  # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 3 over Rational Field
             defined by:
              w - t,
              v - t,
              u - t

            ::

            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
            sage: P2.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: H = Hom(P2, P1)
            sage: f = H([x^2, y*z])
            sage: X = P2.subscheme([x - y])
            sage: f(X)                                                                  # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

            ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(PS)
            sage: f = H([x^3, x*y^2, x*z^2])
            sage: X = PS.subscheme([x - y])
            sage: X._forward_image(f)                                                   # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: P1.<u,v> = ProjectiveSpace(QQ, 1)
            sage: Y = P1.subscheme([u - v])
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: Y._forward_image(f)                                                   # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: subscheme must be in ambient space of domain of map
        """
        if check:
            if not f.is_morphism():
                raise TypeError("map must be a morphism")
        dom = f.domain()
        codom = f.codomain()
        if check:
            if self.ambient_space() != dom:
                raise TypeError("subscheme must be in ambient space of domain of map")
        CR_dom = dom.coordinate_ring()
        CR_codom = codom.coordinate_ring()
        n = CR_dom.ngens()
        m = CR_codom.ngens()
        # can't call eliminate if the base ring is polynomial so we do it ourselves
        # with a lex ordering
        R = PolynomialRing(f.base_ring(), n + m, 'tempvar', order='lex')
        Rvars = R.gens()[0:n]
        phi = CR_dom.hom(Rvars, R)
        zero = n * [0]
        psi = R.hom(zero + list(CR_codom.gens()), CR_codom)
        # set up ideal
        L = R.ideal([phi(t) for t in self.defining_polynomials()] + [R.gen(n + i) - phi(f[i]) for i in range(m)])
        G = L.groebner_basis()  # eliminate
        newL = []
        # get only the elimination ideal portion
        for i in range(len(G) - 1, 0, -1):
            v = G[i].variables()
            if all(Rvars[j] not in v for j in range(n)):
                newL.append(psi(G[i]))
        return codom.subscheme(newL)

    def preimage(self, f, k=1, check=True):
        r"""
        The subscheme that maps to this scheme by the map `f^k`.

        In particular, `f^{-k}(V(h_1,\ldots,h_t)) = V(h_1 \circ f^k, \ldots, h_t \circ f^k)`.
        Map must be a morphism and also must be an endomorphism for `k > 1`.

        INPUT:

        - ``f`` -- a map whose codomain contains this scheme

        - ``k`` -- positive integer

        - ``check`` -- boolean; if ``False`` no input checking is done

        OUTPUT: a subscheme in the domain of ``f``

        EXAMPLES::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(PS)
            sage: f = H([y^2, x^2, z^2])
            sage: X = PS.subscheme([x - y])
            sage: X.preimage(f)                                                         # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Integer Ring
             defined by:
              -x^2 + y^2

        ::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(QQ, 4)
            sage: H = End(P)
            sage: f = H([x^2 - y^2, y^2, z^2, w^2, t^2 + w^2])
            sage: f.rational_preimages(P.subscheme([x - z, t^2, w - t]))                # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 4 over Rational Field
             defined by:
              x^2 - y^2 - z^2,
              w^4 + 2*w^2*t^2 + t^4,
              -t^2

        ::

            sage: P1.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P3.<u,v,w,t> = ProjectiveSpace(QQ, 3)
            sage: H = Hom(P1, P3)
            sage: X = P3.subscheme([u - v, 2*u - w, u + t])
            sage: f = H([x^2, y^2, x^2 + y^2, x*y])
            sage: X.preimage(f)                                                         # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 1 over Rational Field
             defined by:
              x^2 - y^2,
              x^2 - y^2,
              x^2 + x*y

        ::

            sage: P1.<x,y> = ProjectiveSpace(QQ, 1)
            sage: P3.<u,v,w,t> = ProjectiveSpace(QQ, 3)
            sage: H = Hom(P3, P1)
            sage: X = P1.subscheme([x - y])
            sage: f = H([u^2, v^2])
            sage: X.preimage(f)                                                         # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: H = End(PS)
            sage: f = H([x^2, x^2, x^2])
            sage: X = PS.subscheme([x - y])
            sage: X.preimage(f)                                                         # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: map must be a morphism

        ::

            sage: PS.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: P1.<u,v> = ProjectiveSpace(ZZ, 1)
            sage: Y = P1.subscheme([u^2 - v^2])
            sage: H = End(PS)
            sage: f = H([x^2, y^2, z^2])
            sage: Y.preimage(f)                                                         # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: subscheme must be in ambient space of codomain

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: Y = P.subscheme([x - y])
            sage: H = End(P)
            sage: f = H([x^2, y^2, z^2])
            sage: Y.preimage(f, k=2)                                                    # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field
             defined by:
              x^4 - y^4
        """
        dom = f.domain()
        codom = f.codomain()
        if check:
            if not f.is_morphism():
                raise TypeError("map must be a morphism")
            if self.ambient_space() != codom:
                raise TypeError("subscheme must be in ambient space of codomain")
            k = ZZ(k)
            if k <= 0:
                raise ValueError("k (=%s) must be a positive integer" % (k))
            if k > 1 and not f.is_endomorphism():
                raise TypeError("map must be an endomorphism")
        R = codom.coordinate_ring()
        F = f.as_dynamical_system().nth_iterate_map(k) if k > 1 else f
        dic = {R.gen(i): F[i] for i in range(codom.dimension_relative() + 1)}
        return dom.subscheme([t.subs(dic) for t in self.defining_polynomials()])

    def dual(self):
        r"""
        Return the projective dual of the given subscheme of projective space.

        INPUT:

        - ``X`` -- a subscheme of projective space. At present, ``X`` is
          required to be an irreducible and reduced hypersurface defined
          over `\QQ` or a finite field.

        OUTPUT: the dual of ``X`` as a subscheme of the dual projective space

        EXAMPLES:

        The dual of a smooth conic in the plane is also a smooth conic::

            sage: R.<x, y, z> = QQ[]
            sage: P.<x, y, z> = ProjectiveSpace(2, QQ)
            sage: I = R.ideal(x^2 + y^2 + z^2)
            sage: X = P.subscheme(I)
            sage: X.dual()                                                              # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 2 over Rational Field
             defined by:
              y0^2 + y1^2 + y2^2

        The dual of the twisted cubic curve in projective 3-space is a singular
        quartic surface. In the following example, we compute the dual of this
        surface, which by double duality is equal to the twisted cubic itself.
        The output is the twisted cubic as an intersection of three quadrics::

            sage: R.<x, y, z, w> = QQ[]
            sage: P.<x, y, z, w> = ProjectiveSpace(3, QQ)
            sage: I = R.ideal(y^2*z^2 - 4*x*z^3 - 4*y^3*w + 18*x*y*z*w - 27*x^2*w^2)
            sage: X = P.subscheme(I)
            sage: X.dual()                                                              # needs sage.libs.singular
            Closed subscheme of Projective Space of dimension 3 over
             Rational Field defined by:
              y2^2 - y1*y3,
              y1*y2 - y0*y3,
              y1^2 - y0*y2

        The singular locus of the quartic surface in the last example
        is itself supported on a twisted cubic::

            sage: X.Jacobian().radical()                                                # needs sage.libs.singular
            Ideal (z^2 - 3*y*w, y*z - 9*x*w, y^2 - 3*x*z) of Multivariate
            Polynomial Ring in x, y, z, w over Rational Field

        An example over a finite field::

            sage: R = PolynomialRing(GF(61), 'a,b,c')
            sage: P.<a, b, c> = ProjectiveSpace(2, R.base_ring())
            sage: X = P.subscheme(R.ideal(a*a + 2*b*b + 3*c*c))
            sage: X.dual()                                                              # needs sage.libs.singular sage.rings.finite_rings
            Closed subscheme of Projective Space of dimension 2 over
             Finite Field of size 61 defined by:
              y0^2 - 30*y1^2 - 20*y2^2

        TESTS::

            sage: # needs sage.rings.padics
            sage: R = PolynomialRing(Qp(3), 'a,b,c')
            sage: P.<a, b, c> = ProjectiveSpace(2, R.base_ring())
            sage: X = P.subscheme(R.ideal(a*a + 2*b*b + 3*c*c))
            sage: X.dual()
            Traceback (most recent call last):
            ...
            NotImplementedError: base ring must be QQ or a finite field
        """
        from sage.libs.singular.function_factory import ff

        K = self.base_ring()
        if not (isinstance(K, RationalField) or K in Fields().Finite()):
            raise NotImplementedError("base ring must be QQ or a finite field")
        I = self.defining_ideal()
        m = I.ngens()
        n = I.ring().ngens() - 1
        if (m != 1 or (n < 1) or I.is_zero()
                or I.is_trivial() or not I.is_prime()):
            raise NotImplementedError("At the present, the method is only"
                                      " implemented for irreducible and"
                                      " reduced hypersurfaces and the given"
                                      " list of generators for the ideal must"
                                      " have exactly one element.")
        R = PolynomialRing(K, 'x', n + 1)
        from sage.schemes.projective.projective_space import ProjectiveSpace
        Pd = ProjectiveSpace(n, K, 'y')
        Rd = Pd.coordinate_ring()
        x = R.variable_names()
        y = Rd.variable_names()
        S = PolynomialRing(K, x + y + ('t',))
        if S.has_coerce_map_from(I.ring()):
            T = PolynomialRing(K, 'w', n + 1)
            I_S = (I.change_ring(T)).change_ring(S)
        else:
            I_S = I.change_ring(S)
        f_S = I_S.gens()[0]
        z = S.gens()
        J = I_S
        for i in range(n + 1):
            J = J + S.ideal(z[-1] * f_S.derivative(z[i]) - z[i + n + 1])

        try:
            sat = ff.elim__lib.sat_with_exp
        except NameError:
            sat = ff.elim__lib.sat

        max_ideal = S.ideal(z[n + 1: 2 * n + 2])
        J_sat_gens = sat(J, max_ideal)[0]
        J_sat = S.ideal(J_sat_gens)
        L = J_sat.elimination_ideal(z[0: n + 1] + (z[-1],))
        return Pd.subscheme(L.change_ring(Rd))

    def degree(self):
        r"""
        Return the degree of this projective subscheme.

        If `P(t) = a_{m}t^m + \ldots + a_{0}` is the Hilbert
        polynomial of this subscheme, then the degree is `a_{m} m!`.

        OUTPUT: integer

        EXAMPLES::

            sage: P.<x,y,z,w,t,u> = ProjectiveSpace(QQ, 5)
            sage: X = P.subscheme([x^7 + x*y*z*t^4 - u^7])
            sage: X.degree()                                                            # needs sage.libs.singular
            7

            sage: P.<x,y,z,w> = ProjectiveSpace(GF(13), 3)
            sage: X = P.subscheme([y^3 - w^3, x + 7*z])
            sage: X.degree()                                                            # needs sage.libs.singular
            3

            sage: # needs sage.libs.singular sage.schemes
            sage: P.<x,y,z,w,u> = ProjectiveSpace(QQ, 4)
            sage: C = P.curve([x^7 - y*z^3*w^2*u, w*x^2 - y*u^2, z^3 + y^3])
            sage: C.degree()
            63
        """
        P = self.defining_ideal().hilbert_polynomial()
        return P.leading_coefficient() * P.degree().factorial()

    def intersection_multiplicity(self, X, P):
        r"""
        Return the intersection multiplicity of this subscheme and the subscheme ``X`` at the point ``P``.

        This uses the intersection_multiplicity function for affine subschemes on affine patches of this subscheme
        and ``X`` that contain ``P``.

        INPUT:

        - ``X`` -- subscheme in the same ambient space as this subscheme

        - ``P`` -- a point in the intersection of this subscheme with ``X``

        OUTPUT: integer

        EXAMPLES::

            sage: # needs sage.schemes
            sage: P.<x,y,z> = ProjectiveSpace(GF(5), 2)
            sage: C = Curve([x^4 - z^2*y^2], P)
            sage: D = Curve([y^4*z - x^5 - x^3*z^2], P)
            sage: Q1 = P([0,1,0])
            sage: C.intersection_multiplicity(D, Q1)                                    # needs sage.libs.singular
            4
            sage: Q2 = P([0,0,1])
            sage: C.intersection_multiplicity(D, Q2)                                    # needs sage.libs.singular
            6

        ::

            sage: # needs sage.rings.number_field
            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^4 + 1)
            sage: P.<x,y,z,w> = ProjectiveSpace(K, 3)
            sage: X = P.subscheme([x^2 + y^2 - z*w])
            sage: Y = P.subscheme([y*z - x*w, z - w])
            sage: Q1 = P([b^2,1,0,0])
            sage: X.intersection_multiplicity(Y, Q1)                                    # needs sage.libs.singular
            1
            sage: Q2 = P([1/2*b^3 - 1/2*b, 1/2*b^3 - 1/2*b, 1, 1])
            sage: X.intersection_multiplicity(Y, Q2)                                    # needs sage.libs.singular
            1

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: X = P.subscheme([x^2 - z^2, y^3 - w*x^2])
            sage: Y = P.subscheme([w^2 - 2*x*y + z^2, y^2 - w^2])
            sage: Q = P([1,1,-1,1])
            sage: X.intersection_multiplicity(Y, Q)                                     # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 3
            over Rational Field defined by: z^2 + w^2 - 2*y, y^2 - w^2)
            must be proper and finite
        """
        try:
            self.ambient_space()(P)
        except TypeError:
            raise TypeError("(={}) must be a point in the ambient space of this subscheme and (={})".format(P, X))
        # find an affine chart of the ambient space of this curve that contains P
        n = self.ambient_space().dimension_relative()
        for i in range(n + 1):
            if P[i] != 0:
                break
        X1 = self.affine_patch(i)
        X2 = X.affine_patch(i)
        return X1.intersection_multiplicity(X2, X1(P.dehomogenize(i)))

    def multiplicity(self, P):
        r"""
        Return the multiplicity of ``P`` on this subscheme.

        This is computed as the multiplicity of the corresponding point on an affine patch of this subscheme
        that contains ``P``. This subscheme must be defined over a field. An error is returned if ``P``
        not a point on this subscheme.

        INPUT:

        - ``P`` -- a point on this subscheme

        OUTPUT: integer

        EXAMPLES::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(QQ, 4)
            sage: X = P.subscheme([y^2 - x*t, w^7 - t*w*x^5 - z^7])
            sage: Q1 = P([0,0,1,1,1])
            sage: X.multiplicity(Q1)                                                    # needs sage.libs.singular
            1
            sage: Q2 = P([1,0,0,0,0])
            sage: X.multiplicity(Q2)                                                    # needs sage.libs.singular
            3
            sage: Q3 = P([0,0,0,0,1])
            sage: X.multiplicity(Q3)                                                    # needs sage.libs.singular
            7

        ::

            sage: # needs sage.rings.real_mpfr
            sage: P.<x,y,z,w> = ProjectiveSpace(CC, 3)
            sage: X = P.subscheme([z^5*x^2*w - y^8])
            sage: Q = P([2,0,0,1])
            sage: X.multiplicity(Q)                                                     # needs sage.libs.singular
            5

        ::

            sage: # needs sage.libs.singular sage.schemes
            sage: P.<x,y,z,w> = ProjectiveSpace(GF(29), 3)
            sage: C = Curve([y^17 - x^5*w^4*z^8, x*y - z^2], P)
            sage: Q = P([3,0,0,1])
            sage: C.multiplicity(Q)
            8
        """
        if self.base_ring() not in Fields():
            raise TypeError("subscheme must be defined over a field")

        # check whether P is a point on this subscheme
        try:
            P = self(P)
        except TypeError:
            raise TypeError(f"(={P}) is not a point on (={self})")

        # find an affine chart of the ambient space of self that contains P
        i = 0
        while P[i] == 0:
            i += 1
        X = self.affine_patch(i)
        return X.multiplicity(X(P.dehomogenize(i)))

    def veronese_embedding(self, d, CS=None, order='lex'):
        r"""
        Return the degree ``d`` Veronese embedding of this projective subscheme.

        INPUT:

        - ``d`` -- positive integer

        - ``CS`` -- (default: ``None``) a projective ambient space to embed
          into. If the projective ambient space of this subscheme is of
          dimension `N`, the dimension of ``CS`` must be
          `\binom{N + d}{d} - 1`. This is constructed if not specified.

        - ``order`` -- string (default: ``'lex'``); a monomial order to use to
          arrange the monomials defining the embedding. The monomials will be
          arranged from greatest to least with respect to this order.

        OUTPUT: a scheme morphism from this subscheme to its image by the
        degree ``d`` Veronese embedding

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: L = P.subscheme([y - x])
            sage: v = L.veronese_embedding(2); v                                        # needs sage.libs.singular
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 2
                    over Rational Field defined by: -x + y
              To:   Closed subscheme of Projective Space of dimension 5
                    over Rational Field defined by:
                      -x4^2 + x3*x5,
                      x2 - x4,
                      x1 - x3,
                      x0 - x3
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^2 : x*y : x*z : y^2 : y*z : z^2)
            sage: v.codomain().degree()                                                 # needs sage.libs.singular
            2
            sage: C = P.subscheme([y*z - x^2])
            sage: C.veronese_embedding(2).codomain().degree()                           # needs sage.libs.singular
            4

        twisted cubic::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: Q.<u,v,s,t> = ProjectiveSpace(QQ, 3)
            sage: P.subscheme([]).veronese_embedding(3, Q)                              # needs sage.libs.singular
            Scheme morphism:
              From: Closed subscheme of Projective Space of dimension 1
                    over Rational Field defined by: (no polynomials)
              To:   Closed subscheme of Projective Space of dimension 3
                    over Rational Field defined by:
                      -s^2 + v*t,
                      -v*s + u*t,
                      -v^2 + u*s
              Defn: Defined on coordinates by sending (x : y) to
                    (x^3 : x^2*y : x*y^2 : y^3)
        """
        # construct map between projective spaces
        v = self.ambient_space().veronese_embedding(d, CS, order)
        # return this map restricted to self and its image
        return Hom(self, v(self))(v.defining_polynomials())


class AlgebraicScheme_subscheme_projective_field(AlgebraicScheme_subscheme_projective):
    """
    Algebraic subschemes of projective spaces defined over fields.
    """
    def _morphism(self, *args, **kwds):
        r"""
        Construct a morphism determined by action on points of ``self``.

        For internal use only.

        INPUT:

        - same as for
          :class:`~sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space`.

        OUTPUT: :class:`~sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space`

        TESTS::

            sage: P1.<x,y> = ProjectiveSpace(1, QQ)
            sage: P2 = ProjectiveSpace(2, QQ)
            sage: H12 = P1.Hom(P2)
            sage: H12([x^2, x*y, y^2])    # indirect doctest
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 : x*y : y^2)
            sage: P1._morphism(H12, [x^2, x*y, y^2])
            Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x : y) to (x^2 : x*y : y^2)
        """
        return SchemeMorphism_polynomial_projective_subscheme_field(*args, **kwds)

    def Chow_form(self):
        r"""
        Return the Chow form associated to this subscheme.

        For a `k`-dimensional subvariety of `\mathbb{P}^N` of degree `D`.
        The `(N-k-1)`-dimensional projective linear subspaces of `\mathbb{P}^N`
        meeting `X` form a hypersurface in the Grassmannian `G(N-k-1,N)`.
        The homogeneous form of degree `D` defining this hypersurface in Plucker
        coordinates is called the Chow form of `X`.

        The base ring needs to be a number field, finite field, or `\QQbar`.

        ALGORITHM:

        For a `k`-dimension subscheme `X` consider the `k+1` linear forms
        `l_i = u_{i0}x_0 + \cdots + u_{in}x_n`. Let `J` be the ideal in the
        polynomial ring `K[x_i,u_{ij}]` defined by the equations of `X` and the `l_i`.
        Let `J'` be the saturation of `J` with respect to the irrelevant ideal of
        the ambient projective space of `X`. The elimination ideal `I = J' \cap K[u_{ij}]`
        is a principal ideal, let `R` be its generator. The Chow form is obtained by
        writing `R` as a polynomial in Plucker coordinates (i.e. bracket polynomials).
        [DS1994]_.

        OUTPUT: a homogeneous polynomial

        EXAMPLES::

            sage: P.<x0,x1,x2,x3> = ProjectiveSpace(GF(17), 3)
            sage: X = P.subscheme([x3 + x1, x2 - x0, x2 - x3])
            sage: X.Chow_form()                                                         # needs sage.libs.singular
            t0 - t1 + t2 + t3

        ::

            sage: P.<x0,x1,x2,x3> = ProjectiveSpace(QQ, 3)
            sage: X = P.subscheme([x3^2 - 101*x1^2 - 3*x2*x0])
            sage: X.Chow_form()                                                         # needs sage.libs.singular
            t0^2 - 101*t2^2 - 3*t1*t3

        ::

            sage: # needs sage.libs.singular
            sage: P.<x0,x1,x2,x3> = ProjectiveSpace(QQ, 3)
            sage: X = P.subscheme([x0*x2 - x1^2, x0*x3 - x1*x2, x1*x3 - x2^2])
            sage: Ch = X.Chow_form(); Ch
            t2^3 + 2*t2^2*t3 + t2*t3^2 - 3*t1*t2*t4 - t1*t3*t4 + t0*t4^2 + t1^2*t5
            sage: Y = P.subscheme_from_Chow_form(Ch, 1); Y
            Closed subscheme of Projective Space of dimension 3 over Rational Field
             defined by:
              x2^2*x3 - x1*x3^2,                    -x2^3 + x0*x3^2,
              -x2^2*x3 + x1*x3^2,                   x1*x2*x3 - x0*x3^2,
              3*x1*x2^2 - 3*x0*x2*x3,               -2*x1^2*x3 + 2*x0*x2*x3,
              -3*x1^2*x2 + 3*x0*x1*x3,              x1^3 - x0^2*x3,
              x2^3 - x1*x2*x3,                      -3*x1*x2^2 + 2*x1^2*x3 + x0*x2*x3,
              2*x0*x2^2 - 2*x0*x1*x3,               3*x1^2*x2 - 2*x0*x2^2 - x0*x1*x3,
              -x0*x1*x2 + x0^2*x3,                  -x0*x1^2 + x0^2*x2,
              -x1^3 + x0*x1*x2,                     x0*x1^2 - x0^2*x2
            sage: I = Y.defining_ideal()
            sage: I.saturation(I.ring().ideal(list(I.ring().gens())))[0]
            Ideal (x2^2 - x1*x3, x1*x2 - x0*x3, x1^2 - x0*x2)
             of Multivariate Polynomial Ring in x0, x1, x2, x3 over Rational Field
        """
        I = self.defining_ideal()
        P = self.ambient_space()
        R = P.coordinate_ring()
        N = P.dimension() + 1
        d = self.dimension()
        # create the ring for the generic linear hyperplanes
        # u0x0 + u1x1 + ...
        SS = PolynomialRing(R.base_ring(), 'u', N * (d + 1), order='lex')
        vars = SS.variable_names() + R.variable_names()
        S = PolynomialRing(R.base_ring(), vars, order='lex')
        n = S.ngens()
        newcoords = [S.gen(n - N + t) for t in range(N)]
        # map the generators of the subscheme into the ring with the hyperplane variables
        phi = R.hom(newcoords, S)
        phi(self.defining_polynomials()[0])
        # create the dim(X)+1 linear hyperplanes
        l = []
        for i in range(d + 1):
            t = 0
            for j in range(N):
                t += S.gen(N * i + j) * newcoords[j]
            l.append(t)
        # intersect the hyperplanes with X
        J = phi(I) + S.ideal(l)
        # saturate the ideal with respect to the irrelevant ideal
        J2 = J.saturation(S.ideal([phi(u) for u in R.gens()]))[0]
        # eliminate the original variables to be left with the hyperplane coefficients 'u'
        E = J2.elimination_ideal(newcoords)
        # create the plucker coordinates
        D = binomial(N, N - d - 1)  # number of plucker coordinates
        tvars = [f't{i}' for i in range(D)]  # plucker coordinates
        T = PolynomialRing(R.base_ring(), tvars + list(S.variable_names()), order='lex')
        L = []
        coeffs = [T.gen(i) for i in range(len(tvars), N * (d + 1) + len(tvars))]
        M = matrix(T, d + 1, N, coeffs)
        i = 0
        for c in M.minors(d + 1):
            L.append(T.gen(i) - c)
            i += 1
        # create the ideal that we can use for eliminating to get a polynomial
        # in the plucker coordinates (brackets)
        br = T.ideal(L)
        # create a mapping into a polynomial ring over the plucker coordinates
        # and the hyperplane coefficients
        psi = S.hom(coeffs + [0 for _ in range(N)], T)
        E2 = T.ideal([psi(u) for u in E.gens()] + br)
        # eliminate the hyperplane coefficients
        CH = E2.elimination_ideal(coeffs)
        # CH should be a principal ideal, but because of the relations among
        # the plucker coordinates, the elimination will probably have several generators

        # get the relations among the plucker coordinates
        rel = br.elimination_ideal(coeffs)
        # reduce CH with respect to the relations
        reduced = [f.reduce(rel) for f in CH.gens()]
        # find the principal generator

        # polynomial ring in just the plucker coordinates
        T2 = PolynomialRing(R.base_ring(), tvars)
        alp = T.hom(tvars + (N * (d + 1) + N) * [0], T2)
        # get the degrees of the reduced generators of CH
        degs = [u.degree() for u in reduced]
        mind = max(degs)
        # need the smallest degree form that did not reduce to 0
        for d in degs:
            if d < mind and d > 0:
                mind = d
        ind = degs.index(mind)
        CF = reduced[ind]  # this should be the Chow form of X
        # check that it is correct (i.e., it is a principal generator for CH + the relations)
        rel2 = rel + [CF]
        assert all(f in rel2 for f in CH.gens()), "did not find a principal generator"
        return alp(CF)

    def global_height(self, prec=None):
        """
        Return the (projective) global height of the subscheme.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default ``RealField`` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: NF.<a> = NumberField(x^2 - 5)
            sage: P.<x,y,z> = ProjectiveSpace(NF, 2)
            sage: X = P.subscheme([x^2 + y*z, 2*y*z, 3*x*y])
            sage: X.global_height()                                                     # needs sage.libs.singular
            0.000000000000000

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme([z^2 - 101*y^2 - 3*x*z])
            sage: X.global_height()             # long time                             # needs sage.libs.singular
            4.61512051684126
        """
        return self.Chow_form().global_height(prec)

    def local_height(self, v, prec=None):
        """
        Return the (projective) local height of the subscheme.

        INPUT:

        - ``v`` -- a prime or prime ideal of the base ring

        - ``prec`` -- desired floating point precision (default:
          default ``RealField`` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: NF.<a> = NumberField(x^2 - 5)
            sage: I = NF.ideal(3)
            sage: P.<x,y,z> = ProjectiveSpace(NF, 2)
            sage: X = P.subscheme([3*x*y - 5*x*z, y^2])
            sage: X.local_height(I)                                                     # needs sage.libs.singular
            0.000000000000000

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme([z^2 - 101*y^2 - 3*x*z])
            sage: X.local_height(2)                                                     # needs sage.libs.singular
            0.000000000000000
        """
        return self.Chow_form().local_height(v, prec)

    def local_height_arch(self, i, prec=None):
        """
        Return the local height at the ``i``-th infinite place of the subscheme.

        INPUT:

        - ``i`` -- integer

        - ``prec`` -- desired floating point precision (default:
          default ``RealField`` precision)

        OUTPUT: a real number

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: NF.<a> = NumberField(x^2 - 5)
            sage: P.<x,y,z> = ProjectiveSpace(NF, 2)
            sage: X = P.subscheme([x^2 + y*z, 3*x*y])
            sage: X.local_height_arch(1)                                                # needs sage.libs.singular
            0.0000000000000000000000000000000

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: X = P.subscheme([z^2 - 101*y^2 - 3*x*z])
            sage: X.local_height_arch(1)                                                # needs sage.libs.singular
            4.61512051684126
        """
        return self.Chow_form().local_height_arch(i, prec)
