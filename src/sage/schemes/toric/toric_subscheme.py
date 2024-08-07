# sage.doctest: needs sage.geometry.polyhedron sage.graphs
r"""
Subschemes of toric space

AUTHORS:

- David Kohel (2005): initial version.

- William Stein (2005): initial version.

- Andrey Novoseltsev (2010-05-17): subschemes of toric varieties.
"""

# ****************************************************************************
# Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme

lazy_import("sage.calculus.functions", "jacobian")


class AlgebraicScheme_subscheme_toric(AlgebraicScheme_subscheme):
    r"""
    Construct an algebraic subscheme of a toric variety.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`~ToricVariety_field.subscheme` method of :class:`toric
        varieties <sage.schemes.toric.variety.ToricVariety_field>`.

    INPUT:

    - ``toric_variety`` -- ambient :class:`toric variety
      <ToricVariety_field>`

    - ``polynomials`` -- single polynomial, list, or ideal of defining
      polynomials in the coordinate ring of ``toric_variety``

    OUTPUT: an :class:`algebraic subscheme of a toric variety
    <AlgebraicScheme_subscheme_toric>`.

    TESTS::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1.inject_variables()
        Defining s, t, x, y
        sage: import sage.schemes.toric.toric_subscheme as SCM
        sage: X = SCM.AlgebraicScheme_subscheme_toric(
        ....:       P1xP1, [x*s + y*t, x^3 + y^3])
        sage: X
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3

    A better way to construct the same scheme as above::

        sage: P1xP1.subscheme([x*s + y*t, x^3 + y^3])
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3
    """

    # Implementation note: if the toric variety is affine you should
    # construct instances of the derived class
    # AlgebraicScheme_subscheme_affine_toric instead.

    def __init__(self, toric_variety, polynomials):
        r"""
        See :class:`AlgebraicScheme_subscheme_toric` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: import sage.schemes.toric.toric_subscheme as SCM
            sage: X = SCM.AlgebraicScheme_subscheme_toric(
            ....:       P1xP1, [x*s + y*t, x^3 + y^3])
            sage: X
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 4 affine patches defined by:
              s*x + t*y,
              x^3 + y^3
        """
        # Just to make sure that keyword arguments will be passed correctly
        super().__init__(toric_variety, polynomials)

    def _morphism(self, *args, **kwds):
        r"""
        Construct a morphism determined by action on points of ``self``.

        INPUT:

        - same as for
          :class:`~sage.schemes.toric.morphism.SchemeMorphism_polynomial_toric_variety`.

        OUTPUT: a :class:`~sage.schemes.toric.morphism.SchemeMorphism_polynomial_toric_variety`

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: P1 = P1xP1.subscheme(s - t)
            sage: H = P1.Hom(P1xP1)
            sage: H([s, s, x, y])
            Scheme morphism:
              From: Closed subscheme of 2-d CPR-Fano toric variety
              covered by 4 affine patches defined by:
              s - t
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [s : t : x : y] to
                    [s : s : x : y]

            sage: sbar, tbar, xbar, ybar = P1.coordinate_ring().gens()                  # needs sage.libs.singular
            sage: P1._morphism(H, [sbar, sbar, xbar, ybar])                             # needs sage.libs.singular
            Scheme morphism:
              From: Closed subscheme of 2-d CPR-Fano toric variety
              covered by 4 affine patches defined by:
              s - t
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [s : t : x : y] to
                    [t : t : x : y]
        """
        from sage.schemes.toric.morphism import SchemeMorphism_polynomial_toric_variety
        return SchemeMorphism_polynomial_toric_variety(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        r"""
        Construct a Hom-set for ``self``.

        INPUT:

        - same as for
          :class:`~sage.schemes.toric.homset.SchemeHomset_points_toric_field`.

        OUTPUT: a :class:`~sage.schemes.toric.homset.SchemeHomset_points_subscheme_toric_field`

        TESTS::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: quadric = P2.subscheme([x^2 + y^2 + z^2])
            sage: quadric._point_homset(Spec(QQ), quadric)
            Set of rational points of Closed subscheme of 2-d CPR-Fano
            toric variety covered by 3 affine patches defined by:
              x^2 + y^2 + z^2
            sage: type(quadric.point_set())
            <class 'sage.schemes.toric.homset.SchemeHomset_points_subscheme_toric_field_with_category'>
        """
        from sage.schemes.toric.homset import SchemeHomset_points_subscheme_toric_field
        return SchemeHomset_points_subscheme_toric_field(*args, **kwds)

    def fan(self):
        """
        Return the fan of the ambient space.

        OUTPUT: a fan

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P(2)
            sage: E = P2.subscheme([x^2 + y^2 + z^2])
            sage: E.fan()
            Rational polyhedral fan in 2-d lattice N
        """
        return self.ambient_space().fan()

    def affine_patch(self, i):
        r"""
        Return the ``i``-th affine patch of ``self`` as an affine
        toric algebraic scheme.

        INPUT:

        - ``i`` -- integer; index of a generating cone of the fan of the
          ambient space of ``self``

        OUTPUT:

        - subscheme of an affine :class:`toric variety
          <sage.schemes.toric.variety.ToricVariety_field>`
          corresponding to the pull-back of ``self`` by the embedding
          morphism of the ``i``-th :meth:`affine patch of the ambient
          space
          <sage.schemes.toric.variety.ToricVariety_field.affine_patch>`
          of ``self``.

        The result is cached, so the ``i``-th patch is always the same object
        in memory.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: patch1 = P1xP1.affine_patch(1)
            sage: patch1.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d CPR-Fano toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [t : x] to [1 : t : x : 1]
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: P1 = P1xP1.subscheme(x - y)
            sage: subpatch = P1.affine_patch(1)
            sage: subpatch
            Closed subscheme of 2-d affine toric variety defined by:
              x - 1
        """
        i = int(i)   # implicit type checking
        try:
            return self._affine_patches[i]
        except AttributeError:
            self._affine_patches = {}
        except KeyError:
            pass
        ambient_patch = self.ambient_space().affine_patch(i)
        phi_p = ambient_patch.embedding_morphism().defining_polynomials()
        patch = ambient_patch.subscheme(
                            [p(phi_p) for p in self.defining_polynomials()])
        patch._embedding_morphism = patch.hom(phi_p, self, check=False)
        self._affine_patches[i] = patch
        return patch

    def affine_algebraic_patch(self, cone=None, names=None):
        r"""
        Return the affine patch corresponding to ``cone`` as an affine
        algebraic scheme.

        INPUT:

        - ``cone`` -- a :class:`Cone
          <sage.geometry.cone.ConvexRationalPolyhedralCone>` `\sigma`
          of the fan. It can be omitted for an affine toric variety,
          in which case the single generating cone is used.

        OUTPUT:

        An :class:`affine algebraic subscheme
        <sage.schemes.affine.affine_subscheme.AlgebraicScheme_subscheme_affine>`
        corresponding to the patch `\mathop{Spec}(\sigma^\vee \cap M)`
        associated to the cone `\sigma`.

        See also :meth:`affine_patch`, which expresses the patches as
        subvarieties of affine toric varieties instead.

        REFERENCES:

        ..

            David A. Cox, "The Homogeneous Coordinate Ring of a Toric
            Variety", Lemma 2.2.
            :arxiv:`alg-geom/9210008v2`

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: cone = P2.fan().generating_cone(0)
            sage: V = P2.subscheme(x^3 + y^3 + z^3)
            sage: V.affine_algebraic_patch(cone)
            Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              z0^3 + z1^3 + 1

            sage: # needs fpylll sage.libs.singular
            sage: cone = Cone([(0,1), (2,1)])
            sage: A2Z2.<x,y> = AffineToricVariety(cone)
            sage: A2Z2.affine_algebraic_patch()
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              -z0*z1 + z2^2
            sage: V = A2Z2.subscheme(x^2 + y^2 - 1)
            sage: patch = V.affine_algebraic_patch();  patch
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              -z0*z1 + z2^2,
              z0 + z1 - 1
            sage: nbhd_patch = V.neighborhood([1,0]).affine_algebraic_patch();  nbhd_patch
            Closed subscheme of Affine Space of dimension 3 over Rational Field defined by:
              -z0*z1 + z2^2,
              z0 + z1 - 1
            sage: nbhd_patch.embedding_center()
            (0, 1, 0)

        Here we got two defining equations. The first one describes
        the singularity of the ambient space and the second is the
        pull-back of `x^2+y^2-1` ::

            sage: lp = LatticePolytope([(1,0,0), (1,1,0), (1,1,1), (1,0,1), (-2,-1,-1)],
            ....:                      lattice=ToricLattice(3))
            sage: X.<x,y,u,v,t> = CPRFanoToricVariety(Delta_polar=lp)
            sage: Y = X.subscheme(x*v + y*u + t)
            sage: cone = Cone([(1,0,0), (1,1,0), (1,1,1), (1,0,1)])
            sage: Y.affine_algebraic_patch(cone)
            Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              z0*z2 - z1*z3,
              z1 + z3 + 1
        """
        from sage.modules.free_module_element import vector
        from sage.misc.misc_c import prod
        ambient = self.ambient_space()
        fan = ambient.fan()
        if cone is None:
            assert ambient.is_affine()
            cone = fan.generating_cone(0)
        else:
            cone = fan.embed(cone)
        # R/I = C[sigma^dual cap M]
        R, I, dualcone = ambient._semigroup_ring(cone, names)

        # inhomogenize the Cox homogeneous polynomial with respect to the given cone
        inhomogenize = {ambient.coordinate_ring().gen(i): 1
                        for i in range(fan.nrays())
                        if i not in cone.ambient_ray_indices()}
        polynomials = [p.subs(inhomogenize) for p in self.defining_polynomials()]

        # map the monomial x^{D_m} to m, see reference.
        n_rho_matrix = cone.rays().matrix()

        def pullback_polynomial(p):
            result = R.zero()
            for coefficient, monomial in p:
                exponent = monomial.exponents()[0]
                exponent = [ exponent[i] for i in cone.ambient_ray_indices() ]
                exponent = vector(ZZ,exponent)
                m = n_rho_matrix.solve_right(exponent)
                assert all(x in ZZ for x in m), \
                    'The polynomial '+str(p)+' does not define a ZZ-divisor!'
                m_coeffs = dualcone.Hilbert_coefficients(m)
                result += coefficient * prod(R.gen(i)**m_coeffs[i]
                                             for i in range(0,R.ngens()))
            return result

        # construct the affine algebraic scheme to use as patch
        polynomials = [pullback_polynomial(_) for _ in polynomials]
        from sage.schemes.affine.affine_space import AffineSpace
        patch_cover = AffineSpace(R)
        polynomials = list(I.gens()) + polynomials
        polynomials = [x for x in polynomials if not x.is_zero()]
        patch = patch_cover.subscheme(polynomials)

        # TODO: If the cone is not smooth, then the coordinate_ring()
        # of the affine toric variety is wrong; it should be the
        # G-invariant part. So we can't construct the embedding
        # morphism in that case.
        if cone.is_smooth():
            x = ambient.coordinate_ring().gens()
            phi = []
            for i in range(0,fan.nrays()):
                if i in cone.ambient_ray_indices():
                    phi.append(pullback_polynomial(x[i]))
                else:
                    phi.append(1)
            patch._embedding_morphism = patch.hom(phi, self)
        else:
            patch._embedding_morphism = (NotImplementedError,
               'I only know how to construct embedding morphisms for smooth patches')

        try:
            point = self.embedding_center()
        except AttributeError:
            return patch

        # it remains to find the preimage of point
        # map m to the monomial x^{D_m}, see reference.
        F = ambient.coordinate_ring().fraction_field()
        image = []
        for m in dualcone.Hilbert_basis():
            x_Dm = prod([ F.gen(i)**(m*n) for i,n in enumerate(fan.rays()) ])
            image.append(x_Dm)
        patch._embedding_center = tuple( f(list(point)) for f in image )
        return patch

    def _best_affine_patch(self, point):
        r"""
        Return the best affine patch of the ambient toric variety.

        INPUT:

        - ``point`` -- a point of the algebraic subscheme

        OUTPUT: integer. The index of the patch. See :meth:`affine_patch`

        EXAMPLES::

            sage: P.<x,y,z> = toric_varieties.P2()
            sage: S = P.subscheme(x + 2*y + 3*z)
            sage: S._best_affine_patch(P.point([2,-3,0]))
            1
            sage: S._best_affine_patch([2,-3,0])
            1
        """
        # TODO: this method should pick a "best" patch in the sense
        # that it is numerically stable to dehomogenize, see the
        # corresponding method for projective varieties.
        point = list(point)
        zeros = {i for i, coord in enumerate(point) if coord == 0}
        for cone_idx, cone in enumerate(self.ambient_space().fan().generating_cones()):
            if zeros.issubset(cone.ambient_ray_indices()):
                return cone_idx
        assert False, 'The point must not have been a point of the toric variety.'

    def neighborhood(self, point):
        r"""
        Return a toric algebraic scheme isomorphic to neighborhood of
        the ``point``.

        INPUT:

        - ``point`` -- a point of the toric algebraic scheme

        OUTPUT:

        An affine toric algebraic scheme (polynomial equations in an
        affine toric variety) with fixed
        :meth:`~AlgebraicScheme.embedding_morphism` and
        :meth:`~AlgebraicScheme.embedding_center`.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: P.<x,y,z> = toric_varieties.P2()
            sage: S = P.subscheme(x + 2*y + 3*z)
            sage: s = S.point([0,-3,2]); s
            [0 : -3 : 2]
            sage: patch = S.neighborhood(s); patch
            Closed subscheme of 2-d affine toric variety defined by:
              x + 2*y + 6
            sage: patch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d affine toric variety defined by: x + 2*y + 6
              To:   Closed subscheme of 2-d CPR-Fano toric variety
                    covered by 3 affine patches defined by: x + 2*y + 3*z
              Defn: Defined on coordinates by sending [x : y] to [-2*y - 6 : y : 2]
            sage: patch.embedding_center()
            [0 : -3]
            sage: patch.embedding_morphism()(patch.embedding_center())
            [0 : -3 : 2]

        A more complicated example::

            sage: # needs sage.libs.singular
            sage: dP6.<x0,x1,x2,x3,x4,x5> = toric_varieties.dP6()
            sage: twoP1 = dP6.subscheme(x0*x3)
            sage: patch = twoP1.neighborhood([0,1,2, 3,4,5]); patch
            Closed subscheme of 2-d affine toric variety defined by:
              3*x0
            sage: patch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d affine toric variety defined by: 3*x0
              To:   Closed subscheme of 2-d CPR-Fano toric variety
                    covered by 6 affine patches defined by: x0*x3
              Defn: Defined on coordinates by sending [x0 : x1] to [0 : x1 : 2 : 3 : 4 : 5]
            sage: patch.embedding_center()
            [0 : 1]
            sage: patch.embedding_morphism()(patch.embedding_center())
            [0 : 1 : 2 : 3 : 4 : 5]
        """
        point = list(point)
        self._check_satisfies_equations(point)
        PP = self.ambient_space()
        fan = PP.fan()
        cone_idx = self._best_affine_patch(point)
        cone = fan.generating_cone(cone_idx)

        patch_cover = PP.affine_patch(cone_idx)
        R = patch_cover.coordinate_ring()
        phi = []
        point_preimage = []
        for i in range(fan.nrays()):
            try:
                ray_index = cone.ambient_ray_indices().index(i)
                phi.append(R.gen(ray_index))
                point_preimage.append(point[i])
            except ValueError:
                phi.append(point[i])
        pullback_polys = [f(phi) for f in self.defining_polynomials()]
        patch = patch_cover.subscheme(pullback_polys)
        S = patch.coordinate_ring()
        phi_reduced = [S(t) for t in phi]

        patch._embedding_center = patch(point_preimage)
        patch._embedding_morphism = patch.hom(phi_reduced,self)
        return patch

    def dimension(self):
        """
        Return the dimension of ``self``.

        OUTPUT: integer; if ``self`` is empty, `-1` is returned

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: P1 = P1xP1.subscheme(s - t)
            sage: P1.dimension()
            1
            sage: P1xP1.subscheme([s - t, (s-t)^2]).dimension()
            1
            sage: P1xP1.subscheme([s, t]).dimension()
            -1
        """
        if '_dimension' in self.__dict__:
            return self._dimension
        npatches = self.ambient_space().fan().ngenerating_cones()
        dims = [ self.affine_patch(i).dimension() for i in range(0,npatches) ]
        self._dimension = max(dims)
        return self._dimension

    def is_smooth(self, point=None):
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
            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: cuspidal_curve = P2.subscheme([y^2*z - x^3])
            sage: cuspidal_curve
            Closed subscheme of 2-d CPR-Fano toric variety covered by 3 affine patches defined by:
              -x^3 + y^2*z
            sage: cuspidal_curve.is_smooth([1,1,1])
            True
            sage: cuspidal_curve.is_smooth([0,0,1])
            False
            sage: cuspidal_curve.is_smooth()
            False

        Any sufficiently generic cubic hypersurface is smooth::

            sage: P2.subscheme([y^2*z-x^3+z^3+1/10*x*y*z]).is_smooth()                  # needs sage.libs.singular
            True

        A more complicated example::

            sage: # needs sage.libs.singular
            sage: dP6.<x0,x1,x2,x3,x4,x5> = toric_varieties.dP6()
            sage: disjointP1s = dP6.subscheme(x0*x3)
            sage: disjointP1s.is_smooth()
            True
            sage: intersectingP1s = dP6.subscheme(x0*x1)
            sage: intersectingP1s.is_smooth()
            False

        A smooth hypersurface in a compact singular toric variety::

            sage: # needs sage.libs.singular
            sage: lp = LatticePolytope([(1,0,0), (1,1,0), (1,1,1), (1,0,1), (-2,-1,-1)],
            ....:                      lattice=ToricLattice(3))
            sage: X.<x,y,u,v,t> = CPRFanoToricVariety(Delta_polar=lp)
            sage: Y = X.subscheme(x*v + y*u + t)
            sage: cone = Cone([(1,0,0), (1,1,0), (1,1,1), (1,0,1)])
            sage: Y.is_smooth()
            True
        """
        if point is not None:
            toric_patch = self.neighborhood(point)
            return toric_patch.is_smooth(toric_patch.embedding_center())

        # testing smoothness everywhere tends to be expensive
        if '_smooth' in self.__dict__:
            return self._smooth
        npatches = self.ambient_space().fan().ngenerating_cones()
        self._smooth = all(self.affine_patch(i).is_smooth() for i in range(0,npatches))
        return self._smooth

    def is_nondegenerate(self):
        r"""
        Check if ``self`` is nondegenerate.

        OUTPUT:

        Whether the variety is nondegenerate, that is, the intersection
        with every open torus orbit is smooth and transversal.

        EXAMPLES::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: P2.subscheme([x^3 + y^3 + z^3]).is_nondegenerate()
            True
            sage: P2.subscheme([x*y*z]).is_nondegenerate()
            False
            sage: X = P2.subscheme([(x-y)^2*(x+y) + x*y*z + z^3])
            sage: X.is_smooth()
            True
            sage: X.is_nondegenerate()
            False

        A K3 surface in `\mathbf{P}^1 \times \mathbf{P}^1 \times \mathbf{P}^1`::

            sage: diamond = lattice_polytope.cross_polytope(3)
            sage: fan = FaceFan(diamond)
            sage: P1xP1xP1 = ToricVariety(fan)
            sage: z0, z1, z2, z3, z4, z5 = P1xP1xP1.gens()
            sage: t = 5
            sage: F = z0^2*z1^2*z2^2 + z1^2*z2^2*z3^2 + z0^2*z2^2*z4^2\
            ....: + z2^2*z3^2*z4^2 + t*z0*z1*z2*z3*z4*z5 + z0^2*z1^2*z5^2\
            ....: + z1^2*z3^2*z5^2 + z0^2*z4^2*z5^2 + z3^2*z4^2*z5^2
            sage: X = P1xP1xP1.subscheme([F])
            sage: X.is_smooth()
            True
            sage: X.is_nondegenerate()
            False

        Taking a random change of variables breaks the symmetry, but
        makes the surface nondegenerate::

            sage: F1 = F.subs(z0=1*z0 + 1*z3, z3=1*z0 + 2*z3,
            ....:             z1=-2*z1 + -1*z4, z4=1*z1 + 2*z4,
            ....:             z2=-3*z2 + -1*z5, z5=-3*z2 + 2*z5)
            sage: Y = P1xP1xP1.subscheme([F1])
            sage: Y.is_smooth()
            True
            sage: Y.is_nondegenerate()
            True

        This example is from Hamm, :arxiv:`1106.1826v1`. It addresses
        an issue raised at :issue:`15239`::

            sage: X = toric_varieties.WP([1,4,2,3], names='z0 z1 z2 z3')
            sage: X.inject_variables()
            Defining z0, z1, z2, z3
            sage: g0 = z1^3 + z2^6 + z3^4
            sage: g = g0 - 2*z3^2*z0^6 + z2*z0^10 + z0^12
            sage: Y = X.subscheme([g])
            sage: Y.is_nondegenerate()
            False

        It handles nonzero characteristic::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: f = x^5 + 2*x*y^4 + y^5 - 2*y^3*z^2 + x*z^4 - 2*z^5
            sage: P2.change_ring(GF(5)).subscheme([f]).is_nondegenerate()
            True
            sage: P2.change_ring(GF(7)).subscheme([f]).is_nondegenerate()
            False

        TESTS:

        Some corner cases discussed at :issue:`15239`::

            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: P2.subscheme([]).is_nondegenerate()
            False
            sage: P2.subscheme([x]).is_nondegenerate()
            False
        """
        X = self.ambient_space()
        fan = X.fan()
        SR = X.Stanley_Reisner_ideal()
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(X.base_ring(), fan.nrays() + SR.ngens(), 't')
        slack = R.gens()[fan.nrays():]
        SR = SR.change_ring(R)

        def restrict(cone):
            patch = {}
            divide = {}
            for i in cone.ambient_ray_indices():
                patch[R.gen(i)] = R.zero()   # restrict to torus orbit
                # divide out highest power of R.gen(i)
                divide[R.gen(i)] = R.one()
            ideal = self.defining_ideal().change_ring(R)
            ideal = ideal.subs(patch)
            mat = jacobian(ideal.gens(), R.gens()[:fan.nrays()])
            minors = mat.minors(self.codimension())
            minors = tuple([ideal.reduce(m) for m in minors])
            Jac_patch = R.ideal(ideal.gens() + minors)
            SR_patch = R.ideal([monomial * slack[i] - R.one()
                                for i, monomial in
                                enumerate(SR.subs(divide).gens())])
            return ideal, Jac_patch + SR_patch

        for dim in range(0, fan.dim() + 1):
            for cone in fan(dim):
                ideal1, ideal2 = restrict(cone)
                if ideal1.is_zero() or ideal2.dimension() != -1:
                    return False

        return True

    def is_schon(self):
        r"""
        Check if ``self`` is schon (nondegenerate).

        See :meth:`is_nondegenerate` for further documentation.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: P2.<x,y,z> = toric_varieties.P2()
            sage: X = P2.subscheme([(x-y)^2*(x+y) + x*y*z + z^3])
            sage: X.is_smooth()
            True
            sage: X.is_schon()
            False
        """
        return self.is_nondegenerate()


class AlgebraicScheme_subscheme_affine_toric(AlgebraicScheme_subscheme_toric):
    r"""
    Construct an algebraic subscheme of an affine toric variety.

    .. WARNING::

        You should not create objects of this class directly. The preferred
        method to construct such subschemes is to use
        :meth:`~ToricVariety_field.subscheme` method of
        :class:`toric varieties <ToricVariety_field>`.

    INPUT:

    - ``toric_variety`` -- ambient :class:`affine toric variety
      <ToricVariety_field>`

    - ``polynomials`` -- single polynomial, list, or ideal of defining
      polynomials in the coordinate ring of ``toric_variety``

    OUTPUT:

    A :class:`algebraic subscheme of an affine toric variety
    <AlgebraicScheme_subscheme_affine_toric>`.

    TESTS::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1.inject_variables()
        Defining s, t, x, y
        sage: import sage.schemes.toric.toric_subscheme as SCM
        sage: X = SCM.AlgebraicScheme_subscheme_toric(
        ....:       P1xP1, [x*s + y*t, x^3 + y^3])
        sage: X
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3

    A better way to construct the same scheme as above::

        sage: P1xP1.subscheme([x*s + y*t, x^3 + y^3])
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          s*x + t*y,
          x^3 + y^3
    """

    def __init__(self, toric_variety, polynomials):
        r"""
        See :class:`AlgebraicScheme_subscheme_toric` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.inject_variables()
            Defining s, t, x, y
            sage: import sage.schemes.toric.toric_subscheme as SCM
            sage: X = SCM.AlgebraicScheme_subscheme_toric(
            ....:       P1xP1, [x*s + y*t, x^3 + y^3])
            sage: X
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 4 affine patches defined by:
              s*x + t*y,
              x^3 + y^3
        """
        assert toric_variety.is_affine(), 'The toric variety must be affine!'
        # Just to make sure that keyword arguments will be passed correctly
        super().__init__(toric_variety, polynomials)

    def dimension(self):
        """
        Return the dimension of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: P1xP1.<s0,s1,t0,t1> = toric_varieties.P1xP1()
            sage: P1 = P1xP1.subscheme(s0 - s1)
            sage: P1.dimension()
            1

        A more complicated example where the ambient toric variety is
        not smooth::

            sage: # needs sage.libs.singular
            sage: X.<x,y> = toric_varieties.A2_Z2()
            sage: X.is_smooth()
            False
            sage: Y = X.subscheme([x*y, x^2]); Y
            Closed subscheme of 2-d affine toric variety defined by:
              x*y,
              x^2
            sage: Y.dimension()
            1
        """
        if '_dimension' in self.__dict__:
            return self._dimension

        if self.ambient_space().is_smooth():
            self._dimension = self.defining_ideal().dimension()
        else:
            self._dimension = self.affine_algebraic_patch().dimension()
        return self._dimension

    def is_smooth(self, point=None):
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
            sage: A2.<x,y> = toric_varieties.A2()
            sage: cuspidal_curve = A2.subscheme([y^2 - x^3])
            sage: cuspidal_curve
            Closed subscheme of 2-d affine toric variety defined by:
              -x^3 + y^2
            sage: cuspidal_curve.is_smooth([1,1])
            True
            sage: cuspidal_curve.is_smooth([0,0])
            False
            sage: cuspidal_curve.is_smooth()
            False
            sage: circle = A2.subscheme(x^2 + y^2 - 1)
            sage: circle.is_smooth([1,0])
            True
            sage: circle.is_smooth()
            True

        A more complicated example where the ambient toric variety is
        not smooth::

            sage: # needs sage.libs.singular
            sage: X.<x,y> = toric_varieties.A2_Z2()    # 2-d affine space mod Z/2
            sage: X.is_smooth()
            False
            sage: Y = X.subscheme([x*y, x^2])   # (twice the x=0 curve) mod Z/2
            sage: Y
            Closed subscheme of 2-d affine toric variety defined by:
              x*y,
              x^2
            sage: Y.dimension()   # Y is a Weil divisor but not Cartier
            1
            sage: Y.is_smooth()
            True
            sage: Y.is_smooth([0,0])
            True
        """
        if point is not None:
            self._check_satisfies_equations(point)
            if self.ambient_space().is_smooth():
                R = self.ambient_space().coordinate_ring()
                point_subs = dict(zip(R.gens(), point))
                Jac = self.Jacobian().subs(point_subs)
                return not Jac.is_zero()
            else:
                self._embedding_center = self.point(point)
                affine = self.affine_algebraic_patch()
                return affine.is_smooth(affine.embedding_center())

        # testing smoothness everywhere tends to be expensive
        if '_smooth' in self.__dict__:
            return self._smooth

        if self.ambient_space().is_smooth():
            sing_dim = self.Jacobian().dimension()
            self._smooth = (sing_dim == -1)
        else:
            self._smooth = self.affine_algebraic_patch().is_smooth()

        return self._smooth
