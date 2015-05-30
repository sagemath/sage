r"""
Library of commonly used, famous, or interesting polytopes

This module gathers several constructors of polytopes that can be reached
through ``polytopes.<tab>``. For example, here is the hypercube in dimension 5::

    sage: polytopes.hypercube(5)
    A 5-dimensional polyhedron in ZZ^5 defined as the convex hull of 32 vertices

The following constructions are available

.. csv-table::
    :class: contentstable
    :widths: 30
    :delim: |

    :meth:`~sage.geometry.polyhedron.library.Polytopes.Birkhoff_polytope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.buckyball`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cross_polytope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cube`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cuboctahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cyclic_polytope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.flow_polytope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.great_rhombicuboctahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.hypercube`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.hypersimplex`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.icosahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.Kirkman_icosahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.parallelotope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.pentakis_dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.permutahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.regular_polygon`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.rhombic_dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.simplex`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.small_rhombicuboctahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.twenty_four_cell`

REFERENCES:

..  [Fetter2012]
    Hans L. Fetter,
    "A Polyhedron Full of Surprises",
    Mathematics Magazine 85 (2012), no. 5, 334-342.
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#                     2011 Volker Braun <vbraun.name@gmail.com>
#                     2015 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import itertools

from sage.rings.all import Integer, ZZ, RR, QQ, RDF, AA, QQbar
from sage.combinat.permutation import Permutations
from sage.groups.perm_gps.permgroup_named import AlternatingGroup
from sage.misc.decorators import rename_keyword
from sage.misc.superseded import deprecated_function_alias

from constructor import Polyhedron

from sage.graphs.digraph import DiGraph

def zero_sum_projection(d):
    r"""
    Return a matrix corresponding to the projection on the orthogonal of
    `(1,1,\ldots,1)` in dimension `d`.

    The projection maps the orthonormal basis `(1,-1,0,\ldots,0) / \sqrt(2)`,
    `(1,1,-1,0,\ldots,0) / \sqrt(3)`, \ldots, `(1,1,\ldots,1,-1) / \sqrt(d)` to
    the canonical basis in `\RR^{d-1}`.

    OUTPUT:
    
    A matrix of dimensions `(d-1)\times d` defined over :class:`RDF
    <sage.rings.real_double.RealDoubleField_class>`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.library import zero_sum_projection
        sage: zero_sum_projection(2)
        [ 0.7071067811865475 -0.7071067811865475]
        sage: zero_sum_projection(3)
        [ 0.7071067811865475 -0.7071067811865475                 0.0]
        [ 0.4082482904638631  0.4082482904638631 -0.8164965809277261]
    """
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
    basis = [vector(RDF,[1]*i + [-i] + [0]*(d-i-1)) for i in range(1,d)]
    return matrix(RDF, [v / v.norm() for v in basis])

def project_points(*points):
    """
    Projects a set of points into a vector space of dimension one less.

    The projection is isometric to the orthogonal projection on the hyperplane
    made of zero sum vector. Hence, if the set of points have all equal sums,
    then their projection is isometric (as a set of points).

    The projection used is the matrix given by :func:`zero_sum_projection`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.library import project_points
        sage: project_points([2,-1,3,2])          # abs tol 1e-15
        [(2.1213203435596424, -2.041241452319315, -0.577350269189626)]
        sage: project_points([1,2,3],[3,3,5])     # abs tol 1e-15
        [(-0.7071067811865475, -1.2247448713915892), (0.0, -1.6329931618554523)]

    These projections are compatible with the restriction. More precisely, given
    a vector `v`, the projection of `v` restricted to the first `i` coordinates
    will be equal to the projection of the first `i+1` coordinates of `v`::

        sage: project_points([1,2])    # abs tol 1e-15
        [(-0.7071067811865475)]
        sage: project_points([1,2,3])  # abs tol 1e-15
        [(-0.7071067811865475, -1.2247448713915892)]
        sage: project_points([1,2,3,4])     # abs tol 1e-15
        [(-0.7071067811865475, -1.2247448713915892, -1.7320508075688776)]
        sage: project_points([1,2,3,4,0])   # abs tol 1e-15
        [(-0.7071067811865475, -1.2247448713915892, -1.7320508075688776, 2.23606797749979)]

    Check that it is (almost) an isometry::

        sage: V = map(vector, IntegerVectors(n=5,length=3))
        sage: P = project_points(*V)
        sage: for i in range(21):
        ....:     for j in range(21):
        ....:         assert abs((V[i]-V[j]).norm() - (P[i]-P[j]).norm()) < 0.00001
    """
    if not points:
        return []
    from sage.modules.free_module_element import vector
    vecs = [vector(RDF,p) for p in points]
    m = zero_sum_projection(len(vecs[0]))
    return [m*v for v in vecs]

class Polytopes():
    """
    A class of constructors for commonly used, famous, or interesting
    polytopes.
    """

    flow_polytope = staticmethod(DiGraph.flow_polytope)

    def regular_polygon(self, n, exact=True, base_ring=None):
        """
        Return a regular polygon with `n` vertices.

        INPUT:

        - ``n`` -- a positive integer, the number of vertices.

        - ``exact`` -- (boolean, default ``True``) if ``False`` floating point
          numbers are used for coordinates.

        - ``base_ring`` -- a ring in which the coordinates will lie. It is
          ``None`` by default. If it is not provided and ``exact`` is ``True``
          then it will be the field of real algebraic number, if ``exact`` is
          ``False`` it will be the real double field.

        EXAMPLES::

            sage: octagon = polytopes.regular_polygon(8)
            sage: octagon
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 8 vertices
            sage: octagon.n_vertices()
            8
            sage: v = octagon.volume()
            sage: v
            2.828427124746190?
            sage: v == 2*QQbar(2).sqrt()
            True

        Its non exact version::

            sage: polytopes.regular_polygon(3, exact=False).vertices()
            (A vertex at (0.0, 1.0),
             A vertex at (0.8660254038, -0.5),
             A vertex at (-0.8660254038, -0.5))
            sage: polytopes.regular_polygon(25, exact=False).n_vertices()
            25
        """
        n = ZZ(n)
        if n <= 2:
            raise ValueError("n (={}) must be an integer greater than 2".format(n))

        if base_ring is None:
            if exact:
                base_ring = AA
            else:
                base_ring = RDF

        try:
            omega = 2*base_ring.pi() / n
            verts = [((i*omega).sin(), (i*omega).cos()) for i in range(n)]
        except AttributeError:
            z = QQbar.zeta(n)
            verts = [(base_ring((z**k).imag()), base_ring((z**k).real())) for k in range(n)]

        return Polyhedron(vertices=verts, base_ring=base_ring)

    def Birkhoff_polytope(self, n):
        """
        Return the Birkhoff polytope with `n!` vertices.

        The vertices of this polyhedron are the (flattened) `n` by `n`
        permutation matrices. So the ambient vector space has dimension `n^2`
        but the dimension of the polyhedron is `(n-1)^2`.

        INPUT:

        - ``n`` -- a positive integer giving the size of the permutation matrices.

        .. SEEALSO::

            :meth:`sage.matrix.matrix2.Matrix.as_sum_of_permutations` -- return
            the current matrix as a sum of permutation matrices

        EXAMPLES::

            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: b3.f_vector()
            (1, 6, 15, 18, 9, 1)
            sage: print b3.ambient_dim(), b3.dim()
            9 4
            sage: b3.is_lattice_polytope()
            True
            sage: p3 = b3.ehrhart_polynomial()     # optional - latte_int
            sage: p3                               # optional - latte_int
            1/8*t^4 + 3/4*t^3 + 15/8*t^2 + 9/4*t + 1
            sage: [p3(i) for i in [1,2,3,4]]       # optional - latte_int
            [6, 21, 55, 120]
            sage: [len((i*b3).integral_points()) for i in [1,2,3,4]]
            [6, 21, 55, 120]

            sage: b4 = polytopes.Birkhoff_polytope(4)
            sage: print b4.n_vertices(), b4.ambient_dim(), b4.dim()
            24 16 9
        """
        from itertools import permutations
        verts = []
        for p in permutations(range(n)):
            verts.append( [ZZ.one() if p[i]==j else ZZ.zero() for j in range(n) for i in range(n) ] )
        return Polyhedron(vertices=verts, base_ring=ZZ)

    @rename_keyword(deprecation=18213, dim_n='dim')
    def simplex(self, dim=3, project=False):
        """
        Return the ``dim`` dimensional simplex.

        The `d`-simplex is the convex hull in `\RR^{d+1}` of the standard basis
        `(1,0,\ldots,0)`, `(0,1,\ldots,0)`, \ldots, `(0,0,\ldots,1)`. For more
        information, see the :wikipedia:`Simplex`.

        INPUT:

        - ``dim`` -- The dimension of the simplex, a positive
          integer.

        - ``project`` -- (boolean, default ``False``) if ``True``, the polytope
          is (isometrically) projected to a vector space of dimension ``dim-1``.
          This operation turns the coordinates into floating point
          approximations and corresponds to the projection given by the matrix
          from :func:`zero_sum_projection`.

        EXAMPLES::

            sage: s5 = polytopes.simplex(5)
            sage: s5
            A 5-dimensional polyhedron in ZZ^6 defined as the convex hull of 6 vertices
            sage: s5.f_vector()
            (1, 6, 15, 20, 15, 6, 1)

            sage: s5 = polytopes.simplex(5, project=True)
            sage: s5
            A 5-dimensional polyhedron in RDF^5 defined as the convex hull of 6 vertices

        Its volume is `\sqrt{d+1} / d!`::

            sage: s5 = polytopes.simplex(5, project=True)
            sage: s5.volume()      # abs tol 1e-10
            0.0204124145231931
            sage: sqrt(6.) / factorial(5)
            0.0204124145231931

            sage: s6 = polytopes.simplex(6, project=True)
            sage: s6.volume()      # abs tol 1e-10
            0.00367465459870082
            sage: sqrt(7.) / factorial(6)
            0.00367465459870082
        """
        verts = list((ZZ ** (dim+1)).basis())
        if project: verts = project_points(*verts)
        return Polyhedron(vertices=verts)

    n_simplex = deprecated_function_alias(18213, simplex)

    def icosahedron(self, exact=True, base_ring=None):
        """
        Return an icosahedron with edge length 1.

        The icosahedron is one of the Platonic sold. It has 20 faces and is dual
        to the :meth:`dodecahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- (optional) the ring in which the coordinates will
          belong to.  Note that this ring must contain `\sqrt(5)`. If it is not
          provided and ``exact=True`` it will be the number field
          `\QQ[\sqrt(5)]` and if ``exact=False`` it will be the real double
          field.

        EXAMPLES::

            sage: ico = polytopes.icosahedron()
            sage: ico.f_vector()
            (1, 12, 30, 20, 1)
            sage: ico.volume()
            5/12*sqrt5 + 5/4

        Its non exact version::

            sage: ico = polytopes.icosahedron(exact=False)
            sage: ico.base_ring()
            Real Double Field
            sage: ico.volume()
            2.1816949907715726

        A version using `AA <sage.rings.qqbar.AlgebraicRealField>`::

            sage: ico = polytopes.icosahedron(base_ring=AA)   # long time
            sage: ico.base_ring()                             # long time
            Algebraic Real Field
            sage: ico.volume()                                # long time
            2.181694990624913?

        Note that if base ring is provided it must contain the square root of
        `5`. Otherwise you will get an error::

            sage: polytopes.icosahedron(base_ring=QQ)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 1/4*sqrt(5) + 1/4 to a rational
        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            if base_ring is None:
                base_ring = RDF
            g = (1 + base_ring(5).sqrt()) / 2

        r12 = base_ring.one() / 2
        z = base_ring.zero()
        pts = [[z, s1*r12, s2*g/2] for s1,s2 in itertools.product([1,-1],repeat=2)]
        verts = [p(v) for p in AlternatingGroup(3) for v in pts]
        return Polyhedron(vertices=verts, base_ring=base_ring)


    def dodecahedron(self, exact=True, base_ring=None):
        """
        Return a dodecahedron.

        The dodecahedron is the Platonic solid dual to the :meth:`icosahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- (optional) the ring in which the coordinates will
          belong to.  Note that this ring must contain `\sqrt(5)`. If it is not
          provided and ``exact=True`` it will be the number field
          `\QQ[\sqrt(5)]` and if ``exact=False`` it will be the real double
          field.

        EXAMPLES::

            sage: d12 = polytopes.dodecahedron()
            sage: d12.f_vector()
            (1, 20, 30, 12, 1)
            sage: d12.volume()
            -176*sqrt5 + 400
            sage: numerical_approx(_)
            6.45203596003699

            sage: d12 = polytopes.dodecahedron(exact=False)
            sage: d12.base_ring()
            Real Double Field

        Here is an error with a field that does not contain `\sqrt(5)`::

            sage: polytopes.dodecahedron(base_ring=QQ)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 1/4*sqrt(5) + 1/4 to a rational
        """
        return self.icosahedron(exact=exact, base_ring=base_ring).polar()

    def small_rhombicuboctahedron(self, exact=True, base_ring=None):
        """
        Return the (small) rhombicuboctahedron.

        The rhombicuboctahedron is an Archimedean solid with 24 vertices and 26
        faces. See the :wikipedia:`Rhombicuboctahedron` for more information.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False`` it
          will be the real double field.

        EXAMPLES::

            sage: sr = polytopes.small_rhombicuboctahedron()
            sage: sr.f_vector()
            (1, 24, 48, 26, 1)
            sage: sr.volume()
            80/3*sqrt2 + 32

        The faces are `8` equilateral triangles and `18` squares::

            sage: sum(1 for f in sr.faces(2) if len(f.vertices()) == 3)
            8
            sage: sum(1 for f in sr.faces(2) if len(f.vertices()) == 4)
            18

        Its non exact version::

            sage: sr = polytopes.small_rhombicuboctahedron(False)
            sage: sr
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 24
            vertices
            sage: sr.f_vector()
            (1, 24, 48, 26, 1)
        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(2, 'sqrt2')
            sqrt2 = K.gen()
            base_ring = K
        else:
            if base_ring is None:
                base_ring = RDF
            sqrt2 = base_ring(2).sqrt()

        one = base_ring.one()
        a = sqrt2 + one
        verts = []
        verts.extend([s1*one, s2*one, s3*a] for s1,s2,s3 in itertools.product([1,-1], repeat=3))
        verts.extend([s1*one, s3*a, s2*one] for s1,s2,s3 in itertools.product([1,-1], repeat=3))
        verts.extend([s1*a, s2*one, s3*one] for s1,s2,s3 in itertools.product([1,-1], repeat=3))
        return Polyhedron(vertices=verts)


    def great_rhombicuboctahedron(self, exact=True, base_ring=None):
        """
        Return the great rhombicuboctahedron.

        The great rohombicuboctahedron (or truncated cuboctahedron) is an
        Archimedean solid with 48 vertices and 26 faces. For more information
        see the :wikipedia:`Truncated_cuboctahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False`` it
          will be the real double field.

        EXAMPLES::

            sage: gr = polytopes.great_rhombicuboctahedron()  # long time ~ 3sec
            sage: gr.f_vector()                               # long time
            (1, 48, 72, 26, 1)

        A faster implementation is obtained by setting ``exact=False``::

            sage: gr = polytopes.great_rhombicuboctahedron(exact=False)
            sage: gr.f_vector()
            (1, 48, 72, 26, 1)

        Its faces are 4 squares, 8 regular hexagons and 6 regular octagons::

            sage: sum(1 for f in gr.faces(2) if len(f.vertices()) == 4)
            12
            sage: sum(1 for f in gr.faces(2) if len(f.vertices()) == 6)
            8
            sage: sum(1 for f in gr.faces(2) if len(f.vertices()) == 8)
            6
        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(2, 'sqrt2')
            sqrt2 = K.gen()
            base_ring = K
        else:
            if base_ring is None:
                base_ring = RDF
            sqrt2 = base_ring(2).sqrt()

        one = base_ring.one()
        v1 = sqrt2 + 1
        v2 = 2*sqrt2 + 1
        verts = [ [s1*z1, s2*z2, s3*z3]
                       for z1,z2,z3 in itertools.permutations([one,v1,v2])
                       for s1,s2,s3 in itertools.product([1,-1], repeat=3)]
        return Polyhedron(vertices=verts, base_ring=base_ring)

    def rhombic_dodecahedron(self):
        """
        Return the rhombic dodecahedron.

        The rhombic dodecahedron is a a polytope  dual to the cuboctahedron. It
        has 14 vertices and 12 faces. For more information see
        the :wikipedia:`Rhombic_dodecahedron`.

        .. SEEALSO::

            :meth:`cuboctahedron`

        EXAMPLES::

            sage: rd = polytopes.rhombic_dodecahedron()
            sage: rd.f_vector()
            (1, 14, 24, 12, 1)

        Its faces are 12 quadrilaterals (not all identical)::

            sage: sum(1 for f in rd.faces(2) if len(f.vertices()) == 4)
            12

        Some more computations::

            sage: p = rd.ehrhart_polynomial()    # optional - latte_int
            sage: p                              # optional - latte_int
            16*t^3 + 12*t^2 + 4*t + 1
            sage: [p(i) for i in [1,2,3,4]]      # optional - latte_int
            [33, 185, 553, 1233]
            sage: [len((i*rd).integral_points()) for i in [1,2,3,4]]
            [33, 185, 553, 1233]
        """
        v = [[2,0,0],[-2,0,0],[0,2,0],[0,-2,0],[0,0,2],[0,0,-2]]
        v.extend((itertools.product([1,-1], repeat=3)))
        return Polyhedron(vertices=v, base_ring=ZZ)

    def cuboctahedron(self):
        """
        Return the cuboctahedron.

        The cuboctahedron is an Archimedean solid with 12 vertices and 14 faces
        dual to the rhombic dodecahedron. It can be defined as the convex hull
        of the twelve vertices `(0, \pm 1, \pm 1)`, `(\pm 1, 0, \pm 1)` and
        `(\pm 1, \pm 1, 0)`. For more information, see the
        :wikipedia:`Cuboctahedron`.

        .. SEEALSO::

            :meth:`rhombic_dodecahedron`

        EXAMPLES::

            sage: co = polytopes.cuboctahedron()
            sage: co.f_vector()
            (1, 12, 24, 14, 1)

        Its faces are 8 triangles and 6 squares::

            sage: sum(1 for f in co.faces(2) if len(f.vertices()) == 3)
            8
            sage: sum(1 for f in co.faces(2) if len(f.vertices()) == 4)
            6

        Some more computation::

            sage: co.volume()
            20/3
            sage: co.ehrhart_polynomial()      # optional - latte_int
            20/3*t^3 + 8*t^2 + 10/3*t + 1
        """
        v = [ [ 0, -1, -1], [ 0, 1,-1], [ 1,-1, 0],
              [ 1,  1,  0], [ 1, 0, 1], [ 1, 0,-1],
              [ 0,  1,  1], [ 0,-1, 1], [-1, 0, 1],
              [-1,  1,  0], [-1, 0,-1], [-1,-1, 0] ]
        return Polyhedron(vertices=v, base_ring=ZZ)

    def buckyball(self, exact=True, base_ring=None):
        """
        Return the bucky ball.

        The bucky ball, also known as the truncated icosahedron is an Archimedean solid.
        It has 32 faces and 60 vertices.

        .. SEEALSO::

            :meth:`icosahedron`

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False`` it
          will be the real double field.

        EXAMPLES::

            sage: bb = polytopes.buckyball()   # long time - 6secs
            sage: bb.f_vector()                # long time
            (1, 60, 90, 32, 1)
            sage: bb.base_ring()               # long time
            Number Field in sqrt5 with defining polynomial x^2 - 5

        A much faster implementation using floating point approximations::

            sage: bb = polytopes.buckyball(exact=False)
            sage: bb.f_vector()
            (1, 60, 90, 32, 1)
            sage: bb.base_ring()
            Real Double Field

        Its faces are 5 regular pentagons and 6 regular hexagons::

            sage: sum(1 for f in bb.faces(2) if len(f.vertices()) == 5)
            12
            sage: sum(1 for f in bb.faces(2) if len(f.vertices()) == 6)
            20
        """
        return self.icosahedron(exact=exact, base_ring=base_ring).edge_truncation()

    def pentakis_dodecahedron(self, exact=True, base_ring=None):
        """
        Return the pentakis dodecahedron.

        The pentakis dodecahedron (orkisdodecahedron) is a face-regular,
        vertex-uniform polytope dual to the truncated icosahedron.  It has 60
        faces and 32 vertices. See the :wikipedia:`Pentakis_dodecahedron` for more
        information.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False`` it
          will be the real double field.

        EXAMPLES::

            sage: pd = polytopes.pentakis_dodecahedron()    # long time - ~10 sec
            sage: pd.n_vertices()                           # long time
            32
            sage: pd.n_inequalities()                       # long time
            60

        A much faster implementation is obtained when setting ``exact=False``::

            sage: pd = polytopes.pentakis_dodecahedron(exact=False)
            sage: pd.n_vertices()
            32
            sage: pd.n_inequalities()
            60

        The 60 are triangles::

            sage: all(len(f.vertices()) == 3 for f in pd.faces(2))
            True
        """
        return self.buckyball(exact=exact, base_ring=base_ring).polar()

    def Kirkman_icosahedron(self):
        """
        Return the Kirkman icosahedron.

        The Kirkman icosahedron is a 3-polytope with integer coordinates: `(\pm
        9, \pm 6, \pm 6)`, `(\pm 12, \pm 4, 0)`, `(0, \pm 12, \pm 8)`, `(\pm 6,
        0, \pm 12)`. See [Fetter2012]_ for more information.

        EXAMPLES::

            sage: ki = polytopes.Kirkman_icosahedron()
            sage: ki.f_vector()
            (1, 20, 38, 20, 1)

            sage: ki.volume()
            6528

            sage: vertices = ki.vertices()
            sage: edges = [[vector(edge[0]),vector(edge[1])] for edge in ki.bounded_edges()]
            sage: edge_lengths = [norm(edge[0]-edge[1]) for edge in edges]
            sage: union(edge_lengths)
            [7, 8, 9, 11, 12, 14, 16]
        """
        vertices = [[9, 6, 6], [-9, 6, 6], [9, -6, 6], [9, 6, -6],
                    [-9, -6, 6], [-9, 6, -6], [9, -6, -6], [-9, -6, -6],
                    [12, 4, 0], [-12, 4, 0], [12, -4, 0], [-12, -4, 0],
                    [0, 12, 8], [0, -12, 8], [0, 12, -8], [0, -12, -8],
                    [6, 0, 12], [-6, 0, 12], [6, 0, -12], [-6, 0, -12]]
        return Polyhedron(vertices=vertices, base_ring=ZZ)

    def twenty_four_cell(self):
        """
        Return the standard 24-cell polytope.

        The 24-cell polyhedron (also called icositetrachoron or octaplex) is a
        regular polyhedron in 4-dimension. For more information see
        the :wikipedia:`24-cell`.

        EXAMPLES::

            sage: p24 = polytopes.twenty_four_cell()
            sage: p24.f_vector()
            (1, 24, 96, 96, 24, 1)
            sage: v = next(p24.vertex_generator())
            sage: for adj in v.neighbors(): print adj
            A vertex at (-1/2, -1/2, -1/2, 1/2)
            A vertex at (-1/2, -1/2, 1/2, -1/2)
            A vertex at (-1, 0, 0, 0)
            A vertex at (-1/2, 1/2, -1/2, -1/2)
            A vertex at (0, -1, 0, 0)
            A vertex at (0, 0, -1, 0)
            A vertex at (0, 0, 0, -1)
            A vertex at (1/2, -1/2, -1/2, -1/2)

            sage: p24.volume()
            2
        """
        q12 = QQ((1,2))
        verts = list(itertools.product([q12,-q12], repeat=4))
        B4 = (ZZ**4).basis()
        verts.extend(v for v in B4)
        verts.extend(-v for v in B4)
        return Polyhedron(vertices=verts)

    def six_hundred_cell(self, exact=False):
        """
        Return the standard 600-cell polytope.

        The 600-cell is a 4-dimensional regular polytope. In many ways this is
        an analogue of the icosahedron.

        .. WARNING::

            The coordinates are not exact by default. The computation with exact
            coordinates takes a huge amount of time.

        INPUT:

        - ``exact`` - (boolean, default ``False``) if ``True`` use exact
          coordinates instead of floating point approximations

        EXAMPLES::

            sage: p600 = polytopes.six_hundred_cell()
            sage: p600
            A 4-dimensional polyhedron in RDF^4 defined as the convex hull of 120 vertices
            sage: p600.f_vector()
            (1, 120, 720, 1200, 600, 1)

        Computation with exact coordinates is currently too long to be useful::

            sage: p600 = polytopes.six_hundred_cell(exact=True) # not tested - very long time
            sage: len(list(p600.bounded_edges()))               # not tested - very long time
            120
        """
        if exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            g = (1 + RDF(5).sqrt()) / 2
            base_ring = RDF

        q12 = base_ring(1) / base_ring(2)
        z   = base_ring.zero()
        verts = [[s1*q12, s2*q12, s3*q12, s4*q12] for s1,s2,s3,s4 in itertools.product([1,-1], repeat=4)]
        V = (base_ring)**4
        verts.extend(V.basis())
        verts.extend(-v for v in V.basis())
        pts = [[s1 * q12, s2*g/2, s3/(2*g), z] for (s1,s2,s3) in itertools.product([1,-1], repeat=3)]
        for p in AlternatingGroup(4):
            verts.extend(p(x) for x in pts)
        return Polyhedron(vertices=verts, base_ring=base_ring)

    @rename_keyword(deprecation=18213, points_n='n', dim_n='dim')
    def cyclic_polytope(self, dim, n, base_ring=QQ):
        """
        Return a cyclic polytope.

        A cyclic polytope of dimension ``dim`` with ``n`` vertices is the convex
        hull of the points  ``(t,t^2,...,t^dim)`` with `t \in \{0,1,...,n-1\}` .
        For more information, see the :wikipedia:`Cyclic_polytope`.

        INPUT:

        - ``dim`` -- positive integer. the dimension of the polytope.

        - ``n`` -- positive integer. the number of vertices.

        - ``base_ring`` -- either ``QQ`` (default) or ``RDF``.

        EXAMPLES::

            sage: c = polytopes.cyclic_polytope(4,10)
            sage: c.f_vector()
            (1, 10, 45, 70, 35, 1)
        """
        verts = [[t**i for i in range(1,dim+1)] for t in range(n)]
        return Polyhedron(vertices=verts, base_ring=base_ring)

    @rename_keyword(deprecation=18213, dim_n='dim')
    def hypersimplex(self, dim, k, project=False):
        """
        Return the hypersimplex in dimension ``dim`` and parameter ``k``.

        The hypersimplex `\Delta_{d,k}` is the convex hull of the vertices made
        of `k` ones and `d-k` zeros. It lies in the `d-1` hyperplane of vectors
        of sum `k`. If you want a projected version to `\RR^{d-1}` (with
        floating point coordinates) then set ``project=True`` in the options.

        .. SEEALSO::

            :meth:`simplex`

        INPUT:

        - ``dim`` -- the dimension

        - ``n`` -- the numbers ``(1,...,n)`` are permuted

        - ``project`` -- (boolean, default ``False``) if ``True``, the polytope
          is (isometrically) projected to a vector space of dimension ``dim-1``.
          This operation turns the coordinates into floating point
          approximations and corresponds to the projection given by the matrix
          from :func:`zero_sum_projection`.

        EXAMPLES::

            sage: h_4_2 = polytopes.hypersimplex(4, 2)
            sage: h_4_2
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 6 vertices
            sage: h_4_2.f_vector()
            (1, 6, 12, 8, 1)
            sage: h_4_2.ehrhart_polynomial()    # optional - latte_int
            2/3*t^3 + 2*t^2 + 7/3*t + 1

            sage: h_7_3 = polytopes.hypersimplex(7, 3, project=True)
            sage: h_7_3
            A 6-dimensional polyhedron in RDF^6 defined as the convex hull of 35 vertices
            sage: h_7_3.f_vector()
            (1, 35, 210, 350, 245, 84, 14, 1)
        """
        verts = Permutations([0]*(dim-k) + [1]*k).list()
        if project: verts = project_points(*verts)
        return Polyhedron(vertices=verts)

    def permutahedron(self, n, project=False):
        """Return the standard permutahedron of (1,...,n)

        The permutahedron (or permutohedron) is the convex hull of the
        permutations of `\{1,\ldots,n\}` seen as vectors. The edges
        between the permutations correspond to multiplication on the
        right by an elementary transposition in the
        :class:`~sage.groups.perm_gps.permgroup_named.SymmetricGroup`.

        If we take the graph in which the vertices correspond to
        vertices of the polyhedron, and edges to edges, we get the
        :meth:`~sage.graphs.graph_generators.GraphGenerators.BubbleSortGraph`.

        INPUT:

        - ``n`` -- integer

        - ``project`` -- (boolean, default ``False``) if ``True``, the polytope
          is (isometrically) projected to a vector space of dimension ``dim-1``.
          This operation turns the coordinates into floating point
          approximations and corresponds to the projection given by the matrix
          from :func:`zero_sum_projection`.

        EXAMPLES::

            sage: perm4 = polytopes.permutahedron(4)
            sage: perm4
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 24 vertices
            sage: perm4.is_lattice_polytope()
            True
            sage: perm4.ehrhart_polynomial()   # optional - latte_int
            16*t^3 + 15*t^2 + 6*t + 1

            sage: perm4 = polytopes.permutahedron(4, project=True)
            sage: perm4
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 24 vertices
            sage: perm4.plot()
            Graphics3d Object
            sage: perm4.graph().is_isomorphic(graphs.BubbleSortGraph(4))
            True

        .. SEEALSO::

            * :meth:`~sage.graphs.graph_generators.GraphGenerators.BubbleSortGraph`
        """
        verts = list(itertools.permutations(range(1,n+1)))
        if project: verts = project_points(*verts)
        return Polyhedron(vertices=verts)

    def hypercube(self, dim):
        """
        Return a hypercube in the given dimension

        The `d` dimensional hypercube is the convex hull of the points `(\pm 1,
        \pm 1, \ldots, \pm 1)` in `\RR^d`. For more information see
        the :wikipedia:`Hypercube`.

        INPUT:

        - ``dim`` -- integer. The dimension of the cube.

        EXAMPLES::

            sage: four_cube = polytopes.hypercube(4)
            sage: four_cube.is_simple()
            True
            sage: four_cube.base_ring()
            Integer Ring
            sage: four_cube.volume()
            16
            sage: four_cube.ehrhart_polynomial()    # optional - latte_int
            16*t^4 + 32*t^3 + 24*t^2 + 8*t + 1
        """
        return Polyhedron(vertices = list(itertools.product([1,-1], repeat=dim)))

    def cube(self):
        r"""
        Return the cube.

        The cube is the Platonic solid that is obtained as the convex hull of
        the points `(\pm 1, \pm 1, \pm 1)`. It generalizes into several
        dimension into hypercubes.

        .. SEEALSO::

            :meth:`hypercube`

        EXAMPLES::

            sage: c = polytopes.cube()
            sage: c
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: c.f_vector()
            (1, 8, 12, 6, 1)
            sage: c.volume()
            8
            sage: c.plot()
            Graphics3d Object
        """
        return self.hypercube(3)

    n_cube = deprecated_function_alias(18213, hypercube)

    def cross_polytope(self, dim):
        """
        Return a cross-polytope in dimension ``dim``.

        A cross-polytope is a higher dimensional generalization of the
        octahedron. It is the convex hull of the `2d` points `(\pm 1, 0, \ldots,
        0)`, `(0, \pm 1, \ldots, 0)`, \ldots, `(0, 0, \ldots, \pm 1)`.
        See the :wikipedia:`Cross-polytope` for more information.

        INPUT:

        - ``dim`` -- integer. The dimension of the cross-polytope.

        EXAMPLES::

            sage: four_cross = polytopes.cross_polytope(4)
            sage: four_cross.f_vector()
            (1, 8, 24, 32, 16, 1)
            sage: four_cross.is_simple()
            False
        """
        verts = list((ZZ**dim).basis())
        verts.extend([-v for v in verts])
        return Polyhedron(vertices=verts)

    def parallelotope(self, generators):
        r"""
        Return the parallelotope spanned by the generators.

        The parallelotope is the multi-dimensional generalization of a
        parallelogram (2 generators) and a parallelepiped (3 generators).

        INPUT:

        - ``generators`` -- a list vector of vectors of same dimension

        EXAMPLES::

            sage: polytopes.parallelotope([ (1,0), (0,1) ])
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: polytopes.parallelotope([[1,2,3,4],[0,1,0,7],[3,1,0,2],[0,0,1,0]])
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 16 vertices

            sage: K = QuadraticField(2, 'sqrt2')
            sage: sqrt2 = K.gen()
            sage: polytopes.parallelotope([ (1,sqrt2), (1,-1) ])
            A 2-dimensional polyhedron in (Number Field in sqrt2 with defining
            polynomial x^2 - 2)^2 defined as the convex hull of 4 vertices
        """
        from sage.modules.free_module_element import vector
        from sage.structure.sequence import Sequence
        generators = map(vector,generators)
        V = Sequence(generators).universe()
        R = V.base_ring()

        from itertools import combinations
        par =  [ V.zero() ]
        par.extend(sum(c) for k in range(1,len(generators)+1) for c in combinations(generators,k))
        return Polyhedron(vertices=par, base_ring=R)

polytopes = Polytopes()
