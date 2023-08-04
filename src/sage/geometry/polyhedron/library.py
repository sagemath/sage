# -*- coding: utf-8 -*-
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
    :meth:`~sage.geometry.polyhedron.library.Polytopes.associahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.bitruncated_six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.buckyball`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cantellated_one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cantellated_six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cantitruncated_one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cantitruncated_six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cross_polytope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cube`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cuboctahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.cyclic_polytope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.flow_polytope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.Gosset_3_21`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.grand_antiprism`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.great_rhombicuboctahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.hypercube`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.hypersimplex`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.icosahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.icosidodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.Kirkman_icosahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.octahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.omnitruncated_one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.omnitruncated_six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.parallelotope`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.pentakis_dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.permutahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.generalized_permutahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.rectified_one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.rectified_six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.regular_polygon`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.rhombic_dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.rhombicosidodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.runcinated_one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.runcitruncated_one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.runcitruncated_six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.simplex`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.small_rhombicuboctahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.snub_cube`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.snub_dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.tetrahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.truncated_cube`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.truncated_dodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.truncated_icosidodecahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.truncated_tetrahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.truncated_octahedron`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.truncated_one_hundred_twenty_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.truncated_six_hundred_cell`
    :meth:`~sage.geometry.polyhedron.library.Polytopes.twenty_four_cell`
"""
########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#                     2011 Volker Braun <vbraun.name@gmail.com>
#                     2015 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2019 Jean-Philippe Labbé <labbe@math.fu-berlin.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

import itertools

from sage.rings.integer_ring import ZZ
from sage.misc.lazy_import import lazy_import
from sage.rings.rational_field import QQ
lazy_import('sage.combinat.permutation', 'Permutations')
lazy_import('sage.groups.perm_gps.permgroup_named', 'AlternatingGroup')
from .constructor import Polyhedron
from .parent import Polyhedra
lazy_import('sage.graphs.digraph', 'DiGraph')
lazy_import('sage.graphs.graph', 'Graph')
lazy_import('sage.combinat.root_system.associahedron', 'Associahedron')

def zero_sum_projection(d, base_ring=None):
    r"""
    Return a matrix corresponding to the projection on the orthogonal of
    `(1,1,\ldots,1)` in dimension `d`.

    The projection maps the orthonormal basis `(1,-1,0,\ldots,0) / \sqrt(2)`,
    `(1,1,-1,0,\ldots,0) / \sqrt(3)`, \ldots, `(1,1,\ldots,1,-1) / \sqrt(d)` to
    the canonical basis in `\RR^{d-1}`.

    OUTPUT:

    A matrix of dimensions `(d-1)\times d` defined over ``base_ring`` (default:
    :class:`RDF <sage.rings.real_double.RealDoubleField_class>`).

    EXAMPLES::

        sage: from sage.geometry.polyhedron.library import zero_sum_projection
        sage: zero_sum_projection(2)
        [ 0.7071067811865475 -0.7071067811865475]
        sage: zero_sum_projection(3)
        [ 0.7071067811865475 -0.7071067811865475                 0.0]
        [ 0.4082482904638631  0.4082482904638631 -0.8164965809277261]

    Exact computation in :class:`AA <sage.rings.qqbar.AlgebraicRealField>`::

        sage: zero_sum_projection(3, base_ring=AA)                                      # needs sage.rings.number_field
        [ 0.7071067811865475? -0.7071067811865475?                    0]
        [ 0.4082482904638630?  0.4082482904638630? -0.8164965809277260?]

    """
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
    if base_ring is None:
        from sage.rings.real_double import RDF as base_ring
    basis = [vector(base_ring, [1]*i + [-i] + [0]*(d-i-1)) for i in range(1, d)]
    return matrix(base_ring, [v / v.norm() for v in basis])


def project_points(*points, **kwds):
    """
    Projects a set of points into a vector space of dimension one less.

    INPUT:

    - ``points``... -- the points to project.

    - ``base_ring`` -- (defaults to ``RDF`` if keyword is ``None`` or not
      provided in ``kwds``) the base ring to use.

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

    These projections are compatible with the restriction. More precisely,
    given a vector `v`, the projection of `v` restricted to the first `i`
    coordinates will be equal to the projection of the first `i+1` coordinates
    of `v`::

        sage: project_points([1,2])    # abs tol 1e-15
        [(-0.7071067811865475)]
        sage: project_points([1,2,3])  # abs tol 1e-15
        [(-0.7071067811865475, -1.2247448713915892)]
        sage: project_points([1,2,3,4])     # abs tol 1e-15
        [(-0.7071067811865475, -1.2247448713915892, -1.7320508075688776)]
        sage: project_points([1,2,3,4,0])   # abs tol 1e-15
        [(-0.7071067811865475, -1.2247448713915892, -1.7320508075688776, 2.23606797749979)]

    Check that it is (almost) an isometry::

        sage: V = list(map(vector, IntegerVectors(n=5, length=3)))
        sage: P = project_points(*V)
        sage: for i in range(21):                                                       # needs sage.combinat
        ....:     for j in range(21):
        ....:         assert abs((V[i]-V[j]).norm() - (P[i]-P[j]).norm()) < 0.00001

    Example with exact computation::

        sage: V = [ vector(v) for v in IntegerVectors(n=4, length=2) ]
        sage: P = project_points(*V, base_ring=AA)                                      # needs sage.combinat sage.rings.number_field
        sage: for i in range(len(V)):                                                   # needs sage.combinat sage.rings.number_field
        ....:     for j in range(len(V)):
        ....:         assert (V[i]-V[j]).norm() == (P[i]-P[j]).norm()

    """
    if not points:
        return []
    base_ring = kwds.pop('base_ring', None)
    if base_ring is None:
        from sage.rings.real_double import RDF as base_ring
    from sage.modules.free_module_element import vector
    vecs = [vector(base_ring, p) for p in points]
    m = zero_sum_projection(len(vecs[0]), base_ring=base_ring)
    return [m * v for v in vecs]


def gale_transform_to_polytope(vectors, base_ring=None, backend=None):
    r"""
    Return the polytope associated to the list of vectors forming a Gale transform.

    This function is the inverse of
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.gale_transform`
    up to projective transformation.

    INPUT:

    - ``vectors`` -- the vectors of the Gale transform

    - ``base_ring`` -- string (default: `None`);
      the base ring to be used for the construction

    - ``backend`` -- string (default: `None`);
      the backend to use to create the polytope

    .. NOTE::

        The order of the input vectors will not be preserved.

        If the center of the (input) vectors is the origin,
        the function is much faster and might give a nicer representation
        of the polytope.

        If this is not the case, the vectors will be scaled
        (each by a positive scalar) accordingly to obtain the polytope.

    .. SEEALSO::

        :func`~sage.geometry.polyhedron.library.gale_transform_to_primal`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.library import gale_transform_to_polytope
        sage: points = polytopes.octahedron().gale_transform()
        sage: points
        ((0, -1), (-1, 0), (1, 1), (1, 1), (-1, 0), (0, -1))
        sage: P = gale_transform_to_polytope(points); P
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
        sage: P.vertices()
        (A vertex at (-1, 0, 0),
         A vertex at (0, -1, 0),
         A vertex at (0, 0, -1),
         A vertex at (0, 0, 1),
         A vertex at (0, 1, 0),
         A vertex at (1, 0, 0))

    One can specify the base ring::

        sage: gale_transform_to_polytope(
        ....:     [(1,1), (-1,-1), (1,0),
        ....:      (-1,0), (1,-1), (-2,1)]).vertices()
        (A vertex at (-25, 0, 0),
         A vertex at (-15, 50, -60),
         A vertex at (0, -25, 0),
         A vertex at (0, 0, -25),
         A vertex at (16, -35, 54),
         A vertex at (24, 10, 31))
        sage: gale_transform_to_polytope(
        ....:     [(1,1), (-1,-1), (1,0),
        ....:      (-1,0), (1,-1), (-2,1)],
        ....:     base_ring=RDF).vertices()
        (A vertex at (-0.64, 1.4, -2.16),
         A vertex at (-0.96, -0.4, -1.24),
         A vertex at (0.6, -2.0, 2.4),
         A vertex at (1.0, 0.0, 0.0),
         A vertex at (0.0, 1.0, 0.0),
         A vertex at (0.0, 0.0, 1.0))

    One can also specify the backend::

        sage: gale_transform_to_polytope(
        ....:     [(1,1), (-1,-1), (1,0),
        ....:      (-1,0), (1,-1), (-2,1)],
        ....:     backend='field').backend()
        'field'
        sage: gale_transform_to_polytope(
        ....:     [(1,1), (-1,-1), (1,0),
        ....:      (-1,0), (1,-1), (-2,1)],
        ....:     backend='cdd', base_ring=RDF).backend()
        'cdd'

    A gale transform corresponds to a polytope if and only if
    every oriented (linear) hyperplane
    has at least two vectors on each side.
    See Theorem 6.19 of [Zie2007]_.
    If this is not the case, one of two errors is raised.

    If there is such a hyperplane with no vector on one side,
    the vectors are not totally cyclic::

        sage: gale_transform_to_polytope([(0,1), (1,1), (1,0), (-1,0)])
        Traceback (most recent call last):
        ...
        ValueError: input vectors not totally cyclic

    If every hyperplane has at least one vector on each side, then the gale
    transform corresponds to a point configuration.
    It corresponds to a polytope if and only if this point configuration is
    convex and if and only if every hyperplane contains at least two vectors of
    the gale transform on each side.

    If this is not the case, an error is raised::

        sage: gale_transform_to_polytope([(0,1), (1,1), (1,0), (-1,-1)])
        Traceback (most recent call last):
        ...
        ValueError: the gale transform does not correspond to a polytope
    """
    vertices = gale_transform_to_primal(vectors, base_ring, backend)
    P = Polyhedron(vertices=vertices, base_ring=base_ring, backend=backend)

    if not P.n_vertices() == len(vertices):
        # If the input vectors are not totally cyclic, ``gale_transform_to_primal``
        # raises an error.
        # As no error was raised so far, the gale transform corresponds to
        # to a point configuration.
        # It corresponds to a polytope if and only if
        # ``vertices`` are in convex position.
        raise ValueError("the gale transform does not correspond to a polytope")

    return P

def gale_transform_to_primal(vectors, base_ring=None, backend=None):
    r"""
    Return a point configuration dual to a totally cyclic vector configuration.

    This is the dehomogenized vector configuration dual to the input.
    The dual vector configuration is acyclic and can therefore
    be dehomogenized as the input is totally cyclic.

    INPUT:

    - ``vectors`` -- the ordered vectors of the Gale transform

    - ``base_ring`` -- string (default: `None`);
      the base ring to be used for the construction

    - ``backend`` -- string (default: `None`);
      the backend to be use to construct a polyhedral,
      used internally in case the center is not the origin,
      see :func:`~sage.geometry.polyhedron.constructor.Polyhedron`

    OUTPUT: An ordered point configuration as list of vectors.

    .. NOTE::

        If the center of the (input) vectors is the origin,
        the function is much faster and might give a nicer representation
        of the point configuration.

        If this is not the case, the vectors will be scaled
        (each by a positive scalar) accordingly.

    ALGORITHM:

    Step 1: If the center of the (input) vectors is not the origin,
    we do an appropriate transformation to make it so.

    Step 2: We add a row of ones on top of ``Matrix(vectors)``.
    The right kernel of this larger matrix is the dual configuration space,
    and a basis of this space provides the dual point configuration.

    More concretely, the dual vector configuration (inhomogeneous)
    is obtained by taking a basis of the right kernel of ``Matrix(vectors)``.
    If the center of the (input) vectors is the origin,
    there exists a basis of the right kernel of the form
    ``[[1], [V]]``, where ``[1]`` represents a row of ones.
    Then, ``V`` is a dehomogenization and thus the dual point configuration.

    To extend ``[1]`` to a basis of ``Matrix(vectors)``, we add
    a row of ones to ``Matrix(vectors)`` and calculate a basis of the
    right kernel of the obtained matrix.

    REFERENCES:

        For more information, see Section 6.4 of [Zie2007]_
        or Definition 2.5.1 and Definition 4.1.35 of [DLRS2010]_.

    .. SEEALSO::

        :func`~sage.geometry.polyhedron.library.gale_transform_to_polytope`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.library import gale_transform_to_primal
        sage: points = ((0, -1), (-1, 0), (1, 1), (1, 1), (-1, 0), (0, -1))
        sage: gale_transform_to_primal(points)
        [(0, 0, 1), (0, 1, 0), (1, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]

    One can specify the base ring::

        sage: p = [(1,1), (-1,-1), (1,0), (-1,0), (1,-1), (-2,1)]
        sage: gtpp = gale_transform_to_primal(p); gtpp
        [(16, -35, 54),
         (24, 10, 31),
         (-15, 50, -60),
         (-25, 0, 0),
         (0, -25, 0),
         (0, 0, -25)]
        sage: (matrix(RDF, gtpp)/25 +
        ....:  matrix(gale_transform_to_primal(p, base_ring=RDF))).norm() < 1e-15
        True

    One can also specify the backend to be used internally::

        sage: gale_transform_to_primal(p, backend='field')
        [(48, -71, 88),
         (84, -28, 99),
         (-77, 154, -132),
         (-55, 0, 0),
         (0, -55, 0),
         (0, 0, -55)]
        sage: gale_transform_to_primal(p, backend='normaliz')           # optional - pynormaliz
        [(16, -35, 54),
         (24, 10, 31),
         (-15, 50, -60),
         (-25, 0, 0),
         (0, -25, 0),
         (0, 0, -25)]

    The input vectors should be totally cyclic::

        sage: gale_transform_to_primal([(0,1), (1,0), (1,1), (-1,0)])
        Traceback (most recent call last):
        ...
        ValueError: input vectors not totally cyclic

        sage: gale_transform_to_primal(
        ....:     [(1,1,0), (-1,-1,0), (1,0,0),
        ....:      (-1,0,0), (1,-1,0), (-2,1,0)], backend='field')
        Traceback (most recent call last):
        ...
        ValueError: input vectors not totally cyclic
    """
    from sage.modules.free_module_element import vector
    from sage.matrix.constructor import Matrix
    if base_ring:
        vectors = tuple(vector(base_ring, x) for x in vectors)
    else:
        vectors = tuple(vector(x) for x in vectors)

    if not sum(vectors).is_zero():
        # The center of the input vectors shall be the origin.
        # If this is not the case, we scale them accordingly.
        # This has the adventage that right kernel of ``vectors`` can be
        # presented in the form ``[[1], [V]]``, where ``V`` are the points
        # in the dual point configuration.
        # (Dehomogenization is straightforward.)

        # Scaling of the vectors is equivalent to finding a hyperplane that intersects
        # all vectors of the dual point configuration. But if the input is already provided
        # such that the vectors add up to zero, the coordinates might be nicer.
        # (And this is faster.)

        if base_ring:
            ker = Matrix(base_ring, vectors).left_kernel()
        else:
            ker = Matrix(vectors).left_kernel()
        solutions = Polyhedron(lines=tuple(ker.basis_matrix()), base_ring=base_ring, backend=backend)

        from sage.matrix.special import identity_matrix
        pos_orthant = Polyhedron(rays=identity_matrix(len(vectors)), base_ring=base_ring, backend=backend)
        pos_solutions = solutions.intersection(pos_orthant)
        if base_ring is ZZ:
            pos_solutions = pos_solutions.change_ring(ZZ)

        # Any integral point in ``pos_solutions`` will correspond to scaling-factors
        # that make ``sum(vectors)`` zero.
        x = pos_solutions.representative_point()
        if not all(y > 0 for y in x):
            raise ValueError("input vectors not totally cyclic")
        vectors = tuple(vec*x[i] for i,vec in enumerate(vectors))

    # The right kernel of ``vectors`` has a basis of the form ``[[1], [V]]``,
    # where ``V`` is the dehomogenized dual point configuration.
    # If we append a row of ones to ``vectors``, ``V`` is just the right kernel.
    if base_ring:
        m = Matrix(base_ring, vectors).transpose().stack(Matrix(base_ring, [[1]*len(vectors)]))
    else:
        m = Matrix(vectors).transpose().stack(Matrix([[1]*len(vectors)]))

    if m.rank() != len(vectors[0]) + 1:
        # The given vectors do not span the ambient space,
        # then there exists a nonnegative value vector.
        raise ValueError("input vectors not totally cyclic")

    return m.right_kernel_matrix(basis='computed').columns()


class Polytopes():
    """
    A class of constructors for commonly used, famous, or interesting
    polytopes.
    """

    def regular_polygon(self, n, exact=True, base_ring=None, backend=None):
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

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: # needs sage.rings.number_field
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

        TESTS::

            sage: # optional - pynormaliz, needs sage.rings.number_field
            sage: octagon = polytopes.regular_polygon(8, backend='normaliz')
            sage: octagon
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 8 vertices
            sage: octagon.n_vertices()
            8
            sage: octagon.volume()
            2*a
            sage: TestSuite(octagon).run()      # long time

            sage: TestSuite(polytopes.regular_polygon(5, exact=False)).run()
        """
        n = ZZ(n)
        if n <= 2:
            raise ValueError("n (={}) must be an integer greater than 2".format(n))

        if base_ring is None:
            if exact:
                from sage.rings.qqbar import AA as base_ring
            else:
                from sage.rings.real_double import RDF as base_ring

        try:
            omega = 2*base_ring.pi() / n
            verts = [((i*omega).sin(), (i*omega).cos()) for i in range(n)]
        except AttributeError:
            from sage.rings.qqbar import QQbar
            z = QQbar.zeta(n)
            verts = [(base_ring((z**k).imag()), base_ring((z**k).real())) for k in range(n)]

        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def Birkhoff_polytope(self, n, backend=None):
        """
        Return the Birkhoff polytope with `n!` vertices.

        The vertices of this polyhedron are the (flattened) `n` by `n`
        permutation matrices. So the ambient vector space has dimension `n^2`
        but the dimension of the polyhedron is `(n-1)^2`.

        INPUT:

        - ``n`` -- a positive integer giving the size of the permutation matrices.

        - ``backend`` -- the backend to use to create the polytope.

        .. SEEALSO::

            :meth:`sage.matrix.matrix2.Matrix.as_sum_of_permutations` -- return
            the current matrix as a sum of permutation matrices

        EXAMPLES::

            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: b3.f_vector()
            (1, 6, 15, 18, 9, 1)
            sage: b3.ambient_dim(), b3.dim()
            (9, 4)
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
            sage: b4.n_vertices(), b4.ambient_dim(), b4.dim()
            (24, 16, 9)

        TESTS::

            sage: b4norm = polytopes.Birkhoff_polytope(4,backend='normaliz')    # optional - pynormaliz
            sage: TestSuite(b4norm).run()                                       # optional - pynormaliz
            sage: TestSuite(polytopes.Birkhoff_polytope(3)).run()
        """
        from itertools import permutations
        verts = []
        for p in permutations(range(n)):
            verts.append([ZZ.one() if p[i] == j else ZZ.zero()
                          for j in range(n) for i in range(n)])
        return Polyhedron(vertices=verts, base_ring=ZZ, backend=backend)

    def simplex(self, dim=3, project=False, base_ring=None, backend=None):
        r"""
        Return the ``dim`` dimensional simplex.

        The `d`-simplex is the convex hull in `\RR^{d+1}` of the standard basis
        `(1,0,\ldots,0)`, `(0,1,\ldots,0)`, \ldots, `(0,0,\ldots,1)`. For more
        information, see the :wikipedia:`Simplex`.

        INPUT:

        - ``dim`` -- The dimension of the simplex, a positive
          integer.

        - ``project`` -- (boolean, default ``False``) if ``True``, the polytope
          is (isometrically) projected to a vector space of dimension
          ``dim-1``.  This corresponds to the projection given by the matrix
          from :func:`zero_sum_projection`.  By default, this operation turns
          the coordinates into floating point approximations (see
          ``base_ring``).

        - ``base_ring`` -- the base ring to use to create the polytope.
          If ``project`` is ``False``, this defaults to `\ZZ`.
          Otherwise, it defaults to ``RDF``.

        - ``backend`` -- the backend to use to create the polytope.

        .. SEEALSO::

            :meth:`tetrahedron`

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

        Computation in algebraic reals::

            sage: s3 = polytopes.simplex(3, project=True, base_ring=AA)                 # needs sage.rings.number_field
            sage: s3.volume() == sqrt(3+1) / factorial(3)                               # needs sage.rings.number_field
            True

        TESTS::

            sage: s6norm = polytopes.simplex(6,backend='normaliz')      # optional - pynormaliz
            sage: TestSuite(s6norm).run()                               # optional - pynormaliz
            sage: TestSuite(polytopes.simplex(5)).run()
        """
        verts = list((ZZ**(dim + 1)).basis())
        if project:
            # Handling of default in base_ring is delegated to project_points
            verts = project_points(*verts, base_ring=base_ring)
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def icosahedron(self, exact=True, base_ring=None, backend=None):
        r"""
        Return an icosahedron with edge length 1.

        The icosahedron is one of the Platonic solids. It has 20 faces
        and is dual to the :meth:`dodecahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- (optional) the ring in which the coordinates will
          belong to.  Note that this ring must contain `\sqrt(5)`. If it is not
          provided and ``exact=True`` it will be the number field
          `\QQ[\sqrt(5)]` and if ``exact=False`` it will be the real double
          field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: ico = polytopes.icosahedron()                                         # needs sage.rings.number_field
            sage: ico.f_vector()                                                        # needs sage.rings.number_field
            (1, 12, 30, 20, 1)
            sage: ico.volume()                                                          # needs sage.rings.number_field
            5/12*sqrt5 + 5/4

        Its non exact version::

            sage: ico = polytopes.icosahedron(exact=False)                              # needs sage.groups
            sage: ico.base_ring()                                                       # needs sage.groups
            Real Double Field
            sage: ico.volume()                  # known bug                             # needs sage.groups
            2.181694990...

        A version using `AA <sage.rings.qqbar.AlgebraicRealField>`::

            sage: ico = polytopes.icosahedron(base_ring=AA)     # long time             # needs sage.groups sage.rings.number_field
            sage: ico.base_ring()                               # long time             # needs sage.groups sage.rings.number_field
            Algebraic Real Field
            sage: ico.volume()                                  # long time             # needs sage.groups sage.rings.number_field
            2.181694990624913?

        Note that if base ring is provided it must contain the square root of
        `5`. Otherwise you will get an error::

            sage: polytopes.icosahedron(base_ring=QQ)                                   # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 1/4*sqrt(5) + 1/4 to a rational

        TESTS::

            sage: # optional - pynormaliz, needs sage.groups sage.rings.number_field
            sage: ico = polytopes.icosahedron(backend='normaliz')
            sage: ico.f_vector()
            (1, 12, 30, 20, 1)
            sage: ico.volume()
            5/12*sqrt5 + 5/4
            sage: TestSuite(ico).run()

            sage: ico = polytopes.icosahedron(exact=False)                              # needs sage.groups
            sage: TestSuite(ico).run(skip="_test_lawrence")                             # needs sage.groups

        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            g = (1 + base_ring(5).sqrt()) / 2

        r12 = base_ring.one() / 2
        z = base_ring.zero()
        pts = [[z, s1 * r12, s2 * g / 2]
               for s1, s2 in itertools.product([1, -1], repeat=2)]
        verts = [p(v) for p in AlternatingGroup(3) for v in pts]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def dodecahedron(self, exact=True, base_ring=None, backend=None):
        r"""
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

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: # needs sage.groups sage.rings.number_field
            sage: d12 = polytopes.dodecahedron()
            sage: d12.f_vector()
            (1, 20, 30, 12, 1)
            sage: d12.volume()
            -176*sqrt5 + 400
            sage: numerical_approx(_)
            6.45203596003699

            sage: d12 = polytopes.dodecahedron(exact=False)                             # needs sage.groups
            sage: d12.base_ring()                                                       # needs sage.groups
            Real Double Field

        Here is an error with a field that does not contain `\sqrt(5)`::

            sage: polytopes.dodecahedron(base_ring=QQ)                                  # needs sage.groups sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 1/4*sqrt(5) + 1/4 to a rational

        TESTS::

            sage: d12 = polytopes.dodecahedron(backend='normaliz')      # optional - pynormaliz, needs sage.groups sage.rings.number_field
            sage: d12.f_vector()                                        # optional - pynormaliz, needs sage.groups sage.rings.number_field
            (1, 20, 30, 12, 1)
            sage: TestSuite(d12).run()                                  # optional - pynormaliz, needs sage.groups sage.rings.number_field
        """
        return self.icosahedron(exact=exact, base_ring=base_ring, backend=backend).polar()

    def small_rhombicuboctahedron(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the (small) rhombicuboctahedron.

        The rhombicuboctahedron is an Archimedean solid with 24 vertices and 26
        faces. See the :wikipedia:`Rhombicuboctahedron` for more information.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False``
          it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: sr = polytopes.small_rhombicuboctahedron()                            # needs sage.rings.number_field
            sage: sr.f_vector()                                                         # needs sage.rings.number_field
            (1, 24, 48, 26, 1)
            sage: sr.volume()                                                           # needs sage.rings.number_field
            80/3*sqrt2 + 32

        The faces are `8` equilateral triangles and `18` squares::

            sage: sum(1 for f in sr.facets() if len(f.vertices()) == 3)                 # needs sage.rings.number_field
            8
            sage: sum(1 for f in sr.facets() if len(f.vertices()) == 4)                 # needs sage.rings.number_field
            18

        Its non exact version::

            sage: sr = polytopes.small_rhombicuboctahedron(False)
            sage: sr
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of
            24 vertices
            sage: sr.f_vector()
            (1, 24, 48, 26, 1)

        TESTS::

            sage: # optional - pynormaliz, needs sage.rings.number_field
            sage: sr = polytopes.small_rhombicuboctahedron(backend='normaliz')
            sage: sr.f_vector()
            (1, 24, 48, 26, 1)
            sage: sr.volume()
            80/3*sqrt2 + 32
            sage: TestSuite(sr).run()           # long time
        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(2, 'sqrt2')
            sqrt2 = K.gen()
            base_ring = K
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            sqrt2 = base_ring(2).sqrt()

        one = base_ring.one()
        a = sqrt2 + one
        verts = []
        verts.extend([s1*one, s2*one, s3*a] for s1, s2, s3 in itertools.product([1, -1], repeat=3))
        verts.extend([s1*one, s3*a, s2*one] for s1, s2, s3 in itertools.product([1, -1], repeat=3))
        verts.extend([s1*a, s2*one, s3*one] for s1, s2, s3 in itertools.product([1, -1], repeat=3))
        return Polyhedron(vertices=verts, backend=backend)

    def great_rhombicuboctahedron(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the great rhombicuboctahedron.

        The great rhombicuboctahedron (or truncated cuboctahedron) is an
        Archimedean solid with 48 vertices and 26 faces. For more information
        see the :wikipedia:`Truncated_cuboctahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False``
          it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: gr = polytopes.great_rhombicuboctahedron()    # long time             # needs sage.rings.number_field
            sage: gr.f_vector()                 # long time                             # needs sage.rings.number_field
            (1, 48, 72, 26, 1)

        A faster implementation is obtained by setting ``exact=False``::

            sage: gr = polytopes.great_rhombicuboctahedron(exact=False)
            sage: gr.f_vector()
            (1, 48, 72, 26, 1)

        Its facets are 4 squares, 8 regular hexagons and 6 regular octagons::

            sage: sum(1 for f in gr.facets() if len(f.vertices()) == 4)
            12
            sage: sum(1 for f in gr.facets() if len(f.vertices()) == 6)
            8
            sage: sum(1 for f in gr.facets() if len(f.vertices()) == 8)
            6
        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            base_ring = QuadraticField(2, 'sqrt2')
            sqrt2 = base_ring.gen()
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            sqrt2 = base_ring(2).sqrt()

        one = base_ring.one()
        v1 = sqrt2 + 1
        v2 = 2 * sqrt2 + 1
        verts = [[s1 * z1, s2 * z2, s3 * z3]
                 for z1, z2, z3 in itertools.permutations([one, v1, v2])
                 for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def rhombic_dodecahedron(self, backend=None):
        """
        Return the rhombic dodecahedron.

        The rhombic dodecahedron is a polytope dual to the cuboctahedron. It
        has 14 vertices and 12 faces. For more information see
        the :wikipedia:`Rhombic_dodecahedron`.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        .. SEEALSO::

            :meth:`cuboctahedron`

        EXAMPLES::

            sage: rd = polytopes.rhombic_dodecahedron()
            sage: rd.f_vector()
            (1, 14, 24, 12, 1)

        Its facets are 12 quadrilaterals (not all identical)::

            sage: sum(1 for f in rd.facets() if len(f.vertices()) == 4)
            12

        Some more computations::

            sage: p = rd.ehrhart_polynomial()    # optional - latte_int
            sage: p                              # optional - latte_int
            16*t^3 + 12*t^2 + 4*t + 1
            sage: [p(i) for i in [1,2,3,4]]      # optional - latte_int
            [33, 185, 553, 1233]
            sage: [len((i*rd).integral_points()) for i in [1,2,3,4]]
            [33, 185, 553, 1233]

        TESTS::

            sage: rd_norm = polytopes.rhombic_dodecahedron(backend='normaliz')  # optional - pynormaliz
            sage: TestSuite(rd_norm).run()                                      # optional - pynormaliz
        """
        v = [[2,0,0],[-2,0,0],[0,2,0],[0,-2,0],[0,0,2],[0,0,-2]]
        v.extend((itertools.product([1, -1], repeat=3)))
        return Polyhedron(vertices=v, base_ring=ZZ, backend=backend)

    def cuboctahedron(self, backend=None):
        r"""
        Return the cuboctahedron.

        The cuboctahedron is an Archimedean solid with 12 vertices and 14 faces
        dual to the rhombic dodecahedron. It can be defined as the convex hull
        of the twelve vertices `(0, \pm 1, \pm 1)`, `(\pm 1, 0, \pm 1)` and
        `(\pm 1, \pm 1, 0)`. For more information, see the
        :wikipedia:`Cuboctahedron`.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        .. SEEALSO::

            :meth:`rhombic_dodecahedron`

        EXAMPLES::

            sage: co = polytopes.cuboctahedron()
            sage: co.f_vector()
            (1, 12, 24, 14, 1)

        Its facets are 8 triangles and 6 squares::

            sage: sum(1 for f in co.facets() if len(f.vertices()) == 3)
            8
            sage: sum(1 for f in co.facets() if len(f.vertices()) == 4)
            6

        Some more computation::

            sage: co.volume()
            20/3
            sage: co.ehrhart_polynomial()      # optional - latte_int
            20/3*t^3 + 8*t^2 + 10/3*t + 1

        TESTS::

            sage: co_norm = polytopes.cuboctahedron(backend='normaliz')      # optional - pynormaliz
            sage: TestSuite(co_norm).run()                                   # optional - pynormaliz
        """
        v = [[0, -1, -1], [0, 1, -1], [0, -1, 1], [0, 1, 1],
             [-1, -1, 0], [1, -1, 0], [-1, 1, 0], [1, 1, 0],
             [-1, 0, -1], [1, 0, -1], [-1, 0, 1], [1, 0, 1]]
        return Polyhedron(vertices=v, base_ring=ZZ, backend=backend)

    def truncated_cube(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the truncated cube.

        The truncated cube is an Archimedean solid with 24 vertices
        and 14 faces. It can be defined as the convex hull of the 24 vertices
        `(\pm x, \pm 1, \pm 1), (\pm 1, \pm x, \pm 1), (\pm 1, \pm 1, \pm x)`
        where `x = \sqrt(2) - 1`. For more information, see the
        :wikipedia:`Truncated_cube`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\sqrt{2}]` and if ``exact=False`` it
          will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: co = polytopes.truncated_cube()                                       # needs sage.rings.number_field
            sage: co.f_vector()                                                         # needs sage.rings.number_field
            (1, 24, 36, 14, 1)

        Its facets are 8 triangles and 6 octogons::

            sage: sum(1 for f in co.facets() if len(f.vertices()) == 3)                 # needs sage.rings.number_field
            8
            sage: sum(1 for f in co.facets() if len(f.vertices()) == 8)                 # needs sage.rings.number_field
            6

        Some more computation::

            sage: co.volume()                                                           # needs sage.rings.number_field
            56/3*sqrt2 - 56/3

        TESTS::

            sage: co = polytopes.truncated_cube(backend='normaliz')     # optional - pynormaliz, needs sage.rings.number_field
            sage: co.f_vector()                 # optional - pynormaliz                 # needs sage.rings.number_field
            (1, 24, 36, 14, 1)
            sage: TestSuite(co).run()           # optional - pynormaliz                 # needs sage.rings.number_field

        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(2, 'sqrt2')
            sqrt2 = K.gen()
            g = sqrt2 - 1
            base_ring = K
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            g = base_ring(2).sqrt() - 1

        v = [[a * g, b, c] for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
        v += [[a, b * g, c] for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
        v += [[a, b, c * g] for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
        return Polyhedron(vertices=v, base_ring=base_ring, backend=backend)

    def tetrahedron(self, backend=None):
        """
        Return the tetrahedron.

        The tetrahedron is a Platonic solid with 4 vertices and 4 faces
        dual to itself. It can be defined as the convex hull
        of the 4 vertices `(0, 0, 0)`, `(1, 1, 0)`, `(1, 0, 1)` and
        `(0, 1, 1)`. For more information, see the
        :wikipedia:`Tetrahedron`.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        .. SEEALSO::

            :meth:`simplex`

        EXAMPLES::

            sage: co = polytopes.tetrahedron()
            sage: co.f_vector()
            (1, 4, 6, 4, 1)

        Its facets are 4 triangles::

            sage: sum(1 for f in co.facets() if len(f.vertices()) == 3)
            4

        Some more computation::

            sage: co.volume()
            1/3
            sage: co.ehrhart_polynomial()      # optional - latte_int
            1/3*t^3 + t^2 + 5/3*t + 1

        TESTS::

            sage: t_norm = polytopes.tetrahedron(backend='normaliz')    # optional - pynormaliz
            sage: TestSuite(t_norm).run()                               # optional - pynormaliz
        """
        v = [[0, 0, 0], [1, 0, 1], [1, 1, 0], [0, 1, 1]]
        return Polyhedron(vertices=v, base_ring=ZZ, backend=backend)

    def truncated_tetrahedron(self, backend=None):
        r"""
        Return the truncated tetrahedron.

        The truncated tetrahedron is an Archimedean solid with 12
        vertices and 8 faces. It can be defined as the convex hull off
        all the permutations of `(\pm 1, \pm 1, \pm 3)` with an even
        number of minus signs. For more information, see the
        :wikipedia:`Truncated_tetrahedron`.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: co = polytopes.truncated_tetrahedron()
            sage: co.f_vector()
            (1, 12, 18, 8, 1)

        Its facets are 4 triangles and 4 hexagons::

            sage: sum(1 for f in co.facets() if len(f.vertices()) == 3)
            4
            sage: sum(1 for f in co.facets() if len(f.vertices()) == 6)
            4

        Some more computation::

            sage: co.volume()
            184/3
            sage: co.ehrhart_polynomial()      # optional - latte_int
            184/3*t^3 + 28*t^2 + 26/3*t + 1

        TESTS::

            sage: tt_norm = polytopes.truncated_tetrahedron(backend='normaliz')     # optional - pynormaliz
            sage: TestSuite(tt_norm).run()                                          # optional - pynormaliz
        """
        v = [(3,1,1), (1,3,1), (1,1,3),
             (-3,-1,1), (-1,-3,1), (-1,-1,3),
             (-3,1,-1), (-1,3,-1), (-1,1,-3),
             (3,-1,-1), (1,-3,-1), (1,-1,-3)]
        return Polyhedron(vertices=v, base_ring=ZZ, backend=backend)

    def truncated_octahedron(self, backend=None):
        r"""
        Return the truncated octahedron.

        The truncated octahedron is an Archimedean solid with 24
        vertices and 14 faces. It can be defined as the convex hull
        off all the permutations of `(0, \pm 1, \pm 2)`. For more
        information, see the :wikipedia:`Truncated_octahedron`.

        This is also known as the permutohedron of dimension 3.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: co = polytopes.truncated_octahedron()
            sage: co.f_vector()
            (1, 24, 36, 14, 1)

        Its facets are 6 squares and 8 hexagons::

            sage: sum(1 for f in co.facets() if len(f.vertices()) == 4)
            6
            sage: sum(1 for f in co.facets() if len(f.vertices()) == 6)
            8

        Some more computation::

            sage: co.volume()
            32
            sage: co.ehrhart_polynomial()       # optional - latte_int                  # needs sage.combinat
            32*t^3 + 18*t^2 + 6*t + 1

        TESTS::

            sage: to_norm = polytopes.truncated_octahedron(backend='normaliz')  # optional - pynormaliz, needs sage.combinat
            sage: TestSuite(to_norm).run()                                      # optional - pynormaliz, needs sage.combinat
        """
        v = [(0, e, f) for e in [-1, 1] for f in [-2, 2]]
        v = [(xyz[sigma(1) - 1], xyz[sigma(2) - 1], xyz[sigma(3) - 1])
             for sigma in Permutations(3) for xyz in v]
        return Polyhedron(vertices=v, base_ring=ZZ, backend=backend)

    def octahedron(self, backend=None):
        r"""
        Return the octahedron.

        The octahedron is a Platonic solid with 6 vertices and 8 faces
        dual to the cube. It can be defined as the convex hull
        of the six vertices `(0, 0, \pm 1)`, `(\pm 1, 0, 0)` and
        `(0, \pm 1, 0)`. For more information, see the
        :wikipedia:`Octahedron`.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: co = polytopes.octahedron()
            sage: co.f_vector()
            (1, 6, 12, 8, 1)

        Its facets are 8 triangles::

            sage: sum(1 for f in co.facets() if len(f.vertices()) == 3)
            8

        Some more computation::

            sage: co.volume()
            4/3
            sage: co.ehrhart_polynomial()      # optional - latte_int
            4/3*t^3 + 2*t^2 + 8/3*t + 1

        TESTS::

            sage: o_norm = polytopes.octahedron(backend='normaliz')     # optional - pynormaliz
            sage: TestSuite(o_norm).run()                               # optional - pynormaliz
        """
        v = [[0, 0, -1], [0, 0, 1], [1, 0, 0],
             [-1, 0, 0], [0, 1, 0], [0, -1, 0]]
        return Polyhedron(vertices=v, base_ring=ZZ, backend=backend)

    def snub_cube(self, exact=False, base_ring=None, backend=None, verbose=False):
        """
        Return a snub cube.

        The snub cube is an Archimedean solid. It has 24 vertices and 38 faces.
        For more information see the :wikipedia:`Snub_cube`.

        The constant `z` used in constructing this polytope is the reciprocal
        of the tribonacci constant, that is, the solution of the equation
        `x^3 + x^2 + x - 1 = 0`.
        See :wikipedia:`Generalizations_of_Fibonacci_numbers#Tribonacci_numbers`.

        INPUT:

        - ``exact`` -- (boolean, default ``False``) if ``True`` use exact
          coordinates instead of floating point approximations

        - ``base_ring`` -- the field to use. If ``None`` (the default),
          construct the exact number field needed (if ``exact`` is ``True``) or
          default to ``RDF`` (if ``exact`` is ``True``).

        - ``backend`` -- the backend to use to create the polytope.  If
          ``None`` (the default), the backend will be selected automatically.

        EXAMPLES::

            sage: # needs sage.groups
            sage: sc_inexact = polytopes.snub_cube(exact=False); sc_inexact
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 24 vertices
            sage: sc_inexact.f_vector()
            (1, 24, 60, 38, 1)

            sage: # long time, needs sage.groups sage.rings.number_field
            sage: sc_exact = polytopes.snub_cube(exact=True)
            sage: sc_exact.f_vector()
            (1, 24, 60, 38, 1)
            sage: sorted(sc_exact.vertices())
            [A vertex at (-1, -z, -z^2),
             A vertex at (-1, -z^2, z),
             A vertex at (-1, z^2, -z),
             A vertex at (-1, z, z^2),
             A vertex at (-z, -1, z^2),
             A vertex at (-z, -z^2, -1),
             A vertex at (-z, z^2, 1),
             A vertex at (-z, 1, -z^2),
             A vertex at (-z^2, -1, -z),
             A vertex at (-z^2, -z, 1),
             A vertex at (-z^2, z, -1),
             A vertex at (-z^2, 1, z),
             A vertex at (z^2, -1, z),
             A vertex at (z^2, -z, -1),
             A vertex at (z^2, z, 1),
             A vertex at (z^2, 1, -z),
             A vertex at (z, -1, -z^2),
             A vertex at (z, -z^2, 1),
             A vertex at (z, z^2, -1),
             A vertex at (z, 1, z^2),
             A vertex at (1, -z, z^2),
             A vertex at (1, -z^2, -z),
             A vertex at (1, z^2, z),
             A vertex at (1, z, -z^2)]
            sage: sc_exact.is_combinatorially_isomorphic(sc_inexact)
            True

        TESTS::

            sage: sc = polytopes.snub_cube(exact=True, backend='normaliz')      # optional - pynormaliz, needs sage.groups sage.rings.number_field
            sage: sc.f_vector()                                                 # optional - pynormaliz, needs sage.groups sage.rings.number_field
            (1, 24, 60, 38, 1)

        """
        def construct_z(field):
            # z here is the reciprocal of the tribonacci constant, that is, the
            # solution of the equation x^3 + x^2 + x - 1 = 0.
            tsqr33 = 3 * field(33).sqrt()
            return ((17 + tsqr33)**QQ((1, 3)) - (-17 + tsqr33)**QQ((1, 3)) - 1) / 3

        if exact and base_ring is None:
            # construct the exact number field
            from sage.rings.qqbar import AA
            from sage.rings.number_field.number_field import NumberField
            R = QQ['x']
            f = R([-1, 1, 1, 1])
            embedding = construct_z(AA)
            base_ring = NumberField(f, name='z', embedding=embedding)
            z = base_ring.gen()
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            z = construct_z(base_ring)

        verts = []
        z2 = z ** 2
        A3 = AlternatingGroup(3)
        for e in [-1, 1]:
            for f in [-1, 1]:
                for g in [-1, 1]:
                    if e * f * g == -1:
                        v = [e, f * z, g * z2]
                        for p in A3:
                            verts += [p(v)]
                    else:
                        v = [f * z, e, g * z2]
                        for p in A3:
                            verts += [p(v)]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def buckyball(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the bucky ball.

        The bucky ball, also known as the truncated icosahedron is an
        Archimedean solid.  It has 32 faces and 60 vertices.

        .. SEEALSO::

            :meth:`icosahedron`

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False``
          it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: bb = polytopes.buckyball()    # long time                             # needs sage.groups sage.rings.number_field
            sage: bb.f_vector()                 # long time                             # needs sage.groups sage.rings.number_field
            (1, 60, 90, 32, 1)
            sage: bb.base_ring()                # long time                             # needs sage.groups sage.rings.number_field
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        A much faster implementation using floating point approximations::

            sage: bb = polytopes.buckyball(exact=False)                                 # needs sage.groups
            sage: bb.f_vector()                                                         # needs sage.groups
            (1, 60, 90, 32, 1)
            sage: bb.base_ring()                                                        # needs sage.groups
            Real Double Field

        Its facets are 5 regular pentagons and 6 regular hexagons::

            sage: sum(1 for f in bb.facets() if len(f.vertices()) == 5)                 # needs sage.groups
            12
            sage: sum(1 for f in bb.facets() if len(f.vertices()) == 6)                 # needs sage.groups
            20

        TESTS::

            sage: # optional - pynormaliz, needs sage.groups sage.rings.number_field
            sage: bb = polytopes.buckyball(backend='normaliz')
            sage: bb.f_vector()
            (1, 60, 90, 32, 1)
            sage: bb.base_ring()
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        """
        return self.icosahedron(exact=exact, base_ring=base_ring, backend=backend).truncation()

    def icosidodecahedron(self, exact=True, backend=None):
        """
        Return the icosidodecahedron.

        The Icosidodecahedron is a polyhedron with twenty triangular faces and
        twelve pentagonal faces. For more information see the
        :wikipedia:`Icosidodecahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: id = polytopes.icosidodecahedron()                                    # needs sage.groups sage.rings.number_field
            sage: id.f_vector()                                                         # needs sage.groups sage.rings.number_field
            (1, 30, 60, 32, 1)

        TESTS::

            sage: id = polytopes.icosidodecahedron(exact=False); id                     # needs sage.groups sage.rings.number_field
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 30 vertices
            sage: TestSuite(id).run(skip=["_test_is_combinatorially_isomorphic",        # needs sage.groups sage.rings.number_field
            ....:                         "_test_product",
            ....:                         "_test_pyramid",
            ....:                         "_test_lawrence"])

            sage: # optional - pynormaliz, needs sage.groups sage.rings.number_field
            sage: id = polytopes.icosidodecahedron(backend='normaliz')
            sage: id.f_vector()
            (1, 30, 60, 32, 1)
            sage: id.base_ring()
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?
            sage: TestSuite(id).run()           # long time
        """
        from sage.rings.number_field.number_field import QuadraticField
        from itertools import product

        K = QuadraticField(5, 'sqrt5')
        one = K.one()
        phi = (one+K.gen())/2

        gens = [((-1)**a*one/2, (-1)**b*phi/2, (-1)**c*(one+phi)/2)
                for a, b, c in product([0, 1], repeat=3)]
        gens.extend([(0, 0, phi), (0, 0, -phi)])

        verts = []
        for p in AlternatingGroup(3):
            verts.extend(p(x) for x in gens)

        if exact:
            return Polyhedron(vertices=verts, base_ring=K, backend=backend)
        else:
            from sage.rings.real_mpfr import RR
            verts = [(RR(x), RR(y), RR(z)) for x, y, z in verts]
            return Polyhedron(vertices=verts, backend=backend)

    def icosidodecahedron_V2(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the icosidodecahedron.

        The icosidodecahedron is an Archimedean solid.
        It has 32 faces and 30 vertices. For more information, see the
        :wikipedia:`Icosidodecahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to.
          If it is not provided and ``exact=True`` it will be a the number
          field `\QQ[\phi]` where `\phi` is the golden ratio and if
          ``exact=False`` it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: id = polytopes.icosidodecahedron_V2()   # long time - 6secs
            sage: id.f_vector()                           # long time
            (1, 30, 60, 32, 1)
            sage: id.base_ring()                          # long time
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        A much faster implementation using floating point approximations::

            sage: id = polytopes.icosidodecahedron_V2(exact=False)
            sage: id.f_vector()
            (1, 30, 60, 32, 1)
            sage: id.base_ring()
            Real Double Field

        Its facets are 20 triangles and 12 regular pentagons::

            sage: sum(1 for f in id.facets() if len(f.vertices()) == 3)
            20
            sage: sum(1 for f in id.facets() if len(f.vertices()) == 5)
            12

        TESTS::

            sage: # optional - pynormaliz
            sage: id = polytopes.icosidodecahedron_V2(backend='normaliz')
            sage: id.f_vector()
            (1, 30, 60, 32, 1)
            sage: id.base_ring()
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?
            sage: TestSuite(id).run()           # long time
        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            g = (1 + base_ring(5).sqrt()) / 2

        pts = [[g, 0, 0], [-g, 0, 0]]
        pts += [[s1 * base_ring.one() / 2, s2 * g / 2, s3 * (1 + g)/2]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        verts = pts
        verts += [[v[1], v[2], v[0]] for v in pts]
        verts += [[v[2], v[0], v[1]] for v in pts]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def truncated_dodecahedron(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the truncated dodecahedron.

        The truncated dodecahedron is an Archimedean solid.
        It has 32 faces and 60 vertices. For more information, see the
        :wikipedia:`Truncated dodecahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False``
          it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: td = polytopes.truncated_dodecahedron()                               # needs sage.rings.number_field
            sage: td.f_vector()                                                         # needs sage.rings.number_field
            (1, 60, 90, 32, 1)
            sage: td.base_ring()                                                        # needs sage.rings.number_field
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        Its facets are 20 triangles and 12 regular decagons::

            sage: sum(1 for f in td.facets() if len(f.vertices()) == 3)                 # needs sage.rings.number_field
            20
            sage: sum(1 for f in td.facets() if len(f.vertices()) == 10)                # needs sage.rings.number_field
            12

        The faster implementation using floating point approximations does not
        fully work unfortunately, see https://github.com/cddlib/cddlib/pull/7
        for a detailed discussion of this case::

            sage: td = polytopes.truncated_dodecahedron(exact=False)  # random
            doctest:warning
            ...
            UserWarning: This polyhedron data is numerically complicated; cdd
            could not convert between the inexact V and H representation
            without loss of data. The resulting object might show
            inconsistencies.
            sage: td.f_vector()
            Traceback (most recent call last):
            ...
            ValueError: not all vertices are intersections of facets
            sage: td.base_ring()
            Real Double Field

        TESTS::

            sage: # optional - pynormaliz, needs sage.rings.number_field
            sage: td = polytopes.truncated_dodecahedron(backend='normaliz')
            sage: td.f_vector()
            (1, 60, 90, 32, 1)
            sage: td.base_ring()
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            g = (1 + base_ring(5).sqrt()) / 2

        z = base_ring.zero()
        pts = [[z, s1 * base_ring.one() / g, s2 * (2 + g)]
                for s1, s2 in itertools.product([1, -1], repeat=2)]
        pts += [[s1 * base_ring.one() / g, s2 * g, s3 * (2 * g)]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        pts += [[s1 * g, s2 * base_ring(2), s3 * (g ** 2)]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        verts = pts
        verts += [[v[1], v[2], v[0]] for v in pts]
        verts += [[v[2], v[0], v[1]] for v in pts]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def pentakis_dodecahedron(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the pentakis dodecahedron.

        The pentakis dodecahedron (orkisdodecahedron) is a face-regular,
        vertex-uniform polytope dual to the truncated icosahedron.  It has 60
        facets and 32 vertices. See the :wikipedia:`Pentakis_dodecahedron` for more
        information.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False``
          it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: pd = polytopes.pentakis_dodecahedron()    # long time - ~10 sec
            sage: pd.n_vertices()                           # long time
            32
            sage: pd.n_inequalities()                       # long time
            60

        A much faster implementation is obtained when setting ``exact=False``::

            sage: pd = polytopes.pentakis_dodecahedron(exact=False)                     # needs sage.groups
            sage: pd.n_vertices()                                                       # needs sage.groups
            32
            sage: pd.n_inequalities()                                                   # needs sage.groups
            60

        The 60 are triangles::

            sage: all(len(f.vertices()) == 3 for f in pd.facets())                      # needs sage.groups
            True
        """
        return self.buckyball(exact=exact, base_ring=base_ring, backend=backend).polar()

    def Kirkman_icosahedron(self, backend=None):
        r"""
        Return the Kirkman icosahedron.

        The Kirkman icosahedron is a 3-polytope with integer coordinates: `(\pm
        9, \pm 6, \pm 6)`, `(\pm 12, \pm 4, 0)`, `(0, \pm 12, \pm 8)`, `(\pm 6,
        0, \pm 12)`. See [Fe2012]_ for more information.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: ki = polytopes.Kirkman_icosahedron()
            sage: ki.f_vector()
            (1, 20, 38, 20, 1)

            sage: ki.volume()
            6528

            sage: vertices = ki.vertices()
            sage: edges = [[vector(edge[0]),vector(edge[1])] for edge in ki.bounded_edges()]
            sage: edge_lengths = [norm(edge[0]-edge[1]) for edge in edges]
            sage: sorted(set(edge_lengths))
            [7, 8, 9, 11, 12, 14, 16]

        TESTS::

            sage: ki_norm = polytopes.Kirkman_icosahedron(backend='normaliz')   # optional - pynormaliz
            sage: TestSuite(ki_norm).run()                                      # optional - pynormaliz
        """
        vertices = [[9, 6, 6], [-9, 6, 6], [9, -6, 6], [9, 6, -6],
                    [-9, -6, 6], [-9, 6, -6], [9, -6, -6], [-9, -6, -6],
                    [12, 4, 0], [-12, 4, 0], [12, -4, 0], [-12, -4, 0],
                    [0, 12, 8], [0, -12, 8], [0, 12, -8], [0, -12, -8],
                    [6, 0, 12], [-6, 0, 12], [6, 0, -12], [-6, 0, -12]]
        return Polyhedron(vertices=vertices, base_ring=ZZ, backend=backend)

    def rhombicosidodecahedron(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the rhombicosidodecahedron.

        The rhombicosidodecahedron is an Archimedean solid.
        It has 62 faces and 60 vertices. For more information, see the
        :wikipedia:`Rhombicosidodecahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False``
          it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: rid = polytopes.rhombicosidodecahedron()   # long time (6secs)
            sage: rid.f_vector()                             # long time
            (1, 60, 120, 62, 1)
            sage: rid.base_ring()                            # long time
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        A much faster implementation using floating point approximations::

            sage: rid = polytopes.rhombicosidodecahedron(exact=False)
            sage: rid.f_vector()
            (1, 60, 120, 62, 1)
            sage: rid.base_ring()
            Real Double Field

        Its facets are 20 triangles, 30 squares and 12 pentagons::

            sage: sum(1 for f in rid.facets() if len(f.vertices()) == 3)
            20
            sage: sum(1 for f in rid.facets() if len(f.vertices()) == 4)
            30
            sage: sum(1 for f in rid.facets() if len(f.vertices()) == 5)
            12

        TESTS::

            sage: # optional - pynormaliz
            sage: rid = polytopes.rhombicosidodecahedron(backend='normaliz')
            sage: rid.f_vector()
            (1, 60, 120, 62, 1)
            sage: rid.base_ring()
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            g = (1 + base_ring(5).sqrt()) / 2

        pts = [[s1 * base_ring.one(), s2 * base_ring.one(), s3 * (g**3)]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        pts += [[s1 * (g**2), s2 * g, s3 * 2 * g]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        pts += [[s1 * (2 + g), 0, s2 * (g**2)]
                for s1, s2 in itertools.product([1, -1], repeat=2)]
        # the vertices are all even permutations of the lists in pts
        verts = pts
        verts += [[v[1], v[2], v[0]] for v in pts]
        verts += [[v[2], v[0], v[1]] for v in pts]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def truncated_icosidodecahedron(self, exact=True, base_ring=None, backend=None):
        r"""
        Return the truncated icosidodecahedron.

        The truncated icosidodecahedron is an Archimedean solid.
        It has 62 faces and 120 vertices. For more information, see the
        :wikipedia:`Truncated_icosidodecahedron`.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) If ``False`` use an
          approximate ring for the coordinates.

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided and ``exact=True`` it will be a the number field
          `\QQ[\phi]` where `\phi` is the golden ratio and if ``exact=False``
          it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: ti = polytopes.truncated_icosidodecahedron()  # long time
            sage: ti.f_vector()                                 # long time
            (1, 120, 180, 62, 1)
            sage: ti.base_ring()                                # long time
            Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?

        The implementation using floating point approximations is much faster::

            sage: ti = polytopes.truncated_icosidodecahedron(exact=False) # random
            sage: ti.f_vector()
            (1, 120, 180, 62, 1)
            sage: ti.base_ring()
            Real Double Field

        Its facets are 30 squares, 20 hexagons and 12 decagons::

            sage: sum(1 for f in ti.facets() if len(f.vertices()) == 4)
            30
            sage: sum(1 for f in ti.facets() if len(f.vertices()) == 6)
            20
            sage: sum(1 for f in ti.facets() if len(f.vertices()) == 10)
            12

        TESTS::

            sage: ti = polytopes.truncated_icosidodecahedron(backend='normaliz')  # optional - pynormaliz
            sage: ti.f_vector()
            (1, 120, 180, 62, 1)
            sage: ti.base_ring()                                                  # optional - pynormaliz
            Number Field in sqrt5 with defining polynomial x^2 - 5 with sqrt5 = 2.236067977499790?

        """
        if base_ring is None and exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            if base_ring is None:
                from sage.rings.real_double import RDF as base_ring
            g = (1 + base_ring(5).sqrt()) / 2

        pts = [[s1 * 1 / g, s2 * 1 / g, s3 * (3 + g)]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        pts += [[s1 * 2 / g, s2 * g, s3 * (1 + 2 * g)]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        pts += [[s1 * 1 / g, s2 * (g**2), s3 * (-1 + 3 * g)]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        pts += [[s1 * (-1 + 2 * g), s2 * 2 * base_ring.one(), s3 * (2 + g)]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        pts += [[s1 * g, s2 * 3 * base_ring.one(), s3 * 2 * g]
                for s1, s2, s3 in itertools.product([1, -1], repeat=3)]
        # the vertices are all ever permutations of the lists in pts
        verts = pts
        verts += [[v[1], v[2], v[0]] for v in pts]
        verts += [[v[2], v[0], v[1]] for v in pts]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def snub_dodecahedron(self, base_ring=None, backend=None, verbose=False):
        """
        Return the snub dodecahedron.

        The snub dodecahedron is an Archimedean solid.
        It has 92 faces and 60 vertices. For more information, see the
        :wikipedia:`Snub_dodecahedron`.

        INPUT:

        - ``base_ring`` -- the ring in which the coordinates will belong to. If
          it is not provided it will be the real double field.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES:

        Only the backend using the optional normaliz package can construct
        the snub dodecahedron in reasonable time::

            sage: sd = polytopes.snub_dodecahedron(base_ring=AA,        # optional - pynormaliz, long time
            ....:                                  backend='normaliz')
            sage: sd.f_vector()                                         # optional - pynormaliz, long time
            (1, 60, 150, 92, 1)
            sage: sd.base_ring()                                        # optional - pynormaliz, long time
            Algebraic Real Field

        Its facets are 80 triangles and 12 pentagons::

            sage: sum(1 for f in sd.facets()                            # optional - pynormaliz, long time
            ....:       if len(f.vertices()) == 3)
            80
            sage: sum(1 for f in sd.facets()                            # optional - pynormaliz, long time
            ....:       if len(f.vertices()) == 5)
            12

        TESTS:

        The cdd backend with floating point arithmetic fails for this
        polytope::

            sage: sd = polytopes.snub_dodecahedron()        # not tested
            sage: sd.f_vector()                             # not tested
            (1, 60, 150, 92, 1)
            sage: sd.base_ring()                            # not tested
            Real Double Field

        """
        if base_ring is None:
            from sage.rings.real_double import RDF as base_ring
        phi = (1 + base_ring(5).sqrt()) / 2
        xi = ((phi/2 + (phi - ZZ(5)/27).sqrt()/2)**(~ZZ(3)) +
              (phi/2 - (phi - ZZ(5)/27).sqrt()/2)**(~ZZ(3)))

        alpha = xi - 1 / xi
        beta = xi * phi + phi**2 + phi / xi
        signs = [[-1,-1,-1], [-1,1,1], [1,-1,1], [1,1,-1]]

        pts = [[s1 * 2 * alpha, s2 * 2 * base_ring.one(), s3 * 2 * beta]
                for s1, s2, s3 in signs]
        pts += [[s1 * (alpha + beta/phi + phi), s2 * (-alpha * phi + beta + 1/phi), s3 * (alpha/phi + beta * phi - 1)]
                for s1, s2, s3 in signs]
        pts += [[s1 * (alpha + beta/phi - phi), s2 * (alpha * phi - beta + 1/phi), s3 * (alpha/phi + beta * phi + 1)]
                for s1, s2, s3 in signs]
        pts += [[s1 * (-alpha/phi + beta * phi + 1), s2 * (-alpha + beta/phi - phi), s3 * (alpha * phi + beta - 1/phi)]
                for s1, s2, s3 in signs]
        pts += [[s1 * (-alpha/phi + beta * phi - 1), s2 * (alpha - beta/phi - phi), s3 * (alpha * phi + beta + 1/phi)]
                for s1, s2, s3 in signs]

        # the vertices are all even permutations of the lists in pts
        verts = pts
        verts += [[v[1], v[2], v[0]] for v in pts]
        verts += [[v[2], v[0], v[1]] for v in pts]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend, verbose=verbose)

    def twenty_four_cell(self, backend=None):
        """
        Return the standard 24-cell polytope.

        The 24-cell polyhedron (also called icositetrachoron or octaplex) is a
        regular polyhedron in 4-dimension. For more information see
        the :wikipedia:`24-cell`.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: p24 = polytopes.twenty_four_cell()
            sage: p24.f_vector()
            (1, 24, 96, 96, 24, 1)
            sage: v = next(p24.vertex_generator())
            sage: for adj in v.neighbors(): print(adj)
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

        TESTS::

            sage: tfcell = polytopes.twenty_four_cell(backend='normaliz')    # optional - pynormaliz
            sage: TestSuite(tfcell).run()                                    # optional - pynormaliz
        """
        q12 = QQ((1, 2))
        verts = list(itertools.product([q12, -q12], repeat=4))
        B4 = (ZZ**4).basis()
        verts.extend(B4)
        verts.extend(-v for v in B4)
        return Polyhedron(vertices=verts, backend=backend)

    def runcitruncated_six_hundred_cell(self, exact=True, backend=None):
        """
        Return the runcitruncated 600-cell.

        The runcitruncated 600-cell is a 4-dimensional 4-uniform polytope in
        the `H_4` family. It has 7200 vertices. For more information see
        :wikipedia:`Runcitruncated 600-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.runcitruncated_six_hundred_cell(backend='normaliz')  # not tested - very long time
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of
            7200 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[1, 1, 0, 1], exact=exact, backend=backend, regular=True)

    def cantitruncated_six_hundred_cell(self, exact=True, backend=None):
        """
        Return the cantitruncated 600-cell.

        The cantitruncated 600-cell is a 4-dimensional 4-uniform polytope in
        the `H_4` family. It has 7200 vertices. For more information see
        :wikipedia:`Cantitruncated 600-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.cantitruncated_six_hundred_cell(exact=True,          # not tested - very long time
            ....:                                           backend='normaliz')
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 7200 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[1, 1, 1, 0], exact=exact, backend=backend, regular=True)

    def bitruncated_six_hundred_cell(self, exact=True, backend=None):
        """
        Return the bitruncated 600-cell.

        The bitruncated 600-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 3600 vertices. For more information see
        :wikipedia:`Bitruncated 600-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.runcinated_six_hundred_cell(exact=True,              # not tested - very long time
            ....:                                       backend='normaliz')
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 3600 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[0, 1, 1, 0], exact=exact, backend=backend, regular=True)

    def cantellated_six_hundred_cell(self, exact=False, backend=None):
        """
        Return the cantellated 600-cell.

        The cantellated 600-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 3600 vertices. For more information see
        :wikipedia:`Cantellated 600-cell`.

        .. WARNING::

            The coordinates are inexact by default. The computation with
            inexact coordinates (using the backend ``'cdd'``) issues a
            UserWarning on inconsistencies.

        INPUT:

        - ``exact`` - (boolean, default ``False``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.cantellated_six_hundred_cell()       # not tested - very long time
            doctest:warning
            ...
            UserWarning: This polyhedron data is numerically complicated; cdd
            could not convert between the inexact V and H representation
            without loss of data. The resulting object might show
            inconsistencies.
            A 4-dimensional polyhedron in RDF^4 defined as the convex hull of 3600 vertices

        It is possible to use the backend ``'normaliz'`` to get an exact
        representation::

            sage: polytopes.cantellated_six_hundred_cell(exact=True,  # not tested - long time
            ....:                                        backend='normaliz')
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 3600 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[1, 0, 1, 0], exact=exact, backend=backend, regular=True)

    def truncated_six_hundred_cell(self, exact=False, backend=None):
        """
        Return the truncated 600-cell.

        The truncated 600-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 1440 vertices. For more information see
        :wikipedia:`Truncated 600-cell`.

        .. WARNING::

            The coordinates are not exact by default. The computation with
            exact coordinates takes a huge amount of time.

        INPUT:

        - ``exact`` - (boolean, default ``False``) if ``True`` use exact
          coordinates instead of floating point approximations

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.truncated_six_hundred_cell()  # not tested - long time
            A 4-dimensional polyhedron in RDF^4 defined as the convex hull of 1440 vertices

        It is possible to use the backend ``'normaliz'`` to get an exact
        representation::

            sage: polytopes.truncated_six_hundred_cell(exact=True,backend='normaliz') # not tested - long time ~16sec
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 1440 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[1, 1, 0, 0], exact=exact, backend=backend, regular=True)

    def rectified_six_hundred_cell(self, exact=True, backend=None):
        """
        Return the rectified 600-cell.

        The rectified 600-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 720 vertices. For more information see
        :wikipedia:`Rectified 600-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.rectified_six_hundred_cell(backend='normaliz')  # not tested - long time ~14sec
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 720 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[0, 1, 0, 0], exact=exact, backend=backend, regular=True)

    def six_hundred_cell(self, exact=False, backend=None):
        """
        Return the standard 600-cell polytope.

        The 600-cell is a 4-dimensional regular polytope. In many ways this is
        an analogue of the icosahedron.

        .. WARNING::

            The coordinates are not exact by default. The computation with
            exact coordinates takes a huge amount of time.

        INPUT:

        - ``exact`` - (boolean, default ``False``) if ``True`` use exact
          coordinates instead of floating point approximations

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: p600 = polytopes.six_hundred_cell(); p600                             # needs sage.groups
            A 4-dimensional polyhedron in RDF^4 defined as the convex hull of 120 vertices
            sage: p600.f_vector()               # long time (~2sec)                     # needs sage.groups
            (1, 120, 720, 1200, 600, 1)

        Computation with exact coordinates is currently too long to be useful::

            sage: p600 = polytopes.six_hundred_cell(exact=True)         # long time, not tested, needs sage.groups
            sage: len(list(p600.bounded_edges()))                       # long time, not tested, needs sage.groups
            720

        TESTS::

            sage: # optional - pynormaliz, needs sage.groups sage.rings.number_field
            sage: p600 = polytopes.six_hundred_cell(exact=True,
            ....:                                   backend='normaliz')
            sage: len(list(p600.bounded_edges()))                       # long time
            720
        """
        if exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            from sage.rings.real_double import RDF as base_ring
            g = (1 + base_ring(5).sqrt()) / 2

        q12 = base_ring(1) / base_ring(2)
        z = base_ring.zero()
        verts = [[s1*q12, s2*q12, s3*q12, s4*q12] for s1,s2,s3,s4 in itertools.product([1,-1], repeat=4)]
        V = (base_ring)**4
        verts.extend(V.basis())
        verts.extend(-v for v in V.basis())
        pts = [[s1 * q12, s2*g/2, s3/(2*g), z] for (s1,s2,s3) in itertools.product([1,-1], repeat=3)]
        for p in AlternatingGroup(4):
            verts.extend(p(x) for x in pts)
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def grand_antiprism(self, exact=True, backend=None, verbose=False):
        """
        Return the grand antiprism.

        The grand antiprism is a 4-dimensional non-Wythoffian uniform polytope.
        The coordinates were taken from http://eusebeia.dyndns.org/4d/gap. For
        more information, see the :wikipedia:`Grand_antiprism`.

        .. WARNING::

            The coordinates are exact by default. The computation with exact
            coordinates is not as fast as with floating point approximations.
            If you find this method to be too slow, consider using floating
            point approximations

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``False`` use floating
          point approximations instead of exact coordinates

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: gap = polytopes.grand_antiprism()  # not tested - very long time
            sage: gap                                # not tested - very long time
            A 4-dimensional polyhedron in (Number Field in sqrt5 with defining
            polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^4 defined as
            the convex hull of 100 vertices

        Computation with the backend ``'normaliz'`` is instantaneous::

            sage: gap_norm = polytopes.grand_antiprism(backend='normaliz')      # optional - pynormaliz, needs sage.rings.number_field
            sage: gap_norm                                                      # optional - pynormaliz, needs sage.rings.number_field
            A 4-dimensional polyhedron in (Number Field in sqrt5 with defining
            polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^4 defined as
            the convex hull of 100 vertices

        Computation with approximated coordinates is also faster, but inexact::

            sage: gap = polytopes.grand_antiprism(exact=False) # random
            sage: gap
            A 4-dimensional polyhedron in RDF^4 defined as the convex hull of 100 vertices
            sage: gap.f_vector()
            (1, 100, 500, 720, 320, 1)
            sage: len(list(gap.bounded_edges()))
            500
        """
        from itertools import product

        if exact:
            from sage.rings.number_field.number_field import QuadraticField
            K = QuadraticField(5, 'sqrt5')
            sqrt5 = K.gen()
            g = (1 + sqrt5) / 2
            base_ring = K
        else:
            from sage.rings.real_double import RDF as base_ring
            g = (1 + base_ring(5).sqrt()) / 2

        q12 = base_ring(1) / base_ring(2)
        z = base_ring.zero()
        verts = [[s1*q12, s2*q12, s3*q12, s4*q12] for s1,s2,s3,s4 in product([1,-1], repeat=4)]
        V = (base_ring)**4
        verts.extend(V.basis()[2:])
        verts.extend(-v for v in V.basis()[2:])

        verts.extend([s1 * q12, s2/(2*g), s3*g/2, z] for (s1,s2,s3) in product([1,-1], repeat=3))
        verts.extend([s3*g/2, s1 * q12, s2/(2*g), z] for (s1,s2,s3) in product([1,-1], repeat=3))
        verts.extend([s2/(2*g), s3*g/2, s1 * q12, z] for (s1,s2,s3) in product([1,-1], repeat=3))

        verts.extend([s1 * q12, s2*g/2, z, s3/(2*g)] for (s1,s2,s3) in product([1,-1], repeat=3))
        verts.extend([s3/(2*g), s1 * q12, z, s2*g/2] for (s1,s2,s3) in product([1,-1], repeat=3))
        verts.extend([s2*g/2, s3/(2*g), z, s1 * q12] for (s1,s2,s3) in product([1,-1], repeat=3))

        verts.extend([s1 * q12, z, s2/(2*g), s3*g/2] for (s1,s2,s3) in product([1,-1], repeat=3))

        verts.extend([z, s1 * q12, s2*g/2, s3/(2*g)] for (s1,s2,s3) in product([1,-1], repeat=3))

        verts.extend([z, s1/(2*g), q12, g/2] for s1 in [1, -1])
        verts.extend([z, s1/(2*g), -q12, -g/2] for s1 in [1, -1])

        verts.extend([z, s1*g/2, 1/(2*g), q12] for s1 in [1, -1])
        verts.extend([z, s1*g/2, -1/(2*g), -q12] for s1 in [1, -1])

        verts.extend([s1*g/2, z, q12, -1/(2*g)] for s1 in [1, -1])
        verts.extend([s1*g/2, z, -q12, 1/(2*g)] for s1 in [1, -1])

        verts.extend([s1/(2*g), z, g/2, -q12] for s1 in [1, -1])
        verts.extend([s1/(2*g), z, -g/2, q12] for s1 in [1, -1])

        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend, verbose=verbose)

    def Gosset_3_21(self, backend=None):
        r"""
        Return the Gosset `3_{21}` polytope.

        The Gosset `3_{21}` polytope is a uniform 7-polytope. It has 56
        vertices, and 702 facets: `126` `3_{11}` and `576` `6`-simplex. For
        more information, see the :wikipedia:`3_21_polytope`.

        INPUT:

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: g = polytopes.Gosset_3_21(); g
            A 7-dimensional polyhedron in ZZ^8 defined as the convex hull of 56 vertices
            sage: g.f_vector() # not tested (~16s)
            (1, 56, 756, 4032, 10080, 12096, 6048, 702, 1)

        TESTS::

            sage: G321 = polytopes.Gosset_3_21(backend='normaliz')   # optional - pynormaliz
            sage: TestSuite(G321).run()                              # optional - pynormaliz, long time
        """
        from itertools import combinations
        verts = []
        for i, j in combinations(range(8), 2):
            x = [1]*8
            x[i] = x[j] = -3
            verts.append(x)
            verts.append([-xx for xx in x])

        return Polyhedron(vertices=verts, base_ring=ZZ, backend=backend)

    def cyclic_polytope(self, dim, n, base_ring=QQ, backend=None):
        r"""
        Return a cyclic polytope.

        A cyclic polytope of dimension ``dim`` with ``n`` vertices is the
        convex hull of the points  ``(t,t^2,...,t^dim)`` with `t \in
        \{0,1,...,n-1\}` .  For more information, see the
        :wikipedia:`Cyclic_polytope`.

        INPUT:

        - ``dim`` -- positive integer. the dimension of the polytope.

        - ``n`` -- positive integer. the number of vertices.

        - ``base_ring`` -- either ``QQ`` (default) or ``RDF``.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: c = polytopes.cyclic_polytope(4,10)
            sage: c.f_vector()
            (1, 10, 45, 70, 35, 1)

        TESTS::

            sage: cp = polytopes.cyclic_polytope(4,10,backend='normaliz')  # optional - pynormaliz
            sage: TestSuite(cp).run()                                      # optional - pynormaliz
        """
        verts = [[t**i for i in range(1, dim+1)] for t in range(n)]
        return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

    def hypersimplex(self, dim, k, project=False, backend=None):
        r"""
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
          is (isometrically) projected to a vector space of dimension
          ``dim-1``.  This operation turns the coordinates into floating point
          approximations and corresponds to the projection given by the matrix
          from :func:`zero_sum_projection`.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: h_4_2 = polytopes.hypersimplex(4, 2)
            sage: h_4_2
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 6 vertices
            sage: h_4_2.f_vector()
            (1, 6, 12, 8, 1)
            sage: h_4_2.ehrhart_polynomial()                            # optional - latte_int
            2/3*t^3 + 2*t^2 + 7/3*t + 1
            sage: TestSuite(h_4_2).run()

            sage: # needs sage.combinat
            sage: h_7_3 = polytopes.hypersimplex(7, 3, project=True)
            sage: h_7_3
            A 6-dimensional polyhedron in RDF^6 defined as the convex hull of 35 vertices
            sage: h_7_3.f_vector()
            (1, 35, 210, 350, 245, 84, 14, 1)
            sage: TestSuite(h_7_3).run(skip=["_test_pyramid", "_test_lawrence"])
        """
        verts = Permutations([0] * (dim - k) + [1] * k).list()
        if project:
            verts = project_points(*verts)
        return Polyhedron(vertices=verts, backend=backend)

    def permutahedron(self, n, project=False, backend=None):
        r"""
        Return the standard permutahedron of (1,...,n).

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
          is (isometrically) projected to a vector space of dimension
          ``dim-1``.  This operation turns the coordinates into floating point
          approximations and corresponds to the projection given by the matrix
          from :func:`zero_sum_projection`.

        - ``backend`` -- the backend to use to create the polytope.

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
            sage: perm4.plot()                                                          # needs sage.plot
            Graphics3d Object
            sage: perm4.graph().is_isomorphic(graphs.BubbleSortGraph(4))                # needs sage.graphs
            True

        As both Hrepresentation an Vrepresentation are known, the permutahedron can be set
        up with both using the backend ``field``. The following takes very very long time
        to recompute, e.g. with backend ``ppl``::

            sage: polytopes.permutahedron(8, backend='field')  # (~1s)
            A 7-dimensional polyhedron in QQ^8 defined as the convex hull of 40320 vertices
            sage: polytopes.permutahedron(9, backend='field')  # not tested (memory consumption)  # (~5s)
            A 8-dimensional polyhedron in QQ^9 defined as the convex hull of 362880 vertices

        .. SEEALSO::

            * :meth:`~sage.graphs.graph_generators.GraphGenerators.BubbleSortGraph`

        TESTS::

            sage: p4 = polytopes.permutahedron(4,backend='normaliz')   # optional - pynormaliz
            sage: TestSuite(p4).run()                                  # optional - pynormaliz

        Check that precomputed data is correct::

            sage: P = polytopes.permutahedron(5, backend='field')
            sage: TestSuite(P).run()  # long time
        """
        verts = itertools.permutations(range(1, n + 1))
        if project:
            verts = project_points(*verts)
            return Polyhedron(vertices=verts, backend=backend)
        else:
            parent = Polyhedra(ZZ, n, backend=backend)

            def tri(m):
                return (m * (m + 1)) // 2

            # Each proper `S \subset [n]` corresponds exactly to
            # a facet that minimizes the coordinates in `S`.
            # The minimal sum for `m` coordinates is `(m*(m+1))/2`.
            ieqs = ((-tri(sum(x)),) + x
                    for x in itertools.product([0,1], repeat=n)
                    if 0 < sum(x) < n)

            # Adding the defining equality.
            eqns = ((-tri(n),) + tuple(1 for _ in range(n)),)

            return parent([verts, [], []], [ieqs, eqns],
                          Vrep_minimal=True, Hrep_minimal=True, pref_rep="Hrep")

    def generalized_permutahedron(self, coxeter_type, point=None, exact=True, regular=False, backend=None):
        r"""
        Return the generalized permutahedron of type ``coxeter_type`` as the
        convex hull of the orbit of ``point`` in the fundamental cone.

        This generalized permutahedron lies in the vector space used in the
        geometric representation, that is, in the default case, the dimension
        of generalized permutahedron equals the dimension of the space.

        INPUT:

        - ``coxeter_type`` -- a Coxeter type; given as a pair [type,rank],
          where type is a letter and rank is the number of generators.

        - ``point`` -- a list (default: ``None``); a point given by its
          coordinates in the weight basis. If ``None`` is given, the point
          `(1, 1, 1, \ldots)` is used.

        - ``exact`` - (boolean, default ``True``) if ``False`` use floating
          point approximations instead of exact coordinates

        - ``regular`` -- boolean (default: ``False``); whether to apply a
          linear transformation making the vertex figures isometric.

        - ``backend`` -- backend to use to create the polytope; (default:
          ``None``)

        EXAMPLES::

            sage: perm_a3 = polytopes.generalized_permutahedron(['A',3]); perm_a3       # needs sage.combinat
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices

        You can put the starting point along the hyperplane of the first
        generator::

            sage: # needs sage.combinat
            sage: perm_a3_011 = polytopes.generalized_permutahedron(['A',3], [0,1,1])
            sage: perm_a3_011
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 12 vertices
            sage: perm_a3_110 = polytopes.generalized_permutahedron(['A',3], [1,1,0])
            sage: perm_a3_110
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 12 vertices
            sage: perm_a3_110.is_combinatorially_isomorphic(perm_a3_011)
            True
            sage: perm_a3_101 = polytopes.generalized_permutahedron(['A',3], [1,0,1])
            sage: perm_a3_101
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 12 vertices
            sage: perm_a3_110.is_combinatorially_isomorphic(perm_a3_101)
            False
            sage: perm_a3_011.f_vector()
            (1, 12, 18, 8, 1)
            sage: perm_a3_101.f_vector()
            (1, 12, 24, 14, 1)

        The usual output does not necessarily give a polyhedron with isometric
        vertex figures::

            sage: perm_a2 = polytopes.generalized_permutahedron(['A',2])                # needs sage.combinat
            sage: perm_a2.vertices()                                                    # needs sage.combinat
            (A vertex at (-1, -1),
             A vertex at (-1, 0),
             A vertex at (0, -1),
             A vertex at (0, 1),
             A vertex at (1, 0),
             A vertex at (1, 1))

        It works also with Coxeter types that lead to non-rational coordinates::

            sage: perm_b3 = polytopes.generalized_permutahedron(['B',3])        # long time, needs sage.combinat sage.rings.number_field
            sage: perm_b3                                                       # long time, needs sage.combinat sage.rings.number_field
            A 3-dimensional polyhedron in
             (Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?)^3
             defined as the convex hull of 48 vertices

        Setting ``regular=True`` applies a linear transformation to get
        isometric vertex figures and the result is inscribed. This cannot be done using
        rational coordinates. We first do the computations using floating point
        approximations (``RDF``)::

            sage: perm_a2_inexact = polytopes.generalized_permutahedron(                # needs sage.combinat
            ....:     ['A',2], exact=False)
            sage: sorted(perm_a2_inexact.vertices())                                    # needs sage.combinat
            [A vertex at (-1.0, -1.0),
             A vertex at (-1.0, 0.0),
             A vertex at (0.0, -1.0),
             A vertex at (0.0, 1.0),
             A vertex at (1.0, 0.0),
             A vertex at (1.0, 1.0)]

            sage: perm_a2_inexact_reg = polytopes.generalized_permutahedron(            # needs sage.combinat
            ....:     ['A',2], exact=False, regular=True)
            sage: sorted(perm_a2_inexact_reg.vertices())                                # needs sage.combinat
            [A vertex at (-1.0, 0.0),
             A vertex at (-0.5, -0.8660254038),
             A vertex at (-0.5, 0.8660254038),
             A vertex at (0.5, -0.8660254038),
             A vertex at (0.5, 0.8660254038),
             A vertex at (1.0, 0.0)]

        We can do the same computation using exact arithmetic with the field ``AA``::

            sage: perm_a2_reg = polytopes.generalized_permutahedron(                    # needs sage.combinat sage.rings.number_field
            ....:     ['A',2], regular=True)
            sage: V = sorted(perm_a2_reg.vertices()); V  # random                       # needs sage.combinat sage.rings.number_field
            [A vertex at (-1, 0),
             A vertex at (-1/2, -0.866025403784439?),
             A vertex at (-1/2, 0.866025403784439?),
             A vertex at (1/2, -0.866025403784439?),
             A vertex at (1/2, 0.866025403784439?),
             A vertex at (1.000000000000000?, 0.?e-18)]

        Even though the numbers look like floating point approximations, the computation is
        actually exact. We can clean up the display a bit using ``exactify``::

            sage: for v in V:                                                           # needs sage.combinat sage.rings.number_field
            ....:     for x in v:
            ....:         x.exactify()
            sage: V                                                                     # needs sage.combinat sage.rings.number_field
            [A vertex at (-1, 0),
             A vertex at (-1/2, -0.866025403784439?),
             A vertex at (-1/2, 0.866025403784439?),
             A vertex at (1/2, -0.866025403784439?),
             A vertex at (1/2, 0.866025403784439?),
             A vertex at (1, 0)]
            sage: perm_a2_reg.is_inscribed()                                            # needs sage.combinat sage.rings.number_field
            True

        Larger examples take longer::

            sage: # needs sage.combinat sage.rings.number_field
            sage: perm_a3_reg = polytopes.generalized_permutahedron(    # long time
            ....:     ['A',3], regular=True); perm_a3_reg
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 24 vertices
            sage: perm_a3_reg.is_inscribed()    # long time
            True
            sage: perm_b3_reg = polytopes.generalized_permutahedron(    # long time (12sec on 64 bits), not tested
            ....:     ['B',3], regular=True); perm_b3_reg
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 48 vertices

        It is faster with the backend ``'number_field'``, which internally uses an embedded
        number field instead of doing the computations directly with the base ring (``AA``)::

            sage: # needs sage.combinat sage.rings.number_field
            sage: perm_a3_reg_nf = polytopes.generalized_permutahedron(
            ....:    ['A',3], regular=True, backend='number_field'); perm_a3_reg_nf
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 24 vertices
            sage: perm_a3_reg_nf.is_inscribed()
            True
            sage: perm_b3_reg_nf = polytopes.generalized_permutahedron(         # long time
            ....:    ['B',3], regular=True, backend='number_field'); perm_b3_reg_nf
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 48 vertices

        It is even faster with the backend ``'normaliz'``::

            sage: # optional - pynormaliz, needs sage.combinat sage.rings.number_field
            sage: perm_a3_reg_norm = polytopes.generalized_permutahedron(
            ....:    ['A',3], regular=True, backend='normaliz'); perm_a3_reg_norm
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 24 vertices
            sage: perm_a3_reg_norm.is_inscribed()
            True
            sage: perm_b3_reg_norm = polytopes.generalized_permutahedron(
            ....:    ['B',3], regular=True, backend='normaliz'); perm_b3_reg_norm
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 48 vertices

        The speedups from using backend ``'normaliz'`` allow us to go even further::

            sage: # optional - pynormaliz, needs sage.combinat sage.rings.number_field
            sage: perm_h3 = polytopes.generalized_permutahedron(
            ....:     ['H',3], backend='normaliz'); perm_h3
            A 3-dimensional polyhedron in
             (Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?)^3
             defined as the convex hull of 120 vertices
            sage: perm_f4 = polytopes.generalized_permutahedron(        # long time
            ....:     ['F',4], backend='normaliz'); perm_f4
            A 4-dimensional polyhedron
             in (Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?)^4
             defined as the convex hull of 1152 vertices

        .. SEEALSO::

            * :meth:`~sage.combinat.root_system.reflection_group_real.permutahedron`
            * :meth:`~sage.categories.finite_coxeter_groups.permutahedron`

        TESTS::

            sage: TestSuite(perm_h3).run()      # optional - pynormaliz                 # needs sage.combinat sage.rings.number_field
        """
        from sage.combinat.root_system.coxeter_group import CoxeterGroup
        try:
            W = CoxeterGroup(coxeter_type)
        except (TypeError, ValueError):
            raise ValueError("cannot build a Coxeter group from {}".format(coxeter_type))
        n = W.one().canonical_matrix().rank()
        weights = W.fundamental_weights()
        if point is None:
            point = [ZZ.one()] * n
        apex = sum(point[i-1] * weights[i] for i in weights.keys())
        # Try to rationalize the starting point
        non_zero_index = list(apex).index([x for x in apex if x != 0][0])
        apex = (QQ(1)/apex[non_zero_index]) * apex
        apex.set_immutable()
        vertices = set()
        # This does not work well with UCF, so we set it to None:
        # br = apex.base_ring()
        br = None
        for w in W:
            # The apex is considered in the space on which it acts and not in
            # the weight space.
            new_point = w * apex
            new_point.set_immutable()
            vertices.add(new_point)
        if regular:
            from sage.rings.qqbar import AA
            from sage.matrix.constructor import matrix
            from sage.modules.free_module_element import vector
            # This transformation fixes the first root and adjust the other
            # roots to have the correct angles
            bf = W.bilinear_form()
            transf_col = [[1] + [0]*(n-1)]
            for i in range(1, n):
                new_col = [0]*i + [1] + [0]*(n-i-1)
                transf_col += [new_col]
                m = matrix(AA, transf_col)
                col = bf.column(i)
                rhs = vector(AA, list(col[:i+1]))
                adjusted_col = m.solve_right(rhs)
                # Then scales the images so that the polytope is inscribed
                c = 1 - sum(adjusted_col[j]**2 for j in range(n) if j != i)
                c = c.sqrt()
                adjusted_col[i] = c
                transf_col[-1] = adjusted_col
            # TODO: Make this matrix into the cyclotomics, the value of c is an
            # algebraic number not anymore in the cyclotomic field.
            transf = matrix(transf_col).transpose()
            vertices = [transf * v.change_ring(AA) for v in vertices]
            br = AA
        if not exact:
            from sage.rings.real_double import RDF
            vertices = [v.change_ring(RDF) for v in vertices]
            br = RDF
        return Polyhedron(vertices=vertices, backend=backend, base_ring=br)

    def omnitruncated_one_hundred_twenty_cell(self, exact=True, backend=None):
        """
        Return the omnitruncated 120-cell.

        The omnitruncated 120-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 14400 vertices. For more information see
        :wikipedia:`Omnitruncated 120-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.omnitruncated_one_hundred_twenty_cell(backend='normaliz')  # not tested - very long time ~10min
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 14400 vertices
        """
        if not exact:
            # cdd finds a numerical inconsistency.
            raise NotImplementedError("cannot compute the convex hull using floating points")
        return self.generalized_permutahedron(['H', 4], exact=exact, backend=backend, regular=True)

    omnitruncated_six_hundred_cell = omnitruncated_one_hundred_twenty_cell

    def runcitruncated_one_hundred_twenty_cell(self, exact=False, backend=None):
        """
        Return the runcitruncated 120-cell.

        The runcitruncated 120-cell is a 4-dimensional 4-uniform polytope in
        the `H_4` family. It has 7200 vertices. For more information see
        :wikipedia:`Runcitruncated 120-cell`.

        .. WARNING::

            The coordinates are inexact by default. The computation with
            inexact coordinates (using the backend ``'cdd'``) issues a
            UserWarning on inconsistencies.

        INPUT:

        - ``exact`` - (boolean, default ``False``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.runcitruncated_one_hundred_twenty_cell(exact=False)  # not tested - very long time
            doctest:warning
            ...
            UserWarning: This polyhedron data is numerically complicated; cdd
            could not convert between the inexact V and H representation
            without loss of data. The resulting object might show
            inconsistencies.

        It is possible to use the backend ``'normaliz'`` to get an exact
        representation::

            sage: polytopes.runcitruncated_one_hundred_twenty_cell(exact=True,   # not tested - very long time
            ....:                                                  backend='normaliz')
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 7200 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[1, 0, 1, 1], exact=exact, backend=backend, regular=True)

    def cantitruncated_one_hundred_twenty_cell(self, exact=True, backend=None):
        """
        Return the cantitruncated 120-cell.

        The cantitruncated 120-cell is a 4-dimensional 4-uniform polytope in
        the `H_4` family. It has 7200 vertices. For more information see
        :wikipedia:`Cantitruncated 120-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.cantitruncated_one_hundred_twenty_cell(exact=True, backend='normaliz')  # not tested - very long time
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 7200 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[0, 1, 1, 1], exact=exact, backend=backend, regular=True)

    def runcinated_one_hundred_twenty_cell(self, exact=False, backend=None):
        """
        Return the runcinated 120-cell.

        The runcinated 120-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 2400 vertices. For more information see
        :wikipedia:`Runcinated 120-cell`.

        .. WARNING::

            The coordinates are inexact by default. The computation with
            inexact coordinates (using the backend ``'cdd'``) issues a
            UserWarning on inconsistencies.

        INPUT:

        - ``exact`` - (boolean, default ``False``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.runcinated_one_hundred_twenty_cell(exact=False)  # not tested - very long time
            doctest:warning ...  UserWarning: This polyhedron data is
            numerically complicated; cdd could not convert between the inexact
            V and H representation without loss of data. The resulting object
            might show inconsistencies.
            A 4-dimensional polyhedron in RDF^4 defined as the convex hull of 2400 vertices

        It is possible to use the backend ``'normaliz'`` to get an exact
        representation::

            sage: polytopes.runcinated_one_hundred_twenty_cell(exact=True,   # not tested - very long time
            ....:                                              backend='normaliz')
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 2400 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[1, 0, 0, 1], exact=exact, backend=backend, regular=True)

    def cantellated_one_hundred_twenty_cell(self, exact=True, backend=None):
        """
        Return the cantellated 120-cell.

        The cantellated 120-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 3600 vertices. For more information see
        :wikipedia:`Cantellated 120-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.cantellated_one_hundred_twenty_cell(backend='normaliz')  # not tested - long time
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 3600 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[0, 1, 0, 1], exact=exact, backend=backend, regular=True)

    def truncated_one_hundred_twenty_cell(self, exact=True, backend=None):
        """
        Return the truncated 120-cell.

        The truncated 120-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 2400 vertices. For more information see
        :wikipedia:`Truncated 120-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.truncated_one_hundred_twenty_cell(backend='normaliz')  # not tested - long time
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 2400 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[0, 0, 1, 1], exact=exact, backend=backend, regular=True)

    def rectified_one_hundred_twenty_cell(self, exact=True, backend=None):
        """
        Return the rectified 120-cell.

        The rectified 120-cell is a 4-dimensional 4-uniform polytope in the
        `H_4` family. It has 1200 vertices. For more information see
        :wikipedia:`Rectified 120-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` - (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.rectified_one_hundred_twenty_cell(backend='normaliz')  # not tested - long time
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 1200 vertices
        """
        return self.generalized_permutahedron(['H', 4], point=[0, 0, 1, 0], exact=exact, backend=backend, regular=True)

    def one_hundred_twenty_cell(self, exact=True, backend=None, construction='coxeter'):
        """
        Return the 120-cell.

        The 120-cell is a 4-dimensional 4-uniform polytope in the `H_4` family.
        It has 600 vertices and 120 facets. For more information see
        :wikipedia:`120-cell`.

        .. WARNING::

            The coordinates are exact by default. The computation with inexact
            coordinates (using the backend ``'cdd'``) returns a numerical
            inconsistency error, and thus cannot be computed.

        INPUT:

        - ``exact`` -- (boolean, default ``True``) if ``True`` use exact
          coordinates instead of floating point approximations.

        - ``backend`` -- the backend to use to create the polytope.

        - ``construction`` -- the construction to use (string, default 'coxeter');
          the other possibility is 'as_permutahedron'.

`       EXAMPLES:

        The classical construction given by Coxeter in [Cox1969]_ is given by::

            sage: polytopes.one_hundred_twenty_cell()                    # not tested - long time ~15 sec.
            A 4-dimensional polyhedron in (Number Field in sqrt5 with defining
             polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^4 defined as
             the convex hull of 600 vertices

        The ``'normaliz'`` is faster::

            sage: P = polytopes.one_hundred_twenty_cell(backend='normaliz'); P  # optional - pynormaliz
            A 4-dimensional polyhedron in (Number Field in sqrt5 with defining
             polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^4 defined as
             the convex hull of 600 vertices

        It is also possible to realize it using the generalized permutahedron
        of type `H_4`::

            sage: polytopes.one_hundred_twenty_cell(backend='normaliz',  # not tested - long time
            ....:                                   construction='as_permutahedron')
            A 4-dimensional polyhedron in AA^4 defined as the convex hull of 600 vertices

        TESTS::

            sage: TestSuite(P).run()                                 # long time, optional - pynormaliz
        """
        if construction == 'coxeter':
            if not exact:
                raise ValueError("The 'cdd' backend produces numerical inconsistencies, use 'exact=True'.")
            from sage.rings.number_field.number_field import QuadraticField
            base_ring = QuadraticField(5, 'sqrt5')
            sqrt5 = base_ring.gen()
            phi = (1 + sqrt5) / 2
            phi_inv = base_ring.one() / phi

            # The 24 permutations of [0,0,±2,±2] (the ± are independent)
            verts = Permutations([0,0,2,2]).list() + Permutations([0,0,-2,-2]).list() + Permutations([0,0,2,-2]).list()

            # The 64 permutations of the following vectors:
            # [±1,±1,±1,±sqrt(5)]
            # [±1/phi^2,±phi,±phi,±phi]
            # [±1/phi,±1/phi,±1/phi,±phi^2]
            from sage.categories.cartesian_product import cartesian_product
            full_perm_vectors = [[[1,-1],[1,-1],[1,-1],[-sqrt5,sqrt5]],
                                 [[phi_inv**2,-phi_inv**2],[phi,-phi],[phi,-phi],[-phi,phi]],
                                 [[phi_inv,-phi_inv],[phi_inv,-phi_inv],[phi_inv,-phi_inv],[-(phi**2),phi**2]]]
            for vect in full_perm_vectors:
                cp = cartesian_product(vect)
                # The group action creates duplicates, so we reduce it:
                verts += list(set([tuple(p) for c in cp for p in Permutations(list(c))]))

            # The 96 even permutations of [0,±1/phi^2,±1,±phi^2]
            # The 96 even permutations of [0,±1/phi,±phi,±sqrt(5)]
            # The 192 even permutations of [±1/phi,±1,±phi,±2]
            even_perm_vectors = [[[0],[phi_inv**2,-phi_inv**2],[1,-1],[-(phi**2),phi**2]],
                                 [[0],[phi_inv,-phi_inv],[phi,-phi],[-sqrt5,sqrt5]],
                                 [[phi_inv,-phi_inv],[1,-1],[phi,-phi],[-2,2]]]
            even_perm = AlternatingGroup(4)
            for vect in even_perm_vectors:
                cp = cartesian_product(vect)
                verts += [p(tuple(c)) for p in even_perm for c in cp]

            return Polyhedron(vertices=verts, base_ring=base_ring, backend=backend)

        elif construction == 'as_permutahedron':
            return self.generalized_permutahedron(['H', 4], point=[0, 0, 0, 1], exact=exact, backend=backend, regular=True)
        else:
            raise ValueError("construction (={}) must be either 'coxeter' or 'as_permutahedron' ".format(construction))

    def hypercube(self, dim, intervals=None, backend=None):
        r"""
        Return a hypercube of the given dimension.

        The ``dim``-dimensional hypercube is by default the convex hull of the
        `2^{\text{dim}}` `\pm 1` vectors of length ``dim``. Alternatively,
        it is the product of ``dim`` line segments given in the ``intervals``.
        For more information see the wikipedia article :wikipedia:`Hypercube`.

        INPUT:

        - ``dim`` -- integer. The dimension of the hypercube.

        - ``intervals`` -- (default = None). It takes the following
          possible inputs:

          - If ``None`` (the default), it returns the `\pm 1`-cube of
            dimension ``dim``.

          - ``'zero_one'`` -- (string). Return the `0/1`-cube.

          - a list of length ``dim``. Its elements are pairs of
            numbers `(a,b)` with `a < b`. The cube will be the product of
            these intervals.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES:

        Create the `\pm 1`-hypercube of dimension 4::

            sage: four_cube = polytopes.hypercube(4)
            sage: four_cube.is_simple()
            True
            sage: four_cube.base_ring()
            Integer Ring
            sage: four_cube.volume()
            16
            sage: four_cube.ehrhart_polynomial()    # optional - latte_int
            16*t^4 + 32*t^3 + 24*t^2 + 8*t + 1

        Return the `0/1`-hypercube of dimension 4::

            sage: z_cube = polytopes.hypercube(4, intervals='zero_one')
            sage: z_cube.vertices()[0]
            A vertex at (1, 0, 1, 1)
            sage: z_cube.is_simple()
            True
            sage: z_cube.base_ring()
            Integer Ring
            sage: z_cube.volume()
            1
            sage: z_cube.ehrhart_polynomial()       # optional - latte_int
            t^4 + 4*t^3 + 6*t^2 + 4*t + 1

        Return the 4-dimensional combinatorial cube that is the product of
        [0,3]^4::

            sage: t_cube = polytopes.hypercube(4, intervals=[[0,3]]*4)

        Checking that t_cube is three times the previous `0/1`-cube::

            sage: t_cube == 3 * z_cube
            True

        TESTS::

            sage: fc = polytopes.hypercube(4,backend='normaliz')   # optional - pynormaliz
            sage: TestSuite(fc).run()                              # optional - pynormaliz

        ::

            sage: ls = [randint(-100,100) for _ in range(4)]
            sage: intervals = [[x, x+randint(1,50)] for x in ls]
            sage: P = polytopes.hypercube(4, intervals, backend='field')
            sage: TestSuite(P).run()

        Check that :trac:`29904` is fixed::

            sage: intervals = [[-2,2]]
            sage: P = polytopes.hypercube(1, intervals, 'field')
            sage: TestSuite(P).run()

        If the dimension ``dim`` is not equal to the length of intervals, an
        error is raised::

            sage: u_cube = polytopes.hypercube(2, intervals=[[0,1],[0,2],[0,3]])
            Traceback (most recent call last):
            ...
            ValueError: the dimension of the hypercube must match the number of intervals

        The intervals must be pairs `(a, b)` with `a < b`::

            sage: w_cube = polytopes.hypercube(3, intervals=[[0,1],[3,2],[0,3]])
            Traceback (most recent call last):
            ...
            ValueError: each interval must be a pair `(a, b)` with `a < b`

        If a string besides 'zero_one' is passed to ``intervals``, return an
        error::

            sage: v_cube = polytopes.hypercube(3, intervals='a_string')
            Traceback (most recent call last):
            ...
            ValueError: the only allowed string is 'zero_one'

        Check that we set up the hypercube correctly::

            sage: ls = [randint(-100,100) for _ in range(4)]
            sage: intervals = [[x, x+randint(1,50)] for x in ls]
            sage: P = polytopes.hypercube(4, intervals, backend='field')
            sage: P1 = polytopes.hypercube(4, intervals, backend='ppl')                 # needs pplpy
            sage: assert P == P1                                                        # needs pplpy

        Check that coercion for input invervals is handled correctly::

            sage: P = polytopes.hypercube(2, [[1/2, 2], [0, 1]])
            sage: P = polytopes.hypercube(2, [[1/2, 2], [0, 1.0]])
            sage: P = polytopes.hypercube(2, [[1/2, 2], [0, AA(2).sqrt()]])             # needs sage.rings.number_field
            sage: P = polytopes.hypercube(2, [[1/2, 2], [0, 1.0]], backend='ppl')       # needs pplpy
            Traceback (most recent call last):
            ...
            ValueError: specified backend ppl cannot handle the intervals
        """
        parent = Polyhedra(ZZ, dim, backend=backend)
        convert = False

        # If the intervals are (a_1,b_1), ..., (a_dim, b_dim),
        # then the inequalites correspond to
        # b_1,b_2,...,b_dim, a_1,a_2,...,a_dim
        # in that order.

        if intervals is None:
            cp = itertools.product((-1,1), repeat=dim)

            # An inequality -x_i       + 1 >= 0 for i <  dim
            # resp.          x_{dim-i} + 1 >= 0 for i >= dim
            ieq_b = lambda i: 1

        elif isinstance(intervals, str):
            if intervals == 'zero_one':
                cp = itertools.product((0,1), repeat=dim)

                # An inequality -x_i       + 1 >= 0 for i <  dim
                # resp.          x_{dim-i} + 0 >= 0 for i >= dim
                ieq_b = lambda i: 1 if i < dim else 0
            else:
                raise ValueError("the only allowed string is 'zero_one'")
        elif len(intervals) == dim:
            if not all(a < b for a,b in intervals):
                raise ValueError("each interval must be a pair `(a, b)` with `a < b`")
            parent = parent.base_extend(sum(a + b for a,b in intervals))
            if parent.base_ring() not in (ZZ, QQ):
                convert = True
            if backend and parent.backend() is not backend:
                # If the parent changed backends, but a backend was specified,
                # the specified backend cannot handle the intervals.
                raise ValueError("specified backend {} cannot handle the intervals".format(backend))

            cp = itertools.product(*intervals)

            # An inequality -x_i       + b_i >= 0 for i <  dim
            # resp.          x_{dim-i} - a_i >= 0 for i >= dim
            ieq_b = lambda i: intervals[i][1] if i < dim \
                              else -intervals[i-dim][0]
        else:
            raise ValueError("the dimension of the hypercube must match the number of intervals")

        # An inequality -x_i       + ieq_b(i)     >= 0 for i <  dim
        # resp.          x_{dim-i} + ieq_b(i-dim) >= 0 for i >= dim
        ieq_A = lambda i, pos: -1 if i == pos           \
                               else 1 if i == pos + dim \
                               else 0
        ieqs = (tuple(ieq_b(i) if pos == 0 else ieq_A(i, pos-1)
                      for pos in range(dim+1))
                for i in range(2*dim))

        return parent([cp, [], []], [ieqs, []], convert=convert, Vrep_minimal=True, Hrep_minimal=True, pref_rep='Hrep')

    def cube(self, intervals=None, backend=None):
        r"""
        Return the cube.

        The cube is the Platonic solid that is obtained as the convex hull of
        the eight `\pm 1` vectors of length 3 (by default). Alternatively, the
        cube is the product of three intervals from ``intervals``.

        .. SEEALSO::

            :meth:`hypercube`

        INPUT:

        - ``intervals`` -- list (default=None). It takes the following
          possible inputs:

          - If the input is ``None`` (the default), returns the convex hull of
            the eight `\pm 1` vectors of length three.

          - ``'zero_one'`` -- (string). Return the `0/1`-cube.

          - a list of 3 lists of length 2. The cube will be a product of
            these three intervals.

        - ``backend`` -- the backend to use to create the polytope.

        OUTPUT:

        A cube as a polyhedron object.

        EXAMPLES:

        Return the `\pm 1`-cube::

            sage: c = polytopes.cube()
            sage: c
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: c.f_vector()
            (1, 8, 12, 6, 1)
            sage: c.volume()
            8
            sage: c.plot()                                                              # needs sage.plot
            Graphics3d Object

        Return the `0/1`-cube::

            sage: cc = polytopes.cube(intervals ='zero_one')
            sage: cc.vertices_list()
            [[1, 0, 0],
             [1, 1, 0],
             [1, 1, 1],
             [1, 0, 1],
             [0, 0, 1],
             [0, 0, 0],
             [0, 1, 0],
             [0, 1, 1]]
        """
        return self.hypercube(3, backend=backend, intervals=intervals)

    def cross_polytope(self, dim, backend=None):
        r"""
        Return a cross-polytope in dimension ``dim``.

        A cross-polytope is a higher dimensional generalization of the
        octahedron. It is the convex hull of the `2d` points `(\pm 1, 0,
        \ldots, 0)`, `(0, \pm 1, \ldots, 0)`, \ldots, `(0, 0, \ldots, \pm 1)`.
        See the :wikipedia:`Cross-polytope` for more information.

        INPUT:

        - ``dim`` -- integer. The dimension of the cross-polytope.

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: four_cross = polytopes.cross_polytope(4)
            sage: four_cross.f_vector()
            (1, 8, 24, 32, 16, 1)
            sage: four_cross.is_simple()
            False

        TESTS::

            sage: cp = polytopes.cross_polytope(4, backend='normaliz')  # optional - pynormaliz
            sage: TestSuite(cp).run()                                   # optional - pynormaliz

        ::

            sage: P = polytopes.cross_polytope(6, backend='field')
            sage: TestSuite(P).run()  # long time

        Check that double description is set up correctly::

            sage: P = polytopes.cross_polytope(6, backend='ppl')
            sage: Q = polytopes.cross_polytope(6, backend='ppl')
            sage: P == Q
            True
        """
        verts = tuple((ZZ**dim).basis())
        verts += tuple(-v for v in verts)
        ieqs = ((1,) + x for x in itertools.product((-1,1), repeat=dim))
        parent = Polyhedra(ZZ, dim, backend=backend)
        return parent([verts, [], []], [ieqs, []], Vrep_minimal=True, Hrep_minimal=True, pref_rep='Vrep')

    def parallelotope(self, generators, backend=None):
        r"""
        Return the zonotope, or parallelotope, spanned by the generators.

        The parallelotope is the multi-dimensional generalization of a
        parallelogram (2 generators) and a parallelepiped (3 generators).

        INPUT:

        - ``generators`` -- a list of vectors of same dimension

        - ``backend`` -- the backend to use to create the polytope.

        EXAMPLES::

            sage: polytopes.parallelotope([ (1,0), (0,1) ])
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: polytopes.parallelotope([[1,2,3,4], [0,1,0,7], [3,1,0,2], [0,0,1,0]])
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 16 vertices

            sage: K = QuadraticField(2, 'sqrt2')                                        # needs sage.rings.number_field
            sage: sqrt2 = K.gen()                                                       # needs sage.rings.number_field
            sage: P = polytopes.parallelotope([(1, sqrt2), (1, -1)]); P                 # needs sage.rings.number_field
            A 2-dimensional polyhedron in (Number Field in sqrt2 with defining
            polynomial x^2 - 2 with sqrt2 = 1.414213562373095?)^2 defined as
            the convex hull of 4 vertices

        TESTS::

            sage: TestSuite(P).run()                                                    # needs sage.rings.number_field
        """
        from sage.modules.free_module_element import vector
        generators = [vector(v) for v in generators]
        if not generators:
            return Polyhedron(backend=backend)

        zero = generators[0] - generators[0]
        intervals = [Polyhedron([zero, gen], backend=backend) for gen in generators]
        return sum(intervals)

    zonotope = parallelotope

    # --------------------------------------------------------
    # imports from other files
    # --------------------------------------------------------
    associahedron = staticmethod(Associahedron)

    try:
        flow_polytope = staticmethod(DiGraph.flow_polytope)
        edge_polytope = staticmethod(Graph.edge_polytope)
        symmetric_edge_polytope = staticmethod(Graph.symmetric_edge_polytope)
    except ImportError:
        pass


polytopes = Polytopes()
