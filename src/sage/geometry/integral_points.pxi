# sage_setup: distribution = sagemath-linbox
r"""
Cython helper methods to compute integral points in polyhedra.
"""

#*****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_check
import copy
import itertools

from sage.matrix.constructor import matrix, column_matrix, vector, diagonal_matrix
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from sage.misc.misc_c import prod
from sage.modules.free_module import FreeModule

##############################################################################
# The basic idea to enumerate the lattice points in the parallelotope
# is to use the Smith normal form of the ray coordinate matrix to
# bring the lattice into a "nice" basis. Here is the straightforward
# implementation. Note that you do not need to reduce to the
# full-dimensional case, the Smith normal form takes care of that for
# you.
#
## def parallelotope_points(spanning_points, lattice):
##     # compute points in the open parallelotope, see [BK2001]
##     R = matrix(spanning_points).transpose()
##     D,U,V = R.smith_form()
##     e = D.diagonal()          # the elementary divisors
##     d = prod(e)               # the determinant
##     u = U.inverse().columns() # generators for gp(semigroup)
##
##     # "inverse" of the ray matrix as far as possible over ZZ
##     # R*Rinv == diagonal_matrix([d]*D.ncols() + [0]*(D.nrows()-D.ncols()))
##     # If R is full rank, this is Rinv = matrix(ZZ, R.inverse() * d)
##     Dinv = D.transpose()
##     for i in range(D.ncols()):
##         Dinv[i,i] = d/D[i,i]
##     Rinv = V * Dinv * U
##
##     gens = []
##     for b in CartesianProduct(*[ range(i) for i in e ]):
##         # this is our generator modulo the lattice spanned by the rays
##         gen_mod_rays = sum( b_i*u_i for b_i, u_i in zip(b,u) )
##         q_times_d = Rinv * gen_mod_rays
##         q_times_d = vector(ZZ,[ q_i % d  for q_i in q_times_d ])
##         gen = lattice(R*q_times_d / d)
##         gen.set_immutable()
##         gens.append(gen)
##     assert(len(gens) == d)
##     return tuple(gens)
#
# The problem with the naive implementation is that it is slow:
#
#   1. You can simplify some of the matrix multiplications
#
#   2. The inner loop keeps creating new matrices and vectors, which
#      is slow. It needs to recycle objects. Instead of creating a new
#      lattice point again and again, change the entries of an
#      existing lattice point and then copy it!


cpdef tuple parallelotope_points(spanning_points, lattice):
    r"""
    Return integral points in the parallelotope starting at the origin
    and spanned by the ``spanning_points``.

    See :meth:`~ConvexRationalPolyhedralCone.semigroup_generators` for a description of the
    algorithm.

    INPUT:

    - ``spanning_points`` -- a non-empty list of linearly independent
      rays (`\ZZ`-vectors or :class:`toric lattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>` elements),
      not necessarily primitive lattice points.

    OUTPUT:

    The tuple of all lattice points in the half-open parallelotope
    spanned by the rays `r_i`,

    .. MATH::

        \mathop{par}(\{r_i\}) =
        \sum_{0\leq a_i < 1} a_i r_i

    By half-open parallelotope, we mean that the
    points in the facets not meeting the origin are omitted.

    EXAMPLES:

    Note how the points on the outward-facing factes are omitted::

        sage: from sage.geometry.integral_points import parallelotope_points
        sage: rays = list(map(vector, [(2,0), (0,2)]))
        sage: parallelotope_points(rays, ZZ^2)
        ((0, 0), (0, 1), (1, 0), (1, 1))

    The rays can also be toric lattice points::

        sage: rays = list(map(ToricLattice(2), [(2,0), (0,2)]))
        sage: parallelotope_points(rays, ToricLattice(2))
        (N(0, 0), N(0, 1), N(1, 0), N(1, 1))

    A non-smooth cone::

        sage: c = Cone([ (1,0), (1,2) ])
        sage: parallelotope_points(c.rays(), c.lattice())
        (N(0, 0), N(1, 1))

    A :exc:`ValueError` is raised if the ``spanning_points`` are not
    linearly independent::

        sage: rays = list(map(ToricLattice(2), [(1,1)]*2))
        sage: parallelotope_points(rays, ToricLattice(2))
        Traceback (most recent call last):
        ...
        ValueError: The spanning points are not linearly independent!

    TESTS::

        sage: rays = list(map(vector,[(-3, -2, -3, -2), (-2, -1, -8, 5), (1, 9, -7, -4), (-3, -1, -2, 2)]))
        sage: len(parallelotope_points(rays, ZZ^4))
        967
    """
    cdef MatrixClass VDinv
    cdef MatrixClass R = matrix(spanning_points).transpose()
    e, d, VDinv = ray_matrix_normal_form(R)
    cdef tuple points = loop_over_parallelotope_points(e, d, VDinv, R, lattice)
    for p in points:
        p.set_immutable()
    return points


cpdef tuple ray_matrix_normal_form(R):
    r"""
    Compute the Smith normal form of the ray matrix for
    :func:`parallelotope_points`.

    INPUT:

    - ``R`` -- `\ZZ`-matrix whose columns are the rays spanning the
      parallelotope

    OUTPUT: a tuple containing ``e``, ``d``, and ``VDinv``

    EXAMPLES::

        sage: from sage.geometry.integral_points import ray_matrix_normal_form
        sage: R = column_matrix(ZZ,[3,3,3])
        sage: ray_matrix_normal_form(R)
        ([3], 3, [1])
    """
    D,U,V = R.smith_form()
    e = D.diagonal()            # the elementary divisors
    cdef Integer d = prod(e)                 # the determinant
    if d == ZZ.zero():
        raise ValueError('The spanning points are not linearly independent!')
    cdef int i
    Dinv = diagonal_matrix(ZZ, [ d // e[i] for i in range(D.ncols()) ])
    VDinv = V * Dinv
    return (e, d, VDinv)


# The optimized version avoids constructing new matrices, vectors, and lattice points
cpdef tuple loop_over_parallelotope_points(e, d, MatrixClass VDinv,
                                           MatrixClass R, lattice,
                                           A=None, b=None):
    r"""
    The inner loop of :func:`parallelotope_points`.

    INPUT:

    See :meth:`parallelotope_points` for ``e``, ``d``, ``VDinv``, ``R``, ``lattice``.

    - ``A``, ``b`` -- either both ``None`` or a vector and number. If
      present, only the parallelotope points satisfying `A x \leq b`
      are returned.

    OUTPUT:

    The points of the half-open parallelotope as a tuple of lattice
    points.

    EXAMPLES::

        sage: e = [3]
        sage: d = prod(e)
        sage: VDinv = matrix(ZZ, [[1]])
        sage: R = column_matrix(ZZ,[3,3,3])
        sage: lattice = ZZ^3
        sage: from sage.geometry.integral_points import loop_over_parallelotope_points
        sage: loop_over_parallelotope_points(e, d, VDinv, R, lattice)
        ((0, 0, 0), (1, 1, 1), (2, 2, 2))

        sage: A = vector(ZZ, [1,0,0])
        sage: b = 1
        sage: loop_over_parallelotope_points(e, d, VDinv, R, lattice, A, b)
        ((0, 0, 0), (1, 1, 1))
    """
    cdef int i, j
    cdef int dim = VDinv.nrows()
    cdef int ambient_dim = R.nrows()
    s = ZZ.zero() # summation variable
    cdef list gens = []
    gen = lattice(ZZ.zero())
    cdef VectorClass q_times_d = vector(ZZ, dim)
    for base in itertools.product(*[ range(i) for i in e ]):
        for i in range(dim):
            s = ZZ.zero()
            for j in range(dim):
                s += VDinv.get_unsafe(i,j) * base[j]
            q_times_d.set_unsafe(i, s % d)
        for i in range(ambient_dim):
            s = ZZ.zero()
            for j in range(dim):
                s += R.get_unsafe(i,j) * q_times_d.get_unsafe(j)
            gen[i] = s / d
        if A is not None:
            s = ZZ.zero()
            for i in range(ambient_dim):
                s += A[i] * gen[i]
            if s > b:
                continue
        gens.append(copy.copy(gen))
    if A is None:
        assert(len(gens) == d)
    return tuple(gens)


##############################################################################
cpdef tuple simplex_points(vertices):
    r"""
    Return the integral points in a lattice simplex.

    INPUT:

    - ``vertices`` -- an iterable of integer coordinate vectors. The
      indices of vertices that span the simplex under
      consideration.

    OUTPUT:

    A tuple containing the integral point coordinates as `\ZZ`-vectors.

    EXAMPLES::

        sage: from sage.geometry.integral_points import simplex_points
        sage: simplex_points([(1,2,3), (2,3,7), (-2,-3,-11)])
        ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

    The simplex need not be full-dimensional::

        sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)])
        sage: simplex_points(simplex.Vrepresentation())
        ((2, 3, 7, 5), (0, 0, -2, 5), (-2, -3, -11, 5), (1, 2, 3, 5))

        sage: simplex_points([(2,3,7)])
        ((2, 3, 7),)

    TESTS::

        sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
        sage: simplex = Polyhedron(v); simplex
        A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
        sage: pts = simplex_points(simplex.Vrepresentation())
        sage: len(pts)
        49
        sage: for p in pts: p.set_immutable()
        sage: len(set(pts))
        49

        sage: all(simplex.contains(p) for p in pts)
        True

        sage: v = [(4,-1,-1,-1), (-1,4,-1,-1), (-1,-1,4,-1), (-1,-1,-1,4), (-1,-1,-1,-1)]
        sage: P4mirror = Polyhedron(v); P4mirror
        A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
        sage: len(simplex_points(P4mirror.Vrepresentation()))
        126

        sage: vertices = list(map(vector, [(1,2,3), (2,3,7), (-2,-3,-11)]))
        sage: for v in vertices: v.set_immutable()
        sage: simplex_points(vertices)
        ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))
    """
    cdef list rays = [vector(ZZ, list(v)) for v in vertices]
    if not rays:
        return ()
    cdef VectorClass origin = rays.pop()
    origin.set_immutable()
    if not rays:
        return (origin,)
    translate_points(rays, origin)

    # Find equation Ax<=b that cuts out simplex from parallelotope
    cdef MatrixClass Rt = matrix(ZZ, rays)
    cdef MatrixClass R = Rt.transpose()
    cdef int nrays = len(rays)
    cdef Integer b
    M = FreeModule(ZZ, nrays)
    if R.is_square():
        b = abs(R.det())
        A = R.solve_left(M([b]*nrays))
    else:
        RtR = Rt * R
        b = abs(RtR.det())
        A = RtR.solve_left(M([b]*nrays)) * Rt

    e, d, VDinv = ray_matrix_normal_form(R)
    lattice = origin.parent()
    points = loop_over_parallelotope_points(e, d, VDinv, R, lattice, A, b) + tuple(rays)
    translate_points(points, -origin)
    return points


cdef translate_points(v_list, VectorClass delta):
    r"""
    Add ``delta`` to each vector in ``v_list``.
    """
    cdef int dim = delta.degree()
    cdef int i
    for v in v_list:
        for i in range(dim):
            v[i] -= delta.get_unsafe(i)


##############################################################################
# For points with "small" coordinates (that is, fitting into a small
# rectangular bounding box) it is faster to naively enumerate the
# points. This saves the overhead of triangulating the polytope etc.

cpdef rectangular_box_points(list box_min, list box_max,
                             polyhedron=None, count_only=False,
                             return_saturated=False):
    r"""
    Return the integral points in the lattice bounding box that are
    also contained in the given polyhedron.

    INPUT:

    - ``box_min`` -- list of integers; the minimal value for each
      coordinate of the rectangular bounding box

    - ``box_max`` -- list of integers; the maximal value for each
      coordinate of the rectangular bounding box

    - ``polyhedron`` -- a
      :class:`~sage.geometry.polyhedron.base.Polyhedron_base`, a PPL
      :class:`~ppl.polyhedron.C_Polyhedron`, or ``None`` (default)

    - ``count_only`` -- boolean (default: ``False``); whether to
      return only the total number of vertices, and not their
      coordinates. Enabling this option speeds up the
      enumeration. Cannot be combined with the ``return_saturated``
      option.

    - ``return_saturated`` -- boolean (default: ``False``); whether to
      also return which inequalities are saturated for each point of
      the polyhedron. Enabling this slows down the enumeration. Cannot
      be combined with the ``count_only`` option.

    OUTPUT:

    By default, this function returns a tuple containing the integral
    points of the rectangular box spanned by ``box_min`` and ``box_max``
    and that lie inside the ``polyhedron``. For sufficiently large
    bounding boxes, this are all integral points of the polyhedron.

    If no polyhedron is specified, all integral points of the
    rectangular box are returned.

    If ``count_only`` is specified, only the total number (an integer)
    of found lattice points is returned.

    If ``return_saturated`` is enabled, then for each integral point a
    pair ``(point, Hrep)`` is returned where ``point`` is the point
    and ``Hrep`` is the set of indices of the H-representation objects
    that are saturated at the point.

    ALGORITHM:

    This function implements the naive algorithm towards counting
    integral points. Given min and max of vertex coordinates, it
    iterates over all points in the bounding box and checks whether
    they lie in the polyhedron. The following optimizations are
    implemented:

    * Cython: Use machine integers and optimizing C/C++ compiler
      where possible, arbitrary precision integers where necessary.
      Bounds checking, no compile time limits.

    * Unwind inner loop (and next-to-inner loop):

      .. MATH::

          Ax \leq b
          \quad \Leftrightarrow \quad
          a_1 x_1 ~\leq~ b - \sum_{i=2}^d a_i x_i

      so we only have to evaluate `a_1 * x_1` in the inner loop.

    * Coordinates are permuted to make the longest box edge the
      inner loop. The inner loop is optimized to run very fast, so
      its best to do as much work as possible there.

    * Continuously reorder inequalities and test the most
      restrictive inequalities first.

    * Use convexity and only find first and last allowed point in
      the inner loop. The points in-between must be points of the
      polyhedron, too.

    EXAMPLES::

        sage: from sage.geometry.integral_points import rectangular_box_points
        sage: rectangular_box_points([0,0,0],[1,2,3])
        ((0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3),
         (0, 1, 0), (0, 1, 1), (0, 1, 2), (0, 1, 3),
         (0, 2, 0), (0, 2, 1), (0, 2, 2), (0, 2, 3),
         (1, 0, 0), (1, 0, 1), (1, 0, 2), (1, 0, 3),
         (1, 1, 0), (1, 1, 1), (1, 1, 2), (1, 1, 3),
         (1, 2, 0), (1, 2, 1), (1, 2, 2), (1, 2, 3))

        sage: from sage.geometry.integral_points import rectangular_box_points
        sage: rectangular_box_points([0,0,0],[1,2,3], count_only=True)
        24

        sage: cell24 = polytopes.twenty_four_cell()
        sage: rectangular_box_points([-1]*4, [1]*4, cell24)
        ((-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1),
         (0, 0, 0, 0),
         (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, 0))
        sage: d = 3
        sage: dilated_cell24 = d*cell24
        sage: len( rectangular_box_points([-d]*4, [d]*4, dilated_cell24) )
        305

        sage: d = 6
        sage: dilated_cell24 = d*cell24
        sage: len( rectangular_box_points([-d]*4, [d]*4, dilated_cell24) )
        3625

        sage: rectangular_box_points([-d]*4, [d]*4, dilated_cell24, count_only=True)
        3625

        sage: polytope = Polyhedron([(-4,-3,-2,-1),(3,1,1,1),(1,2,1,1),(1,1,3,0),(1,3,2,4)])
        sage: pts = rectangular_box_points([-4]*4, [4]*4, polytope); pts
        ((-4, -3, -2, -1), (-1, 0, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1), (1, 1, 3, 0),
         (1, 2, 1, 1), (1, 2, 2, 2), (1, 3, 2, 4), (2, 1, 1, 1), (3, 1, 1, 1))
        sage: all(polytope.contains(p) for p in pts)
        True

        sage: set(map(tuple,pts)) == \
        ....: set([(-4,-3,-2,-1),(3,1,1,1),(1,2,1,1),(1,1,3,0),(1,3,2,4),
        ....:      (0,1,1,1),(1,2,2,2),(-1,0,0,1),(1,1,1,1),(2,1,1,1)])   # computed with PALP
        True

    Long ints and non-integral polyhedra are explicitly allowed::

        sage: polytope = Polyhedron([[1], [10*pi.n()]], base_ring=RDF)                  # needs sage.symbolic
        sage: len(rectangular_box_points([-100], [100], polytope))                      # needs sage.symbolic
        31

        sage: halfplane = Polyhedron(ieqs=[(-1,1,0)])
        sage: rectangular_box_points([0,-1+10^50], [0,1+10^50])
        ((0, 99999999999999999999999999999999999999999999999999),
         (0, 100000000000000000000000000000000000000000000000000),
         (0, 100000000000000000000000000000000000000000000000001))
        sage: len( rectangular_box_points([0,-100+10^50], [1,100+10^50], halfplane) )
        201

    Using a PPL polyhedron::

        sage: # needs pplpy
        sage: from ppl import Variable, Generator_System, C_Polyhedron, point
        sage: gs = Generator_System()
        sage: x = Variable(0); y = Variable(1); z = Variable(2)
        sage: gs.insert(point(0*x + 1*y + 0*z))
        sage: gs.insert(point(0*x + 1*y + 3*z))
        sage: gs.insert(point(3*x + 1*y + 0*z))
        sage: gs.insert(point(3*x + 1*y + 3*z))
        sage: poly = C_Polyhedron(gs)
        sage: rectangular_box_points([0]*3, [3]*3, poly)
        ((0, 1, 0), (0, 1, 1), (0, 1, 2), (0, 1, 3), (1, 1, 0), (1, 1, 1), (1, 1, 2), (1, 1, 3),
         (2, 1, 0), (2, 1, 1), (2, 1, 2), (2, 1, 3), (3, 1, 0), (3, 1, 1), (3, 1, 2), (3, 1, 3))

    Optionally, return the information about the saturated inequalities as well::

        sage: cube = polytopes.cube()
        sage: cube.Hrepresentation(0)
        An inequality (-1, 0, 0) x + 1 >= 0
        sage: cube.Hrepresentation(1)
        An inequality (0, -1, 0) x + 1 >= 0
        sage: cube.Hrepresentation(2)
        An inequality (0, 0, -1) x + 1 >= 0
        sage: rectangular_box_points([0]*3, [1]*3, cube, return_saturated=True)
        (((0, 0, 0), frozenset()),
         ((0, 0, 1), frozenset({2})),
         ((0, 1, 0), frozenset({1})),
         ((0, 1, 1), frozenset({1, 2})),
         ((1, 0, 0), frozenset({0})),
         ((1, 0, 1), frozenset({0, 2})),
         ((1, 1, 0), frozenset({0, 1})),
         ((1, 1, 1), frozenset({0, 1, 2})))

    TESTS:

    Check that this can be interrupted, see :issue:`20781`::

        sage: ieqs = [(-1, -1, -1, -1, -1, -1, -1, -1, -1),
        ....:         (0, -1, 0, 0, 0, 0, 0, 0, 0),
        ....:         (0, -1, 0, 2, -1, 0, 0, 0, 0),
        ....:         (0, 0, -1, -1, 2, -1, 0, 0, 0),
        ....:         (0, 2, 0, -1, 0, 0, 0, 0, 0),
        ....:         (0, 0, 0, 0, 0, 0, 0, -1, 2),
        ....:         (1, 0, 2, 0, -1, 0, 0, 0, 0),
        ....:         (0, 0, 0, 0, -1, 2, -1, 0, 0),
        ....:         (0, 0, 0, 0, 0, 0, 0, 0, -1),
        ....:         (0, 0, 0, 0, 0, -1, 2, -1, 0),
        ....:         (0, 0, 0, 0, 0, 0, -1, 2, -1)]
        sage: P = Polyhedron(ieqs=ieqs)
        sage: alarm(0.5); P.integral_points()
        Traceback (most recent call last):
        ...
        AlarmInterrupt
    """
    assert len(box_min) == len(box_max)
    assert not (count_only and return_saturated)
    cdef int d = len(box_min)
    cdef int i, j
    cdef list diameter = sorted([(box_max[i]-box_min[i], i) for i in range(d)],
                                reverse=True)
    cdef list diameter_value = [x[0] for x in diameter]
    cdef list diameter_index = [x[1] for x in diameter]

    # Construct the inverse permutation
    cdef list orig_perm = list(range(len(diameter_index)))
    for i, j in enumerate(diameter_index):
        orig_perm[j] = i

    box_min = perm_action(diameter_index, box_min)
    box_max = perm_action(diameter_index, box_max)
    cdef InequalityCollection inequalities = InequalityCollection(polyhedron, diameter_index, box_min, box_max)

    if count_only:
        return loop_over_rectangular_box_points(box_min, box_max, inequalities, d, count_only)

    cdef list points = []
    cdef VectorClass v = vector(ZZ, d)
    if not return_saturated:
        for p in loop_over_rectangular_box_points(box_min, box_max, inequalities, d, count_only):
            #  v = vector(ZZ, perm_action(orig_perm, p))   # too slow
            for i in range(d):
                v.set_unsafe(i, Integer(p[orig_perm[i]]))
            v_copy = copy.copy(v)
            v_copy.set_immutable()
            points.append(v_copy)
    else:
        for p, saturated in loop_over_rectangular_box_points_saturated(box_min, box_max,
                                                                       inequalities, d):
            for i in range(d):
                v.set_unsafe(i, Integer(p[orig_perm[i]]))
            v_copy = copy.copy(v)
            v_copy.set_immutable()
            points.append( (v_copy, saturated) )

    return tuple(points)

cdef list perm_action(list p, list lst):
    """
    Return the action of a permutation ``p`` of `(0, ..., n-1)`
    on a list of length `n`.
    """
    return [lst[i] for i in p]

cdef loop_over_rectangular_box_points(list box_min, list box_max,
                                      InequalityCollection inequalities,
                                      int d, bint count_only):
    """
    The inner loop of :func:`rectangular_box_points`.

    INPUT:

    - ``box_min``, ``box_max`` -- the bounding box

    - ``inequalities`` -- a :class:`InequalityCollection` containing
      the inequalities defining the polyhedron

    - ``d`` -- the ambient space dimension

    - ``count_only`` -- whether to only return the total number of
      lattice points

    OUTPUT: the integral points in the bounding box satisfying all
    inequalities
    """
    cdef int inc
    cdef Integer i_min, i_max
    if count_only:
        points = 0
    else:
        points = []
    cdef list p = list(box_min)
    inequalities.prepare_next_to_inner_loop(p)
    while True:
        sig_check()
        inequalities.prepare_inner_loop(p)
        i_min = box_min[0]
        i_max = box_max[0]
        # Find the lower bound for the allowed region
        while i_min <= i_max:
            if inequalities.are_satisfied(i_min):
                break
            i_min += 1
        # Find the upper bound for the allowed region
        while i_min <= i_max:
            if inequalities.are_satisfied(i_max):
                break
            i_max -= 1
        # The points i_min .. i_max are contained in the polyhedron
        if count_only:
            if i_max >= i_min:
                points += i_max - i_min + 1
        else:
            i = i_min
            while i <= i_max:
                p[0] = i
                points.append(tuple(p))
                i += 1
        # finally increment the other entries in p to move on to next inner loop
        inc = 1
        if d == 1:
            return points
        while True:
            if p[inc] == box_max[inc]:
                p[inc] = box_min[inc]
                inc += 1
                if inc == d:
                    return points
            else:
                p[inc] += 1
                break
        if inc > 1:
            inequalities.prepare_next_to_inner_loop(p)


cdef loop_over_rectangular_box_points_saturated(list box_min, list box_max,
                                                InequalityCollection inequalities,
                                                int d):
    """
    The analog of :func:`rectangular_box_points` except that it keeps
    track of which inequalities are saturated.

    INPUT:

    See :func:`rectangular_box_points`.

    OUTPUT:

    The integral points in the bounding box satisfying all
    inequalities, each point being returned by a coordinate vector and
    a set of H-representation object indices.
    """
    cdef int inc
    cdef list points = []
    cdef list p = list(box_min)
    inequalities.prepare_next_to_inner_loop(p)
    while True:
        inequalities.prepare_inner_loop(p)
        i_min = box_min[0]
        i_max = box_max[0]
        # Find the lower bound for the allowed region
        while i_min <= i_max:
            if inequalities.are_satisfied(i_min):
                break
            i_min += 1
        # Find the upper bound for the allowed region
        while i_min <= i_max:
            if inequalities.are_satisfied(i_max):
                break
            i_max -= 1
        # The points i_min .. i_max are contained in the polyhedron
        i = i_min
        while i <= i_max:
            p[0] = i
            saturated = inequalities.satisfied_as_equalities(i)
            points.append( (tuple(p), saturated) )
            i += 1
        # finally increment the other entries in p to move on to next inner loop
        inc = 1
        if d == 1:
            return points
        while True:
            if p[inc] == box_max[inc]:
                p[inc] = box_min[inc]
                inc += 1
                if inc == d:
                    return points
            else:
                p[inc] += 1
                break
        if inc > 1:
            inequalities.prepare_next_to_inner_loop(p)


cdef class Inequality_generic:
    r"""
    An inequality whose coefficients are arbitrary Python/Sage objects

    INPUT:

    - ``A`` -- list of coefficients

    - ``b`` -- element

    OUTPUT: inequality `A x + b \geq 0`

    EXAMPLES::

        sage: from sage.geometry.integral_points import Inequality_generic
        sage: Inequality_generic([2 * pi, sqrt(3), 7/2], -5.5)                          # needs sage.symbolic
        generic: (2*pi, sqrt(3), 7/2) x + -5.50000000000000 >= 0
    """

    cdef list A
    cdef object b
    cdef object coeff
    cdef object cache
    # The index of the inequality in the polyhedron H-representation
    cdef int index

    def __cinit__(self, list A, b, int index=-1):
        """
        The Cython constructor.

        INPUT:

        See :class:`Inequality_generic`.

        EXAMPLES::

            sage: from sage.geometry.integral_points import Inequality_generic
            sage: Inequality_generic([2 * pi, sqrt(3), 7/2], -5.5)                      # needs sage.symbolic
            generic: (2*pi, sqrt(3), 7/2) x + -5.50000000000000 >= 0
        """
        self.A = A
        self.b = b
        self.coeff = 0
        self.cache = 0
        self.index = int(index)

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.geometry.integral_points import Inequality_generic
            sage: Inequality_generic([2,3,7], -5).__repr__()
            'generic: (2, 3, 7) x + -5 >= 0'
        """
        s = 'generic: ('
        s += ', '.join(str(self.A[i]) for i in range(len(self.A)))
        s += ') x + ' + str(self.b) + ' >= 0'
        return s

    cdef prepare_next_to_inner_loop(self, p):
        """
        In :class:`Inequality_int` this method is used to peel of the
        next-to-inner loop.

        See :meth:`InequalityCollection.prepare_inner_loop` for more details.
        """
        pass

    cdef prepare_inner_loop(self, p):
        """
        Peel off the inner loop.

        See :meth:`InequalityCollection.prepare_inner_loop` for more details.
        """
        cdef int j
        self.coeff = self.A[0]
        self.cache = self.b
        for j in range(1, len(self.A)):
            self.cache += self.A[j] * p[j]

    cdef bint is_not_satisfied(self, inner_loop_variable) except -1:
        r"""
        Test the inequality, using the cached value from :meth:`prepare_inner_loop`

        OUTPUT: boolean; whether the inequality is not satisfied
        """
        return inner_loop_variable * self.coeff + self.cache < 0

    cdef bint is_equality(Inequality_generic self, int inner_loop_variable) except -1:
        r"""
        Test the inequality, using the cached value from :meth:`prepare_inner_loop`

        OUTPUT:

        boolean. Given the inequality `Ax + b \geq 0`, this method
        returns whether the equality `Ax + b = 0` is satisfied.
        """
        return inner_loop_variable * self.coeff + self.cache == 0


# if dim>20 then we always use the generic inequalities (Inequality_generic)
DEF INEQ_INT_MAX_DIM = 20

cdef class Inequality_int:
    r"""
    Fast version of inequality in the case that all coefficients fit
    into machine ints.

    INPUT:

    - ``A`` -- list of integers

    - ``b`` -- integer

    - ``max_abs_coordinates`` -- the maximum of the coordinates that
      one wants to evaluate the coordinates on; used for overflow
      checking

    OUTPUT:

    Inequality `A x + b \geq 0`. A :exc:`OverflowError` is raised if a
    machine integer is not long enough to hold the results. A
    :exc:`ValueError` is raised if some of the input is not integral.

    EXAMPLES::

        sage: from sage.geometry.integral_points import Inequality_int
        sage: Inequality_int([2,3,7], -5, [10]*3)
        integer: (2, 3, 7) x + -5 >= 0

        sage: Inequality_int([1]*21, -5, [10]*21)
        Traceback (most recent call last):
        ...
        OverflowError: Dimension limit exceeded.

        sage: Inequality_int([2,3/2,7], -5, [10]*3)
        Traceback (most recent call last):
        ...
        ValueError: Not integral.

        sage: Inequality_int([2,3,7], -5.2, [10]*3)
        Traceback (most recent call last):
        ...
        ValueError: Not integral.

        sage: Inequality_int([2,3,7], -5*10^50, [10]*3)  # actual error message can differ between 32 and 64 bit
        Traceback (most recent call last):
        ...
        OverflowError: ...

    TESTS:

    Check that :issue:`21993` is fixed::

        sage: Inequality_int([18560500, -89466500], 108027, [178933, 37121])
        Traceback (most recent call last):
        ...
        OverflowError: ...
    """
    cdef int A[INEQ_INT_MAX_DIM]
    cdef int b
    cdef int dim
    # the innermost coefficient
    cdef int coeff
    cdef int cache
    # the next-to-innermost coefficient
    cdef int coeff_next
    cdef int cache_next
    # The index of the inequality in the polyhedron H-representation
    cdef int index

    cdef int _to_int(self, x) except? -999:
        if x not in ZZ:
            raise ValueError('Not integral.')
        return int(x)  # raises OverflowError in Cython if necessary

    def __cinit__(self, list A, b, list max_abs_coordinates, int index=-1):
        """
        The Cython constructor.

        See :class:`Inequality_int` for input.

        EXAMPLES::

            sage: from sage.geometry.integral_points import Inequality_int
            sage: Inequality_int([2,3,7], -5, [10]*3)
            integer: (2, 3, 7) x + -5 >= 0
        """
        cdef int i
        self.dim = self._to_int(len(A))
        self.index = int(index)
        if self.dim < 1 or self.dim > INEQ_INT_MAX_DIM:
            raise OverflowError('Dimension limit exceeded.')
        for i in range(self.dim):
            self.A[i] = self._to_int(A[i])
        self.b = self._to_int(b)
        self.coeff = self.A[0]
        if self.dim > 0:
            self.coeff_next = self.A[1]
        # finally, make sure that there cannot be any overflow during the enumeration
        self._to_int(abs(ZZ(b)) + sum( abs(ZZ(A[i])) * ZZ(max_abs_coordinates[i])
                                       for i in range(self.dim) ))

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.geometry.integral_points import Inequality_int
            sage: Inequality_int([2,3,7], -5, [10]*3).__repr__()
            'integer: (2, 3, 7) x + -5 >= 0'
        """
        s = 'integer: ('
        s += ', '.join(str(self.A[i]) for i in range(self.dim))
        s += ') x + ' + str(self.b) + ' >= 0'
        return s

    cdef prepare_next_to_inner_loop(Inequality_int self, p):
        """
        Peel off the next-to-inner loop.

        See :meth:`InequalityCollection.prepare_inner_loop` for more details.
        """
        cdef int j
        self.cache_next = self.b
        for j in range(2, self.dim):
            self.cache_next += self.A[j] * p[j]

    cdef prepare_inner_loop(Inequality_int self, p):
        """
        Peel off the inner loop.

        See :meth:`InequalityCollection.prepare_inner_loop` for more details.
        """
        cdef int j
        if self.dim > 1:
            self.cache = self.cache_next + self.coeff_next * p[1]
        else:
            self.cache = self.cache_next

    cdef bint is_not_satisfied(Inequality_int self, int inner_loop_variable) noexcept:
        return inner_loop_variable * self.coeff + self.cache < 0

    cdef bint is_equality(Inequality_int self, int inner_loop_variable) noexcept:
        return inner_loop_variable * self.coeff + self.cache == 0


cdef class InequalityCollection:
    """
    A collection of inequalities.

    INPUT:

    - ``polyhedron`` -- a polyhedron defining the inequalities

    - ``permutation`` -- list; a 0-based permutation of the coordinates
      Will be used to permute the coordinates of the inequality

    - ``box_min``, ``box_max`` -- the (not permuted) minimal and maximal
      coordinates of the bounding box; used for bounds checking

    EXAMPLES::

        sage: from sage.geometry.integral_points import InequalityCollection
        sage: P_QQ = Polyhedron(identity_matrix(3).columns() + [(-2, -1,-1)], base_ring=QQ)
        sage: ieq = InequalityCollection(P_QQ, [0,1,2], [0]*3,[1]*3); ieq
        The collection of inequalities
        integer: (3, -2, -2) x + 2 >= 0
        integer: (-1, 4, -1) x + 1 >= 0
        integer: (-1, -1, 4) x + 1 >= 0
        integer: (-1, -1, -1) x + 1 >= 0

        sage: P_RR = Polyhedron(identity_matrix(2).columns() + [(-2.7, -1)], base_ring=RDF)
        sage: InequalityCollection(P_RR, [0,1], [0]*2, [1]*2)
        The collection of inequalities
        integer: (-1, -1) x + 1 >= 0
        generic: (-1.0, 3.7) x + 1.0 >= 0
        generic: (1.0, -1.35) x + 1.35 >= 0

        sage: line = Polyhedron(eqns=[(2,3,7)])
        sage: InequalityCollection(line, [0,1], [0]*2, [1]*2 )
        The collection of inequalities
        integer: (3, 7) x + 2 >= 0
        integer: (-3, -7) x + -2 >= 0

    TESTS::

        sage: TestSuite(ieq).run(skip='_test_pickling')
    """
    cdef list ineqs_int
    cdef list ineqs_generic

    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.geometry.integral_points import InequalityCollection
            sage: line = Polyhedron(eqns=[(2,3,7)])
            sage: InequalityCollection(line, [0,1], [0]*2, [1]*2 ).__repr__()
            'The collection of inequalities\ninteger: (3, 7) x + 2 >= 0\ninteger: (-3, -7) x + -2 >= 0'
        """
        s = 'The collection of inequalities\n'
        for ineq in self.ineqs_int:
            s += str(<Inequality_int>ineq) + '\n'
        for ineq in self.ineqs_generic:
            s += str(<Inequality_generic>ineq) + '\n'
        return s.strip()

    cpdef tuple _make_A_b(self, Hrep_obj, list permutation):
        r"""
        Return the coefficients and constant of the H-representation
        object.

        INPUT:

        - ``Hrep_obj`` -- a H-representation object of the polyhedron
        - ``permutation`` -- the permutation of the coordinates to
          apply to ``A``

        OUTPUT:

        A pair ``(A,b)``.

        EXAMPLES::

            sage: from sage.geometry.integral_points import InequalityCollection
            sage: line = Polyhedron(eqns=[(2,3,7)])
            sage: ieq = InequalityCollection(line, [0,1], [0]*2, [1]*2 )
            sage: ieq._make_A_b(line.Hrepresentation(0), [0,1])
            ([3, 7], 2)
            sage: ieq._make_A_b(line.Hrepresentation(0), [1,0])
            ([7, 3], 2)
        """
        cdef list v = list(Hrep_obj)
        cdef list A = perm_action(permutation, v[1:])
        b = v[0]
        try:
            x = lcm([a.denominator() for a in A] + [b.denominator()])
            A = [a * x for a in A]
            b = b * x
        except AttributeError:
            pass
        return (A, b)

    def __cinit__(self, polyhedron, list permutation, box_min, box_max):
        """
        The Cython constructor.

        See the class documentation for the description of the arguments.

        EXAMPLES::

            sage: from sage.geometry.integral_points import InequalityCollection
            sage: line = Polyhedron(eqns=[(2,3,7)])
            sage: InequalityCollection(line, [0,1], [0]*2, [1]*2 )
            The collection of inequalities
            integer: (3, 7) x + 2 >= 0
            integer: (-3, -7) x + -2 >= 0
        """
        cdef list max_abs_coordinates = [max(abs(c_min), abs(c_max))
                                         for c_min, c_max in zip(box_min, box_max)]
        max_abs_coordinates = perm_action(permutation, max_abs_coordinates)
        self.ineqs_int = []
        self.ineqs_generic = []
        if polyhedron is None:
            return

        try:
            # polyhedron is a PPL C_Polyhedron class?
            self._cinit_from_PPL(max_abs_coordinates, permutation, polyhedron)
        except AttributeError:
            try:
                # polyhedron is a Polyhedron class?
                self._cinit_from_Polyhedron(max_abs_coordinates, permutation, polyhedron)
            except AttributeError:
                raise TypeError('Cannot extract Hrepresentation data from polyhedron.')

    cdef _cinit_from_PPL(self, list max_abs_coordinates, list permutation,
                         polyhedron):
        """
        Initialize the inequalities from a PPL C_Polyhedron.

        See __cinit__ for a description of the arguments.

        EXAMPLES::

            sage: # needs pplpy
            sage: from ppl import Variable, Generator_System, C_Polyhedron, point
            sage: gs = Generator_System()
            sage: x = Variable(0); y = Variable(1); z = Variable(2)
            sage: gs.insert(point(0*x + 0*y + 1*z))
            sage: gs.insert(point(0*x + 3*y + 1*z))
            sage: gs.insert(point(3*x + 0*y + 1*z))
            sage: gs.insert(point(3*x + 3*y + 1*z))
            sage: poly = C_Polyhedron(gs)
            sage: from sage.geometry.integral_points import InequalityCollection
            sage: InequalityCollection(poly, [0,2,1], [0]*3, [3]*3 )
            The collection of inequalities
            integer: (0, 1, 0) x + -1 >= 0
            integer: (0, -1, 0) x + 1 >= 0
            integer: (1, 0, 0) x + 0 >= 0
            integer: (0, 0, 1) x + 0 >= 0
            integer: (-1, 0, 0) x + 3 >= 0
            integer: (0, 0, -1) x + 3 >= 0
        """
        cdef list A
        cdef int index
        for index,c in enumerate(polyhedron.minimized_constraints()):
            A = perm_action(permutation, [Integer(mpz) for mpz in c.coefficients()])
            b = Integer(c.inhomogeneous_term())
            try:
                H = Inequality_int(A, b, max_abs_coordinates, index)
                self.ineqs_int.append(H)
            except (OverflowError, ValueError):
                H = Inequality_generic(A, b, index)
                self.ineqs_generic.append(H)
            if c.is_equality():
                A = [ -a for a in A ]
                b = -b
                try:
                    H = Inequality_int(A, b, max_abs_coordinates, index)
                    self.ineqs_int.append(H)
                except (OverflowError, ValueError):
                    H = Inequality_generic(A, b, index)
                    self.ineqs_generic.append(H)

    cdef _cinit_from_Polyhedron(self, list max_abs_coordinates,
                                list permutation, polyhedron):
        """
        Initialize the inequalities from a Sage Polyhedron.

        See __cinit__ for a description of the arguments.

        EXAMPLES::

            sage: from sage.geometry.integral_points import InequalityCollection
            sage: line = Polyhedron(eqns=[(2,3,7)])
            sage: InequalityCollection(line, [0,1], [0]*2, [1]*2 )
            The collection of inequalities
            integer: (3, 7) x + 2 >= 0
            integer: (-3, -7) x + -2 >= 0

        TESTS:

        Check that :issue:`21037` is fixed::

            sage: P = Polyhedron(vertices=((0, 0), (17,3)))
            sage: P += 1/1000*polytopes.regular_polygon(5)                              # needs sage.rings.number_field
            sage: P.integral_points()                                                   # needs sage.rings.number_field
            ((0, 0), (17, 3))
        """
        cdef list A
        for Hrep_obj in polyhedron.inequality_generator():
            A, b = self._make_A_b(Hrep_obj, permutation)
            try:
                H = Inequality_int(A, b, max_abs_coordinates, Hrep_obj.index())
                self.ineqs_int.append(H)
            except (OverflowError, ValueError, TypeError):
                H = Inequality_generic(A, b, Hrep_obj.index())
                self.ineqs_generic.append(H)
        for Hrep_obj in polyhedron.equation_generator():
            A, b = self._make_A_b(Hrep_obj, permutation)
            # add inequality
            try:
                H = Inequality_int(A, b, max_abs_coordinates, Hrep_obj.index())
                self.ineqs_int.append(H)
            except (OverflowError, ValueError, TypeError):
                H = Inequality_generic(A, b, Hrep_obj.index())
                self.ineqs_generic.append(H)
            # add sign-reversed inequality
            A = [ -a for a in A ]
            b = -b
            try:
                H = Inequality_int(A, b, max_abs_coordinates, Hrep_obj.index())
                self.ineqs_int.append(H)
            except (OverflowError, ValueError, TypeError):
                H = Inequality_generic(A, b, Hrep_obj.index())
                self.ineqs_generic.append(H)

    cpdef prepare_next_to_inner_loop(self, p):
        r"""
        Peel off the next-to-inner loop.

        In the next-to-inner loop of :func:`rectangular_box_points`,
        we have to repeatedly evaluate `A x-A_0 x_0+b`. To speed up
        computation, we pre-evaluate

        .. MATH::

            c = b + \sum_{i=2} A_i x_i

        and only compute `A x-A_0 x_0+b = A_1 x_1 +c \geq 0` in the
        next-to-inner loop.

        INPUT:

        - ``p`` -- the point coordinates. Only ``p[2:]`` coordinates
          are potentially used by this method

        EXAMPLES::

             sage: from sage.geometry.integral_points import InequalityCollection, print_cache
             sage: P = Polyhedron(ieqs=[(2,3,7,11)])
             sage: ieq = InequalityCollection(P, [0,1,2], [0]*3,[1]*3); ieq
             The collection of inequalities
             integer: (3, 7, 11) x + 2 >= 0
             sage: ieq.prepare_next_to_inner_loop([2,1,3])
             sage: ieq.prepare_inner_loop([2,1,3])
             sage: print_cache(ieq)
             Cached inner loop: 3 * x_0 + 42 >= 0
             Cached next-to-inner loop: 3 * x_0 + 7 * x_1 + 35 >= 0
        """
        for ineq in self.ineqs_int:
            (<Inequality_int>ineq).prepare_next_to_inner_loop(p)
        for ineq in self.ineqs_generic:
            (<Inequality_generic>ineq).prepare_next_to_inner_loop(p)

    cpdef prepare_inner_loop(self, p):
        r"""
        Peel off the inner loop.

        In the inner loop of :func:`rectangular_box_points`, we have
        to repeatedly evaluate `A x+b\geq 0`. To speed up computation, we pre-evaluate

        .. MATH::

            c = A x - A_0 x_0 +b = b + \sum_{i=1} A_i x_i

        and only test `A_0 x_0 +c \geq 0` in the inner loop.

        You must call :meth:`prepare_next_to_inner_loop` before
        calling this method.

        INPUT:

        - ``p`` -- the coordinates of the point to loop over. Only the
          ``p[1:]`` entries are used

        EXAMPLES::

             sage: from sage.geometry.integral_points import InequalityCollection, print_cache
             sage: P = Polyhedron(ieqs=[(2,3,7,11)])
             sage: ieq = InequalityCollection(P, [0,1,2], [0]*3,[1]*3); ieq
             The collection of inequalities
             integer: (3, 7, 11) x + 2 >= 0
             sage: ieq.prepare_next_to_inner_loop([2,1,3])
             sage: ieq.prepare_inner_loop([2,1,3])
             sage: print_cache(ieq)
             Cached inner loop: 3 * x_0 + 42 >= 0
             Cached next-to-inner loop: 3 * x_0 + 7 * x_1 + 35 >= 0
        """
        for ineq in self.ineqs_int:
            (<Inequality_int>ineq).prepare_inner_loop(p)
        for ineq in self.ineqs_generic:
            (<Inequality_generic>ineq).prepare_inner_loop(p)

    cpdef swap_ineq_to_front(self, int i):
        r"""
        Swap the ``i``-th entry of the list to the front of the list of inequalities.

        INPUT:

        - ``i`` -- integer; the :class:`Inequality_int` to swap to the
          beginning of the list of integral inequalities

        EXAMPLES::

            sage: from sage.geometry.integral_points import InequalityCollection
            sage: P_QQ = Polyhedron(identity_matrix(3).columns() + [(-2, -1,-1)], base_ring=QQ)
            sage: iec = InequalityCollection(P_QQ, [0,1,2], [0]*3,[1]*3)
            sage: iec
            The collection of inequalities
            integer: (3, -2, -2) x + 2 >= 0
            integer: (-1, 4, -1) x + 1 >= 0
            integer: (-1, -1, 4) x + 1 >= 0
            integer: (-1, -1, -1) x + 1 >= 0
            sage: iec.swap_ineq_to_front(3)
            sage: iec
            The collection of inequalities
            integer: (-1, -1, -1) x + 1 >= 0
            integer: (3, -2, -2) x + 2 >= 0
            integer: (-1, 4, -1) x + 1 >= 0
            integer: (-1, -1, 4) x + 1 >= 0
        """
        i_th_entry = self.ineqs_int[i]
        cdef int j
        for j in range(i-1,-1,-1):
            self.ineqs_int[j+1] = self.ineqs_int[j]
        self.ineqs_int[0] = i_th_entry

    cpdef bint are_satisfied(self, inner_loop_variable) except -1:
        r"""
        Return whether all inequalities are satisfied.

        You must call :meth:`prepare_inner_loop` before calling this
        method.

        INPUT:

        - ``inner_loop_variable`` -- integer; the 0th coordinate of
          the lattice point

        OUTPUT: boolean; whether the lattice point is in the polyhedron

        EXAMPLES::

            sage: from sage.geometry.integral_points import InequalityCollection
            sage: line = Polyhedron(eqns=[(2,3,7)])
            sage: ieq = InequalityCollection(line, [0,1], [0]*2, [1]*2 )
            sage: ieq.prepare_next_to_inner_loop([3,4])
            sage: ieq.prepare_inner_loop([3,4])
            sage: ieq.are_satisfied(3)
            False
        """
        cdef int i
        for i in range(len(self.ineqs_int)):
            sig_check()
            ineq = self.ineqs_int[i]
            if (<Inequality_int>ineq).is_not_satisfied(inner_loop_variable):
                if i > 0:
                    self.swap_ineq_to_front(i)
                return False
        for i in range(len(self.ineqs_generic)):
            sig_check()
            ineq = self.ineqs_generic[i]
            if (<Inequality_generic>ineq).is_not_satisfied(inner_loop_variable):
                return False
        return True

    cpdef frozenset satisfied_as_equalities(self, inner_loop_variable):
        """
        Return the inequalities (by their index) that are satisfied as
        equalities.

        INPUT:

        - ``inner_loop_variable`` -- integer; the 0th coordinate of
          the lattice point

        OUTPUT:

        A set of integers in ascending order. Each integer is the
        index of a H-representation object of the polyhedron (either a
        inequality or an equation).

        EXAMPLES::

            sage: from sage.geometry.integral_points import InequalityCollection
            sage: quadrant = Polyhedron(rays=[(1,0), (0,1)])
            sage: ieqs = InequalityCollection(quadrant, [0,1], [-1]*2, [1]*2 )
            sage: ieqs.prepare_next_to_inner_loop([-1,0])
            sage: ieqs.prepare_inner_loop([-1,0])
            sage: ieqs.satisfied_as_equalities(-1)
            frozenset({1})
            sage: ieqs.satisfied_as_equalities(0)
            frozenset({0, 1})
            sage: ieqs.satisfied_as_equalities(1)
            frozenset({1})
        """
        cdef int i
        cdef list result = []
        for i in range(len(self.ineqs_int)):
            sig_check()
            ineq = self.ineqs_int[i]
            if (<Inequality_int>ineq).is_equality(inner_loop_variable):
                result.append( (<Inequality_int>ineq).index )
        for i in range(len(self.ineqs_generic)):
            sig_check()
            ineq = self.ineqs_generic[i]
            if (<Inequality_generic>ineq).is_equality(inner_loop_variable):
                result.append( (<Inequality_generic>ineq).index )
        return frozenset(result)


cpdef print_cache(InequalityCollection inequality_collection):
    r"""
    Print the cached values in :class:`Inequality_int` (for
    debugging/doctesting only).

    EXAMPLES::

        sage: from sage.geometry.integral_points import InequalityCollection, print_cache
        sage: P = Polyhedron(ieqs=[(2,3,7)])
        sage: ieq = InequalityCollection(P, [0,1], [0]*2,[1]*2); ieq
        The collection of inequalities
        integer: (3, 7) x + 2 >= 0
        sage: ieq.prepare_next_to_inner_loop([3,5])
        sage: ieq.prepare_inner_loop([3,5])
        sage: print_cache(ieq)
        Cached inner loop: 3 * x_0 + 37 >= 0
        Cached next-to-inner loop: 3 * x_0 + 7 * x_1 + 2 >= 0
    """
    cdef Inequality_int ieq = <Inequality_int>(inequality_collection.ineqs_int[0])
    print('Cached inner loop: ' +
          str(ieq.coeff) + ' * x_0 + ' + str(ieq.cache) + ' >= 0')
    print('Cached next-to-inner loop: ' +
          str(ieq.coeff) + ' * x_0 + ' +
          str(ieq.coeff_next) + ' * x_1 + ' + str(ieq.cache_next) + ' >= 0')
