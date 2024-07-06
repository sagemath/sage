"""
Find isomorphisms between fans
"""


#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import column_matrix, matrix
from sage.geometry.cone import Cone


class FanNotIsomorphicError(Exception):
    """
    Exception to return if there is no fan isomorphism
    """
    pass


def fan_isomorphic_necessary_conditions(fan1, fan2):
    """
    Check necessary (but not sufficient) conditions for the fans to be isomorphic.

    INPUT:

    - ``fan1``, ``fan2`` -- two fans

    OUTPUT: boolean; ``False`` if the two fans cannot be isomorphic. ``True``
    if the two fans may be isomorphic.

    EXAMPLES::

        sage: fan1 = toric_varieties.P2().fan()                                         # needs palp sage.graphs
        sage: fan2 = toric_varieties.dP8().fan()                                        # needs palp sage.graphs
        sage: from sage.geometry.fan_isomorphism import fan_isomorphic_necessary_conditions
        sage: fan_isomorphic_necessary_conditions(fan1, fan2)                           # needs palp sage.graphs
        False
    """
    if fan1.lattice_dim() != fan2.lattice_dim():
        return False
    if fan1.dim() != fan2.dim():
        return False
    if fan1.nrays() != fan2.nrays():
        return False
    if fan1.ngenerating_cones() != fan2.ngenerating_cones():
        return False
    if fan1.is_complete() != fan2.is_complete():
        return False
    return True


def fan_isomorphism_generator(fan1, fan2):
    """
    Iterate over the isomorphisms from ``fan1`` to ``fan2``.

    ALGORITHM:

    The :meth:`sage.geometry.fan.Fan.vertex_graph` of the two fans is
    compared. For each graph isomorphism, we attempt to lift it to an
    actual isomorphism of fans.

    INPUT:

    - ``fan1``, ``fan2`` -- two fans

    OUTPUT:

    Yields the fan isomorphisms as matrices acting from the right on
    rays.

    EXAMPLES::

        sage: fan = toric_varieties.P2().fan()                                          # needs palp sage.graphs
        sage: from sage.geometry.fan_isomorphism import fan_isomorphism_generator
        sage: sorted(fan_isomorphism_generator(fan, fan))                               # needs palp sage.graphs
        [
        [-1 -1]  [-1 -1]  [ 0  1]  [0 1]  [ 1  0]  [1 0]
        [ 0  1], [ 1  0], [-1 -1], [1 0], [-1 -1], [0 1]
        ]
        sage: m1 = matrix([(1, 0), (0, -5), (-3, 4)])
        sage: m2 = matrix([(3, 0), (1, 0), (-2, 1)])
        sage: m1.elementary_divisors() == m2.elementary_divisors() == [1,1,0]
        True
        sage: fan1 = Fan([Cone([m1*vector([23, 14]), m1*vector([   3,100])]),
        ....:             Cone([m1*vector([-1,-14]), m1*vector([-100, -5])])])
        sage: fan2 = Fan([Cone([m2*vector([23, 14]), m2*vector([   3,100])]),
        ....:             Cone([m2*vector([-1,-14]), m2*vector([-100, -5])])])
        sage: sorted(fan_isomorphism_generator(fan1, fan2))                             # needs sage.graphs
        [
        [-12  1 -5]
        [ -4  0 -1]
        [ -5  0 -1]
        ]

        sage: m0 = identity_matrix(ZZ, 2)
        sage: m1 = matrix([(1, 0), (0, -5), (-3, 4)])
        sage: m2 = matrix([(3, 0), (1, 0), (-2, 1)])
        sage: m1.elementary_divisors() == m2.elementary_divisors() == [1,1,0]
        True
        sage: fan0 = Fan([Cone([m0*vector([1,0]), m0*vector([1,1])]),
        ....:             Cone([m0*vector([1,1]), m0*vector([0,1])])])
        sage: fan1 = Fan([Cone([m1*vector([1,0]), m1*vector([1,1])]),
        ....:             Cone([m1*vector([1,1]), m1*vector([0,1])])])
        sage: fan2 = Fan([Cone([m2*vector([1,0]), m2*vector([1,1])]),
        ....:             Cone([m2*vector([1,1]), m2*vector([0,1])])])
        sage: sorted(fan_isomorphism_generator(fan0, fan0))                             # needs sage.graphs
        [
        [0 1]  [1 0]
        [1 0], [0 1]
        ]
        sage: sorted(fan_isomorphism_generator(fan1, fan1))                             # needs sage.graphs
        [
        [ -3 -20  28]  [1 0 0]
        [ -1  -4   7]  [0 1 0]
        [ -1  -5   8], [0 0 1]
        ]
        sage: sorted(fan_isomorphism_generator(fan1, fan2))                             # needs sage.graphs
        [
        [-24  -3   7]  [-12   1  -5]
        [ -7  -1   2]  [ -4   0  -1]
        [ -8  -1   2], [ -5   0  -1]
        ]
        sage: sorted(fan_isomorphism_generator(fan2, fan1))                             # needs sage.graphs
        [
        [  0   1  -1]  [ 0  1 -1]
        [  1 -13   8]  [ 2 -8  1]
        [  0  -5   4], [ 1  0 -3]
        ]
    """
    if not fan_isomorphic_necessary_conditions(fan1, fan2):
        return

    graph1 = fan1.vertex_graph()
    graph2 = fan2.vertex_graph()
    graph_iso = graph1.is_isomorphic(graph2, edge_labels=True, certificate=True)
    if not graph_iso[0]:
        return
    graph_iso = graph_iso[1]

    # Pick a basis of rays in fan1
    max_cone = fan1(fan1.dim())[0]
    fan1_pivot_rays = max_cone.rays()
    fan1_basis = fan1_pivot_rays + fan1.virtual_rays()   # A QQ-basis for N_1
    fan1_pivot_cones = [ fan1.embed(Cone([r])) for r in fan1_pivot_rays ]

    # The fan2 cones as set(set(ray indices))
    fan2_cones = frozenset(
        frozenset(cone.ambient_ray_indices())
        for cone in fan2.generating_cones() )

    # iterate over all graph isomorphisms graph1 -> graph2
    for perm in graph2.automorphism_group(edge_labels=True):
        # find a candidate m that maps fan1_basis to the image rays under the graph isomorphism
        fan2_pivot_cones = [ perm(graph_iso[c]) for c in fan1_pivot_cones ]
        fan2_pivot_rays = fan2.rays([ c.ambient_ray_indices()[0] for c in fan2_pivot_cones  ])
        fan2_basis = fan2_pivot_rays + fan2.virtual_rays()
        try:
            m = matrix(ZZ, fan1_basis).solve_right(matrix(ZZ, fan2_basis))
            m = m.change_ring(ZZ)
        except (ValueError, TypeError):
            continue # no solution

        # check that the candidate m lifts the vertex graph homomorphism
        graph_image_ray_indices = [ perm(graph_iso[c]).ambient_ray_indices()[0] for c in fan1(1) ]
        try:
            matrix_image_ray_indices = [ fan2.rays().index(r*m) for r in fan1.rays() ]
        except ValueError:
            continue
        if graph_image_ray_indices != matrix_image_ray_indices:
            continue

        # check that the candidate m maps generating cone to generating cone
        image_cones = frozenset( # The image(fan1) cones as set(set(integers)
            frozenset(graph_image_ray_indices[i]
                      for i in cone.ambient_ray_indices())
            for cone in fan1.generating_cones() )
        if image_cones == fan2_cones:
            m.set_immutable()
            yield m


def find_isomorphism(fan1, fan2, check=False):
    """
    Find an isomorphism of the two fans.

    INPUT:

    - ``fan1``, ``fan2`` -- two fans

    - ``check`` -- boolean (default: ``False``); passed to the fan
      morphism constructor, see
      :func:`~sage.geometry.fan_morphism.FanMorphism`

    OUTPUT:

    A fan isomorphism. If the fans are not isomorphic, a
    :exc:`FanNotIsomorphicError` is raised.

    EXAMPLES::

        sage: rays = ((1, 1), (0, 1), (-1, -1), (3, 1))
        sage: cones = [(0,1), (1,2), (2,3), (3,0)]
        sage: fan1 = Fan(cones, rays)

        sage: m = matrix([[-2,3],[1,-1]])
        sage: m.det() == -1
        True
        sage: fan2 = Fan(cones, [vector(r)*m for r in rays])

        sage: from sage.geometry.fan_isomorphism import find_isomorphism
        sage: find_isomorphism(fan1, fan2, check=True)                                  # needs sage.graphs
        Fan morphism defined by the matrix
        [-2  3]
        [ 1 -1]
        Domain fan: Rational polyhedral fan in 2-d lattice N
        Codomain fan: Rational polyhedral fan in 2-d lattice N

        sage: find_isomorphism(fan1, toric_varieties.P2().fan())                        # needs palp sage.graphs
        Traceback (most recent call last):
        ...
        FanNotIsomorphicError

        sage: fan1 = Fan(cones=[[1,3,4,5],[0,1,2,3],[2,3,4],[0,1,5]],
        ....:            rays=[(-1,-1,0),(-1,-1,3),(-1,1,-1),(-1,3,-1),(0,2,-1),(1,-1,1)])
        sage: fan2 = Fan(cones=[[0,2,3,5],[0,1,4,5],[0,1,2],[3,4,5]],
        ....:            rays=[(-1,-1,-1),(-1,-1,0),(-1,1,-1),(0,2,-1),(1,-1,1),(3,-1,-1)])
        sage: fan1.is_isomorphic(fan2)                                                  # needs sage.graphs
        True
    """
    generator = fan_isomorphism_generator(fan1, fan2)
    try:
        m = next(generator)
    except StopIteration:
        raise FanNotIsomorphicError

    from sage.geometry.fan_morphism import FanMorphism
    return FanMorphism(m, domain_fan=fan1, codomain=fan2, check=check)


def fan_2d_cyclically_ordered_rays(fan):
    """
    Return the rays of a 2-dimensional ``fan`` in cyclic order.

    INPUT:

    - ``fan`` -- a 2-dimensional fan

    OUTPUT:

    A :class:`~sage.geometry.point_collection.PointCollection`
    containing the rays in one particular cyclic order.

    EXAMPLES::

        sage: rays = ((1, 1), (-1, -1), (-1, 1), (1, -1))
        sage: cones = [(0,2), (2,1), (1,3), (3,0)]
        sage: fan = Fan(cones, rays)
        sage: fan.rays()
        N( 1,  1),
        N(-1, -1),
        N(-1,  1),
        N( 1, -1)
        in 2-d lattice N
        sage: from sage.geometry.fan_isomorphism import fan_2d_cyclically_ordered_rays
        sage: fan_2d_cyclically_ordered_rays(fan)
        N(-1, -1),
        N(-1,  1),
        N( 1,  1),
        N( 1, -1)
        in 2-d lattice N

    TESTS::

        sage: fan = Fan(cones=[], rays=[], lattice=ZZ^2)
        sage: from sage.geometry.fan_isomorphism import fan_2d_cyclically_ordered_rays
        sage: fan_2d_cyclically_ordered_rays(fan)
        Empty collection
        in Ambient free module of rank 2 over the principal ideal domain Integer Ring
    """
    assert fan.lattice_dim() == 2
    import math
    rays = [ (math.atan2(r[0],r[1]), r) for r in fan.rays() ]
    rays = [ r[1] for r in sorted(rays) ]
    from sage.geometry.point_collection import PointCollection
    return PointCollection(rays, fan.lattice())


def fan_2d_echelon_forms(fan):
    """
    Return echelon forms of all cyclically ordered ray matrices.

    Note that the echelon form of the ordered ray matrices are unique
    up to different cyclic orderings.

    INPUT:

    - ``fan`` -- a fan

    OUTPUT:

    A set of matrices. The set of all echelon forms for all different
    cyclic orderings.

    EXAMPLES::

        sage: fan = toric_varieties.P2().fan()                                          # needs palp sage.graphs
        sage: from sage.geometry.fan_isomorphism import fan_2d_echelon_forms
        sage: fan_2d_echelon_forms(fan)                                                 # needs palp sage.graphs
        frozenset({[ 1  0 -1]
                   [ 0  1 -1]})

        sage: fan = toric_varieties.dP7().fan()                                         # needs palp sage.graphs
        sage: sorted(fan_2d_echelon_forms(fan))                                         # needs palp sage.graphs
        [
        [ 1  0 -1 -1  0]  [ 1  0 -1 -1  0]  [ 1  0 -1 -1  1]  [ 1  0 -1  0  1]
        [ 0  1  0 -1 -1], [ 0  1  1  0 -1], [ 0  1  1  0 -1], [ 0  1  0 -1 -1],
        <BLANKLINE>
        [ 1  0 -1  0  1]
        [ 0  1  1 -1 -1]
        ]

    TESTS::

        sage: rays = [(1, 1), (-1, -1), (-1, 1), (1, -1)]
        sage: cones = [(0,2), (2,1), (1,3), (3,0)]
        sage: fan1 = Fan(cones, rays)
        sage: from sage.geometry.fan_isomorphism import fan_2d_echelon_form, fan_2d_echelon_forms
        sage: echelon_forms = fan_2d_echelon_forms(fan1)
        sage: S4 = CyclicPermutationGroup(4)                                            # needs sage.groups
        sage: rays.reverse()
        sage: cones = [(3,1), (1,2), (2,0), (0,3)]
        sage: for i in range(100):                                                      # needs sage.groups
        ....:     m = random_matrix(ZZ,2,2)
        ....:     if abs(det(m)) != 1: continue
        ....:     perm = S4.random_element()
        ....:     perm_cones = [ (perm(c[0]+1)-1, perm(c[1]+1)-1) for c in cones ]
        ....:     perm_rays = [ rays[perm(i+1)-1] for i in range(len(rays)) ]
        ....:     fan2 = Fan(perm_cones, rays=[m*vector(r) for r in perm_rays])
        ....:     assert fan_2d_echelon_form(fan2) in echelon_forms

    The trivial case was fixed in :issue:`18613`::

        sage: fan = Fan([], lattice=ToricLattice(2))
        sage: fan_2d_echelon_forms(fan)
        frozenset({[]})
        sage: parent(list(_)[0])
        Full MatrixSpace of 2 by 0 dense matrices over Integer Ring
    """
    if fan.nrays() == 0:
        return frozenset([fan_2d_echelon_form(fan)])
    rays = list(fan_2d_cyclically_ordered_rays(fan))
    echelon_forms = []
    for i in range(2):
        for j in range(len(rays)):
            echelon_forms.append(column_matrix(rays).echelon_form())
            first = rays.pop(0)
            rays.append(first)
        rays.reverse()
    return frozenset(echelon_forms)


def fan_2d_echelon_form(fan):
    """
    Return echelon form of a cyclically ordered ray matrix.

    INPUT:

    - ``fan`` -- a fan

    OUTPUT:

    A matrix. The echelon form of the rays in one particular cyclic
    order.

    EXAMPLES::

        sage: fan = toric_varieties.P2().fan()                                          # needs palp sage.graphs
        sage: from sage.geometry.fan_isomorphism import fan_2d_echelon_form
        sage: fan_2d_echelon_form(fan)                                                  # needs palp sage.graphs
        [ 1  0 -1]
        [ 0  1 -1]
    """
    ray_matrix = fan_2d_cyclically_ordered_rays(fan).column_matrix()
    return ray_matrix.echelon_form()
