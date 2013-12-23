# -*- coding: utf-8 -*-
"""
Enumerate Points of a Toric Variety

The classes here are not meant to be instatiated manually. Instead,
you should always use the methods of the :class:`point set
<sage.schemes.toric.homset.SchemeHomset_points_toric_field>` of the
variety.

In this module, points are always represented by tuples instead of
Sage's class for points of the toric variety. All Sage library code
must then convert it to proper point objects before returning it to
the user.

EXAMPLES::

    sage: P2 = toric_varieties.P2(base_ring=GF(3))
    sage: point_set = P2.point_set()
    sage: point_set.cardinality()
    13
    sage: iter(point_set).next()
    [0 : 0 : 1]
    sage: list(point_set)[0:5]
    [[0 : 0 : 1], [1 : 0 : 0], [0 : 1 : 0], [0 : 1 : 1], [0 : 1 : 2]]
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



from copy import copy

from sage.combinat.cartesian_product import CartesianProduct
from sage.misc.misc import powerset, prod
from sage.misc.cachefunc import cached_method


class InfinitePointEnumerator(object):
    
    def __init__(self, fan, ring):
        """
        Point enumerator for infinite fields.

        INPUT:
        
        - ``fan`` -- fan of the toric variety.

        - ``ring`` -- infinite base ring over which to enumerate
          points.

        TESTS::

            sage: from sage.schemes.toric.points import InfinitePointEnumerator
            sage: fan = toric_varieties.P2().fan()
            sage: n = InfinitePointEnumerator(fan, QQ)
            sage: ni = iter(n)
            sage: [ni.next() for k in range(10)]
            [(0, 1, 1), (1, 1, 1), (-1, 1, 1), (1/2, 1, 1), (-1/2, 1, 1), 
             (2, 1, 1), (-2, 1, 1), (1/3, 1, 1), (-1/3, 1, 1), (3, 1, 1)]
          
            sage: X = ToricVariety(Fan([], lattice=ZZ^0))
            sage: X.point_set().cardinality()
            1
            sage: X.base_ring().is_finite()
            False
            sage: X.point_set().list()
            ([],)
        """
        self.ring = ring
        self.fan = fan

    def __iter__(self):
        """
        Iterate over the points.

        OUTPUT:

        Iterator over points.

        EXAMPLES::

            sage: from sage.schemes.toric.points import InfinitePointEnumerator
            sage: fan = toric_varieties.P2().fan()
            sage: n = InfinitePointEnumerator(fan, QQ)
            sage: ni = iter(n)
            sage: [ni.next() for k in range(5)]
            [(0, 1, 1), (1, 1, 1), (-1, 1, 1), (1/2, 1, 1), (-1/2, 1, 1)]
        """
        rays = self.fan().rays() + self.fan().virtual_rays()
        n = len(rays)
        if n == 0:
            yield tuple()
        else:
            R = self.ring
            p = [R.one() for k in range(n)]
            for r in R:
                p[0] = r
                yield tuple(p)


class NaiveFinitePointEnumerator(object):
    
    def __init__(self, fan, ring):
        """
        The naive point enumerator.

        This is very slow.
        
        INPUT:
        
        - ``fan`` -- fan of the toric variety.

        - ``ring`` -- finite base ring over which to enumerate points.

        EXAMPLES::

            sage: from sage.schemes.toric.points import NaiveFinitePointEnumerator
            sage: fan = toric_varieties.P2().fan()
            sage: n = NaiveFinitePointEnumerator(fan, GF(3))
            sage: iter(n).next()
            (0, 0, 1)
        """
        assert ring.is_finite()
        self.ring = ring
        self.fan = fan

    @cached_method
    def units(self):
        """
        Return the units in the base field.

        EXAMPLES::

            sage: ne = toric_varieties.P2(base_ring=GF(5)).point_set()._naive_enumerator()
            sage: ne.units()
            (1, 2, 3, 4)
        """
        return tuple(x for x in self.ring if x != 0)
        
    @cached_method
    def roots(self, n):
        """
        Return the n-th roots in the base field

        EXAMPLES::

            sage: ne = toric_varieties.P2(base_ring=GF(5)).point_set()._naive_enumerator()
            sage: ne.roots(2)
            (1, 4)
            sage: ne.roots(3)
            (1,)
            sage: ne.roots(4)
            (1, 2, 3, 4)
        """
        return tuple(x for x in self.ring if x**n == self.ring.one())
        
    def _Chow_group_free(self):
        r"""
        Return the relations coming from the free part of the Chow group

        OUTPUT:

        A tuple containing the elements of $Hom(A_{d-1,\text{free}},
        F^\times)$, including the identity.

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: X.Chow_group().degree(1)
            C3 x Z
            sage: enum = X.point_set()._naive_enumerator()
            sage: enum._Chow_group_free()
            ((1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4), (5, 5, 5), (6, 6, 6))
        """
        units = self.units()
        result = []
        rays = self.fan().rays() + self.fan().virtual_rays()
        ker = rays.matrix().integer_kernel().matrix()
        for phases in CartesianProduct(*([units] * ker.nrows())):
            phases = tuple(prod(mu**exponent for mu, exponent in zip(phases, column))
                           for column in ker.columns())
            result.append(phases)
        return tuple(sorted(result))

    def _Chow_group_torsion(self):
        r"""
        Return the relations coming from the torison part of the Chow group

        OUTPUT:

        A tuple containing the non-identity elements of
        $Hom(A_{d-1,\text{tors}}, F^\times)$

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: X.Chow_group().degree(1)
            C3 x Z
            sage: enum = X.point_set()._naive_enumerator()
            sage: enum._Chow_group_torsion()
            ((1, 2, 4), (1, 4, 2))
        """
        if self.fan.is_smooth():
            return tuple()
        rays = self.fan().rays() + self.fan().virtual_rays()
        image = rays.column_matrix().image()
        torsion = image.saturation().quotient(image)
        result = set()
        for t in torsion:
            t_lift = t.lift()
            for root in self.roots(t.order()):
                phases = tuple(root**exponent for exponent in t_lift)
                result.add(phases)
        result.remove(tuple(self.ring.one() for r in rays))
        return tuple(sorted(result))

    @cached_method
    def rescalings(self):
        """
        Return the rescalings of homogeneous coordinates.

        OUTPUT:

        A tuple containing all points that are equivalent to
        `[1:1:\dots:1]`, the distinguished point of the big torus
        orbit.
        
        EXAMPLES::

            sage: ni = toric_varieties.P2_123(base_ring=GF(5)).point_set()._naive_enumerator()
            sage: ni.rescalings()
            ((1, 1, 1), (1, 4, 4), (4, 2, 3), (4, 3, 2))

            sage: ni = toric_varieties.dP8(base_ring=GF(3)).point_set()._naive_enumerator()
            sage: ni.rescalings()
            ((1, 1, 1, 1), (1, 2, 2, 2), (2, 1, 2, 1), (2, 2, 1, 2))

            sage: ni = toric_varieties.P1xP1(base_ring=GF(3)).point_set()._naive_enumerator()
            sage: ni.rescalings()
            ((1, 1, 1, 1), (1, 1, 2, 2), (2, 2, 1, 1), (2, 2, 2, 2))
        """
        free = self._Chow_group_free()
        tors = self._Chow_group_torsion()
        if len(tors) == 0:  # optimization for smooth fans
            return free
        result = set(free)
        for f, t in CartesianProduct(free, tors):
            phases = tuple(x*y for x, y in zip(f, t))
            result.add(phases)
        return tuple(sorted(result))

    def orbit(self, point):
        """
        Return the orbit of homogeneous coordinates under rescalings.

        OUTPUT:

        The set of all homogeneous coordinates that are equivalent to ``point``.

        EXAMPLES::

            sage: ne = toric_varieties.P2_123(base_ring=GF(7)).point_set()._naive_enumerator()
            sage: sorted(ne.orbit([1, 0, 0]))
            [(1, 0, 0), (2, 0, 0), (4, 0, 0)]
            sage: sorted(ne.orbit([0, 1, 0]))
            [(0, 1, 0), (0, 6, 0)]
            sage: sorted(ne.orbit([0, 0, 1]))
            [(0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4), (0, 0, 5), (0, 0, 6)]
            sage: sorted(ne.orbit([1, 1, 0]))
            [(1, 1, 0), (1, 6, 0), (2, 1, 0), (2, 6, 0), (4, 1, 0), (4, 6, 0)]
        """
        result = set()
        for phases in self.rescalings():
            p = tuple(mu*z for mu, z in zip(point, phases))
            result.add(p)
        return frozenset(result)

    def cone_iter(self):
        """
        Iterate over all cones of the fan

        OUTPUT:

        Iterator over the cones, starting with the high-dimensional
        ones.

        EXAMPLES::

            sage: ne = toric_varieties.dP6(base_ring=GF(11)).point_set()._naive_enumerator()
            sage: for cone in ne.cone_iter(): 
            ....:     print cone.ambient_ray_indices()
            (0, 1)
            (1, 2)
            (2, 3)
            (3, 4)
            (4, 5)
            (0, 5)
            (0,)
            (1,)
            (2,)
            (3,)
            (4,)
            (5,)
            ()
        """
        fan = self.fan
        for d in range(fan.dim(), -1, -1):
            for cone in fan.cones(d):
                yield cone
            
    def coordinate_iter(self):
        """
        Iterate over all distinct homogeneous coordinates.

        This method does NOT identify homogeneous coordinates that are
        equivalent by a homogeneous rescaling.
        
        OUTPUT:

        An iterator over the points.

        EXAMPLES::

            sage: F2 = GF(2)
            sage: ni = toric_varieties.P2(base_ring=F2).point_set()._naive_enumerator()
            sage: list(ni.coordinate_iter())
            [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)]

            sage: ni = toric_varieties.P1xP1(base_ring=F2).point_set()._naive_enumerator()
            sage: list(ni.coordinate_iter())
            [(0, 1, 0, 1), (1, 0, 0, 1), (1, 0, 1, 0),
             (0, 1, 1, 0), (0, 1, 1, 1), (1, 0, 1, 1),
             (1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 1, 1)]

        TESTS::

            sage: V = ToricVariety(Fan([Cone([(1,1)])]), base_ring=GF(3))
            sage: ni = V.point_set()._naive_enumerator()
            sage: list(ni.coordinate_iter())
            [(0, 1), (0, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        """
        units = [x for x in self.ring if x != 0]
        zero = self.ring.zero()
        rays = self.fan.rays() + self.fan.virtual_rays()
        big_torus = [units] * len(rays)
        for cone in self.cone_iter():
            patch = copy(big_torus)
            for i in cone.ambient_ray_indices():
                patch[i] = [zero]
            for p in CartesianProduct(*patch):
                yield tuple(p)

    def __iter__(self):
        """
        Iterate over the distinct points of the toric variety.

        This function does identify orbits under the homogeneous
        rescalings, and returns precisely one representative per
        orbit.

        OUTPUT:

        Iterator over points.

        EXAMPLES:

            sage: ni = toric_varieties.P2(base_ring=GF(2)).point_set()._naive_enumerator()
            sage: list(ni)
            [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)]

            sage: ni = toric_varieties.P1xP1(base_ring=GF(3)).point_set()._naive_enumerator()
            sage: list(ni) 
            [(0, 1, 0, 1), (1, 0, 0, 1), (1, 0, 1, 0), (0, 1, 1, 0), 
             (0, 1, 1, 1), (0, 1, 1, 2), (1, 0, 1, 1), (1, 0, 1, 2), 
             (1, 1, 0, 1), (1, 2, 0, 1), (1, 1, 1, 0), (1, 2, 1, 0), 
             (1, 1, 1, 1), (1, 1, 1, 2), (1, 2, 1, 1), (1, 2, 1, 2)]
        """
        seen = set()
        for p in self.coordinate_iter():
            if p in seen:
                continue
            seen.update(self.orbit(p))
            yield p

