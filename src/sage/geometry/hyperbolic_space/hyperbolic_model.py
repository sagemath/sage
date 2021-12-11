# -*- coding: utf-8 -*-
r"""
Hyperbolic Models

In this module, a hyperbolic model is a collection of data that allow
the user to implement new models of hyperbolic space with minimal effort.
The data include facts about the underlying set (such as whether the
model is bounded), facts about the metric (such as whether the model is
conformal), facts about the isometry group (such as whether it is a
linear or projective group), and more.  Generally speaking, any data
or method that pertains to the model itself -- rather than the points,
geodesics, or isometries of the model -- is implemented in this module.

Abstractly, a model of hyperbolic space is a connected, simply connected
manifold equipped with a complete Riemannian metric of constant curvature
`-1`.  This module records information sufficient to enable computations
in hyperbolic space without explicitly specifying the underlying set or
its Riemannian metric.  Although, see the
`SageManifolds <http://sagemanifolds.obspm.fr/>`_ project if
you would like to take this approach.

This module implements the abstract base class for a model of hyperbolic
space of arbitrary dimension.  It also contains the implementations of
specific models of hyperbolic geometry.

AUTHORS:

- Greg Laun (2013): Initial version.

EXAMPLES:

We illustrate how the classes in this module encode data by comparing
the upper half plane (UHP), Poincaré disk (PD) and hyperboloid (HM)
models.  First we create::

    sage: U = HyperbolicPlane().UHP()
    sage: P = HyperbolicPlane().PD()
    sage: H = HyperbolicPlane().HM()

We note that the UHP and PD models are bounded while the HM model is
not::

   sage: U.is_bounded() and P.is_bounded()
   True
   sage: H.is_bounded()
   False

The isometry groups of UHP and PD are projective, while that of HM is
linear::

    sage: U.is_isometry_group_projective()
    True
    sage: H.is_isometry_group_projective()
    False

The models are responsible for determining if the coordinates of points
and the matrix of linear maps are appropriate for constructing points
and isometries in hyperbolic space::

    sage: U.point_in_model(2 + I)
    True
    sage: U.point_in_model(2 - I)
    False
    sage: U.point_in_model(2)
    False
    sage: U.boundary_point_in_model(2)
    True
"""

#***********************************************************************
#
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.bindable_class import BindableClass
from sage.misc.lazy_import import lazy_import
from sage.functions.other import imag, real
from sage.misc.functional import sqrt
from sage.functions.all import arccosh
from sage.rings.cc import CC
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RR
from sage.rings.infinity import infinity
from sage.symbolic.constants import I
from sage.matrix.constructor import matrix
from sage.categories.homset import Hom

from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON, LORENTZ_GRAM
from sage.geometry.hyperbolic_space.hyperbolic_point import (
            HyperbolicPoint, HyperbolicPointUHP)
from sage.geometry.hyperbolic_space.hyperbolic_isometry import (
            HyperbolicIsometry, HyperbolicIsometryUHP,
            HyperbolicIsometryPD, HyperbolicIsometryKM, moebius_transform)
from sage.geometry.hyperbolic_space.hyperbolic_geodesic import (
            HyperbolicGeodesic, HyperbolicGeodesicUHP, HyperbolicGeodesicPD,
            HyperbolicGeodesicKM, HyperbolicGeodesicHM)
from sage.geometry.hyperbolic_space.hyperbolic_coercion import (
            CoercionUHPtoPD, CoercionUHPtoKM, CoercionUHPtoHM,
            CoercionPDtoUHP, CoercionPDtoKM, CoercionPDtoHM,
            CoercionKMtoUHP, CoercionKMtoPD, CoercionKMtoHM,
            CoercionHMtoUHP, CoercionHMtoPD, CoercionHMtoKM)

lazy_import('sage.modules.free_module_element', 'vector')

#####################################################################
## Abstract model


class HyperbolicModel(Parent, UniqueRepresentation, BindableClass):
    r"""
    Abstract base class for hyperbolic models.
    """
    Element = HyperbolicPoint
    _Geodesic = HyperbolicGeodesic
    _Isometry = HyperbolicIsometry

    def __init__(self, space, name, short_name, bounded, conformal,
                 dimension, isometry_group, isometry_group_is_projective):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: TestSuite(UHP).run()
            sage: PD = HyperbolicPlane().PD()
            sage: TestSuite(PD).run()
            sage: KM = HyperbolicPlane().KM()
            sage: TestSuite(KM).run()
            sage: HM = HyperbolicPlane().HM()
            sage: TestSuite(HM).run()
        """
        self._name = name
        self._short_name = short_name
        self._bounded = bounded
        self._conformal = conformal
        self._dimension = dimension
        self._isometry_group = isometry_group
        self._isometry_group_is_projective = isometry_group_is_projective
        from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
        Parent.__init__(self, category=HyperbolicModels(space))

    def _repr_(self):  # Abstract
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP()
            Hyperbolic plane in the Upper Half Plane Model
        """
        return u'Hyperbolic plane in the {}'.format(self._name)

    def _element_constructor_(self, x, is_boundary=None, **graphics_options): #Abstract
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP(2 + I)
            Point in UHP I + 2
        """
        return self.get_point(x, is_boundary, **graphics_options)

    def name(self):  # Abstract
        """
        Return the name of this model.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.name()
            'Upper Half Plane Model'
        """
        return self._name

    def short_name(self):
        """
        Return the short name of this model.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.short_name()
            'UHP'
        """
        return self._short_name

    def is_bounded(self):
        """
        Return ``True`` if ``self`` is a bounded model.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().is_bounded()
            True
            sage: HyperbolicPlane().PD().is_bounded()
            True
            sage: HyperbolicPlane().KM().is_bounded()
            True
            sage: HyperbolicPlane().HM().is_bounded()
            False
        """
        return self._bounded

    def is_conformal(self):
        """
        Return ``True`` if ``self`` is a conformal model.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.is_conformal()
            True
        """
        return self._conformal

    def is_isometry_group_projective(self):
        """
        Return ``True`` if the isometry group of ``self`` is projective.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.is_isometry_group_projective()
            True
        """
        return self._isometry_group_is_projective

    def point_in_model(self, p):
        r"""
        Return ``True`` if the point ``p`` is in the interior of the
        given model and ``False`` otherwise.

        INPUT:

        - any object that can converted into a complex number

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: HyperbolicPlane().UHP().point_in_model(I)
            True
            sage: HyperbolicPlane().UHP().point_in_model(-I)
            False
        """
        return True

    def point_test(self, p): #Abstract
        r"""
        Test whether a point is in the model.  If the point is in the
        model, do nothing.  Otherwise, raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicPlane().UHP().point_test(2 + I)
            sage: HyperbolicPlane().UHP().point_test(2 - I)
            Traceback (most recent call last):
            ...
            ValueError: -I + 2 is not a valid point in the UHP model
        """
        if not (self.point_in_model(p) or self.boundary_point_in_model(p)):
            error_string = "{0} is not a valid point in the {1} model"
            raise ValueError(error_string.format(p, self._short_name))

    def boundary_point_in_model(self, p): #Abstract
        r"""
        Return ``True`` if the point is on the ideal boundary of hyperbolic
        space and ``False`` otherwise.

        INPUT:

        - any object that can converted into a complex number

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: HyperbolicPlane().UHP().boundary_point_in_model(I)
            False
        """
        return True

    def bdry_point_test(self, p): #Abstract
        r"""
        Test whether a point is in the model.  If the point is in the
        model, do nothing; otherwise raise a ``ValueError``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().bdry_point_test(2)
            sage: HyperbolicPlane().UHP().bdry_point_test(1 + I)
            Traceback (most recent call last):
            ...
            ValueError: I + 1 is not a valid boundary point in the UHP model
        """
        if not self._bounded or not self.boundary_point_in_model(p):
            error_string = "{0} is not a valid boundary point in the {1} model"
            raise ValueError(error_string.format(p, self._short_name))

    def isometry_in_model(self, A): #Abstract
        r"""
        Return ``True`` if the input matrix represents an isometry of the
        given model and ``False`` otherwise.

        INPUT:

        - a matrix that represents an isometry in the appropriate model

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: HyperbolicPlane().UHP().isometry_in_model(identity_matrix(2))
            True

            sage: HyperbolicPlane().UHP().isometry_in_model(identity_matrix(3))
            False
        """
        return True

    def isometry_test(self, A): #Abstract
        r"""
        Test whether an isometry ``A`` is in the model.

        If the isometry is in the model, do nothing. Otherwise, raise
        a ``ValueError``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().isometry_test(identity_matrix(2))
            sage: HyperbolicPlane().UHP().isometry_test(matrix(2, [I,1,2,1]))
            Traceback (most recent call last):
            ...
            ValueError:
            [I 1]
            [2 1] is not a valid isometry in the UHP model
        """
        if not self.isometry_in_model(A):
            error_string = "\n{0} is not a valid isometry in the {1} model"
            raise ValueError(error_string.format(A, self._short_name))

    def get_point(self, coordinates, is_boundary=None, **graphics_options):
        r"""
        Return a point in ``self``.

        Automatically determine the type of point to return given either:

        #. the coordinates of a point in the interior or ideal boundary
           of hyperbolic space, or
        #. a :class:`~sage.geometry.hyperbolic_space.hyperbolic_point.HyperbolicPoint` object.

        INPUT:

        - a point in hyperbolic space or on the ideal boundary

        OUTPUT:

        - a :class:`~sage.geometry.hyperbolic_space.hyperbolic_point.HyperbolicPoint`

        EXAMPLES:

        We can create an interior point via the coordinates::

            sage: HyperbolicPlane().UHP().get_point(2*I)
            Point in UHP 2*I

        Or we can create a boundary point via the coordinates::

            sage: HyperbolicPlane().UHP().get_point(23)
            Boundary point in UHP 23

        However we cannot create points outside of our model::

            sage: HyperbolicPlane().UHP().get_point(12 - I)
            Traceback (most recent call last):
            ...
            ValueError: -I + 12 is not a valid point in the UHP model

        ::

            sage: HyperbolicPlane().UHP().get_point(2 + 3*I)
            Point in UHP 3*I + 2

            sage: HyperbolicPlane().PD().get_point(0)
            Point in PD 0

            sage: HyperbolicPlane().KM().get_point((0,0))
            Point in KM (0, 0)

            sage: HyperbolicPlane().HM().get_point((0,0,1))
            Point in HM (0, 0, 1)

            sage: p = HyperbolicPlane().UHP().get_point(I, color="red")
            sage: p.graphics_options()
            {'color': 'red'}

        ::

            sage: HyperbolicPlane().UHP().get_point(12)
            Boundary point in UHP 12

            sage: HyperbolicPlane().UHP().get_point(infinity)
            Boundary point in UHP +Infinity

            sage: HyperbolicPlane().PD().get_point(I)
            Boundary point in PD I

            sage: HyperbolicPlane().KM().get_point((0,-1))
            Boundary point in KM (0, -1)
        """

        if isinstance(coordinates, HyperbolicPoint):
            if coordinates.parent() is not self:
                coordinates = self(coordinates)
            coordinates.update_graphics(True, **graphics_options)
            return coordinates #both Point and BdryPoint

        if is_boundary is None:
            is_boundary = self.boundary_point_in_model(coordinates)
        return self.element_class(self, coordinates, is_boundary, **graphics_options)

    def get_geodesic(self, start, end=None, **graphics_options): #Abstract
        r"""
        Return a geodesic in the appropriate model.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2*I)
            Geodesic in UHP from I to 2*I

            sage: HyperbolicPlane().PD().get_geodesic(0, I/2)
            Geodesic in PD from 0 to 1/2*I

            sage: HyperbolicPlane().KM().get_geodesic((1/2, 1/2), (0,0))
            Geodesic in KM from (1/2, 1/2) to (0, 0)

            sage: HyperbolicPlane().HM().get_geodesic((0,0,1), (1,0, sqrt(2)))
            Geodesic in HM from (0, 0, 1) to (1, 0, sqrt(2))

        TESTS::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(UHP.get_point(I), UHP.get_point(2 + I))
            sage: h = UHP.get_geodesic(I, 2 + I)
            sage: g == h
            True
        """
        if end is None:
            if isinstance(start, HyperbolicGeodesic):
                G = start
                if G.model() is not self:
                    G = G.to_model(self)
                G.update_graphics(True, **graphics_options)
                return G
            raise ValueError("the start and end points must be specified")
        return self._Geodesic(self, self(start), self(end), **graphics_options)

    def get_isometry(self, A):
        r"""
        Return an isometry in ``self`` from the matrix ``A`` in the
        isometry group of ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_isometry(identity_matrix(2))
            Isometry in UHP
            [1 0]
            [0 1]

            sage: HyperbolicPlane().PD().get_isometry(identity_matrix(2))
            Isometry in PD
            [1 0]
            [0 1]

            sage: HyperbolicPlane().KM().get_isometry(identity_matrix(3))
            Isometry in KM
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: HyperbolicPlane().HM().get_isometry(identity_matrix(3))
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if isinstance(A, HyperbolicIsometry):
            if A.model() is not self:
                return A.to_model(self)
            return A
        return self._Isometry(self, A)

    def random_element(self, **kwargs):
        r"""
        Return a random point in ``self``.

        The points are uniformly distributed over the rectangle
        `[-10, 10] \times [0, 10 i]` in the upper half plane model.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().random_element()
            sage: bool((p.coordinates().imag()) > 0)
            True

            sage: p = HyperbolicPlane().PD().random_element()
            sage: HyperbolicPlane().PD().point_in_model(p.coordinates())
            True

            sage: p = HyperbolicPlane().KM().random_element()
            sage: HyperbolicPlane().KM().point_in_model(p.coordinates())
            True

            sage: p = HyperbolicPlane().HM().random_element().coordinates()
            sage: bool((p[0]**2 + p[1]**2 - p[2]**2 - 1) < 10**-8)
            True
        """
        return self.random_point(**kwargs)

    def random_point(self, **kwargs):
        r"""
        Return a random point of ``self``.

        The points are uniformly distributed over the rectangle
        `[-10, 10] \times [0, 10 i]` in the upper half plane model.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().random_point()
            sage: bool((p.coordinates().imag()) > 0)
            True

            sage: PD = HyperbolicPlane().PD()
            sage: p = PD.random_point()
            sage: PD.point_in_model(p.coordinates())
            True
        """
        R = self.realization_of().a_realization()
        return self(R.random_point(**kwargs))

    def random_geodesic(self, **kwargs):
        r"""
        Return a random hyperbolic geodesic.

        Return the geodesic between two random points.

        EXAMPLES::

            sage: h = HyperbolicPlane().PD().random_geodesic()
            sage: bool((h.endpoints()[0].coordinates()).imag() >= 0)
            True
        """
        R = self.realization_of().a_realization()
        g_ends = [R.random_point(**kwargs) for k in range(2)]
        return self.get_geodesic(self(g_ends[0]), self(g_ends[1]))

    def random_isometry(self, preserve_orientation=True, **kwargs):
        r"""
        Return a random isometry in the model of ``self``.

        INPUT:

        - ``preserve_orientation`` -- if ``True`` return an
          orientation-preserving isometry

        OUTPUT:

        - a hyperbolic isometry

        EXAMPLES::

            sage: A = HyperbolicPlane().PD().random_isometry()
            sage: A.preserves_orientation()
            True
            sage: B = HyperbolicPlane().PD().random_isometry(preserve_orientation=False)
            sage: B.preserves_orientation()
            False
        """
        R = self.realization_of().a_realization()
        A = R.random_isometry(preserve_orientation, **kwargs)
        return A.to_model(self)

    ################
    # Dist methods #
    ################

    def dist(self, a, b):
        r"""
        Calculate the hyperbolic distance between ``a`` and ``b``.

        INPUT:

        - ``a``, ``b`` -- a point or geodesic

        OUTPUT:

        - the hyperbolic distance

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: p1 = UHP.get_point(5 + 7*I)
            sage: p2 = UHP.get_point(1.0 + I)
            sage: UHP.dist(p1, p2)
            2.23230104635820

            sage: PD = HyperbolicPlane().PD()
            sage: p1 = PD.get_point(0)
            sage: p2 = PD.get_point(I/2)
            sage: PD.dist(p1, p2)
            arccosh(5/3)

            sage: UHP(p1).dist(UHP(p2))
            arccosh(5/3)

            sage: KM = HyperbolicPlane().KM()
            sage: p1 = KM.get_point((0, 0))
            sage: p2 = KM.get_point((1/2, 1/2))
            sage: numerical_approx(KM.dist(p1, p2))
            0.881373587019543

            sage: HM = HyperbolicPlane().HM()
            sage: p1 = HM.get_point((0,0,1))
            sage: p2 = HM.get_point((1,0,sqrt(2)))
            sage: numerical_approx(HM.dist(p1, p2))
            0.881373587019543

        Distance between a point and itself is 0::

            sage: p = UHP.get_point(47 + I)
            sage: UHP.dist(p, p)
            0

        Points on the boundary are infinitely far from interior points::

            sage: UHP.get_point(3).dist(UHP.get_point(I))
            +Infinity

        TESTS::

            sage: UHP.dist(UHP.get_point(I), UHP.get_point(2*I))
            arccosh(5/4)
            sage: UHP.dist(I, 2*I)
            arccosh(5/4)
        """
        def coords(x):
            return self(x).coordinates()

        if isinstance(a, HyperbolicGeodesic):
            if isinstance(b, HyperbolicGeodesic):
                if not a.is_parallel(b):
                    return 0

                if a.is_ultra_parallel(b):
                    perp = a.common_perpendicular(b)
                    # Find where a and b intersect the common perp...
                    p = a.intersection(perp)[0]
                    q = b.intersection(perp)[0]
                    # ...and return their distance
                    return self._dist_points(coords(p), coords(q))

                raise NotImplementedError("can only compute distance between"
                                          " ultra-parallel and intersecting geodesics")

            # If only one is a geodesic, make sure it's b to make things easier
            a,b = b,a

        if isinstance(b, HyperbolicGeodesic):
            (p, q) = b.ideal_endpoints()
            return self._dist_geod_point(coords(p), coords(q), coords(a))

        return self._dist_points(coords(a), coords(b))

    def _dist_points(self, p1, p2):
        r"""
        Compute the distance between two points.

        INPUT:

        - ``p1``, ``p2`` -- the coordinates of the points

        EXAMPLES::

           sage: HyperbolicPlane().PD()._dist_points(3/5*I, 0)
           arccosh(17/8)
        """
        R = self.realization_of().a_realization()
        phi = R.coerce_map_from(self)
        return R._dist_points(phi.image_coordinates(p1), phi.image_coordinates(p2))

    def _dist_geod_point(self, start, end, p):
        r"""
        Return the hyperbolic distance from a given hyperbolic geodesic
        and a hyperbolic point.

        INPUT:

        - ``start`` -- the start ideal point coordinates of the geodesic
        - ``end`` -- the end ideal point coordinates of the geodesic
        - ``p`` -- the coordinates of the point

        OUTPUT:

        - the hyperbolic distance

        EXAMPLES::

            sage: HyperbolicPlane().PD()._dist_geod_point(3/5*I + 4/5, I, 0)
            arccosh(1/10*sqrt(5)*((sqrt(5) - 1)^2 + 4) + 1)

        If `p` is a boundary point, the distance is infinity::

            sage: HyperbolicPlane().PD()._dist_geod_point(3/5*I + 4/5, I, 12/13*I + 5/13)
            +Infinity
        """
        R = self.realization_of().a_realization()
        assert R is not self

        def phi(c):
            return R.coerce_map_from(self).image_coordinates(c)
        return R._dist_geod_point(phi(start), phi(end), phi(p))

    ####################
    # Isometry methods #
    ####################

    def isometry_from_fixed_points(self, repel, attract):
        r"""
        Given two fixed points ``repel`` and ``attract`` as hyperbolic
        points return a hyperbolic isometry with ``repel`` as repelling
        fixed point and ``attract`` as attracting fixed point.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: PD.isometry_from_fixed_points(-i, i)
            Isometry in PD
            [   3/4  1/4*I]
            [-1/4*I    3/4]

        ::

            sage: p, q = PD.get_point(1/2 + I/2), PD.get_point(6/13 + 9/13*I)
            sage: PD.isometry_from_fixed_points(p, q)
            Traceback (most recent call last):
            ...
            ValueError: fixed points of hyperbolic elements must be ideal

            sage: p, q = PD.get_point(4/5 + 3/5*I), PD.get_point(-I)
            sage: PD.isometry_from_fixed_points(p, q)
            Isometry in PD
            [ 1/6*I - 2/3 -1/3*I - 1/6]
            [ 1/3*I - 1/6 -1/6*I - 2/3]
        """
        R = self.realization_of().a_realization()
        return R.isometry_from_fixed_points(R(self(repel)), R(self(attract))).to_model(self)


#####################################################################
## Upper half plane model

class HyperbolicModelUHP(HyperbolicModel):
    r"""
    Upper Half Plane model.
    """
    Element = HyperbolicPointUHP
    _Geodesic = HyperbolicGeodesicUHP
    _Isometry = HyperbolicIsometryUHP

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: TestSuite(UHP).run()
        """
        HyperbolicModel.__init__(self, space,
            name="Upper Half Plane Model", short_name="UHP",
            bounded=True, conformal=True, dimension=2,
            isometry_group="PSL(2, \\RR)", isometry_group_is_projective=True)

    def _coerce_map_from_(self, X):
        """
        Return if there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.has_coerce_map_from(HyperbolicPlane().PD())
            True
            sage: UHP.has_coerce_map_from(HyperbolicPlane().KM())
            True
            sage: UHP.has_coerce_map_from(HyperbolicPlane().HM())
            True
            sage: UHP.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelPD):
            return CoercionPDtoUHP(Hom(X, self))
        if isinstance(X, HyperbolicModelKM):
            return CoercionKMtoUHP(Hom(X, self))
        if isinstance(X, HyperbolicModelHM):
            return CoercionHMtoUHP(Hom(X, self))
        return super(HyperbolicModelUHP, self)._coerce_map_from_(X)

    def point_in_model(self, p):
        r"""
        Check whether a complex number lies in the open upper half plane.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.point_in_model(1 + I)
            True
            sage: UHP.point_in_model(infinity)
            False
            sage: UHP.point_in_model(CC(infinity))
            False
            sage: UHP.point_in_model(RR(infinity))
            False
            sage: UHP.point_in_model(1)
            False
            sage: UHP.point_in_model(12)
            False
            sage: UHP.point_in_model(1 - I)
            False
            sage: UHP.point_in_model(-2*I)
            False
            sage: UHP.point_in_model(I)
            True
            sage: UHP.point_in_model(0) # Not interior point
            False
        """
        if isinstance(p, HyperbolicPoint):
            return p.is_boundary()
        return bool(imag(CC(p)) > 0)

    def boundary_point_in_model(self, p):
        r"""
        Check whether a complex number is a real number or ``\infty``.
        In the ``UHP.model_name_name``, this is the ideal boundary of
        hyperbolic space.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.boundary_point_in_model(1 + I)
            False
            sage: UHP.boundary_point_in_model(infinity)
            True
            sage: UHP.boundary_point_in_model(CC(infinity))
            True
            sage: UHP.boundary_point_in_model(RR(infinity))
            True
            sage: UHP.boundary_point_in_model(1)
            True
            sage: UHP.boundary_point_in_model(12)
            True
            sage: UHP.boundary_point_in_model(1 - I)
            False
            sage: UHP.boundary_point_in_model(-2*I)
            False
            sage: UHP.boundary_point_in_model(0)
            True
            sage: UHP.boundary_point_in_model(I)
            False
        """
        if isinstance(p, HyperbolicPoint):
            return p.is_boundary()
        im = abs(imag(CC(p)).n())
        return (im < EPSILON) or bool(p == infinity)

    def isometry_in_model(self, A):
        r"""
        Check that ``A`` acts as an isometry on the upper half plane.
        That is, ``A`` must be an invertible `2 \times 2` matrix with real
        entries.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: A = matrix(2,[1,2,3,4])
            sage: UHP.isometry_in_model(A)
            True
            sage: B = matrix(2,[I,2,4,1])
            sage: UHP.isometry_in_model(B)
            False

        An example of a matrix `A` such that `\det(A) \neq 1`, but the `A`
        acts isometrically::

            sage: C = matrix(2,[10,0,0,10])
            sage: UHP.isometry_in_model(C)
            True
        """
        if isinstance(A, HyperbolicIsometry):
            return True
        return bool(A.ncols() == 2 and A.nrows() == 2 and
                    sum([k in RR for k in A.list()]) == 4 and
                    abs(A.det()) > -EPSILON)

    def get_background_graphic(self, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.
        For the upper half space, the background object is the ideal boundary.

        EXAMPLES::

            sage: hp = HyperbolicPlane().UHP().get_background_graphic()
        """
        from sage.plot.line import line
        bd_min = bdry_options.get('bd_min', -5)
        bd_max = bdry_options.get('bd_max', 5)
        return line(((bd_min, 0), (bd_max, 0)), color='black')

    ################
    # Dist methods #
    ################

    def _dist_points(self, p1, p2):
        r"""
        Compute the distance between two points in the Upper Half Plane
        using the hyperbolic metric.

        INPUT:

        - ``p1``, ``p2`` -- the coordinates of the points

        EXAMPLES::

           sage: HyperbolicPlane().UHP()._dist_points(4.0*I, I)
           1.38629436111989
        """
        num = (real(p2) - real(p1))**2 + (imag(p2) - imag(p1))**2
        denom = 2 * imag(p1) * imag(p2)
        if denom == 0:
            return infinity
        return arccosh(1 + num/denom)

    def _dist_geod_point(self, start, end, p):
        r"""
        Return the hyperbolic distance from a given hyperbolic geodesic
        and a hyperbolic point.

        INPUT:

        - ``start`` -- the start ideal point coordinates of the geodesic
        - ``end`` -- the end ideal point coordinates of the geodesic
        - ``p`` -- the coordinates of the point

        OUTPUT:

        - the hyperbolic distance

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP._dist_geod_point(2, infinity, I)
            arccosh(1/10*sqrt(5)*((sqrt(5) - 1)^2 + 4) + 1)

        If `p` is a boundary point, the distance is infinity::

            sage: HyperbolicPlane().UHP()._dist_geod_point(2, infinity, 5)
            +Infinity
        """
        # Here is the trick for computing distance to a geodesic:
        # find an isometry mapping the geodesic to the geodesic between
        # 0 and infinity (so corresponding to the line imag(z) = 0.
        # then any complex number is r exp(i*theta) in polar coordinates.
        # the mutual perpendicular between this point and imag(z) = 0
        # intersects imag(z) = 0 at ri.  So we calculate the distance
        # between r exp(i*theta) and ri after we transform the original
        # point.
        if start + end != infinity:
            # Not a straight line:
            # Map the endpoints to 0 and infinity and the midpoint to 1.
            T = HyperbolicGeodesicUHP._crossratio_matrix(start, (start + end)/2, end)
        else:
            # Is a straight line:
            # Map the endpoints to 0 and infinity and another endpoint to 1.
            T = HyperbolicGeodesicUHP._crossratio_matrix(start, start + 1, end)
        x = moebius_transform(T, p)
        return self._dist_points(x, abs(x)*I)

    #################
    # Point Methods #
    #################

    def random_point(self, **kwargs):
        r"""
        Return a random point in the upper half plane. The points are
        uniformly distributed over the rectangle `[-10, 10] \times [0, 10i]`.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().random_point().coordinates()
            sage: bool((p.imag()) > 0)
            True
        """
        # TODO: use **kwargs to allow these to be set
        real_min = -10
        real_max = 10
        imag_min = 0
        imag_max = 10
        p = RR.random_element(min=real_min, max=real_max) \
            + I * RR.random_element(min=imag_min, max=imag_max)
        return self.get_point(p)

    ####################
    # Isometry Methods #
    ####################

    def isometry_from_fixed_points(self, repel, attract):
        r"""
        Given two fixed points ``repel`` and ``attract`` as complex
        numbers return a hyperbolic isometry with ``repel`` as repelling
        fixed point and ``attract`` as attracting fixed point.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.isometry_from_fixed_points(2 + I, 3 + I)
            Traceback (most recent call last):
            ...
            ValueError: fixed points of hyperbolic elements must be ideal

            sage: UHP.isometry_from_fixed_points(2, 0)
            Isometry in UHP
            [  -1    0]
            [-1/3 -1/3]

        TESTS::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.isometry_from_fixed_points(0, 4)
            Isometry in UHP
            [  -1    0]
            [-1/5 -1/5]
            sage: UHP.isometry_from_fixed_points(UHP.get_point(0), UHP.get_point(4))
            Isometry in UHP
            [  -1    0]
            [-1/5 -1/5]
        """
        if isinstance(repel, HyperbolicPoint):
            repel = repel._coordinates
        if isinstance(attract, HyperbolicPoint):
            attract = attract._coordinates

        if imag(repel) + imag(attract) > EPSILON:
            raise ValueError("fixed points of hyperbolic elements must be ideal")
        repel = real(repel)
        attract = real(attract)
        if repel == infinity:
            A = self._moebius_sending([infinity, attract, attract + 1],
                                     [infinity, attract, attract + 2])
        elif attract == infinity:
            A = self._moebius_sending([repel, infinity, repel + 1],
                                     [repel, infinity, repel + 2])
        else:
            A = self._moebius_sending([repel, attract, infinity],
                                     [repel, attract, max(repel, attract) + 1])
        return self.get_isometry(A)

    def random_isometry(self, preserve_orientation=True, **kwargs):
        r"""
        Return a random isometry in the Upper Half Plane model.

        INPUT:

        - ``preserve_orientation`` -- if ``True`` return an
          orientation-preserving isometry

        OUTPUT:

        - a hyperbolic isometry

        EXAMPLES::

            sage: A = HyperbolicPlane().UHP().random_isometry()
            sage: B = HyperbolicPlane().UHP().random_isometry(preserve_orientation=False)
            sage: B.preserves_orientation()
            False
        """
        [a,b,c,d] = [RR.random_element() for k in range(4)]
        while abs(a*d - b*c) < EPSILON:
            [a,b,c,d] = [RR.random_element() for k in range(4)]
        M = matrix(RDF, 2,[a,b,c,d])
        M = M / (M.det()).abs().sqrt()
        if M.det() > 0:
            if not preserve_orientation:
                M = M * matrix(2,[0,1,1,0])
        elif preserve_orientation:
            M = M * matrix(2,[0,1,1,0])
        return self._Isometry(self, M, check=False)

    ###################
    # Helping Methods #
    ###################

    @staticmethod
    def _moebius_sending(z, w): #UHP
        r"""
        Given two lists ``z`` and ``w`` of three points each in
        `\mathbb{CP}^1`, return the linear fractional transformation
        taking the points in ``z`` to the points in ``w``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import moebius_transform
            sage: bool(abs(moebius_transform(HyperbolicModelUHP._moebius_sending([1,2,infinity],[3 - I, 5*I,-12]),1) - 3 + I) < 10^-4)
            True
            sage: bool(abs(moebius_transform(HyperbolicModelUHP._moebius_sending([1,2,infinity],[3 - I, 5*I,-12]),2) - 5*I) < 10^-4)
            True
            sage: bool(abs(moebius_transform(HyperbolicModelUHP._moebius_sending([1,2,infinity],[3 - I, 5*I,-12]),infinity) + 12) < 10^-4)
            True
        """
        if len(z) != 3 or len(w) != 3:
            raise TypeError("moebius_sending requires each list to be three points long")
        A = HyperbolicGeodesicUHP._crossratio_matrix(z[0],z[1],z[2])
        B = HyperbolicGeodesicUHP._crossratio_matrix(w[0],w[1],w[2])
        return B.inverse() * A

#####################################################################
## Poincaré disk model


class HyperbolicModelPD(HyperbolicModel):
    r"""
    Poincaré Disk Model.
    """
    _Geodesic = HyperbolicGeodesicPD
    _Isometry = HyperbolicIsometryPD

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: TestSuite(PD).run()
        """
        # name should really be 'Poincaré Disk Model', but utf8 is not
        # accepted by repr
        HyperbolicModel.__init__(self, space,
                                 name=u'Poincare Disk Model', short_name="PD",
                                 bounded=True, conformal=True, dimension=2,
                                 isometry_group="PU(1, 1)",
                                 isometry_group_is_projective=True)

    def _coerce_map_from_(self, X):
        """
        Return if there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: PD.has_coerce_map_from(HyperbolicPlane().UHP())
            True
            sage: PD.has_coerce_map_from(HyperbolicPlane().KM())
            True
            sage: PD.has_coerce_map_from(HyperbolicPlane().HM())
            True
            sage: PD.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelUHP):
            return CoercionUHPtoPD(Hom(X, self))
        if isinstance(X, HyperbolicModelKM):
            return CoercionKMtoPD(Hom(X, self))
        if isinstance(X, HyperbolicModelHM):
            return CoercionHMtoPD(Hom(X, self))
        return super(HyperbolicModelPD, self)._coerce_map_from_(X)

    def point_in_model(self, p):
        r"""
        Check whether a complex number lies in the open unit disk.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: PD.point_in_model(1.00)
            False
            sage: PD.point_in_model(1/2 + I/2)
            True
            sage: PD.point_in_model(1 + .2*I)
            False
        """
        if isinstance(p, HyperbolicPoint):
            return p.is_boundary()
        return bool(abs(CC(p)) < 1)

    def boundary_point_in_model(self, p):
        r"""
        Check whether a complex number lies in the open unit disk.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: PD.boundary_point_in_model(1.00)
            True
            sage: PD.boundary_point_in_model(1/2 + I/2)
            False
            sage: PD.boundary_point_in_model(1 + .2*I)
            False
        """
        if isinstance(p, HyperbolicPoint):
            return p.is_boundary()
        return bool(abs(abs(CC(p)) - 1) < EPSILON)

    def isometry_in_model(self, A):
        r"""
        Check if the given matrix ``A`` is in the group `U(1,1)`.

        EXAMPLES::

            sage: z = [CC.random_element() for k in range(2)]; z.sort(key=abs)
            sage: A = matrix(2,[z[1], z[0],z[0].conjugate(),z[1].conjugate()])
            sage: HyperbolicPlane().PD().isometry_in_model(A)
            True
        """
        if isinstance(A, HyperbolicIsometry):
            return True
        # alpha = A[0][0]
        # beta = A[0][1]
        # Orientation preserving and reversing
        return (HyperbolicIsometryPD._orientation_preserving(A) or
                HyperbolicIsometryPD._orientation_preserving(I * A))

    def get_background_graphic(self, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.

        For the Poincaré disk, the background object is the ideal boundary.

        EXAMPLES::

            sage: circ = HyperbolicPlane().PD().get_background_graphic()
        """
        from sage.plot.circle import circle
        return circle((0, 0), 1, axes=False, color='black')


#####################################################################
## Klein disk model

class HyperbolicModelKM(HyperbolicModel):
    r"""
    Klein Model.
    """
    _Geodesic = HyperbolicGeodesicKM
    _Isometry = HyperbolicIsometryKM

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: TestSuite(KM).run()
        """
        HyperbolicModel.__init__(self, space,
            name="Klein Disk Model", short_name="KM",
            bounded=True, conformal=False, dimension=2,
            isometry_group="PSO(2, 1)", isometry_group_is_projective=True)

    def _coerce_map_from_(self, X):
        """
        Return if there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().UHP()
            sage: KM.has_coerce_map_from(HyperbolicPlane().UHP())
            True
            sage: KM.has_coerce_map_from(HyperbolicPlane().PD())
            True
            sage: KM.has_coerce_map_from(HyperbolicPlane().HM())
            True
            sage: KM.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelUHP):
            return CoercionUHPtoKM(Hom(X, self))
        if isinstance(X, HyperbolicModelPD):
            return CoercionPDtoKM(Hom(X, self))
        if isinstance(X, HyperbolicModelHM):
            return CoercionHMtoKM(Hom(X, self))
        return super(HyperbolicModelKM, self)._coerce_map_from_(X)

    def point_in_model(self, p):
        r"""
        Check whether a point lies in the open unit disk.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: KM.point_in_model((1, 0))
            False
            sage: KM.point_in_model((1/2, 1/2))
            True
            sage: KM.point_in_model((1, .2))
            False
        """
        if isinstance(p, HyperbolicPoint):
            return p.is_boundary()
        return len(p) == 2 and bool(p[0]**2 + p[1]**2 < 1)

    def boundary_point_in_model(self, p):
        r"""
        Check whether a point lies in the unit circle, which corresponds
        to the ideal boundary of the hyperbolic plane in the Klein model.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: KM.boundary_point_in_model((1, 0))
            True
            sage: KM.boundary_point_in_model((1/2, 1/2))
            False
            sage: KM.boundary_point_in_model((1, .2))
            False
        """
        if isinstance(p, HyperbolicPoint):
            return p.is_boundary()
        return len(p) == 2 and bool(abs(p[0]**2 + p[1]**2 - 1) < EPSILON)

    def isometry_in_model(self, A):
        r"""
        Check if the given matrix ``A`` is in the group `SO(2,1)`.

        EXAMPLES::

            sage: A = matrix(3, [[1, 0, 0], [0, 17/8, 15/8], [0, 15/8, 17/8]])
            sage: HyperbolicPlane().KM().isometry_in_model(A)
            True
        """
        if isinstance(A, HyperbolicIsometry):
            return True
        return bool((A*LORENTZ_GRAM*A.transpose() - LORENTZ_GRAM).norm()**2 <
                    EPSILON)

    def get_background_graphic(self, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.

        For the Klein model, the background object is the ideal boundary.

        EXAMPLES::

            sage: circ = HyperbolicPlane().KM().get_background_graphic()
        """
        from sage.plot.circle import circle
        return circle((0, 0), 1, axes=False, color='black')

#####################################################################
## Hyperboloid model


class HyperbolicModelHM(HyperbolicModel):
    r"""
    Hyperboloid Model.
    """
    _Geodesic = HyperbolicGeodesicHM

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: TestSuite(HM).run()
        """
        HyperbolicModel.__init__(self, space,
            name="Hyperboloid Model", short_name="HM",
            bounded=False, conformal=True, dimension=2,
            isometry_group="SO(2, 1)", isometry_group_is_projective=False)

    def _coerce_map_from_(self, X):
        """
        Return if there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().UHP()
            sage: HM.has_coerce_map_from(HyperbolicPlane().UHP())
            True
            sage: HM.has_coerce_map_from(HyperbolicPlane().PD())
            True
            sage: HM.has_coerce_map_from(HyperbolicPlane().KM())
            True
            sage: HM.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelUHP):
            return CoercionUHPtoHM(Hom(X, self))
        if isinstance(X, HyperbolicModelPD):
            return CoercionPDtoHM(Hom(X, self))
        if isinstance(X, HyperbolicModelKM):
            return CoercionKMtoHM(Hom(X, self))
        return super(HyperbolicModelHM, self)._coerce_map_from_(X)

    def point_in_model(self, p):
        r"""
        Check whether a complex number lies in the hyperboloid.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: HM.point_in_model((0,0,1))
            True
            sage: HM.point_in_model((1,0,sqrt(2)))
            True
            sage: HM.point_in_model((1,2,1))
            False
        """
        if isinstance(p, HyperbolicPoint):
            return p.is_boundary()
        return len(p) == 3 and bool(abs(p[0]**2 + p[1]**2 - p[2]**2 + 1) < EPSILON)

    def boundary_point_in_model(self, p):
        r"""
        Return ``False`` since the Hyperboloid model has no boundary points.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: HM.boundary_point_in_model((0,0,1))
            False
            sage: HM.boundary_point_in_model((1,0,sqrt(2)))
            False
            sage: HM.boundary_point_in_model((1,2,1))
            False
        """
        return False

    def isometry_in_model(self, A):
        r"""
        Test that the matrix ``A`` is in the group `SO(2,1)^+`.

        EXAMPLES::

           sage: A = diagonal_matrix([1,1,-1])
           sage: HyperbolicPlane().HM().isometry_in_model(A)
           True
        """
        if isinstance(A, HyperbolicIsometry):
            return True
        return bool((A*LORENTZ_GRAM*A.transpose() - LORENTZ_GRAM).norm()**2 < EPSILON)

    def get_background_graphic(self, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.
        For the hyperboloid model, the background object is the hyperboloid
        itself.

        EXAMPLES::

            sage: H = HyperbolicPlane().HM().get_background_graphic()
        """
        from sage.plot.plot3d.all import plot3d
        from sage.symbolic.ring import SR
        hyperboloid_opacity = bdry_options.get('hyperboloid_opacity', .1)
        z_height = bdry_options.get('z_height', 7.0)
        x_max = sqrt((z_height ** 2 - 1) / 2.0)
        x = SR.var('x')
        y = SR.var('y')
        return plot3d((1 + x ** 2 + y ** 2).sqrt(),
                      (x, -x_max, x_max), (y,-x_max, x_max),
                      opacity=hyperboloid_opacity, **bdry_options)
