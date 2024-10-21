# sage_setup: distribution = sagemath-symbolics
r"""
Hyperbolic Points

This module implements points in hyperbolic space of arbitrary dimension.
It also contains the implementations for specific models of
hyperbolic geometry.

This module also implements ideal points in hyperbolic space of arbitrary
dimension. It also contains the implementations for specific models
of hyperbolic geometry.

Note that not all models of hyperbolic space are bounded, meaning that
the ideal boundary is not the topological boundary of the set underlying
tho model.  For example, the unit disk model is bounded with boundary
given by the unit sphere.  The hyperboloid model is not bounded.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct points in the upper half plane model, abbreviated
UHP for convenience::

    sage: UHP = HyperbolicPlane().UHP()
    sage: UHP.get_point(2 + I)
    Point in UHP I + 2
    sage: g = UHP.get_point(3 + I)
    sage: g.dist(UHP.get_point(I))
    arccosh(11/2)

We can also construct boundary points in the upper half plane model::

    sage: UHP.get_point(3)
    Boundary point in UHP 3

Some more examples::

    sage: HyperbolicPlane().UHP().get_point(0)
    Boundary point in UHP 0

    sage: HyperbolicPlane().PD().get_point(I/2)
    Point in PD 1/2*I

    sage: HyperbolicPlane().KM().get_point((0,1))
    Boundary point in KM (0, 1)

    sage: HyperbolicPlane().HM().get_point((0,0,1))
    Point in HM (0, 0, 1)
"""

#***********************************************************************
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from collections.abc import Iterable
from sage.structure.element import Element
from sage.structure.richcmp import richcmp, op_NE
from sage.symbolic.constants import I
from sage.misc.latex import latex
from sage.structure.element import Matrix
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.infinity import infinity
from sage.rings.cc import CC
from sage.rings.real_mpfr import RR
from sage.functions.other import real, imag

from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometry

class HyperbolicPoint(Element):
    r"""
    Abstract base class for hyperbolic points.  This class should never
    be instantiated.

    INPUT:

    - ``model`` -- the model of the hyperbolic space
    - ``coordinates`` -- the coordinates of a hyperbolic point in the
      appropriate model
    - ``is_boundary`` -- whether the point is a boundary point
    - ``check`` -- boolean (default: ``True``); if ``True``, then check to make sure
      the coordinates give a valid point in the model

    EXAMPLES:

    Comparison between different models is performed via coercion::

        sage: UHP = HyperbolicPlane().UHP()
        sage: p = UHP.get_point(.2 + .3*I); p
        Point in UHP 0.200000000000000 + 0.300000000000000*I

        sage: PD = HyperbolicPlane().PD()
        sage: q = PD.get_point(0.2 + 0.3*I); q
        Point in PD 0.200000000000000 + 0.300000000000000*I

        sage: p == q
        False
        sage: PD(p)
        Point in PD 0.231213872832370 - 0.502890173410405*I

        sage: bool(p.coordinates() == q.coordinates())
        True

    Similarly for boundary points::

        sage: p = UHP.get_point(-1); p
        Boundary point in UHP -1

        sage: q = PD.get_point(-1); q
        Boundary point in PD -1

        sage: p == q
        True
        sage: PD(p)
        Boundary point in PD -1

    It is an error to specify a point that does not lie in the
    appropriate model::

        sage: HyperbolicPlane().UHP().get_point(0.2 - 0.3*I)
        Traceback (most recent call last):
        ...
        ValueError: 0.200000000000000 - 0.300000000000000*I is not a valid point in the UHP model

        sage: HyperbolicPlane().PD().get_point(1.2)
        Traceback (most recent call last):
        ...
        ValueError: 1.20000000000000 is not a valid point in the PD model

        sage: HyperbolicPlane().KM().get_point((1,1))
        Traceback (most recent call last):
        ...
        ValueError: (1, 1) is not a valid point in the KM model

        sage: HyperbolicPlane().HM().get_point((1, 1, 1))
        Traceback (most recent call last):
        ...
        ValueError: (1, 1, 1) is not a valid point in the HM model

    It is an error to specify an interior point of hyperbolic space as a
    boundary point::

        sage: HyperbolicPlane().UHP().get_point(0.2 + 0.3*I, is_boundary=True)
        Traceback (most recent call last):
        ...
        ValueError: 0.200000000000000 + 0.300000000000000*I is not a valid boundary point in the UHP model

    TESTS:

    In the PD model, the coordinates of a point are in the unit disk
    in the complex plane `\CC`::

        sage: HyperbolicPlane().PD().get_point(0)
        Point in PD 0
        sage: HyperbolicPlane().PD().get_point(1)
        Boundary point in PD 1

    In the KM model, the coordinates of a point are in the unit disk
    in the real plane `\RR^2`::

        sage: HyperbolicPlane().KM().get_point((0,0))
        Point in KM (0, 0)
        sage: HyperbolicPlane().KM().get_point((1,0))
        Boundary point in KM (1, 0)

    In the HM model, the coordinates of a point are on the
    hyperboloid given by `x^2 + y^2 - z^2 = -1`::

        sage: HyperbolicPlane().HM().get_point((0,0,1))
        Point in HM (0, 0, 1)
        sage: HyperbolicPlane().HM().get_point((0,0,2))
        Traceback (most recent call last):
        ...
        ValueError: (0, 0, 2) is not a valid point in the HM model
        sage: HyperbolicPlane().HM().get_point((1,0,0), is_boundary=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: boundary points are not implemented in the HM model
    """
    def __init__(self, model, coordinates, is_boundary, check=True, **graphics_options):
        r"""
        See ``HyperbolicPoint`` for full documentation.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().get_point(I)
            sage: TestSuite(p).run()
            sage: p1 = HyperbolicPlane().KM().get_point((0,0))
            sage: p2 = HyperbolicPlane().KM().get_point([0,0])
            sage: p1 == p2
            True
        """
        if is_boundary:
            if not model.is_bounded():
                raise NotImplementedError("boundary points are not implemented in the {0} model".format(model.short_name()))
            if check and not model.boundary_point_in_model(coordinates):
                raise ValueError(
                    "{0} is not a valid".format(coordinates) +
                    " boundary point in the {0} model".format(model.short_name()))
        elif check and not model.point_in_model(coordinates):
            raise ValueError(
                "{0} is not a valid".format(coordinates) +
                " point in the {0} model".format(model.short_name()))

        if isinstance(coordinates, Iterable):
            coordinates = vector(coordinates)
        self._coordinates = coordinates
        self._bdry = is_boundary
        self._graphics_options = graphics_options

        Element.__init__(self, model)

    #####################
    # "Private" Methods #
    #####################

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(3 + 4*I)
            Point in UHP 4*I + 3

            sage: HyperbolicPlane().PD().get_point(1/2 + I/2)
            Point in PD 1/2*I + 1/2

            sage: HyperbolicPlane().KM().get_point((1/2, 1/2))
            Point in KM (1/2, 1/2)

            sage: HyperbolicPlane().HM().get_point((0,0,1))
            Point in HM (0, 0, 1)

            sage: HyperbolicPlane().UHP().get_point(infinity)
            Boundary point in UHP +Infinity

            sage: HyperbolicPlane().PD().get_point(-1)
            Boundary point in PD -1

            sage: HyperbolicPlane().KM().get_point((0, -1))
            Boundary point in KM (0, -1)
        """
        if self._bdry:
            base = "Boundary point"
        else:
            base = "Point"
        return base + " in {0} {1}".format(self.parent().short_name(), self._coordinates)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = HyperbolicPlane().UHP().get_point(0)
            sage: latex(p)
            0
            sage: q = HyperbolicPlane().HM().get_point((0,0,1))
            sage: latex(q)
            \left(0,\,0,\,1\right)
        """
        return latex(self._coordinates)

    def _richcmp_(self, other, op):
        r"""
        Comparison of ``self`` and ``other``.

        EXAMPLES::

            sage: p1 = HyperbolicPlane().UHP().get_point(1 + I)
            sage: p2 = HyperbolicPlane().UHP().get_point(2 + I)
            sage: p1 == p2
            False
            sage: p1 == p1
            True

            sage: p1 = HyperbolicPlane().PD().get_point(0)
            sage: p2 = HyperbolicPlane().PD().get_point(1/2 + 2*I/3)
            sage: p1 == p2
            False
            sage: p1 == p1
            True

            sage: p1 = HyperbolicPlane().KM().get_point((0,0))
            sage: p2 = HyperbolicPlane().KM().get_point((0, 1/2))
            sage: p1 == p2
            False

            sage: p1 = HyperbolicPlane().HM().get_point((0,0,1))
            sage: p2 = HyperbolicPlane().HM().get_point((0,0,1/1))
            sage: p1 == p2
            True
        """
        if not (isinstance(other, HyperbolicPoint)
                or self.parent() is other.parent()):
            return op == op_NE
        # bool is required to convert symbolic (in)equalities
        return bool(richcmp(self._coordinates, other._coordinates, op))

    def __rmul__(self, other):
        r"""
        Implement the action of matrices on points of hyperbolic space.

        EXAMPLES::

            sage: A = matrix(2, [0, 1, 1, 0])
            sage: A = HyperbolicPlane().UHP().get_isometry(A)
            sage: A * HyperbolicPlane().UHP().get_point(2 + I)
            Point in UHP 1/5*I + 2/5

        We also lift matrices into isometries::

            sage: B = diagonal_matrix([-1, -1, 1])
            sage: B = HyperbolicPlane().HM().get_isometry(B)                            # needs scipy
            sage: B * HyperbolicPlane().HM().get_point((0, 1, sqrt(2)))
            Point in HM (0, -1, sqrt(2))
        """
        if isinstance(other, HyperbolicIsometry):
            return other(self)
        elif isinstance(other, Matrix):
            # TODO: Currently the __mul__ from the matrices gets called first
            #    and returns an error instead of calling this method
            A = self.parent().get_isometry(other)
            return A(self)
        else:
            raise TypeError("unsupported operand type(s) for *:"
                            "{0} and {1}".format(self, other))

    #######################
    # Setters and Getters #
    #######################

    def coordinates(self):
        r"""
        Return the coordinates of the point.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(2 + I).coordinates()
            I + 2

            sage: HyperbolicPlane().PD().get_point(1/2 + 1/2*I).coordinates()
            1/2*I + 1/2

            sage: HyperbolicPlane().KM().get_point((1/3, 1/4)).coordinates()
            (1/3, 1/4)

            sage: HyperbolicPlane().HM().get_point((0,0,1)).coordinates()
            (0, 0, 1)
        """
        return self._coordinates

    def model(self):
        r"""
        Return the model to which the :class:`HyperbolicPoint` belongs.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(I).model()
            Hyperbolic plane in the Upper Half Plane Model

            sage: HyperbolicPlane().PD().get_point(0).model()
            Hyperbolic plane in the Poincare Disk Model

            sage: HyperbolicPlane().KM().get_point((0,0)).model()
            Hyperbolic plane in the Klein Disk Model

            sage: HyperbolicPlane().HM().get_point((0,0,1)).model()
            Hyperbolic plane in the Hyperboloid Model
        """
        return self.parent()

    def to_model(self, model):
        """
        Convert ``self`` to the ``model``.

        INPUT:

        - ``other`` -- (a string representing) the image model

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: PD.get_point(1/2+I/2).to_model(UHP)
            Point in UHP I + 2
            sage: PD.get_point(1/2+I/2).to_model('UHP')
            Point in UHP I + 2
        """
        if isinstance(model, str):
            model = getattr(self.parent().realization_of(), model)()
        return model(self)

    def is_boundary(self):
        """
        Return ``True`` if ``self`` is a boundary point.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: p = PD.get_point(0.5+.2*I)
            sage: p.is_boundary()
            False
            sage: p = PD.get_point(I)
            sage: p.is_boundary()
            True
        """
        return self._bdry

    def update_graphics(self, update=False, **options):
        r"""
        Update the graphics options of a :class:`HyperbolicPoint`.
        If ``update`` is ``True``, update rather than overwrite.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().get_point(I); p.graphics_options()
            {}

            sage: p.update_graphics(color = "red"); p.graphics_options()
            {'color': 'red'}

            sage: p.update_graphics(color = "blue"); p.graphics_options()
            {'color': 'blue'}

            sage: p.update_graphics(True, size = 20); p.graphics_options()
            {'color': 'blue', 'size': 20}
        """
        if not update:
            self._graphics_options = {}
        self._graphics_options.update(**options)

    def graphics_options(self):
        r"""
        Return the graphics options of the current point.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().get_point(2 + I, color='red')
            sage: p.graphics_options()
            {'color': 'red'}
        """
        return self._graphics_options

    def symmetry_involution(self):
        r"""
        Return the involutory isometry fixing the given point.

        EXAMPLES::

            sage: z = HyperbolicPlane().UHP().get_point(3 + 2*I)
            sage: z.symmetry_involution()
            Isometry in UHP
            [  3/2 -13/2]
            [  1/2  -3/2]

            sage: HyperbolicPlane().UHP().get_point(I).symmetry_involution()
            Isometry in UHP
            [ 0 -1]
            [ 1  0]

            sage: HyperbolicPlane().PD().get_point(0).symmetry_involution()
            Isometry in PD
            [-I  0]
            [ 0  I]

            sage: HyperbolicPlane().KM().get_point((0, 0)).symmetry_involution()
            Isometry in KM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: HyperbolicPlane().HM().get_point((0,0,1)).symmetry_involution()
            Isometry in HM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: p = HyperbolicPlane().UHP().random_element()
            sage: A = p.symmetry_involution()
            sage: p.dist(A*p)  # abs tol 1e-10
            0

            sage: A.preserves_orientation()
            True

            sage: A*A == HyperbolicPlane().UHP().get_isometry(identity_matrix(2))       # needs scipy
            True
        """
        R = self.parent().realization_of().a_realization()
        A = R(self).symmetry_involution()
        return self.parent().get_isometry(A)

    ###########
    # Display #
    ###########

    def show(self, boundary=True, **options):
        r"""
        Plot ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().PD().get_point(0).show()                            # needs sage.plot
            Graphics object consisting of 2 graphics primitives
            sage: HyperbolicPlane().KM().get_point((0,0)).show()                        # needs sage.plot
            Graphics object consisting of 2 graphics primitives
            sage: HyperbolicPlane().HM().get_point((0,0,1)).show()                      # needs sage.plot
            Graphics3d Object
        """
        p = self.coordinates()
        if p == infinity:
            raise NotImplementedError("can't draw the point infinity")

        opts = {'axes': False, 'aspect_ratio': 1}
        opts.update(self.graphics_options())
        opts.update(options)

        from sage.plot.point import point
        from sage.misc.functional import numerical_approx

        if self._bdry:  # It is a boundary point
            p = numerical_approx(p)
            pic = point((p, 0), **opts)
            if boundary:
                bd_pic = self._model.get_background_graphic(bd_min=p - 1,
                                                            bd_max=p + 1)
                pic = bd_pic + pic
        else:  # It is an interior point
            if p in RR:
                p = CC(p)
            else:
                p = numerical_approx(p)
            pic = point(p, **opts)
            if boundary:
                bd_pic = self.parent().get_background_graphic()
                pic = bd_pic + pic
        return pic


class HyperbolicPointUHP(HyperbolicPoint):
    r"""
    A point in the UHP model.

    INPUT:

    - the coordinates of a point in the unit disk in the complex plane `\CC`

    EXAMPLES::

        sage: HyperbolicPlane().UHP().get_point(2*I)
        Point in UHP 2*I

        sage: HyperbolicPlane().UHP().get_point(1)
        Boundary point in UHP 1
    """
    def symmetry_involution(self):
        r"""
        Return the involutory isometry fixing the given point.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(3 + 2*I).symmetry_involution()
            Isometry in UHP
            [  3/2 -13/2]
            [  1/2  -3/2]
        """
        p = self._coordinates
        x, y = real(p), imag(p)
        if y > 0:
            M = matrix([[x/y, -(x**2/y) - y], [1/y, -(x/y)]])
            return self.parent().get_isometry(M)
        raise ValueError("cannot determine the isometry of a boundary point")

    def show(self, boundary=True, **options):
        r"""
        Plot ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(I).show()
            Graphics object consisting of 2 graphics primitives
            sage: HyperbolicPlane().UHP().get_point(0).show()                           # needs sage.plot
            Graphics object consisting of 2 graphics primitives
            sage: HyperbolicPlane().UHP().get_point(infinity).show()
            Traceback (most recent call last):
            ...
            NotImplementedError: can...t draw the point infinity
        """
        p = self.coordinates()
        if p == infinity:
            raise NotImplementedError("can't draw the point infinity")
        opts = {'axes': False, 'aspect_ratio': 1}
        opts.update(self.graphics_options())
        opts.update(options)
        from sage.misc.functional import numerical_approx
        p = numerical_approx(p + 0 * I)
        from sage.plot.point import point
        if self._bdry:
            pic = point((p, 0), **opts)
            if boundary:
                bd_pic = self.parent().get_background_graphic(bd_min=p - 1,
                                                              bd_max=p + 1)
                pic = bd_pic + pic
        else:
            pic = point(p, **opts)
            if boundary:
                cent = real(p)
                bd_pic = self.parent().get_background_graphic(bd_min=cent - 1,
                                                              bd_max=cent + 1)
                pic = bd_pic + pic
        return pic
