# sage_setup: distribution = sagemath-plot
# sage.doctest: needs sage.symbolic
r"""
Graphics 3D object for representing and triangulating isosurfaces

AUTHORS:

- Robert Hanson (2007): initial Java version, in Jmol.
- Carl Witty (2009-01): first Cython version.
- Bill Cauchois (2009): improvements for inclusion into Sage.
"""

#*****************************************************************************
#      Copyright (C) 2009 Carl Witty <Carl.Witty@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Pieces of this file are strongly based on the marching cubes
# implementation in Jmol located at src/org/jmol/jvxl/calc/MarchingCubes.java.
# The data tables are in fact directly copied from there.

#*****************************************************************************
# This copyright is inherited from MarchingCubes.java in the Jmol source code.
#
#  * Copyright (C) 2007 Miguel, Bob, Jmol Development
#  *
#  * Contact: hansonr@stolaf.edu
#  *
#  *  This library is free software; you can redistribute it and/or
#  *  modify it under the terms of the GNU Lesser General Public
#  *  License as published by the Free Software Foundation; either
#  *  version 2.1 of the License, or (at your option) any later version.
#  *
#  *  This library is distributed in the hope that it will be useful,
#  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  *  Lesser General License for more details.
#  *
#  *  You should have received a copy of the GNU Lesser General Public
#  *  License along with this library; if not, write to the Free Software
#  *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
#*****************************************************************************

# There's a framework in here for computing multiple isosurfaces of a
# single function.  Currently, it's only used for a single
# implicit_plot3d where contour=... is given a list, but it's ready to
# be extended.  I think the best UI would be if animate(...) and
# show(...) had a prepass where they went through their rendering
# trees, found ImplicitSurface objects with the same function,
# bounding box, and plot_points (other parameters, such as contours,
# hole, color, would be allowed to be different),
# and arranged to have them all plotted simultaneously.  These
# prepasses have not been written yet.  Combining multiple
# ImplicitSurface plots would be particularly advantageous for animate(...),
# but for a big animation, watch out for memory usage.

# If you have a reasonably simple surface (not a space-filling fractal),
# then if n is your resolution, we have n^3 evaluations of the main
# function, about n^2 evaluations of auxiliary functions (hole, gradient,
# color), and output of size about n^2.

# With this in mind, we pay particular attention to optimizing the n^3
# function evaluations.  (But keep in mind that n may be as small as 20
# or so, so we shouldn't ignore the efficiency of the n^2 parts.)

# Also, we make sure not to ever allocate O(n^3) memory; we do the
# computation a slice at a time to avoid this.  (Jmol always allocates
# n^3 memory when it reads the JVXL file, but that might be on a different
# computer; Tachyon would only allocate memory proportional to the
# output size.)

cimport numpy as np
import numpy as np

from sage.plot.plot3d.transform cimport point_c, face_c, color_c, Transformation
from sage.plot.plot3d.index_face_set cimport IndexFaceSet
from sage.rings.real_double import RDF
from sage.plot.misc import setup_for_eval_on_grid
from sage.plot.colors import check_color_data

from libc.math cimport isnan

include "point_c.pxi"

# The default value for plot_points (a keyword argument to implicit_plot3d()),
# assumed when plot_points is set to "automatic".
DEFAULT_PLOT_POINTS = 40

cdef double nan = float(RDF('NaN'))

cdef inline bint marching_has_edge(double a, double b, double contour, double *f, bint *has_nan) noexcept:
    if isnan(a) or isnan(b):
        has_nan[0] = True
        return False

    has_nan[0] = False

    if (a >= contour) == (b >= contour):
        return False

    f[0] = (contour - a) / (b - a)
    return True

# Returns 0 or 1
cdef inline int marching_is_inside(double v, double contour) noexcept:
    return isnan(v) or v < contour

cdef void interpolate_point_c(point_c *result, double frac, point_c *inputs) noexcept:
    result[0].x = inputs[0].x + frac*(inputs[1].x - inputs[0].x)
    result[0].y = inputs[0].y + frac*(inputs[1].y - inputs[0].y)
    result[0].z = inputs[0].z + frac*(inputs[1].z - inputs[0].z)


cdef class VertexInfo:
    # The point in "integer space"
    cdef point_c pt

    # The gradient at this point in "evaluation space"
    cdef point_c gradient

    # (R,G,B) color ; a color triple
    cdef color_c color

    # This point in "evaluation space"
    cdef point_c eval_pt

    cdef void update_eval_pt(self, point_c *eval_min, point_c *eval_scale) noexcept:
        """
        Use eval_min and eval_scale to transform self.pt into evaluation space
        and store the result in self.eval_pt.
        """
        self.eval_pt.x = eval_min[0].x + eval_scale[0].x*self.pt.x
        self.eval_pt.y = eval_min[0].y + eval_scale[0].y*self.pt.y
        self.eval_pt.z = eval_min[0].z + eval_scale[0].z*self.pt.z

    def __repr__(self):
        """
        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import VertexInfo
            sage: VertexInfo()
            <0.0, 0.0, 0.0>
        """
        return '<{}, {}, {}>'.format(self.pt.x, self.pt.y, self.pt.z)


cdef mk_VertexInfo(double x, double y, double z, point_c *eval_min, point_c *eval_scale):
    cdef VertexInfo v
    v = VertexInfo.__new__(VertexInfo)
    v.pt.x = x
    v.pt.y = y
    v.pt.z = z

    v.update_eval_pt(eval_min, eval_scale)

    return v


cdef class MarchingCubes:
    r"""
    Handles marching cube rendering.

    Protocol:

    1. Create the class.
    2. Call ``process_slice`` once for each X slice, from self.nx > x >= 0.
    3. Call ``finish()``, which returns a list of strings.

    .. NOTE::

        Actually, only 4 slices ever exist; the caller will re-use old
        storage.
    """

    cdef readonly object xrange
    cdef readonly object yrange
    cdef readonly object zrange
    cdef readonly double contour
    cdef readonly int nx
    cdef readonly int ny
    cdef readonly int nz
    cdef readonly Transformation transform
    cdef readonly object region
    cdef readonly object gradient
    cdef readonly bint smooth
    cdef readonly object color_function
    cdef readonly object colormap
    cdef readonly object results

    # We deal with three coordinate systems.  We do most of our work
    # in an integral coordinate system where x ranges from
    # 0 <= x < self.nx, etc.; we do function evaluations where
    # self.xrange[0] <= x <= self.xrange[1], etc.; and output
    # is in a third coordinate system.
    # (Note that in "integer space", function evaluations of the main
    # function happen on integer coordinates; but function evaluations
    # of the auxiliary functions will have one non-integer coordinate.)

    # These parameters convert from integer space to function-evaluation
    # space: eval_coord = eval_min + int_coord * eval_scale
    cdef point_c eval_min
    cdef point_c eval_scale

    # The componentwise reciprocal of eval_scale; just used to change
    # some divisions into multiplications
    cdef point_c eval_scale_inv

    cdef point_c out_origin, out_plus_x, out_plus_y, out_plus_z

    def __init__(self, xrange, yrange, zrange, contour, plot_points,
                 transform=None, region=None, gradient=None, smooth=True,
                 color_function=None, colormap=None):
        """
        TESTS:

        Marching cubes is an abstract base class, you can't do anything with it::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubes
            sage: cube_marcher = MarchingCubes((0, 1), (0, 1), (0, 1), 1, (10, 10, 10))
        """
        self.xrange = xrange
        self.yrange = yrange
        self.zrange = zrange
        self.contour = contour
        self.nx = plot_points[0]
        self.ny = plot_points[1]
        self.nz = plot_points[2]
        self.transform = transform
        self.region = region
        self.gradient = gradient
        self.smooth = smooth

        self.color_function = color_function
        self.colormap = colormap

        self.eval_min.x = xrange[0]
        self.eval_scale.x = (xrange[1] - xrange[0]) / (self.nx - 1)
        self.eval_min.y = yrange[0]
        self.eval_scale.y = (yrange[1] - yrange[0]) / (self.ny - 1)
        self.eval_min.z = zrange[0]
        self.eval_scale.z = (zrange[1] - zrange[0]) / (self.nz - 1)
        self.eval_scale_inv.x = 1/self.eval_scale.x
        self.eval_scale_inv.y = 1/self.eval_scale.y
        self.eval_scale_inv.z = 1/self.eval_scale.z

        cdef point_c zero_pt, origin, plus_x, plus_y, plus_z
        zero_pt.x = 0
        zero_pt.y = 0
        zero_pt.z = 0
        origin = self.eval_min
        plus_x = zero_pt
        plus_x.x = self.eval_scale.x
        plus_y = zero_pt
        plus_y.y = self.eval_scale.y
        plus_z = zero_pt
        plus_z.z = self.eval_scale.z
        if self.transform is not None:
            self.transform.transform_point_c(&self.out_origin, origin)
            self.transform.transform_point_c(&self.out_plus_x, plus_x)
            self.transform.transform_point_c(&self.out_plus_y, plus_y)
            self.transform.transform_point_c(&self.out_plus_z, plus_z)
        else:
            self.out_origin = origin
            self.out_plus_x = plus_x
            self.out_plus_y = plus_y
            self.out_plus_z = plus_z

        self.results = []

    def finish(self):
        """
        Return the results of the marching cubes algorithm as a list.

        The format is specific to the subclass implementing this
        method.

        TESTS:

        By default, it returns an empty list::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubes
            sage: cube_marcher = MarchingCubes((0, 1), (0, 1), (0, 1), 1, (10, 10, 10), None)
            sage: cube_marcher.finish()
            []
        """
        return self.results


cdef class MarchingCubesTriangles(MarchingCubes):
    """
    A subclass of MarchingCubes that returns its results as a list of triangles,
    including their vertices and normals (if ``smooth=True``).

    And also their vertex colors if a vertex coloring function is given.
    """

    cdef readonly np.ndarray x_vertices
    cdef readonly np.ndarray y_vertices
    cdef readonly np.ndarray z_vertices

    cdef readonly np.ndarray y_vertices_swapped
    cdef readonly np.ndarray z_vertices_swapped

    cdef readonly slices

    def __init__(self, *args, **kwargs):
        """
        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubesTriangles
            sage: cube_marcher = MarchingCubesTriangles((0, 1), (0, 1), (0, 1), 1, (10, 10, 10))
        """

        MarchingCubes.__init__(self, *args, **kwargs)

        self.x_vertices = np.empty((self.ny, self.nz), dtype=object)
        self.y_vertices = np.empty((2, self.ny-1, self.nz), dtype=object)
        self.z_vertices = np.empty((2, self.ny, self.nz-1), dtype=object)

        # Create new arrays that share the same underlying data, but
        # have 0 and 1 reversed for the first coordinate.
        self.y_vertices_swapped = self.y_vertices[::-1, ...]
        self.z_vertices_swapped = self.z_vertices[::-1, ...]

        self.slices = []

    def process_slice(self, unsigned int x, np.ndarray slice):
        """
        Process a single slice of function evaluations at the specified `x`
        coordinate.

        EXAMPLES::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubesTriangles
            sage: import numpy as np
            sage: cube_marcher = MarchingCubesTriangles((-2, 2), (-2, 2), (-2, 2), 4, (10,)*3, smooth=False)
            sage: f = lambda x, y, z: x^2 + y^2 + z^2
            sage: slices = np.zeros((10, 10, 10), dtype=np.double)
            sage: for x in reversed(range(0, 10)):
            ....:     for y in range(0, 10):
            ....:         for z in range(0, 10):
            ....:             slices[x, y, z] = f(*[a * (4 / 9) -2 for a in (x, y, z)])
            ....:     cube_marcher.process_slice(x, slices[x, :, :])
            sage: faces = cube_marcher.finish()
            sage: faces[0][0]
            {'x': 1.555555555555..., 'y': -1.111111111111..., 'z': -0.555555555555...}

        We render the isosurface using IndexFaceSet::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: IndexFaceSet([tuple((p['x'], p['y'], p['z']) for p in face) for face in faces])
            Graphics3d Object
        """
        # We go to a lot of effort to avoid repeating computations.
        # (I hope that the related bookkeeping is actually faster
        # than repeating the computations, but I haven't tested.)

        # We get evaluation points, one slice at a time.  But we
        # don't want to process slices of evaluation points;
        # we want to process slices of cubes.  (There is one fewer slice
        # of cubes than there is of evaluation points.)

        # Without any caching, we would repeat this for each cube:
        # Find which vertices are inside/outside the surface.
        # For edges which cross the surface, interpolate to
        #   find the exact crossing location; these will be the vertices
        #   of triangles.
        # Compute the color and surface normal for each of these vertices.
        # Assemble the vertices into triangles, using the "marching cubes"
        #   tables.

        # The problem with this is that each vertex is actually shared
        # among four neighbor cubes (except on the edge of the array),
        # so we would be repeating work (potentially a lot of work, if
        # the user-specified gradient or color routines are expensive)
        # four times.

        # So we cache the information associated with each vertex.
        # Let's call an edge an X, Y, or Z edge depending on which
        # axis it is parallel to; and a vertex is an X, Y, or Z
        # vertex depending on which kind of edge it lies on.

        # X vertices are shared between cubes that are all within a
        # single cube slice.  However, Y and Z vertices are shared
        # between adjacent slices, so we need to keep those caches
        # around.  But we still need only two slices of vertex cache
        # at a time (the two that are adjacent to the current cube slice).

        # To reduce memory allocations, we allocate these caches at the
        # beginning of the run (they are ndarrays of objects).  We
        # set a VertexInfo object in the cache when we first compute it;
        # then we look in the cache when we need it for the other three
        # cubes.

        # The X vertex cache is a 2-D ndarray.  The Y and Z vertex caches
        # are 3-D, where the first (x) dimension is indexed by 0 or 1.

        # The Cython buffer interface (that we use for fast access to
        # ndarray's) has to reinitialize itself separately in each
        # function.  For this reason, we use fewer, bigger functions;
        # in particular, we don't call any functions per-cube that
        # require buffer access.

        # To compute interpolated gradients using central differencing,
        # we need up to 4 slices.  Consider a gradient computed at
        # an X vertex, between slices 1 and 2.  This is interpolated
        # between the gradient at slice 1 and the gradient at slice 2.
        # Computing the gradient at slice 1 requires slices 0 and 2;
        # computing the gradient at slice 2 requires slices 1 and 3.

        # This means we need to queue up slices and process them
        # in a somewhat strange order.  This function only does
        # the queuing, and then passes all the real work off to
        # other functions.

        self.slices = ([slice] + self.slices)[:4]

        if len(self.slices) >= 2:
            self._update_yz_vertices(x+1,
                                    self.slices[0],
                                    self.slices[1],
                                    self.slices[2] if len(self.slices) >= 3 else None)

        if len(self.slices) >= 3:
            self._update_x_vertices(x+1,
                                   self.slices[0],
                                   self.slices[1],
                                   self.slices[2],
                                   self.slices[3] if len(self.slices) >= 4 else None)
            self.process_cubes(self.slices[1], self.slices[2])

        if x == 0:
            self._update_yz_vertices(x,
                                    None,
                                    self.slices[0],
                                    self.slices[1])

            self._update_x_vertices(x,
                                   None,
                                   self.slices[0],
                                   self.slices[1],
                                   self.slices[2] if len(self.slices) >= 3 else None)

            self.process_cubes(self.slices[0], self.slices[1])

    cpdef _update_yz_vertices(self, int x, np.ndarray _prev, np.ndarray _cur, np.ndarray _next):
        """
        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubesTriangles
            sage: import numpy as np
            sage: cube_marcher = MarchingCubesTriangles((0, 1), (0, 1), (0, 1), 0, (2,)*3, smooth=False)
            sage: prev_slice = next_slice = np.ones((2, 2), dtype=np.double)
            sage: cur_slice = prev_slice.copy()
            sage: cur_slice[0, 0] = -1 # Seed the slice data with an "interesting" value.
            sage: cube_marcher._update_yz_vertices(1, prev_slice, cur_slice, next_slice)
            sage: cube_marcher.z_vertices.tolist()
            [[[<1.0, 0.0, 0.5>], [None]], [[None], [None]]]
            sage: cube_marcher.y_vertices.tolist()
            [[[<1.0, 0.5, 0.0>, None]], [[None, None]]]
            sage: cube_marcher.x_vertices.any() # This shouldn't affect the X vertices.
            ...
        """
        (self.y_vertices, self.y_vertices_swapped) = \
            (self.y_vertices_swapped, self.y_vertices)
        (self.z_vertices, self.z_vertices_swapped) = \
            (self.z_vertices_swapped, self.z_vertices)

        cdef bint has_prev = (_prev is not None)
        cdef bint has_next = (_next is not None)

        cdef np.ndarray[double, ndim=2] prev = _prev
        cdef np.ndarray[double, ndim=2] cur = _cur
        cdef np.ndarray[double, ndim=2] next = _next

        cdef np.ndarray[object, ndim=2] y_vertices = self.y_vertices[0,...]
        cdef np.ndarray[object, ndim=2] z_vertices = self.z_vertices[0,...]

        cdef int ny = self.ny
        cdef int nz = self.nz

        cdef int y
        cdef int z
        cdef VertexInfo v
        cdef double frac
        cdef bint has_nan
        cdef point_c gradients[2]
        cdef int i
        for y in range(ny - 1):
            for z in range(nz):
                if marching_has_edge(cur[y,z], cur[y+1,z], self.contour, &frac, &has_nan):
                    v = mk_VertexInfo(x, y+frac, z, &self.eval_min, &self.eval_scale)
                    if self.region is not None:
                        if not self.in_region(v):
                            y_vertices[y,z] = None
                            continue
                    if self.smooth:
                        # We must compute a gradient.
                        if self.gradient is not None:
                            self.apply_point_func(&v.gradient, self.gradient, v)
                        else:
                            # Use central differencing.
                            for i in range(2):
                                self.get_gradient(&gradients[i],
                                                   x, y+i, z,
                                                   cur[y+i,z],
                                                   prev[y+i,z] if has_prev else 0,
                                                   next[y+i,z] if has_next else 0,
                                                   cur[y+i-1,z] if y+i>0 else 0,
                                                   cur[y+i+1,z] if y+i<ny-1 else 0,
                                                   cur[y+i,z-1] if z>0 else 0,
                                                   cur[y+i,z+1] if z<nz-1 else 0)
                            interpolate_point_c(&v.gradient, frac, gradients)
                    if not(self.color_function is None):
                        self.apply_color_func(&v.color, self.color_function,
                                              self.colormap, v)
                    y_vertices[y,z] = <object>v
                else:
                    y_vertices[y,z] = None

        # OK, that updated the Y vertices.  Now we do almost exactly
        # the same thing to update Z vertices.
        for y in range(ny):
            for z in range(nz - 1):
                if marching_has_edge(cur[y,z], cur[y,z+1], self.contour, &frac, &has_nan):
                    v = mk_VertexInfo(x, y, z+frac, &self.eval_min, &self.eval_scale)
                    if self.region is not None:
                        if not self.in_region(v):
                            z_vertices[y,z] = None
                            continue
                    if self.smooth:
                        # We must compute a gradient.
                        if self.gradient is not None:
                            self.apply_point_func(&v.gradient, self.gradient, v)
                        else:
                            # Use central differencing.
                            for i in range(2):
                                self.get_gradient(&gradients[i],
                                                   x, y, z+i,
                                                   cur[y,z+i],
                                                   prev[y,z+i] if has_prev else 0,
                                                   next[y,z+i] if has_next else 0,
                                                   cur[y-1,z+i] if y>0 else 0,
                                                   cur[y+1,z+i] if y<ny-1 else 0,
                                                   cur[y,z+i-1] if z+i>0 else 0,
                                                   cur[y,z+i+1] if z+i<nz-1 else 0)
                            interpolate_point_c(&v.gradient, frac, gradients)
                    if not(self.color_function is None):
                        self.apply_color_func(&v.color, self.color_function,
                                              self.colormap, v)
                    z_vertices[y,z] = <object>v
                else:
                    z_vertices[y,z] = None

    cpdef _update_x_vertices(self, int x, np.ndarray _prev, np.ndarray _left, np.ndarray _right, np.ndarray _next):
        """
        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubesTriangles
            sage: import numpy as np
            sage: cube_marcher = MarchingCubesTriangles((0, 1), (0, 1), (0, 1), 0, (4, 2, 2), smooth=False)
            sage: prev_slice = right_slice = next_slice = np.ones((2, 2), dtype=np.double)
            sage: left_slice = prev_slice.copy()
            sage: left_slice[1, 1] = -1
            sage: cube_marcher._update_x_vertices(1, prev_slice, left_slice, right_slice, next_slice)
            sage: cube_marcher.x_vertices.tolist()
            [[None, None], [None, <1.5, 1.0, 1.0>]]
            sage: cube_marcher.y_vertices.any() or cube_marcher.z_vertices.any() # This shouldn't affect the Y or Z vertices.
            ...
        """
        cdef bint has_prev = (_prev is not None)
        cdef bint has_next = (_next is not None)

        cdef np.ndarray[double, ndim=2] prev = _prev
        cdef np.ndarray[double, ndim=2] left = _left
        cdef np.ndarray[double, ndim=2] right = _right
        cdef np.ndarray[double, ndim=2] next = _next

        cdef np.ndarray[object, ndim=2] x_vertices = self.x_vertices

        cdef int ny = self.ny
        cdef int nz = self.nz

        cdef int y
        cdef int z

        cdef VertexInfo v
        cdef double frac
        cdef bint has_nan
        cdef point_c gradients[2]
        for y in range(ny):
            for z in range(nz):
                if marching_has_edge(left[y,z], right[y,z], self.contour, &frac, &has_nan):
                    v = mk_VertexInfo(x+frac, y, z, &self.eval_min, &self.eval_scale)
                    if self.region is not None:
                        if not self.in_region(v):
                            x_vertices[y,z] = None
                            continue
                    if self.smooth:
                        # We must compute a gradient.
                        if self.gradient is not None:
                            self.apply_point_func(&v.gradient, self.gradient, v)
                        else:
                            # Use central differencing.
                            self.get_gradient(&gradients[0],
                                               x, y, z,
                                               left[y,z],
                                               prev[y,z] if has_prev else 0,
                                               right[y,z],
                                               left[y-1,z] if y>0 else 0,
                                               left[y+1,z] if y<ny-1 else 0,
                                               left[y,z-1] if z>0 else 0,
                                               left[y,z+1] if z<nz-1 else 0)
                            self.get_gradient(&gradients[1],
                                               x, y, z,
                                               right[y,z],
                                               left[y,z],
                                               next[y,z] if has_next else 0,
                                               right[y-1,z] if y>0 else 0,
                                               right[y+1,z] if y<ny-1 else 0,
                                               right[y,z-1] if z>0 else 0,
                                               right[y,z+1] if z<nz-1 else 0)
                            interpolate_point_c(&v.gradient, frac, gradients)
                    if not(self.color_function is None):
                        self.apply_color_func(&v.color, self.color_function,
                                              self.colormap, v)
                    x_vertices[y,z] = <object>v
                else:
                    x_vertices[y,z] = None

    cdef bint in_region(self, VertexInfo v) noexcept:
        return (self.region(v.eval_pt.x, v.eval_pt.y, v.eval_pt.z) > 0)

    cdef apply_point_func(self, point_c *pt, fn, VertexInfo v):
        if isinstance(fn, tuple):
            pt[0].x = fn[0](v.eval_pt.x, v.eval_pt.y, v.eval_pt.z)
            pt[0].y = fn[1](v.eval_pt.x, v.eval_pt.y, v.eval_pt.z)
            pt[0].z = fn[2](v.eval_pt.x, v.eval_pt.y, v.eval_pt.z)
        else:
            t = fn(v.eval_pt.x, v.eval_pt.y, v.eval_pt.z)
            pt[0].x = t[0]
            pt[0].y = t[1]
            pt[0].z = t[2]

    cdef apply_color_func(self, color_c *pt, fn, cm, VertexInfo v):
        t = fn(v.eval_pt.x, v.eval_pt.y, v.eval_pt.z)
        pt[0].r, pt[0].g, pt[0].b, _ = cm(t)

    cdef get_gradient(self,
                      point_c *g,
                      int x, int y, int z,
                      double center,
                      double lx, double ux,
                      double ly, double uy,
                      double lz, double uz):
        # What a mess!  It would be much nicer-looking code to pass slices
        # in here and do the subscripting in here.  Unfortunately,
        # that would also be much slower, because we'd have to re-initialize
        # the Cython buffer interface on each call.

        cdef double dx = ux - lx
        cdef double dy = uy - ly
        cdef double dz = uz - lz

        cdef double gx = dx * self.eval_scale_inv.x
        cdef double gy = dy * self.eval_scale_inv.y
        cdef double gz = dz * self.eval_scale_inv.z

        if x > 0 and x < self.nx - 1:
            gx *= 0.5
        if y > 0 and y < self.ny - 1:
            gy *= 0.5
        if z > 0 and z < self.nz - 1:
            gz *= 0.5

        g[0].x = gx
        g[0].y = gy
        g[0].z = gz

    cpdef process_cubes(self, np.ndarray _left, np.ndarray _right):
        """
        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubesTriangles
            sage: import numpy as np
            sage: cube_marcher = MarchingCubesTriangles((0, 1), (0, 1), (0, 1), 0, (3, 2, 2), smooth=False)
            sage: slices = [np.ones((2, 2), dtype=np.double) for i in range(3)]
            sage: slices[0][1, 1] = -1
            sage: cube_marcher._update_yz_vertices(0, None, slices[0], slices[1])
            sage: cube_marcher._update_x_vertices(0, None, slices[0], slices[1], slices[2])
            sage: cube_marcher.process_cubes(slices[0], slices[1])
            sage: cube_marcher.finish()
            [({'x': 0.0, 'y': 1.0, 'z': 0.5},
              {'x': 0.25, 'y': 1.0, 'z': 1.0},
              {'x': 0.0, 'y': 0.5, 'z': 1.0})]
        """
        cdef np.ndarray[double, ndim=2] left = _left
        cdef np.ndarray[double, ndim=2] right = _right

        cdef np.ndarray[object, ndim=2] x_vertices = self.x_vertices
        cdef np.ndarray[object, ndim=3] y_vertices = self.y_vertices
        cdef np.ndarray[object, ndim=3] z_vertices = self.z_vertices

        cdef int ny = self.ny
        cdef int nz = self.nz

        cdef int y
        cdef int z

        # based on generateSurfaceData in MarchingCubes.java
        cdef int insideMask

        # Cool ASCII art from MarchingCubes.java:
#      *                     Y
#      *                      4 --------4--------- 5
#      *                     /|                   /|
#      *                    / |                  / |
#      *                   /  |                 /  |
#      *                  7   8                5   |
#      *                 /    |               /    9
#      *                /     |              /     |
#      *               7 --------6--------- 6      |
#      *               |      |             |      |
#      *               |      0 ---------0--|----- 1    X
#      *               |     /              |     /
#      *              11    /               10   /
#      *               |   3                |   1
#      *               |  /                 |  /
#      *               | /                  | /
#      *               3 ---------2-------- 2
#      *              Z

        # We see from the above that vertices are labeled 0 to 7, and
        # edges are labeled 0 to 11.

        cdef list all_vertex_info = [None] * 12
        cdef tuple my_triangles

        cdef int i

        for y in range(ny - 1):
            for z in range(nz - 1):
                # For each vertex (0 to 7), set the corresponding bit
                # of insideMask iff the vertex is inside the surface.
                insideMask = 0
                insideMask |= marching_is_inside(left[y,z], self.contour)<<0
                insideMask |= marching_is_inside(right[y,z], self.contour)<<1
                insideMask |= marching_is_inside(right[y,z+1], self.contour)<<2
                insideMask |= marching_is_inside(left[y,z+1], self.contour)<<3
                insideMask |= marching_is_inside(left[y+1,z], self.contour)<<4
                insideMask |= marching_is_inside(right[y+1,z], self.contour)<<5
                insideMask |= marching_is_inside(right[y+1,z+1], self.contour)<<6
                insideMask |= marching_is_inside(left[y+1,z+1], self.contour)<<7

                if insideMask == 0 or insideMask == 255:
                    continue

                # OK, we have a cube on the surface.  Copy all of the vertex
                # info into an array for easier reference.

                all_vertex_info[0] = x_vertices[y,z]
                all_vertex_info[1] = z_vertices[1,y,z]
                all_vertex_info[2] = x_vertices[y,z+1]
                all_vertex_info[3] = z_vertices[0,y,z]
                all_vertex_info[4] = x_vertices[y+1,z]
                all_vertex_info[5] = z_vertices[1,y+1,z]
                all_vertex_info[6] = x_vertices[y+1,z+1]
                all_vertex_info[7] = z_vertices[0,y+1,z]
                all_vertex_info[8] = y_vertices[0,y,z]
                all_vertex_info[9] = y_vertices[1,y,z]
                all_vertex_info[10]= y_vertices[1,y,z+1]
                all_vertex_info[11]= y_vertices[0,y,z+1]

                my_triangles = triangle_table2[insideMask]

                for i in range(0, len(my_triangles), 4):
                    # In wireframe mode, my_triangles[i+3] specifies
                    # whether or not to draw the corresponding edge.
                    # See MarchingCubes.java for details.
                    self.add_triangle(all_vertex_info[my_triangles[i]],
                                      all_vertex_info[my_triangles[i+1]],
                                      all_vertex_info[my_triangles[i+2]])

    cpdef add_triangle(self, VertexInfo v1, VertexInfo v2, VertexInfo v3):
        """
        Called when a new triangle is generated by the marching cubes algorithm
        to update the results array.

        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import MarchingCubesTriangles, VertexInfo
            sage: cube_marcher = MarchingCubesTriangles((0, 1), (0, 1), (0, 1), 0, (10,)*3, smooth=False)
            sage: cube_marcher.add_triangle(VertexInfo(), VertexInfo(), VertexInfo())
            sage: cube_marcher.finish()
            [({'x': 0.0, 'y': 0.0, 'z': 0.0},
              {'x': 0.0, 'y': 0.0, 'z': 0.0},
              {'x': 0.0, 'y': 0.0, 'z': 0.0})]
        """
        if v1 is None or v2 is None or v3 is None:
            # This happens if there is a NaN nearby, or if a hole was
            # specified here.
            return

        cdef:
            point_c v1_ev_pt, v2_ev_pt, v3_ev_pt
            point_c n1_ev_vec, n2_ev_vec, n3_ev_vec
            color_c v1_col, v2_col, v3_col

        if self.transform is not None:
            self.transform.transform_point_c(&v1_ev_pt, v1.eval_pt)
            self.transform.transform_point_c(&v2_ev_pt, v2.eval_pt)
            self.transform.transform_point_c(&v3_ev_pt, v3.eval_pt)
        else:
            v1_ev_pt = v1.eval_pt
            v2_ev_pt = v2.eval_pt
            v3_ev_pt = v3.eval_pt

        face = (v1_ev_pt, v2_ev_pt, v3_ev_pt)

        if not(self.color_function is None):
            v1_col = v1.color
            v2_col = v2.color
            v3_col = v3.color
            face += (v1_col, v2_col, v3_col)

        if self.smooth:
            # XXX I believe this is wrong for non-uniform transforms
            if self.transform is not None:
                self.transform.transform_vector_c(&n1_ev_vec, v1.gradient)
                self.transform.transform_vector_c(&n2_ev_vec, v2.gradient)
                self.transform.transform_vector_c(&n3_ev_vec, v3.gradient)
            else:
                n1_ev_vec = v1.gradient
                n2_ev_vec = v2.gradient
                n3_ev_vec = v3.gradient
            face += (n1_ev_vec, n2_ev_vec, n3_ev_vec)

        self.results.append(face)


cpdef render_implicit(f, xrange, yrange, zrange, plot_points, cube_marchers):
    """
    INPUT:

    - ``f`` -- a (fast!) callable function

    - ``xrange`` -- a 2-tuple (x_min, x_max)

    - ``yrange`` -- a 2-tuple (y_min, y_may)

    - ``zrange`` -- a 2-tuple (z_min, z_maz)

    - ``plot_points`` -- a triple of integers indicating the number of
      function evaluations in each direction

    - ``cube_marchers`` -- list of cube marchers, one for each contour

    OUTPUT:

    A representation of the isosurface, in the format specified by the
    individual cube marchers.

    TESTS::

        sage: from sage.plot.plot3d.implicit_surface import render_implicit, MarchingCubesTriangles
        sage: plot_points, f = (40,)*3, lambda x, y, z: x + y + z
        sage: cube_marcher = MarchingCubesTriangles((0, 1), (0, 1), (0, 1), 1, (10,)*3)
        sage: results = render_implicit(lambda x, y, z: x + y + z,
        ....:                           (0, 1), (0, 1), (0, 1), (10,)*3, [cube_marcher])
        sage: results[0][0]
        {'x': 1.0, 'y': 0.0, 'z': 0.0}
    """

    cdef int nx = plot_points[0]
    cdef int ny = plot_points[1]
    cdef int nz = plot_points[2]

    cdef point_c eval_min, eval_scale

    eval_min.x = xrange[0]
    eval_scale.x = (xrange[1] - xrange[0]) / (nx - 1)
    eval_min.y = yrange[0]
    eval_scale.y = (yrange[1] - yrange[0]) / (ny - 1)
    eval_min.z = zrange[0]
    eval_scale.z = (zrange[1] - zrange[0]) / (nz - 1)

    # A possible future extension would be to allow passing in "f"
    # as a numpy ndarray.  If that were done, we could slightly modify
    # the following code to just pass slices of f to the renderers
    # (no copying of the data would be required).

    # The current marching cube renderers need only at most four slices at
    # a time.

    cdef np.ndarray[double, ndim=3] data = np.zeros((4, ny, nz), dtype=np.double)
    cdef np.ndarray[double, ndim=2] slice

    cdef unsigned int x, y, z
    cdef int n

    cdef double eval_x, eval_y, eval_z

    cdef MarchingCubes marcher

    for n in range(nx):
        x = nx-1-n
        eval_x = eval_min.x + eval_scale.x * x
        slice = data[n % 4, :, :]
        for y in range(ny):
            eval_y = eval_min.y + eval_scale.y * y
            for z in range(nz):
                eval_z = eval_min.z + eval_scale.z * z
                slice[y, z] = f(eval_x, eval_y, eval_z)

        for marcher in cube_marchers:
            marcher.process_slice(x, slice)

    results = []

    for marcher in cube_marchers:
        results.extend(marcher.finish())

    return results


cdef class ImplicitSurface(IndexFaceSet):
    cdef readonly object f
    cdef readonly object vars
    cdef readonly tuple xrange
    cdef readonly tuple yrange
    cdef readonly tuple zrange
    cdef readonly list contours
    cdef readonly object region
    cdef readonly bint smooth
    cdef readonly object gradient
    cdef readonly object color_function
    cdef readonly object colormap
    cdef readonly tuple plot_points

    def __init__(self, f, xrange, yrange, zrange,
                 contour=0, plot_points='automatic',
                 region=None, smooth=True, gradient=None,
                 **kwds):
        """
        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: var('x,y,z')
            (x, y, z)
            sage: G = ImplicitSurface(x^2 + y^2 + z^2, (x,-2, 2), (y,-2, 2), (z,-2, 2), contour=4)
            sage: show(G)

        A colored case::

            sage: t = (1-sin(2*x*y+3*z)**2).function(x,y,z)
            sage: cm = colormaps.autumn
            sage: G = ImplicitSurface(x^2 + y^2 + z^2, (x,-2, 2), (y,-2, 2), (z,-2, 2), contour=4, color=(t,cm))
            sage: G.show(viewer='tachyon')
        """
        self._extra_kwds = kwds
        color_data = None
        if 'color' in kwds:
            try:
                if len(kwds['color']) == 2 and callable(kwds['color'][0]):
                    color_data = kwds['color']
                    kwds.pop('color')
            except (TypeError, AttributeError):
                pass
        if color_data is None:
            # case of a global color
            self.color_function = None
            IndexFaceSet.__init__(self, [], [], **kwds)
        else:
            # case of a color depending on parameters
            cf, cm = check_color_data(color_data)
            self.color_function = cf
            self.colormap = cm
            IndexFaceSet.__init__(self, [], [], texture_list=[], **kwds)
        from sage.ext.fast_eval import fast_float

        orig_f = f
        self.f, ranges, self.vars = setup_for_eval_on_grid(f, [xrange, yrange, zrange], return_vars=True)
        self.xrange = ranges[0][:2]
        self.yrange = ranges[1][:2]
        self.zrange = ranges[2][:2]
        if isinstance(contour, (list, tuple)):
            contours = contour
        else:
            contours = [contour]
        self.contours = [float(c) for c in contours]
        if region is not None:
            self.region = fast_float(region, *self.vars)

        # Comments from Carl Witty, who first wrote this some of this code
        # See Issue 9483
        # When I first wrote the code, I had the idea to create a
        # direct-to-tachyon backend that would use vertex normals
        # to create much nicer-looking plots with smaller numbers
        # of plot_points, and laid the groundwork for this backend
        # with the gradient and smooth arguments. But I abandoned the
        # code without writing this backend (and leaving many other parts
        # of the code unfinished).
        # When William Cauchois took over and finished the code (thank you,
        # William!), he only wrote an IndexFaceSet backend, that can't
        # (currently) make use of vertex normals. So the gradient code is
        # currently useless.
        # But it's still open for somebody to write a direct-to-tachyon backend,
        # or to extend IndexFaceSet to support vertex normals.

        # Since IndexFaceSet doesn't even support smooth shading, we overwrite
        # the passed-in smooth parameter.
        smooth=False

        self.smooth = smooth
        if smooth and gradient is None:
            try:
                gradient = (orig_f.diff(self.vars[0]),
                            orig_f.diff(self.vars[1]),
                            orig_f.diff(self.vars[2]))
            except Exception:
                # Would be nice to have more nuanced error handling here.

                # If anything goes wrong, we'll just use central differencing.
                pass
        if gradient is not None:
            self.gradient = fast_float(gradient, *self.vars)
        if self.color_function is not None:
            self.color_function = fast_float(self.color_function, *self.vars)
        if plot_points == "automatic":
            plot_points = DEFAULT_PLOT_POINTS
        my_plot_points = []
        for i in range(3):
            if isinstance(plot_points, (list, tuple)):
                n = int(plot_points[i])
            else:
                n = int(plot_points)
            if n < 2:
                raise ValueError
            my_plot_points.append(n)
        self.plot_points = tuple(my_plot_points)

    def bounding_box(self):
        """
        Return a bounding box for the ``ImplicitSurface``, as a tuple of two
        3-dimensional points.

        EXAMPLES:

        Note that the bounding box corresponds exactly to the x-, y-, and z- range::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: G = ImplicitSurface(0, (0, 1), (0, 1), (0, 1))
            sage: G.bounding_box()
            ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        """
        return ((self.xrange[0], self.yrange[0], self.zrange[0]),
                (self.xrange[1], self.yrange[1], self.zrange[1]))

    def obj_repr(self, render_params):
        """
        Return a representation of this object in the .obj format.

        TESTS:

        We graph a simple plane::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: var('x,y,z')
            (x, y, z)
            sage: G = ImplicitSurface(x + y + z, (x,-1, 1), (y,-1, 1), (z,-1, 1))
            sage: obj = G.obj_repr(G.default_render_params())
            sage: vertices = obj[2]

        The number of vertices in the OBJ representation should equal the number
        of vertices in the face set::

            sage: len(vertices) == len(G.vertex_list())
            True

        The vertices in the OBJ representation should also be approximately equal
        to the vertices in the face set -- the small error is due to rounding
        which occurs during output (we test only the first 20 points for the
        sake of speed)::

            sage: def points_equal(a, b, epsilon=(1e-5)):
            ....:     return all(abs(x0-x1) < epsilon for x0, x1 in zip(a, b))
            sage: checklist = []
            sage: assert len(vertices) >= 20 # I should hope so, we're rendering at the default resolution!
            sage: for vertex, surf_vertex in list(zip(vertices, G.vertex_list()))[0:20]:
            ....:     checklist.append(points_equal(map(float, vertex.split(' ')[1:]), surf_vertex))
            sage: all(checklist)
            True
        """
        self.triangulate()
        return IndexFaceSet.obj_repr(self, render_params)

    def tachyon_repr(self, render_params):
        """
        Return a representation of this object suitable for use with the Tachyon
        renderer.

        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: var('x,y,z')
            (x, y, z)
            sage: G = ImplicitSurface(x + y + z, (x,-1, 1), (y,-1, 1), (z,-1, 1))
            sage: G.tachyon_repr(G.default_render_params())[0].startswith('TRI')
            True
        """
        self.triangulate()
        return IndexFaceSet.tachyon_repr(self, render_params)

    def jmol_repr(self, render_params):
        """
        Return a representation of this object suitable for use with the Jmol
        renderer.

        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: var('x,y,z')
            (x, y, z)
            sage: G = ImplicitSurface(x + y + z, (x,-1, 1), (y,-1, 1), (z,-1, 1))
            sage: show(G, viewer='jmol')   # indirect doctest
        """
        self.triangulate()
        return IndexFaceSet.jmol_repr(self, render_params)

    def json_repr(self, render_params):
        """
        Return a representation of this object in JavaScript Object Notation (JSON).

        TESTS::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: var('x,y,z')
            (x, y, z)
            sage: G = ImplicitSurface(x + y + z, (x,-1, 1), (y,-1, 1), (z,-1, 1))
            sage: G.json_repr(G.default_render_params())[0].startswith('{"vertices":')
            True
        """
        self.triangulate()
        return IndexFaceSet.json_repr(self, render_params)

    def threejs_repr(self, render_params):
        r"""
        Return a representation of the surface suitable for plotting with three.js.

        EXAMPLES::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: _ = var('x,y,z')
            sage: G = ImplicitSurface(x + y + z, (x,-1, 1), (y,-1, 1), (z,-1, 1))
            sage: G.threejs_repr(G.default_render_params())
            [('surface',
              {'color': '#6666ff',
               'faces': [[0, 1, 2],
                ...
               'opacity': 1.0,
               'vertices': [{'x': ...,
                 'y': -0.9743589743589...,
                 'z': -0.02564102564102...},
                ...
                {'x': -1.0, 'y': 0.9487179487179..., 'z': 0.05128205128205...}]})]
        """
        self.triangulate()
        return IndexFaceSet.threejs_repr(self, render_params)

    def triangulate(self, force=False):
        """
        The IndexFaceSet will be empty until you call this method,
        which generates the faces and vertices according to the
        parameters specified in the constructor for ImplicitSurface.

        Note that if you call this method more than once, subsequent
        invocations will have no effect (this is an optimization to
        avoid repeated work) unless you specify ``force=True`` in the
        keywords.

        EXAMPLES::

            sage: from sage.plot.plot3d.implicit_surface import ImplicitSurface
            sage: var('x,y,z')
            (x, y, z)
            sage: G = ImplicitSurface(x + y + z, (x,-1, 1), (y,-1, 1), (z,-1, 1))
            sage: len(G.vertex_list()), len(G.face_list())
            (0, 0)
            sage: G.triangulate()
            sage: len(G.vertex_list()) > 0, len(G.face_list()) > 0
            (True, True)
            sage: G.show() # This should be fast, since the mesh is already triangulated.
        """
        if self.fcount != 0 and not force:
            # The mesh is already triangulated
            return

        options = dict(xrange=self.xrange, yrange=self.yrange, zrange=self.zrange,
                       region=self.region, smooth=self.smooth,
                       gradient=self.gradient,
                       color_function=self.color_function,
                       colormap=self.colormap,
                       plot_points=self.plot_points)
        cube_marchers = [MarchingCubesTriangles(contour=x, **options) for x in self.contours]
        results = render_implicit(self.f, self.xrange, self.yrange, self.zrange,
                                  self.plot_points, cube_marchers)
        cdef:
            face_c* dest_face
            point_c* dest_vertex
            int fcount = len(results)

        self.realloc(fcount * 3, fcount, fcount * 3)
        for i in range(fcount):
            dest_face = &self._faces[i]
            src_face = results[i]

            dest_face.n = 3
            dest_face.vertices = &self.face_indices[3 * i]

            if self.color_function is not None:
                clist = [(src_face[j]['r'],
                          src_face[j]['g'],
                          src_face[j]['b'])
                         for j in range(3, 6)]
                ct = (sum(clist[j][i] for j in range(3)) / 3 for i in range(3))
                # ct is the mean of the colors of vertices
                dest_face.color.r, dest_face.color.g, dest_face.color.b = ct

            for j in range(3):
                dest_face.vertices[j] = (3 * i) + j
                dest_vertex = &self.vs[(3 * i) + j]
                dest_vertex.x = src_face[j]['x']
                dest_vertex.y = src_face[j]['y']
                dest_vertex.z = src_face[j]['z']


# Data table (courtesy of MarchingCubes.java)
triangle_table2 = ( None,
      ( 0, 8, 3, 7 ),
      ( 0, 1, 9, 7 ), ( 1, 8, 3, 6, 9, 8, 1, 5 ), ( 1, 2, 10, 7 ),
      ( 0, 8, 3, 7, 1, 2, 10, 7 ), ( 9, 2, 10, 6, 0, 2, 9, 5 ),
      ( 2, 8, 3, 6, 2, 10, 8, 1, 10, 9, 8, 3 ), ( 3, 11, 2, 7 ),
      ( 0, 11, 2, 6, 8, 11, 0, 5 ), ( 1, 9, 0, 7, 2, 3, 11, 7 ),
      ( 1, 11, 2, 6, 1, 9, 11, 1, 9, 8, 11, 3 ), ( 3, 10, 1, 6, 11, 10, 3, 5 ),
      ( 0, 10, 1, 6, 0, 8, 10, 1, 8, 11, 10, 3 ),
      ( 3, 9, 0, 6, 3, 11, 9, 1, 11, 10, 9, 3 ), ( 9, 8, 10, 5, 10, 8, 11, 6 ),
      ( 4, 7, 8, 7 ), ( 4, 3, 0, 6, 7, 3, 4, 5 ), ( 0, 1, 9, 7, 8, 4, 7, 7 ),
      ( 4, 1, 9, 6, 4, 7, 1, 1, 7, 3, 1, 3 ), ( 1, 2, 10, 7, 8, 4, 7, 7 ),
      ( 3, 4, 7, 6, 3, 0, 4, 3, 1, 2, 10, 7 ),
      ( 9, 2, 10, 6, 9, 0, 2, 3, 8, 4, 7, 7 ),
      ( 2, 10, 9, 3, 2, 9, 7, 0, 2, 7, 3, 6, 7, 9, 4, 6 ),
      ( 8, 4, 7, 7, 3, 11, 2, 7 ), ( 11, 4, 7, 6, 11, 2, 4, 1, 2, 0, 4, 3 ),
      ( 9, 0, 1, 7, 8, 4, 7, 7, 2, 3, 11, 7 ),
      ( 4, 7, 11, 3, 9, 4, 11, 1, 9, 11, 2, 2, 9, 2, 1, 6 ),
      ( 3, 10, 1, 6, 3, 11, 10, 3, 7, 8, 4, 7 ),
      ( 1, 11, 10, 6, 1, 4, 11, 0, 1, 0, 4, 3, 7, 11, 4, 5 ),
      ( 4, 7, 8, 7, 9, 0, 11, 1, 9, 11, 10, 6, 11, 0, 3, 6 ),
      ( 4, 7, 11, 3, 4, 11, 9, 4, 9, 11, 10, 6 ), ( 9, 5, 4, 7 ),
      ( 9, 5, 4, 7, 0, 8, 3, 7 ), ( 0, 5, 4, 6, 1, 5, 0, 5 ),
      ( 8, 5, 4, 6, 8, 3, 5, 1, 3, 1, 5, 3 ), ( 1, 2, 10, 7, 9, 5, 4, 7 ),
      ( 3, 0, 8, 7, 1, 2, 10, 7, 4, 9, 5, 7 ),
      ( 5, 2, 10, 6, 5, 4, 2, 1, 4, 0, 2, 3 ),
      ( 2, 10, 5, 3, 3, 2, 5, 1, 3, 5, 4, 2, 3, 4, 8, 6 ),
      ( 9, 5, 4, 7, 2, 3, 11, 7 ), ( 0, 11, 2, 6, 0, 8, 11, 3, 4, 9, 5, 7 ),
      ( 0, 5, 4, 6, 0, 1, 5, 3, 2, 3, 11, 7 ),
      ( 2, 1, 5, 3, 2, 5, 8, 0, 2, 8, 11, 6, 4, 8, 5, 5 ),
      ( 10, 3, 11, 6, 10, 1, 3, 3, 9, 5, 4, 7 ),
      ( 4, 9, 5, 7, 0, 8, 1, 5, 8, 10, 1, 2, 8, 11, 10, 3 ),
      ( 5, 4, 0, 3, 5, 0, 11, 0, 5, 11, 10, 6, 11, 0, 3, 6 ),
      ( 5, 4, 8, 3, 5, 8, 10, 4, 10, 8, 11, 6 ), ( 9, 7, 8, 6, 5, 7, 9, 5 ),
      ( 9, 3, 0, 6, 9, 5, 3, 1, 5, 7, 3, 3 ),
      ( 0, 7, 8, 6, 0, 1, 7, 1, 1, 5, 7, 3 ), ( 1, 5, 3, 5, 3, 5, 7, 6 ),
      ( 9, 7, 8, 6, 9, 5, 7, 3, 10, 1, 2, 7 ),
      ( 10, 1, 2, 7, 9, 5, 0, 5, 5, 3, 0, 2, 5, 7, 3, 3 ),
      ( 8, 0, 2, 3, 8, 2, 5, 0, 8, 5, 7, 6, 10, 5, 2, 5 ),
      ( 2, 10, 5, 3, 2, 5, 3, 4, 3, 5, 7, 6 ),
      ( 7, 9, 5, 6, 7, 8, 9, 3, 3, 11, 2, 7 ),
      ( 9, 5, 7, 3, 9, 7, 2, 0, 9, 2, 0, 6, 2, 7, 11, 6 ),
      ( 2, 3, 11, 7, 0, 1, 8, 5, 1, 7, 8, 2, 1, 5, 7, 3 ),
      ( 11, 2, 1, 3, 11, 1, 7, 4, 7, 1, 5, 6 ),
      ( 9, 5, 8, 5, 8, 5, 7, 6, 10, 1, 3, 3, 10, 3, 11, 6 ),
      ( 5, 7, 0, 1, 5, 0, 9, 6, 7, 11, 0, 1, 1, 0, 10, 5, 11, 10, 0, 1 ),
      ( 11, 10, 0, 1, 11, 0, 3, 6, 10, 5, 0, 1, 8, 0, 7, 5, 5, 7, 0, 1 ),
      ( 11, 10, 5, 3, 7, 11, 5, 5 ), ( 10, 6, 5, 7 ),
      ( 0, 8, 3, 7, 5, 10, 6, 7 ), ( 9, 0, 1, 7, 5, 10, 6, 7 ),
      ( 1, 8, 3, 6, 1, 9, 8, 3, 5, 10, 6, 7 ), ( 1, 6, 5, 6, 2, 6, 1, 5 ),
      ( 1, 6, 5, 6, 1, 2, 6, 3, 3, 0, 8, 7 ),
      ( 9, 6, 5, 6, 9, 0, 6, 1, 0, 2, 6, 3 ),
      ( 5, 9, 8, 3, 5, 8, 2, 0, 5, 2, 6, 6, 3, 2, 8, 5 ),
      ( 2, 3, 11, 7, 10, 6, 5, 7 ), ( 11, 0, 8, 6, 11, 2, 0, 3, 10, 6, 5, 7 ),
      ( 0, 1, 9, 7, 2, 3, 11, 7, 5, 10, 6, 7 ),
      ( 5, 10, 6, 7, 1, 9, 2, 5, 9, 11, 2, 2, 9, 8, 11, 3 ),
      ( 6, 3, 11, 6, 6, 5, 3, 1, 5, 1, 3, 3 ),
      ( 0, 8, 11, 3, 0, 11, 5, 0, 0, 5, 1, 6, 5, 11, 6, 6 ),
      ( 3, 11, 6, 3, 0, 3, 6, 1, 0, 6, 5, 2, 0, 5, 9, 6 ),
      ( 6, 5, 9, 3, 6, 9, 11, 4, 11, 9, 8, 6 ), ( 5, 10, 6, 7, 4, 7, 8, 7 ),
      ( 4, 3, 0, 6, 4, 7, 3, 3, 6, 5, 10, 7 ),
      ( 1, 9, 0, 7, 5, 10, 6, 7, 8, 4, 7, 7 ),
      ( 10, 6, 5, 7, 1, 9, 7, 1, 1, 7, 3, 6, 7, 9, 4, 6 ),
      ( 6, 1, 2, 6, 6, 5, 1, 3, 4, 7, 8, 7 ),
      ( 1, 2, 5, 5, 5, 2, 6, 6, 3, 0, 4, 3, 3, 4, 7, 6 ),
      ( 8, 4, 7, 7, 9, 0, 5, 5, 0, 6, 5, 2, 0, 2, 6, 3 ),
      ( 7, 3, 9, 1, 7, 9, 4, 6, 3, 2, 9, 1, 5, 9, 6, 5, 2, 6, 9, 1 ),
      ( 3, 11, 2, 7, 7, 8, 4, 7, 10, 6, 5, 7 ),
      ( 5, 10, 6, 7, 4, 7, 2, 1, 4, 2, 0, 6, 2, 7, 11, 6 ),
      ( 0, 1, 9, 7, 4, 7, 8, 7, 2, 3, 11, 7, 5, 10, 6, 7 ),
      ( 9, 2, 1, 6, 9, 11, 2, 2, 9, 4, 11, 1, 7, 11, 4, 5, 5, 10, 6, 7 ),
      ( 8, 4, 7, 7, 3, 11, 5, 1, 3, 5, 1, 6, 5, 11, 6, 6 ),
      ( 5, 1, 11, 1, 5, 11, 6, 6, 1, 0, 11, 1, 7, 11, 4, 5, 0, 4, 11, 1 ),
      ( 0, 5, 9, 6, 0, 6, 5, 2, 0, 3, 6, 1, 11, 6, 3, 5, 8, 4, 7, 7 ),
      ( 6, 5, 9, 3, 6, 9, 11, 4, 4, 7, 9, 5, 7, 11, 9, 1 ),
      ( 10, 4, 9, 6, 6, 4, 10, 5 ), ( 4, 10, 6, 6, 4, 9, 10, 3, 0, 8, 3, 7 ),
      ( 10, 0, 1, 6, 10, 6, 0, 1, 6, 4, 0, 3 ),
      ( 8, 3, 1, 3, 8, 1, 6, 0, 8, 6, 4, 6, 6, 1, 10, 6 ),
      ( 1, 4, 9, 6, 1, 2, 4, 1, 2, 6, 4, 3 ),
      ( 3, 0, 8, 7, 1, 2, 9, 5, 2, 4, 9, 2, 2, 6, 4, 3 ),
      ( 0, 2, 4, 5, 4, 2, 6, 6 ), ( 8, 3, 2, 3, 8, 2, 4, 4, 4, 2, 6, 6 ),
      ( 10, 4, 9, 6, 10, 6, 4, 3, 11, 2, 3, 7 ),
      ( 0, 8, 2, 5, 2, 8, 11, 6, 4, 9, 10, 3, 4, 10, 6, 6 ),
      ( 3, 11, 2, 7, 0, 1, 6, 1, 0, 6, 4, 6, 6, 1, 10, 6 ),
      ( 6, 4, 1, 1, 6, 1, 10, 6, 4, 8, 1, 1, 2, 1, 11, 5, 8, 11, 1, 1 ),
      ( 9, 6, 4, 6, 9, 3, 6, 0, 9, 1, 3, 3, 11, 6, 3, 5 ),
      ( 8, 11, 1, 1, 8, 1, 0, 6, 11, 6, 1, 1, 9, 1, 4, 5, 6, 4, 1, 1 ),
      ( 3, 11, 6, 3, 3, 6, 0, 4, 0, 6, 4, 6 ), ( 6, 4, 8, 3, 11, 6, 8, 5 ),
      ( 7, 10, 6, 6, 7, 8, 10, 1, 8, 9, 10, 3 ),
      ( 0, 7, 3, 6, 0, 10, 7, 0, 0, 9, 10, 3, 6, 7, 10, 5 ),
      ( 10, 6, 7, 3, 1, 10, 7, 1, 1, 7, 8, 2, 1, 8, 0, 6 ),
      ( 10, 6, 7, 3, 10, 7, 1, 4, 1, 7, 3, 6 ),
      ( 1, 2, 6, 3, 1, 6, 8, 0, 1, 8, 9, 6, 8, 6, 7, 6 ),
      ( 2, 6, 9, 1, 2, 9, 1, 6, 6, 7, 9, 1, 0, 9, 3, 5, 7, 3, 9, 1 ),
      ( 7, 8, 0, 3, 7, 0, 6, 4, 6, 0, 2, 6 ), ( 7, 3, 2, 3, 6, 7, 2, 5 ),
      ( 2, 3, 11, 7, 10, 6, 8, 1, 10, 8, 9, 6, 8, 6, 7, 6 ),
      ( 2, 0, 7, 1, 2, 7, 11, 6, 0, 9, 7, 1, 6, 7, 10, 5, 9, 10, 7, 1 ),
      ( 1, 8, 0, 6, 1, 7, 8, 2, 1, 10, 7, 1, 6, 7, 10, 5, 2, 3, 11, 7 ),
      ( 11, 2, 1, 3, 11, 1, 7, 4, 10, 6, 1, 5, 6, 7, 1, 1 ),
      ( 8, 9, 6, 1, 8, 6, 7, 6, 9, 1, 6, 1, 11, 6, 3, 5, 1, 3, 6, 1 ),
      ( 0, 9, 1, 7, 11, 6, 7, 7 ),
      ( 7, 8, 0, 3, 7, 0, 6, 4, 3, 11, 0, 5, 11, 6, 0, 1 ), ( 7, 11, 6, 7 ),
      ( 7, 6, 11, 7 ), ( 3, 0, 8, 7, 11, 7, 6, 7 ),
      ( 0, 1, 9, 7, 11, 7, 6, 7 ), ( 8, 1, 9, 6, 8, 3, 1, 3, 11, 7, 6, 7 ),
      ( 10, 1, 2, 7, 6, 11, 7, 7 ), ( 1, 2, 10, 7, 3, 0, 8, 7, 6, 11, 7, 7 ),
      ( 2, 9, 0, 6, 2, 10, 9, 3, 6, 11, 7, 7 ),
      ( 6, 11, 7, 7, 2, 10, 3, 5, 10, 8, 3, 2, 10, 9, 8, 3 ),
      ( 7, 2, 3, 6, 6, 2, 7, 5 ), ( 7, 0, 8, 6, 7, 6, 0, 1, 6, 2, 0, 3 ),
      ( 2, 7, 6, 6, 2, 3, 7, 3, 0, 1, 9, 7 ),
      ( 1, 6, 2, 6, 1, 8, 6, 0, 1, 9, 8, 3, 8, 7, 6, 3 ),
      ( 10, 7, 6, 6, 10, 1, 7, 1, 1, 3, 7, 3 ),
      ( 10, 7, 6, 6, 1, 7, 10, 4, 1, 8, 7, 2, 1, 0, 8, 3 ),
      ( 0, 3, 7, 3, 0, 7, 10, 0, 0, 10, 9, 6, 6, 10, 7, 5 ),
      ( 7, 6, 10, 3, 7, 10, 8, 4, 8, 10, 9, 6 ), ( 6, 8, 4, 6, 11, 8, 6, 5 ),
      ( 3, 6, 11, 6, 3, 0, 6, 1, 0, 4, 6, 3 ),
      ( 8, 6, 11, 6, 8, 4, 6, 3, 9, 0, 1, 7 ),
      ( 9, 4, 6, 3, 9, 6, 3, 0, 9, 3, 1, 6, 11, 3, 6, 5 ),
      ( 6, 8, 4, 6, 6, 11, 8, 3, 2, 10, 1, 7 ),
      ( 1, 2, 10, 7, 3, 0, 11, 5, 0, 6, 11, 2, 0, 4, 6, 3 ),
      ( 4, 11, 8, 6, 4, 6, 11, 3, 0, 2, 9, 5, 2, 10, 9, 3 ),
      ( 10, 9, 3, 1, 10, 3, 2, 6, 9, 4, 3, 1, 11, 3, 6, 5, 4, 6, 3, 1 ),
      ( 8, 2, 3, 6, 8, 4, 2, 1, 4, 6, 2, 3 ), ( 0, 4, 2, 5, 4, 6, 2, 3 ),
      ( 1, 9, 0, 7, 2, 3, 4, 1, 2, 4, 6, 6, 4, 3, 8, 6 ),
      ( 1, 9, 4, 3, 1, 4, 2, 4, 2, 4, 6, 6 ),
      ( 8, 1, 3, 6, 8, 6, 1, 0, 8, 4, 6, 3, 6, 10, 1, 3 ),
      ( 10, 1, 0, 3, 10, 0, 6, 4, 6, 0, 4, 6 ),
      ( 4, 6, 3, 1, 4, 3, 8, 6, 6, 10, 3, 1, 0, 3, 9, 5, 10, 9, 3, 1 ),
      ( 10, 9, 4, 3, 6, 10, 4, 5 ), ( 4, 9, 5, 7, 7, 6, 11, 7 ),
      ( 0, 8, 3, 7, 4, 9, 5, 7, 11, 7, 6, 7 ),
      ( 5, 0, 1, 6, 5, 4, 0, 3, 7, 6, 11, 7 ),
      ( 11, 7, 6, 7, 8, 3, 4, 5, 3, 5, 4, 2, 3, 1, 5, 3 ),
      ( 9, 5, 4, 7, 10, 1, 2, 7, 7, 6, 11, 7 ),
      ( 6, 11, 7, 7, 1, 2, 10, 7, 0, 8, 3, 7, 4, 9, 5, 7 ),
      ( 7, 6, 11, 7, 5, 4, 10, 5, 4, 2, 10, 2, 4, 0, 2, 3 ),
      ( 3, 4, 8, 6, 3, 5, 4, 2, 3, 2, 5, 1, 10, 5, 2, 5, 11, 7, 6, 7 ),
      ( 7, 2, 3, 6, 7, 6, 2, 3, 5, 4, 9, 7 ),
      ( 9, 5, 4, 7, 0, 8, 6, 1, 0, 6, 2, 6, 6, 8, 7, 6 ),
      ( 3, 6, 2, 6, 3, 7, 6, 3, 1, 5, 0, 5, 5, 4, 0, 3 ),
      ( 6, 2, 8, 1, 6, 8, 7, 6, 2, 1, 8, 1, 4, 8, 5, 5, 1, 5, 8, 1 ),
      ( 9, 5, 4, 7, 10, 1, 6, 5, 1, 7, 6, 2, 1, 3, 7, 3 ),
      ( 1, 6, 10, 6, 1, 7, 6, 2, 1, 0, 7, 1, 8, 7, 0, 5, 9, 5, 4, 7 ),
      ( 4, 0, 10, 1, 4, 10, 5, 6, 0, 3, 10, 1, 6, 10, 7, 5, 3, 7, 10, 1 ),
      ( 7, 6, 10, 3, 7, 10, 8, 4, 5, 4, 10, 5, 4, 8, 10, 1 ),
      ( 6, 9, 5, 6, 6, 11, 9, 1, 11, 8, 9, 3 ),
      ( 3, 6, 11, 6, 0, 6, 3, 4, 0, 5, 6, 2, 0, 9, 5, 3 ),
      ( 0, 11, 8, 6, 0, 5, 11, 0, 0, 1, 5, 3, 5, 6, 11, 3 ),
      ( 6, 11, 3, 3, 6, 3, 5, 4, 5, 3, 1, 6 ),
      ( 1, 2, 10, 7, 9, 5, 11, 1, 9, 11, 8, 6, 11, 5, 6, 6 ),
      ( 0, 11, 3, 6, 0, 6, 11, 2, 0, 9, 6, 1, 5, 6, 9, 5, 1, 2, 10, 7 ),
      ( 11, 8, 5, 1, 11, 5, 6, 6, 8, 0, 5, 1, 10, 5, 2, 5, 0, 2, 5, 1 ),
      ( 6, 11, 3, 3, 6, 3, 5, 4, 2, 10, 3, 5, 10, 5, 3, 1 ),
      ( 5, 8, 9, 6, 5, 2, 8, 0, 5, 6, 2, 3, 3, 8, 2, 5 ),
      ( 9, 5, 6, 3, 9, 6, 0, 4, 0, 6, 2, 6 ),
      ( 1, 5, 8, 1, 1, 8, 0, 6, 5, 6, 8, 1, 3, 8, 2, 5, 6, 2, 8, 1 ),
      ( 1, 5, 6, 3, 2, 1, 6, 5 ),
      ( 1, 3, 6, 1, 1, 6, 10, 6, 3, 8, 6, 1, 5, 6, 9, 5, 8, 9, 6, 1 ),
      ( 10, 1, 0, 3, 10, 0, 6, 4, 9, 5, 0, 5, 5, 6, 0, 1 ),
      ( 0, 3, 8, 7, 5, 6, 10, 7 ), ( 10, 5, 6, 7 ),
      ( 11, 5, 10, 6, 7, 5, 11, 5 ), ( 11, 5, 10, 6, 11, 7, 5, 3, 8, 3, 0, 7 ),
      ( 5, 11, 7, 6, 5, 10, 11, 3, 1, 9, 0, 7 ),
      ( 10, 7, 5, 6, 10, 11, 7, 3, 9, 8, 1, 5, 8, 3, 1, 3 ),
      ( 11, 1, 2, 6, 11, 7, 1, 1, 7, 5, 1, 3 ),
      ( 0, 8, 3, 7, 1, 2, 7, 1, 1, 7, 5, 6, 7, 2, 11, 6 ),
      ( 9, 7, 5, 6, 9, 2, 7, 0, 9, 0, 2, 3, 2, 11, 7, 3 ),
      ( 7, 5, 2, 1, 7, 2, 11, 6, 5, 9, 2, 1, 3, 2, 8, 5, 9, 8, 2, 1 ),
      ( 2, 5, 10, 6, 2, 3, 5, 1, 3, 7, 5, 3 ),
      ( 8, 2, 0, 6, 8, 5, 2, 0, 8, 7, 5, 3, 10, 2, 5, 5 ),
      ( 9, 0, 1, 7, 5, 10, 3, 1, 5, 3, 7, 6, 3, 10, 2, 6 ),
      ( 9, 8, 2, 1, 9, 2, 1, 6, 8, 7, 2, 1, 10, 2, 5, 5, 7, 5, 2, 1 ),
      ( 1, 3, 5, 5, 3, 7, 5, 3 ), ( 0, 8, 7, 3, 0, 7, 1, 4, 1, 7, 5, 6 ),
      ( 9, 0, 3, 3, 9, 3, 5, 4, 5, 3, 7, 6 ), ( 9, 8, 7, 3, 5, 9, 7, 5 ),
      ( 5, 8, 4, 6, 5, 10, 8, 1, 10, 11, 8, 3 ),
      ( 5, 0, 4, 6, 5, 11, 0, 0, 5, 10, 11, 3, 11, 3, 0, 3 ),
      ( 0, 1, 9, 7, 8, 4, 10, 1, 8, 10, 11, 6, 10, 4, 5, 6 ),
      ( 10, 11, 4, 1, 10, 4, 5, 6, 11, 3, 4, 1, 9, 4, 1, 5, 3, 1, 4, 1 ),
      ( 2, 5, 1, 6, 2, 8, 5, 0, 2, 11, 8, 3, 4, 5, 8, 5 ),
      ( 0, 4, 11, 1, 0, 11, 3, 6, 4, 5, 11, 1, 2, 11, 1, 5, 5, 1, 11, 1 ),
      ( 0, 2, 5, 1, 0, 5, 9, 6, 2, 11, 5, 1, 4, 5, 8, 5, 11, 8, 5, 1 ),
      ( 9, 4, 5, 7, 2, 11, 3, 7 ),
      ( 2, 5, 10, 6, 3, 5, 2, 4, 3, 4, 5, 2, 3, 8, 4, 3 ),
      ( 5, 10, 2, 3, 5, 2, 4, 4, 4, 2, 0, 6 ),
      ( 3, 10, 2, 6, 3, 5, 10, 2, 3, 8, 5, 1, 4, 5, 8, 5, 0, 1, 9, 7 ),
      ( 5, 10, 2, 3, 5, 2, 4, 4, 1, 9, 2, 5, 9, 4, 2, 1 ),
      ( 8, 4, 5, 3, 8, 5, 3, 4, 3, 5, 1, 6 ), ( 0, 4, 5, 3, 1, 0, 5, 5 ),
      ( 8, 4, 5, 3, 8, 5, 3, 4, 9, 0, 5, 5, 0, 3, 5, 1 ), ( 9, 4, 5, 7 ),
      ( 4, 11, 7, 6, 4, 9, 11, 1, 9, 10, 11, 3 ),
      ( 0, 8, 3, 7, 4, 9, 7, 5, 9, 11, 7, 2, 9, 10, 11, 3 ),
      ( 1, 10, 11, 3, 1, 11, 4, 0, 1, 4, 0, 6, 7, 4, 11, 5 ),
      ( 3, 1, 4, 1, 3, 4, 8, 6, 1, 10, 4, 1, 7, 4, 11, 5, 10, 11, 4, 1 ),
      ( 4, 11, 7, 6, 9, 11, 4, 4, 9, 2, 11, 2, 9, 1, 2, 3 ),
      ( 9, 7, 4, 6, 9, 11, 7, 2, 9, 1, 11, 1, 2, 11, 1, 5, 0, 8, 3, 7 ),
      ( 11, 7, 4, 3, 11, 4, 2, 4, 2, 4, 0, 6 ),
      ( 11, 7, 4, 3, 11, 4, 2, 4, 8, 3, 4, 5, 3, 2, 4, 1 ),
      ( 2, 9, 10, 6, 2, 7, 9, 0, 2, 3, 7, 3, 7, 4, 9, 3 ),
      ( 9, 10, 7, 1, 9, 7, 4, 6, 10, 2, 7, 1, 8, 7, 0, 5, 2, 0, 7, 1 ),
      ( 3, 7, 10, 1, 3, 10, 2, 6, 7, 4, 10, 1, 1, 10, 0, 5, 4, 0, 10, 1 ),
      ( 1, 10, 2, 7, 8, 7, 4, 7 ), ( 4, 9, 1, 3, 4, 1, 7, 4, 7, 1, 3, 6 ),
      ( 4, 9, 1, 3, 4, 1, 7, 4, 0, 8, 1, 5, 8, 7, 1, 1 ),
      ( 4, 0, 3, 3, 7, 4, 3, 5 ), ( 4, 8, 7, 7 ),
      ( 9, 10, 8, 5, 10, 11, 8, 3 ), ( 3, 0, 9, 3, 3, 9, 11, 4, 11, 9, 10, 6 ),
      ( 0, 1, 10, 3, 0, 10, 8, 4, 8, 10, 11, 6 ),
      ( 3, 1, 10, 3, 11, 3, 10, 5 ), ( 1, 2, 11, 3, 1, 11, 9, 4, 9, 11, 8, 6 ),
      ( 3, 0, 9, 3, 3, 9, 11, 4, 1, 2, 9, 5, 2, 11, 9, 1 ),
      ( 0, 2, 11, 3, 8, 0, 11, 5 ), ( 3, 2, 11, 7 ),
      ( 2, 3, 8, 3, 2, 8, 10, 4, 10, 8, 9, 6 ), ( 9, 10, 2, 3, 0, 9, 2, 5 ),
      ( 2, 3, 8, 3, 2, 8, 10, 4, 0, 1, 8, 5, 1, 10, 8, 1 ), ( 1, 10, 2, 7 ),
      ( 1, 3, 8, 3, 9, 1, 8, 5 ), ( 0, 9, 1, 7 ), ( 0, 3, 8, 7 ), None )
