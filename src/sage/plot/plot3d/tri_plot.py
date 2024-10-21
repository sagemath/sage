# sage_setup: distribution = sagemath-plot
r"""
Adaptive refinement code for 3d surface plotting

AUTHOR:

- Tom Boothby -- Algorithm design, code
- Joshua Kantor -- Algorithm design
- Marshall Hampton -- Docstrings and doctests

.. TODO::

    - Parametrizations (cylindrical, spherical)
    - Massive optimization
"""

###########################################################################
#       Copyright (C) 2007 Tom Boothby <boothby@u.washington.edu>
#                     2007 Joshua Kantor <jkantor@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.plot.colors import hue
from math import sqrt
import random


class Triangle:
    """
    A graphical triangle class.
    """
    def __init__(self,a,b,c,color=0):
        """
        a, b, c : triples (x,y,z) representing corners on a triangle in 3-space.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import Triangle
            sage: tri = Triangle([0,0,0],[-1,2,3],[0,2,0])
            sage: tri._a
            [0, 0, 0]
            sage: tri.__init__([0,0,1],[-1,2,3],[0,2,0])
            sage: tri._a
            [0, 0, 1]
        """
        self._a = a
        self._b = b
        self._c = c
        self._color = color

    def str(self):
        """
        Return a string representation of an instance of the Triangle
        class of the form

            a b c color

        where a, b, and c are corner coordinates and color is the color.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import Triangle
            sage: tri = Triangle([0,0,0],[-1,2,3],[0,2,0])
            sage: print(tri.str())
            [0, 0, 0] [-1, 2, 3] [0, 2, 0] 0
        """
        return f"{self._a} {self._b} {self._c} {self._color}"

    def set_color(self, color):
        """
        This method will reset the color of the triangle.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import Triangle
            sage: tri = Triangle([0,0,0],[-1,2,3],[0,2,1])
            sage: tri._color
            0
            sage: tri.set_color(1)
            sage: tri._color
            1
        """
        self._color = color

    def get_vertices(self):
        """
        Return a tuple of vertex coordinates of the triangle.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import Triangle
            sage: tri = Triangle([0,0,0],[-1,2,3],[0,2,1])
            sage: tri.get_vertices()
            ([0, 0, 0], [-1, 2, 3], [0, 2, 1])
        """
        return (self._a, self._b, self._c)


class SmoothTriangle(Triangle):
    """
    A class for smoothed triangles.
    """
    def __init__(self,a,b,c,da,db,dc,color=0):
        """
        a, b, c : triples (x,y,z) representing corners on a triangle in 3-space
        da, db, dc : triples (dx,dy,dz) representing the normal vector at each point a,b,c

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import SmoothTriangle
            sage: t = SmoothTriangle([1,2,3],[2,3,4],[0,0,0],[0,0,1],[0,1,0],[1,0,0])
            sage: t._a
            [1, 2, 3]
        """
        self._a = a
        self._b = b
        self._c = c
        self._da = da
        self._db = db
        self._dc = dc
        self._color = color

    def str(self):
        """
        Return a string representation of the SmoothTriangle of the form.

            a b c color da db dc

        where a, b, and c are the triangle corner coordinates,
        da, db, dc are normals at each corner, and color is the color.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import SmoothTriangle
            sage: t = SmoothTriangle([1,2,3],[2,3,4],[0,0,0],[0,0,1],[0,1,0],[1,0,0])
            sage: print(t.str())
            [1, 2, 3] [2, 3, 4] [0, 0, 0] 0 [0, 0, 1] [0, 1, 0] [1, 0, 0]
        """
        return "{} {} {} {} {} {} {}".format(self._a, self._b, self._c, self._color, self._da, self._db, self._dc)

    def get_normals(self):
        """
        Return the normals to vertices a, b, and c.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import SmoothTriangle
            sage: t = SmoothTriangle([1,2,3],[2,3,4],[0,0,0],[0,0,1],[0,1,0],[2,0,0])
            sage: t.get_normals()
            ([0, 0, 1], [0, 1, 0], [2, 0, 0])
        """
        return (self._da, self._db, self._dc)


class TriangleFactory:
    def triangle(self, a, b, c, color=None):
        """
        Parameters:
        a, b, c : triples (x,y,z) representing corners on a triangle in 3-space

        Returns:
        a Triangle object with the specified coordinates

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TriangleFactory
            sage: factory = TriangleFactory()
            sage: tri = factory.triangle([0,0,0],[0,0,1],[1,1,0])
            sage: tri.get_vertices()
            ([0, 0, 0], [0, 0, 1], [1, 1, 0])
        """
        if color is None:
            return Triangle(a,b,c)
        else:
            return Triangle(a,b,c,color)

    def smooth_triangle(self, a, b, c, da, db, dc, color=None):
        """
        Parameters:

        - a, b, c : triples (x,y,z) representing corners on a triangle in 3-space
        - da, db, dc : triples (dx,dy,dz) representing the normal vector at each point a,b,c

        Returns:
        a SmoothTriangle object with the specified coordinates and normals

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TriangleFactory
            sage: factory = TriangleFactory()
            sage: sm_tri = factory.smooth_triangle([0,0,0],[0,0,1],[1,1,0],[0,0,1],[0,2,0],[1,0,0])
            sage: sm_tri.get_normals()
            ([0, 0, 1], [0, 2, 0], [1, 0, 0])
        """
        if color is None:
            return SmoothTriangle(a,b,c,da,db,dc)
        else:
            return SmoothTriangle(a,b,c,da,db,dc,color)

    def get_colors(self, list):
        """
        Parameters:
        list: an iterable collection of values which can be cast into colors
        -- typically an RGB triple, or an RGBA 4-tuple

        Returns:
        a list of single parameters which can be passed into the set_color
        method of the Triangle or SmoothTriangle objects generated by this
        factory.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TriangleFactory
            sage: factory = TriangleFactory()
            sage: factory.get_colors([1,2,3])
            [1, 2, 3]
        """
        return list


class TrianglePlot:
    """
    Recursively plot a function of two variables by building squares of 4 triangles, checking at
    every stage whether or not each square should be split into four more squares.  This way,
    more planar areas get fewer triangles, and areas with higher curvature get more triangles.
    """

    def str(self):
        """
        Return a string listing the objects in the instance of the TrianglePlot class.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TrianglePlot, TriangleFactory
            sage: tf = TriangleFactory()
            sage: t = TrianglePlot(tf, lambda x,y: x^3+y*x-1, (-1, 3), (-2, 100), max_depth = 4)
            sage: len(t.str())
            68980
        """
        return "".join(o.str() for o in self._objects)

    def __init__(self, triangle_factory, f, min_x__max_x, min_y__max_y, g=None,
                       min_depth=4, max_depth=8, num_colors=None, max_bend=.3):
        """

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TrianglePlot, TriangleFactory
            sage: tf = TriangleFactory()
            sage: t = TrianglePlot(tf, lambda x,y: x^2+y^2, (0, 1), (0, 1))
            sage: t._f(1,1)
            2
        """
        (min_x, max_x) = min_x__max_x
        (min_y, max_y) = min_y__max_y
        self._triangle_factory = triangle_factory
        self._f = f
        self._g = g
        self._min_depth = min_depth
        self._max_depth = max_depth
        self._max_bend = max_bend
        self._objects = []
        if min(max_x - min_x, max_y - min_y) == 0:
            raise ValueError('plot rectangle is really a line; make sure min_x != max_x and min_y != max_y')
        self._num_colors = num_colors
        if g is None:
            def fcn(x,y):
                return [self._f(x,y)]
        else:
            def fcn(x,y):
                return [self._f(x,y), self._g(x,y)]

        self._fcn = fcn

        # generate the necessary data to kick-start the recursion
        mid_x = (min_x + max_x)/2
        mid_y = (min_y + max_y)/2
        sw_z = fcn(min_x,min_y)
        nw_z = fcn(min_x,max_y)
        se_z = fcn(max_x,min_y)
        ne_z = fcn(max_x,max_y)
        mid_z = fcn(mid_x,mid_y)

        self._min = min(sw_z[0], nw_z[0], se_z[0], ne_z[0], mid_z[0])
        self._max = max(sw_z[0], nw_z[0], se_z[0], ne_z[0], mid_z[0])

        # jump in and start building blocks
        outer = self.plot_block(min_x, mid_x, max_x, min_y, mid_y, max_y, sw_z, nw_z, se_z, ne_z, mid_z, 0)

        # build the boundary triangles
        self.triangulate(outer.left, outer.left_c)
        self.triangulate(outer.top, outer.top_c)
        self.triangulate(outer.right, outer.right_c)
        self.triangulate(outer.bottom, outer.bottom_c)

        zrange = self._max - self._min
        if num_colors is not None and zrange != 0:
            colors = triangle_factory.get_colors([hue(float(i/num_colors)) for i in range(num_colors)])

            for o in self._objects:
                vertices = o.get_vertices()
                avg_z = (vertices[0][2] + vertices[1][2] + vertices[2][2])/3
                o.set_color(colors[int(num_colors * (avg_z - self._min) / zrange)])

    def plot_block(self, min_x, mid_x, max_x, min_y, mid_y, max_y, sw_z, nw_z, se_z, ne_z, mid_z, depth):
        """
        Recursive triangulation function for plotting.

        First six inputs are scalars, next 5 are 2-dimensional lists, and the depth argument
        keeps track of the depth of recursion.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TrianglePlot, TriangleFactory
            sage: tf = TriangleFactory()
            sage: t = TrianglePlot(tf, lambda x,y: x^2 + y^2, (-1,1), (-1, 1), max_depth=3)
            sage: q = t.plot_block(0,.5,1,0,.5,1,[0,1],[0,1],[0,1],[0,1],[0,1],2)
            sage: q.left
            [[(0, 0, 0)], [(0, 0.500000000000000, 0.250000000000000)], [(0, 1, 0)]]
        """

        if depth < self._max_depth:
            # recursion is still an option -- step in one last level if we're within tolerance
            # and just keep going if we're not.
            # assumption: it's cheap to build triangles, so we might as well use all the data
            # we calculate

            # big square boundary midpoints
            mid_w_z = self._fcn(min_x, mid_y)
            mid_n_z = self._fcn(mid_x, max_y)
            mid_e_z = self._fcn(max_x, mid_y)
            mid_s_z = self._fcn(mid_x, min_y)

            next_depth = depth+1
            if depth < self._min_depth:
                # midpoints locations of sub_squares
                qtr1_x = (min_x + mid_x)/2
                qtr1_y = (min_y + mid_y)/2
                qtr3_x = (mid_x + max_x)/2
                qtr3_y = (mid_y + max_y)/2

                sw_depth = next_depth
                nw_depth = next_depth
                se_depth = next_depth
                ne_depth = next_depth
            else:
                #compute the midpoint-to-corner vectors
                sw_v = (min_x - mid_x, min_y - mid_y, sw_z[0] - mid_z[0])
                nw_v = (min_x - mid_x, max_y - mid_y, nw_z[0] - mid_z[0])
                se_v = (max_x - mid_x, min_y - mid_y, se_z[0] - mid_z[0])
                ne_v = (max_x - mid_x, max_y - mid_y, ne_z[0] - mid_z[0])

                #compute triangle normal unit vectors by taking the cross-products
                #of the midpoint-to-corner vectors.  always go around clockwise
                #so we're guaranteed to have a positive value near 1 when neighboring
                #triangles are parallel
                #However -- crossunit doesn't really return a unit vector.  It returns
                #the length of the vector to avoid numerical instability when the
                #length is nearly zero -- rather than divide by nearly zero, we multiply
                #the other side of the inequality by nearly zero -- in general, this
                #should work a bit better because of the density of floating-point
                #numbers near zero.
                norm_w = crossunit(sw_v, nw_v)
                norm_n = crossunit(nw_v, ne_v)
                norm_e = crossunit(ne_v, se_v)
                norm_s = crossunit(se_v, sw_v)

                #compute the dot products of the triangle unit norms
                e_sw = norm_w[0]*norm_s[0] + norm_w[1]*norm_s[1] + norm_w[2]*norm_s[2]
                e_nw = norm_w[0]*norm_n[0] + norm_w[1]*norm_n[1] + norm_w[2]*norm_n[2]
                e_se = norm_e[0]*norm_s[0] + norm_e[1]*norm_s[1] + norm_e[2]*norm_s[2]
                e_ne = norm_e[0]*norm_n[0] + norm_e[1]*norm_n[1] + norm_e[2]*norm_n[2]

                if e_sw < self._max_bend*norm_s[3]*norm_w[3]:
                    sw_depth = next_depth
                else:
                    sw_depth = self._max_depth
                if e_nw < self._max_bend*norm_n[3]*norm_w[3]:
                    nw_depth = next_depth
                else:
                    nw_depth = self._max_depth
                if e_se < self._max_bend*norm_s[3]*norm_e[3]:
                    se_depth = next_depth
                else:
                    se_depth = self._max_depth
                if e_ne < self._max_bend*norm_n[3]*norm_e[3]:
                    ne_depth = next_depth
                else:
                    ne_depth = self._max_depth

                qtr1_x = min_x + (.325 + random.random()/4)*(mid_x-min_x)
                qtr3_x = mid_x + (.325 + random.random()/4)*(max_x-mid_x)
                qtr1_y = min_y + (.325 + random.random()/4)*(mid_y-min_y)
                qtr3_y = mid_y + (.325 + random.random()/4)*(max_y-mid_y)

            # function evaluated at the midpoints (possibly random)
            mid_sw_z = self._fcn(qtr1_x,qtr1_y)
            mid_nw_z = self._fcn(qtr1_x,qtr3_y)
            mid_se_z = self._fcn(qtr3_x,qtr1_y)
            mid_ne_z = self._fcn(qtr3_x,qtr3_y)

            self.extrema([mid_w_z[0], mid_n_z[0], mid_e_z[0], mid_s_z[0], mid_sw_z[0], mid_se_z[0], mid_nw_z[0], mid_sw_z[0]])

            # recurse into the sub-squares
            sw = self.plot_block(min_x, qtr1_x, mid_x, min_y, qtr1_y, mid_y, sw_z, mid_w_z, mid_s_z, mid_z, mid_sw_z, sw_depth)
            nw = self.plot_block(min_x, qtr1_x, mid_x, mid_y, qtr3_y, max_y, mid_w_z, nw_z, mid_z, mid_n_z, mid_nw_z, nw_depth)
            se = self.plot_block(mid_x, qtr3_x, max_x, min_y, qtr1_y, mid_y, mid_s_z, mid_z, se_z, mid_e_z, mid_se_z, se_depth)
            ne = self.plot_block(mid_x, qtr3_x, max_x, mid_y, qtr3_y, max_y, mid_z, mid_n_z, mid_e_z, ne_z, mid_ne_z, ne_depth)

            # join the sub-squares
            self.interface(1, sw.right, sw.right_c, se.left, se.left_c)
            self.interface(1, nw.right, nw.right_c, ne.left, ne.left_c)
            self.interface(0, sw.top, sw.top_c, nw.bottom, nw.bottom_c)
            self.interface(0, se.top, se.top_c, ne.bottom, ne.bottom_c)

            # get the boundary information about the subsquares
            left = sw.left + nw.left[1:]
            left_c = sw.left_c + nw.left_c
            right = se.right + ne.right[1:]
            right_c = se.right_c + ne.right_c
            top = nw.top + ne.top[1:]
            top_c = nw.top_c + ne.top_c
            bottom = sw.bottom + se.bottom[1:]
            bottom_c = sw.bottom_c + se.bottom_c

        else:
            # just build the square we're in
            if self._g is None:
                sw = [(min_x, min_y, sw_z[0])]
                nw = [(min_x, max_y, nw_z[0])]
                se = [(max_x, min_y, se_z[0])]
                ne = [(max_x, max_y, ne_z[0])]
                c = [[(mid_x, mid_y, mid_z[0])]]
            else:
                sw = [(min_x, min_y, sw_z[0]), sw_z[1]]
                nw = [(min_x, max_y, nw_z[0]), nw_z[1]]
                se = [(max_x, min_y, se_z[0]), se_z[1]]
                ne = [(max_x, max_y, ne_z[0]), ne_z[1]]
                c = [[(mid_x, mid_y, mid_z[0]), mid_z[1]]]

            left = [sw, nw]
            left_c = c
            top = [nw, ne]
            top_c = c
            right = [se, ne]
            right_c = c
            bottom = [sw, se]
            bottom_c = c

        return PlotBlock(left, left_c, top, top_c, right, right_c, bottom, bottom_c)

    def interface(self, n, p, p_c, q, q_c):
        """
        Take a pair of lists of points, and compares the (n)th coordinate, and
        "zips" the lists together into one.  The "centers", supplied in p_c and
        q_c are matched up such that the lists describe triangles whose sides
        are "perfectly" aligned.  This algorithm assumes that p and q start and
        end at the same point, and are sorted smallest to largest.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TrianglePlot, TriangleFactory
            sage: tf = TriangleFactory()
            sage: t = TrianglePlot(tf, lambda x,y: x^2 - y*x, (0, -2), (0, 2), max_depth=3)
            sage: t.interface(1, [[(-1/4, 0, 1/16)], [(-1/4, 1/4, 1/8)]], [[(-1/8, 1/8, 1/32)]], [[(-1/4, 0, 1/16)], [(-1/4, 1/4, 1/8)]], [[(-3/8, 1/8, 3/16)]])
            sage: t._objects[-1].get_vertices()
            ((-1/4, 0, 1/16), (-1/4, 1/4, 1/8), (-3/8, 1/8, 3/16))
        """
        m = [p[0]]      # a sorted union of p and q
        mpc = [p_c[0]]  # centers from p_c corresponding to m
        mqc = [q_c[0]]  # centers from q_c corresponding to m

        i = 1
        j = 1

        while i < len(p_c) or j < len(q_c):
            if p[i][0][n] == q[j][0][n]:
                m.append(p[i])
                mpc.append(p_c[i])
                mqc.append(q_c[j])
                i += 1
                j += 1
            elif p[i][0][n] < q[j][0][n]:
                m.append(p[i])
                mpc.append(p_c[i])
                mqc.append(mqc[-1])
                i += 1
            else:
                m.append(q[j])
                mpc.append(mpc[-1])
                mqc.append(q_c[j])
                j += 1

        m.append(p[-1])

        self.triangulate(m, mpc)
        self.triangulate(m, mqc)

    def triangulate(self, p, c):
        """
        Pass in a list of edge points (p) and center points (c).
        Triangles will be rendered between consecutive edge points and the
        center point with the same index number as the earlier edge point.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TrianglePlot, TriangleFactory
            sage: tf = TriangleFactory()
            sage: t = TrianglePlot(tf, lambda x,y: x^2 - y*x, (0, -2), (0, 2))
            sage: t.triangulate([[[1,0,0],[0,0,1]],[[0,1,1],[1,1,1]]],[[[0,3,1]]])
            sage: t._objects[-1].get_vertices()
            ([1, 0, 0], [0, 1, 1], [0, 3, 1])
        """

        if self._g is None:
            for i in range(0,len(p)-1):
                self._objects.append(self._triangle_factory.triangle(p[i][0], p[i+1][0], c[i][0]))
        else:
            for i in range(0,len(p)-1):
                self._objects.append(self._triangle_factory.smooth_triangle(p[i][0], p[i+1][0], c[i][0],p[i][1], p[i+1][1], c[i][1]))

    def extrema(self, list):
        """
        If the num_colors option has been set, this expands the TrianglePlot's _min and _max
        attributes to include the minimum and maximum of the argument list.

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import TrianglePlot, TriangleFactory
            sage: tf = TriangleFactory()
            sage: t = TrianglePlot(tf, lambda x,y: x^2+y^2, (0, 1), (0, 1), num_colors = 3)
            sage: t._min, t._max
            (0, 2)
            sage: t.extrema([-1,2,3,4])
            sage: t._min, t._max
            (-1, 4)
        """
        if self._num_colors is not None:
            self._min = min(list+[self._min])
            self._max = max(list+[self._max])


def crossunit(u,v):
    """
    This function computes triangle normal unit vectors by taking the
    cross-products of the midpoint-to-corner vectors.  It always goes
    around clockwise so we're guaranteed to have a positive value near
    1 when neighboring triangles are parallel.  However -- crossunit
    doesn't really return a unit vector.  It returns the length of the
    vector to avoid numerical instability when the length is nearly zero
    -- rather than divide by nearly zero, we multiply the other side of
    the inequality by nearly zero -- in general, this should work a bit
    better because of the density of floating-point numbers near zero.

    TESTS::

        sage: from sage.plot.plot3d.tri_plot import crossunit
        sage: crossunit([0,-1,0],[0,0,1])
        (-1, 0, 0, 1.0)
    """
    p = (u[1]*v[2] - v[1]*u[2], u[0]*v[2] - v[0]*u[2], u[0]*v[1] - u[1]*v[0])
    l = sqrt(p[0]**2 + p[1]**2 + p[2]**2)
    return (p[0], p[1], p[2], l)


class PlotBlock:
    """
    A container class to hold information about spatial blocks.
    """
    def __init__(self, left, left_c, top, top_c, right, right_c, bottom, bottom_c):
        """

        TESTS::

            sage: from sage.plot.plot3d.tri_plot import PlotBlock
            sage: pb = PlotBlock([0,0,0],[0,1,0],[0,0,1],[0,0,.5],[1,0,0],[1,1,0],[-2,-2,0],[0,0,0])
            sage: pb.left
            [0, 0, 0]
        """
        self.left = left
        self.left_c = left_c
        self.top = top
        self.top_c = top_c
        self.right = right
        self.right_c = right_c
        self.bottom = bottom
        self.bottom_c = bottom_c
