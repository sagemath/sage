# sage.doctest: needs sage.symbolic
r"""
Parametric surface

Graphics 3D object for triangulating surfaces, and a base class for many other
objects that can be represented by a 2D parametrization.

It takes great care to turn degenerate quadrilaterals into triangles and
to propagate identified points to all attached polygons. This is not
so much to save space as it is to assist the raytracers/other rendering
systems to better understand the surface (and especially calculate correct
surface normals).

AUTHORS:

- Robert Bradshaw (2007-08-26): initial version

EXAMPLES::

    sage: from sage.plot.plot3d.parametric_surface import ParametricSurface, MoebiusStrip
    sage: def f(x,y): return x+y, sin(x)*sin(y), x*y
    sage: P = ParametricSurface(f, (srange(0,10,0.1), srange(-5,5.0,0.1)))
    sage: show(P)

.. PLOT::

    from sage.plot.plot3d.parametric_surface import ParametricSurface
    def f(x,y): return x+y, sin(x)*sin(y), x*y
    sphinx_plot(ParametricSurface(f, (srange(0,10,0.1), srange(-5,5.0,0.1))))

::

    sage: S = MoebiusStrip(1,.2)
    sage: S.is_enclosed()
    False
    sage: show(S)

.. PLOT::

    from sage.plot.plot3d.parametric_surface import MoebiusStrip
    sphinx_plot(MoebiusStrip(1,.2))



By default, the surface is colored with one single color. ::

    sage: P = ParametricSurface(f, (srange(0,10,0.1), srange(-5,5.0,0.1)),
    ....:  color="red")
    sage: P.show()

.. PLOT::

    from sage.plot.plot3d.parametric_surface import ParametricSurface
    def f(x,y): return x+y, sin(x)*sin(y), x*y
    sphinx_plot(ParametricSurface(f, (srange(0,10,0.1), srange(-5,5.0,0.1)), color="red"))

One can instead provide a coloring function and a colormap::

    sage: def f(x,y): return x+y, x-y, x*y
    sage: def c(x,y): return sin((x+y)/2)**2
    sage: cm = colormaps.RdYlGn
    sage: P = ParametricSurface(f, (srange(-5,5,0.1), srange(-5,5.0,0.1)), color=(c,cm))
    sage: P.show(viewer='tachyon')

.. PLOT::

    from sage.plot.plot3d.parametric_surface import ParametricSurface
    def f(x,y): return x+y, x-y, x*y
    def c(x,y): return sin((x+y)/2)**2
    cm = colormaps.RdYlGn
    sphinx_plot(ParametricSurface(f, (srange(-5,5,0.1), srange(-5,5.0,0.1)), color=(c,cm)))

Note that the coloring function should rather have values between 0 and 1.
This value is passed to the chosen colormap.

Another colored example::

    sage: colm = colormaps.autumn
    sage: def g(x,y): return x, y, x**2 + y**2
    sage: P = ParametricSurface(g, (srange(-10,10,0.1), srange(-5,5.0,0.1)), color=(c,colm))
    sage: P.show(viewer='tachyon')

.. PLOT::

    from sage.plot.plot3d.parametric_surface import ParametricSurface
    colm = colormaps.autumn
    def g(x,y): return x, y, x**2 + y**2
    def c(x,y): return sin((x+y)/2)**2
    sphinx_plot(ParametricSurface(g, (srange(-10,10,0.1), srange(-5,5.0,0.1)), color=(c,colm)))

.. NOTE::

    One may override ``eval()`` or ``eval_c()`` in a subclass
    rather than passing in a function for greater speed.
    One also would want to override get_grid.

.. TODO::

    actually remove unused points, fix the below code::

        S = ParametricSurface(f=lambda xy: (xy[0],xy[1],0), domain=(range(10),range(10)))
"""
# ****************************************************************************
#      Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_check

from math import cos, sin
from sage.rings.real_double import RDF

from sage.plot.colors import check_color_data
from .base import RenderParams
from .transform cimport point_c, face_c
from sage.ext.interpreters.wrapper_rdf cimport Wrapper_rdf

include "point_c.pxi"


cdef inline bint smash_edge(point_c* vs, face_c* f, int a, int b):
    if point_c_eq(vs[f.vertices[a]], vs[f.vertices[b]]):
        f.vertices[b] = f.vertices[a]
        f.n -= 1
        return 1
    else:
        return 0


cdef class ParametricSurface(IndexFaceSet):
    """
    Base class that initializes the ParametricSurface
    graphics type. This sets options, the function to be plotted, and the
    plotting array as attributes.

    INPUT:

    - ``f`` - (default: ``None``) The defining function. Either a tuple of
      three functions, or a single function which returns a tuple, taking
      two python floats as input. To subclass, pass ``None`` for ``f`` and
      override ``eval_c`` or ``eval`` instead.

    - ``domain`` - (default: ``None``) A tuple of two lists, defining the
      grid of `u,v` values. If ``None``, this will be calculated automatically.

    - ``color`` - (default: ``None``) A pair `(h,c)` where `h` is
      a function with values in `[0,1]` and `c` is a colormap. The
      color of a point `p` is then defined as the composition
      `c(h(p))`

    EXAMPLES::

        sage: from sage.plot.plot3d.parametric_surface import ParametricSurface
        sage: def f(x,y): return cos(x)*sin(y), sin(x)*sin(y), cos(y)+log(tan(y/2))+0.2*x
        sage: S = ParametricSurface(f, (srange(0,12.4,0.1), srange(0.1,2,0.1)))
        sage: show(S)

        sage: len(S.face_list())
        2214

    .. PLOT::

        from sage.plot.plot3d.parametric_surface import ParametricSurface
        def f(x,y): return cos(x)*sin(y), sin(x)*sin(y), cos(y)+log(tan(y/2))+0.2*x
        sphinx_plot(ParametricSurface(f, (srange(0,12.4,0.1), srange(0.1,2,0.1))))

    The Hessenberg surface:

    ::

        sage: def f(u,v):
        ....:     a = 1
        ....:     from math import cos, sin, sinh, cosh
        ....:     x = cos(a)*(cos(u)*sinh(v)-cos(3*u)*sinh(3*v)/3) + sin(a)*(
        ....:         sin(u)*cosh(v)-sin(3*u)*cosh(3*v)/3)
        ....:     y = cos(a)*(sin(u)*sinh(v)+sin(3*u)*sinh(3*v)/3) + sin(a)*(
        ....:         -cos(u)*cosh(v)-cos(3*u)*cosh(3*v)/3)
        ....:     z = cos(a)*cos(2*u)*cosh(2*v)+sin(a)*sin(2*u)*sinh(2*v)
        ....:     return (x,y,z)
        sage: v = srange(float(0),float((3/2)*pi),float(0.1))
        sage: S = ParametricSurface(f, (srange(float(0),float(pi),float(0.1)),
        ....:                srange(float(-1),float(1),float(0.1))), color="blue")
        sage: show(S)

    .. PLOT::

        def f(u,v):
            a = 1
            from math import cos, sin, sinh, cosh
            x = cos(a)*(cos(u)*sinh(v)-cos(3*u)*sinh(3*v)/3) + sin(a)*(sin(u)*cosh(v)-sin(3*u)*cosh(3*v)/3)
            y = cos(a)*(sin(u)*sinh(v)+sin(3*u)*sinh(3*v)/3) + sin(a)*(-cos(u)*cosh(v)-cos(3*u)*cosh(3*v)/3)
            z = cos(a)*cos(2*u)*cosh(2*v)+sin(a)*sin(2*u)*sinh(2*v)
            return (x,y,z)
        from sage.plot.plot3d.parametric_surface import ParametricSurface
        v = srange(float(0),float((3/2)*pi),float(0.1))
        sphinx_plot(ParametricSurface(f, (srange(float(0),float(pi),float(0.1)),srange(float(-1),float(1),float(0.1))), color="blue"))

    A colored example using the ``color`` keyword::

        sage: def g(x,y): return x, y, - x**2 + y**2
        sage: def c(x,y): return sin((x-y/2)*y/4)**2
        sage: cm = colormaps.gist_rainbow
        sage: P = ParametricSurface(g, (srange(-10,10,0.1),
        ....:   srange(-5,5.0,0.1)),color=(c,cm))
        sage: P.show(viewer='tachyon')

    .. PLOT::

        from sage.plot.plot3d.parametric_surface import ParametricSurface
        def g(x,y): return x, y, - x**2 + y**2
        def c(x,y): return sin((x-y/2)*y/4)**2
        cm = colormaps.gist_rainbow
        sphinx_plot(ParametricSurface(g, (srange(-10,10,0.1), srange(-5,5.0,0.1)),color=(c,cm)))
    """

    def __init__(self, f=None, domain=None, **kwds):
        """
        Create the graphics primitive :class:`ParametricSurface`.  See the
        docstring of this class for full documentation.

        EXAMPLES::

            sage: from sage.plot.plot3d.parametric_surface import ParametricSurface
            sage: def f(x,y): return x+y, sin(x)*sin(y), x*y
            sage: S = ParametricSurface(f, (srange(0,12.4,0.1), srange(0.1,2,0.1)))
        """
        if isinstance(f, list):
            f = tuple(f)
        self.f = f
        self.render_grid = domain
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

    def default_render_params(self):
        """
        Return an instance of RenderParams suitable for plotting this object.

        TESTS::

            sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
            sage: type(MoebiusStrip(3,3).default_render_params())
            <class 'sage.plot.plot3d.base.RenderParams'>
        """
        return RenderParams(ds=.075, crease_threshold=.35)

    def x3d_geometry(self):
        r"""
        Return XML-like representation of the coordinates of all points
        in a triangulation of the object along with an indexing of those
        points.

        TESTS::

            sage: _ = var('x,y')
            sage: P = plot3d(x^2-y^2, (x, -2, 2), (y, -2, 2))
            sage: s = P.x3d_str()    # indirect doctest
            sage: s[:100]
            "<Shape>\n<IndexedFaceSet coordIndex='0,1,..."
        """
        self.triangulate(self.default_render_params())
        return IndexFaceSet.x3d_geometry(self)

    def tachyon_repr(self, render_params):
        """
        Return representation of the object suitable for plotting
        using Tachyon ray tracer.

        TESTS::

            sage: _ = var('x,y')
            sage: P = plot3d(x^2-y^2, (x, -2, 2), (y, -2, 2))
            sage: s = P.tachyon_repr(P.default_render_params())
            sage: s[:2]
            ['TRI V0 -2 -2 0 V1 -2 -1.89744 0.399737 V2 -1.89744 -1.89744 0', 'texture...']
        """
        self.triangulate(render_params)
        return IndexFaceSet.tachyon_repr(self, render_params)

    def obj_repr(self, render_params):
        """
        Return a complete representation of object with name, texture, and
        lists of vertices, faces, and back-faces.

        TESTS::

            sage: _ = var('x,y')
            sage: P = plot3d(x^2-y^2, (x, -2, 2), (y, -2, 2))
            sage: s = P.obj_repr(P.default_render_params())
            sage: s[:2]+s[2][:3]+s[3][:3]
            ['g obj_1',
             'usemtl texture...',
             'v -2 -2 0',
             'v -2 -1.89744 0.399737',
             'v -1.89744 -1.89744 0',
             'f 1 2 3 4',
             'f 2 5 6 3',
             'f 5 7 8 6']
        """
        self.triangulate(render_params)
        return IndexFaceSet.obj_repr(self, render_params)

    def jmol_repr(self, render_params):
        r"""
        Return a representation of the object suitable for plotting
        using Jmol.

        TESTS::

            sage: _ = var('x,y')
            sage: P = plot3d(x^2-y^2, (x, -2, 2), (y, -2, 2))
            sage: s = P.jmol_repr(P.testing_render_params())
            sage: s[:10]
            ['pmesh obj_1 "obj_1.pmesh"\ncolor pmesh  [102,102,255]']
        """
        self.triangulate(render_params)
        return IndexFaceSet.jmol_repr(self, render_params)

    def json_repr(self, render_params):
        """
        Return a representation of the object in JSON format as
        a list with one element, which is a string of a dictionary
        listing vertices, faces and colors.

        TESTS::

            sage: _ = var('x,y')
            sage: P = plot3d(x^2-y^2, (x, -2, 2), (y, -2, 2))
            sage: s = P.json_repr(P.default_render_params())
            sage: print(s[0][:100])
            {"vertices":[{"x":-2,"y":-2,"z":0},{"x":-2,"y":-1.89744,"z":0.399737},{"x":-1.89744,"y":-1.89744,"z"

        One test for :trac:`22688`::

            sage: P = spherical_plot3d(sqrt(x-pi/2),(x,0,pi),(y,0,2*pi))
            sage: s = P.json_repr(P.default_render_params())
            sage: 'nan' in s or 'NaN' in s
            False
        """
        self.triangulate(render_params)
        return IndexFaceSet.json_repr(self, render_params)

    def threejs_repr(self, render_params):
        r"""
        Return a representation of the surface suitable for plotting with three.js.

        EXAMPLES::

            sage: _ = var('x,y')
            sage: P = plot3d(x^2-y^2, (x, -2, 2), (y, -2, 2))
            sage: P.threejs_repr(P.default_render_params())
            [('surface',
              {'color': '#6666ff',
               'faces': [[0, 1, 2, 3],
                ...
               'opacity': 1.0,
               'vertices': [{'x': -2.0, 'y': -2.0, 'z': 0.0},
                ...
                {'x': 2.0, 'y': 2.0, 'z': 0.0}]})]

        """
        self.triangulate(render_params)
        return IndexFaceSet.threejs_repr(self, render_params)

    def is_enclosed(self):
        """
        Return a boolean telling whether or not it is necessary to
        render the back sides of the polygons (assuming, of course,
        that they have the correct orientation).

        This is calculated in by verifying the opposite edges
        of the rendered domain either line up or are pinched together.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(1).is_enclosed()
            True

            sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
            sage: MoebiusStrip(1,0.2).is_enclosed()
            False
        """
        if self.fcount == 0:
            self.triangulate()
        return self.enclosed

    def dual(self):
        """
        Return an ``IndexFaceSet`` which is the dual of the
        :class:`ParametricSurface` object as a triangulated surface.

        EXAMPLES:

        As one might expect, this gives an icosahedron::

            sage: D = dodecahedron()
            sage: D.dual()
            Graphics3d Object

        But any enclosed surface should work::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: T =  Torus(1, .2)
            sage: T.dual()
            Graphics3d Object
            sage: T.is_enclosed()
            True

        Surfaces which are not enclosed, though, should raise an exception::

            sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
            sage: M = MoebiusStrip(3,1)
            sage: M.is_enclosed()
            False
            sage: M.dual()
            Traceback (most recent call last):
            ...
            NotImplementedError: this is only implemented for enclosed surfaces
        """
        # This doesn't completely make sense...
        if self.fcount == 0:
            self.triangulate()
        if not self.is_enclosed():
            raise NotImplementedError("this is only implemented for enclosed surfaces")
        return IndexFaceSet.dual(self)

    def bounding_box(self):
        """
        Return the lower and upper corners of a 3D bounding box for ``self``.

        This is used for rendering and ``self`` should fit entirely within this
        box.

        Specifically, the first point returned should have x, y, and z
        coordinates should be the respective infimum over all points in
        ``self``, and the second point is the supremum.

        EXAMPLES::

            sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
            sage: M = MoebiusStrip(7,3,2)
            sage: M.bounding_box()
            ((-10.0, -7.53907349250478..., -2.9940801852848145),
             (10.0, 7.53907349250478..., 2.9940801852848145))
        """
        # We must triangulate before computing the bounding box; otherwise
        # we'll get an empty bounding box, as the bounding box is computed
        # using the triangulation, and before triangulating the triangulation
        # is empty.
        self.triangulate()
        return IndexFaceSet.bounding_box(self)

    def triangulate(self, render_params=None):
        r"""
        Call self.eval_grid() for all `(u,v)` in
        `\text{urange} \times \text{vrange}` to construct this surface.

        The most complicated part of this code is identifying shared
        vertices and shrinking trivial edges. This is not done so much
        to save memory, rather it is needed so normals of the triangles
        can be calculated correctly.

        TESTS::

            sage: from sage.plot.plot3d.parametric_surface import ParametricSurface, MoebiusStrip
            sage: def f(x,y): return x+y, sin(x)*sin(y), x*y                        # indirect doctests
            sage: P = ParametricSurface(f, (srange(0,10,0.1), srange(-5,5.0,0.1)))  # indirect doctests
            sage: P.show()                                                          # indirect doctests
            sage: S = MoebiusStrip(1, .2)                                           # indirect doctests
            sage: S.show()                                                          # indirect doctests
        """
        cdef double u, v
        if render_params is None:
            render_params = self.default_render_params()
        ds = render_params.ds
        if render_params.transform is not None:
            ds /= render_params.transform.max_scale()
        urange, vrange = self.get_grid(ds)
        urange = [float(u) for u in urange]
        vrange = [float(v) for v in vrange]
        if self.render_grid == (urange, vrange) and self.fcount:
            # Already triangulated at on this grid.
            return

        cdef Py_ssize_t i, j
        cdef Py_ssize_t n = len(urange) - 1
        cdef Py_ssize_t m = len(vrange) - 1
        cdef Py_ssize_t ix = 0

        try:
            self.realloc((m+1)*(n+1), m*n, 4*m*n)
            self.eval_grid(urange, vrange)
        except BaseException:
            self.fcount = self.vcount = 0
            self.render_grid = None
            raise

        # face_c.vertices:
        #
        #   0 - 1
        #   |   |
        #   3 - 2

        cdef face_c *face

        for i in range(n):
            for j in range(m):
                sig_check()
                ix = i*m + j
                face = &self._faces[ix]
                face.n = 4
                face.vertices = &self.face_indices[4*ix]
                if self.color_function is not None:
                    face.color.r, face.color.g, face.color.b, _ = self.colormap(self.color_function(urange[i], vrange[j]))

                # Connect to the i-1 row
                if i == 0:
                    if j == 0:
                        face.vertices[0] = 0
                    else:
                        face.vertices[0] = self._faces[ix-1].vertices[1]
                    face.vertices[1] = j+1
                    smash_edge(self.vs, face, 0, 1)
                else:
                    face.vertices[0] = self._faces[ix-m].vertices[3]
                    face.vertices[1] = self._faces[ix-m].vertices[2]

                # Connect to the j-1 col
                if j == 0:
                    face.vertices[3] = (i+1)*(m+1)
                    smash_edge(self.vs, face, 0, 3)
                else:
                    face.vertices[3] = self._faces[ix-1].vertices[2]

                # This is the newly-seen vertex, identify if it's a triangle
                face.vertices[2] = (i+1)*(m+1)+j+1
                smash_edge(self.vs, face, 1, 2)
                smash_edge(self.vs, face, 3, 2)

        # Now we see if it wraps around or is otherwise enclosed
        self.enclosed = True

        cdef face_c *first
        cdef face_c *last
        cdef point_c first_v0
        cdef point_c first_v1
        cdef point_c first_v3
        cdef point_c last_v1
        cdef point_c last_v2
        cdef point_c last_v3
        for j in range(m):
            sig_check()
            first = &self._faces[j]
            first_v0 = self.vs[first.vertices[0]]
            first_v1 = self.vs[first.vertices[1]]
            if not (point_c_isfinite(first_v0) and point_c_isfinite(first_v1)):
                continue
            last = &self._faces[(n-1)*m+j]
            last_v3 = self.vs[last.vertices[3]]
            last_v2 = self.vs[last.vertices[2]]
            if not (point_c_isfinite(last_v3) and point_c_isfinite(last_v2)):
                continue
            if point_c_eq(first_v0, last_v3):
                last.vertices[3] = first.vertices[0]
            elif first.vertices[0] != first.vertices[1] or last.vertices[3] != last.vertices[2]:
                self.enclosed = False
            if point_c_eq(first_v1, last_v2):
                last.vertices[2] = first.vertices[1]
            elif first.vertices[0] != first.vertices[1] or last.vertices[3] != last.vertices[2]:
                self.enclosed = False

        for i in range(n):
            sig_check()
            first = &self._faces[i*m]
            first_v0 = self.vs[first.vertices[0]]
            first_v3 = self.vs[first.vertices[3]]
            if not (point_c_isfinite(first_v0) and point_c_isfinite(first_v3)):
                continue
            last = &self._faces[i*m + m-1]
            last_v1 = self.vs[last.vertices[1]]
            last_v2 = self.vs[last.vertices[2]]
            if not (point_c_isfinite(last_v1) and point_c_isfinite(last_v2)):
                continue
            if point_c_eq(first_v0, last_v1):
                last.vertices[1] = first.vertices[0]
            elif first.vertices[0] != first.vertices[3] or last.vertices[1] != last.vertices[2]:
                self.enclosed = False
            if point_c_eq(first_v3, last_v2):
                last.vertices[2] = first.vertices[3]
            elif first.vertices[0] != first.vertices[3] or last.vertices[1] != last.vertices[2]:
                self.enclosed = False

        # make sure we deleted the correct point from the triangles
        # so that the correct vertices are the first 3 ones
        for ix in range(n * m):
            sig_check()
            face = &self._faces[ix]
            if face.n == 3:
                if face.vertices[3] == face.vertices[2] or face.vertices[3] == face.vertices[0]:
                    pass
                else:
                    if face.vertices[0] == face.vertices[1]:
                        face.vertices[1] = face.vertices[2]
                    face.vertices[2] = face.vertices[3]

        self._clean_point_list()

        self.render_grid = urange, vrange

    def get_grid(self, ds):
        """
        TESTS::

            sage: from sage.plot.plot3d.parametric_surface import ParametricSurface
            sage: def f(x,y): return x+y,x-y,x*y
            sage: P = ParametricSurface(f)
            sage: P.get_grid(.1)
            Traceback (most recent call last):
            ...
            NotImplementedError: you must override the get_grid method
        """
        if self.render_grid is None:
            raise NotImplementedError("you must override the get_grid method")
        return self.render_grid

    cdef int eval_grid(self, urange, vrange) except -1:
        r"""
        This fills in the points ``self.vs`` for all
        `u \in \text{urange}, v \in \text{vrange}`.
        We assume enough memory has been allocated.

        We branch outside the loops for efficiency. The options for self.f are:

        - ``None`` -- call self.eval_c() or self.eval()
                        (One of these is presumably overridden.)
        - tuple -- split into fx, fy, fz and call each separately
        - callable -- call f(u,v)

        In addition, branches are taken for efficient calling of fast callables
        """
        cdef Py_ssize_t i, j
        cdef Py_ssize_t m = len(urange)
        cdef Py_ssize_t n = len(vrange)
        cdef double u, v
        cdef double uv[2]
        cdef point_c *res
        cdef double* ulist = NULL
        cdef double* vlist = NULL
        cdef bint fast_x, fast_y, fast_z

        if self.f is None:
            ulist = to_double_array(urange)
            vlist = to_double_array(vrange)

            res = self.vs
            for i in range(m):
                u = ulist[i]
                for j in range(n):
                    sig_check()
                    v = vlist[j]
                    self.eval_c(res, u, v)
                    res += 1
        elif isinstance(self.f, tuple):
            fx, fy, fz = self.f

            # First, deal with the fast functions (if any)
            fast_x = isinstance(fx, Wrapper_rdf)
            fast_y = isinstance(fy, Wrapper_rdf)
            fast_z = isinstance(fz, Wrapper_rdf)
            if fast_x or fast_y or fast_z:
                ulist = to_double_array(urange)
                vlist = to_double_array(vrange)

                res = self.vs
                if fast_x:  # must be Wrapper_rdf
                    for i in range(m):
                        uv[0] = ulist[i]
                        for j in range(n):
                            sig_check()
                            uv[1] = vlist[j]
                            (<Wrapper_rdf>fx).call_c(uv, &res.x)
                            res += 1

                res = self.vs
                if fast_y:  # must be Wrapper_rdf
                    for i in range(m):
                        uv[0] = ulist[i]
                        for j in range(n):
                            sig_check()
                            uv[1] = vlist[j]
                            (<Wrapper_rdf>fy).call_c(uv, &res.y)
                            res += 1

                res = self.vs
                if fast_z:  # must be Wrapper_rdf
                    for i in range(m):
                        uv[0] = ulist[i]
                        for j in range(n):
                            sig_check()
                            uv[1] = vlist[j]
                            (<Wrapper_rdf>fz).call_c(uv, &res.z)
                            res += 1

            # Finally, deal with the slow functions (if any)
            if (not fast_x) or (not fast_y) or (not fast_z):
                res = self.vs
                for uu in urange:
                    for vv in vrange:
                        sig_check()
                        if not fast_x:
                            res.x = fx(uu, vv)
                        if not fast_y:
                            res.y = fy(uu, vv)
                        if not fast_z:
                            res.z = fz(uu, vv)
                        res += 1
        else:
            res = self.vs
            for uu in urange:
                for vv in vrange:
                    sig_check()
                    res.x, res.y, res.z = self.f(uu, vv)
                    res += 1

        sig_free(ulist)
        sig_free(vlist)

    # One of the following two methods should be overridden in
    # derived classes.

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
        # can't do a cpdef because of the point_c* argument
        res.x, res.y, res.z = self.eval(u, v)

    def eval(self, double u, double v):
        """
        TESTS::

            sage: from sage.plot.plot3d.parametric_surface import ParametricSurface
            sage: def f(x,y): return x+y,x-y,x*y
            sage: P = ParametricSurface(f,(srange(0,1,0.1),srange(0,1,0.1)))
            sage: P.eval(0,0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def plot(self):
        """
        Draw a 3D plot of this graphics object, which just returns this
        object since this is already a 3D graphics object.
        Needed to support PLOT in doctrings, see :trac:`17498`

        EXAMPLES::

            sage: S = parametric_plot3d( (sin, cos, lambda u: u/10), (0, 20))
            sage: S.plot() is S
            True

        """
        return self


class MoebiusStrip(ParametricSurface):
    """
    Base class for the :class:`MoebiusStrip` graphics type. This sets the
    basic parameters of the object.

    INPUT:

    - ``r`` -- a number which can be coerced to a float, serving roughly
      as the radius of the object

    - ``width`` -- a number which can be coerced to a float, which gives the
      width of the object

    - ``twists`` -- (default: 1) an integer, giving the number of twists in the
      object (where one twist is the 'traditional' Möbius strip)

    EXAMPLES::

        sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
        sage: M = MoebiusStrip(3,3)
        sage: M.show()

    .. PLOT::

        from sage.plot.plot3d.parametric_surface import MoebiusStrip
        sphinx_plot(MoebiusStrip(3,3))
    """

    def __init__(self, r, width, twists=1, **kwds):
        """
        Create the graphics primitive MoebiusStrip. See the docstring of
        this class for full documentation.

        EXAMPLES:

        ::

            sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
            sage: M = MoebiusStrip(3,3); M # Same width and radius, roughly
            Graphics3d Object
            sage: N = MoebiusStrip(7,3,2); N # two twists, lots of open area in the middle
            Graphics3d Object
            sage: O = MoebiusStrip(5,1,plot_points=200,color='red'); O # keywords get passed to plot3d
            Graphics3d Object

        """
        ParametricSurface.__init__(self, **kwds)
        self.r = float(r)
        self.width = float(width)
        self.twists = int(twists)

    def get_grid(self, ds):
        """
        Return appropriate `u` and `v` ranges for this MoebiusStrip instance.

        This is intended for internal use in creating an actual plot.

        INPUT:

        -  ``ds`` -- A number, typically coming from a RenderParams object,
           which helps determine the increment for the `v` range for the
           MoebiusStrip object.

        EXAMPLES::

            sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
            sage: N = MoebiusStrip(7,3,2) # two twists
            sage: N.get_grid(N.default_render_params().ds)
            ([-1, 1], [0.0, 0.12566370614359174, 0.25132741228718347, 0.37699111843077515, ...])
        """
        twoPi = RDF.pi() * 2
        # Previous code, which doesn't seem to use any of the parameters
        # TODO: figure out how to use it properly.
        # res = max(min(twoPi*(self.r+self.twists*self.width)/ds, 10), 6*self.twists, 50)
        res = max(6 * self.twists, 50)
        return [-1, 1], [twoPi * k / res for k in range(res + 1)]

    def eval(self, u, v):
        """
        Return a tuple for `x,y,z` coordinates for the given ``u`` and ``v``
        for this MoebiusStrip instance.

        EXAMPLES::

            sage: from sage.plot.plot3d.parametric_surface import MoebiusStrip
            sage: N = MoebiusStrip(7,3,2) # two twists
            sage: N.eval(-1,0)
            (4.0, 0.0, -0.0)
        """
        return ( (self.r + u*self.width*cos(self.twists*v/2)) * cos(v),
                 (self.r + u*self.width*cos(self.twists*v/2)) * sin(v),
                 u*self.width*sin(self.twists*v/2) )


cdef double* to_double_array(py_list) except NULL:
    cdef double* c_list = <double *>sig_malloc(sizeof(double) * len(py_list))
    if c_list == NULL:
        raise MemoryError
    cdef Py_ssize_t i = 0
    cdef double a
    for a in py_list:
        c_list[i] = a
        i += 1
    return c_list
