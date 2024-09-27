"""
Basic objects such as Sphere, Box, Cone, etc.

AUTHORS:

- Robert Bradshaw 2007-02: initial version
- Robert Bradshaw 2007-08: obj/tachyon rendering, much updating
- Robert Bradshaw 2007-08: cythonization

EXAMPLES::

    sage: from sage.plot.plot3d.shapes import *
    sage: S = Sphere(.5, color='yellow')
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)
    sage: S += Sphere(.1, color='white').translate(.45,-.1,.15)
    sage: S += Sphere(.05, color='black').translate(.51,-.1,.17)
    sage: S += Sphere(.1, color='white').translate(.45, .1,.15)
    sage: S += Sphere(.05, color='black').translate(.51, .1,.17)
    sage: S += Sphere(.1, color='yellow').translate(.5, 0, -.2)
    sage: S.show()
    sage: S.scale(1,1,2).show()

.. PLOT::

    from sage.plot.plot3d.shapes import *
    S = Sphere(.5, color='yellow')
    S += Cone(.5, .5, color='red').translate(0,0,.3)
    S += Sphere(.1, color='white').translate(.45,-.1,.15) + Sphere(.05, color='black').translate(.51,-.1,.17)
    S += Sphere(.1, color='white').translate(.45, .1,.15) + Sphere(.05, color='black').translate(.51, .1,.17)
    S += Sphere(.1, color='yellow').translate(.5, 0, -.2)
    sphinx_plot(S)

::

    sage: from sage.plot.plot3d.shapes import *
    sage: Torus(.7, .2, color=(0,.3,0)).show()

.. PLOT::

    from sage.plot.plot3d.shapes import *
    sphinx_plot(Torus(.7, .2, color=(0,.3,0)))
"""

# ****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.math cimport sqrt, sin, cos, acos, M_PI
from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.misc.decorators import rename_keyword
from sage.plot.plot3d.base import Graphics3dGroup
from sage.plot.plot3d.index_face_set cimport IndexFaceSet, PrimitiveObject
from sage.plot.plot3d.transform cimport point_c


# Helper function to check that Box input is right
def validate_frame_size(size):
    """
    Check that the input is an iterable of length 3 with all
    elements nonnegative and coercible to floats.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import validate_frame_size
        sage: validate_frame_size([3,2,1])
        [3.0, 2.0, 1.0]

    TESTS::

        sage: from sage.plot.plot3d.shapes import validate_frame_size
        sage: validate_frame_size([3,2,-1])
        Traceback (most recent call last):
        ...
        ValueError: each box dimension must be nonnegative
        sage: validate_frame_size([sqrt(-1),3,2])                                       # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: each box dimension must coerce to a float
    """
    if not isinstance(size, (list, tuple)):
        raise TypeError("size must be a list or tuple")
    if len(size) != 3:
        raise TypeError("size must be of length 3")
    try:
        size = [float(x) for x in size]
    except TypeError:
        raise TypeError("each box dimension must coerce to a float")
    for x in size:
        if x < 0:
            raise ValueError("each box dimension must be nonnegative")
    return size


class Box(IndexFaceSet):
    """
    Return a box.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Box

    A square black box::

        sage: show(Box([1,1,1]), color='black')

    .. PLOT::

        from sage.plot.plot3d.shapes import Box
        sphinx_plot(Box([1,1,1], color='black'))

    A red rectangular box::

        sage: show(Box([2,3,4], color='red'))

    .. PLOT::

        from sage.plot.plot3d.shapes import Box
        sphinx_plot(Box([2,3,4], color='red'))

    A stack of boxes::

        sage: show(sum(Box([2,3,1], color='red').translate((0,0,6*i))
        ....:          for i in [0..3]))

    .. PLOT::

        from sage.plot.plot3d.shapes import Box
        P = sum([Box([2,3,1], color='red').translate((0,0,6*i)) for i in range(0,4)])
        sphinx_plot(P)

    A sinusoidal stack of multicolored boxes::

        sage: B = sum(Box([2,4,1/4],                                                    # needs sage.symbolic
        ....:             color=(i/4,i/5,1)).translate((sin(i),0,5-i))
        ....:         for i in [0..20])
        sage: show(B, figsize=6)                                                        # needs sage.symbolic

    .. PLOT::

        from sage.plot.plot3d.shapes import Box
        B = sum([Box([2,4,1/4], color=(i/4.0,i/5.0,1)).translate((sin(i),0,5-i)) for i in range(0,21)])
        sphinx_plot(B)
    """
    def __init__(self, *size, **kwds):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Box
            sage: Box(10, 1, 1) + Box(1, 10, 1) + Box(1, 1, 10)
            Graphics3d Object
        """
        if isinstance(size[0], (tuple, list)):
            size = validate_frame_size(size[0])
        self.size = size
        x, y, z = self.size
        faces = [[(x, y, z), (-x, y, z), (-x,-y, z), ( x,-y, z)],
                 [(x, y, z), ( x, y,-z), (-x, y,-z), (-x, y, z)],
                 [(x, y, z), ( x,-y, z), ( x,-y,-z), ( x, y,-z)] ]
        faces += [list(reversed([(-x,-y,-z) for x,y,z in face])) for face in faces]
        IndexFaceSet.__init__(self, faces, enclosed=True, **kwds)

    def bounding_box(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Box
            sage: Box([1,2,3]).bounding_box()
            ((-1.0, -2.0, -3.0), (1.0, 2.0, 3.0))
        """
        return tuple([-a for a in self.size]), tuple(self.size)

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Box
            sage: Box([1,2,1/4]).x3d_geometry()
            "<Box size='1.0 2.0 0.25'/>"
        """
        return "<Box size='%s %s %s'/>" % tuple(self.size)


@rename_keyword(alpha='opacity')
def ColorCube(size, colors, opacity=1, **kwds):
    """
    Return a cube with given size and sides with given colors.

    INPUT:

    - ``size`` -- 3-tuple of sizes (same as for box and frame)
    - ``colors`` -- list of either 3 or 6 colors
    - ``opacity`` -- (default: 1) opacity of cube sides
    - ``**kwds`` -- passed to the face constructor

    OUTPUT: a 3d graphics object

    EXAMPLES:

    A color cube with translucent sides::

        sage: from sage.plot.plot3d.shapes import ColorCube
        sage: c = ColorCube([1,2,3], ['red', 'blue', 'green', 'black', 'white', 'orange'], opacity=0.5)
        sage: c.show()

    .. PLOT::

        from sage.plot.plot3d.shapes import ColorCube
        c = ColorCube([1,2,3], ['red', 'blue', 'green', 'black', 'white', 'orange'], opacity=0.5)
        sphinx_plot(c)

    ::

        sage: list(c.texture_set())[0].opacity
        0.5

    If you omit the last 3 colors then the first three are repeated (with
    repeated colors on opposing faces)::

        sage: c = ColorCube([0.5,0.5,0.5], ['red', 'blue', 'green'])

    .. PLOT::

        from sage.plot.plot3d.shapes import ColorCube
        c = ColorCube([0.5,0.5,0.5], ['red', 'blue', 'green'])
        sphinx_plot(c)
    """
    if not isinstance(size, (tuple, list)):
        size = (size, size, size)
    box = Box(size)
    faces = box.face_list()
    if len(colors) == 3:
        colors = colors * 2
    all = []
    kwds['opacity'] = opacity

    from sage.plot.plot3d.texture import Texture
    for k in range(6):
        all.append(IndexFaceSet([faces[k]], enclosed=True,
             texture=Texture(colors[k], opacity=opacity),
             **kwds))
    return Graphics3dGroup(all)


cdef class Cone(ParametricSurface):
    """
    A cone, with base in the xy-plane pointing up the z-axis.

    INPUT:

    - ``radius`` -- positive real number

    - ``height`` -- positive real number

    - ``closed`` -- whether or not to include the base (default: ``True``)

    - ``**kwds`` -- passed to the ParametricSurface constructor

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Cone
        sage: c = Cone(3/2, 1, color='red')
        sage: c += Cone(1, 2, color='yellow').translate(3, 0, 0)
        sage: c.show(aspect_ratio=1)

    .. PLOT::

        from sage.plot.plot3d.shapes import Cone
        c = Cone(3/2, 1, color='red') + Cone(1, 2, color='yellow').translate(3, 0, 0)
        sphinx_plot(c)

    We may omit the base::

        sage: Cone(1, 1, closed=False)
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.shapes import Cone
        sphinx_plot(Cone(1, 1, closed=False))

    A spiky plot of the sine function::

        sage: sum(Cone(.1, sin(n), color='yellow').translate(n, sin(n), 0)              # needs sage.symbolic
        ....:     for n in [0..10, step=.1])
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.shapes import Cone
        sphinx_plot(sum(Cone(.1, sin(n/10.0), color='yellow').translate(n/10.0, sin(n/10.0), 0) for n in range(0,100)))

    A Christmas tree::

        sage: T = sum(Cone(exp(-n/5), 4/3*exp(-n/5),                                    # needs sage.symbolic
        ....:              color=(0, .5, 0)).translate(0, 0, -3*exp(-n/5))
        ....:         for n in [1..7])
        sage: T += Cone(1/8, 1, color='brown').translate(0, 0, -3)                      # needs sage.symbolic
        sage: T.show(aspect_ratio=1, frame=False)                                       # needs sage.symbolic

    .. PLOT::

        from sage.plot.plot3d.shapes import Cone
        T = sum(Cone(exp(-n/5.0), 4/3*exp(-n/5.0), color=(0, .5, 0)).translate(0, 0, -3*exp(-n/5.0)) for n in range(8))
        T += Cone(1/8, 1, color='brown').translate(0, 0, -3)
        sphinx_plot(T)
    """
    def __init__(self, radius, height, closed=True, **kwds):
        """
        TESTS::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: c = Cone(1/2, 1, opacity=.5)
        """
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius
        self.height = height
        self.closed = closed

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: Cone(1, 3).x3d_geometry()
            "<Cone bottomRadius='1.0' height='3.0'/>"
        """
        return "<Cone bottomRadius='%s' height='%s'/>" % (self.radius,
                                                          self.height)

    def get_grid(self, ds):
        """
        Return the grid on which to evaluate this parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: Cone(1, 3, closed=True).get_grid(100)
            ([1, 0, -1], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: Cone(1, 3, closed=False).get_grid(100)
            ([1, 0], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: len(Cone(1, 3).get_grid(.001)[1])
            38
        """
        cdef int k, t_res = min(max(int(2*M_PI*self.radius/ds), 5), 37)
        if self.closed:
            urange = [1,0,-1]
        else:
            urange = [1,0]
        vrange = [2*M_PI*k/t_res for k in range(t_res)] + [0.0]
        return urange, vrange

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
        if u == -1:
            res.x, res.y, res.z = 0, 0, 0
        elif u == 0:
            res.x = self.radius*sin(v)
            res.y = self.radius*cos(v)
            res.z = 0
        else: # u == 1:
            res.x, res.y, res.z = 0, 0, self.height


cdef class Cylinder(ParametricSurface):
    """
    A cone, with base in the xy-plane pointing up the z-axis.

    INPUT:

    - ``radius`` -- positive real number

    - ``height`` -- positive real number

    - ``closed`` -- whether or not to include the ends (default: ``True``)

    - ``**kwds`` -- passed to the :class:`ParametricSurface` constructor

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Cylinder
        sage: c = Cylinder(3/2, 1, color='red')
        sage: c += Cylinder(1, 2, color='yellow').translate(3, 0, 0)
        sage: c.show(aspect_ratio=1)

    .. PLOT::

        from sage.plot.plot3d.shapes import Cylinder
        c = Cylinder(3/2, 1, color='red') + Cylinder(1, 2, color='yellow').translate(3, 0, 0)
        sphinx_plot(c)

    We may omit the base::

        sage: Cylinder(1, 1, closed=False)
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.shapes import Cylinder
        sphinx_plot(Cylinder(1, 1, closed=False))

    Some gears::

        sage: # needs sage.symbolic
        sage: G = Cylinder(1, .5) + Cylinder(.25, 3).translate(0, 0, -3)
        sage: G += sum(Cylinder(.2, 1).translate(cos(2*pi*n/9), sin(2*pi*n/9), 0)
        ....:          for n in [1..9])
        sage: G += G.translate(2.3, 0, -.5)
        sage: G += G.translate(3.5, 2, -1)
        sage: G.show(aspect_ratio=1, frame=False)

    .. PLOT::

        from sage.plot.plot3d.shapes import Cylinder
        G = Cylinder(1, .5) + Cylinder(.25, 3).translate(0, 0, -3)
        G += sum(Cylinder(.2, 1).translate(cos(2*pi*n/9.0), sin(2*pi*n/9.0), 0) for n in range(10))
        G += G.translate(2.3, 0, -.5)
        G += G.translate(3.5, 2, -1)
        sphinx_plot(G)
    """
    def __init__(self, radius, height, closed=True, **kwds):
        """
        TESTS::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 1, color='red')
            Graphics3d Object
        """
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius
        self.height = height
        self.closed = closed

    def bounding_box(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 2).bounding_box()
            ((-1.0, -1.0, 0), (1.0, 1.0, 2.0))
        """
        return ((-self.radius, -self.radius, 0),
                (self.radius, self.radius, self.height))

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 2).x3d_geometry()
            "<Cylinder radius='1.0' height='2.0'/>"
        """
        return "<Cylinder radius='%s' height='%s'/>" % (self.radius,
                                                        self.height)

    def tachyon_repr(self, render_params):
        r"""
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: C = Cylinder(1/2, 4, closed=False)
            sage: C.tachyon_repr(C.default_render_params())
            'FCylinder\n   Base 0 0 0\n   Apex 0 0 4.0\n   Rad 0.5\n   texture...     '
            sage: C = Cylinder(1, 2)
            sage: C.tachyon_repr(C.default_render_params())
                ['Ring Center 0 0 0 Normal 0 0 1 Inner 0 Outer 1.0 texture...',
                 'FCylinder\n   Base 0 0 0\n   Apex 0 0 2.0\n   Rad 1.0\n   texture...     ',
                 'Ring Center 0 0 2.0 Normal 0 0 1 Inner 0 Outer 1.0 texture...']
        """
        transform = render_params.transform
        if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])):
            # Tachyon can't do squashed
            return ParametricSurface.tachyon_repr(self, render_params)

        base, top = self.get_endpoints(transform)
        rad = self.get_radius(transform)
        cyl = """FCylinder
   Base %s %s %s
   Apex %s %s %s
   Rad %s
   %s     """%(base[0], base[1], base[2], top[0], top[1], top[2], rad, self.texture.id)
        if self.closed:
            normal = (0,0,1)
            if transform is not None:
                normal = transform.transform_vector(normal)
            base_cap = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % (base[0], base[1], base[2], normal[0], normal[1], normal[2], rad, self.texture.id)
            top_cap  = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % (top[0], top[1], top[2], normal[0], normal[1], normal[2], rad, self.texture.id)
            return [base_cap, cyl, top_cap]
        else:
            return cyl

    def jmol_repr(self, render_params):
        r"""
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder

        For thin cylinders, lines are used::

            sage: C = Cylinder(.1, 4)
            sage: C.jmol_repr(C.default_render_params())
            ['\ndraw line_1 width 0.1 {0 0 0} {0 0 4.0}\ncolor $line_1  [102,102,255]\n']

        For anything larger, we use a pmesh::

            sage: C = Cylinder(3, 1, closed=False)
            sage: C.jmol_repr(C.testing_render_params())
            ['pmesh obj_1 "obj_1.pmesh"\ncolor pmesh  [102,102,255]']
        """
        transform = render_params.transform
        base, top = self.get_endpoints(transform)
        rad = self.get_radius(transform)

        cdef double ratio = sqrt(rad*rad / ((base[0]-top[0])**2 + (base[1]-top[1])**2 + (base[2]-top[2])**2))

        if ratio > .02:
            if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])) or ratio > .05:
                # Jmol can't do squashed
                return ParametricSurface.jmol_repr(self, render_params)

        name = render_params.unique_name('line')
        return ["""
draw %s width %s {%s %s %s} {%s %s %s}\n%s
""" % (name, rad,
       base[0], base[1], base[2],
       top[0], top[1], top[2],
       self.texture.jmol_str("$" + name))]

    def get_endpoints(self, transform=None):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: from sage.plot.plot3d.transform import Transformation
            sage: Cylinder(1, 5).get_endpoints()
            ((0, 0, 0), (0, 0, 5.0))
            sage: Cylinder(1, 5).get_endpoints(Transformation(trans=(1,2,3),
            ....:                                             scale=(2,2,2)))
            ((1.0, 2.0, 3.0), (1.0, 2.0, 13.0))
        """
        if transform is None:
            return (0,0,0), (0,0,self.height)
        else:
            return transform.transform_point((0,0,0)), transform.transform_point((0,0,self.height))

    def get_radius(self, transform=None):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: from sage.plot.plot3d.transform import Transformation
            sage: Cylinder(3, 1).get_radius()
            3.0
            sage: Cylinder(3, 1).get_radius(Transformation(trans=(1,2,3),
            ....:                                          scale=(2,2,2)))
            6.0
        """
        if transform is None:
            return self.radius
        radv = transform.transform_vector((self.radius, 0, 0))
        return sqrt(sum([x * x for x in radv]))

    def get_grid(self, ds):
        """
        Return the grid on which to evaluate this parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: Cylinder(1, 3, closed=True).get_grid(100)
            ([2, 1, -1, -2], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: Cylinder(1, 3, closed=False).get_grid(100)
            ([1, -1], [0.0, 1.2566..., 2.5132..., 3.7699..., 5.0265..., 0.0])
            sage: len(Cylinder(1, 3).get_grid(.001)[1])
            38
        """
        cdef int k, v_res = min(max(int(2*M_PI*self.radius/ds), 5), 37)
        if self.closed:
            urange = [2,1,-1,-2]
        else:
            urange = [1,-1]
        vrange = [2*M_PI*k/v_res for k in range(v_res)] + [0.0]
        return urange, vrange

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
        if u == -2:
            res.x, res.y, res.z = 0, 0, 0
        elif u == -1:
            res.x = self.radius*sin(v)
            res.y = self.radius*cos(v)
            res.z = 0
        elif u == 1:
            res.x = self.radius*sin(v)
            res.y = self.radius*cos(v)
            res.z = self.height
        else: # u == 2:
            res.x, res.y, res.z = 0, 0, self.height


@rename_keyword(alpha='opacity')
def LineSegment(start, end, thickness=1, radius=None, **kwds):
    """
    Create a line segment, which is drawn as a cylinder from start to
    end with radius ``radius``.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import LineSegment, Sphere
        sage: P = (0,0,0.1)
        sage: Q = (0.5,0.6,0.7)
        sage: S = Sphere(.2, color='red').translate(P)
        sage: S += Sphere(.2, color='blue').translate(Q)
        sage: S += LineSegment(P, Q, .05, color='black')
        sage: S.show()

    .. PLOT::

        from sage.plot.plot3d.shapes import LineSegment, Sphere
        P = (0,0,0.1)
        Q = (0.5,0.6,0.7)
        S = Sphere(.2, color='red').translate(P)
        S += Sphere(.2, color='blue').translate(Q)
        S += LineSegment(P, Q, .05, color='black')
        sphinx_plot(S)

    ::

        sage: S = Sphere(.1, color='red').translate(P)
        sage: S += Sphere(.1, color='blue').translate(Q)
        sage: S += LineSegment(P, Q, .15, color='black')
        sage: S.show()

    .. PLOT::

        from sage.plot.plot3d.shapes import LineSegment, Sphere
        P = (0,0,0.1)
        Q = (0.5,0.6,0.7)
        S = Sphere(.1, color='red').translate(P)
        S += Sphere(.1, color='blue').translate(Q)
        S += LineSegment(P, Q, .15, color='black')
        sphinx_plot(S)

    AUTHOR:

    - Robert Bradshaw
    """
    if radius is None:
        radius = thickness/50.0
    start = vector(RDF, start, sparse=False)
    end   = vector(RDF, end, sparse=False)
    zaxis = vector(RDF, (0,0,1), sparse=False)
    diff  = end - start
    height= sqrt(diff.dot_product(diff))
    cyl   = Cylinder(radius, height, **kwds)
    axis  = zaxis.cross_product(diff)
    if axis == 0:
        if diff[2] < 0:
            return cyl.translate(end)
        else:
            return cyl.translate(start)
    else:
        theta = -acos(diff[2]/height)
        return cyl.rotate(axis, theta).translate(start)


def arrow3d(start, end, width=1, radius=None, head_radius=None, head_len=None, **kwds):
    """
    Create a 3d arrow.

    INPUT:

    - ``start`` -- (x,y,z) point; the starting point of the arrow
    - ``end`` -- (x,y,z) point; the end point
    - ``width`` -- (default: 1) how wide the arrow is
    - ``radius`` -- (default: ``width/50.0``) the radius of the arrow
    - ``head_radius`` -- (default: ``3*radius``) radius of arrow head
    - ``head_len`` -- (default: ``3*head_radius``) len of arrow head

    EXAMPLES:

    The default arrow::

        sage: arrow3d((0,0,0), (1,1,1), 1)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(arrow3d((0,0,0), (1,1,1), 1))

    A fat arrow::

        sage: arrow3d((0,0,0), (1,1,1), radius=0.1)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(arrow3d((0,0,0), (1,1,1), radius=0.1))

    A green arrow::

        sage: arrow3d((0,0,0), (1,1,1), color='green')
        Graphics3d Object

    .. PLOT::

        sphinx_plot(arrow3d((0,0,0), (1,1,1), color='green'))

    A fat arrow head::

        sage: arrow3d((2,1,0), (1,1,1), color='green', head_radius=0.3,
        ....:         aspect_ratio=[1,1,1])
        Graphics3d Object

    .. PLOT::

        sphinx_plot(arrow3d((2,1,0), (1,1,1), color='green', head_radius=0.3, aspect_ratio=[1,1,1]))

    Many arrows arranged in a circle (flying spears?)::

        sage: sum(arrow3d((cos(t),sin(t),0), (cos(t),sin(t),1))                         # needs sage.symbolic
        ....:     for t in [0,0.3,..,2*pi])
        Graphics3d Object

    .. PLOT::

        t=0
        G=Graphics()
        while (t<=2*pi):
            G += arrow3d((cos(t),sin(t),0),(cos(t),sin(t),1))
            t +=0.3
        sphinx_plot(G)

    Change the width of the arrow. (Note: for an arrow that scales with zoom, please consider
    the :func:`line3d` function with the option ``arrow_head=True``)::

        sage: arrow3d((0,0,0), (1,1,1), width=1)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(arrow3d((0,0,0), (1,1,1), width=1))

    TESTS:

    If the arrow is too long, the shaft and part of the head is cut off. ::

        sage: a = arrow3d((0,0,0), (0,0,0.5), head_len=1)
        sage: len(a.all)
        1
        sage: type(a.all[0])
        <... 'sage.plot.plot3d.shapes.Cone'>

    Arrows are always constructed pointing up in the z direction from
    the origin, and then rotated/translated into place. This works for
    every arrow direction except the -z direction.  We take care of the
    anomaly by testing to see if the arrow should point in the -z
    direction, and if it should, just scaling the constructed arrow by
    -1 (i.e., every point is sent to its negative). The scaled arrow
    then points downwards. The doctest just tests that the scale of -1
    is applied to the arrow. ::

        sage: a = arrow3d((0,0,0), (0,0,-1))
        sage: a.all[0].get_transformation().transform_point((0,0,1))
        (0.0, 0.0, -1.0)
    """
    if radius is None:
        radius = width/50.0
    if head_radius is None:
        head_radius = 3*radius
    if head_len is None:
        head_len = 3*head_radius
    start = vector(RDF, start, sparse=False)
    end = vector(RDF, end, sparse=False)
    zaxis = vector(RDF, (0,0,1), sparse=False)
    diff = end - start
    length = sqrt(diff.dot_product(diff))
    if length <= head_len:
        arrow = Cone(head_radius*length/head_len, length, **kwds)
    else:
        arrow = Cylinder(radius, length-head_len, **kwds) \
                + Cone(head_radius, head_len, **kwds).translate(0, 0, length-head_len)
    axis = zaxis.cross_product(diff)
    if axis == 0:
        if diff[2] >= 0:
            return arrow.translate(start)
        else:
            return arrow.scale(-1).translate(start)
    else:
        theta = -acos(diff[2]/length)
        return arrow.rotate(axis, theta).translate(start)


cdef class Sphere(ParametricSurface):
    """
    This class represents a sphere centered at the origin.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Sphere
        sage: Sphere(3)
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.shapes import Sphere
        sphinx_plot(Sphere(3))

    Plot with ``aspect_ratio=1`` to see it unsquashed::

        sage: S = Sphere(3, color='blue') + Sphere(2, color='red').translate(0,3,0)
        sage: S.show(aspect_ratio=1)

    .. PLOT::

        from sage.plot.plot3d.shapes import Sphere
        S = Sphere(3, color='blue') + Sphere(2, color='red').translate(0,3,0)
        sphinx_plot(S)

    Scale to get an ellipsoid::

        sage: S = Sphere(1).scale(1,2,1/2)
        sage: S.show(aspect_ratio=1)

    .. PLOT::

        from sage.plot.plot3d.shapes import Sphere
        S = Sphere(1).scale(1,2,1/2)
        sphinx_plot(S)
    """
    def __init__(self, radius, **kwds):
        """
        TESTS::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(3)
            Graphics3d Object
        """
        ParametricSurface.__init__(self, **kwds)
        self.radius = radius

    def bounding_box(self):
        """
        Return the bounding box that contains this sphere.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(3).bounding_box()
            ((-3.0, -3.0, -3.0), (3.0, 3.0, 3.0))
        """
        return ((-self.radius, -self.radius, -self.radius),
                (self.radius, self.radius, self.radius))

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(12).x3d_geometry()
            "<Sphere radius='12.0'/>"
        """
        return "<Sphere radius='%s'/>"%(self.radius)

    def tachyon_repr(self, render_params):
        r"""
        Tachyon can natively handle spheres. Ellipsoids rendering is done
        as a parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: S = Sphere(2)
            sage: S.tachyon_repr(S.default_render_params())
            'Sphere center 0 0 0 Rad 2.0 texture...'
            sage: S.translate(1, 2, 3).scale(3).tachyon_repr(S.default_render_params())
            [['Sphere center 3.0 6.0 9.0 Rad 6.0 texture...']]
            sage: S.scale(1,1/2,1/4).tachyon_repr(S.default_render_params())
            [['TRI V0 0 0 -0.5 V1 0.308116 0.0271646 -0.493844 V2 0.312869 0 -0.493844',
              'texture...',
               ...
              'TRI V0 0.308116 -0.0271646 0.493844 V1 0.312869 0 0.493844 V2 0 0 0.5',
              'texture...']]
        """
        transform = render_params.transform
        if not (transform is None or transform.is_uniform()):
            return ParametricSurface.tachyon_repr(self, render_params)

        if transform is None:
            cen = (0,0,0)
            rad = self.radius
        else:
            cen = transform.transform_point((0,0,0))
            radv = transform.transform_vector((self.radius,0,0))
            rad = sqrt(sum([x*x for x in radv]))
        return "Sphere center %s %s %s Rad %s %s" % (cen[0], cen[1], cen[2], rad, self.texture.id)

    def jmol_repr(self, render_params):
        r"""
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere

        Jmol has native code for handling spheres::

            sage: S = Sphere(2)
            sage: S.jmol_repr(S.default_render_params())
            ['isosurface sphere_1  center {0 0 0} sphere 2.0\ncolor isosurface  [102,102,255]']
            sage: S.translate(10, 100, 1000).jmol_repr(S.default_render_params())
            [['isosurface sphere_1  center {10.0 100.0 1000.0} sphere 2.0\ncolor isosurface  [102,102,255]']]

        It cannot natively handle ellipsoids::

            sage: Sphere(1).scale(2, 3, 4).jmol_repr(S.testing_render_params())
            [['pmesh obj_2 "obj_2.pmesh"\ncolor pmesh  [102,102,255]']]

        Small spheres need extra hints to render well::

            sage: Sphere(.01).jmol_repr(S.default_render_params())
            ['isosurface sphere_1 resolution 100 center {0 0 0} sphere 0.01\ncolor isosurface  [102,102,255]']
        """
        name = render_params.unique_name('sphere')
        transform = render_params.transform
        if not (transform is None or transform.is_uniform()):
            return ParametricSurface.jmol_repr(self, render_params)

        if transform is None:
            cen = (0,0,0)
            rad = self.radius
        else:
            cen = transform.transform_point((0,0,0))
            radv = transform.transform_vector((self.radius,0,0))
            rad = sqrt(sum([x*x for x in radv]))
        if rad < 0.5:
            res = "resolution %s" % min(int(7/rad), 100)
        else:
            res = ""
        return ["isosurface %s %s center {%s %s %s} sphere %s\n%s" % (name, res, cen[0], cen[1], cen[2], rad, self.texture.jmol_str("isosurface"))]

    def get_grid(self, double ds):
        """
        Return the range of variables to be evaluated on to render as a
        parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: Sphere(1).get_grid(100)
            ([-10.0, ..., 10.0],
             [0.0, ..., 3.141592653589793, ..., 0.0])
        """
        cdef int u_res, v_res
        u_res = min(max(int(M_PI*self.radius/ds), 6), 20)
        v_res = min(max(int(2*M_PI * self.radius/ds), 6), 36)
        urange = [-10.0] + [M_PI * k/u_res - M_PI/2 for k in range(1, u_res)] + [10.0]
        vrange = [2*M_PI * k/v_res for k in range(v_res)] + [0.0]
        return urange, vrange

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
        if u == -10:
            res.x, res.y, res.z = 0, 0, -self.radius
        elif u == 10:
            res.x, res.y, res.z = 0, 0, self.radius
        else:
            res.x = self.radius * cos(v) * cos(u)
            res.y = self.radius * sin(v) * cos(u)
            res.z = self.radius * sin(u)


cdef class Torus(ParametricSurface):
    """
    INPUT:

    - ``R`` -- (default: ``1``) outer radius
    - ``r`` -- (default: ``.3``) inner radius

    OUTPUT: a 3d torus

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Torus
        sage: Torus(1, .2).show(aspect_ratio=1)

    .. PLOT::

        from sage.plot.plot3d.shapes import Torus
        sphinx_plot(Torus(1, .2))

    ::

        sage: Torus(1, .7, color='red').show(aspect_ratio=1)

    .. PLOT::

        from sage.plot.plot3d.shapes import Torus
        sphinx_plot(Torus(1, .7, color='red'))

    A rubberband ball::

        sage: show(sum(Torus(1, .03, color=(1, t/30.0, 0)).rotate((1,1,1), t)
        ....:          for t in range(30)))

    .. PLOT::

        from sage.plot.plot3d.shapes import Torus
        sphinx_plot(sum([Torus(1, .03, color=(1, t/30.0, 0)).rotate((1,1,1),t) for t in range(30)]))

    Mmm... doughnuts::

        sage: D = Torus(1, .4, color=(.5, .3, .2))
        sage: D += Torus(1, .3, color='yellow').translate(0, 0, .15)
        sage: G = sum(D.translate(RDF.random_element(-.2, .2),
        ....:                     RDF.random_element(-.2, .2),
        ....:                     .8*t)
        ....:         for t in range(10))
        sage: G.show(aspect_ratio=1, frame=False)

    .. PLOT::

        from sage.plot.plot3d.shapes import Torus
        D = Torus(1, .4, color=(.5, .3, .2)) + Torus(1, .3, color='yellow').translate(0, 0, .15)
        G = sum(D.translate(RDF.random_element(-.2, .2), RDF.random_element(-.2, .2), .8*t) for t in range(10))
        sphinx_plot(G)
    """
    def __init__(self, R=1, r=.3, **kwds):
        """
        TESTS::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: T = Torus(1, .5)
        """
        ParametricSurface.__init__(self, None, **kwds)
        self.R = R
        self.r = r

    def get_grid(self, ds):
        """
        Return the range of variables to be evaluated on to render as a
        parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: Torus(2, 1).get_grid(100)
            ([0.0, -1.047..., -3.141592653589793, ..., 0.0],
             [0.0, 1.047..., 3.141592653589793, ..., 0.0])
        """
        cdef int k, u_divs, v_divs
        u_divs = min(max(int(4*M_PI * self.R/ds), 6), 37)
        v_divs = min(max(int(4*M_PI * self.r/ds), 6), 37)
        urange = [0.0] + [-2*M_PI * k/u_divs for k in range(1, u_divs)] + [0.0]
        vrange = [ 2*M_PI * k/v_divs for k in range(v_divs)] + [0.0]
        return urange, vrange

    cdef int eval_c(self, point_c *res, double u, double v) except -1:
        res.x = (self.R+self.r*sin(v))*sin(u)
        res.y = (self.R+self.r*sin(v))*cos(u)
        res.z = self.r*cos(v)


class Text(PrimitiveObject):
    """
    A text label attached to a point in 3d space. It always starts at the
    origin, translate it to move it elsewhere.

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import Text
        sage: Text("Just a lonely label.")
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.shapes import Text
        sphinx_plot(Text("Just a lonely label. "))

    ::

        sage: pts = [(RealField(10)^3).random_element() for k in range(20)]
        sage: sum(Text(str(P)).translate(P) for P in pts)
        Graphics3d Object

    .. PLOT::

        from sage.plot.plot3d.shapes import Text
        pts = [(RealField(10)**3).random_element() for k in range(20)]
        sphinx_plot(sum(Text(str(P)).translate(P) for P in pts))
    """
    def __init__(self, string, **kwds):
        """
        TESTS::

            sage: from sage.plot.plot3d.shapes import Text
            sage: T = Text("Hi")
        """
        PrimitiveObject.__init__(self, **kwds)
        self.string = string
        self._set_extra_kwds(kwds)

    def x3d_geometry(self):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").x3d_geometry()
            "<Text string='Hi' solid='true'/>"
        """
        return "<Text string='%s' solid='true'/>"%self.string

    def obj_repr(self, render_params):
        """
        The obj file format does not support text strings::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").obj_repr(None)
            ''
        """
        return ''

    def tachyon_repr(self, render_params):
        """
        Strings are not yet supported in Tachyon, so we ignore them for now::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").tachyon_repr(None)
            ''
        """
        return ''
        # Text in Tachyon not implemented yet.
        # I have no idea what the code below is supposed to do.
##         transform = render_params.transform
##         if not (transform is None or transform.is_uniform()):
##             return ParametricSurface.tachyon_repr(self, render_params)
##
##         if transform is None:
##             cen = (0,0,0)
##             rad = self.radius
##         else:
##             cen = transform.transform_point((0,0,0))
##             radv = transform.transform_vector((self.radius,0,0))
##             rad = sqrt(sum([x*x for x in radv]))
##         return "Sphere center %s %s %s Rad %s %s" % (cen[0], cen[1], cen[2], rad, self.texture.id)

    def jmol_repr(self, render_params):
        """
        Labels in jmol must be attached to atoms.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Text
            sage: T = Text("Hi")
            sage: T.jmol_repr(T.testing_render_params())
            ['select atomno = 1', 'color atom  [102,102,255]', 'label "Hi"']
            sage: T = Text("Hi").translate(-1, 0, 0) + Text("Bye").translate(1, 0, 0)
            sage: T.jmol_repr(T.testing_render_params())
            [[['select atomno = 1', 'color atom  [102,102,255]', 'label "Hi"']],
             [['select atomno = 2', 'color atom  [102,102,255]', 'label "Bye"']]]
        """
        cen = (0,0,0)
        if render_params.transform is not None:
            cen = render_params.transform.transform_point(cen)
        render_params.atom_list.append(cen)
        atom_no = len(render_params.atom_list)
        return ['select atomno = %s' % atom_no,
                self.get_texture().jmol_str("atom"),
                'label "%s"' % self.string] #.replace('\n', '|')]

    def threejs_repr(self, render_params):
        r"""
        Return representation of the text suitable for plotting in three.js.

        EXAMPLES::

            sage: T = text3d("Hi", (1, 2, 3), color='red', fontfamily='serif',
            ....:            fontweight='bold', fontstyle='italic', fontsize=20,
            ....:            opacity=0.5)
            sage: T.threejs_repr(T.default_render_params())
            [('text',
              {'color': '#ff0000',
               'fontFamily': ['serif'],
               'fontSize': 20.0,
               'fontStyle': 'italic',
               'fontWeight': 'bold',
               'opacity': 0.5,
               'text': 'Hi',
               'x': 1.0,
               'y': 2.0,
               'z': 3.0})]

        TESTS:

        When created directly via the ``Text`` constructor instead of ``text3d``,
        the text is located at the origin::

            sage: from sage.plot.plot3d.shapes import Text
            sage: T = Text("Hi")
            sage: T.threejs_repr(T.default_render_params())
            [('text',
              {'color': '#6666ff',
               'fontFamily': ['monospace'],
               'fontSize': 14.0,
               'fontStyle': 'normal',
               'fontWeight': 'normal',
               'opacity': 1.0,
               'text': 'Hi',
               'x': 0.0,
               'y': 0.0,
               'z': 0.0})]
        """
        center = (float(0), float(0), float(0))
        if render_params.transform is not None:
            center = render_params.transform.transform_point(center)
        color = '#' + str(self.texture.hex_rgb())
        string = str(self.string)
        text = _validate_threejs_text_style(self._extra_kwds)
        text.update(dict(text=string, x=center[0], y=center[1], z=center[2], color=color))
        return [('text', text)]

    def bounding_box(self):
        """
        Text labels have no extent::

            sage: from sage.plot.plot3d.shapes import Text
            sage: Text("Hi").bounding_box()
            ((0, 0, 0), (0, 0, 0))
        """
        return (0,0,0), (0,0,0)


def _validate_threejs_text_style(style):
    """
    Validate a combination of text styles for use in the Three.js viewer.

    INPUT:

    - ``style`` -- dictionary optionally containing keys: 'color', 'fontSize',
      'fontFamily', 'fontStyle', 'fontWeight', and 'opacity'

    OUTPUT:

    A corrected version of the dict, having printed out warning messages for
    any problems found

    TESTS:

    Default values::

        sage: from sage.plot.plot3d.shapes import _validate_threejs_text_style as validate
        sage: validate(dict())
        {'color': '#000000',
         'fontFamily': ['monospace'],
         'fontSize': 14.0,
         'fontStyle': 'normal',
         'fontWeight': 'normal',
         'opacity': 1.0}

    Color by name or by HTML hex code::

        sage: validate(dict(color='red'))
        {'color': '#ff0000',...}
        sage: validate(dict(color='#ff0000'))
        {'color': '#ff0000',...}
        sage: validate(dict(color='octarine'))
        ...UserWarning: invalid color: octarine...

    Font size in absolute units or in percentage/keyword relative to default::

        sage: validate(dict(fontsize=20.5))
        {...'fontSize': 20.5...}
        sage: validate(dict(fontsize="200%"))
        {...'fontSize': 28...}
        sage: validate(dict(fontsize="x%"))
        ...UserWarning: invalid fontsize: x%...
        sage: validate(dict(fontsize="large"))
        {...'fontSize': 16.8...}
        sage: validate(dict(fontsize="ginormous"))
        ...UserWarning: unknown fontsize: ginormous...

    Font family as list or comma-delimited string::

        sage: validate(dict(fontfamily=[" Times New Roman ", " Georgia", "serif "]))
        {...'fontFamily': ['Times New Roman', 'Georgia', 'serif']...}
        sage: validate(dict(fontfamily=" Times New Roman , Georgia,serif "))
        {...'fontFamily': ['Times New Roman', 'Georgia', 'serif']...}

    Font style keywords including the special syntax for 'oblique' angle::

        sage: validate(dict(fontstyle='italic'))
        {...'fontStyle': 'italic'...}
        sage: validate(dict(fontstyle='oblique 30deg'))
        {...'fontStyle': 'oblique 30deg'...}
        sage: validate(dict(fontstyle='garrish'))
        ...UserWarning: unknown style: garrish...

    Font weight as keyword or integer::

        sage: validate(dict(fontweight='bold'))
        {...'fontWeight': 'bold'...}
        sage: validate(dict(fontweight=500))
        {...'fontWeight': 500...}
        sage: validate(dict(fontweight='roman'))
        {...'fontWeight': 500...}
        sage: validate(dict(fontweight='bold & beautiful'))
        ...UserWarning: unknown fontweight: bold & beautiful...

    Opacity::

        sage: validate(dict(opacity=0.5))
        {...'opacity': 0.5}
    """
    default_color = '#000000' # black
    color = style.get('color', default_color)
    from sage.plot.plot3d.texture import Texture
    try:
        texture = Texture(color=color)
    except ValueError:
        import warnings
        warnings.warn(f"invalid color: {color}, using: {default_color}")
        color = default_color
    else:
        color = '#' + str(texture.hex_rgb())

    default_size = 14.0
    size = style.get('fontsize', default_size)
    try:
        size = float(size)
    except (TypeError, ValueError):
        scale = str(size).lower()
        if scale.endswith('%'):
            try:
                scale = float(scale[:-1]) / 100.0
                size = default_size * scale
            except ValueError:
                import warnings
                warnings.warn(f"invalid fontsize: {size}, using: {default_size}")
                size = default_size
        else:
            from matplotlib.font_manager import font_scalings
            try:
                size = default_size * font_scalings[scale]
            except KeyError:
                import warnings
                warnings.warn(f"unknown fontsize: {size}, using: {default_size}")
                size = default_size

    font = style.get('fontfamily', ['monospace'])
    if isinstance(font, str):
        font = font.split(',')
    font = [str(f).strip() for f in font]

    default_style = 'normal'
    fontstyle = str(style.get('fontstyle', default_style))
    if fontstyle not in ['normal', 'italic'] and not fontstyle.startswith('oblique'): # ex: oblique 30deg
        import warnings
        warnings.warn(f"unknown style: {fontstyle}, using: {default_style}")
        fontstyle = default_style

    default_weight = 'normal'
    weight = style.get('fontweight', default_weight)
    if weight not in ['normal', 'bold']:
        try:
            weight = int(weight)
        except (TypeError, ValueError):
            from matplotlib.font_manager import weight_dict
            try:
                weight = weight_dict[weight]
            except KeyError:
                import warnings
                warnings.warn(f"unknown fontweight: {weight}, using: {default_weight}")
                weight = default_weight

    opacity = float(style.get('opacity', 1.0))

    return dict(color=color, fontSize=size, fontFamily=font, fontStyle=fontstyle,
                fontWeight=weight, opacity=opacity)
