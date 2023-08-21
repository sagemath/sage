r"""
Platonic solids

EXAMPLES: The five platonic solids in a row:

::

    sage: G = tetrahedron((0,-3.5,0), color='blue') + cube((0,-2,0),color=(.25,0,.5))
    sage: G += octahedron(color='red') + dodecahedron((0,2,0), color='orange')
    sage: G += icosahedron(center=(0,4,0), color='yellow')
    sage: G.show(aspect_ratio=[1,1,1])

.. PLOT::

    G = tetrahedron((0,-3.5,0), color='blue') + cube((0,-2,0),color=(.25,0,.5))
    G += octahedron(color='red') + dodecahedron((0,2,0), color='orange')
    G += icosahedron(center=(0,4,0), color='yellow')
    sphinx_plot(G)

All the platonic solids in the same place::

    sage: G = tetrahedron(color='blue',opacity=0.7)
    sage: G += cube(color=(.25,0,.5), opacity=0.7)
    sage: G += octahedron(color='red', opacity=0.7)
    sage: G += dodecahedron(color='orange', opacity=0.7) + icosahedron(opacity=0.7)
    sage: G.show(aspect_ratio=[1,1,1])

.. PLOT::

    G = tetrahedron(color='blue',opacity=0.7)
    G += cube(color=(.25,0,.5), opacity=0.7)
    G += octahedron(color='red', opacity=0.7)
    G += dodecahedron(color='orange', opacity=0.7) + icosahedron(opacity=0.7)
    sphinx_plot(G)

Display nice faces only::

    sage: icosahedron().stickers(['red','blue'], .075, .1)
    Graphics3d Object

.. PLOT::

    sphinx_plot(icosahedron().stickers(['red','blue'], .075, .1))

AUTHORS:

- Robert Bradshaw (2007, 2008): initial version

- William Stein
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

from sage.rings.real_double import RDF
from sage.matrix.constructor import matrix
from sage.misc.decorators import rename_keyword
from .shapes import Box, ColorCube
from .shapes2 import frame3d
from .index_face_set import IndexFaceSet


def index_face_set(face_list, point_list, enclosed, **kwds):
    """
    Helper function that creates ``IndexFaceSet`` object for the
    tetrahedron, dodecahedron, and icosahedron.

    INPUT:

    - ``face_list`` -- list of faces, given explicitly from the
      solid invocation

    - ``point_list`` -- list of points, given explicitly from the
      solid invocation

    - ``enclosed`` -- boolean (default passed is always ``True``
      for these solids)

    TESTS:

    Verify that these are working and passing on keywords::

        sage: tetrahedron(center=(2,0,0),size=2,color='red')
        Graphics3d Object

    ::

        sage: dodecahedron(center=(2,0,0),size=2,color='red')
        Graphics3d Object

    ::

        sage: icosahedron(center=(2,0,0),size=2,color='red')
        Graphics3d Object
    """
    if 'center' in kwds:
        center = kwds['center']
        del kwds['center']
    else:
        center = (0, 0, 0)
    if 'size' in kwds:
        size = kwds['size']
        del kwds['size']
    else:
        size = 1
    I = IndexFaceSet(face_list, point_list, enclosed=enclosed, **kwds)
    return prep(I, center, size, kwds)


def prep(G, center, size, kwds):
    """
    Helper function that scales and translates the platonic
    solid, and passes extra keywords on.

    INPUT:

    - ``center`` -- 3-tuple indicating the center (default passed
      from :func:`index_face_set` is the origin `(0,0,0)`)

    - ``size`` -- number indicating amount to scale by (default
      passed from :func:`index_face_set` is 1)

    - ``kwds`` -- a dictionary of keywords, passed from solid
      invocation by :func:`index_face_set`

    TESTS:

    Verify that scaling and moving the center work together properly,
    and that keywords are passed (see :trac:`10796`)::

        sage: octahedron(center=(2,0,0),size=2,color='red')
        Graphics3d Object
    """
    kwds['threejs_flat_shading'] = True
    G._set_extra_kwds(kwds)
    if size != 1:
        G = G.scale(size)
    if center != (0, 0, 0):
        G = G.translate(center)
    return G

@rename_keyword(alpha='opacity')
def tetrahedron(center=(0, 0, 0), size=1, **kwds):
    """
    A 3d tetrahedron.

    INPUT:

    - ``center`` -- (default: (0,0,0))

    - ``size`` -- (default: 1)

    - ``color`` -- a string (``"red"``, ``"green"``, etc)
      or a tuple (r, g, b) with r, g, b numbers between 0 and 1

    - ``opacity`` -- (default: 1) if less than 1 then is
      transparent

    EXAMPLES: A default colored tetrahedron at the origin::

        sage: tetrahedron()
        Graphics3d Object

    .. PLOT::

        sphinx_plot(tetrahedron())

    A transparent green tetrahedron in front of a solid red one::

        sage: tetrahedron(opacity=0.8, color='green') + tetrahedron((-2,1,0),color='red')
        Graphics3d Object

    .. PLOT::

        sphinx_plot(tetrahedron(opacity=0.8, color='green') + tetrahedron((-2,1,0),color='red'))

    A translucent tetrahedron sharing space with a sphere::

        sage: tetrahedron(color='yellow',opacity=0.7) + sphere(size=.5, color='red')
        Graphics3d Object

    .. PLOT::

        sphinx_plot(tetrahedron(color='yellow',opacity=0.7) + sphere(size = .5, color='red'))

    A big tetrahedron::

        sage: tetrahedron(size=10)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(tetrahedron(size=10))

    A wide tetrahedron::

        sage: tetrahedron(aspect_ratio=[1,1,1]).scale((4,4,1))
        Graphics3d Object

    .. PLOT::

        sphinx_plot(tetrahedron(aspect_ratio=[1,1,1]).scale((4,4,1)))

    A red and blue tetrahedron touching noses::

        sage: tetrahedron(color='red') + tetrahedron((0,0,-2)).scale([1,1,-1])
        Graphics3d Object

    .. PLOT::

        sphinx_plot(tetrahedron(color='red') + tetrahedron((0,0,-2)).scale([1,1,-1]))

    A Dodecahedral complex of 5 tetrahedra (a more elaborate example
    from Peter Jipsen)::

        sage: from math import pi
        sage: v = (sqrt(5.)/2-5/6, 5/6*sqrt(3.)-sqrt(15.)/2, sqrt(5.)/3)
        sage: t = acos(sqrt(5.)/3)/2
        sage: t1 = tetrahedron(aspect_ratio=(1,1,1), opacity=0.5).rotateZ(t)
        sage: t2 = tetrahedron(color='red', opacity=0.5).rotateZ(t).rotate(v,2*pi/5)
        sage: t3 = tetrahedron(color='green', opacity=0.5).rotateZ(t).rotate(v,4*pi/5)
        sage: t4 = tetrahedron(color='yellow', opacity=0.5).rotateZ(t).rotate(v,6*pi/5)
        sage: t5 = tetrahedron(color='orange', opacity=0.5).rotateZ(t).rotate(v,8*pi/5)
        sage: show(t1+t2+t3+t4+t5, frame=False, zoom=1.3)

    .. PLOT::

        v=(sqrt(5.)/2-5/6, 5/6*sqrt(3.)-sqrt(15.)/2, sqrt(5.)/3)
        t=acos(sqrt(5.)/3)/2
        t1=tetrahedron(aspect_ratio=(1,1,1), opacity=0.5).rotateZ(t)
        t2=tetrahedron(color='red', opacity=0.5).rotateZ(t).rotate(v,2*pi/5)
        t3=tetrahedron(color='green', opacity=0.5).rotateZ(t).rotate(v,4*pi/5)
        t4=tetrahedron(color='yellow', opacity=0.5).rotateZ(t).rotate(v,6*pi/5)
        t5=tetrahedron(color='orange', opacity=0.5).rotateZ(t).rotate(v,8*pi/5)
        sphinx_plot(t1+t2+t3+t4+t5)

    AUTHORS:

    - Robert Bradshaw and William Stein
    """
    RR = RDF
    one = RR.one()
    sqrt2 = RR(2).sqrt()
    sqrt6 = RR(6).sqrt()
    point_list = [(0,0,1),
                  (2*sqrt2/3,        0, -one/3),
                  ( -sqrt2/3,  sqrt6/3, -one/3),
                  ( -sqrt2/3, -sqrt6/3, -one/3)]
    face_list = [[0,1,2],[1,3,2],[0,2,3],[0,3,1]]
    if 'aspect_ratio' not in kwds:
        kwds['aspect_ratio'] = [1, 1, 1]
    return index_face_set(face_list, point_list, enclosed=True, center=center, size=size, **kwds)

@rename_keyword(alpha='opacity')
def cube(center=(0, 0, 0), size=1, color=None, frame_thickness=0,
         frame_color=None, **kwds):
    """
    A 3D cube centered at the origin with default side lengths 1.

    INPUT:

    - ``center`` -- (default: (0,0,0))

    - ``size`` -- (default: 1) the side lengths of the
      cube

    - ``color`` -- a string that describes a color; this
      can also be a list of 3-tuples or strings length 6 or 3, in which
      case the faces (and oppositive faces) are colored.

    - ``frame_thickness`` -- (default: 0) if positive,
      then thickness of the frame

    - ``frame_color`` -- (default: None) if given, gives
      the color of the frame

    - ``opacity`` -- (default: 1) if less than 1 then it's
      transparent

    EXAMPLES:

    A simple cube::

        sage: cube()
        Graphics3d Object

    .. PLOT::

        sphinx_plot(cube())

    A red cube::

        sage: cube(color="red")
        Graphics3d Object

    .. PLOT::

        sphinx_plot(cube(color="red"))

    A transparent grey cube that contains a red cube::

        sage: cube(opacity=0.8, color='grey') + cube(size=3/4)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(cube(opacity=0.8, color='grey') + cube(size=3/4))

    A transparent colored cube::

        sage: cube(color=['red', 'green', 'blue'], opacity=0.5)
        Graphics3d Object

    .. PLOT::

        sphinx_plot(cube(color=['red', 'green', 'blue'], opacity=0.5))

    A bunch of random cubes::

        sage: v = [(random(), random(), random()) for _ in [1..30]]
        sage: sum([cube((10*a,10*b,10*c), size=random()/3, color=(a,b,c)) for a,b,c in v])
        Graphics3d Object

    .. PLOT::

        v = [(random(), random(), random()) for _ in range(30)]
        sphinx_plot(sum([cube((10*a,10*b,10*c), size=random()/3, color=(a,b,c)) for a,b,c in v]))

    Non-square cubes (boxes)::

        sage: cube(aspect_ratio=[1,1,1]).scale([1,2,3])
        Graphics3d Object

    .. PLOT::

        sphinx_plot(cube(aspect_ratio=[1,1,1]).scale([1,2,3]))

    ::

        sage: cube(color=['red', 'blue', 'green'],aspect_ratio=[1,1,1]).scale([1,2,3])
        Graphics3d Object

    .. PLOT::

        sphinx_plot(cube(color=['red', 'blue', 'green'],aspect_ratio=[1,1,1]).scale([1,2,3]))

    And one that is colored::

        sage: cube(color=['red', 'blue', 'green', 'black', 'white', 'orange'],
        ....:      aspect_ratio=[1,1,1]).scale([1,2,3])
        Graphics3d Object

    .. PLOT::

        sphinx_plot(cube(color=['red', 'blue', 'green', 'black', 'white', 'orange'],aspect_ratio=[1,1,1]).scale([1,2,3]))

    A nice translucent color cube with a frame::

        sage: c = cube(color=['red', 'blue', 'green'], frame=False, frame_thickness=2,
        ....:          frame_color='brown', opacity=0.8)
        sage: c
        Graphics3d Object

    .. PLOT::

        c = cube(color=['red', 'blue', 'green'], frame=False, frame_thickness=2,frame_color='brown', opacity=0.8)
        sphinx_plot(c)

    A raytraced color cube with frame and transparency::

        sage: c.show(viewer='tachyon')

    This shows :trac:`11272` has been fixed::

        sage: cube(center=(10, 10, 10), size=0.5).bounding_box()
        ((9.75, 9.75, 9.75), (10.25, 10.25, 10.25))

    AUTHORS:

    - William Stein
    """
    if isinstance(color, (list, tuple)) and len(color) > 0 and isinstance(color[0], (list, tuple, str)):
        B = ColorCube(size=[0.5,0.5,0.5], colors=color, **kwds)
    else:
        if color is not None:
            kwds['color'] = color
        B = Box(0.5, 0.5, 0.5, **kwds)
    if frame_thickness > 0:
        if frame_color is None:
            B += frame3d((-0.5,-0.5,-0.5),(0.5,0.5,0.5), thickness=frame_thickness)
        else:
            B += frame3d((-0.5,-0.5,-0.5),(0.5,0.5,0.5), thickness=frame_thickness, color=frame_color)
    return prep(B, center, size, kwds)

@rename_keyword(alpha='opacity')
def octahedron(center=(0, 0, 0), size=1, **kwds):
    r"""
    Return an octahedron.

    INPUT:

    - ``center`` -- (default: (0,0,0))

    - ``size`` -- (default: 1)

    - ``color`` -- a string that describes a color; this
      can also be a list of 3-tuples or strings length 6 or 3, in which
      case the faces (and oppositive faces) are colored.

    - ``opacity`` -- (default: 1) if less than 1 then is
      transparent

    EXAMPLES::

        sage: G = octahedron((1,4,3), color='orange')
        sage: G += octahedron((0,2,1), size=2, opacity=0.6)
        sage: G
        Graphics3d Object

    .. PLOT::

        G = octahedron((1,4,3), color='orange')
        G += octahedron((0,2,1), size=2, opacity=0.6)
        sphinx_plot(G)

    """
    kwds['enclosed'] = True
    if 'aspect_ratio' not in kwds:
        kwds['aspect_ratio'] = [1, 1, 1]
    return prep(Box(1,1,1).dual(**kwds), center, size, kwds)

@rename_keyword(alpha='opacity')
def dodecahedron(center=(0, 0, 0), size=1, **kwds):
    r"""
    A dodecahedron.

    INPUT:

    - ``center`` -- (default: (0,0,0))

    - ``size`` -- (default: 1)

    - ``color`` -- a string that describes a color; this
      can also be a list of 3-tuples or strings length 6 or 3, in which
      case the faces (and oppositive faces) are colored.

    - ``opacity`` -- (default: 1) if less than 1 then is transparent

    EXAMPLES: A plain Dodecahedron::

        sage: dodecahedron()
        Graphics3d Object

    .. PLOT::

        sphinx_plot(dodecahedron())

    A translucent dodecahedron that contains a black sphere::

        sage: G = dodecahedron(color='orange', opacity=0.8)
        sage: G += sphere(size=0.5, color='black')
        sage: G
        Graphics3d Object

    .. PLOT::

        G = dodecahedron(color='orange', opacity=0.8)
        G += sphere(size=0.5, color='black')
        sphinx_plot(G)

    CONSTRUCTION: This is how we construct a dodecahedron. We let one
    point be `Q = (0,1,0)`.

    Now there are three points spaced equally on a circle around the
    north pole. The other requirement is that the angle between them be
    the angle of a pentagon, namely `3\pi/5`. This is enough to
    determine them. Placing one on the `xz`-plane we have.

    `P_1 = \left(t, 0, \sqrt{1-t^2}\right)`

    `P_2 = \left(-\frac{1}{2}t, \frac{\sqrt{3}}{2}t, \sqrt{1-t^2}\right)`

    `P_3 = \left(-\frac{1}{2}t, \frac{\sqrt{3}}{2}t, \sqrt{1-t^2}\right)`

    Solving
    `\frac{(P_1-Q) \cdot (P_2-Q)}{|P_1-Q||P_2-Q|} = \cos(3\pi/5)`
    we get `t = 2/3`.

    Now we have 6 points `R_1, ..., R_6` to close the three
    top pentagons. These can be found by mirroring `P_2` and
    `P_3` by the `yz`-plane and rotating around the
    `y`-axis by the angle `\theta` from `Q` to
    `P_1`. Note that `\cos(\theta) = t = 2/3` and so
    `\sin(\theta) = \sqrt{5}/3`. Rotation gives us the other
    four.

    Now we reflect through the origin for the bottom half.

    AUTHORS:

    - Robert Bradshaw, William Stein
    """
    RR = RDF
    one = RR(1)
    sqrt3 = RR(3).sqrt()
    sqrt5 = RR(5).sqrt()
    R3 = RR**3
    rot = matrix(RR, [[  -one/2,-sqrt3/2, 0],
                      [ sqrt3/2,  -one/2, 0],
                      [       0,       0, 1]])
    rot2 = rot*rot

    # The top
    Q = R3([0,0,1])
    # The first ring
    P1 = R3([2*one/3, 0, sqrt5/3])
    # The second ring
    R1 = R3([sqrt5/3,  1/sqrt3, one/3])
    R2 = R3([sqrt5/3, -1/sqrt3, one/3])

    top = [Q, P1, rot*P1, rot2*P1, R1, rot*R2, rot*R1, rot2*R2, rot2*R1, R2]
    point_list = top + [-p for p in reversed(top)]

    top_faces = [[0,1,4,5,2],
                 [0,2,6,7,3],
                 [0,3,8,9,1],
                 [1,9,13,12,4],
                 [2,5,11,10,6],
                 [3,7,15,14,8]]
    face_list = top_faces + [[19-p for p in reversed(f)] for f in top_faces]

    if 'aspect_ratio' not in kwds:
        kwds['aspect_ratio'] = [1,1,1]
    return index_face_set(face_list, point_list, enclosed=True, center=center, size=size, **kwds)

#    if style == 'vertices' or style == 'edges':
#        from sage.plot.colors import rainbow
#        colors = rainbow(len(vs), 'rgbtuple')
#        #vertex_spheres = [Box(.05, .05, .05, color=color).translate(p) for p in vs]
#        vertex_spheres = [Box(.05, .05, .05, color=c).translate(p) for p,c in zip(vs,colors)]
#        faces = IndexFaceSet([[tuple(vs[i]) for i in f] for f in face_list])
#        vertex_spheres += [faces.stickers(['red','yellow','blue','purple','black','orange'], .1, .1)] # [faces]
#        return Graphics3dGroup(vertex_spheres)

@rename_keyword(alpha='opacity')
def icosahedron(center=(0, 0, 0), size=1, **kwds):
    r"""
    An icosahedron.

    INPUT:

    - ``center`` -- (default: (0, 0, 0))

    - ``size`` -- (default: 1)

    - ``color`` -- a string that describes a color; this
      can also be a list of 3-tuples or strings length 6 or 3, in which
      case the faces (and oppositive faces) are colored.

    - ``opacity`` -- (default: 1) if less than 1 then is transparent

    EXAMPLES::

        sage: icosahedron()
        Graphics3d Object

    .. PLOT::

        sphinx_plot(icosahedron())

    Two icosahedra at different positions of different sizes. ::

        sage: p = icosahedron((-1/2,0,1), color='orange')
        sage: p += icosahedron((2,0,1), size=1/2, color='red', aspect_ratio=[1,1,1])
        sage: p
        Graphics3d Object

    .. PLOT::

        p = icosahedron((-1/2,0,1), color='orange')
        p += icosahedron((2,0,1), size = 0.5, color='red', aspect_ratio=[1,1,1])
        sphinx_plot(p)

    """
    kwds['enclosed'] = True
    if 'aspect_ratio' not in kwds:
        kwds['aspect_ratio'] = [1, 1, 1]
    return prep(dodecahedron().dual(**kwds), center, size, kwds)
