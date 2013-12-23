r"""
Base classes for 3D Graphics objects and plotting

AUTHORS:

- Robert Bradshaw (2007-02): initial version

- Robert Bradshaw (2007-08): Cythonization, much optimization

- William Stein (2008)

TODO: - finish integrating tachyon - good default lights, camera
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from cpython.list cimport *

import os
from math import atan2
from random import randint
import zipfile
from cStringIO import StringIO

import sage.misc.misc

from sage.modules.free_module_element import vector

from sage.rings.real_double import RDF
from sage.misc.functional import sqrt, atan, acos
from sage.misc.temporary_file import tmp_filename

from texture import Texture, is_Texture
from transform cimport Transformation, point_c, face_c
include "point_c.pxi"

from sage.interfaces.tachyon import tachyon_rt

# import the double infinity constant
cdef extern from "math.h":
     enum: INFINITY


default_texture = Texture()
pi = RDF.pi()

cdef class Graphics3d(SageObject):
    """
    This is the baseclass for all 3d graphics objects.
    """
    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: S = sphere((0, 0, 0), 1)
            sage: print S
            Graphics3d Object
        """
        return str(self)

    def _graphics_(self):
        """
        Show graphics.

        The presence of this method is used by the displayhook to
        decide that we want to see a graphical output by default.

        OUTPUT:

        Return ``True`` if graphical output was generated (might not
        be shown in doctest mode), otherwise ``False``.

        EXAMPLES::

            sage: S = sphere((0, 0, 0), 1)
            sage: S._graphics_()
            True
            sage: S  # also productes graphics
            sage: [S, S]
            [Graphics3d Object, Graphics3d Object]
        """
        self.show()
        return True

    def __str__(self):
        """
        EXAMPLES::

            sage: S = sphere((0, 0, 0), 1)
            sage: str(S)
            'Graphics3d Object'
        """
        return "Graphics3d Object"

    def __add__(left, right):
        """
        Addition of objects adds them to the same scene.

        EXAMPLES::
            sage: A = sphere((0,0,0), 1, color='red')
            sage: B = dodecahedron((2, 0, 0), color='yellow')
            sage: A+B

        For convenience, we take 0 and None to be the additive identity::

            sage: A + 0 is A
            True
            sage: A + None is A, 0 + A is A, None + A is A
            (True, True, True)

        In particular, this allows us to use the sum() function without
        having to provide an empty starting object::

            sage: sum(point3d((cos(n), sin(n), n)) for n in [0..10, step=.1])

        A Graphics 3d object can also be added a 2d graphic object::

            sage: A = sphere((0, 0, 0), 1) + circle((0, 0), 1.5)
            sage: A.show(aspect_ratio=1)
        """
        if right == 0 or right is None:
            return left
        elif left == 0 or left is None:
            return right
        elif not isinstance(left, Graphics3d):
            left = left.plot3d()
        elif not isinstance(right, Graphics3d):
            right = right.plot3d()
        return Graphics3dGroup([left, right])

    def _set_extra_kwds(self, kwds):
        """
        Allows one to pass rendering arguments on as if they were set in the constructor.

        EXAMPLES::

            sage: S = sphere((0, 0, 0), 1)
            sage: S._set_extra_kwds({'aspect_ratio': [1, 2, 2]})
            sage: S
        """
        self._extra_kwds = kwds

    def aspect_ratio(self, v=None):
        """
        Sets or gets the preferred aspect ratio of self.

        INPUT:

        - ``v`` -- (default: None) must be a list or tuple of length three,
          or the integer ``1``. If no arguments are provided then the
          default aspect ratio is returned.

        EXAMPLES::

            sage: D = dodecahedron()
            sage: D.aspect_ratio()
            [1.0, 1.0, 1.0]
            sage: D.aspect_ratio([1,2,3])
            sage: D.aspect_ratio()
            [1.0, 2.0, 3.0]
            sage: D.aspect_ratio(1)
            sage: D.aspect_ratio()
            [1.0, 1.0, 1.0]
        """
        if not v is None:
            if v == 1:
                v = (1,1,1)
            if not isinstance(v, (tuple, list)):
                raise TypeError("aspect_ratio must be a list or tuple of "
                                "length 3 or the integer 1")
            self._aspect_ratio = map(float, v)
        else:
            if self._aspect_ratio is None:
                self._aspect_ratio = [1.0,1.0,1.0]
            return self._aspect_ratio

    def frame_aspect_ratio(self, v=None):
        """
        Sets or gets the preferred frame aspect ratio of self.

        INPUT:

        - ``v`` -- (default: None) must be a list or tuple of length three,
          or the integer ``1``. If no arguments are provided then the
          default frame aspect ratio is returned.

        EXAMPLES::

            sage: D = dodecahedron()
            sage: D.frame_aspect_ratio()
            [1.0, 1.0, 1.0]
            sage: D.frame_aspect_ratio([2,2,1])
            sage: D.frame_aspect_ratio()
            [2.0, 2.0, 1.0]
            sage: D.frame_aspect_ratio(1)
            sage: D.frame_aspect_ratio()
            [1.0, 1.0, 1.0]
        """
        if not v is None:
            if v == 1:
                v = (1,1,1)
            if not isinstance(v, (tuple, list)):
                raise TypeError("frame_aspect_ratio must be a list or tuple of "
                                "length 3 or the integer 1")
            self._frame_aspect_ratio = map(float, v)
        else:
            if self._frame_aspect_ratio is None:
                self._frame_aspect_ratio = [1.0,1.0,1.0]
            return self._frame_aspect_ratio

    def _determine_frame_aspect_ratio(self, aspect_ratio):
        a_min, a_max = self._safe_bounding_box()
        return [(a_max[i] - a_min[i])*aspect_ratio[i] for i in range(3)]

    def _safe_bounding_box(self):
        """
        Returns a bounding box but where no side length is 0. This is used
        to avoid zero-division errors for pathological plots.

        EXAMPLES::

            sage: G = line3d([(0, 0, 0), (0, 0, 1)])
            sage: G.bounding_box()
            ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
            sage: G._safe_bounding_box()
            ([-1.0, -1.0, 0.0], [1.0, 1.0, 1.0])
        """
        a_min, a_max = self.bounding_box()
        a_min = list(a_min); a_max = list(a_max)
        for i in range(3):
            if a_min[i] == a_max[i]:
                a_min[i] = a_min[i] - 1
                a_max[i] = a_max[i] + 1
        return a_min, a_max


    def bounding_box(self):
        """
        Returns the lower and upper corners of a 3d bounding box for self.
        This is used for rendering and self should fit entirely within this
        box.

        Specifically, the first point returned should have x, y, and z
        coordinates should be the respective infimum over all points in self,
        and the second point is the supremum.

        The default return value is simply the box containing the origin.

        EXAMPLES::

            sage: sphere((1,1,1), 2).bounding_box()
            ((-1.0, -1.0, -1.0), (3.0, 3.0, 3.0))
            sage: G = line3d([(1, 2, 3), (-1,-2,-3)])
            sage: G.bounding_box()
            ((-1.0, -2.0, -3.0), (1.0, 2.0, 3.0))
        """
        return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))

    def transform(self, **kwds):
        """
        Apply a transformation to self, where the inputs are passed onto a
        TransformGroup object. Mostly for internal use; see the translate,
        scale, and rotate methods for more details.

        EXAMPLES::

            sage: sphere((0,0,0), 1).transform(trans=(1, 0, 0), scale=(2,3,4)).bounding_box()
            ((-1.0, -3.0, -4.0), (3.0, 3.0, 4.0))
        """
        return TransformGroup([self], **kwds)

    def translate(self, *x):
        """
        Return self translated by the given vector (which can be given either
        as a 3-iterable or via positional arguments).

        EXAMPLES::

            sage: icosahedron() + sum(icosahedron(opacity=0.25).translate(2*n, 0, 0) for n in [1..4])
            sage: icosahedron() + sum(icosahedron(opacity=0.25).translate([-2*n, n, n^2]) for n in [1..4])

        TESTS::

            sage: G = sphere((0, 0, 0), 1)
            sage: G.bounding_box()
            ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
            sage: G.translate(0, 0, 1).bounding_box()
            ((-1.0, -1.0, 0.0), (1.0, 1.0, 2.0))
            sage: G.translate(-1, 5, 0).bounding_box()
            ((-2.0, 4.0, -1.0), (0.0, 6.0, 1.0))
        """
        if len(x)==1:
            x = x[0]
        return self.transform(trans=x)

    def scale(self, *x):
        """
        Returns self scaled in the x, y, and z directions.

        EXAMPLES::

            sage: G = dodecahedron() + dodecahedron(opacity=.5).scale(2)
            sage: G.show(aspect_ratio=1)
            sage: G = icosahedron() + icosahedron(opacity=.5).scale([1, 1/2, 2])
            sage: G.show(aspect_ratio=1)

        TESTS::

            sage: G = sphere((0, 0, 0), 1)
            sage: G.scale(2)
            sage: G.scale(1, 2, 1/2).show(aspect_ratio=1)
            sage: G.scale(2).bounding_box()
            ((-2.0, -2.0, -2.0), (2.0, 2.0, 2.0))
        """
        if isinstance(x[0], (tuple, list)):
            x = x[0]
        return self.transform(scale=x)

    def rotate(self, v, theta):
        """
        Returns self rotated about the vector `v` by `theta` radians.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: v = (1,2,3)
            sage: G = arrow3d((0, 0, 0), v)
            sage: G += Cone(1/5, 1).translate((0, 0, 2))
            sage: C = Cone(1/5, 1, opacity=.25).translate((0, 0, 2))
            sage: G += sum(C.rotate(v, pi*t/4) for t in [1..7])
            sage: G.show(aspect_ratio=1)

            sage: from sage.plot.plot3d.shapes import Box
            sage: Box(1/3, 1/5, 1/7).rotate((1, 1, 1), pi/3).show(aspect_ratio=1)
        """
        vx, vy, vz = v
        return self.transform(rot=[vx, vy, vz, theta])

    def rotateX(self, theta):
        """
        Returns self rotated about the `x`-axis by the given angle.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: G = Cone(1/5, 1) + Cone(1/5, 1, opacity=.25).rotateX(pi/2)
            sage: G.show(aspect_ratio=1)
        """
        return self.rotate((1,0,0), theta)

    def rotateY(self, theta):
        """
        Returns self rotated about the `y`-axis by the given angle.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cone
            sage: G = Cone(1/5, 1) + Cone(1/5, 1, opacity=.25).rotateY(pi/3)
            sage: G.show(aspect_ratio=1)
        """
        return self.rotate((0,1,0), theta)

    def rotateZ(self, theta):
        """
        Returns self rotated about the `z`-axis by the given angle.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Box
            sage: G = Box(1/2, 1/3, 1/5) + Box(1/2, 1/3, 1/5, opacity=.25).rotateZ(pi/5)
            sage: G.show(aspect_ratio=1)
        """
        return self.rotate((0,0,1), theta)


    def viewpoint(self):
        """
        Returns the viewpoint of this plot. Currently only a stub for x3d.

        EXAMPLES::

            sage: type(dodecahedron().viewpoint())
            <class 'sage.plot.plot3d.base.Viewpoint'>
        """
        # This should probably be reworked somehow.
        return Viewpoint(0,0,6)

    def default_render_params(self):
        """
        Returns an instance of RenderParams suitable for plotting this object.

        EXAMPLES::

            sage: type(dodecahedron().default_render_params())
            <class 'sage.plot.plot3d.base.RenderParams'>
        """
        return RenderParams(ds=.075)

    def testing_render_params(self):
        """
        Returns an instance of RenderParams suitable for testing this object.
        In particular, it opens up '/dev/null' as an auxiliary zip file for jmol.

        EXAMPLES::

            sage: type(dodecahedron().testing_render_params())
            <class 'sage.plot.plot3d.base.RenderParams'>
        """
        params = RenderParams(ds=.075)
        params.output_archive = zipfile.ZipFile('/dev/null', 'w', zipfile.ZIP_STORED, True)
        return params

    def x3d(self):
        """
        An x3d scene file (as a string) containing the this object.

        EXAMPLES::

            sage: print sphere((1, 2, 3), 5).x3d()
            <X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
            <head>
            <meta name='title' content='sage3d'/>
            </head>
            <Scene>
            <Viewpoint position='0 0 6'/>
            <Transform translation='1 2 3'>
            <Shape><Sphere radius='5.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
            </Scene>
            </X3D>

            sage: G = icosahedron() + sphere((0,0,0), 0.5, color='red')
            sage: print G.x3d()
            <X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
            <head>
            <meta name='title' content='sage3d'/>
            </head>
            <Scene>
            <Viewpoint position='0 0 6'/>
            <Shape>
            <IndexedFaceSet coordIndex='...'>
              <Coordinate point='...'/>
            </IndexedFaceSet>
            <Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            <Transform translation='0 0 0'>
            <Shape><Sphere radius='0.5'/><Appearance><Material diffuseColor='1.0 0.0 0.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
            </Scene>
            </X3D>

        """
        return """
<X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
<head>
<meta name='title' content='sage3d'/>
</head>
<Scene>
%s
%s
</Scene>
</X3D>
"""%(self.viewpoint().x3d_str(), self.x3d_str())

    def tachyon(self):
        """
        An tachyon input file (as a string) containing the this object.

        EXAMPLES::

            sage: print sphere((1, 2, 3), 5, color='yellow').tachyon()
            begin_scene
            resolution 400 400
                     camera
                    ...
                  plane
                    center -2000 -1000 -500
                    normal 2.3 2.4 2.0
                    TEXTURE
                        AMBIENT 1.0 DIFFUSE 1.0 SPECULAR 1.0 OPACITY 1.0
                        COLOR 1.0 1.0 1.0
                        TEXFUNC 0
                Texdef texture...
              Ambient 0.333333333333 Diffuse 0.666666666667 Specular 0.0 Opacity 1
               Color 1.0 1.0 0.0
               TexFunc 0
                Sphere center 1.0 -2.0 3.0 Rad 5.0 texture...
            end_scene

            sage: G = icosahedron(color='red') + sphere((1,2,3), 0.5, color='yellow')
            sage: G.show(viewer='tachyon', frame=false)
            sage: print G.tachyon()
            begin_scene
            ...
            Texdef texture...
              Ambient 0.333333333333 Diffuse 0.666666666667 Specular 0.0 Opacity 1
               Color 1.0 1.0 0.0
               TexFunc 0
            TRI V0 ...
            Sphere center 1.0 -2.0 3.0 Rad 0.5 texture...
            end_scene
        """

        render_params = self.default_render_params()
        # switch from LH to RH coords to be consistent with java rendition
        render_params.push_transform(Transformation(scale=[1,-1,1]))
        return """
begin_scene
resolution 400 400

         camera
            zoom 1.0
            aspectratio 1.0
            antialiasing %s
            raydepth 8
            center  2.3 2.4 2.0
            viewdir  -2.3 -2.4 -2.0
            updir  0.0 0.0 1.0
         end_camera


      light center  4.0 3.0 2.0
            rad 0.2
            color  1.0 1.0 1.0

      plane
        center -2000 -1000 -500
        normal 2.3 2.4 2.0
        TEXTURE
            AMBIENT 1.0 DIFFUSE 1.0 SPECULAR 1.0 OPACITY 1.0
            COLOR 1.0 1.0 1.0
            TEXFUNC 0

    %s

    %s

end_scene""" % (render_params.antialiasing,
               "\n".join(sorted([t.tachyon_str() for t in self.texture_set()])),
               "\n".join(flatten_list(self.tachyon_repr(render_params))))

    def obj(self):
        """
        An .obj scene file (as a string) containing the this object. A
        .mtl file of the same name must also be produced for coloring.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import ColorCube
            sage: print ColorCube(1, ['red', 'yellow', 'blue']).obj()
            g obj_1
            usemtl ...
            v 1 1 1
            v -1 1 1
            v -1 -1 1
            v 1 -1 1
            f 1 2 3 4
            ...
            g obj_6
            usemtl ...
            v -1 -1 1
            v -1 1 1
            v -1 1 -1
            v -1 -1 -1
            f 21 22 23 24
        """
        return "\n".join(flatten_list([self.obj_repr(self.default_render_params()), ""]))

    def export_jmol(self, filename='jmol_shape.jmol', force_reload=False,
                    zoom=100, spin=False, background=(1,1,1), stereo=False,
                    mesh=False, dots=False,
                    perspective_depth = True,
                    orientation = (-764,-346,-545,76.39), **ignored_kwds):
                    # orientation chosen to look same as tachyon
        """
        A jmol scene consists of a script which refers to external files.
        Fortunately, we are able to put all of them in a single zip archive,
        which is the output of this call.

        EXAMPLES::

            sage: out_file = tmp_filename(ext=".jmol")
            sage: G = sphere((1, 2, 3), 5) + cube() + sage.plot.plot3d.shapes.Text("hi")
            sage: G.export_jmol(out_file)
            sage: import zipfile
            sage: z = zipfile.ZipFile(out_file)
            sage: z.namelist()
            ['obj_...pmesh', 'SCRIPT']

            sage: print z.read('SCRIPT')
            data "model list"
            2
            empty
            Xx 0 0 0
            Xx 5.5 5.5 5.5
            end "model list"; show data
            select *
            wireframe off; spacefill off
            set labelOffset 0 0
            background [255,255,255]
            spin OFF
            moveto 0 -764 -346 -545 76.39
            centerAt absolute {0 0 0}
            zoom 100
            frank OFF
            set perspectivedepth ON
            isosurface sphere_1  center {1.0 2.0 3.0} sphere 5.0
            color isosurface  [102,102,255]
            pmesh obj_... "obj_...pmesh"
            color pmesh  [102,102,255]
            select atomno = 1
            color atom  [102,102,255]
            label "hi"
            isosurface fullylit; pmesh o* fullylit; set antialiasdisplay on;

            sage: print z.read(z.namelist()[0])
            24
            0.5 0.5 0.5
            -0.5 0.5 0.5
            ...
            -0.5 -0.5 -0.5
            6
            5
            0
            1
            ...
        """
        render_params = self.default_render_params()
        render_params.mesh = mesh
        render_params.dots = dots
        render_params.output_file = filename
        render_params.force_reload = render_params.randomize_counter = force_reload
        render_params.output_archive = zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED, True)
        # Render the data
        all = flatten_list([self.jmol_repr(render_params), ""])

        f = StringIO()

        if render_params.atom_list:
            # Load the atom model
            f.write('data "model list"\n')
            f.write('%s\nempty\n' % (len(render_params.atom_list) + 1))
            for atom in render_params.atom_list:
                f.write('Xx %s %s %s\n' % atom)
            f.write('Xx 5.5 5.5 5.5\n') # so the zoom fits the box
            f.write('end "model list"; show data\n')
            f.write('select *\n')
            f.write('wireframe off; spacefill off\n')
            f.write('set labelOffset 0 0\n')


        # Set the scene background color
        f.write('background [%s,%s,%s]\n'%tuple([int(a*255) for a in background]))
        if spin:
            f.write('spin ON\n')
        else:
            f.write('spin OFF\n')
        if stereo:
            if stereo is True: stereo = "redblue"
            f.write('stereo %s\n' % stereo)
        if orientation:
            f.write('moveto 0 %s %s %s %s\n'%tuple(orientation))

        f.write('centerAt absolute {0 0 0}\n')
        f.write('zoom %s\n'%zoom)
        f.write('frank OFF\n') # jmol logo

        if perspective_depth:
            f.write('set perspectivedepth ON\n')
        else:
            f.write('set perspectivedepth OFF\n')

        # Put the rest of the object in
        f.write("\n".join(all))
        # Make sure the lighting is correct
        f.write("isosurface fullylit; pmesh o* fullylit; set antialiasdisplay on;\n")

        render_params.output_archive.writestr('SCRIPT', f.getvalue())
        render_params.output_archive.close()

    def json_repr(self, render_params):
        """
        A (possibly nested) list of strings. Each entry is formatted as JSON, so
        that a JavaScript client could eval it and get an object. Each object
        has fields to encapsulate the faces and vertices of self. This
        representation is intended to be consumed by the canvas3d viewer backend.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.json_repr(G.default_render_params())
            []
        """
        return []

    def jmol_repr(self, render_params):
        r"""
        A (possibly nested) list of strings which will be concatenated and
        used by jmol to render self. (Nested lists of strings are used
        because otherwise all the intermediate concatenations can kill
        performance). This may refer to several remove files, which
        are stored in render_parames.output_archive.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.jmol_repr(G.default_render_params())
            []
            sage: G = sphere((1, 2, 3))
            sage: G.jmol_repr(G.default_render_params())
            [['isosurface sphere_1  center {1.0 2.0 3.0} sphere 1.0\ncolor isosurface  [102,102,255]']]
        """
        return []

    def tachyon_repr(self, render_params):
        r"""
        A (possibly nested) list of strings which will be concatenated and
        used by tachyon to render self. (Nested lists of strings are used
        because otherwise all the intermediate concatenations can kill
        performance). This may include a reference to color information which
        is stored elsewhere.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.tachyon_repr(G.default_render_params())
            []
            sage: G = sphere((1, 2, 3))
            sage: G.tachyon_repr(G.default_render_params())
            ['Sphere center 1.0 2.0 3.0 Rad 1.0 texture...']
        """
        return []

    def obj_repr(self, render_params):
        """
        A (possibly nested) list of strings which will be concatenated and
        used to construct an .obj file of self. (Nested lists of strings are
        used because otherwise all the intermediate concatenations can kill
        performance). This may include a reference to color information which
        is stored elsewhere.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.obj_repr(G.default_render_params())
            []
            sage: G = cube()
            sage: G.obj_repr(G.default_render_params())
            ['g obj_1',
             'usemtl ...',
             ['v 0.5 0.5 0.5',
              'v -0.5 0.5 0.5',
              'v -0.5 -0.5 0.5',
              'v 0.5 -0.5 0.5',
              'v 0.5 0.5 -0.5',
              'v -0.5 0.5 -0.5',
              'v 0.5 -0.5 -0.5',
              'v -0.5 -0.5 -0.5'],
             ['f 1 2 3 4',
              'f 1 5 6 2',
              'f 1 4 7 5',
              'f 6 5 7 8',
              'f 7 4 3 8',
              'f 3 2 6 8'],
             []]
        """
        return []

    def texture_set(self):
        """
        Often the textures of a 3d file format are kept separate from the
        objects themselves. This function returns the set of textures used,
        so they can be defined in a preamble or separate file.

        EXAMPLES::

            sage: sage.plot.plot3d.base.Graphics3d().texture_set()
            set([])

            sage: G = tetrahedron(color='red') + tetrahedron(color='yellow') + tetrahedron(color='red', opacity=0.5)
            sage: [t for t in G.texture_set() if t.color == colors.red] # we should have two red textures
            [Texture(texture..., red, ff0000), Texture(texture..., red, ff0000)]
            sage: [t for t in G.texture_set() if t.color == colors.yellow] # ...and one yellow
            [Texture(texture..., yellow, ffff00)]
        """
        return set()

    def mtl_str(self):
        """
        Returns the contents of a .mtl file, to be used to provide coloring
        information for an .obj file.

        EXAMPLES::
            sage: G = tetrahedron(color='red') + tetrahedron(color='yellow', opacity=0.5)
            sage: print G.mtl_str()
            newmtl ...
            Ka 0.5 5e-06 5e-06
            Kd 1.0 1e-05 1e-05
            Ks 0.0 0.0 0.0
            illum 1
            Ns 1
            d 1
            newmtl ...
            Ka 0.5 0.5 5e-06
            Kd 1.0 1.0 1e-05
            Ks 0.0 0.0 0.0
            illum 1
            Ns 1
            d 0.500000000000000
        """
        return "\n\n".join(sorted([t.mtl_str() for t in self.texture_set()])) + "\n"

    def flatten(self):
        """
        Try to reduce the depth of the scene tree by consolidating groups
        and transformations.

        The generic Graphics3d object can't be made flatter.

        EXAMPLES::

            sage: G = sage.plot.plot3d.base.Graphics3d()
            sage: G.flatten() is G
            True
        """
        return self

    def _rescale_for_frame_aspect_ratio_and_zoom(self, b, frame_aspect_ratio, zoom):
        if frame_aspect_ratio is None:
            return (b*zoom,b*zoom,b*zoom), (-b*zoom,-b*zoom,-b*zoom)
        box = [b*w for w in frame_aspect_ratio]
        # Now take the maximum length in box and rescale to b.
        s = b / max(box)
        box_max = tuple([s*w*zoom for w in box])
        box_min = tuple([-w*zoom for w in box_max])
        return box_min, box_max

    def _prepare_for_jmol(self, frame, axes, frame_aspect_ratio, aspect_ratio, zoom):
        from sage.plot.plot import EMBEDDED_MODE
        if EMBEDDED_MODE:
            s = 6
        else:
            s = 3
        box_min, box_max = self._rescale_for_frame_aspect_ratio_and_zoom(s, frame_aspect_ratio, zoom)
        a_min, a_max = self._box_for_aspect_ratio(aspect_ratio, box_min, box_max)
        return self._transform_to_bounding_box(box_min, box_max, a_min, a_max, frame=frame,
                                               axes=axes, thickness=1,
                                               labels = True)   # jmol labels are implemented

    def _prepare_for_tachyon(self, frame, axes, frame_aspect_ratio, aspect_ratio, zoom):
        box_min, box_max = self._rescale_for_frame_aspect_ratio_and_zoom(1.0, frame_aspect_ratio, zoom)
        a_min, a_max = self._box_for_aspect_ratio(aspect_ratio, box_min, box_max)
        return self._transform_to_bounding_box(box_min, box_max, a_min, a_max,
                                               frame=frame, axes=axes, thickness=.75,
                                               labels = False)  # no tachyon text implemented yet

    def _box_for_aspect_ratio(self, aspect_ratio, box_min, box_max):
        # 1. Find a box around self so that when self gets rescaled into the
        # box defined by box_min, box_max, it has the right aspect ratio
        a_min, a_max = self._safe_bounding_box()

        if aspect_ratio == "automatic" or aspect_ratio == [1.0]*3:
            return a_min, a_max

        longest_side = 0; longest_length = a_max[0] - a_min[0]
        shortest_side = 0; shortest_length = a_max[0] - a_min[0]

        for i in range(3):
            s = a_max[i] - a_min[i]
            if s > longest_length:
                longest_length = s
                longest_side = i
            if s < shortest_length:
                shortest_length = s
                shortest_side = i

        # 2. Rescale aspect_ratio so the shortest side is 1.
        r = float(aspect_ratio[shortest_side])
        aspect_ratio = [a/r for a in aspect_ratio]

        # 3. Extend the bounding box of self by rescaling so the sides
        # have the same ratio as aspect_ratio, and without changing
        # the longest side.
        long_box_side = box_max[longest_side] - box_min[longest_side]
        sc = [1.0,1.0,1.0]
        for i in range(3):
            # compute the length we want:
            new_length = longest_length / aspect_ratio[i]
            # change the side length by a_min and a_max so
            # that a_max[i] - a_min[i] = new_length

            # We have to take into account the ratio of the
            # sides after transforming to the bounding box.
            z = long_box_side / (box_max[i] - box_min[i])
            w = new_length / ((a_max[i] - a_min[i]) * z)
            sc[i] = w

        w = min(sc)
        sc = [z/w for z in sc]
        for i in range(3):
            a_min[i] *= sc[i]
            a_max[i] *= sc[i]

        return a_min, a_max

    def _transform_to_bounding_box(self, xyz_min, xyz_max, a_min, a_max, frame, axes, thickness, labels):

        a_min_orig = a_min; a_max_orig = a_max

        # Rescale in each direction
        scale = [float(xyz_max[i] - xyz_min[i]) / (a_max[i] - a_min[i]) for i in range(3)]
        X = self.scale(scale)
        a_min = [scale[i]*a_min[i] for i in range(3)]
        a_max = [scale[i]*a_max[i] for i in range(3)]

        # Translate so lower left corner of original bounding box
        # is in the right spot
        T = [xyz_min[i] - a_min[i] for i in range(3)]
        X = X.translate(T)
        if frame:
            from shapes2 import frame3d, frame_labels
            F = frame3d(xyz_min, xyz_max, opacity=0.5, color=(0,0,0), thickness=thickness)
            if labels:
                F += frame_labels(xyz_min, xyz_max, a_min_orig, a_max_orig)

            X += F

        if axes:
            # draw axes
            from shapes import arrow3d
            A = (arrow3d((min(0,a_min[0]),0, 0), (max(0,a_max[0]), 0,0),
                             thickness, color="blue"),
                 arrow3d((0,min(0,a_min[1]), 0), (0, max(0,a_max[1]), 0),
                             thickness, color="blue"),
                 arrow3d((0, 0, min(0,a_min[2])), (0, 0, max(0,a_max[2])),
                             thickness, color="blue"))
            X += sum(A).translate([-z for z in T])

        return X

    def _process_viewing_options(self, kwds):
        """
        Process viewing options (the keywords passed to show()) and return a new
        dictionary. Defaults will be filled in for missing options and taken from
        self._extra_kwds as well. Options that have the value "automatic" will be
        automatically determined. Finally, the provided dictionary is modified
        to remove all of the keys that were used -- so that the unused keywords
        can be used elsewhere.
        """
        opts = {}
        opts.update(SHOW_DEFAULTS)
        if self._extra_kwds is not None:
            opts.update(self._extra_kwds)
        opts.update(kwds)

        # Remove all of the keys that are viewing options, since the remaining
        # kwds might be passed on.
        for key_to_remove in SHOW_DEFAULTS.keys():
            kwds.pop(key_to_remove, None)

        # deal with any aspect_ratio instances passed from the default options to plot
        if opts['aspect_ratio'] == 'auto':
            opts['aspect_ratio'] = 'automatic'
        if opts['aspect_ratio'] != 'automatic':
            # We need this round about way to make sure that we do not
            # store the aspect ratio that was passed on to show() by the
            # user. We let the .aspect_ratio() method take care of the
            # validity of the arguments that was passed on to show()
            original_aspect_ratio = self.aspect_ratio()
            self.aspect_ratio(opts['aspect_ratio'])
            opts['aspect_ratio'] = self.aspect_ratio()
            self.aspect_ratio(original_aspect_ratio)

        if opts['frame_aspect_ratio'] == 'automatic':
            if opts['aspect_ratio'] != 'automatic':
                # Set the aspect_ratio of the frame to be the same as that
                # of the object we are rendering given the aspect_ratio
                # we'll use for it.
                opts['frame_aspect_ratio'] = \
                    self._determine_frame_aspect_ratio(opts['aspect_ratio'])
            else:
                opts['frame_aspect_ratio'] = self.frame_aspect_ratio()
        else:
            # We need this round about way to make sure that we do not
            # store the frame aspect ratio that was passed on to show() by
            # the user. We let the .frame_aspect_ratio() method take care
            # of the validity of the arguments that was passed on to show()
            original_aspect_ratio = self.frame_aspect_ratio()
            self.frame_aspect_ratio(opts['frame_aspect_ratio'])
            opts['frame_aspect_ratio'] = self.frame_aspect_ratio()
            self.frame_aspect_ratio(original_aspect_ratio)

        if opts['aspect_ratio'] == 'automatic':
            opts['aspect_ratio'] = self.aspect_ratio()

        if not isinstance(opts['figsize'], (list,tuple)):
            opts['figsize'] = [opts['figsize'], opts['figsize']]

        return opts

    def show(self, **kwds):
        """
        INPUT:


        -  ``viewer`` - string (default: 'jmol'), how to view
           the plot

           * 'jmol': Interactive 3D viewer using Java

           * 'tachyon': Ray tracer generates a static PNG image

           * 'java3d': Interactive OpenGL based 3D

           * 'canvas3d': Web-based 3D viewer powered by JavaScript and
             <canvas> (notebook only)

        -  ``filename`` - string (default: a temp file); file
           to save the image to

        -  ``verbosity`` - display information about rendering
           the figure

        -  ``figsize`` - (default: 5); x or pair [x,y] for
           numbers, e.g., [5,5]; controls the size of the output figure. E.g.,
           with Tachyon the number of pixels in each direction is 100 times
           figsize[0]. This is ignored for the jmol embedded renderer.

        -  ``aspect_ratio`` - (default: "automatic") - aspect
           ratio of the coordinate system itself. Give [1,1,1] to make spheres
           look round.

        -  ``frame_aspect_ratio`` - (default: "automatic")
           aspect ratio of frame that contains the 3d scene.

        -  ``zoom`` - (default: 1) how zoomed in

        -  ``frame`` - (default: True) if True, draw a
           bounding frame with labels

        -  ``axes`` - (default: False) if True, draw coordinate
           axes

        -  ``**kwds`` - other options, which make sense for particular
           rendering engines

        CHANGING DEFAULTS: Defaults can be uniformly changed by importing a
        dictionary and changing it. For example, here we change the default
        so images display without a frame instead of with one::

            sage: from sage.plot.plot3d.base import SHOW_DEFAULTS
            sage: SHOW_DEFAULTS['frame'] = False

        This sphere will not have a frame around it::

            sage: sphere((0,0,0))

        We change the default back::

            sage: SHOW_DEFAULTS['frame'] = True

        Now this sphere is enclosed in a frame::

            sage: sphere((0,0,0))

        EXAMPLES: We illustrate use of the ``aspect_ratio`` option::

            sage: x, y = var('x,y')
            sage: p = plot3d(2*sin(x*y), (x, -pi, pi), (y, -pi, pi))
            sage: p.show(aspect_ratio=[1,1,1])

        This looks flattened, but filled with the plot::

            sage: p.show(frame_aspect_ratio=[1,1,1/16])

        This looks flattened, but the plot is square and smaller::

            sage: p.show(aspect_ratio=[1,1,1], frame_aspect_ratio=[1,1,1/8])

        This example shows indirectly that the defaults
        from :func:`~sage.plot.plot.plot` are dealt with properly::

            sage: plot(vector([1,2,3]))

        We use the 'canvas3d' backend from inside the notebook to get a view of
        the plot rendered inline using HTML canvas::

            sage: p.show(viewer='canvas3d')
        """

        opts = self._process_viewing_options(kwds)

        viewer = opts['viewer']
        verbosity = opts['verbosity']
        figsize = opts['figsize']
        aspect_ratio = opts['aspect_ratio'] # this necessarily has a value now
        frame_aspect_ratio = opts['frame_aspect_ratio']
        zoom = opts['zoom']
        frame = opts['frame']
        axes = opts['axes']

        import sage.misc.misc
        try:
            filename = kwds.pop('filename')
        except KeyError:
            filename = tmp_filename()

        from sage.plot.plot import EMBEDDED_MODE
        from sage.doctest import DOCTEST_MODE
        ext = None

        # Tachyon resolution options
        if DOCTEST_MODE:
            opts = '-res 10 10'
            filename = os.path.join(sage.misc.misc.SAGE_TMP, "tmp")
        elif EMBEDDED_MODE:
            opts = '-res %s %s'%(figsize[0]*100, figsize[1]*100)
            filename = sage.misc.temporary_file.graphics_filename()[:-4]
        else:
            opts = '-res %s %s'%(figsize[0]*100, figsize[1]*100)

        if DOCTEST_MODE or viewer=='tachyon' or (viewer=='java3d' and EMBEDDED_MODE):
            T = self._prepare_for_tachyon(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
            tachyon_rt(T.tachyon(), filename+".png", verbosity, True, opts)
            ext = "png"
            import sage.misc.viewer
            viewer_app = sage.misc.viewer.png_viewer()

        if DOCTEST_MODE or viewer=='java3d':
            f = open(filename+".obj", "w")
            f.write("mtllib %s.mtl\n" % filename)
            f.write(self.obj())
            f.close()
            f = open(filename+".mtl", "w")
            f.write(self.mtl_str())
            f.close()
            ext = "obj"
            viewer_app = os.path.join(sage.misc.misc.SAGE_LOCAL, "bin/sage3d")

        if DOCTEST_MODE or viewer=='jmol':
            # Temporary hack: encode the desired applet size in the end of the filename:
            # (This will be removed once we have dynamic resizing of applets in the browser.)
            base, ext = os.path.splitext(filename)
            fg = figsize[0]
            filename = '%s-size%s%s'%(base, fg*100, ext)

            if EMBEDDED_MODE:
                ext = "jmol"
                # jmol doesn't seem to correctly parse the ?params part of a URL
                archive_name = "%s-%s.%s.zip" % (filename, randint(0, 1 << 30), ext)
            else:
                ext = "spt"
                archive_name = "%s.%s.zip" % (filename, ext)

            T = self._prepare_for_jmol(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
            T.export_jmol(archive_name, force_reload=EMBEDDED_MODE, zoom=zoom*100, **kwds)
            viewer_app = os.path.join(sage.misc.misc.SAGE_LOCAL, "bin", "jmol")

            # If the server has a Java installation we can make better static images with Jmol
            # Test for Java then make image with Jmol or Tachyon if no JavaVM
            if EMBEDDED_MODE:
                # We need a script for the Notebook.
                # When the notebook sees this file, it will know to
                # display the static file and the "Make Interactive"
                # button.
                import sagenb
                path = "cells/%s/%s" %(sagenb.notebook.interact.SAGE_CELL_ID, archive_name)
                with open(filename + '.' + ext, 'w') as f:
                    f.write('set defaultdirectory "%s"\n' % path)
                    f.write('script SCRIPT\n')

                # Filename for the static image
                png_path = '.jmol_images'
                sage.misc.misc.sage_makedirs(png_path)
                png_name = os.path.join(png_path, filename + ".jmol.png")

                from sage.interfaces.jmoldata import JmolData
                jdata = JmolData()
                if jdata.is_jvm_available():
                    # Java needs absolute paths
                    archive_name = os.path.abspath(archive_name)
                    png_name = os.path.abspath(png_name)
                    script = '''set defaultdirectory "%s"\nscript SCRIPT\n''' % archive_name
                    jdata.export_image(targetfile=png_name, datafile=script, image_type="PNG", figsize=fg)
                else:
                    # Render the image with tachyon
                    T = self._prepare_for_tachyon(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
                    tachyon_rt(T.tachyon(), png_name, verbosity, True, opts)

        if viewer == 'canvas3d':
            T = self._prepare_for_tachyon(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
            data = flatten_list(T.json_repr(T.default_render_params()))
            f = open(filename + '.canvas3d', 'w')
            f.write('[%s]' % ','.join(data))
            f.close()
            ext = 'canvas3d'

        if ext is None:
            raise ValueError, "Unknown 3d plot type: %s" % viewer

        if not DOCTEST_MODE and not EMBEDDED_MODE:
            if verbosity:
                pipes = "2>&1"
            else:
                pipes = "2>/dev/null 1>/dev/null &"
            os.system('%s "%s.%s" %s' % (viewer_app, filename, ext, pipes))

    def save(self, filename, **kwds):
        """
        Save the graphic to an image file (of type: PNG, BMP, GIF, PPM, or TIFF)
        rendered using Tachyon, or pickle it (stored as an SOBJ so you can load it
        later) depending on the file extension you give the filename.

        INPUT:

        - ``filename`` - Specify where to save the image or object.

        - ``**kwds`` - When specifying an image file to be rendered by Tachyon,
          any of the viewing options accepted by show() are valid as keyword
          arguments to this function and they will behave in the same way.
          Accepted keywords include: ``viewer``, ``verbosity``, ``figsize``,
          ``aspect_ratio``, ``frame_aspect_ratio``, ``zoom``, ``frame``, and
          ``axes``. Default values are provided.

        EXAMPLES::

            sage: f = tmp_filename() + '.png'
            sage: G = sphere()
            sage: G.save(f)

        We demonstrate using keyword arguments to control the appearance of the
        output image::

            sage: G.save(f, zoom=2, figsize=[5, 10])

        But some extra parameters don't make since (like ``viewer``, since
        rendering is done using Tachyon only). They will be ignored::

            sage: G.save(f, viewer='jmol') # Looks the same

        Since Tachyon only outputs PNG images, PIL will be used to convert to
        alternate formats::

            sage: cube().save(tmp_filename(ext='.gif'))
        """
        ext = os.path.splitext(filename)[1].lower()
        if ext == '' or ext == '.sobj':
            SageObject.save(self, filename)
            return
        elif ext in ['.bmp', '.png', '.gif', '.ppm', '.tiff', '.tif', '.jpg', '.jpeg']:
            opts = self._process_viewing_options(kwds)
            T = self._prepare_for_tachyon(
                opts['frame'], opts['axes'], opts['frame_aspect_ratio'],
                opts['aspect_ratio'], opts['zoom']
            )

            if ext == 'png':
                # No conversion is necessary
                out_filename = filename
            else:
                # Save to a temporary file, and then convert using PIL
                out_filename = sage.misc.temporary_file.tmp_filename(ext=ext)
            tachyon_rt(T.tachyon(), out_filename, opts['verbosity'], True,
                '-res %s %s' % (opts['figsize'][0]*100, opts['figsize'][1]*100))
            if ext != 'png':
                import Image
                Image.open(out_filename).save(filename)
        else:
            raise ValueError, 'filetype not supported by save()'



# if you add any default parameters you must update some code below
SHOW_DEFAULTS = {'viewer':'jmol',
                 'verbosity':0,
                 'figsize':5,
                 'aspect_ratio':"automatic",
                 'frame_aspect_ratio':"automatic",
                 'zoom':1,
                 'frame':True,
                 'axes':False}




class Graphics3dGroup(Graphics3d):
    """
    This class represents a collection of 3d objects. Usually they are formed
    implicitly by summing.
    """
    def __init__(self, all=(), rot=None, trans=None, scale=None, T=None):
        """
        EXAMPLES::

            sage: sage.plot.plot3d.base.Graphics3dGroup([icosahedron(), dodecahedron(opacity=.5)])
            sage: type(icosahedron() + dodecahedron(opacity=.5))
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        self.all = list(all)
        self.frame_aspect_ratio(optimal_aspect_ratios([a.frame_aspect_ratio() for a in all]))
        self.aspect_ratio(optimal_aspect_ratios([a.aspect_ratio() for a in all]))
        self._set_extra_kwds(optimal_extra_kwds([a._extra_kwds for a in all if a._extra_kwds is not None]))

    def __add__(self, other):
        """
        We override this here to make large sums more efficient.

        EXAMPLES::
            sage: G = sum(tetrahedron(opacity=1-t/11).translate(t, 0, 0) for t in range(10))
            sage: G
            sage: len(G.all)
            10
        """
        if type(self) is Graphics3dGroup and isinstance(other, Graphics3d):
            self.all.append(other)
            return self
        else:
            return Graphics3d.__add__(self, other)

    def bounding_box(self):
        """
        Box that contains the bounding boxes of
        all the objects that make up self.

        EXAMPLES::

            sage: A = sphere((0,0,0), 5)
            sage: B = sphere((1, 5, 10), 1)
            sage: A.bounding_box()
            ((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0))
            sage: B.bounding_box()
            ((0.0, 4.0, 9.0), (2.0, 6.0, 11.0))
            sage: (A+B).bounding_box()
            ((-5.0, -5.0, -5.0), (5.0, 6.0, 11.0))
            sage: (A+B).show(aspect_ratio=1, frame=True)

            sage: sage.plot.plot3d.base.Graphics3dGroup([]).bounding_box()
            ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        """
        if len(self.all) == 0:
            return Graphics3d.bounding_box(self)
        v = [obj.bounding_box() for obj in self.all]
        return min3([a[0] for a in v]), max3([a[1] for a in v])

    def transform(self, **kwds):
        """
        Transforming this entire group simply makes a transform group with
        the same contents.

        EXAMPLES::

            sage: G = dodecahedron(color='red', opacity=.5) + icosahedron(color='blue')
            sage: G
            sage: G.transform(scale=(2,1/2,1))
            sage: G.transform(trans=(1,1,3))
        """
        T = TransformGroup(self.all, **kwds)
        T._set_extra_kwds(self._extra_kwds)
        return T

    def set_texture(self, **kwds):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red', opacity=.5) + icosahedron((3, 0, 0), color='blue')
            sage: G
            sage: G.set_texture(color='yellow')
            sage: G
        """
        for g in self.all:
            g.set_texture(**kwds)

    def json_repr(self, render_params):
        """
        The JSON representation of a group is simply the concatenation of the
        representations of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1, 2, 3))
            sage: G.json_repr(G.default_render_params())
            [[["{vertices:..."]], [["{vertices:..."]]]
        """
        return [g.json_repr(render_params) for g in self.all]

    def tachyon_repr(self, render_params):
        """
        The tachyon representation of a group is simply the concatenation of
        the representations of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: G.tachyon_repr(G.default_render_params())
            [['Sphere center 0.0 0.0 0.0 Rad 1.0 texture...'],
             ['Sphere center 1.0 2.0 3.0 Rad 1.0 texture...']]
        """
        return [g.tachyon_repr(render_params) for g in self.all]

    def x3d_str(self):
        """
        The x3d representation of a group is simply the concatenation of
        the representation of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: print G.x3d_str()
            <Transform translation='0 0 0'>
            <Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
            <Transform translation='1 2 3'>
            <Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
        """
        return "\n".join([g.x3d_str() for g in self.all])

    def obj_repr(self, render_params):
        """
        The obj representation of a group is simply the concatenation of
        the representation of its objects.

        EXAMPLES::

            sage: G = tetrahedron() + tetrahedron().translate(10, 10, 10)
            sage: G.obj_repr(G.default_render_params())
            [['g obj_1',
              'usemtl ...',
              ['v 0 0 1',
               'v 0.942809 0 -0.333333',
               'v -0.471405 0.816497 -0.333333',
               'v -0.471405 -0.816497 -0.333333'],
              ['f 1 2 3', 'f 2 4 3', 'f 1 3 4', 'f 1 4 2'],
              []],
             [['g obj_2',
               'usemtl ...',
               ['v 10 10 11',
                'v 10.9428 10 9.66667',
                'v 9.5286 10.8165 9.66667',
                'v 9.5286 9.1835 9.66667'],
               ['f 5 6 7', 'f 6 8 7', 'f 5 7 8', 'f 5 8 6'],
               []]]]
        """
        return [g.obj_repr(render_params) for g in self.all]

    def jmol_repr(self, render_params):
        r"""
        The jmol representation of a group is simply the concatenation of
        the representation of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: G.jmol_repr(G.default_render_params())
            [[['isosurface sphere_1  center {0.0 0.0 0.0} sphere 1.0\ncolor isosurface  [102,102,255]']],
             [['isosurface sphere_2  center {1.0 2.0 3.0} sphere 1.0\ncolor isosurface  [102,102,255]']]]
        """
        return [g.jmol_repr(render_params) for g in self.all]

    def texture_set(self):
        """
        The texture set of a group is simply the union of the textures of
        all its objects.

        EXAMPLES::

            sage: G = sphere(color='red') + sphere(color='yellow')
            sage: [t for t in G.texture_set() if t.color == colors.red] # one red texture
            [Texture(texture..., red, ff0000)]
            sage: [t for t in G.texture_set() if t.color == colors.yellow] # one yellow texture
            [Texture(texture..., yellow, ffff00)]

            sage: T = sage.plot.plot3d.texture.Texture('blue'); T
            Texture(texture..., blue, 0000ff)
            sage: G = sphere(texture=T) + sphere((1, 1, 1), texture=T)
            sage: len(G.texture_set())
            1
        """
        return reduce(set.union, [g.texture_set() for g in self.all])

    def flatten(self):
        """
        Try to reduce the depth of the scene tree by consolidating groups
        and transformations.

        EXAMPLES::

            sage: G = sum([circle((0, 0), t) for t in [1..10]], sphere()); G
            sage: G.flatten()
            sage: len(G.all)
            2
            sage: len(G.flatten().all)
            11
        """
        if len(self.all) == 1:
            return self.all[0].flatten()
        all = []
        for g in self.all:
            g = g.flatten()
            if type(g) is Graphics3dGroup:
                all += g.all
            else:
                all.append(g)
        return Graphics3dGroup(all)



class TransformGroup(Graphics3dGroup):
    """
    This class is a container for a group of objects with a common transformation.
    """
    def __init__(self, all=[], rot=None, trans=None, scale=None, T=None):
        """
        EXAMPLES::

            sage: sage.plot.plot3d.base.TransformGroup([sphere()], trans=(1,2,3)) + point3d((0,0,0))

        The are usually constructed implicitly::

            sage: type(sphere((1,2,3)))
            <class 'sage.plot.plot3d.base.TransformGroup'>
            sage: type(dodecahedron().scale(2))
            <class 'sage.plot.plot3d.base.TransformGroup'>

        """
        Graphics3dGroup.__init__(self, all)
        self._rot = rot
        self._trans = trans
        if scale is not None and len(scale) == 1:
            if isinstance(scale, (tuple, list)):
                scale = scale[0]
            scale = (scale, scale, scale)
        self._scale = scale
        if T is not None:
            self.T = T
        self.frame_aspect_ratio(optimal_aspect_ratios([a.frame_aspect_ratio() for a in all]))
        self.aspect_ratio(optimal_aspect_ratios([a.aspect_ratio() for a in all]))
        self._set_extra_kwds(optimal_extra_kwds([a._extra_kwds for a in all if a._extra_kwds is not None]))

    def bounding_box(self):
        """
        Returns the bounding box of self, i.e. the box containing the
        contents of self after applying the transformation.

        EXAMPLES::

            sage: G = cube()
            sage: G.bounding_box()
            ((-0.5, -0.5, -0.5), (0.5, 0.5, 0.5))
            sage: G.scale(4).bounding_box()
            ((-2.0, -2.0, -2.0), (2.0, 2.0, 2.0))
            sage: G.rotateZ(pi/4).bounding_box()
            ((-0.7071067811865475, -0.7071067811865475, -0.5),
             (0.7071067811865475, 0.7071067811865475, 0.5))
        """
        try:
            return self._bounding_box
        except AttributeError:
            pass

        cdef Transformation T = self.get_transformation()
        w = sum([T.transform_bounding_box(obj.bounding_box()) for obj in self.all], ())
        self._bounding_box = point_list_bounding_box(w)
        return self._bounding_box

    def x3d_str(self):
        r"""
        To apply a transformation to a set of objects in x3d, simply make them
        all children of an x3d Transform node.

        EXAMPLES::

            sage: sphere((1,2,3)).x3d_str()
            "<Transform translation='1 2 3'>\n<Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance></Shape>\n\n</Transform>"
        """
        s = "<Transform"
        if self._rot is not None:
            s += " rotation='%s %s %s %s'"%tuple(self._rot)
        if self._trans is not None:
            s += " translation='%s %s %s'"%tuple(self._trans)
        if self._scale is not None:
            s += " scale='%s %s %s'"%tuple(self._scale)
        s += ">\n"
        s += Graphics3dGroup.x3d_str(self)
        s += "\n</Transform>"
        return s

    def json_repr(self, render_params):
        """
        Transformations are applied at the leaf nodes.

        EXAMPLES::

            sage: G = cube().rotateX(0.2)
            sage: G.json_repr(G.default_render_params())
            [["{vertices:[{x:0.5,y:0.589368,z:0.390699},..."]]
        """

        render_params.push_transform(self.get_transformation())
        rep = [g.json_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def tachyon_repr(self, render_params):
        """
        Transformations for Tachyon are applied at the leaf nodes.

        EXAMPLES::

            sage: G = sphere((1,2,3)).scale(2)
            sage: G.tachyon_repr(G.default_render_params())
            [['Sphere center 2.0 4.0 6.0 Rad 2.0 texture...']]
        """
        render_params.push_transform(self.get_transformation())
        rep = [g.tachyon_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def obj_repr(self, render_params):
        """
        Transformations for .obj files are applied at the leaf nodes.

        EXAMPLES::

            sage: G = cube().scale(4).translate(1, 2, 3)
            sage: G.obj_repr(G.default_render_params())
            [[['g obj_1',
               'usemtl ...',
               ['v 3 4 5',
                'v -1 4 5',
                'v -1 0 5',
                'v 3 0 5',
                'v 3 4 1',
                'v -1 4 1',
                'v 3 0 1',
                'v -1 0 1'],
               ['f 1 2 3 4',
                'f 1 5 6 2',
                'f 1 4 7 5',
                'f 6 5 7 8',
                'f 7 4 3 8',
                'f 3 2 6 8'],
               []]]]
        """
        render_params.push_transform(self.get_transformation())
        rep = [g.obj_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def jmol_repr(self, render_params):
        r"""
        Transformations for jmol are applied at the leaf nodes.

        EXAMPLES::

            sage: G = sphere((1,2,3)).scale(2)
            sage: G.jmol_repr(G.default_render_params())
            [[['isosurface sphere_1  center {2.0 4.0 6.0} sphere 2.0\ncolor isosurface  [102,102,255]']]]
        """
        render_params.push_transform(self.get_transformation())
        rep = [g.jmol_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def get_transformation(self):
        """
        Returns the actual transformation object associated with self.

        EXAMPLES::

            sage: G = sphere().scale(100)
            sage: T = G.get_transformation()
            sage: T.get_matrix()
            [100.0   0.0   0.0   0.0]
            [  0.0 100.0   0.0   0.0]
            [  0.0   0.0 100.0   0.0]
            [  0.0   0.0   0.0   1.0]
        """
        try:
            return self.T
        except AttributeError:
            self.T = Transformation(self._scale, self._rot, self._trans)
            return self.T

    def flatten(self):
        """
        Try to reduce the depth of the scene tree by consolidating groups
        and transformations.

        EXAMPLES::

            sage: G = sphere((1,2,3)).scale(100)
            sage: T = G.get_transformation()
            sage: T.get_matrix()
            [100.0   0.0   0.0   0.0]
            [  0.0 100.0   0.0   0.0]
            [  0.0   0.0 100.0   0.0]
            [  0.0   0.0   0.0   1.0]

            sage: G.flatten().get_transformation().get_matrix()
            [100.0   0.0   0.0 100.0]
            [  0.0 100.0   0.0 200.0]
            [  0.0   0.0 100.0 300.0]
            [  0.0   0.0   0.0   1.0]
        """
        G = Graphics3dGroup.flatten(self)
        if isinstance(G, TransformGroup):
            return TransformGroup(G.all, T=self.get_transformation() * G.get_transformation())
        elif isinstance(G, Graphics3dGroup):
            return TransformGroup(G.all, T=self.get_transformation())
        else:
            return TransformGroup([G], T=self.get_transformation())

    def transform(self, **kwds):
        """
        Transforming this entire group can be done by composing transformations.

        EXAMPLES::

            sage: G = dodecahedron(color='red', opacity=.5) + icosahedron(color='blue')
            sage: G
            sage: G.transform(scale=(2,1/2,1))
            sage: G.transform(trans=(1,1,3))
        """
        return Graphics3d.transform(self, **kwds)

class Viewpoint(Graphics3d):
    """
    This class represents a viewpoint, necessary for x3d.

    In the future, there could be multiple viewpoints, and they could have
    more properties. (Currently they only hold a position).
    """
    def __init__(self, *x):
        """
        EXAMPLES::

            sage: sage.plot.plot3d.base.Viewpoint(1, 2, 4).x3d_str()
            "<Viewpoint position='1 2 4'/>"
        """
        if isinstance(x[0], (tuple, list)):
            x = tuple(x[0])
        self.pos = x

    def x3d_str(self):
        """
        EXAMPLES::

            sage: sphere((0,0,0), 100).viewpoint().x3d_str()
            "<Viewpoint position='0 0 6'/>"
        """
        return "<Viewpoint position='%s %s %s'/>"%self.pos



cdef class PrimitiveObject(Graphics3d):
    """
    This is the base class for the non-container 3d objects.
    """
    def __init__(self, **kwds):
        if kwds.has_key('texture'):
            self.texture = kwds['texture']
            if not is_Texture(self.texture):
                self.texture = Texture(self.texture)
        else:
            self.texture = Texture(kwds)

    def set_texture(self, texture=None, **kwds):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red'); G
            sage: G.set_texture(color='yellow'); G
        """
        if not is_Texture(texture):
            texture = Texture(texture, **kwds)
        self.texture = texture

    def get_texture(self):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red')
            sage: G.get_texture()
            Texture(texture..., red, ff0000)
        """
        return self.texture

    def texture_set(self):
        """
        EXAMPLES::

            sage: G = dodecahedron(color='red')
            sage: G.texture_set()
            set([Texture(texture..., red, ff0000)])
        """
        return set([self.texture])

    def x3d_str(self):
        r"""
        EXAMPLES::

            sage: sphere().flatten().x3d_str()
            "<Transform>\n<Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1' specularColor='0.0 0.0 0.0'/></Appearance></Shape>\n\n</Transform>"
        """
        return "<Shape>" + self.x3d_geometry() + self.texture.x3d_str() + "</Shape>\n"

    def tachyon_repr(self, render_params):
        """
        Default behavior is to render the triangulation.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: G = Torus(1, .5)
            sage: G.tachyon_repr(G.default_render_params())
            ['TRI V0 0 1 0.5
            ...
            'texture...']
        """
        return self.triangulation().tachyon_repr(render_params)

    def obj_repr(self, render_params):
        """
        Default behavior is to render the triangulation.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: G = Torus(1, .5)
            sage: G.obj_repr(G.default_render_params())
            ['g obj_1',
             'usemtl ...',
             ['v 0 1 0.5',
             ...
              'f ...'],
             []]
        """
        return self.triangulation().obj_repr(render_params)

    def jmol_repr(self, render_params):
        r"""
        Default behavior is to render the triangulation. The actual polygon
        data is stored in a separate file.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Torus
            sage: G = Torus(1, .5)
            sage: G.jmol_repr(G.testing_render_params())
            ['pmesh obj_1 "obj_1.pmesh"\ncolor pmesh  [102,102,255]']
        """
        return self.triangulation().jmol_repr(render_params)



class BoundingSphere(SageObject):
    """
    A bounding sphere is like a bounding box, but is simpler to deal with and
    behaves better under rotations.
    """
    def __init__(self, cen, r):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,0,0), 1)
            Center (0.0, 0.0, 0.0) radius 1
            sage: BoundingSphere((0,-1,5), 2)
            Center (0.0, -1.0, 5.0) radius 2
            """
        self.cen = vector(RDF, cen)
        self.r = r

    def __repr__(self):
        """
        TESTS::

            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,-1,10), 2)
            Center (0.0, -1.0, 10.0) radius 2
        """
        return "Center %s radius %s" % (self.cen, self.r)

    def __add__(self, other):
        """
        Returns the bounding sphere containing both terms.

        EXAMPLES::

            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,0,0), 1) + BoundingSphere((0,0,0), 2)
            Center (0.0, 0.0, 0.0) radius 2
            sage: BoundingSphere((0,0,0), 1) + BoundingSphere((0,0,100), 1)
            Center (0.0, 0.0, 50.0) radius 51.0
            sage: BoundingSphere((0,0,0), 1) + BoundingSphere((1,1,1), 2)
            Center (0.788675134595, 0.788675134595, 0.788675134595) radius 2.36602540378

        Treat None and 0 as the identity::

            sage: BoundingSphere((1,2,3), 10) + None + 0
            Center (1.0, 2.0, 3.0) radius 10

        """
        if other == 0 or other is None:
            return self
        elif self == 0 or self is None:
            return other
        if self.cen == other.cen:
            return self if self.r > other.r else other
        diff = other.cen - self.cen
        dist = (diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]).sqrt()
        diam = dist + self.r + other.r
        off  = diam/2 - self.r
        return BoundingSphere(self.cen + (off/dist)*diff, diam/2)

    def transform(self, T):
        """
        Returns the bounding sphere of this sphere acted on by T. This always
        returns a new sphere, even if the resulting object is an ellipsoid.

        EXAMPLES::

            sage: from sage.plot.plot3d.transform import Transformation
            sage: from sage.plot.plot3d.base import BoundingSphere
            sage: BoundingSphere((0,0,0), 10).transform(Transformation(trans=(1,2,3)))
            Center (1.0, 2.0, 3.0) radius 10.0
            sage: BoundingSphere((0,0,0), 10).transform(Transformation(scale=(1/2, 1, 2)))
            Center (0.0, 0.0, 0.0) radius 20.0
            sage: BoundingSphere((0,0,3), 10).transform(Transformation(scale=(2, 2, 2)))
            Center (0.0, 0.0, 6.0) radius 20.0
        """
        return BoundingSphere(T.transform_point(self.cen), self.r * T.max_scale())


class RenderParams(SageObject):
    """
    This class is a container for all parameters that may be needed to
    render triangulate/render an object to a certain format. It can
    contain both cumulative and global parameters.

    Of particular note is the transformation object, which holds the
    cumulative transformation from the root of the scene graph to this
    node in the tree.
    """

    _uniq_counter = 0
    randomize_counter = 0
    force_reload = False
    mesh = False
    dots = False
    antialiasing = 8

    def __init__(self, **kwds):
        """
        EXAMPLES::

            sage: params = sage.plot.plot3d.base.RenderParams(foo='x')
            sage: params.transform_list
            []
            sage: params.foo
            'x'
        """
        self.output_file = sage.misc.misc.tmp_filename()
        self.obj_vertex_offset = 1
        self.transform_list = []
        self.transform = None
        self.ds = 1
        self.crease_threshold = .8
        self.__dict__.update(kwds)
        # for jmol, some things (such as labels) must be attached to atoms
        self.atom_list = []

    def push_transform(self, T):
        """
        Push a transformation onto the stack, updating self.transform.

        EXAMPLES::

            sage: from sage.plot.plot3d.transform import Transformation
            sage: params = sage.plot.plot3d.base.RenderParams()
            sage: params.transform is None
            True
            sage: T = Transformation(scale=(10,20,30))
            sage: params.push_transform(T)
            sage: params.transform.get_matrix()
            [10.0  0.0  0.0  0.0]
            [ 0.0 20.0  0.0  0.0]
            [ 0.0  0.0 30.0  0.0]
            [ 0.0  0.0  0.0  1.0]
            sage: params.push_transform(T)  # scale again
            sage: params.transform.get_matrix()
            [100.0   0.0   0.0   0.0]
            [  0.0 400.0   0.0   0.0]
            [  0.0   0.0 900.0   0.0]
            [  0.0   0.0   0.0   1.0]
        """
        self.transform_list.append(self.transform)
        if self.transform is None:
            self.transform = T
        else:
            self.transform = self.transform * T

    def pop_transform(self):
        """
        Remove the last transformation off the stack, resetting self.transform
        to the previous value.

        EXAMPLES::

            sage: from sage.plot.plot3d.transform import Transformation
            sage: params = sage.plot.plot3d.base.RenderParams()
            sage: T = Transformation(trans=(100, 500, 0))
            sage: params.push_transform(T)
            sage: params.transform.get_matrix()
            [  1.0   0.0   0.0 100.0]
            [  0.0   1.0   0.0 500.0]
            [  0.0   0.0   1.0   0.0]
            [  0.0   0.0   0.0   1.0]
            sage: params.push_transform(Transformation(trans=(-100, 500, 200)))
            sage: params.transform.get_matrix()
            [   1.0    0.0    0.0    0.0]
            [   0.0    1.0    0.0 1000.0]
            [   0.0    0.0    1.0  200.0]
            [   0.0    0.0    0.0    1.0]
            sage: params.pop_transform()
            sage: params.transform.get_matrix()
            [  1.0   0.0   0.0 100.0]
            [  0.0   1.0   0.0 500.0]
            [  0.0   0.0   1.0   0.0]
            [  0.0   0.0   0.0   1.0]

        """
        self.transform = self.transform_list.pop()

    def unique_name(self, desc="name"):
        """
        Returns a unique identifier starting with desc.

        EXAMPLES::

            sage: params = sage.plot.plot3d.base.RenderParams()
            sage: params.unique_name()
            'name_1'
            sage: params.unique_name()
            'name_2'
            sage: params.unique_name('texture')
            'texture_3'
        """
        if self.randomize_counter:
            self._uniq_counter = randint(1,1000000)
        else:
            self._uniq_counter += 1
        return "%s_%s" % (desc, self._uniq_counter)

def flatten_list(L):
    """
    This is an optimized routine to turn a list of lists (of lists ...)
    into a single list. We generate data in a non-flat format to avoid
    multiple data copying, and then concatenate it all at the end.

    This is NOT recursive, otherwise there would be a lot of redundant
    copying (which we are trying to avoid in the first place, though at
    least it would be just the pointers).

    EXAMPLES::

        sage: from sage.plot.plot3d.base import flatten_list
        sage: flatten_list([])
        []
        sage: flatten_list([[[[]]]])
        []
        sage: flatten_list([['a', 'b'], 'c'])
        ['a', 'b', 'c']
        sage: flatten_list([['a'], [[['b'], 'c'], ['d'], [[['e', 'f', 'g']]]]])
        ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    """
    if not PyList_CheckExact(L):
        return [L]
    flat = []
    L_stack = []; L_pop = L_stack.pop
    i_stack = []; i_pop = i_stack.pop
    cdef Py_ssize_t i = 0
    while i < PyList_GET_SIZE(L) or PyList_GET_SIZE(L_stack) > 0:
        while i < PyList_GET_SIZE(L):
            tmp = <object>PyList_GET_ITEM(L, i)
            if PyList_CheckExact(tmp):
                PyList_Append(L_stack, L)
                L = tmp
                PyList_Append(i_stack, i)
                i = 0
            else:
                PyList_Append(flat, tmp)
                i += 1
        if PyList_GET_SIZE(L_stack) > 0:
            L = L_pop()
            i = i_pop()
            i += 1
    return flat


def min3(v):
    """
    Return the componentwise minimum of a list of 3-tuples.

    EXAMPLES::

        sage: from sage.plot.plot3d.base import min3, max3
        sage: min3([(-1,2,5), (-3, 4, 2)])
        (-3, 2, 2)
    """
    return tuple([min([a[i] for a in v]) for i in range(3)])

def max3(v):
    """
    Return the componentwise maximum of a list of 3-tuples.

    EXAMPLES::

        sage: from sage.plot.plot3d.base import min3, max3
        sage: max3([(-1,2,5), (-3, 4, 2)])
        (-1, 4, 5)
    """
    return tuple([max([a[i] for a in v]) for i in range(3)])

def point_list_bounding_box(v):
    """
    EXAMPLES::

        sage: from sage.plot.plot3d.base import point_list_bounding_box
        sage: point_list_bounding_box([(1,2,3),(4,5,6),(-10,0,10)])
        ((-10.0, 0.0, 3.0), (4.0, 5.0, 10.0))
        sage: point_list_bounding_box([(float('nan'), float('inf'), float('-inf')), (10,0,10)])
        ((10.0, 0.0, 10.0), (10.0, 0.0, 10.0))
    """
    cdef point_c low, high, cur
    low.x, low.y, low.z = INFINITY, INFINITY, INFINITY
    high.x, high.y, high.z = -INFINITY, -INFINITY, -INFINITY

    for P in v:
        cur.x, cur.y, cur.z = P
        point_c_update_finite_lower_bound(&low, cur)
        point_c_update_finite_upper_bound(&high, cur)
    return ((low.x, low.y, low.z), (high.x, high.y, high.z))

def optimal_aspect_ratios(ratios):
    # average the aspect ratios
    n = len(ratios)
    if n > 0:
        return [max([z[i] for z in ratios]) for i in range(3)]
    else:
        return [1.0,1.0,1.0]

def optimal_extra_kwds(v):
    """
    Given a list v of dictionaries, this function merges them such that
    later dictionaries have precedence.
    """
    if len(v) == 0:
        return {}
    a = dict(v[0])   # make a copy!
    for b in v[1:]:
        for k,w in b.iteritems():
            a[k] = w
    return a

