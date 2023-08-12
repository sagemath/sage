r"""
The Tachyon 3D Ray Tracer

Given any 3D graphics object ``M`` one can compute a raytraced
representation by typing ``M.show(viewer='tachyon')``.
For example, we draw two translucent spheres that contain a red
tube, and render the result using Tachyon.

::

    sage: S = sphere(opacity=0.8, aspect_ratio=[1,1,1])
    sage: L = line3d([(0,0,0),(2,0,0)], thickness=10, color='red')
    sage: M = S + S.translate((2,0,0)) + L
    sage: M.show(viewer='tachyon')

A number of options can be given to the
:meth:`~sage.plot.plot3d.base.Graphics3d.show` method and
correspondingly to the
:meth:`~sage.plot.plot3d.base.Graphics3d.save` method for saving
the generated image to a file::

    sage: M.show(viewer='tachyon',
    ....:    antialiasing=True, raydepth=3,
    ....:    figsize=[12,8], # the image resolution is 100*figsize
    ....:    camera_position=[4, 4.4, 1], # a distant camera position combined with
    ....:    zoom=3, # a large zoom factor will decrease perspective distortion.
    ....:    updir=(0, -0.1, 1), # the camera is slightly tilted
    ....:    viewdir=(-2.,-2.,-0.5), # slightly off-center
    ....:    light_position=(4.0, -3.0, 2.0),
    ....:   )

One can also directly control Tachyon by creating a ``Tachyon`` object
and adding elements of the scene one by one, which gives a huge amount of
flexibility. For example, here we directly use Tachyon to draw 3
spheres on the coordinate axes::

    sage: t = Tachyon(xres=500, yres=500, camera_position=(2,0,0))
    sage: t.light((4,3,2), 0.2, (1,1,1))
    sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
    sage: t.texture('t3', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
    sage: t.texture('t4', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
    sage: t.sphere((0,0.5,0), 0.2, 't2')
    sage: t.sphere((0.5,0,0), 0.2, 't3')
    sage: t.sphere((0,0,0.5), 0.2, 't4')
    sage: t.show()

For scenes with many reflections it is helpful to increase the raydepth option, and turn on antialiasing.  The following scene is an extreme case with many reflections between four cotangent spheres::

    sage: t = Tachyon(camera_position=(0,-4,1), xres=800, yres=600, raydepth=12,
    ....:             aspectratio=.75, antialiasing=4)
    sage: t.light((0.02,0.012,0.001), 0.01, (1,0,0))
    sage: t.light((0,0,10), 0.01, (0,0,1))
    sage: t.texture('s', color=(.8,1,1), opacity=.9, specular=.95, diffuse=.3, ambient=0.05)
    sage: t.texture('p', color=(0,0,1), opacity=1, specular=.2)
    sage: t.sphere((-1,-.57735,-0.7071),1,'s')
    sage: t.sphere((1,-.57735,-0.7071),1,'s')
    sage: t.sphere((0,1.15465,-0.7071),1,'s')
    sage: t.sphere((0,0,0.9259),1,'s')
    sage: t.plane((0,0,-1.9259),(0,0,1),'p')
    sage: t.show() # long time

Different projection options are available. The following examples all
use a sphere and cube::

    sage: cedges = [[[1, 1, 1], [-1, 1, 1]], [[1, 1, 1], [1, -1, 1]],
    ....: [[1, 1, 1], [1, 1, -1]], [[-1, 1, 1], [-1, -1, 1]], [[-1, 1, 1],
    ....: [-1, 1, -1]], [[1, -1, 1], [-1, -1, 1]], [[1, -1, 1], [1, -1, -1]],
    ....: [[-1, -1, 1], [-1, -1, -1]], [[1, 1, -1], [-1, 1, -1]],
    ....: [[1, 1, -1], [1, -1, -1]], [[-1, 1, -1], [-1, -1, -1]],
    ....: [[1, -1, -1], [-1, -1, -1]]]

The default projection is ``'perspective'``::

    sage: t = Tachyon(xres=800, yres=600, camera_position=(-1.5,0.0,0.0), zoom=.2)
    sage: t.texture('t1', color=(0,0,1))
    sage: for ed in cedges:
    ....:     t.fcylinder(ed[0], ed[1], .05, 't1')
    sage: t.light((-4,-4,4), .1, (1,1,1))
    sage: t.show()

Another option is ``projection='fisheye'``, which requires frustum
information. The frustum data is (bottom angle, top angle, left
angle, right angle)::

    sage: t = Tachyon(xres=800, yres=600, camera_position=(-1.5,0.0,0.0),
    ....: projection='fisheye', frustum=(-1.2, 1.2, -1.2, 1.2))
    sage: t.texture('t1', color=(0,0,1))
    sage: for ed in cedges:
    ....:     t.fcylinder(ed[0], ed[1], .05, 't1')
    sage: t.light((-4,-4,4), .1, (1,1,1))
    sage: t.show()

Finally there is the ``projection='perspective_dof'`` option. ::

    sage: T = Tachyon(xres=800, antialiasing=4, raydepth=10,
    ....:             projection='perspective_dof', focallength='1.0', aperture='.0025')
    sage: T.light((0,5,7), 1.0, (1,1,1))
    sage: T.texture('t1', opacity=1, specular=.3)
    sage: T.texture('t2', opacity=1, specular=.3, color=(0,0,1))
    sage: T.texture('t3', opacity=1, specular=1, color=(1,.8,1), diffuse=0.2)
    sage: T.plane((0,0,-1), (0,0,1), 't3')
    sage: ttlist = ['t1', 't2']
    sage: tt = 't1'
    sage: T.cylinder((0,0,.1), (1,1/3,0), .05, 't3')
    sage: for q in srange(-3, 100, .15):
    ....:     if tt == 't1':
    ....:         tt = 't2'
    ....:     else:
    ....:         tt = 't1'
    ....:     T.sphere((q, q/3+.3*sin(3*q), .1+.3*cos(3*q)), .1, tt)
    sage: T.show()

Image files in the ``ppm`` format can be used to tile planes or cover
cylinders or spheres. In this example an image is created and then
used to tile the plane::

    sage: T = Tachyon(xres=800, yres=600, camera_position=(-2.0,-.1,.3),
    ....:             projection='fisheye', frustum=(-1.0, 1.0, -1.0, 1.0))
    sage: T.texture('t1',color=(0,0,1))
    sage: for ed in cedges:
    ....:     T.fcylinder(ed[0], ed[1], .05, 't1')
    sage: T.light((-4,-4,4),.1,(1,1,1))
    sage: fname_png = tmp_filename(ext='.png')
    sage: fname_ppm = tmp_filename(ext='.ppm')
    sage: T.save(fname_png)
    sage: r2 = os.system('convert '+fname_png+' '+fname_ppm)    # optional -- ImageMagick

    sage: # optional - imagemagick
    sage: T = Tachyon(xres=800, yres=600,
    ....:             camera_position=(-2.0,-.1,.3),
    ....:             projection='fisheye',
    ....:             frustum=(-1.0, 1.0, -1.0, 1.0))
    sage: T.texture('t1', color=(1,0,0), specular=.9)
    sage: T.texture('p1', color=(1,1,1), opacity=.1,
    ....:           imagefile=fname_ppm, texfunc=9)
    sage: T.sphere((0,0,0), .5, 't1')
    sage: T.plane((0,0,-1), (0,0,1), 'p1')
    sage: T.light((-4,-4,4), .1, (1,1,1))
    sage: T.show()

AUTHOR:

- John E. Stone (johns@megapixel.com): wrote tachyon ray tracer

- William Stein: sage-tachyon interface

- Joshua Kantor: 3d function plotting

- Tom Boothby: 3d function plotting n'stuff

- Leif Hille: key idea for bugfix for texfunc issue (:trac:`799`)

- Marshall Hampton: improved doctests, rings, axis-aligned boxes.

- Paul Graham: Respect global verbosity settings (:trac:`16228`)

.. TODO::

    - clean up trianglefactory stuff
"""
from .tri_plot import Triangle, SmoothTriangle, TriangleFactory, TrianglePlot

from sage.interfaces.tachyon import tachyon_rt

from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject

from sage.misc.verbose import get_verbose
from sage.misc.temporary_file import tmp_filename

# from sage.ext import fast_tachyon_routines

from math import sqrt


class Tachyon(WithEqualityById, SageObject):
    r"""
    Create a scene the can be rendered using the Tachyon ray tracer.

    INPUT:

    - ``xres`` -- (default 350)
    - ``yres`` -- (default 350)
    - ``zoom`` -- (default 1.0)
    - ``antialiasing`` -- (default ``False``)
    - ``aspectratio``  -- (default 1.0)
    - ``raydepth`` -- (default 8)
    - ``camera_position`` -- (default (-3, 0, 0))
    - ``updir`` -- (default (0, 0, 1))
    - ``look_at`` -- (default (0,0,0))
    - ``viewdir`` -- (default ``None``), otherwise list of three numbers
    - ``projection`` -- ``'PERSPECTIVE'`` (default), ``'perspective_dof'``
      or ``'fisheye'``.
    - ``frustum`` -- (default ``''``), otherwise list of four numbers. Only
      used with ``projection='fisheye'``.
    - ``focallength`` -- (default ''), otherwise a number. Only used
      with ``projection='perspective_dof'``.
    - ``aperture`` -- (default ''), otherwise a number.  Only used
      with ``projection='perspective_dof'``.

    OUTPUT: A Tachyon 3d scene.

    Note that the coordinates are by default such that `z` is
    up, positive `y` is to the {left} and `x` is toward
    you. This is not oriented according to the right hand rule.

    EXAMPLES: Spheres along the twisted cubic.

    ::

        sage: t = Tachyon(xres=512,yres=512, camera_position=(3,0.3,0))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2, diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in srange(-1,1,0.05):
        ....:    k += 1
        ....:    t.sphere((i,i^2-0.5,i^3), 0.1, 't%s'%(k%3))
        sage: t.show()

    Another twisted cubic, but with a white background, got by putting
    infinite planes around the scene.

    ::

        sage: t = Tachyon(xres=512,yres=512, camera_position=(3,0.3,0), raydepth=8)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2,diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: t.texture('white', color=(1,1,1))
        sage: t.plane((0,0,-1), (0,0,1), 'white')
        sage: t.plane((0,-20,0), (0,1,0), 'white')
        sage: t.plane((-20,0,0), (1,0,0), 'white')

    ::

        sage: k=0
        sage: for i in srange(-1,1,0.05):
        ....:    k += 1
        ....:    t.sphere((i,i^2 - 0.5,i^3), 0.1, 't%s'%(k%3))
        ....:    t.cylinder((0,0,0), (0,0,1), 0.05,'t1')
        sage: t.show()

    Many random spheres::

        sage: t = Tachyon(xres=512,yres=512, camera_position=(2,0.5,0.5), look_at=(0.5,0.5,0.5), raydepth=4)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2, diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in range(100):
        ....:    k += 1
        ....:    t.sphere((random(),random(), random()), random()/10, 't%s'%(k%3))
        sage: t.show()

    Points on an elliptic curve, their height indicated by their height
    above the axis::

        sage: # needs sage.schemes
        sage: t = Tachyon(camera_position=(5,2,2), look_at=(0,1,0))
        sage: t.light((10,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
        sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
        sage: E = EllipticCurve('37a')                                                  # needs sage.schemes
        sage: P = E([0,0])                                                              # needs sage.schemes
        sage: Q = P                                                                     # needs sage.schemes
        sage: n = 100
        sage: for i in range(n):   # increase 20 for a better plot                      # needs sage.schemes
        ....:    Q = Q + P
        ....:    t.sphere((Q[1], Q[0], ZZ(i)/n), 0.1, 't%s'%(i%3))
        sage: t.show()                                                                  # needs sage.schemes

    A beautiful picture of rational points on a rank 1 elliptic curve.

    ::

        sage: # needs sage.schemes
        sage: t = Tachyon(xres=1000, yres=800, camera_position=(2,7,4),
        ....:             look_at=(2,0,0), raydepth=4)
        sage: t.light((10,3,2), 1, (1,1,1))
        sage: t.light((10,-3,2), 1, (1,1,1))
        sage: t.texture('black', color=(0,0,0))
        sage: t.texture('red', color=(1,0,0))
        sage: t.texture('grey', color=(.9,.9,.9))
        sage: t.plane((0,0,0),(0,0,1),'grey')
        sage: t.cylinder((0,0,0),(1,0,0),.01,'black')
        sage: t.cylinder((0,0,0),(0,1,0),.01,'black')
        sage: E = EllipticCurve('37a')                                                  # needs sage.schemes
        sage: P = E([0,0])                                                              # needs sage.schemes
        sage: Q = P                                                                     # needs sage.schemes
        sage: n = 100
        sage: for i in range(n):                                                        # needs sage.schemes
        ....:    Q = Q + P
        ....:    c = i/n + .1
        ....:    t.texture('r%s'%i,color=(float(i/n),0,0))
        ....:    t.sphere((Q[0], -Q[1], .01), .04, 'r%s'%i)
        sage: t.show()                          # long time                             # needs sage.schemes

    A beautiful spiral.

    ::

        sage: t = Tachyon(xres=800, yres=800, camera_position=(2,5,2), look_at=(2.5,0,0))
        sage: t.light((0,0,100), 1, (1,1,1))
        sage: t.texture('r', ambient=0.1, diffuse=0.9, specular=0.5,
        ....:           opacity=1.0, color=(1,0,0))
        sage: for i in srange(0,50,0.1):
        ....:    t.sphere((i/10.0,sin(i),cos(i)), 0.05, 'r')
        sage: t.texture('white', color=(1,1,1), opacity=1, specular=1, diffuse=1)
        sage: t.plane((0,0,-100), (0,0,-100), 'white')
        sage: t.show()

    If the optional parameter ``viewdir`` is not set, the camera
    center should not coincide with the point which
    is looked at (see :trac:`7232`)::

        sage: t = Tachyon(xres=80,yres=80, camera_position=(2,5,2), look_at=(2,5,2))
        Traceback (most recent call last):
        ...
        ValueError: camera_position and look_at coincide

    Use of a fisheye lens perspective. ::

        sage: T = Tachyon(xres=800, yres=600, camera_position=(-1.5,-1.5,.3),
        ....:             projection='fisheye', frustum=(-1.0, 1.0, -1.0, 1.0))
        sage: T.texture('t1', color=(0,0,1))
        sage: cedges = [[[1, 1, 1], [-1, 1, 1]], [[1, 1, 1], [1, -1, 1]],
        ....:           [[1, 1, 1], [1, 1, -1]], [[-1, 1, 1], [-1, -1, 1]],
        ....:           [[-1, 1, 1], [-1, 1, -1]], [[1, -1, 1], [-1, -1, 1]],
        ....:           [[1, -1, 1], [1, -1, -1]], [[-1, -1, 1], [-1, -1, -1]],
        ....:           [[1, 1, -1], [-1, 1, -1]], [[1, 1, -1], [1, -1, -1]],
        ....:           [[-1, 1, -1], [-1, -1, -1]], [[1, -1, -1], [-1, -1, -1]]]
        sage: for ed in cedges:
        ....:     T.fcylinder(ed[0], ed[1], .05, 't1')
        sage: T.light((-4,-4,4), .1, (1,1,1))
        sage: T.show()

    Use of the ``projection='perspective_dof'`` option.  This may not be
    implemented correctly. ::

        sage: T = Tachyon(xres=800, antialiasing=4, raydepth=10,
        ....:             projection='perspective_dof', focallength='1.0', aperture='.0025')
        sage: T.light((0,5,7), 1.0, (1,1,1))
        sage: T.texture('t1', opacity=1, specular=.3)
        sage: T.texture('t2', opacity=1, specular=.3, color=(0,0,1))
        sage: T.texture('t3', opacity=1, specular=1, color=(1,.8,1), diffuse=0.2)
        sage: T.plane((0,0,-1), (0,0,1), 't3')
        sage: ttlist = ['t1', 't2']
        sage: tt = 't1'
        sage: T.cylinder((0,0,.1), (1,1/3,0), .05, 't3')
        sage: for q in srange(-3, 100, .15):
        ....:     if tt == 't1':
        ....:         tt = 't2'
        ....:     else:
        ....:         tt = 't1'
        ....:     T.sphere((q, q/3+.3*sin(3*q), .1+.3*cos(3*q)), .1, tt)
        sage: T.show()

    TESTS::

        sage: hash(Tachyon()) # random
        140658972348064
    """
    def __init__(self,
                 xres=350, yres=350,
                 zoom=1.0,
                 antialiasing=False,
                 aspectratio=1.0,
                 raydepth=8,
                 camera_position=None, # default value (-3, 0, 0),
                 camera_center=None, # alternative equivalent name
                 updir=[0, 0, 1],
                 look_at=[0, 0, 0],
                 viewdir=None,
                 projection='PERSPECTIVE',
                 focallength='',
                 aperture='',
                 frustum=''):
        r"""
        Create an instance of the Tachyon class.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t._xres
            350
        """
        self._xres = xres
        self._yres = yres
        self._zoom = zoom
        self._aspectratio = aspectratio
        self._antialiasing = antialiasing
        self._raydepth = raydepth
        if camera_position is not None:
            self._camera_position = camera_position
        elif camera_center is not None: # make sure that old programs continue to work
            self._camera_position = camera_center
        else:
            self._camera_position = (-3, 0, 0) # default value
        self._updir = updir
        self._projection = projection
        self._focallength = focallength
        self._aperture = aperture
        self._frustum = frustum
        self._objects = []
        if viewdir is None:
            if look_at != self._camera_position:
                self._viewdir = [look_at[i] - self._camera_position[i]
                                 for i in range(3)]
            else:
                raise ValueError('camera_position and look_at coincide')
        else:
            self._viewdir = viewdir

    def save_image(self, filename=None, *args, **kwds):
        r"""
        Save an image representation of ``self``.

        The image type is
        determined by the extension of the filename.  For example,
        this could be ``.png``, ``.jpg``, ``.gif``, ``.pdf``,
        ``.svg``.  Currently this is implemented by calling the
        :meth:`save` method of self, passing along all arguments and
        keywords.

        .. NOTE::

            Not all image types are necessarily implemented for all
            graphics types.  See :meth:`save` for more details.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,11), 1,(1,1,1))
            sage: q.texture('s')
            sage: q.sphere((0,-1,1),1,'s')
            sage: tempname = tmp_filename()
            sage: q.save_image(tempname)

        TESTS:

        :meth:`save_image` is used for generating animations::

            sage: def tw_cubic(t):
            ....:     q = Tachyon()
            ....:     q.light((1,1,11), 1,(1,1,1))
            ....:     q.texture('s')
            ....:     for i in srange(-1,t,0.05):
            ....:         q.sphere((i,i^2-0.5,i^3), 0.1, 's')
            ....:     return q

            sage: a = animate([tw_cubic(t) for t in srange(-1,1,.3)])
            sage: a        # optional -- ImageMagick
            Animation with 7 frames
            sage: a.show() # optional -- ImageMagick
        """
        self.save(filename, *args, **kwds)

    def save(self, filename='sage.png', verbose=None, extra_opts=''):
        r"""
        Save rendering of the tachyon scene

        INPUT:

        -  ``filename`` - (default: 'sage.png') output
           filename; the extension of the filename determines the type.
           Supported types include:

        -  ``tga`` - 24-bit (uncompressed)

        -  ``bmp`` - 24-bit Windows BMP (uncompressed)

        -  ``ppm`` - 24-bit PPM (uncompressed)

        -  ``rgb`` - 24-bit SGI RGB (uncompressed)

        -  ``png`` - 24-bit PNG (compressed, lossless)

        -  ``verbose`` - integer (default: ``None``); if no verbosity setting
           is supplied, the verbosity level set by
           ``sage.misc.verbose.set_verbose`` is used.

        -  ``0`` - silent

        -  ``1`` - some output

        -  ``2`` - very verbose output

        -  ``extra_opts`` - passed directly to tachyon command
           line. Use tachyon_rt.usage() to see some of the possibilities.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,11), 1,(1,1,1))
            sage: q.texture('s')
            sage: q.sphere((0,0,0),1,'s')
            sage: tempname = tmp_filename()
            sage: q.save(tempname)
        """
        if verbose is None:
            verbose = get_verbose()
        tachyon_rt(self.str(), filename, verbose, extra_opts)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,11), 1,(1,1,1))
            sage: q.texture('s')
            sage: q.sphere((0,0,0),1,'s')
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: q._rich_repr_(dm)
            OutputImagePng container
        """
        OutputImagePng = display_manager.types.OutputImagePng
        if OutputImagePng not in display_manager.supported_output():
            return
        filename = tmp_filename(ext='.png')
        self.save(filename, **kwds)
        from sage.repl.rich_output.buffer import OutputBuffer
        buf = OutputBuffer.from_file(filename)
        return OutputImagePng(buf)

    def show(self, **kwds):
        r"""
        Create a PNG file of the scene.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        EXAMPLES:

        This example demonstrates how the global Sage verbosity setting
        is used if none is supplied. Firstly, using a global verbosity
        setting of 0 means no extra technical information is displayed,
        and we are simply shown the plot.

        ::

            sage: h = Tachyon(xres=512, yres=512, camera_position=(4,-4,3),
            ....:             viewdir=(-4,4,-3), raydepth=4)
            sage: h.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return float(sin(x*y))
            sage: h.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,
            ....:           opacity=1.0, color=(1.0,0,0))
            sage: h.plot(f, (-4,4), (-4,4), "t0", max_depth=5, initial_depth=3,         # needs sage.symbolic
            ....:        num_colors=60)  # increase min_depth for better picture
            sage: from sage.misc.verbose import set_verbose, get_verbose
            sage: set_verbose(0)
            sage: h.show()                                                              # needs sage.symbolic

        This second example, using a "medium" global verbosity
        setting of 1, displays some extra technical information then
        displays our graph.

        ::

            sage: s = Tachyon(xres=512, yres=512, camera_position=(4,-4,3),
            ....:             viewdir=(-4,4,-3), raydepth=4)
            sage: s.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return float(sin(x*y))
            sage: s.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,
            ....:           opacity=1.0, color=(1.0,0,0))
            sage: s.plot(f, (-4,4), (-4,4), "t0", max_depth=5, initial_depth=3,         # needs sage.symbolic
            ....:        num_colors=60)  # increase min_depth for better picture
            sage: set_verbose(1)
            sage: s.show()                                                              # needs sage.symbolic
            tachyon ...
            Scene contains 2713 objects.
            ...

        The last example shows how you can override the global Sage
        verbosity setting, my supplying a setting level as an argument.
        In this case we chose the highest verbosity setting level, 2,
        so much more extra technical information is shown, along with
        the plot.

        ::

            sage: set_verbose(0)
            sage: d = Tachyon(xres=512, yres=512, camera_position=(4,-4,3),
            ....:             viewdir=(-4,4,-3), raydepth=4)
            sage: d.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return float(sin(x*y))
            sage: d.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,
            ....:           opacity=1.0, color=(1.0,0,0))
            sage: d.plot(f,(-4,4),(-4,4),"t0",max_depth=5,initial_depth=3,              # needs sage.symbolic
            ....:        num_colors=60)  # increase min_depth for better picture
            sage: get_verbose()
            0
            sage: d.show(verbose=2)                                                     # needs sage.symbolic
            tachyon ...
            Scene contains 2713 objects.
            ...
            Scene contains 1 non-gridded objects
            ...
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, **kwds)

    def _res(self):
        r"""
        An internal function that writes the tachyon string for the
        resolution (x and y size of the image).

        EXAMPLES::

            sage: t = Tachyon(xres=300, yres=700)
            sage: t._res()
            '\nresolution 300 700\n'
        """
        return '\nresolution %s %s\n' % (self._xres, self._yres)

    def _camera(self):
        r"""
        An internal function that writes the tachyon string for the
        camera and other rendering information (ray depth, antialiasing).

        EXAMPLES::

            sage: t = Tachyon(raydepth=16, zoom=2, antialiasing=True)
            sage: t._camera().split()[3:10]
            ['zoom', '2.0', 'aspectratio', '1.0', 'antialiasing', '1', 'raydepth']
        """
        camera_out = r"""
           camera
              projection %s""" % (tostr(self._projection))
        if self._focallength != '':
            camera_out = camera_out + r"""
              focallength %s""" % (float(self._focallength))
        if self._aperture != '':
            camera_out = camera_out + r"""
              aperture %s""" % (float(self._aperture))
        camera_out = camera_out + r"""
              zoom %s
              aspectratio %s
              antialiasing %s
              raydepth %s
              center %s
              viewdir %s
              updir %s""" % (float(self._zoom),
                             float(self._aspectratio),
                             int(self._antialiasing),
                             int(self._raydepth),
                             tostr(self._camera_position),
                             tostr(self._viewdir),
                             tostr(self._updir))
        if self._frustum != '':
            camera_out = camera_out + r"""
              frustum %s""" % (tostr(self._frustum))
        camera_out = camera_out + r"""
           end_camera"""
        return camera_out

    def str(self):
        r"""
        Return the complete tachyon scene file as a string.

        EXAMPLES::

            sage: t = Tachyon(xres=500,yres=500, camera_position=(2,0,0))
            sage: t.light((4,3,2), 0.2, (1,1,1))
            sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
            sage: t.texture('t3', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
            sage: t.texture('t4', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
            sage: t.sphere((0,0.5,0), 0.2, 't2')
            sage: t.sphere((0.5,0,0), 0.2, 't3')
            sage: t.sphere((0,0,0.5), 0.2, 't4')
            sage: 'PLASTIC' in t.str()
            True
        """
        return r"""
        begin_scene
        %s
        %s
        %s
        end_scene""" % (self._res(),
                        self._camera(),
                        '\n'.join(x.str() for x in self._objects))

    def light(self, center, radius, color):
        r"""
        Create a light source of the given center, radius, and color.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,1),1.0,(.2,0,.8))
            sage: q.str().split('\n')[17]
            '        light center  1.0 1.0 1.0 '
        """
        self._objects.append(Light(center, radius, color))

    def texfunc(self, type=0, center=(0, 0, 0), rotate=(0, 0, 0),
                scale=(1, 1, 1),
                imagefile=''):
        r"""
        INPUT:

        -  ``type`` - (default: 0)

           0. No special texture, plain shading
           1. 3D checkerboard function, like a rubik's cube
           2. Grit Texture, randomized surface color
           3. 3D marble texture, uses object's base color
           4. 3D wood texture, light and dark brown, not very good yet
           5. 3D gradient noise function (can't remember what it looks
              like)
           6. Don't remember
           7. Cylindrical Image Map, requires ppm filename (with path)
           8. Spherical Image Map, requires ppm filename (with path)
           9. Planar Image Map, requires ppm filename (with path)

        -  ``center`` - (default: (0,0,0))
        -  ``rotate`` - (default: (0,0,0))
        -  ``scale`` - (default: (1,1,1))


        EXAMPLES: We draw an infinite checkerboard::

            sage: t = Tachyon(camera_position=(2,7,4), look_at=(2,0,0))
            sage: t.texture('black', color=(0,0,0), texfunc=1)
            sage: t.plane((0,0,0),(0,0,1),'black')
            sage: t.show()
        """
        type = int(type)
        if type < 0 or type > 9:
            raise ValueError("type must be an integer between 0 and 9")
        return Texfunc(type, center, rotate, scale, imagefile=imagefile).str()

    def texture(self, name, ambient=0.2, diffuse=0.8,
                specular=0.0, opacity=1.0,
                color=(1.0, 0.0, 0.5), texfunc=0, phong=0, phongsize=.5,
                phongtype="PLASTIC", imagefile=''):
        r"""
        INPUT:

        -  ``name`` - string; the name of the texture (to be
           used later)

        -  ``ambient`` - (default: 0.2)

        -  ``diffuse`` - (default: 0.8)

        -  ``specular`` - (default: 0.0)

        -  ``opacity`` - (default: 1.0)

        -  ``color`` - (default: (1.0,0.0,0.5))

        -  ``texfunc`` - (default: 0); a texture function; this
           is either the output of self.texfunc, or a number between 0 and 9,
           inclusive. See the docs for self.texfunc.

        -  ``phong`` - (default: 0)

        -  ``phongsize`` - (default: 0.5)

        -  ``phongtype`` - (default: "PLASTIC")

        EXAMPLES:

        We draw a scene with 4 spheres that illustrates various uses of
        the texture command::

            sage: t = Tachyon(camera_position=(2,5,4), look_at=(2,0,0), raydepth=6)
            sage: t.light((10,3,4), 1, (1,1,1))
            sage: t.texture('mirror', ambient=0.05, diffuse=0.05, specular=.9,
            ....:           opacity=0.9, color=(.8,.8,.8))
            sage: t.texture('grey', color=(.8,.8,.8), texfunc=3)
            sage: t.plane((0,0,0),(0,0,1),'grey')
            sage: t.sphere((4,-1,1), 1, 'mirror')
            sage: t.sphere((0,-1,1), 1, 'mirror')
            sage: t.sphere((2,-1,1), 0.5, 'mirror')
            sage: t.sphere((2,1,1), 0.5, 'mirror')
            sage: show(t)  # known bug (trac #7232)
        """
        if texfunc and not isinstance(texfunc, Texfunc):
            texfunc = self.texfunc(int(texfunc), imagefile=imagefile)
        self._objects.append(Texture(name, ambient, diffuse,
                                     specular, opacity, color, texfunc,
                                     phong, phongsize, phongtype,
                                     imagefile=imagefile))

    def texture_recolor(self, name, colors):
        r"""
        Recolor default textures.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('s')
            sage: q = t.texture_recolor('s',[(0,0,1)])
            sage: t._objects[1]._color
            (0.0, 0.0, 1.0)
        """
        base_tex = None
        names = []
        ident = "SAGETEX%d" % len(self._objects)  # don't collide with other texture names

        for o in self._objects:
            if isinstance(o, Texture) and o._name == name:
                base_tex = o
                break
        if base_tex is None:
            base_tex = Texture(name)

        for i in range(len(colors)):
            n = "%s_%d" % (ident, i)
            self._objects.append(base_tex.recolor(n, colors[i]))
            names.append(n)

        return names

    def sphere(self, center, radius, texture):
        r"""
        Create the scene information for a sphere with the given
        center, radius, and texture.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('sphere_texture')
            sage: t.sphere((1,2,3), .1, 'sphere_texture')
            sage: t._objects[1].str()
            '\n        sphere center  1.0 2.0 3.0  rad 0.1 sphere_texture\n        '
        """
        self._objects.append(Sphere(center, radius, texture))

    def ring(self, center, normal, inner, outer, texture):
        r"""
        Create the scene information for a ring with the given parameters.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.ring([0,0,0], [0,0,1], 1.0, 2.0, 's')
            sage: t._objects[0]._center
            (0.0, 0.0, 0.0)
        """
        self._objects.append(Ring(center, normal, inner, outer, texture))

    def cylinder(self, center, axis, radius, texture):
        r"""
        Create the scene information for a infinite cylinder with the
        given center, axis direction, radius, and texture.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('c')
            sage: t.cylinder((0,0,0),(-1,-1,-1),.1,'c')
        """
        self._objects.append(Cylinder(center, axis, radius, texture))

    def plane(self, center, normal, texture):
        r"""
        Create an infinite plane with the given center and normal.

        TESTS::

            sage: t = Tachyon()
            sage: t.plane((0,0,0),(1,1,1),'s')
            sage: plane_pos = t.str().index('plane')
            sage: t.str()[plane_pos:plane_pos+42]
            'plane center  0.0 0.0 0.0  normal  1.0 1.0'
        """
        self._objects.append(Plane(center, normal, texture))

    def axis_aligned_box(self, min_p, max_p, texture):
        r"""
        Create an axis-aligned box with minimal point ``min_p`` and
        maximum point ``max_p``.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.axis_aligned_box((0,0,0),(2,2,2),'s')
        """
        self._objects.append(Axis_aligned_box(min_p, max_p, texture))

    def fcylinder(self, base, apex, radius, texture):
        r"""
        Finite cylinders are almost the same as infinite ones, but the
        center and length of the axis determine the extents of the
        cylinder.

        The finite cylinder is also really a shell, it
        does not have any caps. If you need to close off the ends of
        the cylinder, use two ring objects, with the inner radius set
        to 0.0 and the normal set to be the axis of the cylinder.
        Finite cylinders are built this way to enhance speed.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.fcylinder((1,1,1),(1,2,3),.01,'s')
            sage: len(t.str())
            451
        """
        self._objects.append(FCylinder(base, apex, radius, texture))

    def triangle(self, vertex_1, vertex_2, vertex_3, texture):
        r"""
        Create a triangle with the given vertices and texture.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('s')
            sage: t.triangle([1,2,3],[4,5,6],[7,8,10],'s')
            sage: t._objects[1].get_vertices()
            ([1, 2, 3], [4, 5, 6], [7, 8, 10])

        """
        self._objects.append(TachyonTriangle(vertex_1, vertex_2, vertex_3,
                                             texture))

    def smooth_triangle(self, vertex_1, vertex_2, vertex_3, normal_1, normal_2, normal_3, texture):
        r"""
        Create a triangle along with a normal vector for smoothing.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.light((1,1,1),.1,(1,1,1))
            sage: t.texture('s')
            sage: t.smooth_triangle([0,0,0],[0,0,1],[0,1,0],[0,1,1],[-1,1,2],[3,0,0],'s')
            sage: t._objects[2].get_vertices()
            ([0, 0, 0], [0, 0, 1], [0, 1, 0])
            sage: t._objects[2].get_normals()
            ([0, 1, 1], [-1, 1, 2], [3, 0, 0])
        """
        self._objects.append(TachyonSmoothTriangle(vertex_1, vertex_2, vertex_3, normal_1, normal_2, normal_3, texture))

    def fractal_landscape(self, res, scale, center, texture):
        r"""
        Axis-aligned fractal landscape.

        Not very useful at the moment.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('s')
            sage: t.fractal_landscape([30,30],[80,80],[0,0,0],'s')
            sage: len(t._objects)
            2
        """
        self._objects.append(FractalLandscape(res, scale, center, texture))

    def plot(self, f, xmin_xmax, ymin_ymax, texture, grad_f=None,
             max_bend=.7, max_depth=5, initial_depth=3, num_colors=None):
        r"""
        INPUT:

        -  ``f`` - Function of two variables, which returns a
           float (or coercible to a float) (xmin,xmax)

        -  ``(ymin,ymax)`` - defines the rectangle to plot over
           texture: Name of texture to be used Optional arguments:

        -  ``grad_f`` - gradient function. If specified,
           smooth triangles will be used.

        -  ``max_bend`` - Cosine of the threshold angle
           between triangles used to determine whether or not to recurse after
           the minimum depth

        -  ``max_depth`` - maximum recursion depth. Maximum
           triangles plotted = `2^{2*max_depth}`

        -  ``initial_depth`` - minimum recursion depth. No
           error-tolerance checking is performed below this depth. Minimum
           triangles plotted: `2^{2*min_depth}`

        -  ``num_colors`` - Number of rainbow bands to color
           the plot with. Texture supplied will be cloned (with different
           colors) using the texture_recolor method of the Tachyon object.


        Plots a function by constructing a mesh with nonstandard sampling
        density without gaps. At very high resolutions (depths 10) it
        becomes very slow. Cython may help. Complexity is approx.
        `O(2^{2*maxdepth})`. This algorithm has been optimized for
        speed, not memory - values from f(x,y) are recycled rather than
        calling the function multiple times. At high recursion depth, this
        may cause problems for some machines.

        Flat Triangles::

            sage: t = Tachyon(xres=512, yres=512, camera_position=(4,-4,3),
            ....:             viewdir=(-4,4,-3), raydepth=4)
            sage: t.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return float(sin(x*y))
            sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,
            ....:           opacity=1.0, color=(1.0,0,0))
            sage: t.plot(f, (-4,4), (-4,4), "t0", max_depth=5, initial_depth=3,         # needs sage.symbolic
            ....:        num_colors=60)  # increase min_depth for better picture
            sage: t.show(verbose=1)                                                     # needs sage.symbolic
            tachyon ...
            Scene contains 2713 objects.
            ...

        Plotting with Smooth Triangles (requires explicit gradient
        function)::

            sage: t = Tachyon(xres=512, yres=512, camera_position=(4,-4,3),
            ....:             viewdir=(-4,4,-3), raydepth=4)
            sage: t.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return float(sin(x*y))
            sage: def g(x,y): return (float(y*cos(x*y)), float(x*cos(x*y)), 1)
            sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,
            ....:           opacity=1.0, color=(1.0,0,0))
            sage: t.plot(f, (-4,4), (-4,4), "t0", max_depth=5, initial_depth=3,         # needs sage.symbolic
            ....:        grad_f=g)  # increase min_depth for better picture
            sage: t.show(verbose=1)                                                     # needs sage.symbolic
            tachyon ...
            Scene contains 2713 objects.
            ...

        Preconditions: f is a scalar function of two variables, grad_f is
        None or a triple-valued function of two variables, min_x !=
        max_x, min_y != max_y

        ::

            sage: f = lambda x,y: x*y
            sage: t = Tachyon()
            sage: t.plot(f,(2.,2.),(-2.,2.),'')
            Traceback (most recent call last):
            ...
            ValueError: plot rectangle is really a line; make sure min_x != max_x and min_y != max_y
        """
        (xmin, xmax) = xmin_xmax
        (ymin, ymax) = ymin_ymax
        factory = TachyonTriangleFactory(self, texture)
        plot = TrianglePlot(factory, f, (xmin, xmax), (ymin, ymax), g=grad_f,
                            min_depth=initial_depth, max_depth=max_depth,
                            max_bend=max_bend, num_colors=num_colors)
        self._objects.append(plot)

    def parametric_plot(self, f, t_0, t_f, tex, r=.1, cylinders=True,
                        min_depth=4, max_depth=8, e_rel=.01, e_abs=.01):
        r"""
        Plot a space curve as a series of spheres and finite cylinders.

        Example (twisted cubic) ::

            sage: f = lambda t: (t,t^2,t^3)
            sage: t = Tachyon(camera_position=(5,0,4))
            sage: t.texture('t')
            sage: t.light((-20,-20,40), 0.2, (1,1,1))
            sage: t.parametric_plot(f,-5,5,'t',min_depth=6)
            sage: t.show(verbose=1)
            tachyon ...
            Scene contains 482 objects.
            ...
        """
        self._objects.append(
            ParametricPlot(f, t_0, t_f, tex, r=r, cylinders=cylinders,
                           min_depth=min_depth, max_depth=max_depth,
                           e_rel=.01, e_abs=.01))


class Light():
    r"""
    Represent lighting objects.

    EXAMPLES::

        sage: from sage.plot.plot3d.tachyon import Light
        sage: q = Light((1,1,1), 1, (1,1,1))
        sage: q._center
        (1.0, 1.0, 1.0)
    """
    def __init__(self, center, radius, color):
        r"""
        Store the center, radius and color.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Light
            sage: q = Light((1,1,1), 1, (1,1,1))
            sage: print(q._center, q._color, q._radius)
            (1.0, 1.0, 1.0) (1.0, 1.0, 1.0) 1.0
        """
        x, y, z = center
        self._center = (float(x), float(y), float(z))
        self._radius = float(radius)
        r, g, b = color
        self._color = (float(r), float(g), float(b))

    def str(self):
        r"""
        Return the tachyon string defining the light source.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Light
            sage: q = Light((1,1,1), 1, (1,1,1))
            sage: print(q.str())
                    light center  1.0 1.0 1.0
                          rad 1.0
                          color  1.0 1.0 1.0

        """
        return r"""
        light center %s
              rad %s
              color %s
        """ % (tostr(self._center), self._radius,
               tostr(self._color))


class Texfunc():

    def __init__(self, ttype=0, center=(0, 0, 0), rotate=(0, 0, 0),
                 scale=(1, 1, 1), imagefile=''):
        r"""
        Create a texture function.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texfunc
            sage: t = Texfunc()
            sage: t._ttype
            0
        """
        self._ttype = ttype
        x, y, z = center
        self._center = (float(x), float(y), float(z))
        x, y, z = rotate
        self._rotate = (float(x), float(y), float(z))
        x, y, z = scale
        self._scale = (float(x), float(y), float(z))
        self._imagefile = imagefile

    def str(self):
        r"""
        Return the scene string for this texture function.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texfunc
            sage: t = Texfunc()
            sage: t.str()
            '0'
        """
        if self._ttype == 0:
            return "0"
        elif self._ttype < 7 and self._ttype > 0:
            return r"""%d center %s rotate %s scale %s""" % (
                self._ttype,
                tostr(self._center),
                tostr(self._rotate),
                tostr(self._scale))
        elif self._ttype < 9:
            return r"""%d %s center %s rotate %s scale %s""" % (
                self._ttype,
                self._imagefile,
                tostr(self._center),
                tostr(self._rotate),
                tostr(self._scale))
        elif self._ttype == 9:
            return r"""%d %s center %s rotate %s scale %s
            uaxis 1.0 0.0 0.0
            vaxis 0.0 1.0 0.0""" % (
                self._ttype,
                self._imagefile,
                tostr(self._center),
                tostr(self._rotate),
                tostr(self._scale))
        else:
            raise ValueError


class Texture():

    def __init__(self, name, ambient=0.2, diffuse=0.8,
                 specular=0.0, opacity=1.0,
                 color=(1.0, 0.0, 0.5), texfunc=0,
                 phong=0, phongsize=0, phongtype="PLASTIC", imagefile=''):
        r"""
        Store texture information.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texture
            sage: t = Texture('w')
            sage: t.str().split()[2:6]
            ['ambient', '0.2', 'diffuse', '0.8']
        """
        self._name = str(name)
        self._ambient = float(ambient)
        self._diffuse = float(diffuse)
        self._specular = float(specular)
        self._opacity = float(opacity)
        r, g, b = color
        self._color = (float(r), float(g), float(b))
        self._texfunc = texfunc
        self._phong = float(phong)
        self._phongsize = float(phongsize)
        self._phongtype = str(phongtype)
        self._imagefile = str(imagefile)

    def recolor(self, name, color):
        r"""
        Return a texture with the new given color.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texture
            sage: t2 = Texture('w')
            sage: t2w = t2.recolor('w2', (.1,.2,.3))
            sage: t2ws = t2w.str()
            sage: color_index = t2ws.find('color')
            sage: t2ws[color_index:color_index+20]
            'color  0.1 0.2 0.3  '
        """
        return Texture(name, self._ambient, self._diffuse, self._specular,
                       self._opacity,
                       color, self._texfunc, self._phong, self._phongsize,
                       self._phongtype, self._imagefile)

    def str(self):
        r"""
        Return the scene string for this texture.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texture
            sage: t = Texture('w')
            sage: t.str().split()[2:6]
            ['ambient', '0.2', 'diffuse', '0.8']
        """
        return r"""
        texdef %s ambient %s diffuse %s specular %s opacity %s
        phong %s %s phong_size %s
        color %s texfunc %s
        """ % (self._name,
               self._ambient,
               self._diffuse,
               self._specular,
               self._opacity,
               self._phongtype,
               self._phong,
               self._phongsize,
               tostr(self._color),
               self._texfunc)


class Sphere():
    r"""
    A class for creating spheres in tachyon.
    """
    def __init__(self, center, radius, texture):
        r"""
        Store the center, radius, and texture information in a class.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Sphere
            sage: t.texture('r', color=(.8,0,0), ambient=.1)
            sage: s = Sphere((1,1,1), 1, 'r')
            sage: s._radius
            1.0
        """
        x, y, z = center
        self._center = (float(x), float(y), float(z))
        self._radius = float(radius)
        self._texture = texture

    def str(self):
        r"""
        Return the scene string for the sphere.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Sphere
            sage: t.texture('r', color=(.8,0,0), ambient = .1)
            sage: s = Sphere((1,1,1), 1, 'r')
            sage: s.str()
            '\n        sphere center  1.0 1.0 1.0  rad 1.0 r\n        '
        """
        return r"""
        sphere center %s rad %s %s
        """ % (tostr(self._center), self._radius, self._texture)


class Ring():
    r"""
    An annulus of zero thickness.
    """
    def __init__(self, center, normal, inner, outer, texture):
        r"""
        Create a ring with the given center, normal, inner radius,
        outer radius, and texture.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Ring
            sage: r = Ring((1,1,1), (1,1,0), 1.0, 2.0, 's')
            sage: r._center
            (1.0, 1.0, 1.0)
        """
        x, y, z = center
        self._center = (float(x), float(y), float(z))
        x, y, z = normal
        self._normal = (float(x), float(y), float(z))
        self._inner = float(inner)
        self._outer = float(outer)
        self._texture = texture

    def str(self):
        r"""
        Return the scene string of the ring.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Ring
            sage: r = Ring((0,0,0), (1,1,0), 1.0, 2.0, 's')
            sage: r.str()
            '\n        ring center  0.0 0.0 0.0  normal  1.0 1.0 0.0  inner 1.0 outer 2.0 s\n        '
        """
        return r"""
        ring center %s normal %s inner %s outer %s %s
        """ % (tostr(self._center), tostr(self._normal),
               self._inner, self._outer, self._texture)


class FractalLandscape():
    r"""
    Axis-aligned fractal landscape.

    Does not seem very useful at the moment, but perhaps will be improved in the future.
    """
    def __init__(self, res, scale, center, texture):
        r"""
        Create a fractal landscape in tachyon.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FractalLandscape
            sage: fl = FractalLandscape([20,20],[30,30],[1,2,3],'s')
            sage: fl._center
            (1.0, 2.0, 3.0)
        """
        x, y = res
        self._res = (int(x), int(y))
        x, y = scale
        self._scale = (int(x), int(y))
        x, y, z = center
        self._center = (float(x), float(y), float(z))
        self._texture = texture

    def str(self):
        r"""
        Return the scene string of the fractal landscape.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FractalLandscape
            sage: fl = FractalLandscape([20,20],[30,30],[1,2,3],'s')
            sage: fl.str()
            '\n        scape res  20 20  scale  30 30  center  1.0 2.0 3.0  s\n        '
        """
        return r"""
        scape res %s scale %s center %s %s
        """ % (tostr(self._res, 2, int), tostr(self._scale, 2, int),
               tostr(self._center), self._texture)


class Cylinder():
    r"""
    An infinite cylinder.
    """
    def __init__(self, center, axis, radius, texture):
        r"""
        Create a cylinder with the given parameters.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Cylinder
            sage: c = Cylinder((0,0,0),(1,1,1),.1,'s')
            sage: c.str()
            '\n        cylinder center  0.0 0.0 0.0  axis  1.0 1.0 1.0  rad 0.1 s\n        '
        """
        x, y, z = center
        self._center = (float(x), float(y), float(z))
        x, y, z = axis
        self._axis = (float(x), float(y), float(z))
        self._radius = float(radius)
        self._texture = texture

    def str(self):
        r"""
        Return the scene string of the cylinder.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Cylinder
            sage: c = Cylinder((0,0,0),(1,1,1),.1,'s')
            sage: c.str()
            '\n        cylinder center  0.0 0.0 0.0  axis  1.0 1.0 1.0  rad 0.1 s\n        '
            """
        return r"""
        cylinder center %s axis %s rad %s %s
        """ % (tostr(self._center), tostr(self._axis), self._radius, self._texture)


class Plane():
    r"""
    An infinite plane.
    """
    def __init__(self, center, normal, texture):
        r"""
        Create the plane object.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Plane
            sage: p = Plane((1,2,3), (1,2,4), 's')
            sage: p.str()
            '\n        plane center  1.0 2.0 3.0  normal  1.0 2.0 4.0  s\n        '
        """
        x, y, z = center
        self._center = (float(x), float(y), float(z))
        x, y, z = normal
        self._normal = (float(x), float(y), float(z))
        self._texture = texture

    def str(self):
        r"""
        Return the scene string of the plane.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Plane
            sage: p = Plane((1,2,3),(1,2,4),'s')
            sage: p.str()
            '\n        plane center  1.0 2.0 3.0  normal  1.0 2.0 4.0  s\n        '
        """
        return r"""
        plane center %s normal %s %s
        """ % (tostr(self._center), tostr(self._normal), self._texture)


class FCylinder():
    r"""
    A finite cylinder.
    """
    def __init__(self, base, apex, radius, texture):
        r"""
        Create a finite cylinder object.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FCylinder
            sage: fc = FCylinder((0,0,0),(1,1,1),.1,'s')
            sage: fc.str()
            '\n        fcylinder base  0.0 0.0 0.0  apex  1.0 1.0 1.0  rad 0.1 s\n        '
        """
        x, y, z = base
        self._center = (float(x), float(y), float(z))
        x, y, z = apex
        self._axis = (float(x), float(y), float(z))
        self._radius = float(radius)
        self._texture = texture

    def str(self):
        r"""
        Return the scene string of the finite cylinder.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FCylinder
            sage: fc = FCylinder((0,0,0),(1,1,1),.1,'s')
            sage: fc.str()
            '\n        fcylinder base  0.0 0.0 0.0  apex  1.0 1.0 1.0  rad 0.1 s\n        '
        """
        return r"""
        fcylinder base %s apex %s rad %s %s
        """ % (tostr(self._center), tostr(self._axis), self._radius, self._texture)


class Axis_aligned_box():
    r"""
    Box with axis-aligned edges with the given min and max coordinates.
    """
    def __init__(self, min_p, max_p, texture):
        r"""
        Create the axis-aligned box object.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Axis_aligned_box
            sage: aab = Axis_aligned_box((0,0,0),(1,1,1),'s')
            sage: aab.str()
            '\n        box min  0.0 0.0 0.0  max  1.0 1.0 1.0  s\n        '
        """
        x, y, z = min_p
        self._min_p = (float(x), float(y), float(z))
        x, y, z = max_p
        self._max_p = (float(x), float(y), float(z))
        self._texture = texture

    def str(self):
        r"""
        Return the scene string of the axis-aligned box.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Axis_aligned_box
            sage: aab = Axis_aligned_box((0,0,0),(1,1,1),'s')
            sage: aab.str()
            '\n        box min  0.0 0.0 0.0  max  1.0 1.0 1.0  s\n        '
        """
        return r"""
        box min %s max %s %s
        """ % (tostr(self._min_p), tostr(self._max_p), self._texture)


class TachyonTriangle(Triangle):
    r"""
    Basic triangle class.
    """
    def str(self):
        r"""
        Return the scene string for a triangle.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangle
            sage: t = TachyonTriangle([-1,-1,-1],[0,0,0],[1,2,3])
            sage: t.str()
            '\n        TRI V0  -1.0 -1.0 -1.0   V1  0.0 0.0 0.0    V2  1.0 2.0 3.0 \n            0\n        '
        """
        return r"""
        TRI V0 %s  V1 %s   V2 %s
            %s
        """ % (tostr(self._a), tostr(self._b), tostr(self._c), self._color)


class TachyonSmoothTriangle(SmoothTriangle):
    r"""
    A triangle along with a normal vector, which is used for smoothing.
    """
    def str(self):
        r"""
        Return the scene string for a smoothed triangle.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonSmoothTriangle
            sage: t = TachyonSmoothTriangle([-1,-1,-1],[0,0,0],[1,2,3],[1,0,0],[0,1,0],[0,0,1])
            sage: t.str()
            '\n        STRI V0  ...  1.0 0.0 0.0  N1  0.0 1.0 0.0   N2  0.0 0.0 1.0 \n             0\n        '
        """
        return r"""
        STRI V0 %s V1 %s  V2 %s
             N0 %s N1 %s  N2 %s
             %s
        """ % (tostr(self._a), tostr(self._b), tostr(self._c),
               tostr(self._da), tostr(self._db), tostr(self._dc), self._color)


class TachyonTriangleFactory(TriangleFactory):
    r"""
    A class to produce triangles of various rendering types.
    """
    def __init__(self, tach, tex):
        r"""
        Initialize with tachyon instance and texture.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttf._texture
            's'
        """
        self._tachyon = tach
        self._texture = tex

    def triangle(self, a, b, c, color=None):
        r"""
        Create a TachyonTriangle with vertices a, b, and c.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttft = ttf.triangle([1,2,3],[3,2,1],[0,2,1])
            sage: ttft.str()
            '\n        TRI V0  1.0 2.0 3.0   V1  3.0 2.0 1.0    V2  0.0 2.0 1.0 \n            s\n        '
        """
        if color is None:
            return TachyonTriangle(a, b, c, self._texture)
        else:
            return TachyonTriangle(a, b, c, color)

    def smooth_triangle(self, a, b, c, da, db, dc, color=None):
        r"""
        Create a TachyonSmoothTriangle.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttfst = ttf.smooth_triangle([0,0,0],[1,0,0],[0,0,1],[1,1,1],[1,2,3],[-1,-1,2])
            sage: ttfst.str()
            '\n        STRI V0  0.0 0.0 0.0  ...'
        """
        if color is None:
            return TachyonSmoothTriangle(a, b, c, da, db, dc, self._texture)
        else:
            return TachyonSmoothTriangle(a, b, c, da, db, dc, color)

    def get_colors(self, list):
        r"""
        Return a list of color labels.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttf.get_colors([(1,1,1)])
            ['SAGETEX1_0']
        """
        return self._tachyon.texture_recolor(self._texture, list)


class ParametricPlot():
    r"""
    Parametric plotting routines.
    """
    def str(self):
        r"""
        Return the tachyon string representation of the parameterized curve.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q.str()[9:69]
            'sphere center  0.0 0.0 0.0  rad 0.1 s\n        \n        fcyli'
        """
        return "".join(o.str() for o in self._objects)

    def __init__(self, f, t_0, t_f, tex, r=.1, cylinders=True,
                 min_depth=4, max_depth=8, e_rel=.01, e_abs=.01):
        r"""
        Create the parametric plotting class.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q._e_rel
            0.01
        """
        self._e_rel = e_rel
        self._e_abs = e_abs
        self._r = r
        self._f = f
        self._tex = tex
        self._cylinders = cylinders
        self._max_depth = max_depth
        self._min_depth = min_depth

        f_0 = f(t_0)
        f_f = f(t_f)
        self._objects = [Sphere(f_0, r, texture=tex)]

        self._plot_step(0, t_0, t_f, f_0, f_f)

    def _plot_step(self, depth, t_0, t_f, f_0, f_f):
        r"""
        Recursively subdivide interval, eventually plotting with cylinders and spheres.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q._plot_step(8,0,1,[0,0,0],[1,1,1])
            sage: len(q._objects)
            515
        """
        if depth < self._max_depth:
            t_mid = (t_f + t_0) / 2
            f_mid = ((f_f[0] + f_0[0]) / 2, (f_f[1] + f_0[1]) / 2, (f_f[2] + f_0[2]) / 2)
            f_val = self._f(t_mid)
            if depth < self._min_depth or self.tol(f_mid, f_val):
                new_depth = depth + 1
            else:
                new_depth = self._max_depth

            self._plot_step(new_depth, t_0, t_mid, f_0, f_val)
            self._plot_step(new_depth, t_mid, t_f, f_val, f_f)
        else:
            if self._cylinders:
                self._objects.append(FCylinder(f_0, f_f, self._r, self._tex))
            self._objects.append(Sphere(f_f, self._r, self._tex))

    def tol(self, est, val):
        r"""
        Check relative, then absolute tolerance.

        If both fail, return ``False``.

        This is a zero-safe error checker.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q.tol([0,0,0],[1,0,0])
            False
            sage: q.tol([0,0,0],[.0001,0,0])
            True
        """
        a, b, c = val
        delta = sqrt((a - est[0])**2 + (b - est[1])**2 + (c - est[2])**2)
        if delta < self._e_abs:
            return True

        r = sqrt(a**2 + b**2 + c**2)
        if delta < self._e_rel * r:
            return True

        return False


def tostr(s, length=3, out_type=float):
    r"""
    Convert vector information to a space-separated string.

    EXAMPLES::

        sage: from sage.plot.plot3d.tachyon import tostr
        sage: tostr((1,1,1))
        ' 1.0 1.0 1.0 '
        sage: tostr('2 3 2')
        '2 3 2'
    """
    if isinstance(s, str):
        return s
    output = ' '
    for an_item in s:
        output = output + str(out_type(an_item)) + ' '
    return output
