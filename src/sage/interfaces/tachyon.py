r"""
The Tachyon Ray Tracer

AUTHOR:

- John E. Stone
"""

#*****************************************************************************
#       Copyright (C) 2006 John E. Stone
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

from sage.cpython.string import bytes_to_str
from sage.misc.pager import pager
from sage.misc.temporary_file import tmp_filename
from sage.structure.sage_object import SageObject


class TachyonRT(SageObject):
    """
    The Tachyon Ray Tracer

    tachyon_rt(model, outfile='sage.png', verbose=1, block=True, extra_opts='')

    INPUT:

    -  ``model`` - a string that describes a 3d model in
       the Tachyon modeling format. Type tachyon_rt.help() for a
       description of this format.

    -  ``outfile`` - (default: 'sage.png') output filename;
       the extension of the filename determines the type. Supported types
       include:

       -  ``tga`` - 24-bit (uncompressed)

       -  ``bmp`` - 24-bit Windows BMP (uncompressed)

       -  ``ppm`` - 24-bit PPM (uncompressed)

       -  ``rgb`` - 24-bit SGI RGB (uncompressed)

       -  ``png`` - 24-bit PNG (compressed, lossless)

    -  ``verbose`` - integer; (default: 1)

       -  ``0`` - silent

       -  ``1`` - some output

       -  ``2`` - very verbose output

    -  ``block`` - bool (default: True); if False, run the
       rendering command in the background.

    -  ``extra_opts`` - passed directly to tachyon command
       line. Use tachyon_rt.usage() to see some of the possibilities.


    OUTPUT:

    - Some text may be displayed onscreen.

    - The file outfile is created.


    EXAMPLES:


    .. automethod:: __call__
    """
    def _repr_(self):
        """
        Returns a brief description of this interface object (the Tachyon raytracer written by John Stone).

        TESTS::

            sage: from sage.interfaces.tachyon import TachyonRT
            sage: t = TachyonRT()
            sage: print(t.__repr__())
            John Stone's Tachyon Ray Tracer
        """
        return "John Stone's Tachyon Ray Tracer"

    def __call__(self, model, outfile='sage.png', verbose=1, extra_opts=''):
        """
        This executes the tachyon program, given a scene file input.

        INPUT:

        - ``model`` -- string. The tachyon model.

        - ``outfile`` -- string, default ``'sage.png'``. The filename
          to save the model to.

        - ``verbose`` -- 0, 1, (default) or 2. The verbosity level.

        - ``extra_opts`` -- string (default: empty string). Extra
          options that will be appended to the tachyon commandline.

        EXAMPLES::

            sage: from sage.interfaces.tachyon import TachyonRT
            sage: tgen = Tachyon()
            sage: tgen.texture('t1')
            sage: tgen.sphere((0,0,0),1,'t1')
            sage: tgen.str()[30:40]
            'resolution'
            sage: t = TachyonRT()
            sage: import os
            sage: t(tgen.str(), outfile=os.devnull)
            tachyon ...
            Tachyon Parallel/Multiprocessor Ray Tracer...

        TESTS::

            sage: from sage.env import SAGE_EXTCODE
            sage: filename = os.path.join(SAGE_EXTCODE, 'doctest', 'invalid', 'syntax_error.tachyon')
            sage: with open(filename, 'r') as f:
            ....:     syntax_error = f.read()
            sage: t(syntax_error, outfile=os.devnull)
            Traceback (most recent call last):
            ...
            RuntimeError: Tachyon Parallel/Multiprocessor Ray Tracer...
            ...
            Parser failed due to an input file syntax error.
            Aborting render.
        """
        modelfile = tmp_filename(ext='.dat')
        with open(modelfile, 'w') as file:
            file.write(model)
        cmd = ['tachyon', modelfile]
        ext = outfile[-4:].lower()
        if ext == '.png':
            cmd += ['-format', 'PNG']
        elif ext == '.tga':
            cmd += ['-format', 'TARGA']
        elif ext == '.bmp':
            cmd += ['-format', 'BMP']
        elif ext == '.ppm':
            cmd += ['-format', 'PPM']
        elif ext == '.rgb':
            cmd += ['-format', 'RGB']
        cmd += ['-o', outfile]
        cmd += extra_opts.split()
        if verbose >= 2:
            cmd += ['+V']
        if verbose:
            print(' '.join(cmd))
        import subprocess
        out = bytes_to_str(subprocess.check_output(cmd))
        if verbose >= 1:
            print(out)
        if out.rstrip().endswith('Aborting render.'):
            raise RuntimeError(out)
        if outfile != os.devnull and os.stat(outfile).st_size == 0:
            raise RuntimeError('tachyon did not abort but output file is empty')

    def usage(self, use_pager=True):
        """
        Returns the basic description of using the Tachyon raytracer (simply what is returned by running tachyon with no input).  The output is paged unless use_pager=False.

        TESTS::

            sage: from sage.interfaces.tachyon import TachyonRT
            sage: t = TachyonRT()
            sage: t.usage(use_pager=False)
            ...
              tachyon modelfile [options]...
            <BLANKLINE>
            Model file formats supported:
              filename.dat ...
        """
        with os.popen('tachyon') as f:
            r = f.read()
        if use_pager:
            pager()(r)
        else:
            print(r)

    def help(self, use_pager=True):
        """
        Prints (pages) the help file written by John Stone describing scene files for Tachyon.  The output is paged unless use_pager=False.

        TESTS::

            sage: from sage.interfaces.tachyon import TachyonRT
            sage: t = TachyonRT()
            sage: t.help(use_pager=False)
            This help, which was written by John Stone, describes ...
        """
        s = r"""
This help, which was written by John Stone, describes how to create
scene files.

At the present time, scene description files are very simple.
The parser can't handle multiple file scene descriptions, although they
may be added in the future.  Most of the objects and their scene description
are closely related to the RAY API
\emph{(See the API docs for additional info.)}

\subsection{Basic Scene Requirements}
  Unlike some other ray tracers out there, RAY requires that you
specify most of the scene parameters in the scene description file itself.
If users would rather specify some of these parameters at the command line,
then I may add that feature in the future.
A scene description file contains keywords, and values associated or grouped
with a keyword.  All keywords can be in caps, lower case, or mixed case
for the convenience of the user.  File names and texture names are
normally case-sensitive, although the behavior for file names is
operating system-dependent.  All values are either character strings, or
floating point numbers.  In some cases, the presence of one keyword will
require additional keyword / value pairs.

  At the moment there are several keywords with values,
that must appear in every scene description file.
Every scene description file must begin with the
{\bf BEGIN\_SCENE} keyword, and end with the {\bf END\_SCENE} keyword.
All definitions and declarations of any kind must be inside the
{\bf BEGIN\_SCENE}, {\bf END\_SCENE} pair.
The {\bf RESOLUTION} keyword is followed by an x resolution
and a y resolution in terms of pixels on each axis.  There are currently
no limits placed on the resolution of an output image other than the
computer's available memory and reasonable execution time.
An example of a simple scene description skeleton is show below:
\begin{verbatim}
BEGIN_SCENE
  RESOLUTION 1024 1024
...
...  Camera definition..
...
...  Other objects, etc..
...

END_SCENE
\end{verbatim}

\subsection{Camera and viewing parameters}
  One of the most important parts of any scene, is the camera position and
orientation.  Having a good angle on a scene can make the difference between
an average looking scene and a strikingly interesting one.  There may be
multiple camera definitions in a scene file, but the last camera definition
overrides all previous definitions.
There are several parameters that control the camera in \RAY,
{\bf PROJECTION}, {\bf ZOOM}, {\bf ASPECTRATIO}, {\bf ANTIALIASING},
 {\bf CENTER}, {\bf RAYDEPTH}, {\bf VIEWDIR}, and {\bf UPDIR}.

The first and last keywords required in the definition of a camera are the
{\bf CAMERA} and {\bf END\_CAMERA} keywords.  The {\bf PROJECTION} keyword
is optional, the remaining camera keywords are required, and must be
written in the sequence they are listed in the examples in this section.

\subsubsection{Camera projection modes}
  The {\bf PROJECTION} keyword must be followed by one of the supported
camera projection mode identifiers {\bf PERSPECTIVE}, {\bf PERSPECTIVE_DOF},
{\bf ORTHOGRAPHIC}, or {\bf FISHEYE}.  The {\bf FISHEYE} projection mode
requires two extra parameters {\bf FOCALLENGTH} and {\bf APERTURE}
which precede the regular camera options.

\begin{verbatim}
Camera
  projection perspective_dof
  focallength 0.75
  aperture 0.02
  Zoom 0.666667
  Aspectratio 1.000000
  Antialiasing 128
  Raydepth 30
  Center  0.000000 0.000000 -2.000000
  Viewdir -0.000000 -0.000000 2.000000
  Updir   0.000000 1.000000 -0.000000
End_Camera
\end{verbatim}

\subsubsection{Common camera parameters}
  The {\bf ZOOM} parameter controls the camera in a way similar to a
telephoto lens on a 35mm camera.  A zoom value of 1.0 is standard,
with a 90 degree field of view.  By changing the zoom factor to 2.0,
the relative size of any feature in the frame is twice as big, while
the field of view is decreased slightly.  The zoom effect is
implemented as a scaling factor on the height and width of the image
plane relative to the world.

  The {\bf ASPECTRATIO} parameter controls the aspect ratio of the resulting
image.  By using the aspect ratio parameter, one can produce images which
look correct on any screen.  Aspect ratio alters the relative width of the
image plane, while keeping the height of the image plane constant.  In
general, most workstation displays have an aspect ratio of 1.0.  To see
what aspect ratio your display has, you can render a simple sphere, at
a resolution of 512x512 and measure the ratio of its width to its height.

The {\bf ANTIALIASING} parameter controls the maximum level of supersampling
used to obtain higher image quality.  The parameter given sets the number of
additional rays to trace per-pixel to attain higher image quality.

  The {\bf RAYDEPTH} parameter tells RAY what the maximum
level of reflections, refraction, or in general the maximum recursion
depth to trace rays to.  A value between 4 and 12 is usually good.  A
value of 1 will disable rendering of reflective or transmissive
objects (they'll be black).

  The remaining three camera parameters are the most important, because
they define the coordinate system of the camera, and its position in the
scene.  The {\bf CENTER} parameter is an X, Y, Z coordinate defining the
center of the camera \emph{(also known as the Center of Projection)}.
Once you have determined where the camera will be placed in the scene, you
need to tell RAY what the camera should be looking at.  The
{\bf VIEWDIR} parameter is a vector indicating the direction the camera
is facing.  It may be useful for me to add a "Look At" type keyword in
the future to make camera aiming easier.  If people want or need the
"Look At" style camera, let me know.  The last parameter needed to completely
define a camera is the "up" direction.  The {\bf UPDIR} parameter is a vector
which points in the direction of the "sky".  I wrote the camera so that
{\bf VIEWDIR} and {\bf UPDIR} don't have to be perpendicular, and there
shouldn't be a need for a "right" vector although some other ray tracers
require it.  Here's a snippet of a camera definition:
\begin{verbatim}
CAMERA
  ZOOM 1.0
  ASPECTRATIO 1.0
  ANTIALIASING 0
  RAYDEPTH 12
  CENTER 0.0 0.0 2.0
  VIEWDIR 0 0 -1
  UPDIR 0 1 0
END_CAMERA
\end{verbatim}


\subsubsection{Viewing frustum}
An optional {\bf FRUSTUM} parameter provides a means for rendering sub-images
in a larger frame, and correct stereoscopic images.  The {\bf FRUSTUM}
keyword must be followed by four floating parameters, which indicate
the top, bottom, left and right coordinates of the image plane in
eye coordinates.   When the projection mode is set to {\bf FISHEYE},
the frustum parameters correspond to spherical coordinates specified
in radians.

\begin{verbatim}
CAMERA
  ZOOM  1.0
  ASPECTRATIO 1.0
  ANTIALIASING 0
  RAYDEPTH  4
  CENTER    0.0 0.0 -6.0
  VIEWDIR   0.0 0.0 1.0
  UPDIR     0.0 1.0 0.0
  FRUSTUM   -0.5 0.5 -0.5 0.5
END_CAMERA
\end{verbatim}


\subsection{Including Files}
The {\bf INCLUDE} keyword is used anywhere after the camera description,
and is immediately followed by a valid filename, for a file containing
additional scene description information.  The included file is opened,
and processing continues as if it were part of the current file, until
the end of the included file is reached.  Parsing of the current file
continues from where it left off prior to the included file.

\subsection{Scene File Comments}
The {\bf $\#$} keyword is used anywhere after the camera description, and
will cause RAY to ignore all characters from the {\bf $\#$} to the end
of the input line.  The {\bf $\#$} character must be surrounded by whitespace
in order to be recognized.  A sequence such as {\bf $\#\#\#$} will not be
recognized as a comment.

\subsection{Lights}
The most frequently used type of lights provided by RAY are positional
point light sources.  The lights are actually small spheres, which are
visible.  A point light is composed of three pieces of
information, a center, a radius (since its a sphere), and a color.
To define a light, simply write the {\bf LIGHT} keyword, followed by
its {\bf CENTER} (a X, Y, Z coordinate), its {\bf RAD} (radius, a scalar),
and its {\bf COLOR} (a Red Green Blue triple).  The radius parameter will
accept any value of 0.0 or greater.  Lights of radius 0.0 will not be
directly visible in the rendered scene, but contribute light to the scene
normally.
For a light, the color values
range from 0.0 to 1.0, any values outside this range may yield unpredictable
results.  A simple light definition looks like this:
\begin{verbatim}
  LIGHT CENTER 4.0 3.0 2.0
        RAD    0.2
        COLOR  0.5 0.5 0.5
\end{verbatim}
This light would be gray colored if seen directly, and would be 50\%
intensity in each RGB color component.


RAY supports simple directional lighting, commonly used in
CAD and scientific visualization programs for its performance
advantages over positional lights.  Directional lights cannot be
seen directly in scenes rendered by \RAY, only their illumination
contributes to the final image.

\begin{verbatim}
DIRECTIONAL_LIGHT
  DIRECTION 0.0 -1.0 0.0
  COLOR   1.0 0.0  0.0
\end{verbatim}

RAY supports spotlights, which are described very similarly to a
point light, but they are attenuated by angle from the direction vector,
based on a  ``falloff start'' angle and ``falloff end''angle.  Between
the starting and ending angles, the illumination is attenuated linearly.
The syntax for a spotlight description in a scene file is as follows.
\begin{verbatim}
SPOTLIGHT
  CENTER  0.0 3.0  17.0
  RAD     0.2
  DIRECTION 0.0 -1.0 0.0
    FALLOFF_START 20.0
    FALLOFF_END   45.0
  COLOR   1.0 0.0  0.0
\end{verbatim}

The lighting system implemented by RAY provides various levels of
distance-based lighting attenuation.  By default, a light is not attenuated
by distance.  If the \emph{attenuation} keywords is present immediately
prior to the light's color, RAY will accept coefficients which are used
to calculate distance-based attenuation, which is applied the light by
multiplying with the resulting value.  The attenuation factor is calculated
from the equation
$$
  1/(K_c + K_l d + k_q d^2)
$$

This attenuation equation should be familiar to some as it
is the same lighting attenuation equation used by OpenGL.
The constant, linear, and quadratic terms are specified in a scene file
as shown in the following example.
\begin{verbatim}
LIGHT
  CENTER  -5.0 0.0 10.0
  RAD     1.0
  ATTENUATION CONSTANT 1.0 LINEAR 0.2 QUADRATIC 0.05
  COLOR   1.0 0.0 0.0
\end{verbatim}



\subsection{Atmospheric effects}
RAY currently only implements one atmospheric effect,
simple distance-based fog.

\subsubsection{Fog}
RAY provides a simple distance-based fog effect intended to provide
functionality similar to that found in OpenGL, for compatibility with
software that requires an OpenGL-like fog implementation.  Much like
OpenGL, RAY provides linear, exponential, and exponential-squared fog.

\begin{verbatim}
  FOG
    LINEAR START 0.0  END 50.0  DENSITY 1.0  COLOR 1.0 1.0 1.0
\end{verbatim}

\begin{verbatim}
  FOG
    EXP START 0.0  END 50.0  DENSITY 1.0  COLOR 1.0 1.0 1.0
\end{verbatim}

\begin{verbatim}
  FOG
    EXP2 START 0.0  END 50.0  DENSITY 1.0  COLOR 1.0 1.0 1.0
\end{verbatim}


\subsection{Objects}

\subsubsection{Spheres}
  Spheres are the simplest object supported by RAY and they are
also the fastest object to render.  Spheres are defined as one would expect,
with a {\bf CENTER}, {\bf RAD} (radius), and a texture.  The texture may
be defined along with the object as discussed earlier, or it may be declared
and assigned a name.
Here's a sphere definition using a previously defined "NitrogenAtom" texture:
\begin{verbatim}
 SPHERE  CENTER 26.4 27.4 -2.4   RAD 1.0   NitrogenAtom
\end{verbatim}
A sphere with an inline texture definition is declared like this:
\begin{verbatim}
 Sphere center 1.0 0.0 10.0
           Rad 1.0
        Texture  Ambient 0.2  Diffuse 0.8  Specular 0.0  Opacity 1.0
                 Color   1.0 0.0 0.5
                 TexFunc 0
\end{verbatim}
Notice that in this example I used mixed case for the keywords, this is
allowable...
Review the section on textures if the texture definitions are confusing.

\subsubsection{Triangles}
  Triangles are also fairly simple objects, constructed by listing the
three vertices of the triangle, and its texture.  The order of the
vertices isn't important, the triangle object is "double sided", so the
surface normal is always pointing back in the direction of the incident ray.
The triangle vertices are listed as {\bf V1}, {\bf V2}, and {\bf V3} each one
is an X, Y, Z coordinate.  An example of a triangle is shown below:
\begin{verbatim}
TRI
  V0  0.0 -4.0 12.0
  V1  4.0 -4.0 8.0
  V2 -4.0 -4.0 8.0
  TEXTURE
    AMBIENT  0.1 DIFFUSE  0.2 SPECULAR 0.7 OPACITY 1.0
    COLOR 1.0 1.0 1.0
    TEXFUNC 0
\end{verbatim}

\subsubsection{Smoothed Triangles}
  Smoothed triangles are just like regular triangles, except that the
  surface normal for each of the three vertices is used to determine the
  surface normal across the triangle by linear interpolation.
  Smoothed triangles yield curved looking objects and have nice
  reflections.
\begin{verbatim}
STRI
  V0 1.4   0.0   2.4
  V1 1.35 -0.37  2.4
  V2 1.36 -0.32  2.45
  N0 -0.9 -0.0  -0.4
  N1 -0.8  0.23 -0.4
  N2 -0.9  0.27 -0.15
  TEXTURE
    AMBIENT  0.1 DIFFUSE  0.2 SPECULAR 0.7 OPACITY 1.0
    COLOR 1.0 1.0 1.0
    TEXFUNC 0
\end{verbatim}

\subsubsection{Infinite Planes}

  Useful for things like desert floors, backgrounds, skies etc, the infinite
plane is pretty easy to use.  An infinite plane only consists of two pieces
of information, the {\bf CENTER} of the plane, and a {\bf NORMAL} to the plane.
The center of the plane is just any point on the plane such that the point
combined with the surface normal define the equation for the plane.
As with triangles, planes are double sided.  Here is an example of an
infinite plane:
\begin{verbatim}
PLANE
  CENTER 0.0 -5.0 0.0
  NORMAL 0.0  1.0 0.0
  TEXTURE
    AMBIENT  0.1 DIFFUSE  0.9 SPECULAR 0.0  OPACITY 1.0
    COLOR  1.0 1.0 1.0
    TEXFUNC  1
      CENTER 0.0 -5.0 0.0
      ROTATE 0. 0.0 0.0
      SCALE  1.0 1.0 1.0
\end{verbatim}

\subsubsection{Rings}
  Rings are a simple object, they are really a not-so-infinite plane.
Rings are simply an infinite plane cut into a washer shaped ring, infinitely
thing just like a plane.  A ring only requires two more pieces of information
than an infinite plane does, an inner and outer radius.  Here's an example
of a ring:
\begin{verbatim}
  Ring
    Center 1.0 1.0 1.0
    Normal 0.0 1.0 0.0
    Inner  1.0
    Outer  5.0
    MyNewRedTexture
\end{verbatim}

\subsubsection{Infinite Cylinders}
  Infinite cylinders are quite simple.  They are defined by a center, an
axis, and a radius.  An example of an infinite cylinder is:
\begin{verbatim}
  Cylinder
    Center 0.0 0.0 0.0
    Axis   0.0 1.0 0.0
    Rad    1.0
    SomeRandomTexture
\end{verbatim}

\subsubsection{Finite Cylinders}
  Finite cylinders are almost the same as infinite ones, but the
  center and length of the axis determine the extents of the cylinder.
  The finite cylinder is also really a shell, it doesn't have any
  caps.  If you need to close off the ends of the cylinder, use two
  ring objects, with the inner radius set to 0.0 and the normal set
  to be the axis of the cylinder.  Finite cylinders are built this
  way to enhance speed.

\begin{verbatim}
  FCylinder
    Center 0.0 0.0 0.0
    Axis   0.0 9.0 0.0
    Rad    1.0
    SomeRandomTexture
\end{verbatim}
This defines a finite cylinder with radius 1.0, going from 0.0 0.0 0.0, to
0.0 9.0 0.0 along the Y axis.  The main difference between an infinite cylinder
and a finite cylinder is in the interpretation of the {\bf AXIS} parameter.
In the case of the infinite cylinder, the length of the axis vector is
ignored.  In the case of the finite cylinder, the axis parameter is used
to determine the length of the overall cylinder.

\subsubsection{Axis Aligned Boxes}
  Axis aligned boxes are fast, but of limited usefulness.  As such, I'm
not going to waste much time explaining 'em.  An axis aligned box is
defined by a {\bf MIN} point, and a {\bf MAX} point.  The volume between
the min and max points is the box.  Here's a simple box:
\begin{verbatim}
  BOX
    MIN -1.0 -1.0 -1.0
    MAX  1.0  1.0  1.0
    Boxtexture1
\end{verbatim}

\subsubsection{Fractal Landscapes}
  Currently fractal landscapes are a built-in function.  In the near future
I'll allow the user to load an image map for use as a heightfield.
Fractal landscapes are currently forced to be axis aligned.  Any suggestion
on how to make them more appealing to users is welcome.  A fractal landscape
is defined by its "resolution" which is the number of grid points along
each axis, and by its scale and center.  The "scale" is how large the
landscape is along the X, and Y axes in world coordinates.  Here's a simple
landscape:
\begin{verbatim}
SCAPE
  RES 30 30
  SCALE 80.0 80.0
  CENTER 0.0 -4.0 20.0
  TEXTURE
    AMBIENT 0.1 DIFFUSE 0.9 SPECULAR 0.0 OPACITY 1.0
    COLOR 1.0 1.0 1.0
    TEXFUNC 0
\end{verbatim}
The landscape shown above generates a square landscape made of 1,800 triangles.
When time permits, the heightfield code will be rewritten to be more
general and to increase rendering speed.

\subsubsection{Arbitrary Quadric Surfaces}
  Docs soon. I need to add these into the parser, must have forgotten
before ;-)

\subsubsection{Volume Rendered Scalar Voxels}
These are a little trickier than the average object :-)
These are likely to change substantially in the very near future so I'm not
going to get too detailed yet.
A volume rendered data set is described by its axis aligned bounding box, and
its resolution along each axis.  The final parameter is the voxel data
file.  If you are seriously interested in messing with these, get hold of
me and I'll give you more info.  Here's a quick example:
\begin{verbatim}
SCALARVOL
  MIN -1.0 -1.0 -0.4
  MAX  1.0  1.0  0.4
  DIM 256 256 100
  FILE /cfs/johns/vol/engine.256x256x110
  TEXTURE
        AMBIENT 1.0 DIFFUSE 0.0 SPECULAR 0.0 OPACITY 8.1
        COLOR 1.0 1.0 1.0
        TEXFUNC 0
\end{verbatim}

\subsection{Texture and Color}
\subsubsection{Simple Texture Characteristics}
  The surface textures applied to an object drastically alter its overall
appearance, making textures and color one of the most important topics in
this manual.  As with many other renderers, textures can be declared and
associated with a name so that they may be used over and over again in
a scene definition with less typing.  If a texture is only need once, or it
is unique to a particular object in the scene, then it may be declared along
with the object it is applied to, and does not need a name.

  The simplest texture definition is a solid color with no image mapping
or procedural texture mapping.  A solid color texture is defined by the
{\bf AMBIENT}, {\bf DIFFUSE}, {\bf SPECULAR}, {\bf OPACITY} and {\bf COLOR}
parameters.  The {\bf AMBIENT} parameter defines the ambient lighting
coefficient to be used when shading the object.  Similarly, the {\bf DIFFUSE}
parameter is the relative contribution of the diffuse shading to the surface
appearance.  The {\bf SPECULAR} parameter is the contribution from perfectly
reflected rays, as if on a mirrored surface.  {\bf OPACITY} defines how
transparent a surface is.  An {\bf OPACITY} value of 0.0 renders the object
completely invisible.  An {\bf OPACITY} value of 1.0 makes the object
completely solid, and non-transmissive.  In general, the values for the
ambient, diffuse, and specular parameters should add up to 1.0, if they don't
then pixels may be over or underexposed quite easily.  These parameters
function in a manner similar to that of other ray tracers.  The {\bf COLOR}
parameter is an RGB triple with each value ranging from 0.0 to 1.0 inclusive.
If the RGB values stray from 0.0 to 1.0, results are undefined.
In the case of solid textures, a final parameter, {\bf TEXFUNC} is set to
zero (integer).

\subsubsection{Texture Declaration and Aliasing}
  To define a simple texture for use on several objects in a scene, the
{\bf TEXDEF} keyword is used.  The {\bf TEXDEF} keyword is followed by
a case sensitive texture name, which will subsequently be used while
defining objects.  If many objects in a scene use the same texture through
texture definition, a significant amount of memory may be saved since only
one copy of the texture is present in memory, and its shared by all
of the objects.  Here is an example of a solid texture definition:
\begin{verbatim}
 TEXDEF MyNewRedTexture
    AMBIENT 0.1 DIFFUSE 0.9 SPECULAR 0.0 OPACITY 1.0
    COLOR 1.0 0.0 0.0  TEXFUNC 0
\end{verbatim}
When this texture is used in an object definition, it is referenced only by
name.  Be careful not to use one of the other keywords as a defined texture,
this will probably cause the parser to explode, as I don't check for use
of keywords as texture names.

  When a texture is declared within an object definition, it appears in
an identical format to the {\bf TEXDEF} declaration, but the {\bf TEXTURE}
keyword is used instead of {\bf TEXDEF}.  If it is useful to have several
names for the same texture (when you are too lazy to actually finish defining
different variations of a wood texture for example, and just want to be
approximately correct for example) aliases can be constructed using the
{\bf TEXALIAS} keyword, along with the alias name, and the original name.
An example of a texture alias is:
\begin{verbatim}
  TEXALIAS MyNewestRedTexture  MyNewRedTexture
\end{verbatim}
This line would alias MyNewestRedTexture to be the same thing as the
previously declared MyNewRedTexture.  Note that the source texture must
be declared before any aliases that use it.

\subsubsection{Image Maps and Procedural Textures} Image maps and
procedural textures very useful in making realistic looking scenes.  A
good image map can do as much for the realism of a wooden table as any
amount of sophisticated geometry or lighting.  Image maps are made by
wrapping an image on to an object in one of three ways, a spherical
map, a cylindrical map, and a planar map.  Procedural textures are
used in a way similar to the image maps, but they are on the fly and
do not use much memory compared to the image maps.  The main
disadvantage of the procedural maps is that they must be hard-coded
into RAY when it is compiled.

  The syntax used for all texture maps is fairly simple to learn.  The biggest
problem with the way that the parser is written now is that the different
mappings are selected by an integer, which is not very user friendly.  I
expect to rewrite this section of the parser sometime in the near future to
alleviate this problem.  When I rewrite the parser, I may also end up altering
the parameters that are used to describe a texture map, and some of them may
become optional rather than required.

\begin{center}
\begin{tabular}{|c|c|}
\multicolumn{2}{c}{Texture Mapping Functions} \\
\hline
{Value for TEXFUNC} & {Mapping and Texture Description}\\
\hline
{0} & {No special texture, plain shading}  \\
{1} & {3D checkerboard function, like a Rubik's cube}  \\
{2} & {Grit Texture, randomized surface color}  \\
{3} & {3D marble texture, uses object's base color}  \\
{4} & {3D wood texture, light and dark brown, not very good yet}  \\
{5} & {3D gradient noise function (can't remember what it look like}  \\
{6} & {Don't remember}  \\
{7} & {Cylindrical Image Map, requires ppm filename}  \\
{8} & {Spherical Image Map, requires ppm filename}  \\
{9} & {Planar Image Map, requires ppm filename}  \\
\hline
\end{tabular}
\end{center}

Here's an example of a sphere, with a spherical image map applied to its
surface:
\begin{verbatim}
SPHERE
  CENTER 2.0  0.0 5.0
  RAD 2.0
  TEXTURE
    AMBIENT  0.4 DIFFUSE  0.8 SPECULAR 0.0  OPACITY 1.0
    COLOR 1.0 1.0 1.0
    TEXFUNC 7 /cfs/johns/imaps/fire644.ppm
      CENTER 2.0 0.0 5.0
      ROTATE 0.0 0.0 0.0
      SCALE  2.0 -2.0 1.0
\end{verbatim}

Basically, the image maps require the center, rotate and scale
parameters so that you can position the image map on the object
properly.
"""
        from sage.misc.sagedoc import format
        f = format(s)
        f = f.replace('{ ','').replace('}','').replace('{','')
        if use_pager:
            pager()(f)
        else:
            print(f)

tachyon_rt = TachyonRT()


