# sage.doctest: needs sage.symbolic
r"""
Animated plots

Animations are generated from a list (or other iterable) of graphics
objects.
Images are produced by calling the ``save_image`` method on each input
object, creating a sequence of PNG files.
These are then assembled to various target formats using different
tools.
In particular, the ``magick/convert`` program from ImageMagick_ can be used to
generate an animated GIF file.
FFmpeg_ (with the command line program ``ffmpeg``) provides support for
various video formats, but also an alternative method of generating
animated GIFs.
For `browsers which support it`_, APNG_ can be used as another
alternative which works without any extra dependencies.

.. WARNING::

    Note that ``ImageMagick`` and ``FFmpeg`` are not included with Sage, and
    must be installed by the user.  On unix systems, type ``which
    magick`` at a command prompt to see if ``magick`` (part of the
    ``ImageMagick`` suite) is installed.  If it is, you will be given its
    location.  Similarly, you can check for ``ffmpeg`` with ``which
    ffmpeg``.  See the websites of ImageMagick_ or FFmpeg_ for
    installation instructions.

EXAMPLES:

The sine function::

    sage: x = SR.var("x")
    sage: sines = [plot(c*sin(x), (-2*pi,2*pi), color=Color(c,0,0), ymin=-1, ymax=1)
    ....:          for c in sxrange(0,1,.2)]
    sage: a = animate(sines)
    sage: print(a)
    Animation with 5 frames
    sage: a.show()                  # long time  # optional -- ImageMagick

Animate using FFmpeg_ instead of ImageMagick::

    sage: a.show(use_ffmpeg=True)   # long time  # optional -- FFmpeg

Animate as an APNG_::

    sage: a.apng(show_path=True)    # long time
    Animation saved to ....png.

An animated :class:`sage.plot.multigraphics.GraphicsArray` of rotating ellipses::

    sage: E = animate((graphics_array([[ellipse((0,0), a, b, angle=t, xmin=-3, xmax=3)
    ....:                                + circle((0,0), 3, color='blue')
    ....:                               for a in range(1,3)]
    ....:                              for b in range(2,4)])
    ....:              for t in sxrange(0, pi/4, .15)))
    sage: str(E)    # animations produced from a generator do not have a known length
    'Animation with unknown number of frames'
    sage: E.show()                  # long time  # optional -- ImageMagick

A simple animation of a circle shooting up to the right::

    sage: c = animate([circle((i,i), 1 - 1/(i+1), hue=i/10)
    ....:              for i in srange(0, 2, 0.2)],
    ....:             xmin=0, ymin=0, xmax=2, ymax=2, figsize=[2,2])
    sage: c.show()                  # long time  # optional -- ImageMagick


Animations of 3d objects::

    sage: s,t = SR.var("s,t")
    sage: def sphere_and_plane(x):
    ....:     return (sphere((0,0,0), 1, color='red', opacity=.5)
    ....:              + parametric_plot3d([t,x,s], (s,-1,1), (t,-1,1),
    ....:                                  color='green', opacity=.7))
    sage: sp = animate([sphere_and_plane(x)
    ....:               for x in sxrange(-1, 1, .3)])
    sage: sp[0]      # first frame
    Graphics3d Object
    sage: sp[-1]     # last frame
    Graphics3d Object
    sage: sp.show()                 # long time  # optional -- ImageMagick

    sage: (x,y,z) = SR.var("x,y,z")
    sage: def frame(t):
    ....:     return implicit_plot3d((x^2 + y^2 + z^2),
    ....:                            (x, -2, 2), (y, -2, 2), (z, -2, 2),
    ....:                            plot_points=60, contour=[1,3,5],
    ....:                            region=lambda x,y,z: x<=t or y>=t or z<=t)
    sage: a = animate([frame(t) for t in srange(.01, 1.5, .2)])
    sage: a[0]                      # long time
    Graphics3d Object
    sage: a.show()                  # long time  # optional -- ImageMagick

If the input objects do not have a ``save_image`` method, then the
animation object attempts to make an image by calling its internal
method :meth:`sage.plot.animate.Animation.make_image`.  This is
illustrated by the following example::

    sage: t = SR.var("t")
    sage: a = animate((sin(c*pi*t) for c in sxrange(1, 2, .2)))
    sage: a.show()                  # long time  # optional -- ImageMagick


AUTHORS:

- William Stein
- John Palmieri
- Niles Johnson (2013-12): Expand to animate more graphics objects
- Martin von Gagern (2014-12): Added APNG support
- Joshua Campbell (2020): interactive animation via Three.js viewer

REFERENCES:

- `ImageMagick <https://www.imagemagick.org>`_
- `FFmpeg <https://www.ffmpeg.org>`_
- `APNG <https://wiki.mozilla.org/APNG_Specification>`_
- `browsers which support it <https://caniuse.com/#feat=apng>`_
"""

############################################################################
#  Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
############################################################################

import builtins
import os
import shlex
import struct
import zlib

from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject
from sage.misc.temporary_file import tmp_dir, tmp_filename
from . import plot
import sage.misc.misc
import sage.misc.viewer


def animate(frames, **kwds):
    r"""
    Animate a list of frames by creating a
    :class:`sage.plot.animate.Animation` object.

    EXAMPLES::

        sage: t = SR.var("t")
        sage: a = animate((cos(c*pi*t) for c in sxrange(1, 2, .2)))
        sage: a.show()              # long time  # optional -- ImageMagick

    See also :mod:`sage.plot.animate` for more examples.
    """
    return Animation(frames, **kwds)


class Animation(WithEqualityById, SageObject):
    r"""
    Return an animation of a sequence of plots of objects.

    INPUT:

    - ``v`` -- iterable of Sage objects; these should preferably be
      graphics objects, but if they aren't, then :meth:`make_image` is
      called on them.

    - ``xmin``, ``xmax``, ``ymin``, ``ymax`` -- the ranges of the x and y axes

    - ``**kwds`` -- all additional inputs are passed onto the rendering
      command. E.g., use ``figsize`` to adjust the resolution and aspect
      ratio.

    EXAMPLES::

        sage: x = SR.var("x")
        sage: a = animate([sin(x + float(k)) for k in srange(0, 2*pi, 0.3)],
        ....:             xmin=0, xmax=2*pi, figsize=[2,1])
        sage: print(a)
        Animation with 21 frames
        sage: print(a[:5])
        Animation with 5 frames
        sage: a.show()              # long time  # optional -- ImageMagick
        sage: a[:5].show()          # long time  # optional -- ImageMagick

    The :meth:`show` method takes arguments to specify the
    delay between frames (measured in hundredths of a second, default
    value 20) and the number of iterations (default: 0, which
    means to iterate forever). To iterate 4 times with half a second
    between each frame::

        sage: a.show(delay=50, iterations=4)  # long time  # optional -- ImageMagick

    An animation of drawing a parabola::

        sage: step = 0.1
        sage: L = Graphics()
        sage: v = []
        sage: for i in srange(0, 1, step):
        ....:     L += line([(i,i^2),(i+step,(i+step)^2)], rgbcolor=(1,0,0), thickness=2)
        ....:     v.append(L)
        sage: a = animate(v, xmin=0, ymin=0)
        sage: a.show()              # long time  # optional -- ImageMagick
        sage: show(L)

    TESTS:

    This illustrates that :issue:`2066` is fixed (setting axes
    ranges when an endpoint is 0)::

        sage: animate([plot(sin, -1,1)], xmin=0, ymin=0)._kwds['xmin']
        0

    We check that :issue:`7981` is fixed::

        sage: x = SR.var("x")
        sage: a = animate([plot(sin(x + float(k)), (0, 2*pi), ymin=-5, ymax=5)
        ....:              for k in srange(0,2*pi,0.3)])
        sage: a.show()              # long time  # optional -- ImageMagick

    Do not convert input iterator to a list, but ensure that
    the frame count is known after rendering the frames::

        sage: x = SR.var("x")
        sage: a = animate((plot(x^p, (x,0,2)) for p in sxrange(1,2,.1)))
        sage: str(a)
        'Animation with unknown number of frames'
        sage: a.png()    # long time
        '.../'
        sage: len(a)     # long time
        10
        sage: a._frames
        <generator object ...

        sage: from sage.plot.animate import Animation
        sage: hash(Animation())  # random
        140658972348064
    """
    def __init__(self, v=None, **kwds):
        r"""
        Return an animation of a sequence of plots of objects.  See
        documentation of :func:`animate` for more details and
        examples.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.3)],  # indirect doctest
            ....:                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: print(a)
            Animation with 21 frames
            sage: a.show()          # long time  # optional -- ImageMagick
        """
        self._frames = v
        self._kwds = kwds

    def _combine_kwds(self, *kwds_tuple):
        """
        Return a dictionary which is a combination of the all the
        dictionaries in kwds_tuple. This also does the appropriate thing
        for taking the mins and maxes of all of the x/y mins/maxes.

        EXAMPLES::

            sage: a = animate([plot(sin, -1,1)], xmin=0, ymin=0)
            sage: kwds1 = {'a':1, 'b':2, 'xmin':2, 'xmax':5}
            sage: kwds2 = {'b':3, 'xmin':0, 'xmax':4}
            sage: kwds = a._combine_kwds(kwds1, kwds2)
            sage: list(sorted(kwds.items()))
            [('a', 1), ('b', 3), ('xmax', 5), ('xmin', 0)]

        Test that the bug reported in :issue:`12107` has been fixed::

            sage: kwds3 = {}
            sage: kwds4 = {'b':3, 'xmin':0, 'xmax':4}
            sage: a._combine_kwds(kwds3, kwds4)['xmin']
            0
        """
        new_kwds = {}

        for kwds in kwds_tuple:
            new_kwds.update(kwds)

        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            values = [v for v in [kwds.get(name, None) for kwds in kwds_tuple] if v is not None]
            if values:
                new_kwds[name] = getattr(builtins, name[1:])(values)
        return new_kwds

    def __getitem__(self, i):
        """
        Get a frame from an animation or
        slice this animation returning a subanimation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ....:               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: print(a)
            Animation with 10 frames
            sage: a.show()          # long time  # optional -- ImageMagick
            sage: frame2 = a[2]     # indirect doctest
            sage: frame2.show()
            sage: print(a[3:7])     # indirect doctest
            Animation with 4 frames
            sage: a[3:7].show()     # long time  # optional -- ImageMagick
        """
        if isinstance(i, slice):
            return Animation(self._frames[i], **self._kwds)
        else:
            return self._frames[i]

    def _repr_(self):
        """
        Print representation for an animation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ....:               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: print(a)
            Animation with 10 frames
            sage: a._repr_()
            'Animation with 10 frames'
        """
        try:
            num = len(self)
        except TypeError:
            num = "unknown number of"
        return "Animation with %s frames" % num

    def __add__(self, other):
        """
        Add two animations. This has the effect of superimposing the two
        animations frame-by-frame.

        EXAMPLES::

            sage: a = animate([line([(0,0),(1,i)],hue=0/3) for i in range(2)],
            ....:             xmin=0, ymin=0, xmax=1, ymax=1,
            ....:             figsize=[1,1], axes=False)
            sage: print(a)
            Animation with 2 frames
            sage: a.show()        # optional -- ImageMagick
            sage: b = animate([line([(0,0),(i,1)],hue=2/3) for i in range(2)],
            ....:             xmin=0, ymin=0, xmax=1, ymax=1,
            ....:             figsize=[1,1], axes=False)
            sage: print(b)
            Animation with 2 frames
            sage: b.show()        # optional -- ImageMagick
            sage: s = a+b         # indirect doctest
            sage: print(s)
            Animation with 2 frames
            sage: len(a), len(b)
            (2, 2)
            sage: len(s)
            2
            sage: s.show()        # optional -- ImageMagick
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self._kwds, other._kwds)

        #Combine the frames
        m = max(len(self), len(other))
        frames = [a+b for a,b in zip(self._frames, other._frames)]
        frames += self._frames[m:] + other._frames[m:]

        return Animation(frames, **kwds)

    def __mul__(self, other):
        """
        Multiply two animations. This has the effect of appending the two
        animations (the second comes after the first).

        EXAMPLES::

            sage: a = animate([line([(0,0),(1,i)],hue=0/3) for i in range(2)],
            ....:             xmin=0, ymin=0, xmax=1, ymax=1,
            ....:             figsize=[1,1], axes=False)
            sage: print(a)
            Animation with 2 frames
            sage: a.show()             # optional -- ImageMagick
            sage: b = animate([line([(0,0),(i,1)],hue=2/3) for i in range(2)],
            ....:             xmin=0, ymin=0, xmax=1, ymax=1,
            ....:             figsize=[1,1], axes=False)
            sage: print(b)
            Animation with 2 frames
            sage: b.show()             # optional -- ImageMagick
            sage: p = a*b              # indirect doctest
            sage: len(a), len(b)
            (2, 2)
            sage: len(p)
            4
            sage: print(p)
            Animation with 4 frames
            sage: p.show()             # optional -- ImageMagick
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self._kwds, other._kwds)

        return Animation(self._frames + other._frames, **kwds)

    def __len__(self):
        """
        Length of ``self``.

        EXAMPLES::

            sage: a = animate([circle((i,0),1,thickness=20*i) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1], axes=False)
            sage: len(a)
            5
        """
        try:
            return self._num_frames
        except AttributeError:
            return len(self._frames)

    def make_image(self, frame, filename, **kwds):
        r"""
        Given a frame which has no ``save_image()`` method, make a graphics
        object and save it as an image with the given filename.  By default, this is
        :meth:`sage.plot.plot.plot`.  To make animations of other objects,
        override this method in a subclass.

        EXAMPLES::

            sage: from sage.plot.animate import Animation
            sage: class MyAnimation(Animation):
            ....:    def make_image(self, frame, filename, **kwds):
            ....:        P = parametric_plot(frame[0], frame[1], **frame[2])
            ....:        P.save_image(filename, **kwds)

            sage: t = SR.var("t")
            sage: x = lambda t: cos(t)
            sage: y = lambda n,t: sin(t)/n
            sage: B = MyAnimation([([x(t), y(i+1,t)], (t,0,1),
            ....:                   {'color':Color((1,0,i/4)), 'aspect_ratio':1, 'ymax':1})
            ....:                  for i in range(4)])

            sage: d = B.png(); v = os.listdir(d); v.sort(); v  # long time
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
            sage: B.show()  # not tested

            sage: class MyAnimation(Animation):
            ....:    def make_image(self, frame, filename, **kwds):
            ....:        G = frame.plot()
            ....:        G.set_axes_range(floor(G.xmin()), ceil(G.xmax()),
            ....:                         floor(G.ymin()), ceil(G.ymax()))
            ....:        G.save_image(filename, **kwds)

            sage: B = MyAnimation([graphs.CompleteGraph(n)
            ....:                  for n in range(7,11)], figsize=5)
            sage: d = B.png()
            sage: v = os.listdir(d); v.sort(); v
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
            sage: B.show()  # not tested
        """
        p = plot.plot(frame, **kwds)
        p.save_image(filename)

    def png(self, dir=None):
        r"""
        Render PNG images of the frames in this animation, saving them
        in ``dir``.  Return the absolute path to that directory.  If
        the frames have been previously rendered and ``dir`` is
        ``None``, just return the directory in which they are stored.

        When ``dir`` is other than ``None``, force re-rendering of
        frames.

        INPUT:

        - ``dir`` -- (default: ``None``) directory in which to store frames; in
          this case, a temporary directory will be created for storing the frames

        OUTPUT: absolute path to the directory containing the PNG images

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([plot(x^2 + n) for n in range(4)], ymin=0, ymax=4)
            sage: d = a.png(); v = os.listdir(d); v.sort(); v  # long time
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
        """
        if dir is None:
            try:
                return self._png_dir
            except AttributeError:
                pass
            dir = tmp_dir()
        i = 0
        for frame in self._frames:
            filename = '%s/%08d.png' % (dir,i)
            try:
                save_image = frame.save_image
            except AttributeError:
                self.make_image(frame, filename, **self._kwds)
            else:
                save_image(filename, **self._kwds)
            i += 1
        self._num_frames = i
        self._png_dir = dir
        return dir

    def graphics_array(self, ncols=3):
        r"""
        Return a :class:`sage.plot.multigraphics.GraphicsArray` with plots of the
        frames of this animation, using the given number of columns.
        The frames must be acceptable inputs for
        :class:`sage.plot.multigraphics.GraphicsArray`.

        EXAMPLES::

            sage: # needs sage.schemes
            sage: E = EllipticCurve('37a')
            sage: v = [E.change_ring(GF(p)).plot(pointsize=30)
            ....:      for p in [97, 101, 103]]
            sage: a = animate(v, xmin=0, ymin=0, axes=False)
            sage: print(a)
            Animation with 3 frames
            sage: a.show()                      # optional -- ImageMagick

        Modify the default arrangement of array::

            sage: g = a.graphics_array(); print(g)                                      # needs sage.schemes
            Graphics Array of size 1 x 3
            sage: g.show(figsize=[6,3])                                                 # needs sage.schemes

        Specify different arrangement of array and save it with a given file name::

            sage: g = a.graphics_array(ncols=2); print(g)                               # needs sage.schemes
            Graphics Array of size 2 x 2
            sage: f = tmp_filename(ext='.png'); print(f)                                # needs sage.schemes
            ...png
            sage: g.save(f)                                                             # needs sage.schemes

        Frames can be specified as a generator too; it is internally converted to a list::

            sage: t = SR.var("t")
            sage: b = animate((plot(sin(c*pi*t)) for c in sxrange(1,2,.2)))
            sage: g = b.graphics_array()
            sage: print(g)
            Graphics Array of size 2 x 3
        """
        ncols = int(ncols)
        frame_list = list(self._frames)
        n = len(frame_list)
        nrows, rem = divmod(n,ncols)
        if rem > 0:
            nrows += 1
        return plot.graphics_array(frame_list, nrows, ncols)

    def gif(self, delay=20, savefile=None, iterations=0, show_path=False,
            use_ffmpeg=False):
        r"""
        Return an animated gif composed from rendering the graphics
        objects in ``self``.

        This method will only work if either (a) the ImageMagick
        software suite is installed, i.e., you have the ``magick/convert``
        command or (b) ``ffmpeg`` is installed.  See the web sites of
        ImageMagick_ and FFmpeg_ for more details.  By default, this
        produces the gif using Imagemagick if it is present.  If this
        can't find ImageMagick or if ``use_ffmpeg`` is True, then it
        uses ``ffmpeg`` instead.

        INPUT:

        - ``delay`` -- (default: 20) delay in hundredths of a
          second between frames

        - ``savefile`` -- file that the animated gif gets saved
          to

        - ``iterations`` -- integer (default: 0); number of
          iterations of animation. If 0, loop forever.

        - ``show_path`` -- boolean (default: ``False``); if True,
          print the path to the saved file

        - ``use_ffmpeg`` -- boolean (default: ``False``); if ``True``, use
          'ffmpeg' by default instead of ImageMagick

        If ``savefile`` is not specified: in notebook mode, display the
        animation; otherwise, save it to a default file name.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([sin(x + float(k))
            ....:              for k in srange(0,2*pi,0.7)],
            ....:             xmin=0, xmax=2*pi, ymin=-1, ymax=1, figsize=[2,1])
            sage: td = tmp_dir()
            sage: a.gif()              # not tested
            sage: a.gif(savefile=td + 'my_animation.gif',               # long time  # optional -- ImageMagick
            ....:       delay=35, iterations=3)
            sage: with open(td + 'my_animation.gif', 'rb') as f:        # long time  # optional -- ImageMagick
            ....:     print(b'GIF8' in f.read())
            True
            sage: a.gif(savefile=td + 'my_animation.gif',               # long time  # optional -- ImageMagick
            ....:       show_path=True)
            Animation saved to .../my_animation.gif.
            sage: a.gif(savefile=td + 'my_animation_2.gif',             # long time  # optional -- FFmpeg
            ....:       show_path=True, use_ffmpeg=True)
            Animation saved to .../my_animation_2.gif.

        .. NOTE::

           If neither ffmpeg nor ImageMagick is installed, you will
           get an error message like this::

              Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
              animation to a GIF file or displaying an animation requires one of these
              packages, so please install one of them and try again.

              See www.imagemagick.org and www.ffmpeg.org for more information.
        """
        from sage.features.imagemagick import ImageMagick
        from sage.features.ffmpeg import FFmpeg

        if not ImageMagick().is_present() and not FFmpeg().is_present():
            raise OSError("Error: Neither ImageMagick nor ffmpeg appear to "
                    "be installed. Saving an animation to a GIF file or "
                    "displaying an animation requires one of these "
                    "packages, so please install one of them and try "
                    "again. See www.imagemagick.org and www.ffmpeg.org "
                    "for more information.")

        if use_ffmpeg or not ImageMagick().is_present():
            self.ffmpeg(savefile=savefile, show_path=show_path,
                        output_format='.gif', delay=delay,
                        iterations=iterations)
        else:
            self._gif_from_imagemagick(savefile=savefile, show_path=show_path,
                        delay=delay, iterations=iterations)

    def _gif_from_imagemagick(self, savefile=None, show_path=False,
            delay=20, iterations=0):
        r"""
        Return a movie showing an animation composed from rendering
        the frames in ``self``.

        This method will only work if ``imagemagick`` is installed (command
        ``magick`` or ``convert``). See https://www.imagemagick.org for information
        about ``imagemagick``.

        INPUT:

        - ``savefile`` -- file that the mpeg gets saved to

        .. warning::

            This will overwrite ``savefile`` if it already exists.

        - ``show_path`` -- boolean (default: ``False``); if ``True``,
          print the path to the saved file

        - ``delay`` -- (default: 20) delay in hundredths of a
          second between frames

        - ``iterations`` -- integer (default: 0); number of iterations
          of animation. If 0, loop forever.

        If ``savefile`` is not specified: in notebook mode, display
        the animation; otherwise, save it to a default file name.  Use
        :func:`sage.misc.verbose.set_verbose` with ``level=1`` to see
        additional output.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:             xmin=0, xmax=2*pi, ymin=-1, ymax=1, figsize=[2,1])
            sage: td = tmp_dir()
            sage: a._gif_from_imagemagick(savefile=td + 'new.gif')  # long time  # optional -- ImageMagick

        .. NOTE::

           If imagemagick is not installed, you will get an error message
           like this::

              FeatureNotPresentError: imagemagick is not available.
              Executable 'magick' not found on PATH.
              Further installation instructions might be available at
              https://www.imagemagick.org/.
        """
        from sage.features.imagemagick import ImageMagick, Magick
        ImageMagick().require()

        if not savefile:
            savefile = tmp_filename(ext='.gif')
        if not savefile.endswith('.gif'):
            savefile += '.gif'
        savefile = os.path.abspath(savefile)

        # running the command
        directory = self.png()
        cmd = [Magick().executable, '-dispose', 'Background',
                '-delay', '%s' % int(delay), '-loop', '%s' % int(iterations),
                '*.png', savefile]
        from subprocess import run
        result = run(cmd, cwd=directory, capture_output=True, text=True)

        # If a problem with the command occurs, print the log before
        # raising an error (more verbose than result.check_returncode())
        if result.returncode:
            print('Command "{}" returned nonzero exit status "{}" '
                  '(with stderr "{}" and stdout "{}").'.format(result.args,
                                        result.returncode,
                                        result.stderr.strip(),
                                        result.stdout.strip()))
            raise OSError("Error: Cannot generate GIF animation. "
                    "The magick/convert command (ImageMagick) is present but does "
                    "not seem to be functional. Verify that the objects "
                    "passed to the animate command can be saved in PNG "
                    "image format. "
                    "See www.imagemagick.org more information.")

        if show_path:
            print("Animation saved to file %s." % savefile)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([plot(x^2 + n) for n in range(2)], ymin=0, ymax=4)
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: a._rich_repr_(dm)       # optional -- ImageMagick
            OutputImageGif container
        """

        iterations = kwds.get('iterations', 0)
        loop = (iterations == 0)

        t = display_manager.types
        supported = display_manager.supported_output()
        format = kwds.pop("format", None)
        if format is None:
            if t.OutputImageGif in supported:
                format = "gif"
            else:
                return # No supported format could be guessed
        suffix = None
        outputType = None
        if format == "gif":
            outputType = t.OutputImageGif
            suffix = ".gif"
        if format == "ogg":
            outputType = t.OutputVideoOgg
        if format == "webm":
            outputType = t.OutputVideoWebM
        if format == "mp4":
            outputType = t.OutputVideoMp4
        if format == "flash":
            outputType = t.OutputVideoFlash
        if format == "matroska":
            outputType = t.OutputVideoMatroska
        if format == "avi":
            outputType = t.OutputVideoAvi
        if format == "wmv":
            outputType = t.OutputVideoWmv
        if format == "quicktime":
            outputType = t.OutputVideoQuicktime
        if format is None:
            raise ValueError("Unknown video format")
        if outputType not in supported:
            return # Sorry, requested format is not supported
        if suffix is not None:
            return display_manager.graphics_from_save(
                self.save, kwds, suffix, outputType)

        # Now we save for OutputVideoBase
        filename = tmp_filename(ext=outputType.ext)
        self.save(filename, **kwds)
        from sage.repl.rich_output.buffer import OutputBuffer
        buf = OutputBuffer.from_file(filename)
        return outputType(buf, loop=loop)

    def show(self, delay=None, iterations=None, **kwds):
        r"""
        Show this animation immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        INPUT:

        - ``delay`` -- (default: 20) delay in hundredths of a
          second between frames

        - ``iterations`` -- integer (default: 0); number of
          iterations of animation. If 0, loop forever.

        - ``format`` -- (default: gif) format to use for output.
          Currently supported formats are: gif,
          ogg, webm, mp4, flash, matroska, avi, wmv, quicktime.

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        .. NOTE::

           Currently this is done using an animated gif, though this
           could change in the future. This requires that either
           ffmpeg or the ImageMagick suite (in particular, the
           ``magick/convert`` command) is installed.

        See also the :meth:`ffmpeg` method.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([sin(x + float(k))
            ....:              for k in srange(0,2*pi,0.7)],
            ....:             xmin=0, xmax=2*pi, figsize=[2,1])
            sage: a.show()                              # long time  # optional -- ImageMagick

        The preceding will loop the animation forever. If you want to show
        only three iterations instead::

            sage: a.show(iterations=3)                  # long time  # optional -- ImageMagick

        To put a half-second delay between frames::

            sage: a.show(delay=50)                      # long time  # optional -- ImageMagick

        You can also make use of the HTML5 video element in the Sage notebook::

            sage: # long time, optional -- FFmpeg
            sage: a.show(format='ogg')
            sage: a.show(format='webm')
            sage: a.show(format='mp4')
            sage: a.show(format='webm', iterations=1)

        Other backends may support other file formats as well::

            sage: # long time, optional -- FFmpeg
            sage: a.show(format='flash')
            sage: a.show(format='matroska')
            sage: a.show(format='avi')
            sage: a.show(format='wmv')
            sage: a.show(format='quicktime')

        TESTS:

        Use of positional parameters is discouraged, will likely get
        deprecated, but should still work for the time being::

            sage: a.show(50, 3)                         # long time  # optional -- ImageMagick

        .. NOTE::

           If you don't have ffmpeg or ImageMagick installed, you will
           get an error message like this::

              Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
              animation to a GIF file or displaying an animation requires one of these
              packages, so please install one of them and try again.

              See www.imagemagick.org and www.ffmpeg.org for more information.
        """

        # Positional parameters for the sake of backwards compatibility
        if delay is not None:
            kwds.setdefault("delay", delay)
        if iterations is not None:
            kwds.setdefault("iterations", iterations)

        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, **kwds)

    def ffmpeg(self, savefile=None, show_path=False, output_format=None,
               ffmpeg_options='', delay=None, iterations=0, pix_fmt='rgb24'):
        r"""
        Return a movie showing an animation composed from rendering
        the frames in ``self``.

        This method will only work if ``ffmpeg`` is installed.  See
        https://www.ffmpeg.org for information about ``ffmpeg``.

        INPUT:

        - ``savefile`` -- file that the mpeg gets saved to

        .. warning::

            This will overwrite ``savefile`` if it already exists.

        - ``show_path`` -- boolean (default: ``False``); if ``True``,
          print the path to the saved file

        - ``output_format`` -- string (default: ``None``); format and
          suffix to use for the video.  This may be ``'mpg'``, ``'mpeg'``,
          ``'avi'``, ``'gif'``, or any other format that ``ffmpeg`` can handle.
          If this is ``None`` and the user specifies ``savefile`` with a
          suffix, say ``savefile='animation.avi'``, try to determine the
          format (``'avi'`` in this case) from that file name.  If no file
          is specified or if the suffix cannot be determined, ``'mpg'`` is
          used.

        - ``ffmpeg_options`` -- string (default: ``''``); this string is
          passed directly to ffmpeg

        - ``delay`` -- integer (default: ``None``); delay in hundredths of a
          second between frames.  The framerate is 100/delay.
          This is not supported for mpeg files: for mpegs, the frame
          rate is always 25 fps.

        - ``iterations`` -- integer (default: 0); number of iterations
          of animation. If 0, loop forever.  This is only supported
          for animated gif output and requires ``ffmpeg`` version 0.9 or
          later.  For older versions, set ``iterations=None``.

        - ``pix_fmt`` -- string (default: ``'rgb24'``); used only for gif
          output.  Different values such as 'rgb8' or 'pal8' may be
          necessary depending on how ffmpeg was installed.  Set
          ``pix_fmt=None`` to disable this option.

        If ``savefile`` is not specified: in notebook mode, display
        the animation; otherwise, save it to a default file name.  Use
        :func:`sage.misc.verbose.set_verbose` with ``level=1`` to see
        additional output.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([sin(x + float(k))
            ....:              for k in srange(0, 2*pi, 0.7)],
            ....:             xmin=0, xmax=2*pi, ymin=-1, ymax=1, figsize=[2,1])
            sage: td = tmp_dir()
            sage: a.ffmpeg(savefile=td + 'new.mpg')                  # long time  # optional -- FFmpeg
            sage: a.ffmpeg(savefile=td + 'new.avi')                  # long time  # optional -- FFmpeg
            sage: a.ffmpeg(savefile=td + 'new.gif')                  # long time  # optional -- FFmpeg
            sage: a.ffmpeg(savefile=td + 'new.mpg', show_path=True)  # long time  # optional -- FFmpeg
            Animation saved to .../new.mpg.

        .. NOTE::

           If ffmpeg is not installed, you will get an error message
           like this::

              FeatureNotPresentError: ffmpeg is not available.
              Executable 'ffmpeg' not found on PATH.
              Further installation instructions might be available at https://www.ffmpeg.org/.

        TESTS::

            sage: a.ffmpeg(output_format='gif',delay=30,iterations=5)  # long time  # optional -- FFmpeg
        """
        from sage.features.ffmpeg import FFmpeg
        FFmpeg().require()

        if savefile is None:
            if output_format is None:
                output_format = '.mpg'
            else:
                if output_format[0] != '.':
                    output_format = '.'+output_format
            savefile = tmp_filename(ext=output_format)
        else:
            if output_format is None:
                suffix = os.path.splitext(savefile)[1]
                if len(suffix) > 0:
                    output_format = suffix
                else:
                    output_format = '.mpg'
        if not savefile.endswith(output_format):
            savefile += output_format
        early_options = ''
        if output_format == '.gif':
            # We try to set reasonable options for gif output.
            #
            # Older versions of ffmpeg (before 0.9, summer 2011)
            # use the option -loop_output instead of -loop.
            # Setting iterations=None is a way of preventing sage
            # from adding the -loop option.  A separate
            # -loop_output option can be added with the
            # ffmpeg_options argument.
            if iterations is not None:
                loop_cmd = f'-loop {iterations} '
            else:
                loop_cmd = ''
            # A pix_fmt value is required for some but not all
            # ffmpeg installations.  Setting pix_fmt=None will
            # prevent sage from adding this option, and it may be
            # controlled separately through ffmpeg_options.
            if pix_fmt is not None:
                pix_fmt_cmd = f'-pix_fmt {pix_fmt} '
            else:
                pix_fmt_cmd = ''
            ffmpeg_options += f' {pix_fmt_cmd}{loop_cmd}'
        if delay is not None and output_format != '.mpeg' and output_format != '.mpg':
            early_options += ' -r %s ' % int(100/delay)
        savefile = os.path.abspath(savefile)
        pngdir = self.png()
        pngs = os.path.join(pngdir, "%08d.png")
        # For ffmpeg, it seems that some options, like '-g ... -r
        # ...', need to come before the input file names, while
        # some options, like '-pix_fmt rgb24', need to come
        # afterwards.  Hence 'early_options' and 'ffmpeg_options'
        # The `-nostdin` is needed to avoid the command to hang, see
        # https://stackoverflow.com/questions/16523746/ffmpeg-hangs-when-run-in-background
        cmd = 'cd {}; {} -nostdin -y -f image2 {} -i {} {} {}'.format(
            shlex.quote(pngdir), shlex.quote(FFmpeg().absolute_filename()),
            early_options, shlex.quote(pngs), ffmpeg_options, shlex.quote(savefile))
        from subprocess import check_call, CalledProcessError, PIPE
        try:
            if sage.misc.verbose.get_verbose() > 0:
                set_stderr = None
            else:
                set_stderr = PIPE
            sage.misc.verbose.verbose("Executing '%s'" % cmd,level=1)
            sage.misc.verbose.verbose("\n---- ffmpeg output below ----\n")
            check_call(cmd, shell=True, stderr=set_stderr)
            if show_path:
                print("Animation saved to file %s." % savefile)
        except (CalledProcessError, OSError):
            print("Error running ffmpeg.")
            raise

    def apng(self, savefile=None, show_path=False, delay=20, iterations=0):
        r"""
        Create an animated PNG composed from rendering the graphics
        objects in ``self``. Return the absolute path to that file.

        Notice that not all web browsers are capable of displaying APNG
        files, though they should still present the first frame of the
        animation as a fallback.

        The generated file is not optimized, so it may be quite large.

        Input:

        - ``delay`` -- (default: 20) delay in hundredths of a
          second between frames

        - ``savefile`` -- file that the animated gif gets saved
          to

        - ``iterations`` -- integer (default: 0); number of
          iterations of animation. If 0, loop forever.

        - ``show_path`` -- boolean (default: ``False``); if True,
          print the path to the saved file

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([sin(x + float(k))
            ....:              for k in srange(0,2*pi,0.7)],
            ....:             xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir()
            sage: a.apng(show_path=True)  # long time
            Animation saved to ....png.
            sage: a.apng(savefile=dir + 'my_animation.png', delay=35, iterations=3)  # long time
            sage: a.apng(savefile=dir + 'my_animation.png', show_path=True)  # long time
            Animation saved to .../my_animation.png.

        If the individual frames have different sizes, an error will be raised::

            sage: a = animate([plot(sin(x), (x, 0, k))
            ....:              for k in range(1,4)],
            ....:             ymin=-1, ymax=1, aspect_ratio=1, figsize=[2,1])
            sage: a.apng()  # long time
            Traceback (most recent call last):
            ...
            ValueError: Chunk IHDR mismatch

        TESTS::

            sage: a = animate([])
            sage: a.apng(show_path=True)
            Animation saved to file ....png.
        """
        pngdir = self.png()
        if savefile is None:
            savefile = tmp_filename(ext='.png')
        with open(savefile, "wb") as out:
            apng = APngAssembler(
                out, len(self),
                delay=delay, num_plays=iterations)
            for i in range(len(self)):
                png = os.path.join(pngdir, "%08d.png" % i)
                apng.add_frame(png)
        if show_path:
            print("Animation saved to file %s." % savefile)

    def save(self, filename=None, show_path=False, use_ffmpeg=False, **kwds):
        r"""
        Save this animation.

        INPUT:

        - ``filename`` -- (default: ``None``) name of save file

        - ``show_path`` -- boolean (default: ``False``); if True,
          print the path to the saved file

        - ``use_ffmpeg`` -- boolean (default: ``False``); if ``True``, use
          'ffmpeg' by default instead of ImageMagick when creating GIF files

        If filename is ``None``, then in notebook mode, display the
        animation; otherwise, save the animation to a GIF file. If
        filename ends in '.html', save an :meth:`interactive` version of
        the animation to an HTML file that uses the Three.js viewer.  If
        filename ends in '.sobj', save to an sobj file.  Otherwise,
        try to determine the format from the filename extension
        ('.mpg', '.gif', '.avi', etc.).  If the format cannot be
        determined, default to GIF.

        For GIF files, either ffmpeg or the ImageMagick suite must be
        installed.  For other movie formats, ffmpeg must be installed.
        sobj and HTML files can be saved with no extra software installed.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: a = animate([sin(x + float(k))
            ....:              for k in srange(0, 2*pi, 0.7)],
            ....:             xmin=0, xmax=2*pi, ymin=-1, ymax=1, figsize=[2,1])
            sage: td = tmp_dir()
            sage: a.save()         # not tested
            sage: a.save(td + 'wave.gif')                   # long time  # optional -- ImageMagick
            sage: a.save(td + 'wave.gif', show_path=True)   # long time  # optional -- ImageMagick
            Animation saved to file .../wave.gif.
            sage: a.save(td + 'wave.avi', show_path=True)   # long time  # optional -- FFmpeg
            Animation saved to file .../wave.avi.
            sage: a.save(td + 'wave0.sobj')
            sage: a.save(td + 'wave1.sobj', show_path=True)
            Animation saved to file .../wave1.sobj.
            sage: a.save(td + 'wave0.html', online=True)
            sage: a.save(td + 'wave1.html', show_path=True, online=True)
            Animation saved to file .../wave1.html.

        TESTS:

        Ensure that we can pass delay and iteration count to the saved
        GIF image (see :issue:`18176`)::

            sage: # long time, optional -- ImageMagick
            sage: a.save(td + 'wave.gif')
            sage: with open(td + 'wave.gif', 'rb') as f:
            ....:     print(b'GIF8' in f.read())
            True
            sage: with open(td + 'wave.gif', 'rb') as f:
            ....:     print(b'!\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00' in f.read())
            True
            sage: a.save(td + 'wave.gif', delay=35)
            sage: with open(td + 'wave.gif', 'rb') as f:
            ....:     print(b'GIF8' in f.read())
            True
            sage: a.save(td + 'wave.gif', iterations=3)
            sage: with open(td + 'wave.gif', 'rb') as f:
            ....:     print(b'!\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00' in f.read())
            False
            sage: with open(td + 'wave.gif', 'rb') as f:
            ....:      check1 = b'!\xff\x0bNETSCAPE2.0\x03\x01\x02\x00\x00'
            ....:      check2 = b'!\xff\x0bNETSCAPE2.0\x03\x01\x03\x00\x00'
            ....:      data = f.read()
            ....:      print(check1 in data or check2 in data)
            True
        """
        if filename is None:
            suffix = '.gif'
        else:
            suffix = os.path.splitext(filename)[1]
            if not suffix:
                suffix = '.gif'

        if filename is None or suffix == '.gif':
            self.gif(savefile=filename, show_path=show_path,
                     use_ffmpeg=use_ffmpeg, **kwds)
        elif suffix == '.sobj':
            SageObject.save(self, filename)
            if show_path:
                print("Animation saved to file %s." % filename)
        elif suffix == '.html':
            self.interactive(**kwds).save(filename)
            if show_path:
                print("Animation saved to file %s." % filename)
        else:
            self.ffmpeg(savefile=filename, show_path=show_path, **kwds)

    def interactive(self, **kwds):
        r"""
        Create an interactive depiction of the animation.

        INPUT:

        - ``**kwds`` -- any of the viewing options accepted by show() are valid
          as keyword arguments to this function and they will behave in the same
          way. Those that are animation-related and recognized by the Three.js
          viewer are: ``animate``, ``animation_controls``, ``auto_play``,
          ``delay``, and ``loop``.

        OUTPUT: a 3D graphics object which, by default, will use the Three.js viewer

        EXAMPLES::

            sage: x = SR.var("x")
            sage: frames = [point3d((sin(x), cos(x), x))
            ....:           for x in (0, pi/16, .., 2*pi)]
            sage: animate(frames).interactive(online=True)
            Graphics3d Object

        Works with frames that are 2D or 3D graphics objects or convertible to
        2D or 3D graphics objects via a ``plot`` or ``plot3d`` method::

            sage: frames = [dodecahedron(), circle(center=(0, 0), radius=1), x^2]
            sage: animate(frames).interactive(online=True, delay=100)
            Graphics3d Object

        .. SEEALSO::

            :ref:`threejs_viewer`
        """
        from sage.plot.plot3d.base import Graphics3d, KeyframeAnimationGroup
        # Attempt to convert frames to Graphics3d objects.
        g3d_frames = []
        for i, frame in enumerate(self._frames):
            if not isinstance(frame, Graphics3d):
                try:
                    frame = frame.plot3d()
                except (AttributeError, TypeError):
                    try:
                        frame = frame.plot().plot3d()
                    except (AttributeError, TypeError):
                        frame = None
                if not isinstance(frame, Graphics3d):
                    raise TypeError(f"Could not convert frame {i} to Graphics3d")
            g3d_frames.append(frame)
        # Give preference to this method's keyword arguments over those provided
        # to animate or the constructor.
        kwds = dict(self._kwds, **kwds)
        # Three.js is the only viewer that supports animation at present.
        if 'viewer' not in kwds:
            kwds['viewer'] = 'threejs'
        return KeyframeAnimationGroup(g3d_frames, **kwds)


class APngAssembler:
    r"""
    Build an APNG_ (Animated PNG) from a sequence of PNG files.
    This is used by the :meth:`sage.plot.animate.Animation.apng` method.

    This code is quite simple; it does little more than copying chunks
    from input PNG files to the output file. There is no optimization
    involved. This does not depend on external programs or libraries.

    INPUT:

    - ``out`` -- a file opened for binary writing to which the data
      will be written

    - ``num_frames`` -- the number of frames in the animation

    - ``num_plays`` -- how often to iterate, 0 means infinitely

    - ``delay`` -- numerator of the delay fraction in seconds

    - ``delay_denominator`` -- denominator of the delay in seconds

    EXAMPLES::

        sage: from sage.plot.animate import APngAssembler
        sage: x = SR.var("x")
        sage: def assembleAPNG():
        ....:     a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
        ....:                 xmin=0, xmax=2*pi, figsize=[2,1])
        ....:     pngdir = a.png()
        ....:     outfile = sage.misc.temporary_file.tmp_filename(ext='.png')
        ....:     with open(outfile, "wb") as f:
        ....:         apng = APngAssembler(f, len(a))
        ....:         for i in range(len(a)):
        ....:             png = os.path.join(pngdir, "{:08d}.png".format(i))
        ....:             apng.add_frame(png, delay=10*i + 10)
        ....:     return outfile
        sage: assembleAPNG()  # long time
        '...png'
    """
    magic = b"\x89PNG\x0d\x0a\x1a\x0a"
    mustmatch = frozenset([b"IHDR", b"PLTE", b"bKGD", b"cHRM", b"gAMA",
                           b"pHYs", b"sBIT", b"tRNS"])

    def __init__(self, out, num_frames,
                 num_plays=0, delay=200, delay_denominator=100):
        r"""
        Initialize for creation of an APNG file.
        """
        self._last_seqno = -1
        self._idx = 0
        self._first = True
        self.out = out
        self.num_frames = num_frames
        self.num_plays = num_plays
        self.default_delay_numerator = delay
        self.default_delay_denominator = delay_denominator
        self._matchref = {}
        self.out.write(self.magic)

    def add_frame(self, pngfile, delay=None, delay_denominator=None):
        r"""
        Add a single frame to the APNG file.

        INPUT:

        - ``pngfile`` -- file name of the PNG file with data for this frame

        - ``delay`` -- numerator of the delay fraction in seconds

        - ``delay_denominator`` -- denominator of the delay in seconds

        If the delay is not specified, the default from the constructor
        applies.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from io import BytesIO
            sage: buf = BytesIO()
            sage: apng = APngAssembler(buf, 2)
            sage: fn = APngAssembler._testData("input1", True)
            sage: apng.add_frame(fn, delay=0x567, delay_denominator=0x1234)
            sage: fn = APngAssembler._testData("input2", True)
            sage: apng.add_frame(fn)
            sage: len(buf.getvalue())
            217
            sage: buf.getvalue() == APngAssembler._testData("anim12", False)
            True
            sage: apng.add_frame(fn)
            Traceback (most recent call last):
            ...
            RuntimeError: Already reached the declared number of frames
        """
        if self._idx == self.num_frames:
            raise RuntimeError("Already reached the declared number of frames")
        self.delay_numerator = self.default_delay_numerator
        self.delay_denominator = self.default_delay_denominator
        self._actl_written = False
        self._fctl_written = False
        if delay is not None:
            self.delay_numerator = delay
        if delay_denominator is not None:
            self.delay_denominator = delay_denominator
        self._add_png(pngfile)
        self._idx += 1
        if self._idx == self.num_frames:
            self._chunk(b"IEND", b"")

    def set_default(self, pngfile):
        r"""
        Add a default image for the APNG file.

        This image is used as a fallback in case some application does
        not understand the APNG format.  This method must be called
        prior to any calls to the ``add_frame`` method, if it is called
        at all.  If it is not called, then the first frame of the
        animation will be the default.

        INPUT:

        - ``pngfile`` -- file name of the PNG file with data
          for the default image

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from io import BytesIO
            sage: buf = BytesIO()
            sage: apng = APngAssembler(buf, 1)
            sage: fn = APngAssembler._testData("input1", True)
            sage: apng.set_default(fn)
            sage: fn = APngAssembler._testData("input2", True)
            sage: apng.add_frame(fn, delay=0x567, delay_denominator=0x1234)
            sage: len(buf.getvalue())
            179
            sage: buf.getvalue() == APngAssembler._testData("still1anim2", False)
            True
            sage: apng.add_frame(fn)
            Traceback (most recent call last):
            ...
            RuntimeError: Already reached the declared number of frames
        """
        if self._idx != 0:
            raise RuntimeError("Default image must precede all animation frames")
        self._actl_written = False
        self._fctl_written = True
        self._add_png(pngfile)

    def _add_png(self, pngfile):
        r"""
        Add data from one PNG still image.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_add_png", reads=False)
            enter _add_png('...png')
              write _current_chunk = (...'\x00\x00\x00\r',...'IHDR',...'\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00',...'\xb8\x1f9\xc6')
              call _copy() -> None
              call _first_IHDR(...'\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00') -> None
              write _current_chunk = (...'\x00\x00\x00\x04',...'gAMA',...'\x00\x01\x86\xa0',...'1\xe8\x96_')
              call _copy() -> None
              write _current_chunk = (...'\x00\x00\x00\x07',...'tIME',...'\x07\xde\x06\x1b\x0b&$',...'\x1f0z\xd5')
              write _current_chunk = (...'\x00\x00\x00\x08',...'IDAT',...'img1data',...'\xce\x8aI\x99')
              call _first_IDAT(...'img1data') -> None
              write _current_chunk = (...'\x00\x00\x00\x00',...'IEND',...'',...'\xaeB`\x82')
              write _first = False
            exit _add_png -> None
            enter _add_png('...png')
              write _current_chunk = (...'\x00\x00\x00\r',...'IHDR',...'\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00',...'\xb8\x1f9\xc6')
              write _current_chunk = (...'\x00\x00\x00\x04',...'gAMA',...'\x00\x01\x86\xa0',...'1\xe8\x96_')
              write _current_chunk = (...'\x00\x00\x00\x04',...'IDAT',...'img2',...'\x0ei\xab\x1d')
              call _next_IDAT(...'img2') -> None
              write _current_chunk = (...'\x00\x00\x00\x04',...'IDAT',...'data',...'f\x94\xcbx')
              call _next_IDAT(...'data') -> None
              write _current_chunk = (...'\x00\x00\x00\x00',...'IEND',...'',...'\xaeB`\x82')
              write _first = False
            exit _add_png -> None
        """
        with open(pngfile, 'rb') as png:
            if png.read(8) != self.magic:
                raise ValueError(f"{pngfile} is not a PNG file")
            while True:
                chead = png.read(8)
                if len(chead) == 0:
                    break
                clen, ctype = struct.unpack(">L4s", chead)
                cdata = png.read(clen)
                ccrc = png.read(4)
                utype = ctype.decode("ascii")
                self._current_chunk = (chead[:4], ctype, cdata, ccrc)
                if ctype in self.mustmatch:
                    ref = self._matchref.get(ctype)
                    if ref is None:
                        self._matchref[ctype] = cdata
                        self._copy()
                    else:
                        if cdata != ref:
                            raise ValueError(f"Chunk {utype} mismatch")
                met = ("_first_" if self._first else "_next_") + utype
                try:
                    met = getattr(self, met)
                except AttributeError:
                    pass
                else:
                    met(cdata)
        self._first = False

    def _seqno(self):
        r"""
        Generate next sequence number.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from io import BytesIO
            sage: buf = BytesIO()
            sage: apng = APngAssembler(buf, 1)
            sage: apng._seqno() == b'\x00\x00\x00\x00'
            True
            sage: apng._seqno() == b'\x00\x00\x00\x01'
            True
            sage: apng._seqno() == b'\x00\x00\x00\x02'
            True
        """
        self._last_seqno += 1
        return struct.pack(">L", self._last_seqno)

    def _first_IHDR(self, data):
        r"""
        Remember image size.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_first_IHDR")
            enter _first_IHDR(...'\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00')
              write width = 3
              write height = 2
            exit _first_IHDR -> None
        """
        w, h, d, ctype, comp, filt, ilace = struct.unpack(">2L5B", data)
        self.width = w
        self.height = h

    def _first_IDAT(self, data):
        r"""
        Write acTL and fcTL, then copy as IDAT.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_first_IDAT")
            enter _first_IDAT(...'img1data')
              call _actl() -> None
              call _fctl() -> None
              call _copy() -> None
            exit _first_IDAT -> None
        """
        self._actl()
        self._fctl()
        self._copy()

    def _next_IDAT(self, data):
        r"""
        Write fcTL, then convert to fdAT.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_next_IDAT")
            enter _next_IDAT(...'img2')
              call _fctl() -> None
              call _seqno() -> ...'\x00\x00\x00\x02'
              call _chunk(...'fdAT', ...'\x00\x00\x00\x02img2') -> None
            exit _next_IDAT -> None
            enter _next_IDAT(...'data')
              call _fctl() -> None
              call _seqno() -> ...'\x00\x00\x00\x03'
              call _chunk(...'fdAT', ...'\x00\x00\x00\x03data') -> None
            exit _next_IDAT -> None
        """
        self._fctl()
        maxlen = 0x7ffffffb
        while len(data) > maxlen:
            self._chunk(b"fdAT", self._seqno() + data[:maxlen])
            data = data[maxlen:]
        self._chunk(b"fdAT", self._seqno() + data)

    def _copy(self):
        r"""
        Copy an existing chunk without modification.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_copy")
            enter _copy()
              read _current_chunk = (...'\x00\x00\x00\r',...'IHDR',...'\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00',...'\xb8\x1f9\xc6')
              read out = <_io.BytesIO... at ...
              read out = <_io.BytesIO... at ...
              read out = <_io.BytesIO... at ...
              read out = <_io.BytesIO... at ...
            exit _copy -> None
            enter _copy()
              read _current_chunk = (...'\x00\x00\x00\x04',...'gAMA',...'\x00\x01\x86\xa0',...'1\xe8\x96_')
            ...
              read _current_chunk = (...'\x00\x00\x00\x08',...'IDAT',...'img1data',...'\xce\x8aI\x99')
            ...
            exit _copy -> None
        """
        for d in self._current_chunk:
            self.out.write(d)

    def _actl(self):
        r"""
        Write animation control data (acTL).

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_actl")
            enter _actl()
              read _actl_written = False
              read num_frames = 2
              read num_plays = 0
              call _chunk(...'acTL',...'\x00\x00\x00\x02\x00\x00\x00\x00') -> None
              write _actl_written = True
            exit _actl -> None
        """
        if self._actl_written:
            return
        data = struct.pack(">2L", self.num_frames, self.num_plays)
        self._chunk(b"acTL", data)
        self._actl_written = True

    def _fctl(self):
        r"""
        Write frame control data (fcTL).

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_fctl")
            enter _fctl()
              read _fctl_written = False
              read width = 3
              read height = 2
              read delay_numerator = 1383
              read delay_denominator = 4660
              call _seqno() -> ...'\x00\x00\x00\x00'
              call _chunk(...'fcTL', ...'\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x00\x05g\x124\x01\x00') -> None
              write _fctl_written = True
            exit _fctl -> None
            enter _fctl()
              read _fctl_written = False
              read width = 3
              read height = 2
              read delay_numerator = 200
              read delay_denominator = 100
              call _seqno() -> ...'\x00\x00\x00\x01'
              call _chunk(...'fcTL', ...'\x00\x00\x00\x01\x00\x00\x00\x03\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00\xc8\x00d\x01\x00') -> None
              write _fctl_written = True
            exit _fctl -> None
            enter _fctl()
              read _fctl_written = True
            exit _fctl -> None
        """
        if self._fctl_written:
            return
        data = struct.pack(
            ">4L2H2B",
            self.width, self.height, 0, 0,
            self.delay_numerator, self.delay_denominator,
            1, 0)
        self._chunk(b"fcTL", self._seqno() + data)
        self._fctl_written = True

    def _chunk(self, ctype, cdata):
        r"""
        Write a new (or modified) chunk of data.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from io import BytesIO
            sage: buf = BytesIO()
            sage: apng = APngAssembler(buf, 1)
            sage: buf.getvalue() == b'\x89PNG\r\n\x1a\n'
            True
            sage: apng._chunk(b"abcd", b"efgh")
            sage: buf.getvalue() == b'\x89PNG\r\n\x1a\n\x00\x00\x00\x04abcdefgh\xae\xef*P'
            True
        """
        ccrc = struct.pack(">L", zlib.crc32(ctype + cdata) & 0xffffffff)
        clen = struct.pack(">L", len(cdata))
        for d in [clen, ctype, cdata, ccrc]:
            self.out.write(d)

    @classmethod
    def _hex2bin(cls, h):
        r"""
        Convert hex data to binary.

        This is a helper method used for testing.
        Most data is given as lower-case hex digits,
        possibly intermixed with whitespace.
        A dot causes the next four bytes to be copied verbatim
        even if they look like hex digits. This is used for chunk types.
        Other characters which are not hex digits are passed verbatim.

        EXAMPLES::

            sage: from sage.plot.animate import APngAssembler
            sage: h2b = APngAssembler._hex2bin
            sage: h2b("0123") == b"\x01\x23"
            True
            sage: h2b(" 01 \n 23 ") == b"\x01\x23"
            True
            sage: h2b(".abcdef") == b"abcd\xef"
            True
            sage: h2b("PNG") == b"PNG"
            True
        """
        b = []
        while h:
            if h[0] in ' \n': # ignore whitespace
                h = h[1:]
            elif h[0] in '0123456789abcdef': # hex byte
                b.append(int(h[:2], 16))
                h = h[2:]
            elif h[0] == '.': # for chunk type
                b.extend(ord(h[i]) for i in range(1, 5))
                h = h[5:]
            else: # for PNG magic
                b.append(ord(h[0]))
                h = h[1:]

        return bytes(b)

    @classmethod
    def _testData(cls, name, asFile):
        r"""
        Retrieve data for test cases.

        INPUT:

        - ``name`` -- the name of the file content.

        - ``asFile`` -- whether to return a binary string of the named data
          or the path of a file containing that data

        EXAMPLES::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testData("input1", False).startswith(b'\x89PNG')
            True
            sage: APngAssembler._testData("input2", True)
            '...png'
        """
        data = {

            # Input 1: one PNG image, except the data makes no real sense
            "input1": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000007.tIME 07de061b0b2624 1f307ad5
            00000008.IDAT 696d673164617461 ce8a4999
            00000000.IEND ae426082""",

            # Input 2: slightly different, data in two chunks
            "input2": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000004.IDAT 696d6732 0e69ab1d
            00000004.IDAT 64617461 6694cb78
            00000000.IEND ae426082""",

            # Expected output 1: both images as frames of an animation
            "anim12": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000008.acTL 0000000200000000 f38d9370
            0000001a.fcTL 000000000000000300000002
                          0000000000000000056712340100 b4f729c9
            00000008.IDAT 696d673164617461 ce8a4999
            0000001a.fcTL 000000010000000300000002
                          000000000000000000c800640100 1b92eb4d
            00000008.fdAT 00000002696d6732 9cfb89a3
            00000008.fdAT 0000000364617461 c966c076
            00000000.IEND ae426082""",

            # Expected output 2: first image as fallback, second as animation
            "still1anim2": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000008.acTL 0000000100000000 b42de9a0
            00000008.IDAT 696d673164617461 ce8a4999
            0000001a.fcTL 000000000000000300000002
                          0000000000000000056712340100 b4f729c9
            00000008.fdAT 00000001696d6732 db5bf373
            00000008.fdAT 0000000264617461 f406e9c6
            00000000.IEND ae426082""",

        }
        d = cls._hex2bin(data[name])
        if asFile:
            from sage.misc.temporary_file import tmp_filename
            fn = tmp_filename(ext='.png')
            with open(fn, 'wb') as f:
                f.write(d)
            return fn
        return d

    @classmethod
    def _testCase1(cls, methodToTrace=None, **kwds):
        r"""
        Run common test case.

        This test case is one animation of two frames.
        The named method (if not None) will be traced during execution.
        This will demonstrate the role of each method in the doctests.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1()
        """
        from sage.doctest.fixtures import trace_method
        from io import BytesIO
        buf = BytesIO()
        apng = cls(buf, 2)
        if methodToTrace is not None:
            trace_method(apng, methodToTrace, **kwds)
        apng.add_frame(cls._testData("input1", True),
                       delay=0x567, delay_denominator=0x1234)
        apng.add_frame(cls._testData("input2", True))
        out = buf.getvalue()
        assert len(out) == 217
        assert out == cls._testData("anim12", False)
