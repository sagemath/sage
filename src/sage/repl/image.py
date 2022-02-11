# -*- encoding: utf-8 -*-
"""
Sage Wrapper for Bitmap Images

Some computations in Sage return bitmap images, for example matrices
can be turned into bitmaps directly. Note that this is different from
all plotting functionality, the latter can equally produce vector
graphics. This module is about bitmaps only, and a shallow wrapper
around ``PIL.Image``. The only difference is that :class:`Image`
is displayed as graphics by the Sage if the UI can.

EXAMPLES::

    sage: from sage.repl.image import Image
    sage: img = Image('RGB', (256, 256), 'white')
    sage: pixels = img.pixels()
    sage: for x in range(img.width()):
    ....:     for y in range(img.height()):
    ....:         pixels[x, y] = (x, y, 100)
    sage: img
    256x256px 24-bit RGB image
    sage: type(img)
    <class 'sage.repl.image.Image'>
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import io

import PIL.Image
from sage.structure.sage_object import SageObject


class Image(SageObject):

    def __init__(self, mode, size, color='white'):
        """
        Creates a new image with the given mode and size.

        INPUT:

        - ``mode`` -- string. The mode to use for the new image. Valid
          options are:

              * ``'1'`` (1-bit pixels, black and white, stored with
                one pixel per byte)

              * ``'L'`` (8-bit pixels, black and white)

              * ``'P'`` (8-bit pixels, mapped to any other mode using
                a color palette)

              * ``'RGB'`` (3x8-bit pixels, true color)

              * ``'RGBA'`` (4x8-bit pixels, true color with
                transparency mask)

              * ``'CMYK'`` (4x8-bit pixels, color separation)

              * ``'YCbCr'`` (3x8-bit pixels, color video format)

              * ``'LAB'`` (3x8-bit pixels, the L*a*b color space)

              * ``'HSV'`` (3x8-bit pixels, Hue, Saturation, Value
                color space)

              * ``'I'`` (32-bit signed integer pixels)

              * ``'F'`` (32-bit floating point pixels)

        - ``size`` -- 2-tuple, containing (width, height) in pixels.

        - ``color`` -- string, numeric or tuple of numeric. What colour to use
          for the image. Default is black.  If given, this should be a
          a tuple with one value per band. When creating RGB images,
          you can also use colour strings as supported by the
          ImageColor module.  If the colour is None, the image is not
          initialised.

        OUTPUT:

        A new :class:`Image` object.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: Image('P', (16, 16), 13)
            16x16px 8-bit Color image
        """
        # pillow does not support Sage integers as color
        from sage.rings.integer import Integer
        if isinstance(color, Integer):
            color = int(color)
        elif isinstance(color, tuple):
            color = tuple(int(i) if isinstance(i, Integer) else i for i in color)
        self._pil = PIL.Image.new(mode, size, color)

    @property
    def pil(self):
        """
        Access the wrapped PIL(low) Image

        OUTPUT:

        The underlying ``PIL.Image.Image object``.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('RGB', (16, 16), 'white')
            sage: img.pil
            <PIL.Image.Image image mode=RGB size=16x16 at 0x...>
        """
        return self._pil

    def pixels(self):
        """
        Return the pixel map

        OUTPUT:

        The PIL PixelAccess object that allows you to get/set the
        pixel data.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('RGB', (16, 16), 'white')
            sage: img.pixels()
            <PixelAccess object at 0x...>
        """
        return self._pil.load()

    def _repr_(self):
        """
        Return string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: Image('RGB', (16, 16), 'white')   # indirect doctest
            16x16px 24-bit RGB image
        """
        modestr = {
            '1':     '{0}x{1}px BW image',
            'L':     '{0}x{1}px 8-bit BW image',
            'P':     '{0}x{1}px 8-bit Color image',
            'RGB':   '{0}x{1}px 24-bit RGB image',
            'RGBA':  '{0}x{1}px 32-bit RGBA image',
            'CMYK':  '{0}x{1}px 24-bit CMYK image',
            'YCbCr': '{0}x{1}px 24-bit YCbCr mage',
            'LAB':   '{0}x{1}px 24-bit LAB image',
            'HSV':   '{0}x{1}px 24-bit HSV image',
            'I':     '{0}x{1}px 32-bit signed integer image',
            'F':     '{0}x{1}px 32-bit float image',
        }
        try:
            mode = modestr[self.pil.mode]
        except AttributeError:
            mode = 'Unknown mode'
        width, height = self.pil.size
        return mode.format(width, height)

    def mode(self):
        """
        Return the color mode

        OUTPUT:

        String. As given when constructing the image.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('YCbCr', (16, 16), 'white')
            sage: img.mode()
            'YCbCr'
        """
        return self.pil.mode

    def width(self):
        """
        Return the horizontal dimension in pixels

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('1', (12, 34), 'white')
            sage: img.width()
            12
            sage: img.height()
            34
        """
        return self.pil.size[0]

    def height(self):
        """
        Return the vertical dimension in pixels

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('1', (12, 34), 'white')
            sage: img.width()
            12
            sage: img.height()
            34
        """
        return self.pil.size[1]

    def save(self, filename):
        r"""
        Save the bitmap image

        INPUT:

        - ``filename`` -- string. The filename to save as. The given
          extension automatically determines the image file type.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('P', (12, 34), 13)
            sage: filename = tmp_filename(ext='.png')
            sage: img.save(filename)
            sage: with open(filename, 'rb') as f:
            ....:     f.read(4) == b'\x89PNG'
            True
        """
        self.pil.save(filename)

    def show(self):
        r"""
        Show this image immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('1', (12, 34), 'white')
            sage: img.show()
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.image import Image
            sage: img = Image('1', (16, 16), 'white')
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: img._rich_repr_(dm)
            OutputImagePng container

            sage: img = Image('F', (16, 16), 'white')   # not supported in PNG
            sage: img._rich_repr_(dm)
            OutputImageGif container
        """
        if display_manager.preferences.graphics == 'disable':
            return
        types = display_manager.types
        preferred = (
            ('PNG',  types.OutputImagePng),
            ('JPEG', types.OutputImageJpg),
            ('GIF',  types.OutputImageGif),
        )
        from sage.repl.rich_output.buffer import OutputBuffer
        for format, output_container in preferred:
            if output_container in display_manager.supported_output():
                stream = io.BytesIO()
                try:
                    self.pil.save(stream, format=format)
                except IOError:
                    # not all formats support all modes, e.g. no alpha support in gif
                    continue
                buf = OutputBuffer(stream.getvalue())
                return output_container(buf)
