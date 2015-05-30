# -*- encoding: utf-8 -*-
r"""
Graphics Output Types

This module defines the rich output types for 2-d images, both vector
and raster graphics.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os

from sage.repl.rich_output.output_basic import OutputBase
from sage.repl.rich_output.buffer import OutputBuffer


class OutputImagePng(OutputBase):

    def __init__(self, png):
        """
        PNG Image

        .. NOTE::

            Every backend that is capable of displaying any kind of
            graphics is supposed to support the PNG format at least.

        INPUT:

        - ``png`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. The
          PNG image data.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: OutputImagePng.example()  # indirect doctest
            OutputImagePng container
        """
        self.png = OutputBuffer(png)

    @classmethod
    def example(cls):
        r"""
        Construct a sample PNG output container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputImagePng`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: OutputImagePng.example()
            OutputImagePng container
            sage: OutputImagePng.example().png
            buffer containing 608 bytes
            sage: OutputImagePng.example().png.get().startswith('\x89PNG')
            True
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.png')
        with open(filename) as f:
            return cls(f.read())
        

class OutputImageGif(OutputBase):

    def __init__(self, gif):
        """
        GIF Image (possibly animated)

        INPUT:

        - ``gif`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. The
          GIF image data.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageGif
            sage: OutputImageGif.example()   # indirect doctest
            OutputImageGif container
        """
        self.gif = OutputBuffer(gif)

    @classmethod
    def example(cls):
        r"""
        Construct a sample GIF output container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputImageGif`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageGif
            sage: OutputImageGif.example()
            OutputImageGif container
            sage: OutputImageGif.example().gif
            buffer containing 408 bytes
            sage: OutputImageGif.example().gif.get().startswith('GIF89a')
            True
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.gif')
        with open(filename) as f:
            return cls(f.read())


class OutputImageJpg(OutputBase):

    def __init__(self, jpg):
        """
        JPEG Image

        INPUT:

        - ``jpeg`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. The
          JPEG image data.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageJpg
            sage: OutputImageJpg.example()   # indirect doctest
            OutputImageJpg container
        """
        self.jpg = OutputBuffer(jpg)

    @classmethod
    def example(cls):
        r"""
        Construct a sample JPEG output container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputImageJpg`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageJpg
            sage: OutputImageJpg.example()
            OutputImageJpg container
            sage: OutputImageJpg.example().jpg
            buffer containing 978 bytes
            sage: OutputImageJpg.example().jpg.get().startswith('\xff\xd8\xff\xe0\x00\x10JFIF')
            True
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.jpg')
        with open(filename) as f:
            return cls(f.read())

        
class OutputImageSvg(OutputBase):

    def __init__(self, svg):
        """
        SVG Image

        INPUT:

        - ``SVG`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. The
          SVG image data.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageSvg
            sage: OutputImageSvg.example()   # indirect doctest
            OutputImageSvg container
        """
        self.svg = OutputBuffer(svg)

    @classmethod
    def example(cls):
        r"""
        Construct a sample SVG output container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputImageSvg`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageSvg
            sage: OutputImageSvg.example()
            OutputImageSvg container
            sage: OutputImageSvg.example().svg
            buffer containing 1422 bytes
            sage: '</svg>' in OutputImageSvg.example().svg.get()
            True
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.svg')
        with open(filename) as f:
            return cls(f.read())


class OutputImagePdf(OutputBase):

    def __init__(self, pdf):
        """
        PDF Image

        INPUT:

        - ``pdf`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. The
          PDF data.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImagePdf
            sage: OutputImagePdf.example()   # indirect doctest
            OutputImagePdf container
        """
        self.pdf = OutputBuffer(pdf)

    @classmethod
    def example(cls):
        r"""
        Construct a sample PDF output container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputImagePdf`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImagePdf
            sage: OutputImagePdf.example()
            OutputImagePdf container
            sage: OutputImagePdf.example().pdf
            buffer containing 4285 bytes
            sage: OutputImagePdf.example().pdf.get().startswith('%PDF-1.4')
            True
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.pdf')
        with open(filename) as f:
            return cls(f.read())


class OutputImageDvi(OutputBase):

    def __init__(self, dvi):
        """
        DVI Image

        INPUT:

        - ``dvi`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. The
          DVI data.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageDvi
            sage: OutputImageDvi.example()     # indirect doctest
            OutputImageDvi container
        """
        self.dvi = OutputBuffer(dvi)

    @classmethod
    def example(cls):
        r"""
        Construct a sample DVI output container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputImageDvi`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImageDvi
            sage: OutputImageDvi.example()
            OutputImageDvi container
            sage: OutputImageDvi.example().dvi
            buffer containing 212 bytes
            sage: 'TeX output' in OutputImageDvi.example().dvi.get()
            True
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.dvi')
        with open(filename) as f:
            return cls(f.read())

