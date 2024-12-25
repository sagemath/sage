# sage_setup: distribution = sagemath-repl
"""
The backend used for doctests

This backend is active during doctests. It should mimic the behavior
of the IPython command line as close as possible. Without actually
launching image viewers, of course.

EXAMPLES::

    sage: from sage.repl.rich_output import get_display_manager
    sage: get_display_manager()
    The Sage display manager using the doctest backend
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys

from sage.repl.rich_output.backend_base import BackendBase
from sage.repl.rich_output.output_catalog import *


class BackendDoctest(BackendBase):

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_doctest import BackendDoctest
            sage: backend = BackendDoctest()
            sage: backend._repr_()
            'doctest'
        """
        return 'doctest'

    def default_preferences(self):
        """
        Return the backend's display preferences.

        Matches the IPython command line display preferences to keep
        the differences between that and the doctests to a minimum.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.preferences.DisplayPreferences`.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.default_preferences()
            Display preferences:
            * align_latex is not specified
            * graphics is not specified
            * supplemental_plot = never
            * text is not specified
        """
        from sage.repl.rich_output.preferences import DisplayPreferences
        return DisplayPreferences(supplemental_plot='never')

    def install(self, **kwds):
        """
        Switch to the doctest backend.

        This method is being called from within
        :meth:`~sage.repl.rich_output.display_manager.DisplayManager.switch_backend`. You
        should never call it by hand.

        INPUT:

        - ``**kwds`` -- none of the optional keyword arguments are used in the
          doctest backend

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_doctest import BackendDoctest
            sage: backend = BackendDoctest()
            sage: backend.install()
            sage: backend.uninstall()
        """
        self._old_displayhook = sys.displayhook
        sys.displayhook = self.get_display_manager().displayhook

    def uninstall(self):
        """
        Switch away from the doctest backend.

        This method is being called from within
        :meth:`~sage.repl.rich_output.display_manager.DisplayManager.switch_backend`. You
        should never call it by hand.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_doctest import BackendDoctest
            sage: backend = BackendDoctest()
            sage: backend.install()
            sage: backend.uninstall()
        """
        sys.displayhook = self._old_displayhook

    def supported_output(self):
        """
        Return the supported output types.

        OUTPUT:

        Set of subclasses of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`, the
        supported output container types.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_doctest import BackendDoctest
            sage: from sage.repl.rich_output.output_catalog import *
            sage: backend = BackendDoctest()
            sage: OutputPlainText in backend.supported_output()
            True
            sage: OutputSceneJmol in backend.supported_output()
            True
        """
        return set([
            OutputPlainText, OutputAsciiArt, OutputUnicodeArt,
            OutputImagePng, OutputImageGif, OutputImageJpg,
            OutputImageSvg, OutputImagePdf, OutputImageDvi,
            OutputSceneJmol, OutputSceneCanvas3d, OutputSceneWavefront,
            OutputVideoOgg, OutputVideoWebM, OutputVideoMp4,
            OutputVideoFlash, OutputVideoMatroska, OutputVideoAvi,
            OutputVideoWmv, OutputVideoQuicktime,
        ])

    def displayhook(self, plain_text, rich_output):
        """
        Display object from displayhook.

        INPUT:

        - ``plain_text`` -- instance of
          :class:`~sage.repl.rich_output.output_basic.OutputPlainText`. The
          plain text version of the output.

        - ``rich_output`` -- instance of an output container class
          (subclass of
          :class:`~sage.repl.rich_output.output_basic.OutputBase`). Guaranteed
          to be one of the output containers returned from
          :meth:`supported_output`, possibly the same as
          ``plain_text``.

        EXAMPLES:

        This ends up calling the displayhook::

            sage: plt = plot(sin)                                                       # needs sage.plot sage.symbolic
            sage: plt                                                                   # needs sage.plot sage.symbolic
            Graphics object consisting of 1 graphics primitive
            sage: plt.show()                                                            # needs sage.plot sage.symbolic

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.displayhook(plt)       # indirect doctest                          # needs sage.plot sage.symbolic
            Graphics object consisting of 1 graphics primitive
        """
        self.validate(rich_output)
        if any(isinstance(rich_output, cls)
               for cls in [OutputPlainText, OutputAsciiArt, OutputLatex, OutputHtml]):
            rich_output.print_to_stdout()
        else:
            plain_text.print_to_stdout()

    def display_immediately(self, plain_text, rich_output):
        """
        Display object immediately.

        INPUT:

        Same as :meth:`displayhook`.

        EXAMPLES:

        The following example does not call the displayhook. More
        precisely, the ``show()`` method returns ``None`` which is
        ignored by the displayhook. When running the example on a Sage
        display backend capable of displaying graphics outside of the
        displayhook, the plot is still shown. Nothing is shown during
        doctests::

            sage: plt = plot(sin)                                                       # needs sage.plot sage.symbolic
            sage: plt                                                                   # needs sage.plot sage.symbolic
            Graphics object consisting of 1 graphics primitive
            sage: plt.show()                                                            # needs sage.plot sage.symbolic

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.display_immediately(plt)   # indirect doctest                      # needs sage.plot sage.symbolic
        """
        self.validate(rich_output)
        types_to_print = [OutputPlainText, OutputAsciiArt, OutputUnicodeArt, OutputHtml]
        if any(isinstance(rich_output, cls) for cls in types_to_print):
            rich_output.print_to_stdout()

    def validate(self, rich_output):
        """
        Perform checks on ``rich_output``.

        INPUT:

        - ``rich_output`` -- instance of a subclass of
          :class:`~sage.repl.rich_output.output_basic.OutputBase`

        OUTPUT: an assertion is triggered if ``rich_output`` is invalid

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: invalid = dm.types.OutputImagePng('invalid')
            sage: backend = dm._backend;  backend
            doctest
            sage: backend.validate(invalid)
            Traceback (most recent call last):
            ...
            AssertionError
            sage: backend.validate(dm.types.OutputPlainText.example())
            sage: backend.validate(dm.types.OutputAsciiArt.example())
            sage: backend.validate(dm.types.OutputLatex.example())
            sage: backend.validate(dm.types.OutputImagePng.example())
            sage: backend.validate(dm.types.OutputImageGif.example())
            sage: backend.validate(dm.types.OutputImageJpg.example())
            sage: backend.validate(dm.types.OutputImageSvg.example())
            sage: backend.validate(dm.types.OutputImagePdf.example())
            sage: backend.validate(dm.types.OutputImageDvi.example())
            sage: backend.validate(dm.types.OutputSceneJmol.example())
            sage: backend.validate(dm.types.OutputSceneWavefront.example())
            sage: backend.validate(dm.types.OutputSceneCanvas3d.example())
            sage: backend.validate(dm.types.OutputVideoOgg.example())
            sage: backend.validate(dm.types.OutputVideoWebM.example())
            sage: backend.validate(dm.types.OutputVideoMp4.example())
            sage: backend.validate(dm.types.OutputVideoFlash.example())
            sage: backend.validate(dm.types.OutputVideoMatroska.example())
            sage: backend.validate(dm.types.OutputVideoAvi.example())
            sage: backend.validate(dm.types.OutputVideoWmv.example())
            sage: backend.validate(dm.types.OutputVideoQuicktime.example())
        """
        if isinstance(rich_output, OutputPlainText):
            pass
        elif isinstance(rich_output, OutputAsciiArt):
            pass
        elif isinstance(rich_output, OutputUnicodeArt):
            pass
        elif isinstance(rich_output, OutputLatex):
            pass
        elif isinstance(rich_output, OutputHtml):
            pass
        elif isinstance(rich_output, OutputImagePng):
            assert rich_output.png.get().startswith(b'\x89PNG')
        elif isinstance(rich_output, OutputImageGif):
            assert rich_output.gif.get().startswith(b'GIF89a')
        elif isinstance(rich_output, OutputImageJpg):
            assert rich_output.jpg.get().startswith(b'\xff\xd8\xff\xe0\x00\x10JFIF')
        elif isinstance(rich_output, OutputImageSvg):
            assert b'</svg>' in rich_output.svg.get()
        elif isinstance(rich_output, OutputImagePdf):
            assert rich_output.pdf.get().startswith(b'%PDF-')
        elif isinstance(rich_output, OutputImageDvi):
            assert b'TeX output' in rich_output.dvi.get()
        elif isinstance(rich_output, OutputSceneJmol):
            assert rich_output.preview_png.get().startswith(b'\x89PNG')
            assert rich_output.scene_zip.get().startswith(b'PK')  # zip archive
        elif isinstance(rich_output, OutputSceneWavefront):
            assert rich_output.obj.get().startswith(b'mtllib ')
            assert rich_output.mtl.get().startswith(b'newmtl ')
        elif isinstance(rich_output, OutputSceneCanvas3d):
            assert rich_output.canvas3d.get().startswith(b'[{"vertices":')
        elif isinstance(rich_output, OutputVideoOgg):
            assert rich_output.video.get().startswith(b'OggS')
        elif isinstance(rich_output, OutputVideoWebM):
            data = rich_output.video.get()
            assert data.startswith(b'\x1a\x45\xdf\xa3')
            assert b'\x42\x82\x84webm' in data
        elif isinstance(rich_output, OutputVideoMp4):
            data = rich_output.video.get()
            assert data[4:8] == b'ftyp'
            assert data.startswith(b'\0\0\0')
            # See http://www.ftyps.com/
            ftyps = [data[i:i+4] for i in range(8, data[3], 4)]
            del ftyps[1]  # version number, not an ftyp
            expected = [b'avc1', b'iso2', b'mp41', b'mp42']
            assert any(i in ftyps for i in expected)
        elif isinstance(rich_output, OutputVideoFlash):
            assert rich_output.video.get().startswith(b'FLV\x01')
        elif isinstance(rich_output, OutputVideoMatroska):
            data = rich_output.video.get()
            assert data.startswith(b'\x1a\x45\xdf\xa3')
            assert b'\x42\x82\x88matroska' in data
        elif isinstance(rich_output, OutputVideoAvi):
            data = rich_output.video.get()
            assert data[:4] == b'RIFF' and data[8:12] == b'AVI '
        elif isinstance(rich_output, OutputVideoWmv):
            assert rich_output.video.get().startswith(b'\x30\x26\xb2\x75')
        elif isinstance(rich_output, OutputVideoQuicktime):
            data = rich_output.video.get()
            assert data[4:12] == b'ftypqt  ' or data[4:8] == b'moov'
        else:
            raise TypeError('rich_output type not supported')
