# -*- encoding: utf-8 -*-
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

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sys
from sage.repl.rich_output.backend_base import BackendBase
from sage.repl.rich_output.output_catalog import *


    
class BackendDoctest(BackendBase):

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_doctest import BackendDoctest
            sage: backend = BackendDoctest()
            sage: backend._repr_()
            'doctest'
        """
        return 'doctest'

    def default_preferences(self):
        """
        Return the backend's display preferences

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
            * graphics is not specified
            * supplemental_plot = never
            * text is not specified
        """
        from sage.repl.rich_output.preferences import DisplayPreferences
        return DisplayPreferences(supplemental_plot='never')

    def install(self, **kwds):
        """
        Switch to the the doctest backend

        This method is being called from within
        :meth:`~sage.repl.rich_output.display_manager.DisplayManager.switch_backend`. You
        should never call it by hand.

        INPUT:

        None of the optional keyword arguments are used in the doctest
        backend.

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
        Switch away from the doctest backend

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
        Return the supported output types
        
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
            OutputPlainText, OutputAsciiArt, OutputUnicodeArt, OutputLatex,
            OutputImagePng, OutputImageGif, OutputImageJpg, 
            OutputImageSvg, OutputImagePdf, OutputImageDvi,
            OutputSceneJmol, OutputSceneCanvas3d, OutputSceneWavefront,
        ])

    def displayhook(self, plain_text, rich_output):
        """
        Display object from displayhook

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

            sage: plt = plot(sin)
            sage: plt
            Graphics object consisting of 1 graphics primitive
            sage: plt.show()

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.displayhook(plt)       # indirect doctest
            Graphics object consisting of 1 graphics primitive
        """
        self.validate(rich_output)
        if any(isinstance(rich_output, cls)
               for cls in [OutputPlainText, OutputAsciiArt, OutputLatex]):
            rich_output.print_to_stdout()
        else:
            plain_text.print_to_stdout()

    def display_immediately(self, plain_text, rich_output):
        """
        Display object immediately

        INPUT:

        Same as :meth:`displayhook`.

        EXAMPLES:

        The following example does not call the displayhook. More
        precisely, the ``show()`` method returns ``None`` which is
        ignored by the displayhook. When running the example on a Sage
        display backend capable of displaying graphics outside of the
        displayhook, the plot is still shown. Nothing is shown during
        doctests::

            sage: plt = plot(sin)
            sage: plt
            Graphics object consisting of 1 graphics primitive
            sage: plt.show()

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.display_immediately(plt)   # indirect doctest
        """
        self.validate(rich_output)
        types_to_print = [OutputPlainText, OutputAsciiArt, OutputUnicodeArt]
        if isinstance(rich_output, OutputLatex):
            print(rich_output.mathjax(display=False))
        elif any(isinstance(rich_output, cls) for cls in types_to_print):
            rich_output.print_to_stdout()

    def validate(self, rich_output):
        """
        Perform checks on ``rich_output``
        
        INPUT:

        - ``rich_output`` -- instance of a subclass of
          :class:`~sage.repl.rich_output.output_basic.OutputBase`.

        OUTPUT:

        An assertion is triggered if ``rich_output`` is invalid.

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
        """
        if isinstance(rich_output, OutputPlainText):
            pass
        elif isinstance(rich_output, OutputAsciiArt):
            pass
        elif isinstance(rich_output, OutputUnicodeArt):
            pass
        elif isinstance(rich_output, OutputLatex):
            assert rich_output.mathjax().startswith('<html>')
        elif isinstance(rich_output, OutputImagePng):
            assert rich_output.png.get().startswith('\x89PNG')
        elif isinstance(rich_output, OutputImageGif):
            assert rich_output.gif.get().startswith('GIF89a')
        elif isinstance(rich_output, OutputImageJpg):
            assert rich_output.jpg.get().startswith('\xff\xd8\xff\xe0\x00\x10JFIF')
        elif isinstance(rich_output, OutputImageSvg):
            assert '</svg>' in rich_output.svg.get()
        elif isinstance(rich_output, OutputImagePdf):
            assert rich_output.pdf.get().startswith('%PDF-')
        elif isinstance(rich_output, OutputImageDvi):
            assert 'TeX output' in rich_output.dvi.get()
        elif isinstance(rich_output, OutputSceneJmol):
            assert rich_output.preview_png.get().startswith('\x89PNG')
            assert rich_output.scene_zip.get().startswith('PK')  # zip archive
        elif isinstance(rich_output, OutputSceneWavefront):
            assert rich_output.obj.get().startswith('mtllib ')
            assert rich_output.mtl.get().startswith('newmtl ')
        elif isinstance(rich_output, OutputSceneCanvas3d):
            assert rich_output.canvas3d.get().startswith('[{vertices:')
        else:
            raise TypeError('rich_output type not supported')
