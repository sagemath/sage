# sage_setup: distribution = sagemath-repl
r"""
Emacs sage-mode Backend for the Sage Rich Output System

This module defines the Emacs backend for :mod:`sage.repl.rich_output`
based on the IPython shell version.
"""

# ****************************************************************************
#       Copyright (C) 2015 Ivan Andrus <darthandrus@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
from sage.repl.rich_output.output_catalog import *


class BackendEmacs(BackendIPythonCommandline):
    """
    Emacs Backend.

    This backend is used by Emacs' sage-mode to have typeset output
    and inline images.

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_emacs import BackendEmacs
        sage: BackendEmacs()
        Emacs sage-mode
    """

    def _repr_(self):
        r"""
        Return string representation of the backend.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_emacs import BackendEmacs
            sage: backend = BackendEmacs()
            sage: backend._repr_()
            'Emacs sage-mode'
        """
        return 'Emacs sage-mode'

    def default_preferences(self):
        """
        Return the backend's display preferences.

        Override this method to change the default preferences when
        using your backend.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.preferences.DisplayPreferences`.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_emacs import BackendEmacs
            sage: backend = BackendEmacs()
            sage: backend.default_preferences()
            Display preferences:
            * align_latex is not specified
            * graphics is not specified
            * supplemental_plot is not specified
            * text is not specified
        """
        from sage.repl.rich_output.preferences import DisplayPreferences
        return DisplayPreferences()

    def displayhook(self, plain_text, rich_output):
        r"""
        Backend implementation of the displayhook.

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

        OUTPUT:

        Because this is based on the IPython commandline display hook,
        it returns the IPython display data, a pair of
        dictionaries. The first dictionary contains mime types as keys
        and the respective output as value. The second dictionary is
        metadata.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText, OutputLatex
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_emacs import BackendEmacs
            sage: backend = BackendEmacs()
            sage: backend.displayhook(plain_text, plain_text)
            ({'text/plain': 'Example plain text output'}, {})
            sage: latex_text = OutputLatex.example()
            sage: backend.displayhook(plain_text, latex_text)
            ({'text/plain': 'BEGIN_TEXT:Example plain text output:END_TEXT\nBEGIN_LATEX:\\newcommand{\\Bold}[1]{\\mathbf{#1}}\\int \\sin\\left(x\\right)\\,{d x}:END_LATEX'},
              {})
        """

        if isinstance(rich_output, OutputPlainText):
            return ({'text/plain': rich_output.text.get_str()}, {})
        elif isinstance(rich_output, OutputAsciiArt):
            return ({'text/plain': rich_output.ascii_art.get_str()}, {})
        elif isinstance(rich_output, OutputLatex):
            text = "BEGIN_TEXT:" + plain_text.text.get_str() + ":END_TEXT\nBEGIN_LATEX:" + \
                   rich_output.latex.get_str() + ":END_LATEX"
            return ({'text/plain': text}, {})

        # TODO: perhaps handle these by returning the data inline,
        # e.g. base64 encoded, so that sage-mode can show inline
        # images for remotely running shells.
        elif isinstance(rich_output, OutputImagePng):
            msg = self.launch_viewer(
                rich_output.png.filename(ext='png'), plain_text.text.get())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImageGif):
            msg = self.launch_viewer(
                rich_output.gif.filename(ext='gif'), plain_text.text.get())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImagePdf):
            msg = self.launch_viewer(
                rich_output.pdf.filename(ext='pdf'), plain_text.text.get())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImageDvi):
            msg = self.launch_viewer(
                rich_output.dvi.filename(ext='dvi'), plain_text.text.get())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputSceneJmol):
            msg = self.launch_jmol(rich_output, plain_text.text.get())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputSceneWavefront):
            msg = self.launch_sage3d(rich_output, plain_text.text.get())
            return ({'text/plain': msg}, {})
        else:
            raise TypeError('rich_output type not supported')
