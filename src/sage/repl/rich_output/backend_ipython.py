# sage_setup: distribution = sagemath-repl
"""
IPython Backend for the Sage Rich Output System

This module defines the IPython backends for
:mod:`sage.repl.rich_output`.
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import html
from IPython.display import publish_display_data
from sage.repl.rich_output.backend_base import BackendBase
from sage.repl.rich_output.output_catalog import *


class BackendIPython(BackendBase):
    """
    Common base for the IPython UIs.

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_ipython import BackendIPython
        sage: BackendIPython()._repr_()
        Traceback (most recent call last):
        ...
        NotImplementedError: derived classes must implement this method
    """

    def install(self, **kwds):
        """
        Switch the Sage rich output to the IPython backend.

        INPUT:

        - ``shell`` -- keyword argument; the IPython shell

        No tests since switching away from the doctest rich output
        backend will break the doctests.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.repl.rich_output.backend_ipython import BackendIPython
            sage: backend = BackendIPython()
            sage: shell = get_test_shell()
            sage: backend.install(shell=shell)
            sage: shell.run_cell('1+1')
            2
        """
        shell = kwds['shell']
        from sage.repl.display.formatter import SageDisplayFormatter
        shell.display_formatter = SageDisplayFormatter(parent=shell)
        shell.configurables.append(shell.display_formatter)

    def set_underscore_variable(self, obj):
        """
        Set the ``_`` builtin variable.

        Since IPython handles the history itself, this does nothing.

        INPUT:

        - ``obj`` -- anything

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.repl.rich_output.backend_ipython import BackendIPython
            sage: backend = BackendIPython()
            sage: backend.set_underscore_variable(123)
            sage: _
            0
        """
        pass

    def display_immediately(self, plain_text, rich_output):
        """
        Show output immediately.

        This method is similar to the rich output :meth:`displayhook`,
        except that it can be invoked at any time.

        INPUT:

        Same as :meth:`displayhook`.

        OUTPUT: this method does not return anything

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: _ = backend.display_immediately(plain_text, plain_text)
            Example plain text output
        """
        formatted, metadata = self.displayhook(plain_text, rich_output)
        if not formatted:
            return
        publish_display_data(data=formatted, metadata=metadata)


class BackendIPythonCommandline(BackendIPython):
    """
    Backend for the IPython Command Line.

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
        sage: BackendIPythonCommandline()
        IPython command line
    """

    def default_preferences(self):
        """
        Return the backend's display preferences.

        The default for the commandline is to not plot graphs since
        the launching of an external viewer is considered too
        disruptive.

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

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend._repr_()
            'IPython command line'
        """
        return 'IPython command line'

    def supported_output(self):
        """
        Return the outputs that are supported by the IPython commandline backend.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).
        The order is ignored.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: supp = backend.supported_output();  supp     # random output
            set([<class 'sage.repl.rich_output.output_graphics.OutputImageGif'>,
                 ...,
                 <class 'sage.repl.rich_output.output_graphics.OutputImagePng'>])
            sage: from sage.repl.rich_output.output_basic import OutputLatex
            sage: OutputLatex in supp
            True
        """
        return set([
            OutputPlainText, OutputAsciiArt, OutputUnicodeArt, OutputLatex,
            OutputImagePng, OutputImageGif,
            OutputImagePdf, OutputImageDvi,
            OutputSceneJmol, OutputSceneWavefront, OutputSceneThreejs,
        ])

    def displayhook(self, plain_text, rich_output):
        """
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

        The IPython commandline display hook returns the IPython
        display data, a pair of dictionaries. The first dictionary
        contains mime types as keys and the respective output as
        value. The second dictionary is metadata.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.displayhook(plain_text, plain_text)
            ({'text/plain': 'Example plain text output'}, {})

        TESTS:

        We verify that unicode strings work::

            sage: class Foo(sage.structure.sage_object.SageObject):
            ....:     def _rich_repr_(self, dm):
            ....:         return dm.types.OutputPlainText('Motörhead')
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.displayhook(Foo())
            ({'text/plain': 'Mot\xf6rhead'}, {})
        """
        if isinstance(rich_output, OutputPlainText):
            return ({'text/plain': rich_output.text.get_str()}, {})
        elif isinstance(rich_output, OutputAsciiArt):
            return ({'text/plain': rich_output.ascii_art.get_str()}, {})
        elif isinstance(rich_output, OutputUnicodeArt):
            return ({'text/plain': rich_output.unicode_art.get_str()}, {})
        elif isinstance(rich_output, OutputLatex):
            return ({'text/plain': rich_output.latex.get_str()}, {})
        elif isinstance(rich_output, OutputImagePng):
            msg = self.launch_viewer(
                rich_output.png.filename(ext='png'), plain_text.text.get_str())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImageGif):
            msg = self.launch_viewer(
                rich_output.gif.filename(ext='gif'), plain_text.text.get_str())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImagePdf):
            msg = self.launch_viewer(
                rich_output.pdf.filename(ext='pdf'), plain_text.text.get_str())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImageDvi):
            msg = self.launch_viewer(
                rich_output.dvi.filename(ext='dvi'), plain_text.text.get_str())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputSceneJmol):
            msg = self.launch_jmol(rich_output, plain_text.text.get_str())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputSceneWavefront):
            msg = self.launch_sage3d(rich_output, plain_text.text.get_str())
            return ({'text/plain': msg}, {})
        elif isinstance(rich_output, OutputSceneThreejs):
            msg = self.launch_viewer(
                rich_output.html.filename(ext='html'), plain_text.text.get_str())
            return ({'text/plain': msg}, {})
        else:
            raise TypeError('rich_output type not supported')

    def display_immediately(self, plain_text, rich_output):
        """
        Show output without going back to the command line prompt.

        This method is similar to the rich output :meth:`displayhook`,
        except that it can be invoked at any time. On the Sage command
        line it launches viewers just like :meth:`displayhook`.

        INPUT:

        Same as :meth:`displayhook`.

        OUTPUT: this method does not return anything

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.display_immediately(plain_text, plain_text)
            Example plain text output
        """
        formatdata, metadata = self.displayhook(plain_text, rich_output)
        print(formatdata['text/plain'])

    def launch_viewer(self, image_file, plain_text):
        """
        Launch external viewer for the graphics file.

        INPUT:

        - ``image_file`` -- string; file name of the image file

        - ``plain_text`` -- string; the plain text representation of
          the image file

        OUTPUT:

        String. Human-readable message indicating whether the viewer
        was launched successfully.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.launch_viewer('/path/to/foo.bar', 'Graphics object')
            'Launched bar viewer for Graphics object'
        """
        base, dot_ext = os.path.splitext(image_file)
        ext = dot_ext.lstrip(os.path.extsep)
        from sage.misc.viewer import viewer
        command = viewer(ext)
        if not command:
            command = viewer.browser()
        from sage.doctest import DOCTEST_MODE
        if not DOCTEST_MODE:
            os.system('{0} {1} 2>/dev/null 1>/dev/null &'
                      .format(command, image_file))
        return 'Launched {0} viewer for {1}'.format(ext, plain_text)

    def launch_jmol(self, output_jmol, plain_text):
        """
        Launch the stand-alone jmol viewer.

        INPUT:

        - ``output_jmol`` --
          :class:`~sage.repl.rich_output.output_graphics3d.OutputSceneJmol`; the
          scene to launch Jmol with

        - ``plain_text`` -- string; the plain text representation

        OUTPUT:

        string; human-readable message indicating that the viewer was launched.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: backend.launch_jmol(OutputSceneJmol.example(), 'Graphics3d object')   # needs sage.plot
            'Launched jmol viewer for Graphics3d object'
        """
        from sage.doctest import DOCTEST_MODE
        from sage.interfaces.jmoldata import JmolData
        jdata = JmolData()
        if not jdata.is_jmol_available() and not DOCTEST_MODE:
            raise RuntimeError('jmol cannot run, no suitable java version found')
        launch_script = output_jmol.launch_script_filename()
        jmol_cmd = 'jmol'
        if not DOCTEST_MODE:
            os.system('{0} {1} 2>/dev/null 1>/dev/null &'
                      .format(jmol_cmd, launch_script))
        return 'Launched jmol viewer for {0}'.format(plain_text)

    def is_in_terminal(self):
        """
        Test whether the UI is meant to run in a terminal.

        See
        :meth:`sage.repl.rich_output.display_manager.DisplayManager.is_in_terminal`
        for details.

        OUTPUT: ``True`` for the IPython commandline

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.is_in_terminal()
            True
        """
        return True

    def threejs_offline_scripts(self):
        """
        Three.js script for the IPython command line.

        OUTPUT: string containing script tag

        EXAMPLES::

            sage: # needs threejs
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.threejs_offline_scripts()
            '...<script ...</script>...'
        """
        from sage.features.threejs import Threejs

        if not Threejs().is_present():
            return ''

        script = Threejs().absolute_filename()

        return '\n<script src="{0}"></script>'.format(script)


IFRAME_TEMPLATE = """
<iframe srcdoc="{escaped_html}"
        width="{width}"
        height="{height}"
        style="border: 0;">
</iframe>
"""


class BackendIPythonNotebook(BackendIPython):
    """
    Backend for the IPython Notebook.

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
        sage: BackendIPythonNotebook()
        IPython notebook
    """

    def _repr_(self):
        """
        Return string representation of the backend.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: backend._repr_()
            'IPython notebook'
        """
        return 'IPython notebook'

    def supported_output(self):
        """
        Return the outputs that are supported by the IPython notebook backend.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).
        The order is ignored.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: supp = backend.supported_output();  supp     # random output
            set([<class 'sage.repl.rich_output.output_graphics.OutputPlainText'>,
                 ...,
                 <class 'sage.repl.rich_output.output_graphics.OutputImagePdf'>])
            sage: from sage.repl.rich_output.output_basic import OutputLatex
            sage: OutputLatex in supp
            True
            sage: from sage.repl.rich_output.output_graphics import OutputImageGif
            sage: OutputImageGif in supp
            True
        """
        return set([
            OutputPlainText, OutputAsciiArt, OutputUnicodeArt, OutputLatex,
            OutputHtml,
            OutputImagePng, OutputImageGif, OutputImageJpg,
            OutputImageSvg, OutputImagePdf,
            OutputSceneJmol, OutputSceneThreejs,
        ])

    def displayhook(self, plain_text, rich_output):
        """
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

        The IPython notebook display hook returns the IPython
        display data, a pair of dictionaries. The first dictionary
        contains mime types as keys and the respective output as
        value. The second dictionary is metadata.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: backend.displayhook(plain_text, plain_text)
            ({'text/plain': 'Example plain text output'}, {})
        """
        if isinstance(rich_output, OutputPlainText):
            return ({'text/plain': rich_output.text.get_str()}, {})
        elif isinstance(rich_output, OutputAsciiArt):
            return ({'text/plain': rich_output.ascii_art.get_str()}, {})
        elif isinstance(rich_output, OutputUnicodeArt):
            return ({'text/plain': rich_output.unicode_art.get_str()}, {})
        elif isinstance(rich_output, OutputLatex):
            return ({'text/latex': rich_output.latex.get_str(),
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        elif isinstance(rich_output, OutputHtml):
            data = {'text/html': rich_output.html.get_str(),
                    'text/plain': plain_text.text.get_str()}
            if rich_output.latex:
                data['text/latex'] = rich_output.latex.get_str()
            return (data, {})
        elif isinstance(rich_output, OutputImagePng):
            return ({'image/png': rich_output.png.get(),
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        elif isinstance(rich_output, OutputImageGif):
            return ({'text/html': rich_output.html_fragment(),
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        elif isinstance(rich_output, OutputImageJpg):
            return ({'image/jpeg': rich_output.jpg.get(),
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        elif isinstance(rich_output, OutputImageSvg):
            return ({'image/svg+xml': rich_output.svg.get(),
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        elif isinstance(rich_output, OutputImagePdf):
            return ({'image/png': rich_output.png.get(),
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        elif isinstance(rich_output, OutputSceneJmol):
            from sage.repl.display.jsmol_iframe import JSMolHtml
            jsmol = JSMolHtml(rich_output, height=500)
            return ({'text/html': jsmol.iframe(),
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        elif isinstance(rich_output, OutputSceneThreejs):
            escaped_html = html.escape(rich_output.html.get_str())
            iframe = IFRAME_TEMPLATE.format(
                escaped_html=escaped_html,
                width='100%',
                height=400,
            )
            return ({'text/html': iframe,
                     'text/plain': plain_text.text.get_str(),
                     }, {})
        else:
            raise TypeError('rich_output type not supported')

    def threejs_offline_scripts(self):
        """
        Three.js script for the IPython notebook.

        OUTPUT: string containing script tag

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: backend.threejs_offline_scripts()                                     # needs sage.plot
            '...<script src="/nbextensions/threejs-sage/r.../three.min.js...<\\/script>...'
        """
        from sage.repl.rich_output import get_display_manager
        from sage.features.threejs import Threejs
        CDN_script = get_display_manager().threejs_scripts(online=True)
        CDN_script = CDN_script.replace('</script>', r'<\/script>').replace('\n', ' \\\n')
        return """
<script src="/nbextensions/threejs-sage/{}/three.min.js"></script>
<script>
  if ( !window.THREE ) document.write('{}');
</script>
        """.format(Threejs().required_version(), CDN_script)
