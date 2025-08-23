# sage_setup: distribution = sagemath-environment
r"""
Determination of programs for viewing web pages, etc.

The function :func:`default_viewer` defines reasonable defaults for
these programs.  To use something else, use ``viewer``.  First
import it::

    sage: from sage.misc.viewer import viewer

On OS X, PDFs are opened by default using the 'open' command, which
runs whatever has been designated as the PDF viewer in the OS.  To
change this to use 'Adobe Reader'::

    sage: viewer.pdf_viewer('open -a /Applications/Adobe\ Reader.app') # not tested

Similarly, you can set ``viewer.browser(...)``, ``viewer.dvi_viewer(...)``,
and ``viewer.png_viewer(...)``.  You can make this change permanent by adding
lines like these to your :file:`SAGE_STARTUP_FILE` (which is
:file:`$HOME/.sage/init.sage` by default)::

    from sage.misc.viewer import viewer
    viewer.pdf_viewer('open -a /Applications/Adobe\ Reader.app')

Functions and classes
---------------------
"""

import platform
from sage.structure.sage_object import SageObject


VIEWERS = ['browser', 'dvi_viewer', 'pdf_viewer', 'png_viewer']


def default_viewer(viewer=None):
    """
    Set up default programs for opening web pages, PDFs, PNGs, and DVI files.

    INPUT:

    - ``viewer`` -- ``None`` or a string; one of ``'browser'``, ``'pdf'``,
      ``'png'``, ``'dvi'``. Return the name of the corresponding program.
      ``None`` is treated the same as ``'browser'``.

    EXAMPLES::

        sage: from sage.misc.viewer import default_viewer
        sage: default_viewer(None) # random -- depends on OS, etc.
        'open'
        sage: default_viewer('pdf') # random -- depends on OS, etc.
        'xdg-open'
        sage: default_viewer('jpg')
        Traceback (most recent call last):
        ...
        ValueError: Unknown type of viewer: jpg.
    """
    import os
    from sage.misc.sage_ostools import have_program

    if isinstance(viewer, str):
        viewer = viewer.lower()

    if 'BROWSER' in os.environ:
        BROWSER = os.environ['BROWSER']
        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER

    elif platform.system() == 'Darwin':
        # Simple on OS X, since there is an open command that opens
        # anything, using the user's preferences.
        BROWSER = 'open -W'
        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER

    elif have_program('xdg-open'):
        # On other OS'es try xdg-open if present.
        # See http://portland.freedesktop.org/xdg-utils-1.0.
        BROWSER = 'xdg-open'
        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER

    else:
        # If all fails try to get something from the environment.
        try:
            BROWSER = os.environ['BROWSER']
        except KeyError:
            BROWSER = 'less'  # silly default; lets hope it doesn't come to this!
            for cmd in ['firefox', 'konqueror', 'mozilla', 'mozilla-firefox']:
                if have_program(cmd):
                    BROWSER = cmd
                    break
        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER

        # Alternatives, if they are set in the environment or available.
        try:
            DVI_VIEWER = os.environ['DVI_VIEWER']
        except KeyError:
            for cmd in ['xdvi', 'kdvi']:
                if have_program(cmd):
                    DVI_VIEWER = cmd
                    break
        try:
            PDF_VIEWER = os.environ['PDF_VIEWER']
        except KeyError:
            for cmd in ['acroread', 'xpdf']:
                if have_program(cmd):
                    PDF_VIEWER = cmd
                    break

    if viewer is None or viewer.startswith('browse'):
        return BROWSER
    elif viewer.startswith('dvi'):
        return DVI_VIEWER
    elif viewer.startswith('png'):
        return PNG_VIEWER
    elif viewer.startswith('pdf'):
        return PDF_VIEWER
    else:
        raise ValueError('Unknown type of viewer: {}.'.format(viewer))


# _viewer_prefs: a dictionary holding global preferences for viewers.
_viewer_prefs = {}


class Viewer(SageObject):
    """
    Set defaults for various viewing applications: a web browser, a
    dvi viewer, a pdf viewer, and a png viewer.

    EXAMPLES::

        sage: from sage.misc.viewer import viewer
        sage: old_browser = viewer.browser()  # indirect doctest
        sage: viewer.browser('open -a /Applications/Firefox.app')
        sage: viewer.browser()
        'open -a /Applications/Firefox.app'
        sage: viewer.browser(old_browser) # restore old value
    """
    def _set(self, app=None, TYPE='browser'):
        r"""
        Change the default viewer. Return the current setting if the
        argument ``app`` is ``None``.

        INPUT:

        - ``app`` -- ``None`` or a string; the program to use
        - ``TYPE`` -- string (default: ``'browser'``); must be in the list
          ``VIEWERS`` defined in :mod:`sage.misc.viewer`

        EXAMPLES::

            sage: from sage.misc.viewer import viewer
            sage: old_browser = viewer.browser()
            sage: viewer.browser('open -a /Applications/Firefox.app') # indirect doctest
            sage: viewer.browser()
            'open -a /Applications/Firefox.app'
            sage: viewer.browser(old_browser) # restore old value
        """
        TYPE = TYPE.lower()
        if TYPE not in VIEWERS:
            raise ValueError('Unrecognized type of viewer: {}'.format(TYPE))
        if app is None:
            try:
                return _viewer_prefs[TYPE]
            except KeyError:
                return default_viewer(TYPE)
        elif app == 'default':
            try:
                del _viewer_prefs[TYPE]
            except KeyError:
                pass
        else:
            _viewer_prefs[TYPE] = app

    def browser(self, app=None):
        r"""
        Change the default browser. Return the current setting if arg
        is ``None``, which is the default.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use

        EXAMPLES::

            sage: from sage.misc.viewer import viewer
            sage: old_browser = viewer.browser()
            sage: viewer.browser('open -a /Applications/Firefox.app') # indirect doctest
            sage: viewer.browser()
            'open -a /Applications/Firefox.app'
            sage: viewer.browser(old_browser) # restore old value
        """
        return self._set(app, TYPE='browser')

    def dvi_viewer(self, app=None):
        r"""
        Change the default dvi viewer. Return the current setting if arg
        is ``None``, which is the default.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use

        EXAMPLES::

            sage: from sage.misc.viewer import viewer
            sage: old_dvi_app = viewer.dvi_viewer()
            sage: viewer.dvi_viewer('/usr/bin/xdvi') # indirect doctest
            sage: viewer.dvi_viewer()
            '/usr/bin/xdvi'
            sage: viewer.dvi_viewer(old_dvi_app) # restore old value
        """
        return self._set(app, TYPE='dvi_viewer')

    def pdf_viewer(self, app=None):
        r"""
        Change the default pdf viewer. Return the current setting if arg
        is ``None``, which is the default.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use

        EXAMPLES::

            sage: from sage.misc.viewer import viewer
            sage: old_pdf_app = viewer.pdf_viewer()
            sage: viewer.pdf_viewer('/usr/bin/pdfopen') # indirect doctest
            sage: viewer.pdf_viewer()
            '/usr/bin/pdfopen'
            sage: viewer.pdf_viewer(old_pdf_app) # restore old value
        """
        return self._set(app, TYPE='pdf_viewer')

    def png_viewer(self, app=None):
        r"""
        Change the default png viewer. Return the current setting if arg
        is ``None``, which is the default.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use

        EXAMPLES::

            sage: from sage.misc.viewer import viewer
            sage: old_png_app = viewer.png_viewer()
            sage: viewer.png_viewer('display') # indirect doctest
            sage: viewer.png_viewer()
            'display'
            sage: viewer.png_viewer(old_png_app) # restore old value
        """
        return self._set(app, TYPE='png_viewer')

    def __call__(self, x=None):
        """
        Return the current setting for a viewer program.

        INPUT:

        - ``x`` -- string

        If ``x`` is ``None`` or starts with 'browse', then return the
        current browser app.  If ``x`` starts with 'dvi', return the
        current dvi viewer, and similarly if ``x`` starts with 'png'
        or 'pdf'.

        EXAMPLES::

            sage: from sage.misc.viewer import viewer
            sage: viewer('pdf') # random -- depends on OS, etc.
            'mozilla'
            sage: viewer('browser') == viewer()
            True
        """
        if isinstance(x, str):
            x = x.lower()

        if x is None or x.startswith('browse'):
            return self.browser()
        elif x.startswith('dvi'):
            return self.dvi_viewer()
        elif x.startswith('png'):
            return self.png_viewer()
        elif x.startswith('pdf'):
            return self.pdf_viewer()


viewer = Viewer()


def browser():
    """
    Return the program used to open a web page.  By default, the
    program used depends on the platform and other factors, like
    settings of certain environment variables.  To use a different
    program, call ``viewer.browser('PROG')``, where 'PROG' is the
    desired program.

    EXAMPLES::

        sage: from sage.misc.viewer import browser
        sage: browser() # random -- depends on OS, etc.
        'open'
    """
    return viewer.browser()


def dvi_viewer():
    """
    Return the program used to display a dvi file.  By default, the
    program used depends on the platform and other factors, like
    settings of certain environment variables.  To use a different
    program, call ``viewer.dvi_viewer('PROG')``, where 'PROG' is the
    desired program.

    EXAMPLES::

        sage: from sage.misc.viewer import dvi_viewer
        sage: dvi_viewer() # random -- depends on OS, etc.
        'open'
    """
    viewer()
    return viewer.dvi_viewer()


def pdf_viewer():
    """
    Return the program used to display a pdf file.  By default, the
    program used depends on the platform and other factors, like
    settings of certain environment variables.  To use a different
    program, call ``viewer.pdf_viewer('PROG')``, where 'PROG' is the
    desired program.

    EXAMPLES::

        sage: from sage.misc.viewer import pdf_viewer, viewer
        sage: old_pdf_app = viewer.pdf_viewer()
        sage: viewer.pdf_viewer('acroread')
        sage: pdf_viewer()
        'acroread'
        sage: viewer.pdf_viewer('old_pdf_app')
    """
    viewer()
    return viewer.pdf_viewer()


def png_viewer():
    """
    Return the program used to display a png file. By default, the
    program used depends on the platform and other factors, like
    settings of certain environment variables.  To use a different
    program, call ``viewer.png_viewer('PROG')``, where 'PROG' is the
    desired program.

    EXAMPLES::

        sage: from sage.misc.viewer import png_viewer
        sage: png_viewer() # random -- depends on OS, etc.
        'xdg-open'
    """
    viewer()
    return viewer.png_viewer()
