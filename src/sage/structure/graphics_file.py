# -*- coding: utf-8 -*-
"""
Wrapper for Graphics Files
"""

import os

from sage.misc.temporary_file import tmp_filename
from sage.structure.sage_object import SageObject
import sage.doctest


class Mime(object):
    TEXT = 'text/plain'
    HTML = 'text/html'
    LATEX = 'text/latex'
    JSON = 'application/json'
    JAVASCRIPT = 'application/javascript'
    PDF = 'application/pdf'
    PNG = 'image/png'
    JPG = 'image/jpeg'
    SVG = 'image/svg+xml'

    JMOL = 'application/jmol'

    @classmethod
    def validate(cls, value):
        """
        Check that input is known mime type

        INPUT:

        - ``value`` -- string.

        OUTPUT:

        Unicode string of that mime type. A ``ValueError`` is raised
        if input is incorrect / unknown.

        EXAMPLES::

            sage: from sage.structure.graphics_file import Mime
            sage: Mime.validate('image/png')
            'image/png'
            sage: Mime.validate('foo/bar')
            Traceback (most recent call last):
            ...
            ValueError: unknown mime type
        """
        value = str(value).lower()
        for k, v in cls.__dict__.items():
            if isinstance(v, str) and v == value:
                return v
        raise ValueError('unknown mime type')

    @classmethod
    def extension(cls, mime_type):
        """
        Return file extension.

        INPUT:

        - ``mime_type`` -- mime type as string.

        OUTPUT:

        String containing the usual file extension for that type of
        file. Excludes ``os.extsep``.

        EXAMPLES::

            sage: from sage.structure.graphics_file import Mime
            sage: Mime.extension('image/png')
            'png'
        """
        try:
            return preferred_filename_ext[mime_type]
        except KeyError:
            raise ValueError('no known extension for mime type')


preferred_filename_ext = {
    Mime.TEXT: 'txt',
    Mime.HTML: 'html',
    Mime.LATEX: 'tex',
    Mime.JSON: 'json',
    Mime.JAVASCRIPT: 'js',
    Mime.PDF: 'pdf',
    Mime.PNG: 'png',
    Mime.JPG: 'jpg',
    Mime.SVG: 'svg',
    Mime.JMOL: 'spt.zip',
}


mimetype_for_ext = dict(
    (value, key) for (key, value) in preferred_filename_ext.items()
)


class GraphicsFile(SageObject):

    def __init__(self, filename, mime_type=None):
        """
        Wrapper around a graphics file.
        """
        self._filename = filename
        if mime_type is None:
            mime_type = self._guess_mime_type(filename)
        self._mime = Mime.validate(mime_type)

    def _guess_mime_type(self, filename):
        """
        Guess mime type from file extension
        """
        ext = os.path.splitext(filename)[1]
        ext = ext.lstrip(os.path.extsep)
        try:
            return mimetype_for_ext[ext]
        except KeyError:
            raise ValueError('unknown file extension, please specify mime type')

    def _repr_(self):
        """
        Return a string representation.
        """
        return 'Graphics file {0}'.format(self.mime())

    def filename(self):
        return self._filename

    def save_as(self, filename):
        """
        Make the file available under a new filename.

        INPUT:

        - ``filename`` -- string. The new filename.

        The newly-created ``filename`` will be a hardlink if
        possible. If not, an independent copy is created.
        """
        try:
            os.link(self.filename(), filename)
        except OSError:
            import shutil
            shutil.copy2(self.filename(), filename)

    def mime(self):
        return self._mime

    def data(self):
        """
        Return a byte string containing the image file.
        """
        with open(self._filename, 'rb') as f:
            return f.read()

    def launch_viewer(self):
        """
        Launch external viewer for the graphics file.

        .. note::

            Does not actually launch a new process when doctesting.

        EXAMPLES::

            sage: from sage.structure.graphics_file import GraphicsFile
            sage: g = GraphicsFile('/tmp/test.png', 'image/png')
            sage: g.launch_viewer()
        """
        if sage.doctest.DOCTEST_MODE:
            return
        if self.mime() == Mime.JMOL:
            return self._launch_jmol()
        from sage.misc.viewer import viewer
        command = viewer(preferred_filename_ext[self.mime()])
        os.system('{0} {1} 2>/dev/null 1>/dev/null &'
                  .format(command, self.filename()))
        # TODO: keep track of opened processes...

    def _launch_jmol(self):
        launch_script = tmp_filename(ext='.spt')
        with open(launch_script, 'w') as f:
            f.write('set defaultdirectory "{0}"\n'.format(self.filename()))
            f.write('script SCRIPT\n')
        os.system('jmol {0} 2>/dev/null 1>/dev/null &'
                  .format(launch_script))


def graphics_from_save(save_function, preferred_mime_types,
                       allowed_mime_types=None, figsize=None, dpi=None):
    """
    Helper function to construct a graphics file.

    INPUT:

    - ``save_function`` -- callable that can save graphics to a file
      and accepts options like
      :meth:`sage.plot.graphics.Graphics.save`.

    - ``preferred_mime_types`` -- list of mime types. The graphics
      output mime types in order of preference (i.e. best quality to
      worst).

    - ``allowed_mime_types`` -- set of mime types (as strings). The
      graphics types that we can display. Output, if any, will be one
      of those.

    - ``figsize`` -- pair of integers (optional). The desired graphics
      size in pixels. Suggested, but need not be respected by the
      output.

    - ``dpi`` -- integer (optional). The desired resolution in dots
      per inch. Suggested, but need not be respected by the output.

    OUTPUT:

    Return an instance of
    :class:`sage.structure.graphics_file.GraphicsFile` encapsulating a
    suitable image file. Image is one of the
    ``preferred_mime_types``. If ``allowed_mime_types`` is specified,
    the resulting file format matches one of these.

    Alternatively, this function can return ``None`` to indicate that
    textual representation is preferable and/or no graphics with the
    desired mime type can be generated.
    """
    # Figure out best mime type
    mime = None
    if allowed_mime_types is None:
        mime = Mime.PNG
    else:
        # order of preference
        for m in preferred_mime_types:
            if m in allowed_mime_types:
                mime = m
                break
    if mime is None:
        return None    # don't know how to generate suitable graphics
    # Generate suitable temp file
    filename = tmp_filename(ext=os.path.extsep + Mime.extension(mime))
    # Call the save_function with the right arguments
    kwds = {}
    if figsize is not None:
        kwds['figsize'] = figsize
    if dpi is not None:
        kwds['dpi'] = dpi
    save_function(filename, **kwds)
    return GraphicsFile(filename, mime)
