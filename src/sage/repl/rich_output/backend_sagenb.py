# -*- encoding: utf-8 -*-
"""
SageNB Backend for the Sage Rich Output System

This module defines the IPython backends for
:mod:`sage.repl.rich_output`.

EXAMPLES:

Install the SageNB displayhook while doctesting for the rest of this
file. This somewhat odd incantation is how SageNB installed its
displayhook before::

    sage: from sage.misc.displayhook import DisplayHook
    sage: import sys
    sage: sys.displayhook = DisplayHook()

    sage: from sage.repl.rich_output import get_display_manager
    sage: get_display_manager()
    The Sage display manager using the SageNB backend

We also enable the SageNB magic global variable::

    sage: import sage.plot.plot
    sage: sage.plot.plot.EMBEDDED_MODE = True

And switch to a temporary directory so our current directory will not
get cluttered with temporary files::

    sage: os.chdir(tmp_dir())

The SageNB notebook is based on saving data files with predictable
filenames::

    sage: os.path.exists('sage0.png')
    False
    sage: Graphics()
    sage: os.path.exists('sage0.png')
    True
    sage: os.remove('sage0.png')

Tables are typeset as html in SageNB::

    sage: table([1, 2, 3])
    <html><div class="notruncate">
    <table  class="table_form">
    <tbody>
    <tr class ="row-a">
    <td><script type="math/tex">1</script></td>
    <td><script type="math/tex">2</script></td>
    <td><script type="math/tex">3</script></td>
    </tr>
    </tbody>
    </table>
    </div></html>
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
import stat
from sage.misc.cachefunc import cached_method
from sage.misc.html import html
from sage.misc.temporary_file import graphics_filename
from sage.doctest import DOCTEST_MODE
from sage.repl.rich_output.backend_base import BackendBase
from sage.repl.rich_output.output_catalog import *
from sage.repl.rich_output.output_video import OutputVideoBase


def world_readable(filename):
    """
    All SageNB temporary files must be world-writeable.
    
    Discussion of this design choice can be found at :trac:`17743`.

    EXAMPLES::

        sage: import os, stat
        sage: f = tmp_filename()
    
    At least on a sane system the temporary files are only readable by
    the user, but not by others in the group or total strangers::

        sage: mode = os.stat(f).st_mode
        sage: bool(mode & stat.S_IRUSR), bool(mode & stat.S_IRGRP), bool(mode & stat.S_IROTH)   # random output
        (True, False, False)

    This function disables that protection::

        sage: from sage.repl.rich_output.backend_sagenb import world_readable
        sage: world_readable(f)
        sage: mode = os.stat(f).st_mode
        sage: bool(mode & stat.S_IRUSR), bool(mode & stat.S_IRGRP), bool(mode & stat.S_IROTH)
        (True, True, True)
    """
    os.chmod(filename, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)


class SageNbOutputSceneJmol(OutputSceneJmol):
    """
    Adapt Jmol rich output container for SageNB.

    For legacy reasons, SageNB expects Jmol files saved under strange
    names. This class takes care of that.

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
        sage: SageNbOutputSceneJmol.example()
        SageNbOutputSceneJmol container
    """
    
    @cached_method
    def sagenb_launch_script_filename(self):
        """
        Return the launch script filename used by SageNB

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
            sage: j = SageNbOutputSceneJmol.example()
            sage: j.sagenb_launch_script_filename()
            'sage0-size32.jmol'
        """
        import PIL.Image
        from StringIO import StringIO
        width, height = PIL.Image.open(StringIO(self.preview_png.get())).size
        ext = '-size{0}.jmol'.format(width)
        return graphics_filename(ext=ext)

    @cached_method
    def _base_filename(self):
        """
        Return the common part of the file name

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
            sage: j = SageNbOutputSceneJmol.example()
            sage: j._base_filename()
            'sage0-size32'
            sage: j.sagenb_launch_script_filename().startswith(j._base_filename())
            True
            sage: j.scene_zip_filename().startswith(j._base_filename())
            True
        """
        dot_jmol = '.jmol'
        filename = self.sagenb_launch_script_filename()
        assert filename.endswith(dot_jmol)
        return filename[:-len(dot_jmol)]
    
    @cached_method
    def scene_zip_filename(self):
        """
        Return the filename for the scene zip archive

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
            sage: j = SageNbOutputSceneJmol.example()
            sage: j.scene_zip_filename()
            'sage0-size32-....jmol.zip'
        """
        from random import randint
        return '{0}-{1}.jmol.zip'.format(
            self._base_filename(),
            randint(0, 1 << 30)
        )

    @cached_method
    def preview_filename(self):
        """
        Return the filename for the png preview

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
            sage: j = SageNbOutputSceneJmol.example()
            sage: j.preview_filename()
            './.jmol_images/sage0-size32.jmol.png'
        """
        directory, filename = os.path.split(self._base_filename())
        if not directory:
            directory = '.'
        return '{0}/.jmol_images/{1}.jmol.png'.format(directory, filename)

    def save_launch_script(self):
        """
        Save the Jmol launch script
        
        See :meth:`sagenb_launch_script_filename`.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
            sage: j = SageNbOutputSceneJmol.example()
            sage: os.path.exists('sage0-size32.jmol')
            False
            sage: j.save_launch_script()
            sage: os.path.exists('sage0-size32.jmol')
            True
        """
        from sagenb.notebook.interact import SAGE_CELL_ID
        path = 'cells/{0}/{1}'.format(
            SAGE_CELL_ID,
            self.scene_zip_filename())
        with open(self.sagenb_launch_script_filename(), 'w') as f:
            f.write('set defaultdirectory "{0}"\n'.format(path))
            f.write('script SCRIPT\n')
        world_readable(self.sagenb_launch_script_filename())

    def save_preview(self):
        """
        Save the preview PNG image
        
        See :meth:`preview_filename`.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
            sage: j = SageNbOutputSceneJmol.example()
            sage: import shutil
            sage: shutil.rmtree('.jmol_images', ignore_errors=True)
            sage: j.save_preview()
            sage: os.listdir('.jmol_images')
            ['sage1-size32.jmol.png']
        """
        from sage.misc.misc import sage_makedirs
        sage_makedirs('.jmol_images')
        self.preview_png.save_as(self.preview_filename())
        world_readable(self.preview_filename())
    
    def embed(self):
        """
        Save all files necessary to embed jmol
        
        EXAMPLES:

        Switch to a new empty temporary directory::

            sage: os.chdir(tmp_dir())

            sage: from sage.repl.rich_output.backend_sagenb import SageNbOutputSceneJmol
            sage: j = SageNbOutputSceneJmol.example()
            sage: j.embed()
            sage: sorted(os.listdir('.'))
            ['.jmol_images', 'sage0-size32-....jmol.zip', 'sage0-size32.jmol']
            sage: sorted(os.listdir('.jmol_images'))
            ['sage0-size32.jmol.png']
        """
        self.save_preview()
        self.save_launch_script()
        self.scene_zip.save_as(self.scene_zip_filename())
        world_readable(self.scene_zip_filename())


class BackendSageNB(BackendBase):

    def _repr_(self):
        """
        Return the string representation
        
        OUTPUT:

        String

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import BackendSageNB
            sage: backend = BackendSageNB()
            sage: backend._repr_()
            'SageNB'
        """
        return 'SageNB'
    
    def supported_output(self):
        """
        Return the outputs that are supported by the SageNB backend.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).
        The order is ignored.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_sagenb import BackendSageNB
            sage: backend = BackendSageNB()
            sage: supp = backend.supported_output();  supp     # random output
            set([<class 'sage.repl.rich_output.output_graphics.OutputPlainText'>, 
                 ...,
                 <class 'sage.repl.rich_output.output_graphics.OutputCanvas3d'>])
            sage: from sage.repl.rich_output.output_basic import OutputLatex
            sage: OutputLatex in supp
            True
        """
        return set([
            OutputPlainText, OutputAsciiArt, OutputLatex,
            OutputHtml,
            OutputImagePng, OutputImageGif, OutputImageJpg,
            OutputImagePdf, OutputImageSvg,
            SageNbOutputSceneJmol,
            OutputSceneCanvas3d,
            OutputVideoOgg, OutputVideoWebM, OutputVideoMp4,
        ])

    def display_immediately(self, plain_text, rich_output):
        r"""
        Show output without waiting for the prompt.

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

        This method does not return anything.

        EXAMPLES::

            sage: import sage.repl.rich_output.output_catalog as catalog
            sage: plain_text = catalog.OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_sagenb import BackendSageNB
            sage: backend = BackendSageNB()
            sage: backend.display_immediately(plain_text, plain_text)
            Example plain text output
            sage: latex = catalog.OutputLatex.example()
            sage: backend.display_immediately(plain_text, latex)
            <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\int \sin\left(x\right)\,{d x}</script></html>
        """
        if isinstance(rich_output, (OutputPlainText, OutputAsciiArt)):
            rich_output.print_to_stdout()
        elif isinstance(rich_output, OutputLatex):
            print(rich_output.mathjax())
        elif isinstance(rich_output, OutputHtml):
            print(rich_output.with_html_tag())
        elif isinstance(rich_output, OutputImagePng):
            self.embed_image(rich_output.png, '.png')
        elif isinstance(rich_output, OutputImageGif):
            self.embed_image(rich_output.gif, '.gif')
        elif isinstance(rich_output, OutputImageJpg):
            self.embed_image(rich_output.jpg, '.jpg')
        elif isinstance(rich_output, OutputImagePdf):
            self.embed_image(rich_output.pdf, '.pdf')
        elif isinstance(rich_output, OutputImageSvg):
            self.embed_image(rich_output.svg, '.svg')
        elif isinstance(rich_output, OutputSceneJmol):
            rich_output.embed()
        elif isinstance(rich_output, OutputSceneCanvas3d):
            self.embed_image(rich_output.canvas3d, '.canvas3d')
        elif isinstance(rich_output, OutputVideoBase):
            self.embed_video(rich_output)
        else:
            raise TypeError('rich_output type not supported, got {0}'.format(rich_output))

    def embed_image(self, output_buffer, file_ext):
        """
        Embed Image in the SageNB worksheet

        SageNB scans per-cell directories for image files, so all we
        have to do here is to save the image at the right place.

        INPUT:

        - ``output_buffer`` --
          :class:`~sage.repl.rich_output.buffer.Buffer`. A buffer
          holding the image data.
        
        - ``file_ext`` -- string. The file extension to use for saving
          the image.

        OUTPUT:

        Nothing is returned. The image file is saved in the
        appropriate place for SageNB.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: rich_output = dm.types.OutputImagePng.example()
            sage: os.path.exists('sage0.png')
            False
            sage: dm._backend.embed_image(rich_output.png, '.png')
            sage: os.path.exists('sage0.png')
            True
        """
        filename = graphics_filename(ext=file_ext)
        output_buffer.save_as(filename)
        world_readable(filename)

    def embed_video(self, video_output):
        filename = graphics_filename(ext=video_output.ext)
        video_output.video.save_as(filename)
        world_readable(filename)
        html(video_output.html_fragment(
            url='cell://' + filename,
            link_attrs='class="file_link"',
        ))
