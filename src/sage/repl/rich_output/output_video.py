# sage_setup: distribution = sagemath-repl
r"""
Video Output Types

This module defines the rich output types for video formats.
"""
# ****************************************************************************
#       Copyright (C) 2015 Martin von Gagern <Martin.vGagern@gmx.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import importlib.resources

from sage.repl.rich_output.output_basic import OutputBase
from sage.repl.rich_output.buffer import OutputBuffer


class OutputVideoBase(OutputBase):

    def __init__(self, video, loop=True):
        """
        Abstract base class for rich video output.

        INPUT:

        - ``video`` -- :class:`~sage.repl.rich_output.buffer.OutputBuffer`;
          the video data
        - ``loop`` -- boolean; whether to repeat the video in an endless loop

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputVideoOgg
            sage: OutputVideoOgg.example()  # indirect doctest
            OutputVideoOgg container
        """
        assert isinstance(video, OutputBuffer)
        self.video = video
        self.loop = loop

    @classmethod
    def example(cls):
        r"""
        Construct a sample video output container.

        This static method is meant for doctests, so they can easily
        construct an example.  The method is implemented in the abstract
        :class:`OutputVideoBase` class, but should get invoked on a
        concrete subclass for which an actual example can exist.

        OUTPUT: an instance of the class on which this method is called

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputVideoOgg
            sage: OutputVideoOgg.example()
            OutputVideoOgg container
            sage: OutputVideoOgg.example().video
            buffer containing 5612 bytes
            sage: OutputVideoOgg.example().ext
            '.ogv'
            sage: OutputVideoOgg.example().mimetype
            'video/ogg'
        """
        with importlib.resources.path(__package__, 'example' + cls.ext) as filename:
            return cls(OutputBuffer.from_file(filename),
                       {'controls': True, 'loop': False})

    def html_fragment(self, url, link_attrs=''):
        r"""
        Construct a HTML fragment for embedding this video.

        INPUT:

        - ``url`` -- string; the URL where the data of this video can be found

        - ``link_attrs`` -- string; can be used to style the fallback link
          which is presented to the user if the video is not supported

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputVideoOgg
            sage: print(OutputVideoOgg.example().html_fragment
            ....:       ('foo', 'class="bar"').replace('><','>\n<'))
            <video autoplay="autoplay" controls="controls" loop="loop">
            <source src="foo" type="video/ogg" />
            <p>
            <a target="_new" href="foo" class="bar">Download video/ogg video</a>
            </p>
            </video>
        """
        attrs = {
            'autoplay': 'autoplay',
            'controls': 'controls',
        }
        if self.loop:
            attrs['loop'] = 'loop'
        attrs = ''.join(' {}="{}"'.format(k, v)
                        for k, v in sorted(attrs.items()))
        txt = ('<video{attrs}>'
               '<source src="{url}" type="{mimetype}" /><p>'
               '<a target="_new" href="{url}" {link_attrs}>'
               'Download {mimetype} video</a></p></video>')
        return txt.format(url=url,
                          mimetype=self.mimetype,
                          attrs=attrs,
                          link_attrs=link_attrs)


class OutputVideoOgg(OutputVideoBase):
    """
    Ogg video, Ogg Theora in particular.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoOgg
        sage: OutputVideoOgg.example()
        OutputVideoOgg container
    """

    ext = ".ogv"
    mimetype = "video/ogg"


class OutputVideoWebM(OutputVideoBase):
    """
    WebM video.

    The video can be encoded using VP8, VP9 or an even more recent codec.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoWebM
        sage: OutputVideoWebM.example()
        OutputVideoWebM container
    """

    ext = ".webm"
    mimetype = "video/webm"


class OutputVideoMp4(OutputVideoBase):
    """
    MPEG 4 video.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoMp4
        sage: OutputVideoMp4.example()
        OutputVideoMp4 container
    """

    ext = ".mp4"
    mimetype = "video/mp4"


class OutputVideoFlash(OutputVideoBase):
    """
    Flash video.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoFlash
        sage: OutputVideoFlash.example()
        OutputVideoFlash container
    """

    ext = ".flv"
    mimetype = "video/x-flv"


class OutputVideoMatroska(OutputVideoBase):
    """
    Matroska Video.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoMatroska
        sage: OutputVideoMatroska.example()
        OutputVideoMatroska container
    """

    ext = ".mkv"
    mimetype = "video/x-matroska"


class OutputVideoAvi(OutputVideoBase):
    """
    AVI video.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoAvi
        sage: OutputVideoAvi.example()
        OutputVideoAvi container
    """

    ext = ".avi"
    mimetype = "video/x-msvideo"


class OutputVideoWmv(OutputVideoBase):
    """
    Windows Media Video.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoWmv
        sage: OutputVideoWmv.example()
        OutputVideoWmv container
    """

    ext = ".wmv"
    mimetype = "video/x-ms-wmv"


class OutputVideoQuicktime(OutputVideoBase):
    """
    Quicktime video.

    EXAMPLES::

        sage: from sage.repl.rich_output.output_catalog import OutputVideoQuicktime
        sage: OutputVideoQuicktime.example()
        OutputVideoQuicktime container
    """

    ext = ".mov"
    mimetype = "video/quicktime"
