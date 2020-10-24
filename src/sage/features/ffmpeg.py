# -*- coding: utf-8 -*-
r"""
Check for FFmpeg
"""
# ****************************************************************************
#       Copyright (C) 2018 Sebastien Labbe <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable


class FFmpeg(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``FFmpeg``

    EXAMPLES::

        sage: from sage.features.ffmpeg import FFmpeg
        sage: FFmpeg().is_present()  # optional: ffmpeg
        FeatureTestResult('FFmpeg', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.ffmpeg import FFmpeg
            sage: isinstance(FFmpeg(), FFmpeg)
            True
        """
        Executable.__init__(self, "FFmpeg", executable="ffmpeg",
                            url="https://www.ffmpeg.org/")
