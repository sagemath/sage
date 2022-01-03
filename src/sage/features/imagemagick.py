# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of ``imagemagick``

Currently we only check for the presence of ``convert``. When needed other
commands like ``magick``, ``magick-script``, ``convert``, ``mogrify``,
``identify``, ``composite``, ``montage``, ``compare``, etc. could be also
checked in this module.
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
from .join_feature import JoinFeature


class ImageMagick(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    ``ImageMagick``

    Currently, only the availability of ``convert`` is checked.

    EXAMPLES::

        sage: from sage.features.imagemagick import ImageMagick
        sage: ImageMagick().is_present()  # optional - imagemagick
        FeatureTestResult('imagemagick', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.imagemagick import ImageMagick
            sage: isinstance(ImageMagick(), ImageMagick)
            True
        """
        JoinFeature.__init__(self, "imagemagick",
                             [Executable("convert", executable="convert")],
                             spkg="imagemagick",
                             url="https://www.imagemagick.org/")


def all_features():
    return [ImageMagick()]
