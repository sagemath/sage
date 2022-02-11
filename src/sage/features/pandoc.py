# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of ``pandoc``
"""
# ****************************************************************************
#       Copyright (C) 2018 Thierry Monteil <sage!lma.metelu.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable


class Pandoc(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``pandoc``

    EXAMPLES::

        sage: from sage.features.pandoc import Pandoc
        sage: Pandoc().is_present()  # optional - pandoc
        FeatureTestResult('pandoc', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.pandoc import Pandoc
            sage: isinstance(Pandoc(), Pandoc)
            True
        """
        Executable.__init__(self, "pandoc", executable="pandoc",
                            url="https://pandoc.org/")


def all_features():
    return [Pandoc()]
