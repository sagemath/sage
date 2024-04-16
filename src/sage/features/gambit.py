# sage_setup: distribution = sagemath-environment
r"""
Check for pygambit
"""

# ****************************************************************************
#       Copyright (C) 2024 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import PythonModule


class pygambit(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of the
    Python package :ref:`pygambit <spkg_pygambit>`.

    EXAMPLES::

        sage: from sage.features.gambit import pygambit
        sage: pygambit().is_present()                           # optional - pygambit
        FeatureTestResult('pygambit', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.gambit import pygambit
            sage: isinstance(pygambit(), pygambit)
            True
        """
        PythonModule.__init__(self, 'pygambit', spkg="pygambit")


def all_features():
    return [pygambit()]
