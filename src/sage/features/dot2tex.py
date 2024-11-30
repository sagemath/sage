# sage_setup: distribution = sagemath-environment
r"""
Check for ``dot2tex``
"""

# *****************************************************************************
#       Copyright (C) 2024 Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import PythonModule


class dot2tex(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of :ref:`dot2tex <spkg_dot2tex>`.

    dot2tex is provided by an optional package in the Sage distribution.

    EXAMPLES::

        sage: from sage.features.dot2tex import dot2tex
        sage: dot2tex().is_present()                     # optional - dot2tex
        FeatureTestResult('dot2tex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.dot2tex import dot2tex
            sage: isinstance(dot2tex(), dot2tex)
            True
        """
        PythonModule.__init__(self, 'dot2tex', spkg='dot2tex')


def all_features():
    return [dot2tex()]
