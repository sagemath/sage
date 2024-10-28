# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``pynormaliz``
"""

# *****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import PythonModule
from .join_feature import JoinFeature


class PyNormaliz(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the
    Python package :ref:`PyNormaliz <spkg_pynormaliz>`.

    EXAMPLES::

        sage: from sage.features.normaliz import PyNormaliz
        sage: PyNormaliz().is_present()                    # optional - pynormaliz
        FeatureTestResult('pynormaliz', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.normaliz import PyNormaliz
            sage: isinstance(PyNormaliz(), PyNormaliz)
            True
        """
        JoinFeature.__init__(self, 'pynormaliz',
                             [PythonModule('PyNormaliz', spkg='pynormaliz')])


def all_features():
    return [PyNormaliz()]
