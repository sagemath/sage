# sage_setup: distribution = sagemath-environment
r"""
Check for symengine_py
"""

# ****************************************************************************
#       Copyright (C) 2023 Dima Pasechnik
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import PythonModule
from .join_feature import JoinFeature


class symengine_py(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of the
    Python package :ref:`symengine_py <spkg_symengine_py>`.

    EXAMPLES::

        sage: from sage.features.symengine_py import symengine_py
        sage: symengine_py().is_present()                    # optional - symengine_py
        FeatureTestResult('symengine_py', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.symengine_py import symengine_py
            sage: isinstance(symengine_py(), symengine_py)
            True
        """
        JoinFeature.__init__(self, 'symengine_py',
                             [PythonModule('symengine', spkg='symengine_py',
                                            url='https://pypi.org/project/symengine')])

def all_features():
    return [symengine_py()]
