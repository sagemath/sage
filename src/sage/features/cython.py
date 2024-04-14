# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``cython``
"""

# *****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import CythonFeature


class sage__misc__cython(CythonFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether :mod:`sage.misc.cython`
    is available and functional.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features import CythonFeature
            sage: from sage.features.cython import sage__misc__cython
            sage: isinstance(sage__misc__cython(), CythonFeature)
            True
        """
        # It suffices to use a trivial CythonFeature because CythonFeature
        # is implemented via sage.misc.cython.
        CythonFeature.__init__(self, "sage.misc.cython", test_code="")


def all_features():
    return [sage__misc__cython()]
