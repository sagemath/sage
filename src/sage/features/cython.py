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

from sage.features import Feature

class sage__misc__cython(Feature):
    r"""
    A :class:`~sage.features.Feature` that is always present,
    indicating that :mod:`sage.misc.cython` is always available.

    This class can be removed once all ``needs sage.misc.cython``
    tags have been eliminated from the doctests.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.cython import sage__misc__cython
            sage: sage__misc__cython()
            Feature('sage.misc.cython')
        """
        Feature.__init__(self, "sage.misc.cython")

    def _is_present(self):
        r"""
        TESTS::

            sage: from sage.features.cython import sage__misc__cython
            sage: sage__misc__cython()._is_present()
            FeatureTestResult('sage.misc.cython', True)

        """
        from sage.features import FeatureTestResult
        return FeatureTestResult(self, True)


def all_features():
    return [sage__misc__cython()]
