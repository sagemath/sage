# -*- coding: utf-8 -*-
r"""
Features for testing the presence of ``cython``
"""

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
