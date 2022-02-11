r"""
Features for testing the presence of ``tdlib``
"""

from . import PythonModule
from .join_feature import JoinFeature


class Tdlib(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the ``tdlib``.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.tdlib import Tdlib
            sage: isinstance(Tdlib(), Tdlib)
            True
        """
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_tdlib' later
        JoinFeature.__init__(self, 'tdlib',
                             [PythonModule('sage.graphs.graph_decompositions.tdlib', spkg='tdlib')])


def all_features():
    return [Tdlib()]
