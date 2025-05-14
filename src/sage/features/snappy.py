r"""
Check for SnapPy
"""
from . import PythonModule


class SnapPy(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of SnapPy.

    SnapPy is provided by an optional package in the Sage distribution.

    EXAMPLES::

        sage: from sage.features.snappy import SnapPy
        sage: SnapPy().is_present()                     # optional - snappy
        FeatureTestResult('snappy', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.snappy import SnapPy
            sage: isinstance(SnapPy(), SnapPy)
            True
        """
        PythonModule.__init__(self, 'snappy', spkg='snappy')


def all_features():
    return [SnapPy()]
