r"""
Check for Regina
"""
from . import PythonModule


class Regina(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of Regina.

    Regina is provided by an optional package in the Sage distribution.

    EXAMPLES::

        sage: from sage.features.regina import Regina
        sage: Regina().is_present()                     # optional - regina
        FeatureTestResult('regina', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.regina import Regina
            sage: isinstance(Regina(), Regina)
            True
        """
        PythonModule.__init__(self, 'regina', spkg='regina')


def all_features():
    return [Regina()]
