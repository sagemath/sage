r"""
Check for Khoca
"""
from . import PythonModule


class Khoca(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of Khoca.

    Khoca is provided by an optional package in the Sage distribution.

    EXAMPLES::

        sage: from sage.features.khoca import Khoca
        sage: Khoca().is_present()                     # optional - khoca
        FeatureTestResult('khoca', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.khoca import Khoca
            sage: isinstance(Khoca(), Khoca)
            True
        """
        PythonModule.__init__(self, 'khoca', spkg='khoca')


def all_features():
    return [Khoca()]
