r"""
Check for phitigra
"""
from . import PythonModule


class Phitigra(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of phitigra.

    Phitigra is provided by an optional package in the Sage distribution.

    EXAMPLES::

        sage: from sage.features.phitigra import Phitigra
        sage: Phitigra().is_present()                     # optional - phitigra
        FeatureTestResult('phitigra', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.phitigra import Phitigra
            sage: isinstance(Phitigra(), Phitigra)
            True
        """
        PythonModule.__init__(self, 'phitigra', spkg='phitigra')


def all_features():
    return [Phitigra()]
