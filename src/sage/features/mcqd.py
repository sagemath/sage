r"""
Features for testing the presence of ``mcqd``
"""

from . import PythonModule
from .join_feature import JoinFeature


class Mcqd(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :mod:`~sage.graphs.mcqd`

    EXAMPLES::

        sage: from sage.features.mcqd import Mcqd
        sage: Mcqd().is_present()  # optional - mcqd
        FeatureTestResult('mcqd', True)
    """

    def __init__(self):
        """
        TESTS::

            sage: from sage.features.mcqd import Mcqd
            sage: isinstance(Mcqd(), Mcqd)
            True
        """
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_mcqd' later
        JoinFeature.__init__(self, 'mcqd',
                             [PythonModule('sage.graphs.mcqd', spkg='mcqd')])


def all_features():
    return [Mcqd()]
