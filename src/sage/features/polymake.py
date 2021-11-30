from . import PythonModule
from .join_feature import JoinFeature


class JuPyMake(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of the ``JuPyMake``
    module, a Python interface to the polymake library.

    EXAMPLES::

        sage: from sage.features.polymake import JuPyMake
        sage: JuPyMake().is_present()  # optional: jupymake
        FeatureTestResult('jupymake', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.polymake import JuPyMake
            sage: isinstance(JuPyMake(), JuPyMake)
            True
        """
        JoinFeature.__init__(self, "jupymake",
                             [PythonModule("JuPyMake", spkg="jupymake")])
