r"""
Feature for testing the presence of (lib)brial
"""

from sage.config import brial_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature


class Brial(BuildFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    :mod:`sage.rings.polynomial.pbori`.

    The :mod:`sage.rings.polynomial.pbori` module in turn depends on
    the presence and usability of libbrial -- a slightly more
    modern fork of PolyBoRi, which hopefully explains the name.

    TESTS::

        sage: from sage.features.brial import Brial
        sage: Brial().is_present()  # needs brial
        FeatureTestResult('brial', True)
        sage: Brial().is_present()  # needs !brial
        FeatureTestResult('brial', False)

    """
    _enabled_in_build = brial_enabled

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.brial import Brial
            sage: isinstance(Brial(), Brial)
            True
        """
        super().__init__("brial", spkg="brial", type="standard")


    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.brial import Brial
            sage: brial = Brial()
            sage: result = brial.is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs brial
            FeatureTestResult('brial', True)

        """
        # The build system installs the sage/rings/pbori source code
        # even when brial support is disabled, so a naive import of
        # that module will actually succeed. We check for a
        # conditionally-compiled cython module instead.
        cython_modname = "sage.rings.polynomial.pbori.pbori"
        result = PythonModule(cython_modname)._is_present()
        result.feature = self
        return result

def all_features():
    return [Brial()]
