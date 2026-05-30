r"""
Feature for testing the presence of libhomfly
"""

from sage.config import libhomfly_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature


class Libhomfly(BuildFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    :mod:`sage.libs.homfly`, the interface to libhomfly.

    TESTS::

        sage: from sage.features.libhomfly import Libhomfly
        sage: Libhomfly().is_present()  # needs libhomfly
        FeatureTestResult('libhomfly', True)
        sage: Libhomfly().is_present()  # needs !libhomfly
        FeatureTestResult('libhomfly', False)

    """
    _enabled_in_build = libhomfly_enabled

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.libhomfly import Libhomfly
            sage: isinstance(Libhomfly(), Libhomfly)
            True

        """
        super().__init__('libhomfly', spkg='libhomfly', type='standard')

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.libhomfly import Libhomfly
            sage: lhf = Libhomfly()
            sage: result = lhf.is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs libhomfly
            FeatureTestResult('libhomfly', True)

        """
        result = PythonModule("sage.libs.homfly")._is_present()
        result.feature = self
        return result


def all_features():
    return [Libhomfly()]
