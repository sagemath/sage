r"""
Feature for testing if the machine is 32-bit.
"""
from . import Feature, FeatureTestResult
import sys


class Is32Bit(Feature):
    r"""
    A :class:`~sage.features.Feature`.

    EXAMPLES::

        sage: from sage.features.bitness import Is32Bit
        sage: Is32Bit()
        Feature('32_bit')
        sage: Is32Bit() is Is32Bit()
        True
    """
    def __init__(self):
        Feature.__init__(self, '32_bit')

    def _is_present(self):
        r"""
        Test whether the machine is 32-bit.

        EXAMPLES::

            sage: from sage.features.bitness import Is32Bit
            sage: Is32Bit()._is_present()  # random
            FeatureTestResult('32_bit', False)
        """
        return sys.maxsize <= (1 << 32)  # cf. sage.doctest.sources.bitness_value


def all_features():
    return [Is32Bit()]
