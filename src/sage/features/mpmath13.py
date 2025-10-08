r"""
Feature for test for mpmath-1.3.x. Used to simultaneously
support v1.3 and v1.4 in the doctests.
"""
from . import Feature


class MpMath13(Feature):
    r"""
    A :class:`~sage.features.Feature`.

    EXAMPLES::

        sage: from sage.features.mpmath13 import MpMath13
        sage: MpMath13()
        Feature('mpmath13')
        sage: MpMath13() is MpMath13()
        True
    """
    def __init__(self):
        Feature.__init__(self, 'mpmath13')

    def _is_present(self):
        r"""
        Test the mpmath version.

        EXAMPLES::

            sage: from sage.features.mpmath13 import MpMath13
            sage: MpMath()._is_present()  # needs mpmath13
            FeatureTestResult('mpmath13', True)

        """
        from mpmath import __version__ as mpmath_version
        return mpmath_version.startswith("1.3")


def all_features():
    return [MpMath13()]
