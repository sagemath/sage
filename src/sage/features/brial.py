r"""
Feature for testing the presence of (lib)brial
"""

from sage.features.join_feature import JoinFeature
from sage.features import PythonModule


class Brial(JoinFeature):
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
        FeatureTestResult('sage.rings.polynomial.pbori.pbori', False)

    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.brial import Brial
            sage: isinstance(Brial(), Brial)
            True
        """
        JoinFeature.__init__(self, 'brial',
                             [PythonModule('sage.rings.polynomial.pbori.pbori')],
                             spkg='sagemath_brial', type='standard')


def all_features():
    return [Brial()]
