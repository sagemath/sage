r"""
Feature for testing the presence of libhomfly
"""

from sage.features.join_feature import JoinFeature
from sage.features import PythonModule


class Libhomfly(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    :mod:`sage.libs.homfly`, the interface to libhomfly.

    TESTS::

        sage: from sage.features.libhomfly import Libhomfly
        sage: Libhomfly().is_present()  # needs libhomfly
        FeatureTestResult('libhomfly', True)
        sage: Libhomfly().is_present()  # needs !libhomfly
        FeatureTestResult('sage.libs.homfly', False)

    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.libhomfly import Libhomfly
            sage: isinstance(Libhomfly(), Libhomfly)
            True

        """
        JoinFeature.__init__(self, 'libhomfly',
                             [PythonModule('sage.libs.homfly')],
                             spkg='libhomfly', type='standard')


def all_features():
    return [Libhomfly()]
