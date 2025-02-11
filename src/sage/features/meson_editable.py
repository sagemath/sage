r"""
Feature for testing if Meson editable install is used.
"""

import sys
from . import Feature, FeatureTestResult


class MesonEditable(Feature):
    r"""
    A :class:`~sage.features.Feature` describing if Meson editable install is used.

    EXAMPLES::

        sage: from sage.features.meson_editable import MesonEditable
        sage: MesonEditable()
        Feature('meson_editable')
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.meson_editable import MesonEditable
            sage: MesonEditable() is MesonEditable()
            True
        """
        Feature.__init__(self, 'meson_editable')

    def _is_present(self):
        r"""
        Test whether Meson editable install is used.

        EXAMPLES::

            sage: from sage.features.meson_editable import MesonEditable
            sage: Internet()._is_present()  # random
            FeatureTestResult('meson_editable', True)
        """
        import sage
        import sys
        result = type(sage.__loader__).__module__ == '_sagemath_editable_loader'
        return FeatureTestResult(self, result)


def all_features():
    return [MesonEditable()]
