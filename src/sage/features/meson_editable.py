r"""
Feature for testing if Meson editable install is used.
"""
from sage.config import is_editable_install

from . import Feature, FeatureTestResult


class MesonEditable(Feature):
    r"""
    A :class:`~sage.features.Feature` describing if Meson editable install
    is used.

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
            sage: MesonEditable()._is_present()  # random
            FeatureTestResult('meson_editable', True)
        """
        return FeatureTestResult(self, is_editable_install())


def all_features():
    return [MesonEditable()]
