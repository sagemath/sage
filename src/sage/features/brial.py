r"""
Feature for testing the presence of (lib)brial
"""

from sage.config import brial_enabled
from sage.features.build_feature import BuildModule


class Brial(BuildModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    :mod:`~sage.rings.polynomial.pbori.pbori`.

    The :mod:`~sage.rings.polynomial.pbori.pbori` module in turn depends on
    the presence and usability of libbrial -- a slightly more
    modern fork of PolyBoRi, which hopefully explains the name.

    EXAMPLES::

        sage: from sage.features.brial import Brial
        sage: Brial().is_present()  # needs brial
        FeatureTestResult('brial', True)
        sage: Brial().is_present()  # needs !brial
        FeatureTestResult('brial', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !brial`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.brial import Brial
        sage: Brial().is_present_at_runtime()  # needs brial
        FeatureTestResult('brial', True)

    """
    _enabled_in_build = brial_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.brial import Brial
            sage: Brial()
            Feature('brial')

        """
        module_name = "sage.rings.polynomial.pbori.pbori"
        super().__init__("brial",
                         module_name,
                         type="standard")


def all_features():
    return [Brial()]
