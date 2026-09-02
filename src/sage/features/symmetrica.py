r"""
Features for the symmetrica library
"""

from sage.config import symmetrica_enabled
from sage.features.build_feature import BuildModule

class Symmetrica(BuildModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the SageMath interface to the symmetrica library.

    EXAMPLES::

        sage: from sage.features.symmetrica import Symmetrica
        sage: Symmetrica().is_present()  # needs symmetrica
        FeatureTestResult('symmetrica', True)
        sage: Symmetrica().is_present()  # needs !symmetrica
        FeatureTestResult('symmetrica', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !symmetrica`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.symmetrica import Symmetrica
        sage: Symmetrica().is_present_at_runtime()  # needs symmetrica
        FeatureTestResult('symmetrica', True)

    """
    _enabled_in_build = symmetrica_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.symmetrica import Symmetrica
            sage: Symmetrica()
            Feature('symmetrica')

        """
        module_name = "sage.libs.symmetrica.symmetrica"
        super().__init__("symmetrica",
                         module_name,
                         spkg="symmetrica",
                         type="standard")


def all_features():
    return [Symmetrica()]
