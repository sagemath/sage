r"""
Feature for testing the presence of libbraiding
"""

from sage.config import libbraiding_enabled
from sage.features.build_feature import BuildModule


class Libbraiding(BuildModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    ``sage.libs.braiding``, the interface to libbraiding.

    EXAMPLES::

        sage: from sage.features.libbraiding import Libbraiding
        sage: Libbraiding().is_present()  # needs libbraiding
        FeatureTestResult('libbraiding', True)
        sage: Libbraiding().is_present()  # needs !libbraiding
        FeatureTestResult('libbraiding', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !libbraiding`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.libbraiding import Libbraiding
        sage: Libbraiding().is_present_at_runtime()  # needs libbraiding
        FeatureTestResult('libbraiding', True)

    """
    _enabled_in_build = libbraiding_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.libbraiding import Libbraiding
            sage: Libbraiding()
            Feature('libbraiding')

        """
        module_name = "sage.libs.braiding"
        super().__init__('libbraiding',
                         module_name,
                         type='standard')


def all_features():
    return [Libbraiding()]
