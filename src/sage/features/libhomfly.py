r"""
Feature for testing the presence of libhomfly
"""

from sage.config import libhomfly_enabled
from sage.features.build_feature import BuildModule


class Libhomfly(BuildModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    ``sage.libs.homfly``, the interface to libhomfly.

    EXAMPLES::

        sage: from sage.features.libhomfly import Libhomfly
        sage: Libhomfly().is_present()  # needs libhomfly
        FeatureTestResult('libhomfly', True)
        sage: Libhomfly().is_present()  # needs !libhomfly
        FeatureTestResult('libhomfly', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !libhomfly`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.libhomfly import Libhomfly
        sage: Libhomfly().is_present_at_runtime()  # needs libhomfly
        FeatureTestResult('libhomfly', True)

    """
    _enabled_in_build = libhomfly_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.libhomfly import Libhomfly
            sage: Libhomfly()
            Feature('libhomfly')

        """
        module_name = "sage.libs.homfly"
        super().__init__('libhomfly',
                         module_name,
                         type='standard')


def all_features():
    return [Libhomfly()]
