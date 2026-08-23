r"""
Feature for testing the presence of primesieve
"""

from sage.config import primesieve_enabled
from sage.features.build_feature import BuildModule


class Primesieve(BuildModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    :mod:`sage.libs.primesieve`, the interface to the primesieve
    library.

    EXAMPLES::

        sage: from sage.features.primesieve import Primesieve
        sage: Primesieve().is_present()  # needs primesieve
        FeatureTestResult('primesieve', True)
        sage: Primesieve().is_present()  # needs !primesieve
        FeatureTestResult('primesieve', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !primesieve`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.primesieve import Primesieve
        sage: Primesieve().is_present_at_runtime()  # needs primesieve
        FeatureTestResult('primesieve', True)
    """
    _enabled_in_build = primesieve_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.primesieve import Primesieve
            sage: Primesieve()
            Feature('primesieve')
        """
        module_name = "sage.libs.primesieve"
        super().__init__('primesieve',
                         module_name,
                         type='standard')


def all_features():
    return [Primesieve()]
