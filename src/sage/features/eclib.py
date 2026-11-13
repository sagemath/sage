r"""
Feature for eclib (the :mod:`sage.libs.eclib.mwrank` module)
"""

from sage.config import eclib_enabled
from sage.features.build_feature import BuildModule


class Eclib(BuildModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    :mod:`sage.libs.eclib.mwrank`.

    EXAMPLES::

        sage: from sage.features.eclib import Eclib
        sage: Eclib().is_present()  # needs eclib
        FeatureTestResult('eclib', True)
        sage: Eclib().is_present()  # needs !eclib
        FeatureTestResult('eclib', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !eclib`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.eclib import Eclib
        sage: Eclib().is_present_at_runtime()  # needs eclib
        FeatureTestResult('eclib', True)

    """
    _enabled_in_build = eclib_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.eclib import Eclib
            sage: Eclib()
            Feature('eclib')

        """
        module_name = "sage.libs.eclib.mwrank"
        super().__init__("eclib", module_name, type="standard")


def all_features():
    return [Eclib()]
