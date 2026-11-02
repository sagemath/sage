r"""
Feature for the m4rie library (matrices over `GF(2^e)`)
"""

from sage.config import m4rie_enabled
from sage.features.build_feature import BuildModule


class M4rie(BuildModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    :mod:`sage.matrix.matrix_gf2e_dense`, based on m4rie.

    EXAMPLES::

        sage: from sage.features.m4rie import M4rie
        sage: M4rie().is_present()  # needs m4rie
        FeatureTestResult('m4rie', True)
        sage: M4rie().is_present()  # needs !m4rie
        FeatureTestResult('m4rie', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !m4rie`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.m4rie import M4rie
        sage: M4rie().is_present_at_runtime()  # needs m4rie
        FeatureTestResult('m4rie', True)

    """
    _enabled_in_build = m4rie_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.m4rie import M4rie
            sage: M4rie()
            Feature('m4rie')

        """
        module_name = "sage.matrix.matrix_gf2e_dense"
        super().__init__('m4rie',
                         module_name,
                         type='standard')


def all_features():
    return [M4rie()]
