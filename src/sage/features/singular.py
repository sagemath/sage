r"""
Features for testing the presence of Singular
"""
from . import Executable
from sage.env import SINGULAR_BIN


class Singular(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the Singular executable.

    EXAMPLES::

        sage: from sage.features.singular import Singular
        sage: Singular().is_present()
        FeatureTestResult('singular', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.singular import Singular
            sage: isinstance(Singular(), Singular)
            True
        """
        Executable.__init__(self, "singular", SINGULAR_BIN,
                            spkg='singular')
