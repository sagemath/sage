r"""
Feature for testing the presence of ``cddlib``
"""

import subprocess

from . import Executable, FeatureTestResult


class CddExecutable(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of an executable
    which comes as a part of ``cddlib``.

    EXAMPLES::

        sage: from sage.features.cddlib import CddExecutable
        sage: CddExecutable().is_present()
        FeatureTestResult('cddexec_gmp', True)
    """
    def __init__(self, name='cddexec_gmp'):
        r"""
        TESTS::

            sage: from sage.features.cddlib import CddExecutable
            sage: isinstance(CddExecutable(), CddExecutable)
            True
        """
        Executable.__init__(self, name=name, executable=name, spkg="cddlib",
                            url="https://github.com/cddlib/cddlib")
