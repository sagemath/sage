# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of the ``mwrank`` program (part
of eclib)
"""

from . import Executable, FeatureTestResult


class Mwrank(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the ``mwrank`` program.

    EXAMPLES::

        sage: from sage.features.mwrank import Mwrank
        sage: Mwrank().is_present()  # needs mwrank
        FeatureTestResult('mwrank', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.mwrank import Mwrank
            sage: isinstance(Mwrank(), Mwrank)
            True
        """
        Executable.__init__(self, 'mwrank', executable='mwrank',
                            spkg='eclib', type='standard')


def all_features():
    return [Mwrank()]
