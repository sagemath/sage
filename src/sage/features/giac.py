# sage_setup: distribution = sagemath-environment
# sage.doctest: optional - giac
r"""
Feature for testing the presence of ``giac``
"""

from . import Executable, FeatureTestResult
from sage.env import SAGE_GIAC_ENABLED

class Giac(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`giac <spkg_giac>`.

    EXAMPLES::

        sage: from sage.features.giac import Giac
        sage: Giac().is_present()  # needs giac
        FeatureTestResult('giac', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.giac import Giac
            sage: isinstance(Giac(), Giac)
            True
        """
        if SAGE_GIAC_ENABLED == "no":
            giac_exe = 'fofobar42barfoo'
        else:
            giac_exe = 'giac'
        Executable.__init__(self, 'giac', executable=giac_exe,
                            spkg='giac', type='optional')

def all_features():
    return [Giac()]
