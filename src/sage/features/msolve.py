# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of msolve

`msolve <https://msolve.lip6.fr/>`_ is a multivariate polynomial system solver.

.. SEEALSO::

    - :mod:`sage.rings.polynomial.msolve`
"""

# *****************************************************************************
#       Copyright (C) 2022 Marc Mezzarobba
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import subprocess
from . import Executable
from . import FeatureTestResult

class msolve(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`msolve <spkg_msolve>`.

    EXAMPLES::

        sage: from sage.features.msolve import msolve
        sage: msolve().is_present()  # optional - msolve
        FeatureTestResult('msolve', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.msolve import msolve
            sage: isinstance(msolve(), msolve)
            True
        """
        Executable.__init__(self, "msolve", executable='msolve',
                            url='https://msolve.lip6.fr/')

    def is_functional(self):
        r"""
        Test if our installation of msolve is working.

        TESTS::

            sage: from sage.features.msolve import msolve
            sage: msolve().is_functional()  # optional - msolve
            FeatureTestResult('msolve', True)
        """
        msolve_out = subprocess.run(["msolve", "-h"], capture_output=True)

#        if msolve_out.returncode != 0:
#            return FeatureTestResult(self, False, reason="msolve -h returned "
#                                f"nonzero exit status {msolve_out.returncode}")
        if (msolve_out.stdout[:46] !=
              b'\nmsolve library for polynomial system solving\n'):
            return FeatureTestResult(self, False,
                                     reason="output of msolve -h not recognized")
        return FeatureTestResult(self, True)

def all_features():
    return [msolve()]
