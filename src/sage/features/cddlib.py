# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``cddlib``
"""

# *****************************************************************************
#       Copyright (C) 2022 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import Executable


class CddExecutable(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of an executable
    which comes as a part of :ref:`cddlib <spkg_cddlib>`.

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
        Executable.__init__(self, name=name, executable=name, spkg='cddlib',
                            url='https://github.com/cddlib/cddlib', type='standard')
