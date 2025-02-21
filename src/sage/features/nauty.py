# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of nauty executables
"""

# *****************************************************************************
#       Copyright (C) 2022 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.env import SAGE_NAUTY_BINS_PREFIX

from . import Executable
from .join_feature import JoinFeature


class NautyExecutable(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for executables from the :ref:`nauty <spkg_nauty>` package.

    EXAMPLES::

        sage: from sage.features.nauty import NautyExecutable
        sage: NautyExecutable('converseg').is_present()                                 # needs nauty
        FeatureTestResult('nauty_converseg', True)
    """
    def __init__(self, name):
        r"""
        TESTS::

            sage: from sage.features.nauty import NautyExecutable
            sage: isinstance(NautyExecutable('geng'), NautyExecutable)
            True
        """
        Executable.__init__(self, name=f"nauty_{name}",
                            executable=f"{SAGE_NAUTY_BINS_PREFIX}{name}",
                            spkg='nauty',
                            type='standard')


class Nauty(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the executables
    which comes as a part of :ref:`nauty <spkg_nauty>`.

    EXAMPLES::

        sage: from sage.features.nauty import Nauty
        sage: Nauty().is_present()                                                      # needs nauty
        FeatureTestResult('nauty', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.nauty import Nauty
            sage: isinstance(Nauty(), Nauty)
            True
        """
        JoinFeature.__init__(self, "nauty",
                             [NautyExecutable(name)
                              for name in ('directg', 'gentourng', 'geng', 'genbg', 'gentreeg', 'genktreeg')])


def all_features():
    return [Nauty()]
