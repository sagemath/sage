# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``palp``
"""
# ****************************************************************************
#       Copyright (C) 2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable
from .join_feature import JoinFeature


class PalpExecutable(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of a :ref:`PALP <spkg_palp>` executable.

    INPUT:

    - ``palpprog`` -- string, one of ``'poly'``, ``'class'``, ``'nef'``, ``'cws'``

    - ``suff`` -- string or ``None``
    """
    def __init__(self, palpprog, suff=None):
        r"""
        TESTS::

            sage: from sage.features.palp import PalpExecutable
            sage: isinstance(PalpExecutable("poly", 5), PalpExecutable)
            True
        """
        if suff:
            Executable.__init__(self, f"palp_{palpprog}_{suff}d",
                                executable=f"{palpprog}-{suff}d.x",
                                spkg='palp', type='standard')
        else:
            Executable.__init__(self, f"palp_{palpprog}",
                                executable=f"{palpprog}.x",
                                spkg='palp', type='standard')


class Palp(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`PALP <spkg_palp>`.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.palp import Palp
            sage: isinstance(Palp(), Palp)
            True
        """
        JoinFeature.__init__(self, "palp",
                             [PalpExecutable(palpprog, suff)
                              for palpprog in ("poly", "class", "nef", "cws")
                              for suff in (None, 4, 5, 6, 11)],
                             description='PALP')


def all_features():
    return [Palp()]
