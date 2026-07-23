# sage_setup: distribution = sagemath-environment
r"""
Check for pygambit
"""
# ****************************************************************************
#       Copyright (C) 2025 SageMath Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import PythonModule


class pygambit(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of the
    Python package :ref:`pygambit <spkg_pygambit>`.

    EXAMPLES::

        sage: from sage.features.gambit import pygambit
        sage: pygambit().is_present()                           # optional - pygambit
        FeatureTestResult('pygambit', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.gambit import pygambit
            sage: isinstance(pygambit(), pygambit)
            True
        """
        PythonModule.__init__(self, 'pygambit', spkg='pygambit',
                              url='https://gambitproject.readthedocs.io/en/latest/install.html#install')

    def resolution(self):
        r"""
        Return a suggestion on how to make :meth:`is_present` pass if it did
        not pass.

        The message gives the command for a non-managed installation of Sage
        (plain ``pip``), the command for a managed installation of Sage such as
        one based on conda or on system packages (``sage --pip``), and a pointer
        to the upstream gambit installation instructions.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.features.gambit import pygambit
            sage: print(pygambit().resolution())
            To install pygambit you can run one of the following:
              * in a non-managed installation of Sage:  pip install pygambit
              * in a managed installation of Sage...:  sage --pip install pygambit
            Further installation instructions are available at
            https://gambitproject.readthedocs.io/en/latest/install.html#install.
        """
        if self._hidden:
            return super().resolution()
        return (
            "To install pygambit you can run one of the following:\n"
            "  * in a non-managed installation of Sage:  pip install pygambit\n"
            "  * in a managed installation of Sage (for example based on conda "
            "or on system packages):  sage --pip install pygambit\n"
            "Further installation instructions are available at\n"
            "{url}.".format(url=self.url)
        )


def all_features():
    return [pygambit()]