# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``rubiks``
"""
# ****************************************************************************
#       Copyright (C) 2020      John H. Palmieri
#                     2021-2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.env import RUBIKS_BINS_PREFIX

from . import Executable
from .join_feature import JoinFeature


class cu2(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``cu2``.

    EXAMPLES::

        sage: from sage.features.rubiks import cu2
        sage: cu2().is_present()  # optional - rubiks
        FeatureTestResult('cu2', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.rubiks import cu2
            sage: isinstance(cu2(), cu2)
            True
        """
        Executable.__init__(self, "cu2", executable=RUBIKS_BINS_PREFIX + "cu2",
                            spkg='rubiks')


class size222(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``size222``.

    EXAMPLES::

        sage: from sage.features.rubiks import size222
        sage: size222().is_present()  # optional - rubiks
        FeatureTestResult('size222', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.rubiks import size222
            sage: isinstance(size222(), size222)
            True
        """
        Executable.__init__(self, "size222", executable=RUBIKS_BINS_PREFIX + "size222",
                            spkg='rubiks')


class optimal(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``optimal``.

    EXAMPLES::

        sage: from sage.features.rubiks import optimal
        sage: optimal().is_present()  # optional - rubiks
        FeatureTestResult('optimal', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.rubiks import optimal
            sage: isinstance(optimal(), optimal)
            True
        """
        Executable.__init__(self, "optimal", executable=RUBIKS_BINS_PREFIX + "optimal",
                            spkg='rubiks')


class mcube(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``mcube``.

    EXAMPLES::

        sage: from sage.features.rubiks import mcube
        sage: mcube().is_present()  # optional - rubiks
        FeatureTestResult('mcube', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.rubiks import mcube
            sage: isinstance(mcube(), mcube)
            True
        """
        Executable.__init__(self, "mcube", executable=RUBIKS_BINS_PREFIX + "mcube",
                            spkg='rubiks')


class dikcube(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``dikcube``.

    EXAMPLES::

        sage: from sage.features.rubiks import dikcube
        sage: dikcube().is_present()  # optional - rubiks
        FeatureTestResult('dikcube', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.rubiks import dikcube
            sage: isinstance(dikcube(), dikcube)
            True
        """
        Executable.__init__(self, "dikcube", executable=RUBIKS_BINS_PREFIX + "dikcube",
                            spkg='rubiks')


class cubex(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``cubex``.

    EXAMPLES::

        sage: from sage.features.rubiks import cubex
        sage: cubex().is_present()  # optional - rubiks
        FeatureTestResult('cubex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.rubiks import cubex
            sage: isinstance(cubex(), cubex)
            True
        """
        Executable.__init__(self, "cubex", executable=RUBIKS_BINS_PREFIX + "cubex",
                            spkg='rubiks')


class Rubiks(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the
    :class:`cu2`, :class:`cubex`, :class:`dikcube`, :class:`mcube`, :class:`optimal`, and
    :class:`size222` programs from the :ref:`rubiks <spkg_rubiks>` package.

    EXAMPLES::

        sage: from sage.features.rubiks import Rubiks
        sage: Rubiks().is_present()  # optional - rubiks
        FeatureTestResult('rubiks', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.rubiks import Rubiks
            sage: isinstance(Rubiks(), Rubiks)
            True
        """
        JoinFeature.__init__(self, "rubiks",
                             [cu2(), size222(), optimal(), mcube(), dikcube(), cubex()],
                             spkg='rubiks')


def all_features():
    return [Rubiks()]
