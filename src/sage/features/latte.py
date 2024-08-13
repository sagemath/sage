# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``latte_int``
"""

# ****************************************************************************
#       Copyright (C) 2018 Vincent Delecroix
#                     2019 Frédéric Chapoton
#                     2021 Matthias Koeppe
#                     2021 Kwankyu Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable
from .join_feature import JoinFeature


LATTE_URL = "https://www.math.ucdavis.edu/~latte/software.php"


class Latte_count(Executable):
    r"""
    Feature for the executable ``count`` from :ref:`LattE integrale <spkg_latte_int>`.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latte import Latte_count
            sage: isinstance(Latte_count(), Latte_count)
            True
        """
        Executable.__init__(self, 'count', executable='count',
                            spkg='latte_int',
                            url=LATTE_URL)


class Latte_integrate(Executable):
    r"""
    Feature for the executable ``integrate`` from :ref:`LattE integrale <spkg_latte_int>`.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latte import Latte_integrate
            sage: isinstance(Latte_integrate(), Latte_integrate)
            True
        """
        Executable.__init__(self, 'integrate', executable='integrate',
                            spkg='latte_int',
                            url=LATTE_URL)


class Latte(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of excecutables
    from :ref:`LattE integrale <spkg_latte_int>`.

    EXAMPLES::

        sage: from sage.features.latte import Latte
        sage: Latte().is_present()  # optional - latte_int
        FeatureTestResult('latte_int', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latte import Latte
            sage: isinstance(Latte(), Latte)
            True
        """
        JoinFeature.__init__(self, 'latte_int',
                             (Latte_count(), Latte_integrate()),
                             description='LattE')


def all_features():
    return [Latte()]
