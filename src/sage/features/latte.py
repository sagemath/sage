# -*- coding: utf-8 -*-
r"""
Features for testing the presence of ``latte_int``
"""
from . import Executable
from .join_feature import JoinFeature


LATTE_URL = "https://www.math.ucdavis.edu/~latte/software.php"


class Latte_count(Executable):
    r"""
    Feature for the executable ``count`` from the LattE suite.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latte import Latte_count
            sage: isinstance(Latte_count(), Latte_count)
            True
        """
        Executable.__init__(self, "count", executable="count",
                            spkg="latte_int",
                            url=LATTE_URL)


class Latte_integrate(Executable):
    r"""
    Feature for the executable ``integrate`` from the LattE suite.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latte import Latte_integrate
            sage: isinstance(Latte_integrate(), Latte_integrate)
            True
        """
        Executable.__init__(self, "integrate", executable="integrate",
                            spkg="latte_int",
                            url=LATTE_URL)


class Latte(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the ``LattE``
    binaries which comes as a part of ``latte_int``.

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
        JoinFeature.__init__(self, "latte_int",
                             (Latte_count(), Latte_integrate()),
                             description="LattE")


def all_features():
    return [Latte()]
