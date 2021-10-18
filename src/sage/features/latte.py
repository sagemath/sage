# -*- coding: utf-8 -*-
r"""
Check for LattE
"""
from . import Executable, Feature, FeatureTestResult
from .join_feature import JoinFeature


LATTE_URL = "https://www.math.ucdavis.edu/~latte/software.php"


class Latte_count(Executable):
    r"""
    Feature for the executable ``count`` from the LattE suite.
    """
    def __init__(self):
        Executable.__init__(self, "count", executable="count",
                            spkg="latte_int",
                            url=LATTE_URL)


class Latte_integrate(Executable):
    r"""
    Feature for the executable ``integrate`` from the LattE suite.
    """
    def __init__(self):
        Executable.__init__(self, "integrate", executable="integrate",
                            spkg="latte_int",
                            url=LATTE_URL)


class Latte(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of the ``LattE``
    binaries which comes as a part of ``latte_int``.

    EXAMPLES::

        sage: from sage.features.latte import Latte
        sage: Latte().is_present()  # optional - latte_int
        FeatureTestResult('LattE', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latte import Latte
            sage: isinstance(Latte(), Latte)
            True
        """
        JoinFeature.__init__(self, "LattE",
                             (Latte_count(), Latte_integrate()))
