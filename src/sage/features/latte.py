# -*- coding: utf-8 -*-
r"""
Check for LattE
"""
from . import Executable, Feature, FeatureTestResult

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


class Latte(Feature):
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
        Feature.__init__(self, "LattE")

    def _is_present(self):
        r"""
        Test for the presence of LattE binaries.

        EXAMPLES::

            sage: from sage.features.latte import Latte
            sage: Latte()._is_present()  # optional - latte_int
            FeatureTestResult('LattE', True)
        """

        test = (Latte_count()._is_present() and
                Latte_integrate()._is_present())
        if not test:
            return test

        return FeatureTestResult(self, True)

    def is_functional(self):
        r"""
        Test whether count and integrate are functionals.

        EXAMPLES::

            sage: from sage.features.latte import Latte
            sage: Latte().is_functional()  # optional - latte_int
            FeatureTestResult('LattE', True)
        """
        test = (Latte_count().is_functional() and
                Latte_integrate().is_functional())
        if not test:
            return test

        return FeatureTestResult(self, True)
