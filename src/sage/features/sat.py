# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of SAT solvers
"""

from . import Executable, PythonModule


class Glucose(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of an
    executable from the :ref:`Glucose SAT solver <spkg_glucose>`.

    EXAMPLES::

        sage: from sage.features.sat import Glucose
        sage: Glucose().is_present()                  # optional - glucose
        FeatureTestResult('glucose', True)
    """
    def __init__(self, executable="glucose"):
        r"""
        TESTS::

            sage: from sage.features.sat import Glucose
            sage: isinstance(Glucose(), Glucose)
            True
        """
        Executable.__init__(self, name=executable, executable=executable,
                            spkg="glucose", type="optional")


class Kissat(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the
    :ref:`Kissat SAT solver <spkg_kissat>`.

    EXAMPLES::

        sage: from sage.features.sat import Kissat
        sage: Kissat().is_present()                             # optional - kissat
        FeatureTestResult('kissat', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sat import Kissat
            sage: isinstance(Kissat(), Kissat)
            True
        """
        Executable.__init__(self, name="kissat", executable="kissat",
                            spkg="kissat", type="optional")


class Pycosat(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`spkg_pycosat`.

    EXAMPLES::

        sage: from sage.features.sat import Pycosat
        sage: Pycosat().is_present()                  # optional - pycosat
        FeatureTestResult('pycosat', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sat import Pycosat
            sage: isinstance(Pycosat(), Pycosat)
            True
        """
        PythonModule.__init__(self, "pycosat",
                              spkg="pycosat", type="optional")


class Pycryptosat(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`spkg_pycryptosat`.

    EXAMPLES::

        sage: from sage.features.sat import Pycryptosat
        sage: Pycryptosat().is_present()              # optional - pycryptosat
        FeatureTestResult('pycryptosat', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sat import Pycryptosat
            sage: isinstance(Pycryptosat(), Pycryptosat)
            True
        """
        PythonModule.__init__(self, "pycryptosat",
                              spkg="pycryptosat", type="optional")


def all_features():
    return [Glucose(),
            Kissat(),
            Pycosat(),
            Pycryptosat()]
