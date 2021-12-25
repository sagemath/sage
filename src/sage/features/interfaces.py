r"""
Features for testing whether interpreter interfaces are functional
"""
import importlib

from . import Feature, FeatureTestResult, PythonModule


class InterfaceFeature(Feature):
    r"""
    A :class:`~sage.features.Feature` describing whether an :class:`~sage.interfaces.interface.Interface` is present and functional.

    TESTS::

        sage: from sage.features.interfaces import InterfaceFeature
        sage: broken = InterfaceFeature("broken_interface", "sage.interfaces.does_not_exist")
        sage: broken.is_present()
        FeatureTestResult('sage.interfaces.does_not_exist', False)
        sage: _.reason
        "Failed to import `sage.interfaces.does_not_exist`: No module named 'sage.interfaces.does_not_exist'"

        sage: also_broken = InterfaceFeature("also_broken_interface", "sage.interfaces.interface")
        sage: also_broken.is_present()
        FeatureTestResult('also_broken_interface', False)
        sage: _.reason
        "Interface also_broken_interface cannot be imported: module 'sage.interfaces.interface' has no attribute 'also_broken_interface'"
    """
    @staticmethod
    def __classcall__(cls, name, module, description=None):
        """
        TESTS::

            sage: from sage.features import PythonModule
            sage: from sage.features.interfaces import InterfaceFeature
            sage: f = InterfaceFeature("test_interface", "sage.interfaces.interface")
            sage: f is InterfaceFeature("test_interface", PythonModule("sage.interfaces.interface"))
            True
        """
        if isinstance(module, str):
            module = PythonModule(module)
        return Feature.__classcall__(cls, name, module, description)

    def __init__(self, name, module, description):
        """
        TESTS::

            sage: from sage.features.interfaces import InterfaceFeature
            sage: f = InterfaceFeature("test_interface", "sage.interfaces.interface")
            sage: isinstance(f, InterfaceFeature)
            True
        """
        super().__init__(name, description=description)
        self.module = module

    def _is_present(self):
        """
        TESTS::

            sage: from sage.features.interfaces import InterfaceFeature
            sage: from sage.interfaces.sage0 import Sage
            sage: f = InterfaceFeature("sage0", "sage.interfaces.sage0")
            sage: f.is_present()
            FeatureTestResult('sage0', True)
        """
        result = self.module.is_present()
        if not result:
            return result
        m = importlib.import_module(self.module.name)
        try:
            interface = getattr(m, self.name)
        except Exception as exception:
            return FeatureTestResult(self, False,
                                     reason=f"Interface {self.name} cannot be imported: {exception}")
        try:
            interface('2+3')
            return FeatureTestResult(self, True)
        except Exception as exception:
            return FeatureTestResult(self, False,
                                     reason=f"Interface {interface} is not functional: {exception}")


# The following are provided by external software only (no SPKG)

class Magma(InterfaceFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether :class:`sage.interfaces.magma.Magma`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Magma
        sage: Magma().is_present()  # random
        FeatureTestResult('magma', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature.__classcall__(cls, 'magma', 'sage.interfaces.magma')


class Matlab(InterfaceFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether :class:`sage.interfaces.matlab.Matlab`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Matlab
        sage: Matlab().is_present()  # random
        FeatureTestResult('matlab', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature.__classcall__(cls, 'matlab', 'sage.interfaces.matlab')


class Mathematica(InterfaceFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether :class:`sage.interfaces.mathematica.Mathematica`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Mathematica
        sage: Mathematica().is_present()  # random
        FeatureTestResult('mathematica', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature.__classcall__(cls, 'mathematica', 'sage.interfaces.mathematica')


class Maple(InterfaceFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether :class:`sage.interfaces.maple.Maple`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Maple
        sage: Maple().is_present()  # random
        FeatureTestResult('maple', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature.__classcall__(cls, 'maple', 'sage.interfaces.maple')


class Macaulay2(InterfaceFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether :class:`sage.interfaces.macaulay2.Macaulay2`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Macaulay2
        sage: Macaulay2().is_present()  # random
        FeatureTestResult('macaulay2', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature.__classcall__(cls, 'macaulay2', 'sage.interfaces.macaulay2')


class Octave(InterfaceFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether :class:`sage.interfaces.octave.Octave`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Octave
        sage: Octave().is_present()  # random
        FeatureTestResult('octave', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature.__classcall__(cls, 'octave', 'sage.interfaces.octave')


class Scilab(InterfaceFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether :class:`sage.interfaces.scilab.Scilab`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Scilab
        sage: Scilab().is_present()  # random
        FeatureTestResult('scilab', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature.__classcall__(cls, 'scilab', 'sage.interfaces.scilab')


def all_features():
    r"""
    Return features corresponding to interpreter interfaces.

    EXAMPLES::

        sage: from sage.features.interfaces import all_features
        sage: list(all_features())
        [Feature('magma'),
         Feature('matlab'),
         Feature('mathematica'),
         Feature('maple'),
         Feature('macaulay2'),
         Feature('octave'),
         Feature('scilab')]
    """
    return [Magma(),
            Matlab(),
            Mathematica(),
            Maple(),
            Macaulay2(),
            Octave(),
            Scilab()]
