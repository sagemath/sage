r"""
Feature for testing the presence of ``pynormaliz``
"""
from . import PythonModule
from .join_feature import JoinFeature


class PyNormaliz(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the
    Python package ``PyNormaliz``.

    EXAMPLES::

        sage: from sage.features.normaliz import PyNormaliz
        sage: PyNormaliz().is_present()                    # optional - pynormaliz
        FeatureTestResult('pynormaliz', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.normaliz import PyNormaliz
            sage: isinstance(PyNormaliz(), PyNormaliz)
            True
        """
        JoinFeature.__init__(self, 'pynormaliz',
                             [PythonModule('PyNormaliz', spkg="pynormaliz")])


def all_features():
    return [PyNormaliz()]
