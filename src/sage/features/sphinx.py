r"""
Check for Sphinx
"""
from . import PythonModule


class Sphinx(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of Sphinx.

    Sphinx is provided by a standard package in the Sage distribution,
    but it can be disabled by ``configure --disable-doc``.

    EXAMPLES::

        sage: from sage.features.sphinx import Sphinx
        sage: Sphinx().is_present()                     # optional - sphinx
        FeatureTestResult('sphinx', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sphinx import Sphinx
            sage: isinstance(Sphinx(), Sphinx)
            True
        """
        PythonModule.__init__(self, 'sphinx', spkg='sphinx')


def all_features():
    return [Sphinx()]
