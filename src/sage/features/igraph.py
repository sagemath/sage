r"""
Check for igraph
"""
from . import PythonModule
from .join_feature import JoinFeature


class python_igraph(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of the
    Python package ``igraph``.

    EXAMPLES::

        sage: from sage.features.igraph import python_igraph
        sage: python_igraph().is_present()                    # optional - python_igraph
        FeatureTestResult('python_igraph', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.igraph import python_igraph
            sage: isinstance(python_igraph(), python_igraph)
            True
        """
        JoinFeature.__init__(self, 'python_igraph',
                             [PythonModule('igraph', spkg="python_igraph",
                                            url="http://igraph.org")])

def all_features():
    return [python_igraph()]
