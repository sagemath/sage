from . import PythonModule
from .join_feature import JoinFeature


class Tdlib(JoinFeature):

    def __init__(self):
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_tdlib' later
        JoinFeature.__init__(self, 'tdlib',
                             [PythonModule('sage.graphs.graph_decompositions.tdlib', spkg='tdlib')])
