from . import PythonModule
from .join_feature import JoinFeature


class Meataxe(JoinFeature):

    def __init__(self):
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_meataxe' later
        JoinFeature.__init__(self, 'meataxe',
                             [PythonModule('sage.matrix.matrix_gfpn_dense', spkg='meataxe')])
