from . import PythonModule
from .join_feature import JoinFeature


class Mcqd(JoinFeature):

    def __init__(self):
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_mcqd' later
        JoinFeature.__init__(self, 'mcqd',
                             [PythonModule('sage.graphs.mcqd', spkg='mcqd')])
