r"""
Features for testing the presence of :class:`MixedIntegerLinearProgram` backends
"""

from . import Feature, FeatureTestResult
from .join_feature import JoinFeature


class MIPBackend(Feature):
    r"""
    A :class:`~sage.features.Feature` describing whether a :class:`MixedIntegerLinearProgram` backend is available.
    """
    def _is_present(self):
        r"""
        Test for the presence of a :class:`MixedIntegerLinearProgram` backend.

        EXAMPLES::

            sage: from sage.features.mip_backends import CPLEX
            sage: CPLEX()._is_present()  # optional - cplex
            FeatureTestResult('cplex', True)
        """
        try:
            from sage.numerical.mip import MixedIntegerLinearProgram
            MixedIntegerLinearProgram(solver=self.name)
            return FeatureTestResult(self, True)
        except Exception:
            return FeatureTestResult(self, False)


class CPLEX(MIPBackend):
    r"""
    A :class:`~sage.features.Feature` describing whether the :class:`MixedIntegerLinearProgram` backend ``CPLEX`` is available.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.mip_backends import CPLEX
            sage: CPLEX()._is_present()  # optional - cplex
            FeatureTestResult('cplex', True)
        """
        MIPBackend.__init__(self, 'cplex',
                            spkg='sage_numerical_backends_cplex')


class Gurobi(MIPBackend):
    r"""
    A :class:`~sage.features.Feature` describing whether the :class:`MixedIntegerLinearProgram` backend ``Gurobi`` is available.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.mip_backends import Gurobi
            sage: Gurobi()._is_present()  # optional - gurobi
            FeatureTestResult('gurobi', True)
        """
        MIPBackend.__init__(self, 'gurobi',
                            spkg='sage_numerical_backends_gurobi')


class COIN(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether the :class:`MixedIntegerLinearProgram` backend ``COIN`` is available.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.mip_backends import COIN
            sage: COIN()._is_present()  # optional - sage_numerical_backends_coin
            FeatureTestResult('sage_numerical_backends_coin', True)
        """
        JoinFeature.__init__(self, 'sage_numerical_backends_coin',
                             [MIPBackend('coin')],
                             spkg='sage_numerical_backends_coin')


def all_features():
    return [CPLEX(),
            Gurobi(),
            COIN()]
