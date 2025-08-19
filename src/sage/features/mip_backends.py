# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of :class:`MixedIntegerLinearProgram` backends
"""

# *****************************************************************************
#       Copyright (C) 2021-2022 Matthias Koeppe
#                     2021      Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import Feature, PythonModule, FeatureTestResult
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


class CVXOPT(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing whether the :class:`MixedIntegerLinearProgram` backend ``CVXOPT`` is available.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.mip_backends import CVXOPT
            sage: CVXOPT()._is_present()  # optional - cvxopt
            FeatureTestResult('cvxopt', True)
        """
        JoinFeature.__init__(self, 'cvxopt',
                             [MIPBackend('CVXOPT'),
                              PythonModule('cvxopt')],
                             spkg='cvxopt',
                             type='standard')


def all_features():
    return [CPLEX(),
            Gurobi(),
            COIN(),
            CVXOPT()]
