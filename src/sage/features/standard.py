# sage_setup: distribution = sagemath-environment
r"""
Check for various standard packages (for modularized distributions)

These features are provided by standard packages in the Sage distribution.
"""

# *****************************************************************************
#       Copyright (C) 2023 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import PythonModule
from .join_feature import JoinFeature


def all_features():
    return [PythonModule('cvxopt', spkg='pkg:pypi/cvxopt', type='standard'),
            PythonModule('fpylll', spkg='pkg:pypi/fpylll', type='standard'),
            JoinFeature('ipython', (PythonModule('IPython'),), spkg='pkg:pypi/ipython', type='standard'),
            JoinFeature('lrcalc_python', (PythonModule('lrcalc'),), spkg='pkg:pypi/lrcalc', type='standard'),
            PythonModule('mpmath', spkg='pkg:pypi/mpmath', type='standard'),
            PythonModule('networkx', spkg='pkg:pypi/networkx', type='standard'),
            PythonModule('numpy', spkg='pkg:pypi/numpy', type='standard'),
            PythonModule('pexpect', spkg='pkg:pypi/pexpect', type='standard'),
            JoinFeature('pillow', (PythonModule('PIL'),), spkg='pkg:pypi/pillow', type='standard'),
            JoinFeature('pplpy', (PythonModule('ppl'),), spkg='pkg:pypi/pplpy', type='standard'),
            PythonModule('primecountpy', spkg='pkg:pypi/primecountpy', type='standard'),
            PythonModule('ptyprocess', spkg='pkg:pypi/ptyprocess', type='standard'),
            PythonModule('pyparsing', spkg='pkg:pypi/pyparsing', type='standard'),
            PythonModule('requests', spkg='pkg:pypi/requests', type='standard'),
            PythonModule('rpy2', spkg='pkg:pypi/rpy2', type='standard'),
            PythonModule('scipy', spkg='pkg:pypi/scipy', type='standard'),
            PythonModule('sympy', spkg='pkg:pypi/sympy', type='standard')]
