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
    return [PythonModule('cvxopt', spkg='cvxopt', type='standard'),
            PythonModule('fpylll', spkg='fpylll', type='standard'),
            JoinFeature('ipython', (PythonModule('IPython'),), spkg='ipython', type='standard'),
            JoinFeature('lrcalc_python', (PythonModule('lrcalc'),), spkg='lrcalc_python', type='standard'),
            JoinFeature('mpmath', (PythonModule('sage.libs.mpmath'),
                                   PythonModule('sage.libs.mpmath.all'),),
                        spkg='sagemath_mpmath', type='standard'),
            PythonModule('networkx', spkg='networkx', type='standard'),
            PythonModule('numpy', spkg='numpy', type='standard'),
            PythonModule('pexpect', spkg='pexpect', type='standard'),
            JoinFeature('pillow', (PythonModule('PIL'),), spkg='pillow', type='standard'),
            JoinFeature('pplpy', (PythonModule('ppl'),), spkg='pplpy', type='standard'),
            PythonModule('primecountpy', spkg='primecountpy', type='standard'),
            PythonModule('ptyprocess', spkg='ptyprocess', type='standard'),
            PythonModule('pyparsing', spkg='pyparsing', type='standard'),
            PythonModule('requests', spkg='requests', type='standard'),
            PythonModule('rpy2', spkg='rpy2', type='standard'),
            PythonModule('scipy', spkg='scipy', type='standard'),
            PythonModule('sympy', spkg='sympy', type='standard')]
