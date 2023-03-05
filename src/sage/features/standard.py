r"""
Check for various standard packages (for modularized distributions)

These features are provided by standard packages in the Sage distribution.
"""
from . import PythonModule
from . import JoinFeature


def all_features():
    return [PythonModule('cvxopt', spkg='cvxopt'),
            PythonModule('fpylll', spkg='fpylll'),
            PythonModule('IPython', spkg='ipython'),
            PythonModule('lrcalc', spkg='lrcalc_python'),
            PythonModule('mpmath', spkg='mpmath'),
            PythonModule('networkx', spkg='networkx'),
            PythonModule('numpy', spkg='numpy'),
            PythonModule('pexpect', spkg='pexpect'),
            PythonModule('PIL', spkg='pillow'),
            PythonModule('ppl', spkg='pplpy'),
            PythonModule('primecountpy', spkg='primecountpy'),
            PythonModule('ptyprocess', spkg='ptyprocess'),
            PythonModule('requests', spkg='requests'),
            PythonModule('rpy2', spkg='rpy2'),
            PythonModule('scipy', spkg='scipy'),
            PythonModule('sympy', spkg='sympy')]
