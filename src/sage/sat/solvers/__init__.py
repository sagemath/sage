from __future__ import absolute_import
from .satsolver import SatSolver

from .dimacs import Glucose, RSat

try:
    from .cryptominisat import CryptoMiniSat
except ImportError:
    pass
