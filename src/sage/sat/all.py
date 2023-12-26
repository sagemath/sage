# sage_setup: distribution = sagemath-combinat
from sage.misc.lazy_import import lazy_import
lazy_import('sage.sat.solvers.satsolver', 'SAT')
del lazy_import
