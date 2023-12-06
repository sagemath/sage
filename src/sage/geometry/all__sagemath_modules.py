# sage_setup: distribution = sagemath-modules
from sage.misc.lazy_import import lazy_import

lazy_import('sage.geometry.toric_lattice', 'ToricLattice')
del lazy_import
