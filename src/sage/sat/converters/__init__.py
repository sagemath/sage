# sage_setup: distribution = sagemath-combinat
from sage.misc.lazy_import import lazy_import

from .anf2cnf import ANF2CNFConverter

lazy_import('sage.rings.polynomial.pbori.cnf', 'CNFEncoder', as_='PolyBoRiCNFEncoder')
