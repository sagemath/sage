from sage.misc.lazy_import import lazy_import

lazy_import('sage.algebras.group_algebra', 'GroupAlgebra')

from .algebra import Algebra
from .finite_dimensional_algebras.all import FiniteDimensionalAlgebra
from .clifford_algebra import CliffordAlgebra, ExteriorAlgebra
from .weyl_algebra import DifferentialWeylAlgebra
lazy_import('sage.algebras.octonion_algebra', 'OctonionAlgebra')

import sage.algebras.catalog as algebras
