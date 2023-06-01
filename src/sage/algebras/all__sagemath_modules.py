from sage.misc.lazy_import import lazy_import

lazy_import('sage.algebras.group_algebra', 'GroupAlgebra')

from .algebra import Algebra
from .finite_dimensional_algebras.all import FiniteDimensionalAlgebra
