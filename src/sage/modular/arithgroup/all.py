# Note: the is_xxx functions are imported to here, but not from here up to sage.modular.all, so
# they are invisible to the user but easy to import all in one go by other code that needs them.

from sage.modular.arithgroup.arithgroup_generic import is_ArithmeticSubgroup, ArithmeticSubgroup
from sage.modular.arithgroup.congroup_generic import is_CongruenceSubgroup, CongruenceSubgroupBase, CongruenceSubgroup_constructor as CongruenceSubgroup
from sage.modular.arithgroup.congroup_gammaH import GammaH_constructor as GammaH, is_GammaH, GammaH_class
from sage.modular.arithgroup.congroup_gamma1 import Gamma1_constructor as Gamma1, is_Gamma1, Gamma1_class
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor as Gamma0, is_Gamma0, Gamma0_class
from sage.modular.arithgroup.congroup_gamma import Gamma_constructor as Gamma, is_Gamma, Gamma_class
from sage.modular.arithgroup.congroup_sl2z import SL2Z, is_SL2Z, SL2Z_class

from sage.misc.lazy_import import lazy_import
lazy_import('sage.modular.arithgroup.arithgroup_perm', 'ArithmeticSubgroup_Permutation')

from sage.modular.arithgroup.congroup import (degeneracy_coset_representatives_gamma0,
                                              degeneracy_coset_representatives_gamma1)

from sage.modular.arithgroup.farey_symbol import Farey as FareySymbol
del lazy_import
