from sage.misc.lazy_import import lazy_import

from .quatalg.all import *

from .modsym.all import *

from .modform.all import *

from .ssmod.all import *

from .abvar.all import *

from .dirichlet import (DirichletGroup,
                        kronecker_character, kronecker_character_upside_down,
                        trivial_character)

from .arithgroup.all import (Gamma0, Gamma1, GammaH, Gamma, SL2Z,
                             ArithmeticSubgroup_Permutation,
                             CongruenceSubgroup, FareySymbol)

from .cusps import Cusp, Cusps

lazy_import('sage.modular.dims', ('dimension_cusp_forms',
                                  'dimension_new_cusp_forms',
                                  'dimension_eis',
                                  'dimension_modular_forms',
                                  'sturm_bound'),
            deprecation=(32647, 'removed from main namespace'))

from .etaproducts import (EtaGroup, EtaProduct, EtaGroupElement,
                          AllCusps, CuspFamily)

lazy_import('sage.modular.multiple_zeta', ['Multizeta', 'Multizetas'])

from .overconvergent.all import *

from .local_comp.all import *

from .cusps_nf import NFCusp, NFCusps, Gamma0_NFCusps

from .btquotients.all import *

from .pollack_stevens.all import *

from .quasimodform.all import *
