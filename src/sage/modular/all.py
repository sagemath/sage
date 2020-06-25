"""
Test for deprecations of imports into global namespace::

    sage: buzzard_tpslopes
    doctest:warning...:
    DeprecationWarning:
    Importing buzzard_tpslopes from here is deprecated. If you need to use it, please import it directly from sage.modular.buzzard
    See https://trac.sagemath.org/27066 for details.
    <function buzzard_tpslopes at ...>
"""
from __future__ import absolute_import
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

from .dims import (dimension_cusp_forms,
                   dimension_new_cusp_forms,
                   dimension_eis,
                   dimension_modular_forms,
                   sturm_bound)

lazy_import("sage.modular.buzzard", 'buzzard_tpslopes', deprecation=27066)

from .etaproducts import (EtaGroup, EtaProduct, EtaGroupElement,
                          AllCusps, CuspFamily)

lazy_import('sage.modular.multiple_zeta', ['Multizeta', 'Multizetas'])

from .overconvergent.all import *

from .local_comp.all import *

from .cusps_nf import NFCusp, NFCusps, Gamma0_NFCusps

from .btquotients.all import *

from .pollack_stevens.all import *

del absolute_import
