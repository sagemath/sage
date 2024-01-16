# sage_setup: distribution = sagemath-schemes
###########################################################
# Re-bindings for unpickling
#
# Because we have a large number of pickles which were
# created when Gamma0, Gamma1, etc., were still classes
# as opposed to functions, we can't unpickle these unless
# sage.modular.congroup.Gamma0, etc. are classes. So we
# re-bind these as such below -- note that these bindings
# are *DIFFERENT* than the bindings for Gamma0, etc. that
# get imported in all.py.
#
# See trac #5059
#
###########################################################

from sage.modular.arithgroup.congroup_gamma0 import Gamma0_class
from sage.modular.arithgroup.congroup_gamma1 import Gamma1_class
from sage.modular.arithgroup.congroup_gammaH import GammaH_class

Gamma0 = Gamma0_class
Gamma1 = Gamma1_class
GammaH = GammaH_class

# For unpickling.
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor
from sage.modular.arithgroup.congroup_gamma1 import Gamma1_constructor
from sage.modular.arithgroup.congroup_gammaH import GammaH_constructor
from sage.modular.arithgroup.congroup_sl2z import SL2Z_class, _SL2Z_ref
