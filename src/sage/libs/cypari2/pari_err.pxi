from sage.misc.superseded import deprecation
deprecation(20213, '''the file "pari_err.pxi" is deprecated, use sig_on instead of pari_catch_sig_on''')

from cysignals.signals cimport sig_on as pari_catch_sig_on
from cysignals.signals cimport sig_str as pari_catch_sig_str
from cysignals.signals cimport sig_off as pari_catch_sig_off
