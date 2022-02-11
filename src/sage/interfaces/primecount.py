from sage.misc.superseded import deprecation
deprecation(32894, "the module sage.interfaces.primecount is deprecated - use primecountpy.primecount instead")
from sage.misc.lazy_import import lazy_import
lazy_import("primecountpy.primecount", ['phi', 'nth_prime', 'prime_pi', 'prime_pi_128'], 
    deprecation=(32894, "the module sage.interfaces.primecount is deprecated - use primecountpy.primecount instead"))

