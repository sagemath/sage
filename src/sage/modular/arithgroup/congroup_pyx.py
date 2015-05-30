"""
Deprecated congruence subgroups module

TESTS::

    sage: from sage.modular.arithgroup.congroup_pyx import degeneracy_coset_representatives_gamma0
    sage: len(degeneracy_coset_representatives_gamma0(13, 1, 1))
    doctest:...: DeprecationWarning:
    Importing degeneracy_coset_representatives_gamma0 from here is deprecated. If you need to use it, please import it directly from sage.modular.arithgroup.congroup
    See http://trac.sagemath.org/17824 for details.
    14
"""

from sage.misc.lazy_import import lazy_import
from sage.misc.superseded import deprecation
lazy_import('sage.modular.arithgroup.congroup',
    ('degeneracy_coset_representatives_gamma0',
    'degeneracy_coset_representatives_gamma1'),
    deprecation=17824)

def generators_helper(coset_reps, level, Mat2Z):
    deprecation(17824, "The congroup_pyx module is deprecated, use sage.modular.arithgroup.congroup.generators_helper instead and remove the last Mat2Z argument")
    from sage.modular.arithgroup.congroup import generators_helper
    return generators_helper(coset_reps, level)
