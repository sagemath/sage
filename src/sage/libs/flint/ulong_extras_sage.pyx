from .types cimport n_factor_t
from .ulong_extras cimport n_factor_init, n_factor


def n_factor_to_list(unsigned long n, int proved):
    """
    A wrapper around ``n_factor``.

    EXAMPLES::

        sage: from sage.libs.flint.ulong_extras_sage import n_factor_to_list
        sage: n_factor_to_list(60, 20)
        [(2, 2), (3, 1), (5, 1)]
        sage: n_factor_to_list((10**6).next_prime() + 1, 0)
        [(2, 2), (53, 2), (89, 1)]
    """
    cdef n_factor_t f
    n_factor_init(&f)
    n_factor(&f, n, proved)
    return [(f.p[i], int(f.exp[i])) for i in range(f.num)]
