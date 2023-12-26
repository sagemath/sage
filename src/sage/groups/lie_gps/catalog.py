# sage_setup: distribution = sagemath-groups
r"""
Type ``groups.lie.<tab>`` to access examples of Lie groups.
"""
from sage.misc.lazy_import import lazy_import as _lazy_import

# We use lazy import because the module depends on sage.manifolds and thus on sage.symbolic
_lazy_import('sage.groups.lie_gps.nilpotent_lie_group', 'NilpotentLieGroup', as_='Nilpotent')
