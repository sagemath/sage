"""
Index of combinatorial theories to be used in flag algebras


The :obj:`algebras.flag_theories` object may be used to access the theories already implemented.

{INDEX_OF_FUNCTIONS}

To import these names into the global namespace, use::

    sage: from sage.algebras.flag_theories_catalog import *
"""

# ****************************************************************************
#       Copyright (C) 2023 LEVENTE BODNAR <bodnalev at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import as _lazy_import

_lazy_import('sage.algebras.flag_algebras', 'CombinatorialTheory', 'CombinatorialTheory', at_startup=True)

_lazy_import('sage.algebras.flag_algebras', 'GraphTheory', 'GraphTheory', at_startup=True)
_lazy_import('sage.algebras.flag_algebras', 'DiGraphTheory', 'DiGraphTheory', at_startup=True)
_lazy_import('sage.algebras.flag_algebras', 'ThreeGraphTheory', 'ThreeGraphTheory', at_startup=True)
_lazy_import('sage.algebras.flag_algebras', 'TournamentTheory', 'TournamentTheory', at_startup=True)
_lazy_import('sage.algebras.flag_algebras', 'PermutationTheory', 'PermutationTheory', at_startup=True)
_lazy_import('sage.algebras.flag_algebras', 'OEGraphTheory', 'OEGraphTheory', at_startup=True)
_lazy_import('sage.algebras.flag_algebras', 'OVGraphTheory', 'OVGraphTheory', at_startup=True)
_lazy_import('sage.algebras.flag_algebras', 'RamseyGraphTheory', 'RamseyGraphTheory', at_startup=True)

from sage.misc.rest_index_of_methods import gen_rest_table_index as _gen_rest_table_index
import sys as _sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=_gen_rest_table_index(_sys.modules[__name__], only_local_functions=False))