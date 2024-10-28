# sage_setup: distribution = sagemath-objects
from sage.structure.factorization import Factorization

from sage.structure.sequence import Sequence, seq

from sage.structure.unique_representation import UniqueRepresentation

from sage.structure.sage_object import SageObject

from sage.structure.element import (
    canonical_coercion,
    coercion_model,
    get_coercion_model,
    coercion_traceback,
    parent
)

from sage.structure.parent import Parent

from sage.structure.parent_gens import localvars

from sage.structure.proof import all as proof

from sage.misc.lazy_import import lazy_import
lazy_import('sage.structure.formal_sum', ['FormalSums', 'FormalSum'])
del lazy_import

from sage.structure.mutability import Mutability

from sage.structure.element_wrapper import ElementWrapper
