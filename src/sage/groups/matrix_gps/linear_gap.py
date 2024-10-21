# sage_setup: distribution = sagemath-gap

"""
Linear Groups with GAP
"""

from sage.groups.matrix_gps.linear import LinearMatrixGroup_generic
from sage.groups.matrix_gps.named_group_gap import NamedMatrixGroup_gap
from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap


class LinearMatrixGroup_gap(NamedMatrixGroup_gap, LinearMatrixGroup_generic, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The general or special linear group in GAP.

    TESTS:

    Check that :issue:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
        sage: G = GL(3,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """
    pass
