r"""
Homsets between simplicial complexes: deprecated

The current version is :mod:`sage.topology.simplicial_complex_homset`.
"""
from sage.misc.superseded import deprecated_function_alias
import sage.topology.simplicial_complex_homset

is_SimplicialComplexHomset = deprecated_function_alias(31925,
                               sage.topology.simplicial_complex_homset.is_SimplicialComplexHomset)
SimplicialComplexHomset = deprecated_function_alias(31925,
                            sage.topology.simplicial_complex_homset.SimplicialComplexHomset)

