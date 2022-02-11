r"""
Morphisms of simplicial complexes: deprecated

The current version is :mod:`sage.topology.simplicial_complex_morphism`.
"""

from sage.misc.superseded import deprecated_function_alias
import sage.topology.simplicial_complex_morphism

is_SimplicialComplexMorphism = deprecated_function_alias(31925, 
                                 sage.topology.simplicial_complex_morphism.is_SimplicialComplexMorphism)
SimplicialComplexMorphism = deprecated_function_alias(31925,
                              sage.topology.simplicial_complex_morphism.SimplicialComplexMorphism)
