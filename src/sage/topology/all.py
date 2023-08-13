# sage.doctest: needs sage.graphs
from .simplicial_complex import SimplicialComplex, Simplex

from .simplicial_complex_morphism import SimplicialComplexMorphism

from .delta_complex import DeltaComplex, delta_complexes

from .cubical_complex import CubicalComplex, cubical_complexes

from sage.misc.lazy_import import lazy_import
lazy_import('sage.topology.filtered_simplicial_complex', 'FilteredSimplicialComplex')

lazy_import('sage.topology', 'simplicial_complex_catalog', 'simplicial_complexes')
lazy_import('sage.topology', 'simplicial_set_catalog', 'simplicial_sets')


# # For taking care of old pickles
# from sage.misc.persist import register_unpickle_override
# register_unpickle_override('sage.topology.simplicial_complex_examples', 'SimplicialSurface', SimplicialComplex)
