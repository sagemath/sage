# sage.doctest: needs sage.graphs
from sage.topology.simplicial_complex import SimplicialComplex, Simplex

from sage.topology.simplicial_complex_morphism import SimplicialComplexMorphism

from sage.topology.delta_complex import DeltaComplex, delta_complexes

from sage.topology.cubical_complex import CubicalComplex, cubical_complexes

from sage.misc.lazy_import import lazy_import
lazy_import('sage.topology.filtered_simplicial_complex', 'FilteredSimplicialComplex')

lazy_import('sage.topology', 'simplicial_complex_catalog', 'simplicial_complexes')
lazy_import('sage.topology', 'simplicial_set_catalog', 'simplicial_sets')

lazy_import('sage.topology.moment_angle_complex', 'MomentAngleComplex')

# # For taking care of old pickles
# from sage.misc.persist import register_unpickle_override
# register_unpickle_override('sage.topology.simplicial_complex_examples', 'SimplicialSurface', SimplicialComplex)
del lazy_import
