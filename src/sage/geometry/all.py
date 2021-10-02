from sage.misc.lazy_import import lazy_import

from .polyhedron.all import *
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_interface', 'HyperbolicPlane')
from .polyhedral_complex import PolyhedralComplex

lazy_import('sage.geometry.cone', ['Cone', 'random_cone'])
lazy_import('sage.geometry', 'cone_catalog', 'cones')
lazy_import('sage.geometry.fan', ['Fan', 'FaceFan', 'NormalFan', 'Fan2d'])
lazy_import('sage.geometry.fan_morphism', 'FanMorphism')
lazy_import('sage.geometry.lattice_polytope',
            ['LatticePolytope', 'NefPartition',
             'ReflexivePolytope', 'ReflexivePolytopes'])
lazy_import('sage.geometry', 'lattice_polytope')
lazy_import('sage.geometry.toric_lattice', 'ToricLattice')
lazy_import('sage.geometry', 'toric_plotter')
lazy_import('sage.geometry.voronoi_diagram', 'VoronoiDiagram')
lazy_import('sage.geometry.ribbon_graph', 'RibbonGraph')
lazy_import('sage.geometry.hyperplane_arrangement.arrangement', 'HyperplaneArrangements')
lazy_import('sage.geometry.hyperplane_arrangement.library', 'hyperplane_arrangements')
