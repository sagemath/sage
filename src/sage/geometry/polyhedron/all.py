# sage_setup: distribution = sagemath-polyhedra

from sage.misc.lazy_import import lazy_import
lazy_import('sage.geometry.polyhedron.constructor', 'Polyhedron')
lazy_import('sage.geometry.polyhedron.library', 'polytopes')
lazy_import('sage.geometry.polyhedron.combinatorial_polyhedron.base',
            'CombinatorialPolyhedron')
del lazy_import
