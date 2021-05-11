"""
Combinatorial design features that are imported by default in the interpreter namespace
"""
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.designs.block_design', 'BlockDesign')

lazy_import('sage.combinat.designs.incidence_structures', 'IncidenceStructure')

lazy_import('sage.combinat.designs.incidence_structures',
            'IncidenceStructure', as_='Hypergraph')

lazy_import('sage.combinat.designs.covering_design',
            ['CoveringDesign', 'schonheim', 'trivial_covering_design'])

from . import design_catalog as designs
