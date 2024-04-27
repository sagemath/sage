r"""
Combinatorial designs and incidence structures

All designs can be accessed by ``designs.<tab>`` and are listed in the
design catalog:

- :ref:`sage.combinat.designs.design_catalog`

**Design-related classes**

- :ref:`sage.combinat.designs.incidence_structures`
- :ref:`sage.combinat.designs.covering_design`

**Constructions**

- :ref:`sage.combinat.designs.block_design`
- :ref:`sage.combinat.designs.bibd`
- :ref:`sage.combinat.designs.resolvable_bibd`
- :ref:`sage.combinat.designs.group_divisible_designs`
- :ref:`sage.combinat.designs.latin_squares`
- :ref:`sage.combinat.designs.orthogonal_arrays`
- :ref:`sage.combinat.designs.orthogonal_arrays_build_recursive`
- :ref:`sage.combinat.designs.orthogonal_arrays_find_recursive`
- :ref:`sage.combinat.designs.difference_family`
- :ref:`sage.combinat.designs.difference_matrices`
- :ref:`sage.combinat.designs.steiner_quadruple_systems`
- :ref:`sage.combinat.designs.twographs`
- :ref:`sage.combinat.designs.database`
- :ref:`sage.combinat.designs.gen_quadrangles_with_spread`

**Technical things**

- :ref:`sage.combinat.designs.ext_rep`
- :ref:`sage.combinat.designs.designs_pyx`
- :ref:`sage.combinat.designs.subhypergraph_search`
- :ref:`sage.combinat.designs.evenly_distributed_sets`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import


lazy_import('sage.combinat.designs.incidence_structures', 'IncidenceStructure')
lazy_import('sage.combinat.designs.incidence_structures', 'IncidenceStructure', 'BlockDesign')

lazy_import('sage.combinat.designs.incidence_structures',
            'IncidenceStructure', as_='Hypergraph')

lazy_import('sage.combinat.designs.covering_design',
            ['CoveringDesign', 'schonheim', 'trivial_covering_design'])

from . import design_catalog as designs
