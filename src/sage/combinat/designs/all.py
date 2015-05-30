"""
Combinatorial design features that are imported by default in the interpreter namespace
"""
from block_design import (BlockDesign)

from ext_rep import (designs_from_XML, designs_from_XML_url)

from incidence_structures import (IncidenceStructure,
                          IncidenceStructureFromMatrix)

from incidence_structures import IncidenceStructure as Hypergraph

from covering_design import (CoveringDesign,
                             schonheim,
                             trivial_covering_design)

import design_catalog as designs
