"""
Combinatorial design features that are imported by default in the interpreter namespace
"""
from __future__ import absolute_import

from .block_design import (BlockDesign)

from .ext_rep import (designs_from_XML, designs_from_XML_url)

from .incidence_structures import (IncidenceStructure)

from .incidence_structures import IncidenceStructure as Hypergraph

from .covering_design import (CoveringDesign,
                             schonheim,
                             trivial_covering_design)

from . import design_catalog as designs
