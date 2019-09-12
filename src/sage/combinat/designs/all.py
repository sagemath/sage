"""
Combinatorial design features that are imported by default in the interpreter namespace

Test for deprecations of imports into global namespace::

    sage: designs_from_XML
    doctest:warning...:
    DeprecationWarning: 
    Importing designs_from_XML from here is deprecated. If you need to use it, please import it directly from sage.combinat.designs.ext_rep
    See https://trac.sagemath.org/27066 for details.
    ...

    sage: designs_from_XML_url
    doctest:warning...:
    DeprecationWarning: 
    Importing designs_from_XML_url from here is deprecated. If you need to use it, please import it directly from sage.combinat.designs.ext_rep
    See https://trac.sagemath.org/27066 for details.
    ...
"""
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import("sage.combinat.designs.ext_rep", ['designs_from_XML',
                                              'designs_from_XML_url'],
            deprecation=27066)

lazy_import('sage.combinat.designs.block_design', 'BlockDesign')

lazy_import('sage.combinat.designs.incidence_structures', 'IncidenceStructure')

lazy_import('sage.combinat.designs.incidence_structures',
            'IncidenceStructure', as_='Hypergraph')

lazy_import('sage.combinat.designs.covering_design',
            ['CoveringDesign', 'schonheim', 'trivial_covering_design'])

from . import design_catalog as designs

del absolute_import
