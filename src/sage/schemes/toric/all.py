# sage.doctest: needs sage.geometry.polyhedron sage.graphs
from sage.misc.lazy_import import lazy_import

lazy_import('sage.schemes.toric.weierstrass', 'WeierstrassForm')
lazy_import('sage.schemes.toric.variety', ['AffineToricVariety', 'ToricVariety'])
lazy_import('sage.schemes.toric.library', 'toric_varieties')
lazy_import('sage.schemes.toric.fano_variety', ['FanoToricVariety', 'CPRFanoToricVariety'])
lazy_import('sage.schemes.toric.batyrev', 'SmoothFanoToricVariety')
lazy_import('sage.schemes.toric.batyrev_library', 'BTF')
lazy_import('sage.schemes.toric.ideal', 'ToricIdeal')
del lazy_import
