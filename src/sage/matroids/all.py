"""
Matroids
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import
lazy_import('sage.matroids.matroid', 'Matroid')
lazy_import('sage.matroids', 'matroids_catalog', 'matroids')
del lazy_import
del install_doc
