"""
Oriented Matroids
"""
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import
lazy_import('sage.matroids.oriented_matroids.oriented_matroid', 'OrientedMatroid')
lazy_import('sage.matroids.oriented_matroids.abstract_oriented_matroid', 'AbstractOrientedMatroid')