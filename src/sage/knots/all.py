from sage.misc.lazy_import import lazy_import
from sage.features.databases import DatabaseKnotInfo

lazy_import('sage.knots.knot', ['Knot', 'Knots'])
lazy_import('sage.knots.link', 'Link')
lazy_import('sage.knots.knotinfo', ['KnotInfo', 'KnotInfoSeries'])
