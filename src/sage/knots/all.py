from sage.misc.lazy_import import lazy_import

lazy_import('sage.knots.knot', ['Knot', 'Knots'])
lazy_import('sage.knots.link', 'Link')
lazy_import('sage.knots.mosaic', ['Mosaic', 'MosaicTile', 'random_mosaic',
                                  'rational_tangle', 'tangle_join'])
lazy_import('sage.knots.knotinfo', ['KnotInfo', 'KnotInfoSeries'])
