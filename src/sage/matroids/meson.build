py.install_sources(
  'advanced.py',
  'all.py',
  'catalog.py',
  'constructor.py',
  'database_collections.py',
  'database_matroids.py',
  'dual_matroid.py',
  'matroids_catalog.py',
  'matroids_plot_helpers.py',
  'minor_matroid.py',
  'named_matroids.py',
  'rank_matroid.py',
  'utilities.py',
  subdir: 'sage/matroids',
)

extension_data = {
  'basis_exchange_matroid' : files('basis_exchange_matroid.pyx'),
  'basis_matroid' : files('basis_matroid.pyx'),
  'circuit_closures_matroid' : files('circuit_closures_matroid.pyx'),
  'circuits_matroid' : files('circuits_matroid.pyx'),
  'extension' : files('extension.pyx'),
  'lean_matrix' : files('lean_matrix.pyx'),
  'linear_matroid' : files('linear_matroid.pyx'),
  'matroid' : files('matroid.pyx'),
  'set_system' : files('set_system.pyx'),
  'union_matroid' : files('union_matroid.pyx'),
  'unpickling' : files('unpickling.pyx'),
  'flats_matroid' : files('flats_matroid.pyx'),
  'graphic_matroid' : files('graphic_matroid.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/matroids',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

