py.install_sources(
  '__init__.py',
  'all.py',
  'basic.py',
  'chessboard.py',
  'classical_geometries.py',
  'degree_sequence.py',
  'distance_regular.pyx',
  'families.py',
  'generators_test.py',
  'intersection.py',
  'platonic_solids.py',
  'random.py',
  'shortened_000_111_extended_binary_Golay_code_graph.pickle.xz',
  'smallgraphs.py',
  'trees.pxd',
  'trees.pyx',
  'world_map.py',
  subdir: 'sage/graphs/generators',
)

extension_data = {
  'distance_regular' : files('distance_regular.pyx'),
  'trees' : files('trees.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/generators',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

