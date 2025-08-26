py.install_sources(
  '__init__.py',
  'all.py',
  'basic.py',
  'chessboard.py',
  'classical_geometries.py',
  'degree_sequence.py',
  'families.py',
  'generators_test.py',
  'intersection.py',
  'platonic_solids.py',
  'random.py',
  'shortened_000_111_extended_binary_Golay_code_graph.pickle.xz',
  'smallgraphs.py',
  'world_map.py',
  subdir: 'sage/graphs/generators',
)

extension_data = {'distance_regular' : files('distance_regular.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/generators',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

