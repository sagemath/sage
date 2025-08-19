py.install_sources(
  '__init__.py',
  'all.py',
  'basic.py',
  'chessboard.py',
  'classical_geometries.py',
  'degree_sequence.py',
  'families.py',
  'intersection.py',
  'platonic_solids.py',
  'random.py',
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

