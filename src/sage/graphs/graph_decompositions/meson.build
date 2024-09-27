tdlib = cc.find_library('tdlib', required: false)
# Cannot be found via pkg-config
rw = cc.find_library('rw')

py.install_sources(
  'all.py',
  'all__sagemath_tdlib.py',
  'fast_digraph.pxd',
  'modular_decomposition.py',
  'rankwidth.pxd',
  'tree_decomposition.pxd',
  'vertex_separation.pxd',
  subdir: 'sage/graphs/graph_decompositions',
)

extension_data = {
  'bandwidth' : files('bandwidth.pyx'),
  'cutwidth' : files('cutwidth.pyx'),
  'fast_digraph' : files('fast_digraph.pyx'),
  'graph_products' : files('graph_products.pyx'),
  'rankwidth' : files('rankwidth.pyx'),
  'tree_decomposition' : files('tree_decomposition.pyx'),
  'vertex_separation' : files('vertex_separation.pyx'),
  'slice_decomposition' : files('slice_decomposition.pyx'),
  'tdlib' : files('tdlib.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/graph_decompositions',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, cysignals, gmp, rw],
  )
endforeach

extension_data_cpp = {'clique_separators': files('clique_separators.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/graph_decompositions',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, cysignals, gmp, rw],
  )
endforeach

if tdlib.found()
  py.extension_module(
    'tdlib',
    sources: 'tdlib.pyx',
    subdir: 'sage/graphs/graph_decompositions',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, cysignals, tdlib],
  )
endif

