py.install_sources(
  '__init__.py',
  'all.py',
  'boost_graph.pxd',
  'c_graph.pxd',
  'dense_graph.pxd',
  'graph_backends.pxd',
  'overview.py',
  'sparse_graph.pxd',
  'static_dense_graph.pxd',
  'static_sparse_backend.pxd',
  'static_sparse_graph.pxd',
  subdir: 'sage/graphs/base',
)

extension_data = {
  'dense_graph' : files('dense_graph.pyx'),
  'graph_backends' : files('graph_backends.pyx'),
  'sparse_graph' : files('sparse_graph.pyx'),
  'static_dense_graph' : files('static_dense_graph.pyx'),
  'static_sparse_backend' : files('static_sparse_backend.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/base',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

extension_data_cpp = {
  'boost_graph': files('boost_graph.pyx'),
  'c_graph': files('c_graph.pyx'),
  'static_sparse_graph': files('static_sparse_graph.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/base',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

