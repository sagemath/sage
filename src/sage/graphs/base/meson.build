boost_cpp17_compatible = false
boost = dependency(
  'boost',
  version: '>=1.83.0',
  required: false,
  modules: ['graph', 'tuple'],
)
if boost.found()
  # Need at least v1.83.0 (https://github.com/boostorg/functional/commit/6a573e4b8333ee63ee62ce95558c3667348db233)
  # for C++17 compatibility.
  boost_cpp17_compatible = true
else
  boost = dependency(
    'boost',
    version: '>=1.66.0',
    required: false,
    modules: ['graph', 'tuple'],
  )
endif

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
  deps = [py_dep, cysignals, gmp]
  override_options = ['cython_language=cpp']
  if name == 'boost_graph'
    deps += [boost]
    if not boost_cpp17_compatible
      override_options += ['cpp_std=c++11']
    endif
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/base',
    install: true,
    override_options: override_options,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: deps,
  )
endforeach

