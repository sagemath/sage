if cc.has_header('treedec/combinations.hpp', required: get_option('tdlib'))
  # tdlib is a header-only library
  tdlib = declare_dependency()
else
  tdlib = disabler()
endif
# Cannot be found via pkg-config
rw = cc.find_library('rw', required: get_option('rankwidth'), disabler: true)

py.install_sources(
  '__init__.py',
  'all.py',
  'bandwidth.pyx',
  'clique_separators.pyx',
  'cutwidth.pyx',
  'fast_digraph.pxd',
  'fast_digraph.pyx',
  'graph_products.pyx',
  'modular_decomposition.pxd',
  'modular_decomposition.pyx',
  'rankwidth.pxd',
  'rankwidth.pyx',
  'slice_decomposition.pxd',
  'slice_decomposition.pyx',
  'tdlib.pyx',
  'tree_decomposition.pxd',
  'tree_decomposition.pyx',
  'vertex_separation.pxd',
  'vertex_separation.pyx',
  subdir: 'sage/graphs/graph_decompositions',
)

extension_data = {
  'bandwidth': files('bandwidth.pyx'),
  'cutwidth': files('cutwidth.pyx'),
  'fast_digraph': files('fast_digraph.pyx'),
  'graph_products': files('graph_products.pyx'),
  'rankwidth': files('rankwidth.pyx'),
  'tree_decomposition': files('tree_decomposition.pyx'),
  'vertex_separation': files('vertex_separation.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp]
  if name == 'rankwidth'
    deps += [rw]
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/graph_decompositions',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: deps,
  )

  all_deps_found = true
  foreach dep : deps
    if not dep.found()
      all_deps_found = false
      break
    endif
  endforeach
  if not all_deps_found
    configure_file(
      input: not_available,
      output: name + '.py',
      install_tag: 'python-runtime',
      install_dir: sage_install_dir / 'graphs/graph_decompositions',
      install: true,
      copy: true,
    )
  endif
endforeach

extension_data_cpp = {
  'clique_separators': files('clique_separators.pyx'),
  'slice_decomposition': files('slice_decomposition.pyx'),
  'modular_decomposition': files('modular_decomposition.pyx'),
}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, cysignals, gmp, rw]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/graph_decompositions',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: deps,
  )

  all_deps_found = true
  foreach dep : deps
    if not dep.found()
      all_deps_found = false
      break
    endif
  endforeach
  if not all_deps_found
    configure_file(
      input: not_available,
      output: name + '.py',
      install_tag: 'python-runtime',
      install_dir: sage_install_dir / 'graphs/graph_decompositions',
      install: true,
      copy: true,
    )
  endif
endforeach

deps = [py_dep, cysignals, tdlib]
py.extension_module(
  'tdlib',
  sources: 'tdlib.pyx',
  subdir: 'sage/graphs/graph_decompositions',
  install: true,
  override_options: ['cython_language=cpp'],
  include_directories: [inc_cpython, inc_data_structures],
  dependencies: deps,
)
all_deps_found = true
foreach dep : deps
  if not dep.found()
    all_deps_found = false
    break
  endif
endforeach
if not all_deps_found
  configure_file(
    input: not_available,
    output: name + '.py',
    install_tag: 'python-runtime',
    install_dir: sage_install_dir / 'graphs/graph_decompositions',
    install: true,
    copy: true,
  )
endif
