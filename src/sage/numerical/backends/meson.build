# Cannot be found via pkg-config
glpk = cc.find_library('glpk')

py.install_sources(
  '__init__.py',
  'all.py',
  'all__sagemath_polyhedra.py',
  'cvxopt_backend_test.py',
  'cvxpy_backend.pxd',
  'cvxpy_backend_test.py',
  'generic_backend.pxd',
  'generic_backend_test.py',
  'generic_sdp_backend.pxd',
  'glpk_backend.pxd',
  'glpk_backend_test.py',
  'glpk_exact_backend.pxd',
  'glpk_exact_backend_test.py',
  'glpk_graph_backend.pxd',
  'interactivelp_backend.pxd',
  'interactivelp_backend_test.py',
  'logging_backend.py',
  'matrix_sdp_backend.pxd',
  'ppl_backend_test.py',
  'scip_backend.pxd',
  'scip_backend_test.py',
  subdir: 'sage/numerical/backends',
)

extension_data = {
  'cvxopt_backend' : files('cvxopt_backend.pyx'),
  'cvxopt_sdp_backend' : files('cvxopt_sdp_backend.pyx'),
  'cvxpy_backend' : files('cvxpy_backend.pyx'),
  'generic_backend' : files('generic_backend.pyx'),
  'generic_sdp_backend' : files('generic_sdp_backend.pyx'),
  'glpk_backend' : files('glpk_backend.pyx'),
  'glpk_exact_backend' : files('glpk_exact_backend.pyx'),
  'glpk_graph_backend' : files('glpk_graph_backend.pyx'),
  'interactivelp_backend' : files('interactivelp_backend.pyx'),
  'matrix_sdp_backend' : files('matrix_sdp_backend.pyx'),
  'ppl_backend' : files('ppl_backend.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, glpk, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/numerical/backends',
    install: true,
    include_directories: [inc_cpython, inc_rings],
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
      install_dir: sage_install_dir / 'numerical/backends',
      install: true,
      copy: true,
    )
  endif
endforeach

extension_data_cpp = {'scip_backend': files('scip_backend.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/numerical/backends',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

