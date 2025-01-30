py.install_sources(
  '__init__.py',
  'all.py',
  'all__sagemath_polyhedra.py',
  'interactive_simplex_method.py',
  'knapsack.py',
  'linear_functions.pxd',
  'linear_tensor.py',
  'linear_tensor_constraints.py',
  'linear_tensor_element.pxd',
  'mip.pxd',
  'optimize.py',
  'sdp.pxd',
  subdir: 'sage/numerical',
)

extension_data = {
  'gauss_legendre' : files('gauss_legendre.pyx'),
  'linear_functions' : files('linear_functions.pyx'),
  'linear_tensor_element' : files('linear_tensor_element.pyx'),
  'mip' : files('mip.pyx'),
  'sdp' : files('sdp.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cypari2, gmp, mpfr]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/numerical',
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
      install_dir: sage_install_dir / 'numerical',
      install: true,
      copy: true,
    )
  endif
endforeach

subdir('backends')
