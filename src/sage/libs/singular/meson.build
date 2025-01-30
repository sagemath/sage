py.install_sources(
  '__init__.py',
  'decl.pxd',
  'function.pxd',
  'function_factory.py',
  'groebner_strategy.pxd',
  'polynomial.pxd',
  'ring.pxd',
  'singular.pxd',
  'standard_options.py',
  subdir: 'sage/libs/singular',
)

extension_data_cpp = {
  'function': files('function.pyx'),
  'groebner_strategy': files('groebner_strategy.pyx'),
  'option': files('option.pyx'),
  'polynomial': files('polynomial.pyx'),
  'ring': files('ring.pyx'),
  'singular': files('singular.pyx'),
}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, cysignals, givaro, gmp, singular]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/singular',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ntl, inc_rings, inc_rings_finite],
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
      install_dir: sage_install_dir / 'libs/singular',
      install: true,
      copy: true,
    )
  endif
endforeach

