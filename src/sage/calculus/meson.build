py.install_sources(
  '__init__.py',
  'all.py',
  'calculus.py',
  'desolvers.py',
  'expr.py',
  'functional.py',
  'functions.py',
  'interpolation.pxd',
  'ode.pxd',
  'predefined.py',
  'test_sympy.py',
  'tests.py',
  'wester.py',
  subdir: 'sage/calculus',
)

extension_data = {
  'integration' : files('integration.pyx'),
  'interpolation' : files('interpolation.pyx'),
  'interpolators' : files('interpolators.pyx'),
  'ode' : files('ode.pyx'),
  'riemann' : files('riemann.pyx'),
  'var' : files('var.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp, gsl, interpreters_dep]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/calculus',
    install: true,
    include_directories: [inc_numpy],
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
      install_dir: sage_install_dir / 'calculus',
      install: true,
      copy: true,
    )
  endif
endforeach

subdir('transforms')
