py.install_sources(
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
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/calculus',
    install: true,
    include_directories: [inc_numpy],
    dependencies: [py_dep, cysignals, gmp, gsl, interpreters_dep],
  )
endforeach

subdir('transforms')
