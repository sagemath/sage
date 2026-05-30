py.install_sources(
  '__init__.py',
  'decl.pxd',
  'function.pxd',
  'function.pyx',
  'function_factory.py',
  'groebner_strategy.pxd',
  'groebner_strategy.pyx',
  'option.pyx',
  'polynomial.pxd',
  'polynomial.pyx',
  'ring.pxd',
  'ring.pyx',
  'singular.pxd',
  'singular.pyx',
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
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/singular',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ntl, inc_rings, inc_rings_finite],
    dependencies: [py_dep, cysignals, givaro, gmp, singular],
  )
endforeach

