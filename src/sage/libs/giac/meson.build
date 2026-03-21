py.install_sources(
  '__init__.py',
  'auto-methods.pxi',
  'context.py',
  'gb.py',
  'giac.pxd',
  'giac.pyx',
  'keywords.pxi',
  'misc.h',
  subdir: 'sage/libs/giac',
)

py.extension_module(
  'giac',
  sources: files('giac.pyx'),
  subdir: 'sage/libs/giac',
  install: true,
  override_options: ['cython_language=cpp'],
  include_directories: [inc_cpython, inc_rings],
  dependencies: [py_dep, cysignals, giac_dep, gmpy2, gmpxx, gmp],
)
