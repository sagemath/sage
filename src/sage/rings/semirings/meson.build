py.install_sources(
  '__init__.py',
  'all.py',
  'non_negative_integer_semiring.py',
  'tropical_mpolynomial.py',
  'tropical_polynomial.py',
  'tropical_variety.py',
  subdir: 'sage/rings/semirings',
)

extension_data = {'tropical_semiring' : files('tropical_semiring.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/semirings',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, gmp],
  )
endforeach

