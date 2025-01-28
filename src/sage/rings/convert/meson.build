py.install_sources(
  '__init__.py',
  'all.py',
  'mpfi.pxd',
  subdir: 'sage/rings/convert',
)

extension_data = {'mpfi' : files('mpfi.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/convert',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cypari2, gmp, gsl, mpfi, mpfr, pari],
  )
endforeach

