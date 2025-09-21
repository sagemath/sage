py.install_sources(
  '__init__.py',
  'all.py',
  'utils.pxd',
  'utils.pyx',
  subdir: 'sage/libs/mpmath',
)

extension_data = {
  'utils' : files('utils.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/mpmath',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings],
    dependencies: [py_dep, cysignals, gmp, gmpy2, mpfr],
  )
endforeach

