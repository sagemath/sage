py.install_sources(
  '__init__.py',
  'all.py',
  'convert_flint.pxd',
  'convert_gmp.pxd',
  'convert_sage.pxd',
  'convert_sage_complex_double.pxd',
  'convert_sage_real_double.pxd',
  'convert_sage_real_mpfr.pxd',
  'misc.pxd',
  'tests.py',
  subdir: 'sage/libs/pari',
)

extension_data = {
  'convert_flint' : files('convert_flint.pyx'),
  'convert_gmp' : files('convert_gmp.pyx'),
  'convert_sage' : files('convert_sage.pyx'),
  'convert_sage_complex_double' : files('convert_sage_complex_double.pyx'),
  'convert_sage_matrix' : files('convert_sage_matrix.pyx'),
  'convert_sage_real_double' : files('convert_sage_real_double.pyx'),
  'convert_sage_real_mpfr' : files('convert_sage_real_mpfr.pyx'),
  'misc' : files('misc.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/pari',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_gsl, inc_rings],
    dependencies: [py_dep, cypari2, cysignals, flint, gmp, gsl, mpfr, pari],
  )
endforeach

