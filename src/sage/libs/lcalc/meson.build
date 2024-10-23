lcalc = dependency('lcalc', version: '>= 2.0.0', required: not is_windows, disabler: true)

py.install_sources(
  '__init__.py',
  'lcalc_Lfunction.pxd',
  'lcalc_sage.h',
  subdir: 'sage/libs/lcalc',
)

extension_data_cpp = {'lcalc_Lfunction': files('lcalc_Lfunction.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/lcalc',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cypari2, cysignals, gmp, lcalc, m, mpfr, ntl],
  )
endforeach

