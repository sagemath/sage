# Once https://github.com/mesonbuild/meson/pull/12616 is released, we can simply do the following:
# lcalc = dependency('lcalc', version: '>= 2.0.0')
lcalc = dependency(
  'lcalc',
  version: '>= 2.0.0',
  required: false,
  allow_fallback: false,
)
if not lcalc.found()
  lcalc = subproject('lcalc').get_variable('lcalc_dep')
endif

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

