py.install_sources(
  '__init__.py',
  'arb_wrap.h',
  'arith.pyx',
  subdir: 'sage/libs/arb',
)

extension_data = {'arith': files('arith.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/arb',
    install: true,
    include_directories: [inc_cpython, inc_flint, inc_rings],
    dependencies: [py_dep, flint, gmp, mpfr],
  )
endforeach

