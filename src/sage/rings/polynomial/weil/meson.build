py.install_sources(
  '__init__.py',
  'all.py',
  'power_sums.h',
  subdir: 'sage/rings/polynomial/weil',
)

extension_data = {'weil_polynomials' : files('weil_polynomials.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/polynomial/weil',
    install: true,
    include_directories: [inc_cpython, inc_flint, inc_rings, inc_src],
    dependencies: [py_dep, cysignals, flint, gmp],
  )
endforeach

