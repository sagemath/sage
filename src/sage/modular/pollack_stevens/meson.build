py.install_sources(
  '__init__.py',
  'all.py',
  'dist.pxd',
  'distributions.py',
  'fund_domain.py',
  'manin_map.py',
  'modsym.py',
  'padic_lseries.py',
  'sigma0.py',
  'space.py',
  subdir: 'sage/modular/pollack_stevens',
)

extension_data = {'dist' : files('dist.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modular/pollack_stevens',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, gmp],
  )
endforeach

