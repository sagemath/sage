py.install_sources(
  '__init__.py',
  'cryptominisat.py',
  'dimacs.py',
  'picosat.py',
  'sat_lp.py',
  'satsolver.pxd',
  subdir: 'sage/sat/solvers',
)

extension_data = {'satsolver' : files('satsolver.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/sat/solvers',
    install: true,
    include_directories: [],
    dependencies: [py_dep, gmp],
  )
endforeach

