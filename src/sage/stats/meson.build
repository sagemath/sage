py.install_sources('all.py', 'basic_stats.py', 'r.py', subdir: 'sage/stats')

extension_data = {
  'intlist' : files('intlist.pyx'),
  'time_series' : files('time_series.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/stats',
    install: true,
    include_directories: [inc_cpython, inc_numpy],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

subdir('distributions')
subdir('hmm')
