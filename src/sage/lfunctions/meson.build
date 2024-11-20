py.install_sources(
  'all.py',
  'dokchitser.py',
  'lcalc.py',
  'pari.py',
  'sympow.py',
  subdir: 'sage/lfunctions',
)

extension_data = {'zero_sums' : files('zero_sums.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/lfunctions',
    install: true,
    include_directories: [inc_flint],
    dependencies: [py_dep, flint, gmp],
  )
endforeach

