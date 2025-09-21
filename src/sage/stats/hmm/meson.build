py.install_sources(
  '__init__.py',
  'all.py',
  'distributions.pxd',
  'hmm.pxd',
  'util.pxd',
  subdir: 'sage/stats/hmm',
)

extension_data = {
  'chmm' : files('chmm.pyx'),
  'distributions' : files('distributions.pyx'),
  'hmm' : files('hmm.pyx'),
  'util' : files('util.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/stats/hmm',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

