py.install_sources(
  'all.py',
  'dft.py',
  'dwt.pxd',
  'fft.pxd',
  subdir: 'sage/calculus/transforms',
)

extension_data = {'dwt' : files('dwt.pyx'), 'fft' : files('fft.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/calculus/transforms',
    install: true,
    include_directories: [inc_gsl],
    dependencies: [py_dep, cysignals, gmp, gsl],
  )
endforeach

