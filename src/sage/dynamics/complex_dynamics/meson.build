py.install_sources(
  '__init__.py',
  'all.py',
  'mandel_julia.py',
  subdir: 'sage/dynamics/complex_dynamics',
)

extension_data = {'mandel_julia_helper' : files('mandel_julia_helper.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/dynamics/complex_dynamics',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

