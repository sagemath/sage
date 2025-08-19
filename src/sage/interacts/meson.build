py.install_sources(
  '__init__.py',
  'algebra.py',
  'all.py',
  'calculus.py',
  'fractals.py',
  'geometry.py',
  'library.py',
  'statistics.py',
  subdir: 'sage/interacts',
)

extension_data = {'library_cython' : files('library_cython.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/interacts',
    install: true,
    include_directories: [],
    dependencies: [py_dep, gmp],
  )
endforeach

