py.install_sources(
  'algebra.py',
  'algebra_elements.pxd',
  'all.py',
  'ar_quiver.py',
  'homspace.py',
  'morphism.py',
  'path_semigroup.py',
  'paths.pxd',
  'representation.py',
  subdir: 'sage/quivers',
)

extension_data = {
  'algebra_elements' : files('algebra_elements.pyx'),
  'paths' : files('paths.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/quivers',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

