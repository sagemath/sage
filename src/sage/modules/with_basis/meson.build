py.install_sources(
  'all.py',
  'cell_module.py',
  'invariant.py',
  'morphism.py',
  'representation.py',
  'subquotient.py',
  subdir: 'sage/modules/with_basis',
)

extension_data = {'indexed_element' : files('indexed_element.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modules/with_basis',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, gmp],
  )
endforeach

