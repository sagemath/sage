py.install_sources(
  'all.py',
  'finite_dimensional_algebra.py',
  'finite_dimensional_algebra_element.pxd',
  'finite_dimensional_algebra_ideal.py',
  'finite_dimensional_algebra_morphism.py',
  subdir: 'sage/algebras/finite_dimensional_algebras',
)

extension_data = {
  'finite_dimensional_algebra_element' : files(
    'finite_dimensional_algebra_element.pyx',
  ),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras/finite_dimensional_algebras',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, gmp],
  )
endforeach

