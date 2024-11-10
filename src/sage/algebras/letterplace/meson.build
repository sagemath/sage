py.install_sources(
  'all.py',
  'free_algebra_element_letterplace.pxd',
  'free_algebra_letterplace.pxd',
  subdir: 'sage/algebras/letterplace',
)

extension_data_cpp = {
  'free_algebra_element_letterplace': files(
    'free_algebra_element_letterplace.pyx',
  ),
  'free_algebra_letterplace': files('free_algebra_letterplace.pyx'),
  'letterplace_ideal': files('letterplace_ideal.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras/letterplace',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, gmp, singular],
  )
endforeach

