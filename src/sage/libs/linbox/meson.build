py.install_sources(
  '__init__.py',
  'conversion.pxd',
  'fflas.pxd',
  'givaro.pxd',
  'linbox.pxd',
  'linbox_flint_interface.pxd',
  subdir: 'sage/libs/linbox',
)

extension_data_cpp = {
  'linbox_flint_interface': files('linbox_flint_interface.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/linbox',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_flint],
    dependencies: [py_dep, blas, flint, fplll, givaro, gmp, gmpxx, linbox],
  )
endforeach

