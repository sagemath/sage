py.install_sources(
  '__init__.py',
  'all.py',
  'ext_impl.pxd',
  'ext_main.pxd',
  'utils.pxd',
  subdir: 'sage/libs/mpmath',
)

extension_data = {
  'ext_impl' : files('ext_impl.pyx'),
  'ext_libmp' : files('ext_libmp.pyx'),
  'ext_main' : files('ext_main.pyx'),
  'utils' : files('utils.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/mpmath',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings],
    dependencies: [py_dep, cysignals, gmp, mpfr],
  )
endforeach

