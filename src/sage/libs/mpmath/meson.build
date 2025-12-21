# Determine the mpmath version at build time so that we can support
# both 1.3 and 1.4 simultaneously (for a while).
mpmath14 = run_command(
  py, [
  '-c',
  'import mpmath; print(mpmath.__version__);'
  ],
  check: true,
).stdout().strip().startswith('1.4')


if mpmath14
  py.install_sources(
    '__init__.py',
    'all.py',
    'utils.pxd',
    'utils.pyx',
    subdir: 'sage/libs/mpmath',
  )

  extension_data = {
    'utils' : files('utils.pyx'),
  }

  foreach name, pyx : extension_data
    py.extension_module(
      name,
      sources: pyx,
      subdir: 'sage/libs/mpmath',
      install: true,
      include_directories: [inc_cpython, inc_ext, inc_rings],
      dependencies: [py_dep, cysignals, gmp, gmpy2, mpfr],
    )
  endforeach

else
  # mpmath-1.3.x
  py.install_sources(
    '__init__.py',
    'all.py',
    'ext_impl.pxd',
    'ext_impl.pyx',
    'ext_libmp.pyx',
    'ext_main.pxd',
    'ext_main.pyx',
    'utils.pxd',
    'utils13.pyx',
    subdir: 'sage/libs/mpmath',
  )

  extension_data = {
    'ext_impl' : files('ext_impl.pyx'),
    'ext_libmp' : files('ext_libmp.pyx'),
    'ext_main' : files('ext_main.pyx'),
    'utils' : files('utils13.pyx'),
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

endif
