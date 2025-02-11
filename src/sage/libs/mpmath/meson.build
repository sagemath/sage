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
  deps = [py_dep, cysignals, gmp, mpfr]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/mpmath',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings],
    dependencies: deps,
  )

  all_deps_found = true
  foreach dep : deps
    if not dep.found()
      all_deps_found = false
      break
    endif
  endforeach
  if not all_deps_found
    configure_file(
      input: not_available,
      output: name + '.py',
      install_tag: 'python-runtime',
      install_dir: sage_install_dir / 'libs/mpmath',
      install: true,
      copy: true,
    )
  endif
endforeach

