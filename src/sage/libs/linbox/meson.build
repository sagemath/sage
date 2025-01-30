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
  deps = [py_dep, blas, flint, fplll, givaro, gmp, gmpxx, linbox]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/linbox',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_flint],
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
      install_dir: sage_install_dir / 'libs/linbox',
      install: true,
      copy: true,
    )
  endif
endforeach

