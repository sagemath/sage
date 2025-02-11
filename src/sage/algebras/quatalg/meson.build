py.install_sources(
  '__init__.py',
  'all.py',
  'quaternion_algebra.py',
  'quaternion_algebra_element.pxd',
  subdir: 'sage/algebras/quatalg',
)

extension_data_cpp = {
  'quaternion_algebra_cython': files('quaternion_algebra_cython.pyx'),
  'quaternion_algebra_element': files('quaternion_algebra_element.pyx'),
}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, flint, gmp, m, ntl]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras/quatalg',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_ntl, inc_rings],
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
      install_dir: sage_install_dir / 'algebras/quatalg',
      install: true,
      copy: true,
    )
  endif
endforeach

