py.install_sources(
  '__init__.py',
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
  deps = [py_dep, gmp, singular]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras/letterplace',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
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
      install_dir: sage_install_dir / 'algebras/letterplace',
      install: true,
      copy: true,
    )
  endif
endforeach

