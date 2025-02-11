py.install_sources(
  '__init__.py',
  'all.py',
  'cell_module.py',
  'indexed_element.pxd',
  'invariant.py',
  'morphism.py',
  'representation.py',
  'subquotient.py',
  subdir: 'sage/modules/with_basis',
)

extension_data = {'indexed_element' : files('indexed_element.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modules/with_basis',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
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
      install_dir: sage_install_dir / 'modules/with_basis',
      install: true,
      copy: true,
    )
  endif
endforeach

