py.install_sources(
  '__init__.py',
  'all.py',
  'ambient.py',
  'apply.pxd',
  'boundary.py',
  'element.py',
  'g1list.py',
  'ghlist.py',
  'hecke_operator.py',
  'manin_symbol.pxd',
  'manin_symbol_list.py',
  'modsym.py',
  'modular_symbols.py',
  'p1list.pxd',
  'p1list_nf.py',
  'relation_matrix.py',
  'space.py',
  'subspace.py',
  'tests.py',
  subdir: 'sage/modular/modsym',
)

extension_data = {
  'apply' : files('apply.pyx'),
  'heilbronn' : files('heilbronn.pyx'),
  'manin_symbol' : files('manin_symbol.pyx'),
  'p1list' : files('p1list.pyx'),
  'relation_matrix_pyx' : files('relation_matrix_pyx.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, flint, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modular/modsym',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_rings],
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
      install_dir: sage_install_dir / 'modular/modsym',
      install: true,
      copy: true,
    )
  endif
endforeach

