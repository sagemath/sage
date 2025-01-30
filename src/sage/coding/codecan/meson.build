py.install_sources(
  '__init__.py',
  'all.py',
  'codecan.pxd',
  subdir: 'sage/coding/codecan',
)

extension_data = {
  'autgroup_can_label' : files('autgroup_can_label.pyx'),
  'codecan' : files('codecan.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gap, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/coding/codecan',
    install: true,
    include_directories: [
      inc_cpython,
      inc_data_structures,
      inc_partn_ref2,
      inc_rings,
    ],
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
      install_dir: sage_install_dir / 'coding/codecan',
      install: true,
      copy: true,
    )
  endif
endforeach

