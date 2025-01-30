py.install_sources(
  '__init__.py',
  'all.py',
  'refinement_generic.pxd',
  subdir: 'sage/groups/perm_gps/partn_ref2',
)

extension_data = {'refinement_generic' : files('refinement_generic.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp, gap]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/groups/perm_gps/partn_ref2',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_partn_ref2],
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
      install_dir: sage_install_dir / 'groups/perm_gps/partn_ref2',
      install: true,
      copy: true,
    )
  endif
endforeach

