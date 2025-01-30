py.install_sources(
  '__init__.py',
  'all.py',
  'constructor.py',
  'cubegroup.py',
  'permgroup.py',
  'permgroup_element.pxd',
  'permgroup_morphism.py',
  'permgroup_named.py',
  'permutation_groups_catalog.py',
  'symgp_conjugacy_class.py',
  subdir: 'sage/groups/perm_gps',
)

extension_data = {'permgroup_element' : files('permgroup_element.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gap, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/groups/perm_gps',
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
      install_dir: sage_install_dir / 'groups/perm_gps',
      install: true,
      copy: true,
    )
  endif
endforeach

subdir('partn_ref')
subdir('partn_ref2')
