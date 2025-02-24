py.install_sources(
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
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/groups/perm_gps',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings],
    dependencies: [py_dep, cysignals, gap, gmp],
  )
endforeach

subdir('partn_ref')
subdir('partn_ref2')
