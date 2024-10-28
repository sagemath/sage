py.install_sources(
  'all.py',
  'refinement_generic.pxd',
  subdir: 'sage/groups/perm_gps/partn_ref2',
)

extension_data = {'refinement_generic' : files('refinement_generic.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/groups/perm_gps/partn_ref2',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_partn_ref2],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

