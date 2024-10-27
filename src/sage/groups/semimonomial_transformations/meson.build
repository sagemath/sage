py.install_sources(
  'all.py',
  'semimonomial_transformation.pxd',
  'semimonomial_transformation_group.py',
  subdir: 'sage/groups/semimonomial_transformations',
)

extension_data = {
  'semimonomial_transformation' : files('semimonomial_transformation.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/groups/semimonomial_transformations',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, gmp],
  )
endforeach

