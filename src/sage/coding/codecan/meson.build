py.install_sources('all.py', 'codecan.pxd', subdir: 'sage/coding/codecan')

extension_data = {
  'autgroup_can_label' : files('autgroup_can_label.pyx'),
  'codecan' : files('codecan.pyx'),
}

foreach name, pyx : extension_data
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
    dependencies: [py_dep, cysignals, gap, gmp],
  )
endforeach

