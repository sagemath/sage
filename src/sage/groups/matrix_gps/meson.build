py.install_sources(
  'all.py',
  'binary_dihedral.py',
  'catalog.py',
  'coxeter_group.py',
  'finitely_generated.py',
  'finitely_generated_gap.py',
  'group_element.pxd',
  'group_element_gap.pxd',
  'heisenberg.py',
  'homset.py',
  'isometries.py',
  'linear.py',
  'linear_gap.py',
  'matrix_group.py',
  'matrix_group_gap.py',
  'morphism.py',
  'named_group.py',
  'named_group_gap.py',
  'orthogonal.py',
  'orthogonal_gap.py',
  'pickling_overrides.py',
  'symplectic.py',
  'symplectic_gap.py',
  'unitary.py',
  'unitary_gap.py',
  subdir: 'sage/groups/matrix_gps',
)

extension_data = {
  'group_element' : files('group_element.pyx'),
  'group_element_gap' : files('group_element_gap.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, gmp]
  if name == 'group_element_gap'
    deps += [gap]
  endif    

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/groups/matrix_gps',
    install: true,
    include_directories: [inc_cpython],
    dependencies: deps,
  )
endforeach

