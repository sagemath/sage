py.install_sources(
  'MOLS_handbook_data.py',
  'all.py',
  'bibd.py',
  'block_design.py',
  'covering_array.py',
  'covering_design.py',
  'database.py',
  'design_catalog.py',
  'designs_pyx.pxd',
  'difference_family.py',
  'difference_matrices.py',
  'ext_rep.py',
  'group_divisible_designs.py',
  'incidence_structures.py',
  'latin_squares.py',
  'orthogonal_arrays.py',
  'orthogonal_arrays_build_recursive.py',
  'resolvable_bibd.py',
  'steiner_quadruple_systems.py',
  'twographs.py',
  subdir: 'sage/combinat/designs',
)

extension_data = {
  'designs_pyx' : files('designs_pyx.pyx'),
  'evenly_distributed_sets' : files('evenly_distributed_sets.pyx'),
  'gen_quadrangles_with_spread' : files('gen_quadrangles_with_spread.pyx'),
  'orthogonal_arrays_find_recursive' : files(
    'orthogonal_arrays_find_recursive.pyx',
  ),
  'subhypergraph_search' : files('subhypergraph_search.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/designs',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

