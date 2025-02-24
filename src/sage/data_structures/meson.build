py.install_sources(
  'all.py',
  'binary_matrix.pxd',
  'binary_search.pxd',
  'bitset.pxd',
  'bitset_base.pxd',
  'bitset_intrinsics.h',
  'blas_dict.pxd',
  'bounded_integer_sequences.pxd',
  'list_of_pairs.pxd',
  'mutable_poset.py',
  'sparse_bitset.pxd',
  'stream.py',
  subdir: 'sage/data_structures',
)

extension_data = {
  'binary_search' : files('binary_search.pyx'),
  'bitset' : files('bitset.pyx'),
  'bitset_base' : files('bitset_base.pyx'),
  'blas_dict' : files('blas_dict.pyx'),
  'bounded_integer_sequences' : files('bounded_integer_sequences.pyx'),
  'list_of_pairs' : files('list_of_pairs.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/data_structures',
    install: true,
    include_directories: [
      inc_cpython,
      inc_data_structures,
      inc_flint,
      inc_rings,
    ],
    dependencies: [py_dep, cysignals, flint, gmp],
  )
endforeach

