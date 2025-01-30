py.install_sources(
  '__init__.py',
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
  'pairing_heap.h',
  'pairing_heap.pxd',
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
  deps = [py_dep, cysignals, flint, gmp]
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
      install_dir: sage_install_dir / 'data_structures',
      install: true,
      copy: true,
    )
  endif
endforeach

extension_data_cpp = {'pairing_heap' : files('pairing_heap.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/data_structures/',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach
