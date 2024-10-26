py.install_sources(
  'all.py',
  'all__sagemath_objects.py',
  'cartesian_product.py',
  'condition_set.py',
  'disjoint_set.pxd',
  'disjoint_union_enumerated_sets.py',
  'family.pxd',
  'finite_enumerated_set.py',
  'finite_set_map_cy.pxd',
  'finite_set_maps.py',
  'image_set.py',
  'integer_range.py',
  'non_negative_integers.py',
  'positive_integers.py',
  'primes.py',
  'pythonclass.pxd',
  'real_set.py',
  'recursively_enumerated_set.pxd',
  'set.py',
  'set_from_iterator.py',
  'totally_ordered_finite_set.py',
  subdir: 'sage/sets',
)

extension_data = {
  'disjoint_set' : files('disjoint_set.pyx'),
  'family' : files('family.pyx'),
  'finite_set_map_cy' : files('finite_set_map_cy.pyx'),
  'pythonclass' : files('pythonclass.pyx'),
  'recursively_enumerated_set' : files('recursively_enumerated_set.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/sets',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

