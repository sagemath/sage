py.install_sources(
    'all.py',
    'all__sagemath_objects.py',
    'cartesian_product.py',
    'condition_set.py',
    'disjoint_union_enumerated_sets.py',
    'family.py',
    'finite_enumerated_set.py',
    'finite_set_maps.py',
    'image_set.py',
    'integer_range.py',
    'non_negative_integers.py',
    'positive_integers.py',
    'primes.py',
    'real_set.py',
    'set.py',
    'set_from_iterator.py',
    'totally_ordered_finite_set.py',
    subdir: 'sage/sets',
)

extension_data = {
    #'disjoint_set': files('disjoint_set.pyx'), Needs other imports
    'finite_set_map_cy': files('finite_set_map_cy.pyx'),
    'pythonclass': files('pythonclass.pyx'),
    'recursively_enumerated_set': files('recursively_enumerated_set.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/sets',
        install: true,
        include_directories: [inc_cpython],
        dependencies: [py_dep, cysignals],
    )
endforeach
