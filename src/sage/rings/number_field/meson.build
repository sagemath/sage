py.install_sources(
    'S_unit_solver.py',
    'all.py',
    'bdd_height.py',
    'class_group.py',
    'galois_group.py',
    'homset.py',
    'maps.py',
    'morphism.py',
    'number_field.py',
    'number_field_ideal.py',
    'number_field_ideal_rel.py',
    'number_field_rel.py',
    'order.py',
    'selmer_group.py',
    'small_primes_of_degree_one.py',
    'splitting_field.py',
    'structure.py',
    'totallyreal_phc.py',
    'totallyreal_rel.py',
    'unit_group.py',
    subdir: 'sage/rings/number_field',
)

extension_data_cpp = {
    'number_field_element': files('number_field_element.pyx'),
}

extension_data = {
    'number_field_base': files('number_field_base.pyx'),
    'number_field_element_base': files('number_field_element_base.pyx'),
    #'number_field_element_quadratic': files('number_field_element_quadratic.pyx'), - needs arb
    'number_field_morphisms': files('number_field_morphisms.pyx'),
    'totallyreal': files('totallyreal.pyx'),
    'totallyreal_data': files('totallyreal_data.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/number_field',
        install: true,
        include_directories: [inc_cpython, inc_ntl],
        dependencies: [py_dep, cysignals, cypari2],
    )
endforeach
foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/number_field',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_ntl],
        dependencies: [py_dep, cysignals, cypari2],
    )
endforeach
