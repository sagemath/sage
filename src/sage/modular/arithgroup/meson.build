inc_farey = include_directories('.')
py.install_sources(
    'all.py',
    'arithgroup_generic.py',
    'arithgroup_perm.py',
    'congroup_gamma.py',
    'congroup_gamma0.py',
    'congroup_gamma1.py',
    'congroup_gammaH.py',
    'congroup_generic.py',
    'congroup_sl2z.py',
    'tests.py',
    subdir: 'sage/modular/arithgroup',
)

extension_data = {
    'arithgroup_element': files('arithgroup_element.pyx'),
    'congroup': files('congroup.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/modular/arithgroup',
        install: true,
        include_directories: [inc_cpython, inc_ext, inc_farey, inc_flint, inc_rings],
        dependencies: [py_dep, cysignals, flint, gmp],
    )
endforeach

extension_data_cpp = {
    'farey_symbol': files('farey_symbol.pyx', 'sl2z.cpp', 'farey.cpp'),
}

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/modular/arithgroup',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_ext, inc_farey, inc_flint, inc_rings],
        dependencies: [py_dep, cysignals, flint, gmp],
    )
endforeach

