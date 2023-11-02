py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/lcalc',
)

extension_data_cpp = {
    'lcalc_Lfunction': files('lcalc_Lfunction.pyx'),
}

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/lcalc',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_rings],
        dependencies: [py_dep, cypari2, cysignals, gmp, lfunction, m, mpfr, ntl],
    )
endforeach

