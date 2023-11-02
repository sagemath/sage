py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/glpk',
)

extension_data = {
    'error': files('error.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/glpk',
        install: true,
        include_directories: [inc_cpython],
        dependencies: [py_dep, cysignals, gmp],
    )
endforeach

