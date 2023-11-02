py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/gsl',
)

extension_data = {
    'array': files('array.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/gsl',
        install: true,
        include_directories: [inc_gsl],
        dependencies: [py_dep, cysignals, gmp, gsl],
    )
endforeach

