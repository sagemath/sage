py.install_sources(
    '__init__.py',
    'all.py',
    subdir: 'sage/libs/mpmath',
)

extension_data = {
    'ext_impl': files('ext_impl.pyx'),
    'ext_libmp': files('ext_libmp.pyx'),
    'ext_main': files('ext_main.pyx'),
    'utils': files('utils.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/mpmath',
        install: true,
        include_directories: [inc_cpython],
        dependencies: [py_dep, cysignals, cypari2, gmp, mpfr],
    )
endforeach
