py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/gmp',
)

extension_data = {
    'pylong': files('pylong.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/gmp',
        install: true,
        include_directories: [inc_cpython],
        dependencies: py_dep,
    )
endforeach
