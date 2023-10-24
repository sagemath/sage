py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/gmp',
)

py.extension_module(
    'pylong.pyx',
    subdir: 'sage/libs/gmp',
    install: true,
    dependencies: py_dep,
)
