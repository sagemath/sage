py.install_sources(
    'all.py',
    '__init__.py',
    '_py2_random.py',
    subdir: 'sage/cpython',
)

py.extension_module(
    'debug.pyx',
    'getattr.pyx',
    'atexit.pyx',
    'wrapperdescr.pyx',
    'builtin_types.pyx',
    'cython_metaclass.pyx',
    'dict_del_by_value.pyx',
    'string.pyx',
    'type.pyx',
    subdir: 'sage/cpython',
    install: true,
    dependencies: py_dep,
)
