py.install_sources(
  '__init__.py',
  '_py2_random.py',
  'all.py',
  'cython_metaclass.h',
  'cython_metaclass.pxd',
  'dict_del_by_value.pxd',
  'dict_internal.h',
  'getattr.pxd',
  'pycore_long.h',
  'pycore_long.pxd',
  'python_debug.h',
  'python_debug.pxd',
  'pyx_visit.h',
  'string.pxd',
  'string_impl.h',
  'type.pxd',
  'wrapperdescr.pxd',
  subdir: 'sage/cpython',
)

extension_data = {
  'atexit' : files('atexit.pyx'),
  'builtin_types' : files('builtin_types.pyx'),
  'cython_metaclass' : files('cython_metaclass.pyx'),
  'debug' : files('debug.pyx'),
  'dict_del_by_value' : files('dict_del_by_value.pyx'),
  'getattr' : files('getattr.pyx'),
  'string' : files('string.pyx'),
  'type' : files('type.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/cpython',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, gmp],
  )
endforeach

