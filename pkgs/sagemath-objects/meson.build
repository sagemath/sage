project(
  'sagemath-objects',
  ['c', 'cpp', 'cython'],
  version: files('VERSION.txt'),
  meson_version: '>= 1.2.0',
)

set_variable('distribution_sagemath_objects', true)

# Python module
# https://mesonbuild.com/Python-module.html
py_module = import('python')
py = py_module.find_installation(pure: false)
py_dep = py.dependency()

# Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cython = meson.get_compiler('cython')

py.install_sources(
    'sage/cpython/__init__.py',
    'sage/cpython/_py2_random.py',
    'sage/cpython/all.py',
    'sage/cpython/cython_metaclass.h',
    'sage/cpython/cython_metaclass.pxd',
    'sage/cpython/dict_del_by_value.pxd',
    'sage/cpython/dict_internal.h',
    'sage/cpython/getattr.pxd',
    'sage/cpython/pycore_long.h',
    'sage/cpython/pycore_long.pxd',
    'sage/cpython/python_debug.h',
    'sage/cpython/python_debug.pxd',
    'sage/cpython/pyx_visit.h',
    'sage/cpython/string.pxd',
    'sage/cpython/string_impl.h',
    'sage/cpython/type.pxd',
    'sage/cpython/wrapperdescr.pxd',
    subdir: 'sage/cpython/',
)

extension_data = {
    'atexit': files('sage/cpython/atexit.pyx'),
    'builtin_types': files('sage/cpython/builtin_types.pyx'),
    'cython_metaclass': files('sage/cpython/cython_metaclass.pyx'),
    'debug': files('sage/cpython/debug.pyx'),
    'dict_del_by_value': files('sage/cpython/dict_del_by_value.pyx'),
    'getattr': files('sage/cpython/getattr.pyx'),
    'string': files('sage/cpython/string.pyx'),
    'type': files('sage/cpython/type.pyx'),
    'wrapperdescr': files('sage/cpython/wrapperdescr.pyx'),
}

if get_variable('distribution_sagemath_objects', false)
foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        install: true,
        subdir: 'sage/cpython/',
        include_directories: include_directories('sage/cpython/'),
        dependencies: [py_dep],
    )
endforeach
endif
