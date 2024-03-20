project(
  'sagemath-mpmath',
  ['c', 'cpp', 'cython'],
  version: files('VERSION.txt'),
  meson_version: '>= 1.2.0',
)

set_variable('distribution_sagemath_mpmath', true)

# Python module
# https://mesonbuild.com/Python-module.html
py_module = import('python')
py = py_module.find_installation(pure: false)
py_dep = py.dependency()

# Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cython = meson.get_compiler('cython')

subdir('sage')
