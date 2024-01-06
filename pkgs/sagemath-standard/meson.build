project(
  'sagemath-standard',
  ['c', 'cpp', 'cython'],
  version: files('VERSION.txt'),
  meson_version: '>= 1.2.0',
)

# same list in setup.py
set_variable('distribution_', true)  # catch-all
set_variable('distribution_sagemath_objects', true)
set_variable('distribution_sagemath_brial', true)
set_variable('distribution_sagemath_categories', true)
set_variable('distribution_sagemath_combinat', true)
set_variable('distribution_sagemath_eclib', true)
set_variable('distribution_sagemath_environment', true)
set_variable('distribution_sagemath_flint', true)
set_variable('distribution_sagemath_gap', true)
set_variable('distribution_sagemath_giac', true)
set_variable('distribution_sagemath_glpk', true)
set_variable('distribution_sagemath_graphs', true)
set_variable('distribution_sagemath_groups', true)
set_variable('distribution_sagemath_homfly', true)
set_variable('distribution_sagemath_lcalc', true)
set_variable('distribution_sagemath_libbraiding', true)
set_variable('distribution_sagemath_libecm', true)
set_variable('distribution_sagemath_linbox', true)
set_variable('distribution_sagemath_modules', true)
set_variable('distribution_sagemath_mpmath', true)
set_variable('distribution_sagemath_ntl', true)
set_variable('distribution_sagemath_objects', true)
set_variable('distribution_sagemath_pari', true)
set_variable('distribution_sagemath_plot', true)
set_variable('distribution_sagemath_polyhedra', true)
set_variable('distribution_sagemath_repl', true)
set_variable('distribution_sagemath_schemes', true)
set_variable('distribution_sagemath_singular', true)
set_variable('distribution_sagemath_symbolics', true)

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
