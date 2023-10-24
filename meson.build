project(
  'SageMath',
  ['c', 'cpp', 'cython'],
  version: '10.2.beta8',
  license: 'GPL v3',
  meson_version: '>= 1.2.0',
)

# Python module
# https://mesonbuild.com/Python-module.html
py_module = import('python')
py = py_module.find_installation(pure: false)
py_dep = py.dependency()

# Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cython = meson.get_compiler('cython')

# Generate the configuration file
conf_data = configuration_data()

conf_data.set('PACKAGE_VERSION', '1.2.3')

maxima = find_program('maxima', required: true)
conf_data.set('SAGE_MAXIMA', maxima.path())
# Conda's ecl does not have any problems with Maxima, so nothing needs to be set here:
conf_data.set('SAGE_MAXIMA_FAS', '')

# Kenzo cannot yet be provided by the system, so we always use the SAGE_LOCAL path for now.
conf_data.set('SAGE_KENZO_FAS', '\'${prefix}\'/lib/ecl/kenzo.fas')

conf_data.set('SAGE_ARB_LIBRARY', 'arb')

#ntl = dependency('ntl', required: true) doesn't work, so ask the compiler directly:
ntl = cc.find_library('ntl', required: true)
# It can be found, so we don't have to set anything here:
conf_data.set('NTL_INCDIR', '')
conf_data.set('NTL_LIBDIR', '')

ecl_config = find_program('ecl-config', required: true)
conf_data.set('SAGE_ECL_CONFIG', ecl_config.path())

conf_data.set('SAGE_ARCHFLAGS', 'unset')

# not needed when using conda, as we then don't build any pc files
conf_data.set('SAGE_PKG_CONFIG_PATH', '')

openmp = dependency('openmp', required : false)
if openmp.found()
  conf_data.set('OPENMP_CFLAGS', '-fopenmp')
  conf_data.set('OPENMP_CXXFLAGS', '-fopenmp')
endif

configure_file(input : 'pkgs/sage-conf_conda/_sage_conf/_conf.py.in',
               output : '_conf.py',
               configuration : conf_data)

test(
    'import',
    py,
    args: [
      '-c',
      'from sage.categories.category import Category; print(Category)',
    ],
    )

# py_with_pytest = py_module.find_installation(required: false, modules: ['pytest'])
# if py_with_pytest.found()
#   test(
#     'pytest', 
#     py_with_pytest, 
#     args: [
#       '-m', 
#       'pytest', 
#       '-c', 
#       meson.current_source_dir() / 'tox.ini', 
#       '--doctest-modules', 
#       meson.current_source_dir() / 'src' / 'sage' / 'categories'
#     ],
#     timeout: 0,
#   )
# endif

subdir('src')
