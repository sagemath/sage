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
