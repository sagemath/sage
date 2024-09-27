project(
  'SageMath',
  ['c', 'cpp', 'cython'],
  version: files('src/VERSION.txt'),
  license: 'GPL v3',
)

# Python module
# https://mesonbuild.com/Python-module.html
py_module = import('python')
py = py_module.find_installation(pure: false)
py_dep = py.dependency()

# Additional targets
py_with_pytest = py_module.find_installation(
  required: false,
  modules: ['pytest'],
)
if py_with_pytest.found()
  test(
    'pytest',
    py_with_pytest,
    args: [
      '-m',
      'pytest',
      '-c',
      meson.current_source_dir() / 'tox.ini',
      '--doctest',
      meson.current_source_dir() / 'src' / 'sage' / 'categories',
    ],
    timeout: 0,
  )
endif

subdir('src')
